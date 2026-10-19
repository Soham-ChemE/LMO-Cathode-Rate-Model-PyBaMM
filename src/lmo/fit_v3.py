"""v3: the v2.1 double-layer model plus one parasitic oxidation reaction at the positive
electrode, fitted to the v2.1 training data plus the measured excess charge.

v2.1 matches the rests that start at 3.5 V and misses the ones that start at the 4.5 V charge
cut-off, where the cell also loses charge to a side reaction (first-cycle coulombic efficiency
91 % at C/10). v3 tests whether one anodic Tafel reaction at the LMO surface
(:mod:`lmo.oxidation`) explains the after-charge rests. Two parameters are added to the four of
v2.1: j_ref, the parasitic current density at 4.5 V, and alpha_ox, its transfer coefficient
(Tafel slope 59 mV / alpha_ox per decade). The training data are exactly those of v2.1 (C/5 and
C/2 cycle-1 curves, the six C/10, C/5 and C/2 cycle-1 rests at lambda = 1) plus two constraint
residuals: the model's charge-minus-discharge capacity of C/10 cycle 1 and of C/5 cycle 1
against the measured excess charge. 2C is never used. SPMe first (three starting points), then
the DFN from the SPMe optimum.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import time

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

from . import PROCESSED_DIR, RESULTS_DIR
from .fit import BOUNDS_LOG10, PARAM_NAMES, RateSimulator, covariance_report
from .fit_v2 import HEADLINE_LAMBDA, REST_KINDS, REST_TRAINING, TRAINING, V2Residual, block_label, rest_jumps, rest_key
from .fit_v21 import EXTRA_BOUNDS_LOG10 as V21_BOUNDS
from .fit_v21 import VARIANTS as V21_VARIANTS
from .io import get_step
from .model import DEFAULT_DISC_MM, DEFAULT_FIT, PARTICLE_RADIUS_M, disc_geometry
from .oxidation import OPTION, PARAM_ALPHA, PARAM_J_REF, PARAM_U_REF, U_REF_V, oxidation_parameters

PREFIX = "v3_"
VARIANT = "dl+ox"
MODEL_OPTIONS = {**V21_VARIANTS["dl"]["options"], OPTION: "Tafel"}
EXTRA = ("C_dl", "j_ref", "alpha_ox")
BOUNDS = {**BOUNDS_LOG10, **V21_BOUNDS, "j_ref": (-6.0, 1.0), "alpha_ox": (-1.7, 0.5)}  # j_ref in A/m2 of particle surface
# The excess-charge constraint: cycle 1 of these rates, model (charge - discharge capacity)
# against measured. Weighted like the rests (lambda_q = 1: the two residuals together carry the
# summed squared weight of the curve block when each misses the measured excess by the
# tolerance), with the tolerance set by the data: the C/10 run read its current on a coarse
# range (+-2.5 %, i.e. +-0.04 mAh on a 0.14 mAh excess) and C/5 is good to a few percent.
EXCESS_TRAINING = [("C/10", 1), ("C/5", 1)]
EXCESS_TOLERANCE = 0.3  # relative
EXCESS_LAMBDA = 1.0
CURVE_SCALE_V = 0.025  # the curve RMSE scale that converts the relative excess error into volts
# Three starting points: the v2.1 SPMe optimum for Ds, k, R_contact and C_dl, with the side
# reaction started as a strong leak with a conventional Tafel slope, a weak one, and a leak
# that is almost potential-independent.
STARTS_SIDE = [{"j_ref": 2e-2, "alpha_ox": 0.5}, {"j_ref": 3e-3, "alpha_ox": 0.1}, {"j_ref": 1e-2, "alpha_ox": 0.03}]
DEFAULT_SIDE = {"C_dl": 3.0, "j_ref": 1e-2, "alpha_ox": 0.1}


def param_hook(params: pybamm.ParameterValues):
    params.update(
        {
            "Positive electrode double-layer capacity [F.m-2]": pybamm.InputParameter("C_dl"),
            **oxidation_parameters(pybamm.InputParameter("j_ref"), pybamm.InputParameter("alpha_ox")),
        },
        check_already_exists=False,
    )


def sim_kwargs() -> dict:
    return {"model_options": dict(MODEL_OPTIONS), "param_hook": param_hook, "extra_inputs": EXTRA}


def excess_key(rate: str) -> str:
    return f"{rate} excess charge"


class V3Residual(V2Residual):
    """The v2.1 residual (curves plus weighted rests) plus the excess-charge constraint blocks."""

    def __init__(self, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, lam=HEADLINE_LAMBDA):
        super().__init__(kind, disc_mm, lam, "trace", sim_kwargs=sim_kwargs(), param_names=list(PARAM_NAMES) + list(EXTRA))
        assert all(cyc == 1 and rate != "2C" and rate in self.sims for rate, cyc in EXCESS_TRAINING)
        self.excess_meas = {}
        for rate, cyc in EXCESS_TRAINING:
            chg, dis = get_step(rate, cyc, "charge"), get_step(rate, cyc, "discharge")
            self.excess_meas[rate] = chg.total_capacity_mAh - dis.total_capacity_mAh
        self.n_excess = len(EXCESS_TRAINING)
        self.w_excess = {rate: float(np.sqrt(EXCESS_LAMBDA * self.n_curve / self.n_excess) * CURVE_SCALE_V / (EXCESS_TOLERANCE * q)) for rate, q in self.excess_meas.items()}
        self.blocks = self.blocks + [("excess", rate, "cycle 1") for rate, _ in EXCESS_TRAINING]

    def block_sizes(self) -> dict[str, int]:
        return {**super().block_sizes(), **{excess_key(r): 1 for r in self.excess_meas}}

    def _n(self, b) -> int:
        return 1 if b[0] == "excess" else super()._n(b)

    def _weight(self, b) -> float:
        return {"curve": 1.0, "rest": self.w_rest}.get(b[0]) or self.w_excess[b[1]]

    def evaluate(self, fit: dict, blocks=None) -> dict[tuple, np.ndarray]:
        """Unweighted residuals: model - measured voltage (V) for curves and rests, model -
        measured excess charge (mAh) for the constraint blocks."""
        from .fit import model_on_measured_grid

        blocks = blocks or self.blocks
        rates = {b[1] for b in blocks}
        frames = {r: self.sims[r].run(fit, include_rests=True) for r in rates}  # one simulation per rate
        out = {}
        for b in blocks:
            typ, rate, k = b
            if typ == "curve":
                t, v = self.curve_points[(rate, k)]
                out[b] = model_on_measured_grid(frames[rate][k], t, k) - v
            elif typ == "rest":
                r = self.rests[(rate, k)]
                idx = self.rest_idx[(rate, k)]
                f = frames[rate][f"rest after {k}"]
                out[b] = np.interp(r.time_s[idx], f["time_s"].to_numpy(), f["voltage_V"].to_numpy()) - r.voltage_V[idx]
            else:
                model = frames[rate]["charge"]["capacity_mAh"].iloc[-1] - frames[rate]["discharge"]["capacity_mAh"].iloc[-1]
                out[b] = np.array([model - self.excess_meas[rate]])
        return out

    def weighted(self, theta, blocks=None) -> np.ndarray:
        blocks = blocks or self.blocks
        try:
            res = self.evaluate(self.unpack(theta), blocks)
        except (pybamm.SolverError, KeyError):
            return np.full(sum(self._n(b) for b in blocks), 1.0)
        return np.concatenate([res[b] * self._weight(b) for b in blocks])

    def rmse_by_type(self, fit: dict) -> dict:
        res = self.evaluate(fit)
        curve = np.concatenate([v for b, v in res.items() if b[0] == "curve"])
        rest = np.concatenate([v for b, v in res.items() if b[0] == "rest"])
        rest_chg = np.concatenate([v for b, v in res.items() if b[0] == "rest" and b[2] == "charge"])
        rest_dis = np.concatenate([v for b, v in res.items() if b[0] == "rest" and b[2] == "discharge"])
        rm = lambda x: float(1000 * np.sqrt(np.mean(x**2)))
        return {
            "curves_mV": rm(curve),
            "rests_mV": rm(rest),
            "rest_traces_mV": rm(rest),
            "rests_after_charge_mV": rm(rest_chg),
            "rests_after_discharge_mV": rm(rest_dis),
            "per_block_mV": {block_label(b): rm(v) for b, v in res.items() if b[0] != "excess"},
            "excess_charge_mAh": {b[1]: {"model": float(v[0] + self.excess_meas[b[1]]), "measured": self.excess_meas[b[1]], "error": float(v[0])} for b, v in res.items() if b[0] == "excess"},
        }


def rest_rmse_split(per_block_mV: dict, n_points: dict) -> dict:
    """After-charge and after-discharge rest RMSE from a per-block table (used to read the same
    numbers from the v2.1 result file)."""
    out = {}
    for after in REST_KINDS:
        keys = [rest_key(r, after) for r, _ in REST_TRAINING]
        ss = sum(per_block_mV[k] ** 2 * n_points[k] for k in keys)
        out[f"rests_after_{after}_mV"] = float(np.sqrt(ss / sum(n_points[k] for k in keys)))
    return out


def fit_v3(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, x0: dict | None = None, blocks=None, verbose=True, starts=None, res: V3Residual | None = None):
    res = res or V3Residual(kind, disc_mm)
    lb = np.array([BOUNDS[n][0] for n in res.param_names])
    ub = np.array([BOUNDS[n][1] for n in res.param_names])
    t0 = time.time()
    best, best_cost, tried = None, None, []
    for start in starts or [x0]:
        x = {**DEFAULT_FIT, **DEFAULT_SIDE, **(start or {})}
        theta0 = np.array([np.log10(x[n]) for n in res.param_names])
        r = least_squares(lambda th: res.weighted(th, blocks), theta0, bounds=(lb, ub), diff_step=1e-3, xtol=1e-5, ftol=1e-6, verbose=2 if verbose else 0)
        tried.append({"start": {n: x[n] for n in res.param_names}, "cost": float(np.sum(r.fun**2)), "optimum": {n: float(v) for n, v in res.unpack(r.x).items()}})
        if best is None or tried[-1]["cost"] < best_cost:
            best, best_cost = r, tried[-1]["cost"]
    best.starts_tried = tried
    return best, res, time.time() - t0


def _loo_refit(args):
    kind, disc_mm, best, keep, label = args
    r, res, _ = fit_v3(kind, disc_mm, x0=best, blocks=keep, verbose=False)
    d = res.unpack(r.x)
    d["dropped"] = label
    return d


def _block_label(b) -> str:
    return excess_key(b[1]) if b[0] == "excess" else block_label(b)


def v21_base(kind: str) -> dict | None:
    """The v2.1 headline parameters of the same model (the v3 starting point), if available."""
    for name in (f"v2_1_fit_{kind.lower()}.json", "v2_1_fit_spme.json"):
        p = RESULTS_DIR / name
        if p.exists():
            d = json.load(open(p))["parameters"]
            return {n: d[n]["value"] for n in ("Ds", "k", "R_contact", "C_dl")}
    return None


def starts_for(kind: str) -> list[dict]:
    """SPMe: three starts from the v2.1 optimum. DFN: the v3 SPMe optimum alone."""
    if kind != "SPMe":
        p = RESULTS_DIR / f"{PREFIX}fit_spme.json"
        if p.exists():
            d = json.load(open(p))["parameters"]
            return [{n: d[n]["value"] for n in list(PARAM_NAMES) + list(EXTRA)}]
    base = v21_base(kind) or {}
    return [{**base, **side} for side in STARTS_SIDE]


def leakage_current_table(fit: dict, disc_mm: float = DEFAULT_DISC_MM) -> dict:
    """The fitted parasitic current (per cell, mA) as a function of the particle surface
    potential, for the record and the figure."""
    from .model import COATING_THICKNESS_M, F, T_K

    geo = disc_geometry(disc_mm)
    area = 3 * geo["active_volume_fraction"] / PARTICLE_RADIUS_M * geo["area_m2"] * COATING_THICKNESS_M
    R = 8.314462618
    volts = [3.5, 3.8, 4.0, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8]
    i_mA = {f"{v:.1f}": float(1000 * fit["j_ref"] * area * np.exp(fit["alpha_ox"] * F / (R * T_K) * (v - U_REF_V))) for v in volts}
    return {
        "particle_surface_area_m2": area,
        "tafel_slope_mV_per_decade": float(1000 * np.log(10) * R * T_K / (F * fit["alpha_ox"])),
        "current_at_reference_potential_mA": i_mA[f"{U_REF_V:.1f}"],
        "current_mA_vs_surface_potential_V": i_mA,
    }


def coulombic_efficiency_check(fit: dict, kind: str, disc_mm: float = DEFAULT_DISC_MM, rates=("C/10", "C/5", "C/2")) -> dict:
    """Model vs measured charge, discharge and coulombic efficiency for both cycles of the
    low-rate runs (cycle 1 of C/10 and C/5 is training through the constraint; everything else
    is a check). One two-cycle simulation per rate."""
    out = {}
    for rate in rates:
        sim = RateSimulator(rate, 2, kind, disc_mm, mesh="v2", **sim_kwargs())
        sol = sim.sim.solve(inputs={n: float(fit[n]) for n in sim.input_names})
        from .model import solution_frames

        frames = solution_frames(sol, 2, rate)
        sim.sim._solution = None
        for cyc in (1, 2):
            fr = {f["step"].iloc[0]: f for f in frames if f["cycle"].iloc[0] == cyc}
            qc, qd = float(fr["charge"]["capacity_mAh"].iloc[-1]), float(fr["discharge"]["capacity_mAh"].iloc[-1])
            mc, md = get_step(rate, cyc, "charge"), get_step(rate, cyc, "discharge")
            out[f"{rate} cycle {cyc}"] = {
                "role": "training (excess-charge constraint)" if (rate, cyc) in EXCESS_TRAINING else "check",
                "charge_mAh_model": qc,
                "charge_mAh_measured": mc.total_capacity_mAh,
                "discharge_mAh_model": qd,
                "discharge_mAh_measured": md.total_capacity_mAh,
                "excess_mAh_model": qc - qd,
                "excess_mAh_measured": mc.total_capacity_mAh - md.total_capacity_mAh,
                "coulombic_efficiency_model": qd / qc,
                "coulombic_efficiency_measured": md.total_capacity_mAh / mc.total_capacity_mAh,
            }
        del sol
    return out


def summarise(result, res: V3Residual, elapsed: float, loo=None) -> dict:
    fit = res.unpack(result.x)
    rep = covariance_report(result)
    rm = res.rmse_by_type(fit)
    names = res.param_names
    out = {
        "version": "v3",
        "variant": VARIANT,
        "model_options": dict(MODEL_OPTIONS),
        "side_reaction": f"j_ox = j_ref exp(alpha_ox F (phi_s - phi_e - {U_REF_V}) / RT) on the positive particle surface",
        "model": res.kind,
        "disc_mm": res.disc_mm,
        "training_curves": [f"{r} cycle {c} charge+discharge" for r, c in TRAINING],
        "training_rests": [f"{r} cycle {c} rest after {a}" for r, c in REST_TRAINING for a in REST_KINDS],
        "training_excess_charge": [f"{r} cycle {c} charge minus discharge capacity" for r, c in EXCESS_TRAINING],
        "two_c_used": False,
        "particle_radius_m": PARTICLE_RADIUS_M,
        "weighting": {
            "lambda": res.lam,
            "rest_weight": res.w_rest,
            "rest_points": "every logged sample of each rest",
            "n_curve_points": res.n_curve,
            "n_rest_points": res.n_rest,
            "excess_lambda": EXCESS_LAMBDA,
            "excess_tolerance_relative": EXCESS_TOLERANCE,
            "excess_weight_per_mAh": res.w_excess,
            "excess_rule": "w = sqrt(lambda_q N_curve / N_q) * 25 mV / (tolerance * measured excess): the constraint block carries the curve block's weight when each excess misses by the tolerance",
        },
        "training_rmse_curves_V": rm["curves_mV"] / 1000,
        "training_rmse_rests_V": rm["rests_mV"] / 1000,
        "training_rmse_rest_traces_V": rm["rest_traces_mV"] / 1000,
        "training_rmse_rests_after_charge_V": rm["rests_after_charge_mV"] / 1000,
        "training_rmse_rests_after_discharge_V": rm["rests_after_discharge_mV"] / 1000,
        "training_rmse_per_block_mV": rm["per_block_mV"],
        "excess_charge_mAh": rm["excess_charge_mAh"],
        "weighted_rmse_V": rep["rmse_V"],
        "n_residuals": int(len(result.fun)),
        "fit_time_s": elapsed,
        "converged": bool(result.success),
        "starts_tried": getattr(result, "starts_tried", None),
        "parameters": {},
        "log10_correlation": {f"{a}-{b}": float(rep["corr"][i, j]) for i, a in enumerate(names) for j, b in enumerate(names) if i < j},
    }
    for i, n in enumerate(names):
        sd = rep["sd_log10"][i]
        entry = {"value": fit[n], "log10": float(result.x[i]), "sd_log10_jacobian": float(sd), "ci95_low_jacobian": float(10 ** (result.x[i] - 1.96 * sd)), "ci95_high_jacobian": float(10 ** (result.x[i] + 1.96 * sd))}
        if loo:
            vals = np.array([d[n] for d in loo])
            entry["leave_one_step_out_min"] = float(vals.min())
            entry["leave_one_step_out_max"] = float(vals.max())
        out["parameters"][n] = entry
    tau = PARTICLE_RADIUS_M**2 / fit["Ds"]
    sd = rep["sd_log10"][0]
    out["parameters"]["tau_R2_over_Ds_s"] = {"value": tau, "ci95_low_jacobian": float(tau / 10 ** (1.96 * sd)), "ci95_high_jacobian": float(tau * 10 ** (1.96 * sd))}
    if loo:
        taus = np.array([PARTICLE_RADIUS_M**2 / d["Ds"] for d in loo])
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_min"] = float(taus.min())
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_max"] = float(taus.max())
    out["parameters"]["R_contact"]["unit"] = "Ohm"
    out["parameters"]["Ds"]["unit"] = "m2/s"
    out["parameters"]["k"]["unit"] = "m^2.5/(mol^0.5 s) in j0 = k F sqrt(c_e c_s (c_max - c_s))"
    out["parameters"]["C_dl"]["unit"] = "F/m2 of particle surface (effective; absorbs the unknown real surface area)"
    out["parameters"]["j_ref"]["unit"] = f"A/m2 of particle surface at {U_REF_V} V (effective; absorbs the unknown real surface area)"
    out["parameters"]["alpha_ox"]["unit"] = "transfer coefficient of the oxidation; Tafel slope = 59 mV / alpha_ox per decade"
    out["leakage_current"] = leakage_current_table(fit, res.disc_mm)
    if loo:
        tab = [leakage_current_table({**fit, "j_ref": d["j_ref"], "alpha_ox": d["alpha_ox"]}, res.disc_mm) for d in loo]
        out["leakage_current"]["leave_one_out_current_at_reference_potential_mA"] = [float(t["current_at_reference_potential_mA"]) for t in tab]
        out["leakage_current"]["leave_one_out_tafel_slope_mV_per_decade"] = [float(t["tafel_slope_mV_per_decade"]) for t in tab]
    out["rest_jumps"] = rest_jumps(res, fit)
    return out


def fit_with_loo(kind: str, disc_mm: float = DEFAULT_DISC_MM, starts=None, leave_one_out: bool = True):
    result, res, elapsed = fit_v3(kind, disc_mm, starts=starts)
    best = res.unpack(result.x)
    jobs = [(kind, disc_mm, best, [b for b in res.blocks if b != drop], _block_label(drop)) for drop in res.blocks]
    workers = int(os.environ.get("LMO_WORKERS", "1"))
    if not leave_one_out:
        loo = []
    elif workers > 1:
        with multiprocessing.get_context("spawn").Pool(workers) as pool:
            loo = pool.map(_loo_refit, jobs)
    else:
        loo = [_loo_refit(j) for j in jobs]
    summary = summarise(result, res, elapsed, loo)
    summary["leave_one_step_out"] = loo
    return result, res, summary


def run(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, starts=None, leave_one_out: bool = True) -> dict:
    """Fit (SPMe: three starts; DFN: from the SPMe optimum) with leave-one-block-out refits,
    the coulombic-efficiency check on both cycles, and the files
    results/v3_fit_<model>.json and data/processed/v3_model_rests_<model>.csv."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    result, res, summary = fit_with_loo(kind, disc_mm, starts=starts or starts_for(kind), leave_one_out=leave_one_out)
    best = res.unpack(result.x)
    summary["coulombic_efficiency"] = coulombic_efficiency_check(best, kind, disc_mm)
    summary["v2_1_start"] = v21_base(kind)
    with open(RESULTS_DIR / f"{PREFIX}fit_{kind.lower()}.json", "w") as f:
        json.dump(summary, f, indent=2)
    res.rest_frames(best).to_csv(PROCESSED_DIR / f"{PREFIX}model_rests_{kind.lower()}.csv", index=False)
    return summary


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    out = run(kind)
    print(json.dumps({k: out[k] for k in ("parameters", "log10_correlation", "training_rmse_curves_V", "training_rmse_rests_V", "excess_charge_mAh", "coulombic_efficiency")}, indent=2))
