"""v2 fit: the same three parameters and model as v1, trained on the C/5 and C/2 curves plus
the 30 s rests that follow the C/10, C/5 and C/2 steps.

v1 fitted only the charge/discharge curves of C/5 and C/2, which cannot tell an instantaneous
ohmic drop from a transport loss that builds up over time (k vs R_contact correlation 0.94); it
put the loss into a 96.5 Ohm contact resistance and over-predicted the 2C capacity. The rest
after each step separates the two: the jump in the first tens of milliseconds is ohmic (plus
the fastest part of the charge-transfer relaxation), the relaxation over the following 30 s is
diffusion. The v2 residual therefore adds the measured rest traces (every logged point, time
measured from the last sample under current) to the v1 curve residual. No 2C data of any kind
(curves or rests) is used.

Weighting: a rest has 10 to 25 logged points against 300 per curve step, so unweighted it
would barely move the fit. Each rest residual is multiplied by w = sqrt(lambda * N_curve /
N_rest), so that the rest block carries ``lambda`` times the summed squared weight of the curve
block. lambda = 1 (equal block weight) on the full rest traces is the headline; the pipeline
also refits at the other entries of ``WEIGHTINGS`` as a sensitivity check. One of them,
"jump only", keeps only the first logged sample of each rest (the instantaneous jump): the
model has no double-layer capacitance and a single particle size, so the one rest feature it
can represent is that jump, and this variant is reported in full (prefix ``v2_jump_``).
"""

from __future__ import annotations

import json
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

from . import PROCESSED_DIR, RESULTS_DIR
from .fit import BOUNDS_LOG10, N_POINTS_PER_STEP, PARAM_NAMES, RateSimulator, TRAINING, covariance_report, model_on_measured_grid, unpack
from .io import Rest, get_rest
from .model import DEFAULT_DISC_MM, DEFAULT_FIT, PARTICLE_RADIUS_M

# Rests used for training: cycle 1 of every rate except 2C, after charge and after discharge.
REST_TRAINING = [("C/10", 1), ("C/5", 1), ("C/2", 1)]
REST_KINDS = ("charge", "discharge")
REST_PERIOD_S = 0.05  # model output interval inside the rests
HEADLINE_LAMBDA = 1.0
# name -> (lambda, which rest points enter the residual: the whole trace or the first sample)
WEIGHTINGS = {
    "unweighted": (None, "trace"),
    "lambda=0.25": (0.25, "trace"),
    "lambda=1": (1.0, "trace"),
    "lambda=4": (4.0, "trace"),
    "jump only, lambda=1": (1.0, "jump"),
}
HEADLINE = "lambda=1"
JUMP_VARIANT = "jump only, lambda=1"
PREFIX = {"trace": "v2_", "jump": "v2_jump_"}
FORBIDDEN_RATE = "2C"


def rest_key(rate: str, after: str) -> str:
    return f"{rate} rest after {after}"


class V2Residual:
    """Curve residuals (v1) plus weighted rest residuals, sharing one simulation per rate."""

    def __init__(
        self,
        kind: str = "SPMe",
        disc_mm: float = DEFAULT_DISC_MM,
        lam: float | None = HEADLINE_LAMBDA,
        mode: str = "trace",
        sim_kwargs: dict | None = None,
        param_names=PARAM_NAMES,
    ):
        rates = {r for r, _ in TRAINING} | {r for r, _ in REST_TRAINING}
        assert FORBIDDEN_RATE not in rates, "2C must never enter the fit"
        assert mode in ("trace", "jump")
        cycles = dict(TRAINING + REST_TRAINING)
        self.kind, self.disc_mm, self.lam, self.mode = kind, disc_mm, lam, mode
        self.param_names = list(param_names)
        self.sims = {
            r: RateSimulator(r, cycles[r], kind, disc_mm, rest_period_s=REST_PERIOD_S, mesh="v2", **(sim_kwargs or {})) for r in sorted(rates)
        }
        self.curve_points = {}
        for rate, cyc in TRAINING:
            for k, (t, v) in self.sims[rate].residual_points(N_POINTS_PER_STEP).items():
                self.curve_points[(rate, k)] = (t, v)
        self.rests: dict[tuple[str, str], Rest] = {(r, a): get_rest(r, c, a) for r, c in REST_TRAINING for a in REST_KINDS}
        # Which measured rest samples enter the residual.
        self.rest_idx = {key: (np.arange(len(r.time_s)) if mode == "trace" else np.array([0])) for key, r in self.rests.items()}
        n_curve = sum(len(v[0]) for v in self.curve_points.values())
        n_rest = sum(len(idx) for idx in self.rest_idx.values())
        self.n_curve, self.n_rest = n_curve, n_rest
        self.w_rest = 1.0 if lam is None else float(np.sqrt(lam * n_curve / n_rest))
        self.blocks = [("curve", rate, k) for rate, k in self.curve_points] + [("rest", rate, a) for rate, a in self.rests]

    def block_sizes(self) -> dict[str, int]:
        return {**{f"{r} {k}": len(v[0]) for (r, k), v in self.curve_points.items()}, **{rest_key(r, a): len(idx) for (r, a), idx in self.rest_idx.items()}}

    def _n(self, b) -> int:
        return len(self.curve_points[(b[1], b[2])][0]) if b[0] == "curve" else len(self.rest_idx[(b[1], b[2])])

    def evaluate(self, fit: dict, blocks=None) -> dict[tuple, np.ndarray]:
        """Unweighted residual (model - measured, V) of every block, keyed like ``self.blocks``."""
        blocks = blocks or self.blocks
        rates = {b[1] for b in blocks}
        frames = {r: self.sims[r].run(fit, include_rests=True) for r in rates}
        out = {}
        for b in blocks:
            typ, rate, k = b
            if typ == "curve":
                t, v = self.curve_points[(rate, k)]
                out[b] = model_on_measured_grid(frames[rate][k], t, k) - v
            else:
                r = self.rests[(rate, k)]
                idx = self.rest_idx[(rate, k)]
                f = frames[rate][f"rest after {k}"]
                out[b] = np.interp(r.time_s[idx], f["time_s"].to_numpy(), f["voltage_V"].to_numpy()) - r.voltage_V[idx]
        return out

    def unpack(self, theta_log10) -> dict:
        return {n: 10.0 ** v for n, v in zip(self.param_names, theta_log10)}

    def weighted(self, theta, blocks=None) -> np.ndarray:
        blocks = blocks or self.blocks
        try:
            res = self.evaluate(self.unpack(theta), blocks)
        except (pybamm.SolverError, KeyError):
            # Solver failure, or an experiment that terminated before all its steps ran
            # (the frame for a step is then missing): a flat penalty residual.
            return np.full(sum(self._n(b) for b in blocks), 1.0)
        return np.concatenate([res[b] * (self.w_rest if b[0] == "rest" else 1.0) for b in blocks])

    def rmse_by_type(self, fit: dict) -> dict:
        """Unweighted RMSE of the curves, of the fitted rest points, and of the full rest traces
        (the last is the same number whatever ``mode`` is, so variants can be compared)."""
        res = self.evaluate(fit)
        curve = np.concatenate([v for b, v in res.items() if b[0] == "curve"])
        rest = np.concatenate([v for b, v in res.items() if b[0] == "rest"])
        full = V2Residual.__new__(V2Residual)
        full.__dict__.update(self.__dict__)
        full.rest_idx = {key: np.arange(len(r.time_s)) for key, r in self.rests.items()}
        trace = np.concatenate([v for b, v in full.evaluate(fit).items() if b[0] == "rest"])
        return {
            "curves_mV": float(1000 * np.sqrt(np.mean(curve**2))),
            "rests_mV": float(1000 * np.sqrt(np.mean(rest**2))),
            "rest_traces_mV": float(1000 * np.sqrt(np.mean(trace**2))),
            "per_block_mV": {f"{b[1]} {b[2] if b[0] == 'curve' else 'rest after ' + b[2]}": float(1000 * np.sqrt(np.mean(v**2))) for b, v in res.items()},
        }

    def rest_frames(self, fit: dict) -> pd.DataFrame:
        """Model rest traces (fine grid) and the model on the measured rest times, for figures."""
        rows = []
        for rate, _ in REST_TRAINING:
            frames = self.sims[rate].run(fit, include_rests=True)
            for after in REST_KINDS:
                f = frames[f"rest after {after}"]
                r = self.rests[(rate, after)]
                rows.append(pd.DataFrame({"rate": rate, "after": after, "source": "model", "time_s": f["time_s"], "voltage_V": f["voltage_V"]}))
                rows.append(pd.DataFrame({"rate": rate, "after": after, "source": "measured", "time_s": r.time_s, "voltage_V": r.voltage_V}))
        return pd.concat(rows, ignore_index=True)


def fit_v2(
    kind: str = "SPMe",
    disc_mm: float = DEFAULT_DISC_MM,
    lam: float | None = HEADLINE_LAMBDA,
    x0: dict | None = None,
    blocks=None,
    verbose=True,
    mode: str = "trace",
):
    res = V2Residual(kind, disc_mm, lam, mode)
    x0 = {**DEFAULT_FIT, **(x0 or {})}
    theta0 = np.array([np.log10(x0[n]) for n in PARAM_NAMES])
    lb = np.array([BOUNDS_LOG10[n][0] for n in PARAM_NAMES])
    ub = np.array([BOUNDS_LOG10[n][1] for n in PARAM_NAMES])
    t0 = time.time()
    result = least_squares(
        lambda th: res.weighted(th, blocks), theta0, bounds=(lb, ub), diff_step=1e-3, xtol=1e-5, ftol=1e-6, verbose=2 if verbose else 0
    )
    return result, res, time.time() - t0


def block_label(b) -> str:
    return f"{b[1]} {b[2]}" if b[0] == "curve" else rest_key(b[1], b[2])


def summarise(result, res: V2Residual, elapsed: float, loo: list[dict] | None = None) -> dict:
    fit = unpack(result.x)
    rep = covariance_report(result)
    rm = res.rmse_by_type(fit)
    out = {
        "version": "v2" if res.mode == "trace" else "v2 (jump only)",
        "model": res.kind,
        "disc_mm": res.disc_mm,
        "training_curves": [f"{r} cycle {c} charge+discharge" for r, c in TRAINING],
        "training_rests": [f"{r} cycle {c} rest after {a}" for r, c in REST_TRAINING for a in REST_KINDS],
        "two_c_used": False,
        "particle_radius_m": PARTICLE_RADIUS_M,
        "weighting": {
            "lambda": res.lam,
            "rest_points": "every logged sample of each rest" if res.mode == "trace" else "first logged sample of each rest only (the instantaneous jump)",
            "rest_weight": res.w_rest,
            "rule": "w_rest = sqrt(lambda * N_curve / N_rest): the rest block carries lambda times the summed squared weight of the curve block",
            "n_curve_points": res.n_curve,
            "n_rest_points": res.n_rest,
            "rest_period_s": REST_PERIOD_S,
        },
        "block_sizes": res.block_sizes(),
        "training_rmse_curves_V": rm["curves_mV"] / 1000,
        "training_rmse_rests_V": rm["rests_mV"] / 1000,
        "training_rmse_rest_traces_V": rm["rest_traces_mV"] / 1000,
        "training_rmse_per_block_mV": rm["per_block_mV"],
        "weighted_rmse_V": rep["rmse_V"],
        "n_residuals": int(len(result.fun)),
        "fit_time_s": elapsed,
        "converged": bool(result.success),
        "parameters": {},
        "log10_correlation": {
            f"{a}-{b}": float(rep["corr"][i, j]) for i, a in enumerate(PARAM_NAMES) for j, b in enumerate(PARAM_NAMES) if i < j
        },
    }
    for i, n in enumerate(PARAM_NAMES):
        sd = rep["sd_log10"][i]
        entry = {
            "value": fit[n],
            "log10": float(result.x[i]),
            "sd_log10_jacobian": float(sd),
            "ci95_low_jacobian": float(10 ** (result.x[i] - 1.96 * sd)),
            "ci95_high_jacobian": float(10 ** (result.x[i] + 1.96 * sd)),
        }
        if loo:
            vals = np.array([d[n] for d in loo])
            entry["leave_one_step_out_min"] = float(vals.min())
            entry["leave_one_step_out_max"] = float(vals.max())
        out["parameters"][n] = entry
    tau = PARTICLE_RADIUS_M**2 / fit["Ds"]
    sd = rep["sd_log10"][0]
    out["parameters"]["tau_R2_over_Ds_s"] = {
        "value": tau,
        "ci95_low_jacobian": float(tau / 10 ** (1.96 * sd)),
        "ci95_high_jacobian": float(tau * 10 ** (1.96 * sd)),
    }
    if loo:
        taus = np.array([PARTICLE_RADIUS_M**2 / d["Ds"] for d in loo])
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_min"] = float(taus.min())
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_max"] = float(taus.max())
    out["parameters"]["R_contact"]["unit"] = "Ohm"
    out["parameters"]["Ds"]["unit"] = "m2/s"
    out["parameters"]["k"]["unit"] = "m^2.5/(mol^0.5 s) in j0 = k F sqrt(c_e c_s (c_max - c_s))"
    # Model-side rest jumps at the optimum, for the comparison with the measured ones.
    out["rest_jumps"] = rest_jumps(res, fit)
    return out


def rest_jumps(res: V2Residual, fit: dict) -> dict:
    """Model vs measured instantaneous jump and 30 s relaxation of every training rest."""
    out = {}
    for rate, _ in REST_TRAINING:
        frames = res.sims[rate].run(fit, include_rests=True)
        for after in REST_KINDS:
            r = res.rests[(rate, after)]
            f = frames[f"rest after {after}"]
            v_model = np.interp(r.time_s, f["time_s"].to_numpy(), f["voltage_V"].to_numpy())
            i_mA = abs(r.current_before_mA)
            out[rest_key(rate, after)] = {
                "current_mA": i_mA,
                "first_sample_delay_s": r.first_sample_delay_s,
                "jump_measured_mV": 1000 * abs(r.jump_V),
                "jump_model_mV": float(1000 * abs(v_model[0] - r.v_before_V)),
                "instantaneous_resistance_measured_ohm": r.instantaneous_resistance_ohm,
                "instantaneous_resistance_model_ohm": float(abs(v_model[0] - r.v_before_V) / i_mA * 1000),
                "relaxation_30s_measured_mV": 1000 * abs(r.total_relaxation_V),
                "relaxation_30s_model_mV": float(1000 * abs(v_model[-1] - r.v_before_V)),
                "end_of_rest_error_mV": float(1000 * (v_model[-1] - r.voltage_V[-1])),
            }
    return out


def _loo_refit(args):
    kind, disc_mm, lam, mode, best, keep, label = args
    r, _, _ = fit_v2(kind, disc_mm, lam, x0=best, blocks=keep, verbose=False, mode=mode)
    d = unpack(r.x)
    d["dropped"] = label
    return d


def _fit_with_loo(kind, disc_mm, lam, mode, x0=None):
    """Fit, then leave-one-block-out refits from the optimum. The refits are independent, so
    they run in ``LMO_WORKERS`` processes (default 1); the result does not depend on that."""
    result, res, elapsed = fit_v2(kind, disc_mm, lam, x0=x0, mode=mode)
    best = unpack(result.x)
    jobs = [(kind, disc_mm, lam, mode, best, [b for b in res.blocks if b != drop], block_label(drop)) for drop in res.blocks]
    workers = int(os.environ.get("LMO_WORKERS", "1"))
    if workers > 1:
        with Pool(workers) as pool:
            loo = pool.map(_loo_refit, jobs)
    else:
        loo = [_loo_refit(j) for j in jobs]
    summary = summarise(result, res, elapsed, loo)
    summary["leave_one_step_out"] = loo
    return result, res, summary


def run(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, weighting_sensitivity: bool = True, x0: dict | None = None) -> dict:
    """Headline fit (full rest traces, lambda = 1) with leave-one-block-out refits, the
    weighting sensitivity, and the jump-only variant with its own leave-one-out refits."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    lam, mode = WEIGHTINGS[HEADLINE]
    result, res, summary = _fit_with_loo(kind, disc_mm, lam, mode, x0)
    best = unpack(result.x)
    summary["weighting"]["name"] = HEADLINE
    with open(RESULTS_DIR / f"{PREFIX[mode]}fit_{kind.lower()}.json", "w") as f:
        json.dump(summary, f, indent=2)
    res.rest_frames(best).to_csv(PROCESSED_DIR / f"{PREFIX[mode]}model_rests_{kind.lower()}.csv", index=False)
    if not weighting_sensitivity:
        return summary
    sens = {"model": kind, "headline": HEADLINE, "weightings": {}}
    for name, (lam_i, mode_i) in WEIGHTINGS.items():
        if name == HEADLINE:
            r, rs, el = result, res, summary["fit_time_s"]
        elif name == JUMP_VARIANT:
            r, rs, jump_summary = _fit_with_loo(kind, disc_mm, lam_i, mode_i, x0=best)
            el = jump_summary["fit_time_s"]
            jump_summary["weighting"]["name"] = name
            with open(RESULTS_DIR / f"{PREFIX[mode_i]}fit_{kind.lower()}.json", "w") as f:
                json.dump(jump_summary, f, indent=2)
            rs.rest_frames(unpack(r.x)).to_csv(PROCESSED_DIR / f"{PREFIX[mode_i]}model_rests_{kind.lower()}.csv", index=False)
        else:
            r, rs, el = fit_v2(kind, disc_mm, lam_i, x0=best, verbose=False, mode=mode_i)
        fit = unpack(r.x)
        rep = covariance_report(r)
        rm = rs.rmse_by_type(fit)
        sens["weightings"][name] = {
            "lambda": lam_i,
            "rest_points": "trace" if mode_i == "trace" else "first sample only",
            "rest_weight": rs.w_rest,
            "n_rest_points": rs.n_rest,
            "parameters": fit,
            "tau_R2_over_Ds_s": PARTICLE_RADIUS_M**2 / fit["Ds"],
            "training_rmse_curves_mV": rm["curves_mV"],
            "training_rmse_fitted_rest_points_mV": rm["rests_mV"],
            "training_rmse_rest_traces_mV": rm["rest_traces_mV"],
            "k_R_contact_correlation": float(rep["corr"][1, 2]),
            "rest_jumps": rest_jumps(rs, fit),
            "fit_time_s": el,
        }
    with open(RESULTS_DIR / f"v2_weighting_sensitivity_{kind.lower()}.json", "w") as f:
        json.dump(sens, f, indent=2)
    summary["weighting_sensitivity"] = sens["weightings"]
    return summary


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    out = run(kind, weighting_sensitivity="--no-sens" not in sys.argv)
    print(json.dumps({k: out[k] for k in ("parameters", "log10_correlation", "training_rmse_curves_V", "training_rmse_rests_V")}, indent=2))
