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
block. lambda = 1 (equal block weight) is the headline; the pipeline also refits at the other
values in ``WEIGHTINGS`` as a sensitivity check.
"""

from __future__ import annotations

import json
import time

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
WEIGHTINGS = {"unweighted": None, "lambda=0.25": 0.25, "lambda=1": 1.0, "lambda=4": 4.0}
FORBIDDEN_RATE = "2C"


def rest_key(rate: str, after: str) -> str:
    return f"{rate} rest after {after}"


class V2Residual:
    """Curve residuals (v1) plus weighted rest residuals, sharing one simulation per rate."""

    def __init__(self, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, lam: float | None = HEADLINE_LAMBDA):
        rates = {r for r, _ in TRAINING} | {r for r, _ in REST_TRAINING}
        assert FORBIDDEN_RATE not in rates, "2C must never enter the fit"
        cycles = dict(TRAINING + REST_TRAINING)
        self.kind, self.disc_mm, self.lam = kind, disc_mm, lam
        self.sims = {r: RateSimulator(r, cycles[r], kind, disc_mm, rest_period_s=REST_PERIOD_S, mesh="v2") for r in sorted(rates)}
        self.curve_points = {}
        for rate, cyc in TRAINING:
            for k, (t, v) in self.sims[rate].residual_points(N_POINTS_PER_STEP).items():
                self.curve_points[(rate, k)] = (t, v)
        self.rests: dict[tuple[str, str], Rest] = {(r, a): get_rest(r, c, a) for r, c in REST_TRAINING for a in REST_KINDS}
        n_curve = sum(len(v[0]) for v in self.curve_points.values())
        n_rest = sum(len(r.time_s) for r in self.rests.values())
        self.n_curve, self.n_rest = n_curve, n_rest
        self.w_rest = 1.0 if lam is None else float(np.sqrt(lam * n_curve / n_rest))
        self.blocks = [("curve", rate, k) for rate, k in self.curve_points] + [("rest", rate, a) for rate, a in self.rests]

    def block_sizes(self) -> dict[str, int]:
        return {**{f"{r} {k}": len(v[0]) for (r, k), v in self.curve_points.items()}, **{rest_key(r, a): len(x.time_s) for (r, a), x in self.rests.items()}}

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
                f = frames[rate][f"rest after {k}"]
                out[b] = np.interp(r.time_s, f["time_s"].to_numpy(), f["voltage_V"].to_numpy()) - r.voltage_V
        return out

    def weighted(self, theta, blocks=None) -> np.ndarray:
        blocks = blocks or self.blocks
        try:
            res = self.evaluate(unpack(theta), blocks)
        except pybamm.SolverError:
            n = sum(len(self.curve_points[(b[1], b[2])][0]) if b[0] == "curve" else len(self.rests[(b[1], b[2])].time_s) for b in blocks)
            return np.full(n, 1.0)
        return np.concatenate([res[b] * (self.w_rest if b[0] == "rest" else 1.0) for b in blocks])

    def rmse_by_type(self, fit: dict) -> dict:
        res = self.evaluate(fit)
        curve = np.concatenate([v for b, v in res.items() if b[0] == "curve"])
        rest = np.concatenate([v for b, v in res.items() if b[0] == "rest"])
        return {
            "curves_mV": float(1000 * np.sqrt(np.mean(curve**2))),
            "rests_mV": float(1000 * np.sqrt(np.mean(rest**2))),
            "per_block_mV": {f"{b[1]} {b[2] if b[0] == 'curve' else 'rest after ' + b[2]}": float(1000 * np.sqrt(np.mean(v**2))) for b, v in res.items()},
        }

    def rest_frames(self, fit: dict) -> pd.DataFrame:
        """Model rest traces (fine grid) and the model on the measured rest times, for figures."""
        rows = []
        for rate in sorted({r for r, _ in REST_TRAINING}):
            frames = self.sims[rate].run(fit, include_rests=True)
            for after in REST_KINDS:
                f = frames[f"rest after {after}"]
                r = self.rests[(rate, after)]
                rows.append(pd.DataFrame({"rate": rate, "after": after, "source": "model", "time_s": f["time_s"], "voltage_V": f["voltage_V"]}))
                rows.append(pd.DataFrame({"rate": rate, "after": after, "source": "measured", "time_s": r.time_s, "voltage_V": r.voltage_V}))
        return pd.concat(rows, ignore_index=True)


def fit_v2(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, lam: float | None = HEADLINE_LAMBDA, x0: dict | None = None, blocks=None, verbose=True):
    res = V2Residual(kind, disc_mm, lam)
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
        "version": "v2",
        "model": res.kind,
        "disc_mm": res.disc_mm,
        "training_curves": [f"{r} cycle {c} charge+discharge" for r, c in TRAINING],
        "training_rests": [f"{r} cycle {c} rest after {a}" for r, c in REST_TRAINING for a in REST_KINDS],
        "two_c_used": False,
        "particle_radius_m": PARTICLE_RADIUS_M,
        "weighting": {
            "lambda": res.lam,
            "rest_weight": res.w_rest,
            "rule": "w_rest = sqrt(lambda * N_curve / N_rest): the rest block carries lambda times the summed squared weight of the curve block",
            "n_curve_points": res.n_curve,
            "n_rest_points": res.n_rest,
            "rest_period_s": REST_PERIOD_S,
        },
        "block_sizes": res.block_sizes(),
        "training_rmse_curves_V": rm["curves_mV"] / 1000,
        "training_rmse_rests_V": rm["rests_mV"] / 1000,
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
    for rate in sorted({r for r, _ in REST_TRAINING}):
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


def run(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, leave_one_out: bool = True, weighting_sensitivity: bool = True, x0: dict | None = None) -> dict:
    result, res, elapsed = fit_v2(kind, disc_mm, HEADLINE_LAMBDA, x0=x0)
    best = unpack(result.x)
    loo = []
    if leave_one_out:
        for drop in res.blocks:
            keep = [b for b in res.blocks if b != drop]
            r, _, _ = fit_v2(kind, disc_mm, HEADLINE_LAMBDA, x0=best, blocks=keep, verbose=False)
            d = unpack(r.x)
            d["dropped"] = block_label(drop)
            loo.append(d)
    summary = summarise(result, res, elapsed, loo)
    summary["leave_one_step_out"] = loo
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / f"v2_fit_{kind.lower()}.json", "w") as f:
        json.dump(summary, f, indent=2)
    res.rest_frames(best).to_csv(PROCESSED_DIR / f"v2_model_rests_{kind.lower()}.csv", index=False)
    if weighting_sensitivity:
        sens = {"model": kind, "headline_lambda": HEADLINE_LAMBDA, "weightings": {}}
        for name, lam in WEIGHTINGS.items():
            if lam == HEADLINE_LAMBDA:
                r, rs, el = result, res, elapsed
            else:
                r, rs, el = fit_v2(kind, disc_mm, lam, x0=best, verbose=False)
            fit = unpack(r.x)
            rep = covariance_report(r)
            rm = rs.rmse_by_type(fit)
            sens["weightings"][name] = {
                "lambda": lam,
                "rest_weight": rs.w_rest,
                "parameters": fit,
                "tau_R2_over_Ds_s": PARTICLE_RADIUS_M**2 / fit["Ds"],
                "training_rmse_curves_mV": rm["curves_mV"],
                "training_rmse_rests_mV": rm["rests_mV"],
                "k_R_contact_correlation": float(rep["corr"][1, 2]),
                "fit_time_s": el,
            }
        with open(RESULTS_DIR / f"v2_weighting_sensitivity_{kind.lower()}.json", "w") as f:
            json.dump(sens, f, indent=2)
        summary["weighting_sensitivity"] = sens["weightings"]
    return summary


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    out = run(kind, leave_one_out="--no-loo" not in sys.argv, weighting_sensitivity="--no-sens" not in sys.argv)
    print(json.dumps({k: out[k] for k in ("parameters", "log10_correlation", "training_rmse_curves_V", "training_rmse_rests_V")}, indent=2))
