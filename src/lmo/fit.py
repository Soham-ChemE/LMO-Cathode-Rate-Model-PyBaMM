"""Fit the three effective parameters to C/5 cycle 1 and C/2 cycle 1 (charge and discharge).

Parameters are fitted in log10 space with scipy ``least_squares``:
  * Ds   positive particle diffusivity at the fixed 5 um radius, reported as tau = R^2/Ds
  * k    reaction rate constant of the Butler-Volmer exchange current density
  * R_c  lumped contact resistance

The residual is the model minus measured voltage on each training step at ~300 points per step,
with the model held at the cut-off voltage once it has terminated (so a model that ends a step
too early or too late is penalised through the voltage curve, not through a separate capacity
term). Uncertainties come from the Jacobian covariance at the optimum and, as a second, more
conservative view, from refits that leave one of the four training steps out.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

from . import RATES, RESULTS_DIR
from .io import Step, get_step
from .model import (
    DEFAULT_DISC_MM,
    DEFAULT_FIT,
    PARTICLE_RADIUS_M,
    build_parameters,
    gcd_experiment,
    make_model,
    make_simulation,
    ocv_to_sto,
    solution_frames,
)

TRAINING = [("C/5", 1), ("C/2", 1)]
PARAM_NAMES = ["Ds", "k", "R_contact"]
BOUNDS_LOG10 = {"Ds": (-17.0, -11.0), "k": (-14.0, -8.0), "R_contact": (-1.0, 3.0)}
N_POINTS_PER_STEP = 300


class RateSimulator:
    """One PyBaMM simulation of the full measured protocol at one rate, parameters as inputs."""

    def __init__(self, rate: str, cycle: int, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, n_cycles: int | None = None):
        self.rate, self.cycle, self.kind = rate, cycle, kind
        self.steps = {k: get_step(rate, cycle, k) for k in ("charge", "discharge")}
        # The run is simulated from its first cycle so that the state entering ``cycle`` carries
        # the history of the earlier cycle (matters for 2C, whose cycle 1 started over-lithiated).
        self.n_cycles = n_cycles or cycle
        first = {k: get_step(rate, 1, k) for k in ("charge", "discharge")}
        currents = (first["charge"].mean_current_mA, first["discharge"].mean_current_mA)
        self.start_voltage = float(first["charge"].voltage_V[0])
        params = build_parameters(
            Ds=pybamm.InputParameter("Ds"),
            k=pybamm.InputParameter("k"),
            R_contact=pybamm.InputParameter("R_contact"),
            disc_mm=disc_mm,
            initial_sto=ocv_to_sto(self.start_voltage),
        )
        self.sim = make_simulation(gcd_experiment(currents, n_cycles=self.n_cycles), params, kind, make_model(kind))

    def run(self, fit: dict) -> dict[str, pd.DataFrame]:
        sol = self.sim.solve(inputs={n: float(fit[n]) for n in PARAM_NAMES})
        frames = solution_frames(sol, self.n_cycles, self.rate)
        return {f["step"].iloc[0]: f for f in frames if f["cycle"].iloc[0] == self.cycle}

    def residual_points(self, n: int = N_POINTS_PER_STEP) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Measured (t, V) subsampled to n evenly spaced points on each step."""
        out = {}
        for k, s in self.steps.items():
            idx = np.unique(np.linspace(0, len(s.time_s) - 1, n).astype(int))
            out[k] = (s.time_s[idx], s.voltage_V[idx])
        return out


def model_on_measured_grid(frame: pd.DataFrame, t_meas: np.ndarray, kind: str, v_min=3.5, v_max=4.5) -> np.ndarray:
    """Interpolate the model step onto measured times; past the model's end hold the cut-off."""
    v_end = v_max if kind == "charge" else v_min
    return np.interp(t_meas, frame["time_s"].to_numpy(), frame["voltage_V"].to_numpy(), right=v_end)


def unpack(theta_log10: np.ndarray) -> dict:
    return {n: 10.0 ** v for n, v in zip(PARAM_NAMES, theta_log10)}


def make_residual(sims: list[RateSimulator], steps_to_use=None):
    points = {(s.rate, k): v for s in sims for k, v in s.residual_points().items()}

    def residual(theta):
        fit = unpack(theta)
        res = []
        for s in sims:
            try:
                frames = s.run(fit)
            except pybamm.SolverError:
                return np.full(sum(len(v[0]) for v in points.values()), 1.0)
            for k in ("charge", "discharge"):
                if steps_to_use is not None and (s.rate, k) not in steps_to_use:
                    continue
                t, v = points[(s.rate, k)]
                res.append(model_on_measured_grid(frames[k], t, k) - v)
        return np.concatenate(res)

    return residual


def fit_parameters(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, x0: dict | None = None, steps_to_use=None, verbose=True):
    sims = [RateSimulator(rate, cyc, kind, disc_mm) for rate, cyc in TRAINING]
    residual = make_residual(sims, steps_to_use)
    x0 = {**DEFAULT_FIT, **(x0 or {})}
    theta0 = np.array([np.log10(x0[n]) for n in PARAM_NAMES])
    lb = np.array([BOUNDS_LOG10[n][0] for n in PARAM_NAMES])
    ub = np.array([BOUNDS_LOG10[n][1] for n in PARAM_NAMES])
    t0 = time.time()
    result = least_squares(residual, theta0, bounds=(lb, ub), diff_step=1e-3, xtol=1e-5, ftol=1e-6, verbose=2 if verbose else 0)
    elapsed = time.time() - t0
    return result, sims, elapsed


def covariance_report(result) -> dict:
    """Jacobian covariance of the log10 parameters; 95% intervals as multiplicative factors."""
    J = result.jac
    n, p = J.shape
    s2 = np.sum(result.fun**2) / max(n - p, 1)
    try:
        cov = np.linalg.inv(J.T @ J) * s2
    except np.linalg.LinAlgError:
        cov = np.full((p, p), np.nan)
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    return {"sd_log10": sd, "corr": corr, "rmse_V": float(np.sqrt(np.mean(result.fun**2)))}


def summarise(result, kind: str, disc_mm: float, elapsed: float, loo: list[dict] | None = None) -> dict:
    fit = unpack(result.x)
    rep = covariance_report(result)
    out = {
        "model": kind,
        "disc_mm": disc_mm,
        "training": [f"{r} cycle {c} charge+discharge" for r, c in TRAINING],
        "particle_radius_m": PARTICLE_RADIUS_M,
        "training_rmse_V": rep["rmse_V"],
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
    return out


def run(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, leave_one_out: bool = True, x0: dict | None = None) -> dict:
    result, sims, elapsed = fit_parameters(kind, disc_mm, x0=x0)
    best = unpack(result.x)
    loo = []
    if leave_one_out:
        all_steps = [(r, k) for r, _ in TRAINING for k in ("charge", "discharge")]
        for drop in all_steps:
            keep = [s for s in all_steps if s != drop]
            r, _, _ = fit_parameters(kind, disc_mm, x0=best, steps_to_use=keep, verbose=False)
            d = unpack(r.x)
            d["dropped"] = f"{drop[0]} {drop[1]}"
            loo.append(d)
    summary = summarise(result, kind, disc_mm, elapsed, loo)
    summary["leave_one_step_out"] = loo
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_DIR / f"fit_{kind.lower()}.json", "w") as f:
        json.dump(summary, f, indent=2)
    return summary


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    print(json.dumps(run(kind), indent=2))
