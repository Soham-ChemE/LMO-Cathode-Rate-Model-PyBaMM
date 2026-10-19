"""Blind validation: predict every rate with the fitted parameters, simulate the three CVs,
and run the disc-diameter sensitivity check.

Only C/5 cycle 1 and C/2 cycle 1 were used for fitting. 2C cycle 2 is the blind prediction;
C/10 is a low-rate sanity check; the CVs are a different technique altogether.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pybamm

from . import CVS, PROCESSED_DIR, RATES, REPRESENTATIVE_CYCLE, RESULTS_DIR
from .fit import PARAM_NAMES, RateSimulator, TRAINING, fit_parameters, unpack
from .io import cv_cycle, cv_peaks, get_step, load_cv
from .model import DEFAULT_DISC_MM, build_parameters, default_var_pts, disc_geometry, make_model, ocv_to_sto


def load_fit(kind: str) -> dict:
    with open(RESULTS_DIR / f"fit_{kind.lower()}.json") as f:
        d = json.load(f)
    return {n: d["parameters"][n]["value"] for n in PARAM_NAMES}


def step_metrics(model: pd.DataFrame, meas, kind: str) -> dict:
    """Capacity, mean voltage, energy and V-Q / V-t RMSE of one model step against the data."""
    v_end = 4.5 if kind == "charge" else 3.5
    q_model = model["capacity_mAh"].to_numpy()
    v_model = model["voltage_V"].to_numpy()
    t_model = model["time_s"].to_numpy()
    cap_model = float(q_model[-1])
    e_model = float(np.trapezoid(v_model * np.abs(model["current_A"].to_numpy()), t_model) / 3.6)
    # V-Q RMSE on the measured capacity axis, model held at the cut-off beyond its own capacity.
    v_on_q = np.interp(meas.capacity_mAh, q_model, v_model, right=v_end)
    v_on_t = np.interp(meas.time_s, t_model, v_model, right=v_end)
    common = meas.capacity_mAh <= min(cap_model, meas.total_capacity_mAh)
    return {
        "capacity_mAh_model": cap_model,
        "capacity_mAh_measured": meas.total_capacity_mAh,
        "capacity_error_pct": 100 * (cap_model - meas.total_capacity_mAh) / meas.total_capacity_mAh,
        "mean_voltage_V_model": e_model / cap_model,
        "mean_voltage_V_measured": meas.mean_voltage_V,
        "mean_voltage_error_mV": 1000 * (e_model / cap_model - meas.mean_voltage_V),
        "energy_mWh_model": e_model,
        "energy_mWh_measured": meas.energy_mWh,
        "energy_error_pct": 100 * (e_model - meas.energy_mWh) / meas.energy_mWh,
        "vq_rmse_mV": float(1000 * np.sqrt(np.mean((v_on_q - meas.voltage_V) ** 2))),
        "vq_rmse_common_range_mV": float(1000 * np.sqrt(np.mean((v_on_q[common] - meas.voltage_V[common]) ** 2))),
        "vt_rmse_mV": float(1000 * np.sqrt(np.mean((v_on_t - meas.voltage_V) ** 2))),
        "v_start_model": float(v_model[0]),
        "v_start_measured": float(meas.voltage_V[0]),
    }


def predict_all(kind: str = "SPMe", fit: dict | None = None, disc_mm: float = DEFAULT_DISC_MM, save: bool = True) -> dict:
    fit = fit or load_fit(kind)
    out = {"model": kind, "disc_mm": disc_mm, "parameters": fit, "rates": {}}
    frames_out = []
    for rate in RATES:
        cyc = REPRESENTATIVE_CYCLE[rate]
        sim = RateSimulator(rate, cyc, kind, disc_mm)
        frames = sim.run(fit)
        role = "training" if (rate, cyc) in TRAINING else ("blind prediction" if rate == "2C" else "low-rate check")
        entry = {"cycle": cyc, "role": role}
        for k in ("charge", "discharge"):
            entry[k] = step_metrics(frames[k], get_step(rate, cyc, k), k)
            frames_out.append(frames[k])
        out["rates"][rate] = entry
    two_c = out["rates"]["2C"]["discharge"]
    out["headline"] = {
        "statement": "2C cycle-2 discharge capacity, never used in fitting",
        "predicted_mAh": two_c["capacity_mAh_model"],
        "measured_mAh": two_c["capacity_mAh_measured"],
        "error_pct": two_c["capacity_error_pct"],
        "abs_error_pct": abs(two_c["capacity_error_pct"]),
        "mean_voltage_error_mV": two_c["mean_voltage_error_mV"],
        "energy_error_pct": two_c["energy_error_pct"],
        "vq_rmse_mV": two_c["vq_rmse_mV"],
        "measured_capacity_loss_vs_c10_pct": 100
        * (1 - two_c["capacity_mAh_measured"] / out["rates"]["C/10"]["discharge"]["capacity_mAh_measured"]),
        "predicted_capacity_loss_vs_c10_pct": 100
        * (1 - two_c["capacity_mAh_model"] / out["rates"]["C/10"]["discharge"]["capacity_mAh_model"]),
    }
    if save:
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        pd.concat(frames_out).to_csv(PROCESSED_DIR / f"model_gcd_{kind.lower()}.csv", index=False)
        with open(RESULTS_DIR / f"validation_{kind.lower()}.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


# ----------------------------------------------------------------------------- CV simulation


def simulate_cv(scan: str, fit: dict, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM) -> pd.DataFrame:
    """Drive the model with the measured voltage programme (both cycles) in voltage mode."""
    raw = load_cv(scan)
    t = raw["time_s"].to_numpy()
    v = raw["voltage_V"].to_numpy()
    # Replace the first seconds (instrument inrush at 3.0 V) with the clean ramp.
    params = build_parameters(disc_mm=disc_mm, initial_sto=ocv_to_sto(float(v[0])), v_min=2.0, v_max=5.5, **fit)
    params.update({"Voltage function [V]": lambda tt: pybamm.Interpolant(t, v, tt, name="measured CV programme")}, check_already_exists=False)
    model = make_model(kind, **{"operating mode": "voltage"})
    sim = pybamm.Simulation(model, parameter_values=params, solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8), var_pts=default_var_pts())
    sol = sim.solve(t_eval=[0, t[-1]], t_interp=t)
    df = pd.DataFrame(
        {
            "time_s": sol["Time [s]"].entries,
            "voltage_V": sol["Voltage [V]"].entries,
            "current_mA": -1000 * sol["Current [A]"].entries,  # EC-Lab sign: positive = oxidation (charge)
            "cycle": raw["cycle"].to_numpy()[: len(sol.t)],
        }
    )
    df["scan"] = scan
    return df


def cv_peaks_from_frame(df: pd.DataFrame, cycle: int = 2) -> dict:
    sel = df[df["cycle"] == cycle]
    v = sel["voltage_V"].to_numpy()
    i = sel["current_mA"].to_numpy()
    i_turn = int(np.argmax(v))
    i_ox = int(np.argmax(i[: i_turn + 1]))
    i_red = i_turn + int(np.argmin(i[i_turn:]))
    return {
        "oxidation_peak_V": float(v[i_ox]),
        "oxidation_peak_mA": float(i[i_ox]),
        "reduction_peak_V": float(v[i_red]),
        "reduction_peak_mA": float(i[i_red]),
    }


def b_value(scan_rates, peak_currents) -> dict:
    """Slope of log10|i_p| vs log10(scan rate): 0.5 diffusion-limited, 1.0 capacitive/surface."""
    x = np.log10(np.asarray(scan_rates))
    y = np.log10(np.abs(np.asarray(peak_currents)))
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return {"b": float(slope), "log10_intercept": float(intercept), "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else 1.0}


def randles_sevcik_D(peak_current_A, scan_rate_V_s, area_m2, conc_mol_m3, n: int = 1, T: float = 298.15) -> float:
    """Apparent diffusivity (m2/s) from i_p = 0.4463 n F A C sqrt(n F v D / (R T)) (reversible, planar)."""
    F, R = 96485.33212, 8.314462618
    pref = 0.4463 * n * F * area_m2 * conc_mol_m3 * np.sqrt(n * F * scan_rate_V_s / (R * T))
    return float((abs(peak_current_A) / pref) ** 2)


def validate_cv(kind: str = "SPMe", fit: dict | None = None, disc_mm: float = DEFAULT_DISC_MM, save: bool = True) -> dict:
    fit = fit or load_fit(kind)
    geo = disc_geometry(disc_mm)
    from .ocv import LMO_C_MAX_MOL_M3, sto_window

    sto_min, sto_max = sto_window(get_step("C/10", 2, "discharge").total_capacity_mAh)
    c_active = LMO_C_MAX_MOL_M3 * (sto_max - sto_min)  # lithium that actually cycles
    out = {"model": kind, "disc_mm": disc_mm, "scans": {}}
    frames = []
    rates_m, iox_m, ired_m, iox_s, ired_s = [], [], [], [], []
    for scan, meta in CVS.items():
        meas = cv_peaks(scan, cycle=2)
        df = simulate_cv(scan, fit, kind, disc_mm)
        frames.append(df)
        simp = cv_peaks_from_frame(df, cycle=2)
        out["scans"][scan] = {
            "scan_rate_V_s": meta["scan_rate_V_s"],
            "measured": {k: meas[k] for k in simp},
            "model": simp,
            "oxidation_peak_shift_mV": 1000 * (simp["oxidation_peak_V"] - meas["oxidation_peak_V"]),
            "reduction_peak_shift_mV": 1000 * (simp["reduction_peak_V"] - meas["reduction_peak_V"]),
            "oxidation_peak_current_error_pct": 100 * (simp["oxidation_peak_mA"] / meas["oxidation_peak_mA"] - 1),
            "reduction_peak_current_error_pct": 100 * (simp["reduction_peak_mA"] / meas["reduction_peak_mA"] - 1),
            "peak_separation_measured_mV": 1000 * (meas["oxidation_peak_V"] - meas["reduction_peak_V"]),
            "peak_separation_model_mV": 1000 * (simp["oxidation_peak_V"] - simp["reduction_peak_V"]),
        }
        rates_m.append(meta["scan_rate_V_s"])
        iox_m.append(meas["oxidation_peak_mA"])
        ired_m.append(meas["reduction_peak_mA"])
        iox_s.append(simp["oxidation_peak_mA"])
        ired_s.append(simp["reduction_peak_mA"])
    out["b_value"] = {
        "measured_oxidation": b_value(rates_m, iox_m),
        "measured_reduction": b_value(rates_m, ired_m),
        "model_oxidation": b_value(rates_m, iox_s),
        "model_reduction": b_value(rates_m, ired_s),
    }
    # Apparent Randles-Sevcik diffusivity: geometric disc area, concentration of cycling lithium.
    d_rs = {
        scan: {
            "oxidation_m2_s": randles_sevcik_D(iox_m[j] / 1000, rates_m[j], geo["area_m2"], c_active),
            "reduction_m2_s": randles_sevcik_D(ired_m[j] / 1000, rates_m[j], geo["area_m2"], c_active),
        }
        for j, scan in enumerate(CVS)
    }
    out["randles_sevcik"] = {
        "assumptions": "reversible planar electrode, geometric disc area, C = cycling lithium concentration",
        "area_m2": geo["area_m2"],
        "concentration_mol_m3": c_active,
        "per_scan": d_rs,
        "apparent_D_geometric_mean_m2_s": float(np.exp(np.mean([np.log(v) for d in d_rs.values() for v in d.values()]))),
        "fitted_Ds_m2_s": fit["Ds"],
    }
    out["randles_sevcik"]["ratio_apparent_to_fitted"] = out["randles_sevcik"]["apparent_D_geometric_mean_m2_s"] / fit["Ds"]
    if save:
        pd.concat(frames).to_csv(PROCESSED_DIR / f"model_cv_{kind.lower()}.csv", index=False)
        with open(RESULTS_DIR / f"cv_{kind.lower()}.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


# ------------------------------------------------------------------ disc diameter sensitivity


def disc_sensitivity(kind: str = "SPMe", discs=(12.0, 14.0, 16.0), save: bool = True) -> dict:
    """Refit the three parameters for each assumed disc diameter and re-predict 2C."""
    out = {"model": kind, "discs": {}}
    for d in discs:
        result, _, elapsed = fit_parameters(kind, d, verbose=False)
        fit = unpack(result.x)
        pred = predict_all(kind, fit, d, save=False)
        out["discs"][str(d)] = {
            "geometry": disc_geometry(d),
            "parameters": fit,
            "training_rmse_V": float(np.sqrt(np.mean(result.fun**2))),
            "predicted_2c_capacity_mAh": pred["headline"]["predicted_mAh"],
            "predicted_2c_error_pct": pred["headline"]["error_pct"],
            "predicted_2c_vq_rmse_mV": pred["headline"]["vq_rmse_mV"],
            "fit_time_s": elapsed,
        }
    errs = [v["predicted_2c_error_pct"] for v in out["discs"].values()]
    out["spread_of_2c_error_pct"] = float(max(errs) - min(errs))
    if save:
        with open(RESULTS_DIR / f"sensitivity_disc_{kind.lower()}.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    print(json.dumps(predict_all(kind)["headline"], indent=2))
    print(json.dumps(validate_cv(kind), indent=2))
