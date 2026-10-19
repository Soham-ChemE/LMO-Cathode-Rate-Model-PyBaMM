"""Measured pseudo-OCV of the LMO electrode and dQ/dV analysis, both from the C/10 data.

PyBaMM ships no LiMn2O4 open-circuit potential, so the model uses the one measured here: the
average of the C/10 cycle-2 charge and discharge voltage curves on a common normalised-capacity
axis, smoothed and forced monotonic. The normalised capacity is mapped onto a lithium
stoichiometry window of Li_x Mn2O4 so that PyBaMM can use it as a function of stoichiometry.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, savgol_filter

from . import NOMINAL_CAPACITY_MAH, PROCESSED_DIR, RESULTS_DIR
from .io import get_step

# LiMn2O4: molar mass and density give the maximum lithium concentration in the solid.
LMO_MOLAR_MASS_KG_MOL = 0.18081  # Li1 Mn2 O4
LMO_DENSITY_KG_M3 = 4280.0
LMO_C_MAX_MOL_M3 = LMO_DENSITY_KG_M3 / LMO_MOLAR_MASS_KG_MOL  # 23671 mol/m3
LMO_THEORETICAL_MAH_G = 96485.0 / LMO_MOLAR_MASS_KG_MOL / 3.6 / 1000  # 148.2 mAh/g
RATED_MAH_G = 110.0  # MTI rating used to size the electrode mass from the EC-Lab capacity

# The 3.5 V end of the window is taken as Li_1 Mn2O4 (sto = 1). Discharging below ~3.3 V would
# enter the tetragonal Li2Mn2O4 plateau, which this cell never visited.
STO_MAX = 1.0


def active_mass_g() -> float:
    """Active LMO mass implied by the EC-Lab capacity and the 110 mAh/g rating."""
    return NOMINAL_CAPACITY_MAH / RATED_MAH_G


def theoretical_capacity_mAh() -> float:
    return active_mass_g() * LMO_THEORETICAL_MAH_G


def sto_window(usable_capacity_mAh: float) -> tuple[float, float]:
    """(sto_min, sto_max) covered by ``usable_capacity_mAh`` of the electrode."""
    return STO_MAX - usable_capacity_mAh / theoretical_capacity_mAh(), STO_MAX


def build_pseudo_ocv(cycle: int = 2, n_grid: int = 401, smooth_window: int = 21) -> pd.DataFrame:
    """Average C/10 charge and discharge on a normalised-capacity axis.

    ``x`` runs from 0 (discharged, 3.5 V) to 1 (charged, 4.5 V) and each step is normalised by
    its own capacity, which lines up both voltage cut-offs and the plateau transition near
    4.05-4.08 V (it sits at x = 0.47-0.48 in both steps). The charge step is ~6% longer than the
    discharge because of parasitic oxidation at high voltage; normalising by its own capacity
    spreads that excess over the whole curve. The discharge capacity defines the stoichiometry
    window and the capacity axis of the output.
    """
    chg = get_step("C/10", cycle, "charge")
    dis = get_step("C/10", cycle, "discharge")
    q_ref = dis.total_capacity_mAh
    x_chg = chg.capacity_mAh / chg.total_capacity_mAh
    x_dis = 1.0 - dis.capacity_mAh / q_ref
    x = np.linspace(0.0, 1.0, n_grid)
    v_chg = np.interp(x, x_chg, chg.voltage_V)
    v_dis = np.interp(x, x_dis[::-1], dis.voltage_V[::-1])
    v_avg = 0.5 * (v_chg + v_dis)
    v_smooth = savgol_filter(v_avg, smooth_window, 3)
    # Monotonic (voltage must increase with x): cumulative maximum, then a tiny slope floor.
    v_mono = np.maximum.accumulate(v_smooth)
    v_mono = v_mono + 1e-4 * x  # guarantees strictly increasing for the interpolant
    sto_min, sto_max = sto_window(q_ref)
    sto = sto_max - x * (sto_max - sto_min)
    return pd.DataFrame(
        {
            "x_norm": x,
            "stoichiometry": sto,
            "ocv_V": v_mono,
            "charge_V": v_chg,
            "discharge_V": v_dis,
            "capacity_mAh": x * q_ref,
        }
    )


def ocv_table_for_pybamm(ocv: pd.DataFrame, pad: float = 0.08) -> np.ndarray:
    """(sto, V) table sorted by increasing sto, extended linearly past both ends.

    The extension lets the solver reach the 3.5 and 4.5 V cut-offs with a well-defined slope
    instead of failing at the edge of the measured window.
    """
    sto = ocv["stoichiometry"].to_numpy()[::-1]
    v = ocv["ocv_V"].to_numpy()[::-1]
    # End slopes from the outer 2% of the measured window.
    n = max(3, len(sto) // 50)
    s_lo = np.polyfit(sto[:n], v[:n], 1)[0]
    s_hi = np.polyfit(sto[-n:], v[-n:], 1)[0]
    s_lo = min(s_lo, -1.0)  # at least 1 V per unit sto outside the window
    s_hi = min(s_hi, -1.0)
    lo = np.linspace(sto[0] - pad, sto[0], 9, endpoint=False)
    hi = np.linspace(sto[-1], sto[-1] + pad, 10)[1:]
    v_lo = v[0] + s_lo * (lo - sto[0])
    v_hi = v[-1] + s_hi * (hi - sto[-1])
    return np.column_stack([np.r_[lo, sto, hi], np.r_[v_lo, v, v_hi]])


def dqdv(step, bin_mV: float = 5.0, smooth_bins: int = 5) -> pd.DataFrame:
    """Incremental capacity dQ/dV (mAh/V) of one step on a uniform voltage grid.

    Capacity is binned in ``bin_mV`` voltage bins (robust to the 10 mV / 5 s recording
    condition) and lightly smoothed with a moving average.
    """
    v = step.voltage_V
    q = step.capacity_mAh
    order = np.argsort(v)
    v_s, q_s = v[order], q[order]
    edges = np.arange(np.floor(v.min() * 200) / 200, v.max() + bin_mV / 1000, bin_mV / 1000)
    # Monotonic Q(V) by taking the cumulative capacity reached at each voltage edge.
    q_edges = np.interp(edges, v_s, q_s)
    if step.kind == "discharge":
        q_edges = q_edges[0] - q_edges  # discharge: Q grows as V falls
    dq = np.abs(np.diff(q_edges))
    dv = np.diff(edges)
    centres = 0.5 * (edges[1:] + edges[:-1])
    y = uniform_filter1d(dq / dv, smooth_bins, mode="nearest")
    return pd.DataFrame({"voltage_V": centres, "dqdv_mAh_V": y})


def dqdv_peaks(df: pd.DataFrame, v_min: float = 3.8, v_max: float = 4.3, n: int = 2) -> list[dict]:
    sel = df[(df["voltage_V"] > v_min) & (df["voltage_V"] < v_max)]
    idx, props = find_peaks(sel["dqdv_mAh_V"].to_numpy(), prominence=0.5)
    peaks = sorted(zip(props["prominences"], idx), reverse=True)[:n]
    out = [
        {"voltage_V": float(sel["voltage_V"].iloc[i]), "dqdv_mAh_V": float(sel["dqdv_mAh_V"].iloc[i])}
        for _, i in peaks
    ]
    return sorted(out, key=lambda p: p["voltage_V"])


def run(cycle: int = 2) -> dict:
    """Build the pseudo-OCV and dQ/dV, write data/processed and results/ocv.json."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ocv = build_pseudo_ocv(cycle)
    ocv.to_csv(PROCESSED_DIR / "pseudo_ocv.csv", index=False)
    table = ocv_table_for_pybamm(ocv)
    np.savetxt(PROCESSED_DIR / "lmo_ocp_vs_sto.csv", table, delimiter=",", header="stoichiometry,ocp_V", comments="")

    dis = get_step("C/10", cycle, "discharge")
    chg = get_step("C/10", cycle, "charge")
    ic_dis = dqdv(dis)
    ic_chg = dqdv(chg)
    ic_dis.to_csv(PROCESSED_DIR / "dqdv_c10_discharge.csv", index=False)
    ic_chg.to_csv(PROCESSED_DIR / "dqdv_c10_charge.csv", index=False)
    sto_min, sto_max = sto_window(dis.total_capacity_mAh)
    hyst = ocv["charge_V"] - ocv["discharge_V"]
    summary = {
        "source": f"C/10 cycle {cycle}, average of charge and discharge vs normalised capacity",
        "reference_discharge_capacity_mAh": dis.total_capacity_mAh,
        "charge_capacity_mAh": chg.total_capacity_mAh,
        "active_mass_g": active_mass_g(),
        "lmo_c_max_mol_m3": LMO_C_MAX_MOL_M3,
        "theoretical_capacity_mAh": theoretical_capacity_mAh(),
        "stoichiometry_window": {"sto_min": sto_min, "sto_max": sto_max},
        "ocv_at_x0_V": float(ocv["ocv_V"].iloc[0]),
        "ocv_at_x1_V": float(ocv["ocv_V"].iloc[-1]),
        "mean_charge_discharge_hysteresis_mV": float(hyst.mean() * 1000),
        "median_charge_discharge_hysteresis_mV": float(hyst.median() * 1000),
        "dqdv_peaks_discharge": dqdv_peaks(ic_dis),
        "dqdv_peaks_charge": dqdv_peaks(ic_chg),
    }
    with open(RESULTS_DIR / "ocv.json", "w") as f:
        json.dump(summary, f, indent=2)
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
