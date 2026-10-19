"""Parsers for the EC-Lab text exports in data/raw.

GCD files: ``time/s, Ewe/V, <I>/mA, Capacity/mA.h, cycle number`` (EC-Lab cycles 0 and 1,
renumbered here to 1 and 2). CV files: ``time/s, Ewe/V, <I>/mA, cycle number`` (cycles 1 and 2).
Files are tab separated, latin-1, CRLF.

Capacity is always obtained by integrating the measured current with the trapezoid rule, never
from the EC-Lab capacity column, because the C/10 run was recorded on a current range whose
reading differs from the set point.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import CVS, RATES, RAW_DIR, REPRESENTATIVE_CYCLE

GCD_COLUMNS = ["time_s", "voltage_V", "current_mA", "eclab_capacity_mAh", "cycle"]
CV_COLUMNS = ["time_s", "voltage_V", "current_mA", "cycle"]


def _read_eclab(path, columns):
    df = pd.read_csv(path, sep="\t", header=0, encoding="latin-1", usecols=range(len(columns)))
    df.columns = columns
    df["cycle"] = df["cycle"].round().astype(int)
    return df


def load_gcd(rate: str) -> pd.DataFrame:
    """Load one galvanostatic file and label its steps.

    Returns a DataFrame with the raw columns plus ``cycle`` (1-based), ``step`` in
    {"rest", "charge", "discharge"} and ``step_id`` (a running index of contiguous steps).
    """
    df = _read_eclab(RAW_DIR / f"{RATES[rate]['file']}.txt", GCD_COLUMNS)
    df["cycle"] = df["cycle"] + 1
    threshold = 0.02 * np.abs(df["current_mA"]).max()
    state = np.where(df["current_mA"] > threshold, 1, np.where(df["current_mA"] < -threshold, -1, 0))
    df["step"] = np.array(["rest", "charge", "discharge"])[state]  # -1 indexes "discharge"
    df["step_id"] = np.r_[0, np.cumsum(np.diff(state) != 0)]
    return df


def load_cv(scan: str) -> pd.DataFrame:
    """Load one CV file. ``cycle`` is as recorded by EC-Lab (1 and 2)."""
    return _read_eclab(RAW_DIR / f"{CVS[scan]['file']}.txt", CV_COLUMNS)


def integrate_charge_mAh(time_s, current_mA):
    """Trapezoid integral of |I| dt in mAh."""
    return float(np.trapezoid(np.abs(current_mA), time_s) / 3600.0)


@dataclass
class Step:
    rate: str
    cycle: int
    kind: str  # "charge" or "discharge"
    time_s: np.ndarray
    voltage_V: np.ndarray
    current_mA: np.ndarray
    capacity_mAh: np.ndarray = field(repr=False)  # cumulative, from integrated current

    @property
    def total_capacity_mAh(self) -> float:
        return float(self.capacity_mAh[-1])

    @property
    def energy_mWh(self) -> float:
        return float(np.trapezoid(np.abs(self.current_mA) * self.voltage_V, self.time_s) / 3600.0)

    @property
    def mean_voltage_V(self) -> float:
        return self.energy_mWh / self.total_capacity_mAh

    @property
    def mean_current_mA(self) -> float:
        return float(np.trapezoid(np.abs(self.current_mA), self.time_s) / (self.time_s[-1] - self.time_s[0]))

    @property
    def duration_s(self) -> float:
        return float(self.time_s[-1] - self.time_s[0])


def get_step(rate: str, cycle: int, kind: str, df: pd.DataFrame | None = None) -> Step:
    """Extract one charge or discharge step. Time is reset to zero at the step start."""
    if df is None:
        df = load_gcd(rate)
    sel = df[(df["cycle"] == cycle) & (df["step"] == kind)]
    if sel.empty:
        raise KeyError(f"no {kind} step for {rate} cycle {cycle}")
    # A cycle should contain exactly one contiguous step of each kind.
    ids = sel["step_id"].unique()
    if len(ids) != 1:
        raise ValueError(f"{rate} cycle {cycle} has {len(ids)} {kind} segments")
    t = sel["time_s"].to_numpy()
    i = sel["current_mA"].to_numpy()
    v = sel["voltage_V"].to_numpy()
    q = np.concatenate([[0.0], np.cumsum(0.5 * (np.abs(i[1:]) + np.abs(i[:-1])) * np.diff(t)) / 3600.0])
    return Step(rate, cycle, kind, t - t[0], v, i, q)


def gcd_summary() -> pd.DataFrame:
    """Per rate / cycle / step table of integrated capacity, energy and mean voltage."""
    rows = []
    for rate in RATES:
        df = load_gcd(rate)
        for cycle in sorted(df["cycle"].unique()):
            for kind in ("charge", "discharge"):
                s = get_step(rate, cycle, kind, df)
                rows.append(
                    {
                        "rate": rate,
                        "cycle": int(cycle),
                        "step": kind,
                        "representative": cycle == REPRESENTATIVE_CYCLE[rate],
                        "capacity_mAh": s.total_capacity_mAh,
                        "energy_mWh": s.energy_mWh,
                        "mean_voltage_V": s.mean_voltage_V,
                        "mean_current_mA": s.mean_current_mA,
                        "duration_h": s.duration_s / 3600.0,
                        "v_start_V": float(s.voltage_V[0]),
                        "v_end_V": float(s.voltage_V[-1]),
                    }
                )
    out = pd.DataFrame(rows)
    # Coulombic efficiency per cycle (discharge / charge).
    ce = (
        out.pivot_table(index=["rate", "cycle"], columns="step", values="capacity_mAh")
        .eval("discharge / charge")
        .rename("coulombic_efficiency")
    )
    return out.merge(ce.reset_index(), on=["rate", "cycle"])


def cv_cycle(scan: str, cycle: int = 2, drop_inrush_s: float = 60.0) -> pd.DataFrame:
    """One CV cycle with time reset to zero.

    ``drop_inrush_s`` removes the first seconds of the 0.2 mV/s file, which start with an
    instrument inrush current at 3.0 V (see data/README.md). It only affects cycle 1.
    """
    df = load_cv(scan)
    sel = df[df["cycle"] == cycle].copy()
    if cycle == 1 and drop_inrush_s > 0:
        sel = sel[sel["time_s"] > sel["time_s"].iloc[0] + drop_inrush_s]
    sel["time_s"] -= sel["time_s"].iloc[0]
    return sel.reset_index(drop=True)


def cv_peaks(scan: str, cycle: int = 2) -> dict:
    """Oxidation and reduction peak position and current of one CV cycle.

    The anodic peak is the current maximum on the positive-going sweep, the cathodic peak the
    minimum on the negative-going sweep. Auto-range steps (data/README.md) are plateaus, not
    maxima, so they do not win the argmax.
    """
    sel = cv_cycle(scan, cycle)
    v = sel["voltage_V"].to_numpy()
    i = sel["current_mA"].to_numpy()
    i_turn = int(np.argmax(v))
    i_ox = int(np.argmax(i[: i_turn + 1]))
    i_red = i_turn + int(np.argmin(i[i_turn:]))
    return {
        "scan": scan,
        "cycle": cycle,
        "scan_rate_V_s": CVS[scan]["scan_rate_V_s"],
        "oxidation_peak_V": float(v[i_ox]),
        "oxidation_peak_mA": float(i[i_ox]),
        "reduction_peak_V": float(v[i_red]),
        "reduction_peak_mA": float(i[i_red]),
    }
