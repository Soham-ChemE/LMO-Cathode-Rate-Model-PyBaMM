"""Overpotential breakdown at each C-rate from the model's voltage-component variables.

For a discharge step PyBaMM gives V = OCV_bulk + eta_particle + eta_reaction + eta_conc
+ eta_ohmic_electrolyte + eta_ohmic_solid - eta_contact (all but the last are negative on
discharge). The closure of that sum is checked and reported. The capacity-weighted mean of each
term over the discharge and its value at the cut-off are tabulated, and the lithium still left in
the particles when the cut-off is reached is converted into the capacity the model says the rate
has "locked out".
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import PROCESSED_DIR, RATES, RESULTS_DIR
from .model import ocv_to_sto
from .ocv import sto_window, theoretical_capacity_mAh
from .io import get_step

COMPONENTS = {
    "contact_ohmic": ("eta_contact_V", -1.0),
    "reaction": ("eta_reaction_V", 1.0),
    "particle_diffusion": ("eta_particle_V", 1.0),
    "electrolyte_concentration": ("eta_electrolyte_conc_V", 1.0),
    "electrolyte_ohmic": ("eta_electrolyte_ohmic_V", 1.0),
    "solid_ohmic": ("eta_solid_ohmic_V", 1.0),
}


def breakdown_frame(step: pd.DataFrame) -> pd.DataFrame:
    """Signed voltage losses (V, positive = loss) of one step on its own capacity axis."""
    out = pd.DataFrame({"time_s": step["time_s"], "capacity_mAh": step["capacity_mAh"], "voltage_V": step["voltage_V"], "ocv_V": step["ocv_V"]})
    sign = 1.0 if step["step"].iloc[0] == "discharge" else -1.0
    for name, (col, s) in COMPONENTS.items():
        out[name] = -sign * s * step[col]
    out["closure_V"] = step["ocv_V"] - sign * out[list(COMPONENTS)].sum(axis=1) - step["voltage_V"]
    return out


def summarise_step(step: pd.DataFrame, rate: str) -> dict:
    bd = breakdown_frame(step)
    q = bd["capacity_mAh"].to_numpy()
    w = np.gradient(q)
    w = w / w.sum()
    total = bd[list(COMPONENTS)].sum(axis=1).to_numpy()
    res = {
        "rate": rate,
        "step": step["step"].iloc[0],
        "mean_current_mA": float(np.abs(step["current_A"]).mean() * 1000),
        "capacity_mAh": float(q[-1]),
        "closure_max_abs_V": float(bd["closure_V"].abs().max()),
        "capacity_weighted_mean_mV": {k: float(1000 * np.sum(w * bd[k].to_numpy())) for k in COMPONENTS},
        "at_cutoff_mV": {k: float(1000 * bd[k].iloc[-1]) for k in COMPONENTS},
        "total_mean_mV": float(1000 * np.sum(w * total)),
        "total_at_cutoff_mV": float(1000 * total[-1]),
        "bulk_ocv_at_cutoff_V": float(bd["ocv_V"].iloc[-1]),
    }
    shares = {k: v / res["total_mean_mV"] for k, v in res["capacity_weighted_mean_mV"].items()}
    res["share_of_mean_loss"] = shares
    if res["step"] == "discharge":
        # Lithium left in the particles at the cut-off, as capacity the model could not deliver.
        sto_end = ocv_to_sto(res["bulk_ocv_at_cutoff_V"])
        sto_min, sto_max = sto_window(get_step("C/10", 2, "discharge").total_capacity_mAh)
        res["model_capacity_locked_out_mAh"] = float((sto_max - sto_end) * theoretical_capacity_mAh())
    return res


def run(kind: str = "SPMe", prefix: str = "") -> dict:
    """``prefix="v2_"`` reads and writes the v2 files."""
    df = pd.read_csv(PROCESSED_DIR / f"{prefix}model_gcd_{kind.lower()}.csv")
    out = {"model": kind, "components": list(COMPONENTS), "steps": []}
    if prefix:
        out["version"] = prefix.rstrip("_")
    frames = []
    for rate in RATES:
        for step_kind in ("discharge", "charge"):
            sel = df[(df["rate"] == rate) & (df["step"] == step_kind)]
            if sel.empty:
                continue
            out["steps"].append(summarise_step(sel, rate))
            bd = breakdown_frame(sel)
            bd["rate"] = rate
            bd["step"] = step_kind
            frames.append(bd)
    pd.concat(frames).to_csv(PROCESSED_DIR / f"{prefix}overpotentials_{kind.lower()}.csv", index=False)
    dis = {s["rate"]: s for s in out["steps"] if s["step"] == "discharge"}
    out["discharge_mean_loss_mV"] = {r: s["total_mean_mV"] for r, s in dis.items()}
    out["discharge_dominant_term"] = {r: max(s["capacity_weighted_mean_mV"], key=s["capacity_weighted_mean_mV"].get) for r, s in dis.items()}
    out["two_c"] = {
        "mean_loss_mV": dis["2C"]["total_mean_mV"],
        "shares": dis["2C"]["share_of_mean_loss"],
        "at_cutoff_mV": dis["2C"]["at_cutoff_mV"],
        "bulk_ocv_at_cutoff_V": dis["2C"]["bulk_ocv_at_cutoff_V"],
        "capacity_locked_out_mAh": dis["2C"]["model_capacity_locked_out_mAh"],
    }
    with open(RESULTS_DIR / f"{prefix}overpotentials_{kind.lower()}.json", "w") as f:
        json.dump(out, f, indent=2)
    return out


if __name__ == "__main__":
    import sys

    print(json.dumps(run(sys.argv[1] if len(sys.argv) > 1 else "SPMe")["two_c"], indent=2))
