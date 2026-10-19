"""All figures. Every curve comes from data/raw via the parser or from data/processed."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import CVS, PROCESSED_DIR, RATES, REPRESENTATIVE_CYCLE, RESULTS_DIR
from .diagnostics import COMPONENTS
from .io import cv_cycle, get_step
from .plotting import COLORS, CV_COLORS, MEASURED_STYLE, MODEL_STYLE, savefig

COMPONENT_LABELS = {
    "contact_ohmic": "Contact / ohmic",
    "reaction": "Charge transfer",
    "particle_diffusion": "Particle diffusion",
    "electrolyte_concentration": "Electrolyte concentration",
    "electrolyte_ohmic": "Electrolyte ohmic",
    "solid_ohmic": "Solid-phase ohmic",
}
COMPONENT_COLORS = ["#7f7f7f", "#ff7f0e", "#17becf", "#bcbd22", "#8c564b", "#e377c2"]


def _model_frames(kind, prefix=""):
    return pd.read_csv(PROCESSED_DIR / f"{prefix}model_gcd_{kind.lower()}.csv")


def fig_headline(kind: str = "DFN", prefix: str = ""):
    df = _model_frames(kind, prefix)
    val = json.load(open(RESULTS_DIR / f"{prefix}validation_{kind.lower()}.json"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, step in zip(axes, ("discharge", "charge")):
        for rate in RATES:
            cyc = REPRESENTATIVE_CYCLE[rate]
            m = get_step(rate, cyc, step)
            s = df[(df["rate"] == rate) & (df["step"] == step)]
            role = val["rates"][rate]["role"]
            tag = " (prediction)" if role == "blind prediction" else (" (fit)" if role == "training" else "")
            ax.plot(m.capacity_mAh, m.voltage_V, color=COLORS[rate], label=f"{rate} measured", **MEASURED_STYLE)
            ax.plot(s["capacity_mAh"], s["voltage_V"], color=COLORS[rate], label=f"{rate} {kind}{tag}", **MODEL_STYLE)
        ax.set_xlabel("Capacity (mAh)")
        ax.set_title(f"{step.capitalize()}")
        ax.set_xlim(0, 1.6)
    axes[0].set_ylabel("Voltage vs Li/Li$^+$ (V)")
    axes[0].set_ylim(3.4, 4.6)
    axes[0].legend(fontsize=7.5, ncol=2, loc="lower left")
    h = val["headline"]
    trained = "C/5 and C/2 curves plus the C/10, C/5, C/2 rests (v2)" if prefix else "C/5 and C/2 only"
    fig.suptitle(
        f"LiMn$_2$O$_4$ half cell: measured vs {kind} fitted on {trained}. "
        f"2C prediction: {h['predicted_mAh']:.3f} mAh vs {h['measured_mAh']:.3f} measured ({h['error_pct']:+.1f}%)",
        fontsize=10,
    )
    return savefig(fig, f"{prefix}headline_measured_vs_model")


def fig_fit_time(kind: str):
    df = _model_frames(kind)
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5))
    for row, rate in enumerate(("C/5", "C/2")):
        for col, step in enumerate(("charge", "discharge")):
            ax = axes[row, col]
            m = get_step(rate, 1, step)
            s = df[(df["rate"] == rate) & (df["step"] == step)]
            ax.plot(m.time_s / 3600, m.voltage_V, color=COLORS[rate], label="measured", **MEASURED_STYLE)
            ax.plot(s["time_s"] / 3600, s["voltage_V"], color="k", label=f"{kind} fit", **MODEL_STYLE)
            ax.set_title(f"{rate} cycle 1 {step} (training)")
            ax.set_xlabel("Time (h)")
            ax.set_ylabel("Voltage (V)")
            ax.legend()
    fig.tight_layout()
    return savefig(fig, f"fit_training_{kind.lower()}")


def fig_ocv():
    ocv = pd.read_csv(PROCESSED_DIR / "pseudo_ocv.csv")
    meta = json.load(open(RESULTS_DIR / "ocv.json"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.plot(ocv["capacity_mAh"], ocv["charge_V"], color=COLORS["C/10"], alpha=0.5, lw=1, label="C/10 charge (cycle 2)")
    ax.plot(ocv["capacity_mAh"], ocv["discharge_V"], color=COLORS["C/10"], alpha=0.5, lw=1, ls=":", label="C/10 discharge (cycle 2)")
    ax.plot(ocv["capacity_mAh"], ocv["ocv_V"], color="k", lw=2, label="pseudo-OCV (average, smoothed)")
    ax.set_xlabel("Capacity removed (mAh)")
    ax.set_ylabel("Potential vs Li/Li$^+$ (V)")
    ax.set_title("Measured pseudo-OCV")
    sec = ax.secondary_xaxis("top", functions=(lambda q: 1 - q / meta["theoretical_capacity_mAh"], lambda s: (1 - s) * meta["theoretical_capacity_mAh"]))
    sec.set_xlabel("x in Li$_x$Mn$_2$O$_4$")
    ax.legend(fontsize=8)
    ax = axes[1]
    for kind, ls in (("discharge", "-"), ("charge", "--")):
        ic = pd.read_csv(PROCESSED_DIR / f"dqdv_c10_{kind}.csv")
        ax.plot(ic["voltage_V"], ic["dqdv_mAh_V"], color=COLORS["C/10"], ls=ls, label=f"C/10 {kind}")
        for p in meta[f"dqdv_peaks_{kind}"]:
            ax.annotate(f"{p['voltage_V']:.2f} V", (p["voltage_V"], p["dqdv_mAh_V"]), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=8)
    ax.set_xlim(3.8, 4.4)
    ax.set_xlabel("Voltage vs Li/Li$^+$ (V)")
    ax.set_ylabel("dQ/dV (mAh/V)")
    ax.set_title("Incremental capacity at C/10: two spinel plateaus")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return savefig(fig, "pseudo_ocv_and_dqdv")


def fig_overpotentials(kind: str, prefix: str = ""):
    ov = json.load(open(RESULTS_DIR / f"{prefix}overpotentials_{kind.lower()}.json"))
    bd = pd.read_csv(PROCESSED_DIR / f"{prefix}overpotentials_{kind.lower()}.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    dis = [s for s in ov["steps"] if s["step"] == "discharge"]
    x = np.arange(len(dis))
    bottom = np.zeros(len(dis))
    for comp, color in zip(COMPONENTS, COMPONENT_COLORS):
        vals = np.array([s["capacity_weighted_mean_mV"][comp] for s in dis])
        ax.bar(x, vals, bottom=bottom, color=color, label=COMPONENT_LABELS[comp], width=0.6)
        bottom += vals
    ax.set_xticks(x, [s["rate"] for s in dis])
    for xi, rate in zip(x, [s["rate"] for s in dis]):
        ax.get_xticklabels()[xi].set_color(COLORS[rate])
    ax.set_ylabel("Capacity-weighted mean voltage loss on discharge (mV)")
    ax.set_title(f"Overpotential breakdown ({kind}{', v2' if prefix else ''})")
    ax.legend(fontsize=8)
    ax = axes[1]
    sel = bd[(bd["rate"] == "2C") & (bd["step"] == "discharge")]
    ax.plot(sel["capacity_mAh"], sel["ocv_V"], color="k", lw=1.2, label="bulk OCV (state of charge)")
    cum = sel["ocv_V"].to_numpy().copy()
    for comp, color in zip(COMPONENTS, COMPONENT_COLORS):
        new = cum - sel[comp].to_numpy()
        ax.fill_between(sel["capacity_mAh"], new, cum, color=color, alpha=0.8, label=COMPONENT_LABELS[comp])
        cum = new
    ax.plot(sel["capacity_mAh"], sel["voltage_V"], color=COLORS["2C"], lw=2, label="terminal voltage (model)")
    m = get_step("2C", 2, "discharge")
    ax.plot(m.capacity_mAh, m.voltage_V, color=COLORS["2C"], lw=1.2, ls=":", label="terminal voltage (measured)")
    ax.axhline(3.5, color="grey", lw=0.8, ls="--")
    ax.set_xlabel("Capacity (mAh)")
    ax.set_ylabel("Voltage vs Li/Li$^+$ (V)")
    ax.set_title("2C discharge: where the voltage goes")
    ax.set_ylim(3.3, 4.3)
    ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    return savefig(fig, f"{prefix}overpotentials_{kind.lower()}")


def fig_cv(kind: str, prefix: str = ""):
    sim = pd.read_csv(PROCESSED_DIR / f"{prefix}model_cv_{kind.lower()}.csv")
    cv = json.load(open(RESULTS_DIR / f"{prefix}cv_{kind.lower()}.json"))
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    for ax, scan in zip(axes[:3], CVS):
        m = cv_cycle(scan, 2)
        s = sim[(sim["scan"] == scan) & (sim["cycle"] == 2)]
        ax.plot(m["voltage_V"], m["current_mA"], color=CV_COLORS[scan], label="measured (cycle 2)", **MEASURED_STYLE)
        ax.plot(s["voltage_V"], s["current_mA"], color="k", label=f"{kind}", **MODEL_STYLE)
        ax.set_title(scan)
        ax.set_xlabel("Voltage vs Li/Li$^+$ (V)")
        ax.set_ylabel("Current (mA)")
        ax.legend(fontsize=8)
    ax = axes[3]
    rates = np.array([CVS[s]["scan_rate_V_s"] for s in CVS]) * 1000
    for key, marker, label in (("measured", "o", "measured"), ("model", "s", kind)):
        iox = [abs(cv["scans"][s][key]["oxidation_peak_mA"]) for s in CVS]
        ired = [abs(cv["scans"][s][key]["reduction_peak_mA"]) for s in CVS]
        ax.loglog(rates, iox, marker=marker, ls="none", color="#d62728", label=f"oxidation, {label}")
        ax.loglog(rates, ired, marker=marker, ls="none", color="#1f77b4", label=f"reduction, {label}", mfc="none" if key == "model" else None)
    b = cv["b_value"]
    for name, color, ls in (("measured_oxidation", "#d62728", "-"), ("measured_reduction", "#1f77b4", "-"), ("model_oxidation", "#d62728", "--"), ("model_reduction", "#1f77b4", "--")):
        xx = np.linspace(rates.min() * 0.9, rates.max() * 1.1, 10)
        ax.loglog(xx, 10 ** b[name]["log10_intercept"] * (xx / 1000) ** b[name]["b"], color=color, ls=ls, lw=1, label=f"b = {b[name]['b']:.2f}")
    ax.set_xlabel("Scan rate (mV/s)")
    ax.set_ylabel("|Peak current| (mA)")
    ax.set_title("b-value: log i$_p$ vs log v")
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout()
    return savefig(fig, f"{prefix}cv_{kind.lower()}")


def fig_sensitivity(kind: str = "SPMe"):
    sens = json.load(open(RESULTS_DIR / f"sensitivity_disc_{kind.lower()}.json"))
    val = json.load(open(RESULTS_DIR / f"validation_{kind.lower()}.json"))
    discs = sorted(sens["discs"], key=float)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    d = [float(x) for x in discs]
    axes[0].plot(d, [sens["discs"][x]["predicted_2c_capacity_mAh"] for x in discs], "o-", color=COLORS["2C"], label="predicted")
    axes[0].axhline(val["headline"]["measured_mAh"], color="k", ls="--", label="measured")
    axes[0].set_ylabel("2C discharge capacity (mAh)")
    axes[0].legend()
    from .model import PLAUSIBLE_POROSITY

    axes[1].axhspan(*PLAUSIBLE_POROSITY, color="#2ca02c", alpha=0.12, label="plausible porosity")
    axes[1].plot(d, [sens["discs"][x]["geometry"]["porosity"] for x in discs], "o-", color="grey", label="porosity")
    axes[1].plot(d, [sens["discs"][x]["geometry"]["active_volume_fraction"] for x in discs], "s-", color="k", label="active volume fraction")
    axes[1].set_ylabel("Volume fraction")
    axes[1].legend(fontsize=8)
    ref = sens["discs"]["10.0"]["parameters"] if "10.0" in sens["discs"] else sens["discs"][discs[0]]["parameters"]
    for name, marker, color, label in (("Ds", "o", "#17becf", "D$_s$"), ("k", "s", "#ff7f0e", "k"), ("R_contact", "^", "#7f7f7f", "R$_{contact}$")):
        axes[2].plot(d, [sens["discs"][x]["parameters"][name] / ref[name] for x in discs], marker=marker, ls="-", color=color, label=label)
    axes[2].axhline(1, color="k", lw=0.6)
    axes[2].set_ylabel("Refitted parameter / value at 10 mm")
    axes[2].set_ylim(0.9, 1.1)
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.set_xlabel("Assumed disc diameter (mm)")
        ax.set_xticks(d)
    fig.suptitle(f"Sensitivity to the unknown electrode area ({kind}, refitted per disc)", fontsize=10)
    fig.tight_layout()
    return savefig(fig, f"sensitivity_disc_{kind.lower()}")


def fig_raw_data():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for rate in RATES:
        for cyc, ls in ((1, "-"), (2, ":")):
            m = get_step(rate, cyc, "discharge")
            axes[0].plot(m.capacity_mAh, m.voltage_V, color=COLORS[rate], ls=ls, lw=1.4, label=f"{rate} cycle {cyc}")
    axes[0].set_xlabel("Discharge capacity (mAh)")
    axes[0].set_ylabel("Voltage vs Li/Li$^+$ (V)")
    axes[0].set_title("Measured discharges, both cycles")
    axes[0].legend(fontsize=7, ncol=2)
    for scan in CVS:
        m = cv_cycle(scan, 2)
        axes[1].plot(m["voltage_V"], m["current_mA"], color=CV_COLORS[scan], lw=1.4, label=f"{scan} cycle 2")
    axes[1].set_xlabel("Voltage vs Li/Li$^+$ (V)")
    axes[1].set_ylabel("Current (mA)")
    axes[1].set_title("Measured CVs")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    return savefig(fig, "raw_data_overview")


def make_all(headline_kind: str = "DFN", kinds=("SPMe", "DFN")):
    paths = [fig_raw_data(), fig_ocv(), fig_headline(headline_kind)]
    for kind in kinds:
        if (RESULTS_DIR / f"validation_{kind.lower()}.json").exists():
            paths += [fig_fit_time(kind), fig_overpotentials(kind), fig_cv(kind)]
        if (RESULTS_DIR / f"sensitivity_disc_{kind.lower()}.json").exists():
            paths.append(fig_sensitivity(kind))
    return paths
