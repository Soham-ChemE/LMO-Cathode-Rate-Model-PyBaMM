"""v2: the post-step rests as data (model-free check), the post-hoc 2C rest comparison, the
v1 vs v2 comparison, and the v2 figures.

Nothing here feeds back into the fit. The 2C rests are compared with the v2 model only after
the fit, and the comparison file simply reads the v1 and v2 result files side by side.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import FIGURES_DIR, PROCESSED_DIR, RATES, RESULTS_DIR
from .fit import RateSimulator
from .fit_v2 import REST_KINDS, REST_PERIOD_S, REST_TRAINING, WEIGHTINGS, rest_key
from .io import all_rests, get_rest, get_step, rest_summary
from .plotting import COLORS, MEASURED_STYLE, MODEL_STYLE, savefig
from .validate import load_fit, load_fit_all

TRAINED_RATES = [r for r, _ in REST_TRAINING]


# ----------------------------------------------------------------- 2. model-free check


def measured_rests(save: bool = True) -> dict:
    """Instantaneous resistance and 30 s relaxation of every rest, straight from the raw files.

    The first rest sample is logged ~48 ms after the current is switched off, so its jump holds
    the ohmic drop plus whatever charge-transfer overpotential has relaxed by then: the
    resistance is an upper bound on the pure ohmic resistance. The 2C rows are listed for
    completeness but flagged; they are never used in fitting.
    """
    tab = rest_summary()
    tab["used_in_v2_fit"] = [(r in TRAINED_RATES) and (c == 1) for r, c in zip(tab["rate"], tab["cycle"])]
    tab["post_hoc_only"] = tab["rate"] == "2C"
    by_rate = {}
    for rate in RATES:
        sel = tab[tab["rate"] == rate]
        chg = sel[sel["after"] == "charge"]["instantaneous_resistance_ohm"]
        dis = sel[sel["after"] == "discharge"]["instantaneous_resistance_ohm"]
        by_rate[rate] = {
            "after_charge_ohm": [float(x) for x in chg],
            "after_discharge_ohm": [float(x) for x in dis],
            "after_charge_mean_ohm": float(chg.mean()),
            "after_discharge_mean_ohm": float(dis.mean()),
            "cycle1_after_charge_ohm": float(sel[(sel["after"] == "charge") & (sel["cycle"] == 1)]["instantaneous_resistance_ohm"].iloc[0]),
            "used_in_v2_fit": rate in TRAINED_RATES,
        }
    fitted = tab[tab["rate"] != "2C"]
    out = {
        "note": (
            "jump = first rest sample minus last sample under current; resistance = |jump| / |I|. The first sample "
            "is logged tens of milliseconds after switch-off and already contains the fastest part of the charge-"
            "transfer relaxation, so this is an upper bound on the pure ohmic resistance."
        ),
        "first_sample_delay_s": {"min": float(tab["first_sample_delay_s"].min()), "max": float(tab["first_sample_delay_s"].max())},
        "rests": tab.to_dict(orient="records"),
        "by_rate": by_rate,
        "after_charge_mean_ohm_c10_c5_c2": float(fitted[fitted["after"] == "charge"]["instantaneous_resistance_ohm"].mean()),
        "after_discharge_mean_ohm_c10_c5_c2": float(fitted[fitted["after"] == "discharge"]["instantaneous_resistance_ohm"].mean()),
        "after_charge_mean_ohm_2c_post_hoc": by_rate["2C"]["after_charge_mean_ohm"],
        "after_discharge_mean_ohm_2c_post_hoc": by_rate["2C"]["after_discharge_mean_ohm"],
        "v1_R_contact_ohm": _v1_R_contact(),
    }
    if save:
        tab.to_csv(PROCESSED_DIR / "v2_rests_measured.csv", index=False)
        with open(RESULTS_DIR / "v2_rests_measured.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


def _v1_R_contact(kind: str = "DFN"):
    p = RESULTS_DIR / f"fit_{kind.lower()}.json"
    if not p.exists():
        p = RESULTS_DIR / "fit_spme.json"
    return json.load(open(p))["parameters"]["R_contact"]["value"] if p.exists() else None


# ------------------------------------------------------- 5. post-hoc 2C rests (not fitted)


def posthoc_2c_rests(kind: str = "SPMe", fit: dict | None = None, save: bool = True, prefix: str = "v2_", sim_kwargs: dict | None = None) -> dict:
    """Compare the measured 2C rest jumps with the fitted model's, after the fit. Not used for fitting.
    ``prefix`` "v2_" (default) or "v2_1_" with the variant's ``sim_kwargs``."""
    fit = fit or (load_fit_all(kind, prefix) if sim_kwargs else load_fit(kind, prefix))
    sim = RateSimulator("2C", 2, kind, rest_period_s=REST_PERIOD_S, mesh="v2", **(sim_kwargs or {}))
    frames = sim.run(fit, include_rests=True)
    version = prefix.rstrip("_").replace("v2_1", "v2.1")
    out = {"model": kind, "version": version, "post_hoc_only": True, "note": "2C rests were never used in fitting; comparison made after the fit", "rests": {}}
    traces = []
    for after in REST_KINDS:
        r = get_rest("2C", 2, after)
        f = frames[f"rest after {after}"]
        v_model = np.interp(r.time_s, f["time_s"].to_numpy(), f["voltage_V"].to_numpy())
        i_mA = abs(r.current_before_mA)
        out["rests"][rest_key("2C", after)] = {
            "cycle": 2,
            "current_mA": i_mA,
            "jump_measured_mV": 1000 * abs(r.jump_V),
            "jump_model_mV": float(1000 * abs(v_model[0] - r.v_before_V)),
            "instantaneous_resistance_measured_ohm": r.instantaneous_resistance_ohm,
            "instantaneous_resistance_model_ohm": float(abs(v_model[0] - r.v_before_V) / i_mA * 1000),
            "relaxation_30s_measured_mV": 1000 * abs(r.total_relaxation_V),
            "relaxation_30s_model_mV": float(1000 * abs(v_model[-1] - r.v_before_V)),
            "rest_rmse_mV": float(1000 * np.sqrt(np.mean((v_model - r.voltage_V) ** 2))),
        }
        traces.append(pd.DataFrame({"rate": "2C", "after": after, "source": "model", "time_s": f["time_s"], "voltage_V": f["voltage_V"]}))
        traces.append(pd.DataFrame({"rate": "2C", "after": after, "source": "measured", "time_s": r.time_s, "voltage_V": r.voltage_V}))
    out["fitted_R_contact_ohm"] = fit["R_contact"]
    out["v2_R_contact_ohm"] = fit["R_contact"]  # kept for older readers of the v2 file
    out["measured_2c_after_charge_ohm"] = out["rests"][rest_key("2C", "charge")]["instantaneous_resistance_measured_ohm"]
    out["measured_2c_after_discharge_ohm"] = out["rests"][rest_key("2C", "discharge")]["instantaneous_resistance_measured_ohm"]
    if save:
        pd.concat(traces).to_csv(PROCESSED_DIR / f"{prefix}model_rests_2c_{kind.lower()}.csv", index=False)
        with open(RESULTS_DIR / f"{prefix}rest_posthoc_2c_{kind.lower()}.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


# ------------------------------------------------------------------- 7. v1 vs v2 comparison


def _load(name):
    p = RESULTS_DIR / name
    return json.load(open(p)) if p.exists() else None


def comparison(kinds=("SPMe", "DFN"), save: bool = True) -> dict:
    out = {"models": {}}
    for kind in kinds:
        k = kind.lower()
        versions = {}
        for ver, pre in (("v1", ""), ("v2", "v2_"), ("v2.1", "v2_1_")):
            fj, vj = _load(f"{pre}fit_{k}.json"), _load(f"{pre}validation_{k}.json")
            if fj and vj:
                versions[ver] = (fj, vj, _load(f"{pre}overpotentials_{k}.json"), _load(f"{pre}cv_{k}.json"))
        if "v1" not in versions or "v2" not in versions:
            continue
        v1v = versions["v1"][1]
        entry = {}
        for ver, (fj, vj, oj, cj) in versions.items():
            h = vj["headline"]
            c10 = vj["rates"]["C/10"]["discharge"]
            e = {
                "parameters": {n: fj["parameters"][n]["value"] for n in fj["parameters"] if n != "tau_R2_over_Ds_s"},
                "model_options": fj.get("model_options", {}),
                "variant": fj.get("variant"),
                "tau_R2_over_Ds_s": fj["parameters"]["tau_R2_over_Ds_s"]["value"],
                "k_R_contact_correlation": fj["log10_correlation"]["k-R_contact"],
                "Ds_k_correlation": fj["log10_correlation"]["Ds-k"],
                "Ds_R_contact_correlation": fj["log10_correlation"]["Ds-R_contact"],
                "training_rmse_curves_mV": 1000 * (fj.get("training_rmse_curves_V") or fj["training_rmse_V"]),
                "training_rmse_rests_mV": 1000 * fj["training_rmse_rests_V"] if "training_rmse_rests_V" in fj else None,
                "two_c": {
                    "capacity_mAh": h["predicted_mAh"],
                    "capacity_error_pct": h["error_pct"],
                    "mean_voltage_error_mV": h["mean_voltage_error_mV"],
                    "energy_error_pct": h["energy_error_pct"],
                    "vq_rmse_mV": h["vq_rmse_mV"],
                },
                "c10": {"capacity_error_pct": c10["capacity_error_pct"], "vq_rmse_mV": c10["vq_rmse_mV"]},
                "rates": {r: {"capacity_error_pct": x["discharge"]["capacity_error_pct"], "vq_rmse_mV": x["discharge"]["vq_rmse_mV"]} for r, x in vj["rates"].items()},
            }
            if oj:
                e["two_c_shares"] = oj["two_c"]["shares"]
                e["two_c_mean_loss_mV"] = oj["two_c"]["mean_loss_mV"]
            if cj:
                e["cv_max_peak_shift_mV"] = max(max(abs(x["oxidation_peak_shift_mV"]), abs(x["reduction_peak_shift_mV"])) for x in cj["scans"].values())
                e["b_value_model"] = {"oxidation": cj["b_value"]["model_oxidation"]["b"], "reduction": cj["b_value"]["model_reduction"]["b"]}
            entry[ver] = e
        entry["measured_2c_capacity_mAh"] = v1v["headline"]["measured_mAh"]
        entry["two_c_abs_error_change_pct_points"] = abs(entry["v2"]["two_c"]["capacity_error_pct"]) - abs(entry["v1"]["two_c"]["capacity_error_pct"])
        entry["two_c_improved"] = bool(entry["two_c_abs_error_change_pct_points"] < 0)
        if "v2.1" in entry:
            e21, e1 = entry["v2.1"]["two_c"], entry["v1"]["two_c"]
            entry["v2_1_two_c_abs_error_change_pct_points"] = abs(e21["capacity_error_pct"]) - abs(e1["capacity_error_pct"])
            entry["v2_1_two_c_improved"] = bool(abs(e21["capacity_error_pct"]) < abs(e1["capacity_error_pct"]) and e21["vq_rmse_mV"] < e1["vq_rmse_mV"])
        out["models"][kind] = entry
    out["rests_measured"] = _load("v2_rests_measured.json") and {
        k: v for k, v in _load("v2_rests_measured.json").items() if k.startswith("after_") or k == "v1_R_contact_ohm"
    }
    if save:
        with open(RESULTS_DIR / "v2_comparison.json", "w") as f:
            json.dump(out, f, indent=2)
    return out


# ------------------------------------------------------------------------------- figures


def _logt(t):
    t = np.asarray(t, dtype=float).copy()
    t[t <= 0] = 1e-3
    return t


def fig_rests_measured():
    rests = all_rests()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, after in zip(axes[:2], REST_KINDS):
        for r in rests:
            if r.after != after or r.cycle != 1:
                continue
            dv = 1000 * (r.voltage_V - r.v_before_V)
            ls = ":" if r.rate == "2C" else "-"
            ax.plot(_logt(r.time_s), dv, marker="o", ms=3, ls=ls, color=COLORS[r.rate], label=f"{r.rate}{' (not fitted)' if r.rate == '2C' else ''}")
        ax.set_xscale("log")
        ax.set_xlabel("Time after switch-off (s)")
        ax.set_ylabel("V(rest) - V(last sample under current) (mV)")
        ax.set_title(f"Rest after {after} (cycle 1)")
        ax.legend(fontsize=8)
    ax = axes[2]
    tab = rest_summary()
    x = np.arange(len(RATES))
    for j, (after, marker) in enumerate(zip(REST_KINDS, ("o", "s"))):
        for cyc, mfc in ((1, None), (2, "none")):
            sel = tab[(tab["after"] == after) & (tab["cycle"] == cyc)].set_index("rate").reindex(list(RATES))
            ax.plot(x + (j - 0.5) * 0.15, sel["instantaneous_resistance_ohm"], marker=marker, ls="none", color="k", mfc=mfc, label=f"after {after}, cycle {cyc}")
    v1 = _v1_R_contact()
    if v1:
        ax.axhline(v1, color="#d62728", ls="--", label=f"v1 fitted R_contact ({v1:.0f} Ω)")
    v2 = _load("v2_fit_dfn.json") or _load("v2_fit_spme.json")
    if v2:
        ax.axhline(v2["parameters"]["R_contact"]["value"], color="#2ca02c", ls="--", label=f"v2 fitted R_contact ({v2['parameters']['R_contact']['value']:.0f} Ω)")
    ax.set_xticks(x, list(RATES))
    ax.set_ylabel("Instantaneous jump / current (Ω)")
    ax.set_title("Model-free resistance per rest (upper bound on ohmic)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return savefig(fig, "v2_rests_measured")


def fig_fit_rests(kind: str, prefix: str = "v2_"):
    df = pd.read_csv(PROCESSED_DIR / f"{prefix}model_rests_{kind.lower()}.csv")
    fitj = json.load(open(RESULTS_DIR / f"{prefix}fit_{kind.lower()}.json"))
    label = prefix.rstrip("_").replace("v2_1", "v2.1")
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for j, rate in enumerate(TRAINED_RATES):
        for i, after in enumerate(REST_KINDS):
            ax = axes[i, j]
            r = get_rest(rate, 1, after)
            m = df[(df["rate"] == rate) & (df["after"] == after) & (df["source"] == "model")]
            ax.plot(_logt(np.r_[0.0, r.time_s]), np.r_[r.v_before_V, r.voltage_V], marker="o", ms=3, color=COLORS[rate], label="measured", **MEASURED_STYLE)
            ax.plot(_logt(m["time_s"]), m["voltage_V"], color="k", label=f"{kind} {label} fit", **MODEL_STYLE)
            jr = fitj["rest_jumps"][rest_key(rate, after)]
            ax.set_title(f"{rate} rest after {after}: jump {jr['jump_measured_mV']:.0f} mV meas / {jr['jump_model_mV']:.0f} mV model", fontsize=9)
            ax.set_xscale("log")
            ax.set_xlabel("Time after switch-off (s)")
            ax.set_ylabel("Voltage (V)")
            ax.legend(fontsize=8)
    opts = fitj.get("model_options")
    fig.suptitle(
        f"{label} training rests ({kind}{', ' + ', '.join(f'{a} = {b}' for a, b in opts.items()) if opts else ''}): measured vs fitted; "
        "the point at 1 ms is the last sample under current",
        fontsize=10,
    )
    fig.tight_layout()
    return savefig(fig, f"{prefix}fit_rests_{kind.lower()}")


def fig_comparison(kind: str = "DFN"):
    comp = json.load(open(RESULTS_DIR / "v2_comparison.json"))["models"][kind]
    versions = [("v1", "", "#7f7f7f", ":"), ("v2", "v2_", "k", "--")]
    if "v2.1" in comp:
        versions.append(("v2.1", "v2_1_", "#d62728", "-."))
    curves = {ver: pd.read_csv(PROCESSED_DIR / f"{pre}model_gcd_{kind.lower()}.csv") for ver, pre, _, _ in versions}
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    ax = axes[0]
    m = get_step("2C", 2, "discharge")
    ax.plot(m.capacity_mAh, m.voltage_V, color=COLORS["2C"], label="2C measured (cycle 2)", **MEASURED_STYLE)
    for ver, _, color, _ in versions:
        s = curves[ver][(curves[ver]["rate"] == "2C") & (curves[ver]["step"] == "discharge")]
        e = comp[ver]["two_c"]
        ax.plot(s["capacity_mAh"], s["voltage_V"], color=color, label=f"{ver} {kind}: {e['capacity_mAh']:.3f} mAh ({e['capacity_error_pct']:+.1f} %)", **MODEL_STYLE)
    ax.set_xlabel("Discharge capacity (mAh)")
    ax.set_ylabel("Voltage vs Li/Li$^+$ (V)")
    ax.set_title("Blind 2C discharge: " + " vs ".join(v[0] for v in versions))
    ax.set_ylim(3.4, 4.4)
    ax.legend(fontsize=8)
    ax = axes[1]
    for rate in RATES:
        if rate == "2C":
            continue
        mm = get_step(rate, 1, "discharge")
        ax.plot(mm.capacity_mAh, mm.voltage_V, color=COLORS[rate], label=f"{rate} measured", **MEASURED_STYLE)
        for ver, _, _, ls in versions:
            s = curves[ver][(curves[ver]["rate"] == rate) & (curves[ver]["step"] == "discharge")]
            ax.plot(s["capacity_mAh"], s["voltage_V"], color=COLORS[rate], ls=ls, lw=1.3, label=f"{rate} {ver}")
    ax.set_xlabel("Discharge capacity (mAh)")
    ax.set_title("C/10, C/5, C/2 discharges (v1 dotted, v2 dashed" + (", v2.1 dash-dot)" if "v2.1" in comp else ")"))
    ax.set_ylim(3.4, 4.6)
    ax.legend(fontsize=6, ncol=3)
    ax = axes[2]
    last_pre = versions[-1][1]
    posthoc = _load(f"{last_pre}rest_posthoc_2c_{kind.lower()}.json")
    for after, color in zip(REST_KINDS, (COLORS["2C"], "#17becf")):
        r = get_rest("2C", 2, after)
        ax.plot(_logt(np.r_[0.0, r.time_s]), 1000 * (np.r_[r.v_before_V, r.voltage_V] - r.v_before_V), marker="o", ms=3, color=color, label=f"after {after}, measured", **MEASURED_STYLE)
        for ver, pre, _, ls in versions[1:]:
            path = PROCESSED_DIR / f"{pre}model_rests_2c_{kind.lower()}.csv"
            if path.exists():
                tr = pd.read_csv(path)
                s = tr[(tr["after"] == after) & (tr["source"] == "model")]
                ax.plot(_logt(s["time_s"]), 1000 * (s["voltage_V"] - r.v_before_V), color=color, ls=ls, lw=1.6, label=f"after {after}, {ver} {kind}")
    ax.set_xscale("log")
    ax.set_xlabel("Time after switch-off (s)")
    ax.set_ylabel("V(rest) - V(cut-off) (mV)")
    title = "2C rests, post-hoc (never fitted)"
    if posthoc:
        title += f"\nmeasured {posthoc['measured_2c_after_charge_ohm']:.0f} Ω after charge vs {versions[-1][0]} R_contact {posthoc['fitted_R_contact_ohm']:.0f} Ω"
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    return savefig(fig, "v2_comparison_2c")


def fig_weighting(kind: str):
    sens = json.load(open(RESULTS_DIR / f"v2_weighting_sensitivity_{kind.lower()}.json"))
    names = list(sens["weightings"])
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    w = [sens["weightings"][n] for n in names]
    axes[0].plot(x, [e["parameters"]["R_contact"] for e in w], "o-", color="#7f7f7f")
    v1 = _v1_R_contact(kind)
    if v1:
        axes[0].axhline(v1, color="#d62728", ls="--", label=f"v1 ({v1:.0f} Ω)")
    meas = _load("v2_rests_measured.json")
    if meas:
        axes[0].axhline(meas["after_charge_mean_ohm_c10_c5_c2"], color="k", ls=":", label="measured jump / I, after charge")
    axes[0].set_ylabel("R_contact (Ω)")
    axes[0].legend(fontsize=8)
    axes[1].semilogy(x, [e["parameters"]["Ds"] for e in w], "o-", color="#17becf", label="Ds (m²/s)")
    axes[1].semilogy(x, [e["parameters"]["k"] for e in w], "s-", color="#ff7f0e", label="k")
    axes[1].legend(fontsize=8)
    axes[1].set_ylabel("Fitted value")
    axes[2].plot(x, [e["training_rmse_curves_mV"] for e in w], "o-", color="k", label="curves RMSE (mV)")
    axes[2].plot(x, [e["training_rmse_rest_traces_mV"] for e in w], "s-", color="#2ca02c", label="rest traces RMSE (mV)")
    ax2 = axes[2].twinx()
    ax2.plot(x, [e["k_R_contact_correlation"] for e in w], "^--", color="#9467bd", label="k vs R_contact corr.")
    ax2.set_ylabel("Correlation")
    axes[2].legend(fontsize=8, loc="upper left")
    ax2.legend(fontsize=8, loc="upper right")
    for ax in axes:
        ax.set_xticks(x, names, fontsize=8)
        ax.set_xlabel("Rest weighting")
    fig.suptitle(f"v2 sensitivity to the rest weighting ({kind}); w_rest = sqrt(λ N_curve / N_rest)", fontsize=10)
    fig.tight_layout()
    return savefig(fig, f"v2_weighting_sensitivity_{kind.lower()}")


def make_figures(headline_kind: str = "DFN", kinds=("SPMe", "DFN")):
    from .figures import fig_cv, fig_headline, fig_overpotentials

    paths = [fig_rests_measured()]
    for kind in kinds:
        if (RESULTS_DIR / f"v2_fit_{kind.lower()}.json").exists():
            paths.append(fig_fit_rests(kind))
        if (RESULTS_DIR / f"v2_weighting_sensitivity_{kind.lower()}.json").exists():
            paths.append(fig_weighting(kind))
        if (RESULTS_DIR / f"v2_validation_{kind.lower()}.json").exists():
            paths += [fig_overpotentials(kind, "v2_"), fig_cv(kind, "v2_")]
        if (RESULTS_DIR / f"v2_1_fit_{kind.lower()}.json").exists():
            paths.append(fig_fit_rests(kind, "v2_1_"))
        if (RESULTS_DIR / f"v2_1_validation_{kind.lower()}.json").exists():
            paths += [fig_overpotentials(kind, "v2_1_"), fig_cv(kind, "v2_1_")]
    if (RESULTS_DIR / f"v2_validation_{headline_kind.lower()}.json").exists():
        paths.append(fig_headline(headline_kind, "v2_"))
    if (RESULTS_DIR / f"v2_1_validation_{headline_kind.lower()}.json").exists():
        paths.append(fig_headline(headline_kind, "v2_1_"))
    if (RESULTS_DIR / "v2_comparison.json").exists() and headline_kind in json.load(open(RESULTS_DIR / "v2_comparison.json"))["models"]:
        paths.append(fig_comparison(headline_kind))
    return paths
