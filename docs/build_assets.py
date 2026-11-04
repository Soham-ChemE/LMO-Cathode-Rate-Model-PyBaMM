"""Build every asset used by docs/index.html and the README hero from repository data.

Reads (all produced by `python scripts/run_all.py`, or shipped in the repo):
  data/raw/*.txt                      EC-Lab exports, parsed with src/lmo/io.py
  data/processed/model_gcd_dfn.csv    DFN v1 galvanostatic curves at all four rates
  data/processed/overpotentials_dfn.csv  DFN v1 voltage-loss components along each discharge
  data/processed/model_cv_dfn.csv     DFN v1 cyclic voltammograms
  results/*.json                      every number quoted on the page
Writes:
  docs/assets/data.js                 window.LMO = {...} (classic script, works from file://)
  docs/assets/fig_*_{dark,light}.svg  figures for the page
  docs/assets/social-preview.png      1280x640
Run from the repository root with the package installed:  python docs/build_assets.py
"""
import json, os, sys, subprocess
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from lmo import REPRESENTATIVE_CYCLE, RATES, CVS  # noqa: E402
from lmo.io import get_step, cv_cycle  # noqa: E402

OUT = os.path.join(ROOT, "docs", "assets"); os.makedirs(OUT, exist_ok=True)
R = lambda n: json.load(open(os.path.join(ROOT, "results", n)))
val, fit, cv, ov = R("validation_dfn.json"), R("fit_dfn.json"), R("cv_dfn.json"), R("overpotentials_dfn.json")
meas = R("measured_summary.json"); rests = R("v2_rests_measured.json"); v2 = R("v2_fit_dfn.json"); v21 = R("v2_1_fit_dfn.json")
ORDER = ["C/10", "C/5", "C/2", "2C"]

def thin(x, y, n=140):
    idx = np.unique(np.linspace(0, len(x) - 1, n).astype(int))
    return [[round(float(x[i]), 4), round(float(y[i]), 4)] for i in idx]

# ---------------- measured and modelled discharge curves ----------------
model = pd.read_csv(os.path.join(ROOT, "data", "processed", "model_gcd_dfn.csv"))
curves = []
for rate in ORDER:
    cyc = REPRESENTATIVE_CYCLE[rate]
    s = get_step(rate, cyc, "discharge")
    m = model[(model["rate"] == rate) & (model["cycle"] == cyc) & (model["step"] == "discharge")]
    e = val["rates"][rate]
    curves.append({"rate": rate, "cycle": cyc, "role": e["role"], "measured": thin(s.capacity_mAh, s.voltage_V),
                   "model": thin(m["capacity_mAh"].to_numpy(), m["voltage_V"].to_numpy()),
                   "q_meas": round(e["discharge"]["capacity_mAh_measured"], 3), "q_model": round(e["discharge"]["capacity_mAh_model"], 3),
                   "q_err_pct": round(e["discharge"]["capacity_error_pct"], 1), "v_meas": round(e["discharge"]["mean_voltage_V_measured"], 3),
                   "v_model": round(e["discharge"]["mean_voltage_V_model"], 3), "e_meas": round(e["discharge"]["energy_mWh_measured"], 2),
                   "e_model": round(e["discharge"]["energy_mWh_model"], 2), "vq_rmse_mV": round(e["discharge"]["vq_rmse_mV"]),
                   "current_mA": round(abs(s.mean_current_mA), 3)})

# ---------------- measured efficiencies (representative cycles) ----------------
g = meas["gcd"]; eff = {}
for rate in ORDER:
    cyc = REPRESENTATIVE_CYCLE[rate]
    ch = next(r for r in g if r["rate"] == rate and r["cycle"] == cyc and r["step"] == "charge")
    dc = next(r for r in g if r["rate"] == rate and r["cycle"] == cyc and r["step"] == "discharge")
    eff[rate] = {"rte_pct": round(100 * dc["energy_mWh"] / ch["energy_mWh"], 1), "ce_pct": round(100 * dc["capacity_mAh"] / ch["capacity_mAh"], 1)}
c10 = next(r for r in g if r["rate"] == "C/10" and r["cycle"] == 1 and r["step"] == "discharge")
tc = next(r for r in g if r["rate"] == "2C" and r["cycle"] == 2 and r["step"] == "discharge")
retention = {"capacity_pct": round(100 * tc["capacity_mAh"] / c10["capacity_mAh"], 1), "energy_pct": round(100 * tc["energy_mWh"] / c10["energy_mWh"], 1)}

# ---------------- overpotentials along the 2C discharge ----------------
op = pd.read_csv(os.path.join(ROOT, "data", "processed", "overpotentials_dfn.csv"))
o2 = op[(op["rate"] == "2C") & (op["step"] == "discharge")]
comps = ["contact_ohmic", "reaction", "particle_diffusion", "electrolyte_concentration", "electrolyte_ohmic", "solid_ohmic"]
idx = np.unique(np.linspace(0, len(o2) - 1, 120).astype(int))
op_profile = {"capacity_mAh": [round(float(v), 4) for v in o2["capacity_mAh"].to_numpy()[idx]],
              **{c: [round(1000 * float(v), 1) for v in o2[c].to_numpy()[idx]] for c in comps},
              "ocv_V": [round(float(v), 4) for v in o2["ocv_V"].to_numpy()[idx]], "voltage_V": [round(float(v), 4) for v in o2["voltage_V"].to_numpy()[idx]]}
two = ov["two_c"]
shares = {c: round(100 * two["shares"][c], 1) for c in comps}

# ---------------- CV peaks ----------------
cv_rows = []
for scan, e in cv["scans"].items():
    cv_rows.append({"scan": scan, "ox_meas": round(e["measured"]["oxidation_peak_V"], 2), "ox_model": round(e["model"]["oxidation_peak_V"], 2),
                    "red_meas": round(e["measured"]["reduction_peak_V"], 2), "red_model": round(e["model"]["reduction_peak_V"], 2),
                    "ox_shift_mV": round(abs(e["oxidation_peak_shift_mV"])), "red_shift_mV": round(abs(e["reduction_peak_shift_mV"]))})
cv_max_shift = max(max(r["ox_shift_mV"], r["red_shift_mV"]) for r in cv_rows)
mcv = pd.read_csv(os.path.join(ROOT, "data", "processed", "model_cv_dfn.csv"))
cv_curves = []
for scan in CVS:
    d = cv_cycle(scan, 2); mm = mcv[(mcv["scan"] == scan) & (mcv["cycle"] == 2)]
    cv_curves.append({"scan": scan, "measured": thin(d["voltage_V"].to_numpy(), d["current_mA"].to_numpy(), 400),
                      "model": thin(mm["voltage_V"].to_numpy(), mm["current_mA"].to_numpy(), 400)})

# ---------------- tests ----------------
n_tests = None
try:
    out = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q"], cwd=ROOT, capture_output=True, text=True, timeout=600).stdout
    n_tests = int([l for l in out.splitlines() if "collected" in l][-1].split()[0])
except Exception as ex:  # pragma: no cover
    print("could not count tests:", ex)

p = fit["parameters"]; h = val["headline"]
stats = {"train_rmse_mV": round(1000 * fit["training_rmse_V"]), "n_params": 3,
         "c10_cap_err_pct": round(abs(val["rates"]["C/10"]["discharge"]["capacity_error_pct"]), 1),
         "two_c_v_err_mV": round(abs(h["mean_voltage_error_mV"])), "two_c_cap_err_pct": round(h["error_pct"], 1),
         "two_c_pred_mAh": round(h["predicted_mAh"], 3), "two_c_meas_mAh": round(h["measured_mAh"], 3),
         "cv_max_shift_mV": cv_max_shift, "shares": shares, "two_c_mean_loss_mV": round(two["mean_loss_mV"]),
         "tau_s": round(p["tau_R2_over_Ds_s"]["value"]), "Ds": p["Ds"]["value"], "k": p["k"]["value"], "R_contact": round(p["R_contact"]["value"], 1),
         "R_contact_loo": [round(p["R_contact"]["leave_one_step_out_min"]), round(p["R_contact"]["leave_one_step_out_max"])],
         "rest_after_charge_ohm": round(rests["after_charge_mean_ohm_c10_c5_c2"]), "rest_after_discharge_ohm": round(rests["after_discharge_mean_ohm_c10_c5_c2"]),
         "rest_2c_after_charge_ohm": round(rests["after_charge_mean_ohm_2c_post_hoc"]),
         "v2_rest_rmse_mV": round(1000 * v2["training_rmse_rests_V"]), "v21_rest_rmse_mV": round(1000 * v21["training_rmse_rests_V"]),
         "v21_R_contact": round(v21["parameters"]["R_contact"]["value"]), "n_tests": n_tests, "retention": retention, "eff": eff,
         "b_meas": [round(cv["b_value"]["measured_oxidation"]["b"], 2), round(cv["b_value"]["measured_reduction"]["b"], 2)],
         "b_model": [round(cv["b_value"]["model_oxidation"]["b"], 2), round(cv["b_value"]["model_reduction"]["b"], 2)]}
rest_rows = [{"rate": r, **{k: round(v, 0) if isinstance(v, float) else v for k, v in rests["by_rate"][r].items() if k in ("after_charge_mean_ohm", "after_discharge_mean_ohm", "used_in_v2_fit")}} for r in ORDER]

data = {"curves": curves, "op_profile": op_profile, "cv_rows": cv_rows, "cv_curves": cv_curves, "stats": stats, "rest_rows": rest_rows,
        "cell": {"coating_um": 85, "foil_um": 15, "particle_um": 5, "disc_mm": 10, "active_wt_pct": 94.2, "rated_mAh_g": 110, "capacity_mAh": 1.434}}
with open(os.path.join(OUT, "data.js"), "w") as f:
    f.write("// Generated by docs/build_assets.py from results/*.json and data/. Do not edit by hand.\n")
    f.write("window.LMO = " + json.dumps(data, separators=(",", ":")) + ";\n")
print("data.js:", {k: stats[k] for k in ("train_rmse_mV", "two_c_v_err_mV", "cv_max_shift_mV", "n_tests", "retention")})

# ---------------- figures ----------------
THEMES = {"dark": ("#e6edf3", "#8b98a9"), "light": ("#0b1220", "#5b6676")}; DARK = "#0b1220"
RATE_COL = {"C/10": "#22d3ee", "C/5": "#34d399", "C/2": "#f59e0b", "2C": "#f87171"}
COMP_COL = {"contact_ohmic": "#f87171", "reaction": "#f59e0b", "particle_diffusion": "#34d399", "electrolyte_concentration": "#22d3ee", "electrolyte_ohmic": "#60a5fa", "solid_ohmic": "#a78bfa"}
COMP_LAB = {"contact_ohmic": "contact / ohmic", "reaction": "charge transfer", "particle_diffusion": "particle diffusion", "electrolyte_concentration": "electrolyte conc.", "electrolyte_ohmic": "electrolyte ohmic", "solid_ohmic": "solid ohmic"}
def style(ax, INK, MUTED):
    ax.set_facecolor("none")
    for s in ax.spines.values(): s.set_color(MUTED); s.set_linewidth(0.8)
    ax.tick_params(colors=INK, labelsize=9); ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK); ax.title.set_color(INK); ax.grid(True, color=MUTED, alpha=0.2, lw=0.6)
for THEME, (INK, MUTED) in THEMES.items():
    save = lambda fig, name: (fig.savefig(os.path.join(OUT, f"{name}_{THEME}.svg"), transparent=True, bbox_inches="tight"), plt.close(fig), print("wrote", name, THEME))
    fig, ax = plt.subplots(figsize=(6.6, 4.2)); style(ax, INK, MUTED)
    for c in curves:
        mq, mv = zip(*c["measured"]); oq, ovv = zip(*c["model"])
        ax.plot(mq, mv, color=RATE_COL[c["rate"]], lw=2, label=f"{c['rate']} measured ({c['role']})"); ax.plot(oq, ovv, color=RATE_COL[c["rate"]], lw=1.4, ls="--")
    ax.set_xlabel("Discharge capacity (mAh)"); ax.set_ylabel("Cell voltage (V vs Li/Li⁺)"); ax.set_ylim(3.45, 4.55)
    ax.legend(fontsize=8, frameon=False, labelcolor=INK, loc="lower left"); ax.set_title("Solid: measured. Dashed: DFN (3 parameters fitted on C/5 and C/2 only)", fontsize=9)
    save(fig, "fig_curves")
    fig, ax = plt.subplots(figsize=(6.6, 3.8)); style(ax, INK, MUTED)
    q = np.array(op_profile["capacity_mAh"]); base = np.zeros_like(q)
    for c in comps:
        y = np.array(op_profile[c]); ax.fill_between(q, base, base + y, color=COMP_COL[c], alpha=0.85, lw=0, label=f"{COMP_LAB[c]} ({shares[c]:.0f} %)"); base = base + y
    ax.set_xlabel("2C discharge capacity (mAh)"); ax.set_ylabel("Voltage loss (mV)"); ax.legend(fontsize=8, frameon=False, labelcolor=INK, loc="upper left")
    ax.set_title(f"Where the 2C voltage goes: capacity-weighted mean loss {stats['two_c_mean_loss_mV']} mV", fontsize=9); save(fig, "fig_overpotential")
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.9), sharey=True)
    for ax, c in zip(axs, cv_curves):
        style(ax, INK, MUTED); mx, my = zip(*c["measured"]); ox, oy = zip(*c["model"])
        ax.plot(mx, my, color="#22d3ee", lw=1.6, label="measured"); ax.plot(ox, oy, color="#f59e0b", lw=1.2, ls="--", label="DFN"); ax.set_title(c["scan"], fontsize=9); ax.set_xlabel("V vs Li/Li⁺")
    axs[0].set_ylabel("Current (mA)"); axs[0].legend(fontsize=8, frameon=False, labelcolor=INK); fig.tight_layout(); save(fig, "fig_cv")

fig = plt.figure(figsize=(12.8, 6.4), dpi=100); fig.patch.set_facecolor(DARK)
ax = fig.add_axes([0.53, 0.14, 0.44, 0.74]); ax.set_facecolor(DARK)
for c in curves:
    mq, mv = zip(*c["measured"]); oq, ovv = zip(*c["model"]); ax.plot(mq, mv, color=RATE_COL[c["rate"]], lw=2.4); ax.plot(oq, ovv, color=RATE_COL[c["rate"]], lw=1.4, ls="--")
    ax.text(mq[len(mq)//3], mv[len(mv)//3] + 0.04, c["rate"], color=RATE_COL[c["rate"]], fontsize=11, fontweight="bold")
for s in ax.spines.values(): s.set_color("#1f2a3c")
ax.tick_params(colors="#8b98a9"); ax.set_xlabel("Discharge capacity (mAh)", color="#8b98a9"); ax.set_ylabel("V vs Li/Li⁺", color="#8b98a9"); ax.set_ylim(3.45, 4.55)
fig.text(0.05, 0.74, "LiMn₂O₄ rate model", color="#e6edf3", fontsize=40, fontweight="bold")
fig.text(0.05, 0.62, "PyBaMM DFN half cell fitted to my own\nBio-Logic data, validated blind at 2C", color="#e6edf3", fontsize=19, va="top")
fig.text(0.05, 0.36, f"{stats['train_rmse_mV']} mV training RMSE, 3 fitted parameters\n2C mean voltage within {stats['two_c_v_err_mV']} mV, never fitted\nsix CV peaks within {stats['cv_max_shift_mV']} mV", color="#22d3ee", fontsize=14, va="top")
fig.text(0.05, 0.10, "Soham Kavathekar  |  github.com/Soham-ChemE", color="#8b98a9", fontsize=12)
fig.savefig(os.path.join(OUT, "social-preview.png"), dpi=100, facecolor=DARK); plt.close(fig); print("wrote social-preview.png")
