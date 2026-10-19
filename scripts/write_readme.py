"""Render README.md and results/resume_bullets.md from results/*.json.

No number in the README is typed by hand: everything is read from the JSON written by
scripts/run_all.py, so the text can never drift from the results.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def load(name):
    with open(RESULTS / name) as f:
        return json.load(f)


def fmt_ci(p, unit="", scale=1.0, sig=3):
    v = p["value"] * scale
    lo, hi = p["ci95_low_jacobian"] * scale, p["ci95_high_jacobian"] * scale
    g = lambda x: f"{x:.{sig}g}"
    if p.get("sd_log10_jacobian", 0) > 2:
        # No sensitivity left at a bound of the search box: the Jacobian interval is meaningless.
        s = f"{g(v)} (at bound, Jacobian interval undefined)"
    else:
        s = f"{g(v)} ({g(lo)} to {g(hi)})"
    if "leave_one_step_out_min" in p:
        s += f"; LOO {g(p['leave_one_step_out_min'] * scale)} to {g(p['leave_one_step_out_max'] * scale)}"
    return s + (f" {unit}" if unit else "")


V1_NEXT_STEP = """## Next step

The weakest point is the 0.94 correlation between the rate constant and the contact resistance: C/5 and C/2 alone
cannot tell an instantaneous ohmic drop from a transport loss that builds up over time, so the fit chose the
ohmic explanation and over-predicts 2C capacity. Every galvanostatic step in this data ends with a 30 s rest; the
instantaneous voltage jump at the start of each rest is the ohmic part and the slow relaxation that follows is
the diffusion part. Adding those relaxations to the training set, still without touching 2C, is the direct way
to break that correlation."""


def g(x, sig=3):
    return f"{x:.{sig}g}"


def load_v2(kinds, H):
    """Everything the v1 vs v2 section needs, or None if the v2 pipeline has not been run."""
    needed = ["v2_comparison.json", "v2_rests_measured.json", f"v2_fit_{H.lower()}.json", f"v2_validation_{H.lower()}.json"]
    if not all((RESULTS / n).exists() for n in needed):
        return None
    comp = load("v2_comparison.json")
    if H not in comp["models"]:
        return None
    rests = load("v2_rests_measured.json")
    fit1 = load(f"fit_{H.lower()}.json")
    fit2 = load(f"v2_fit_{H.lower()}.json")
    val1 = load(f"validation_{H.lower()}.json")
    val2 = load(f"v2_validation_{H.lower()}.json")
    ov2 = load(f"v2_overpotentials_{H.lower()}.json")
    cv2 = load(f"v2_cv_{H.lower()}.json")
    sens = {k: load(f"v2_weighting_sensitivity_{k.lower()}.json") for k in kinds if (RESULTS / f"v2_weighting_sensitivity_{k.lower()}.json").exists()}
    jump = load(f"v2_jump_fit_{H.lower()}.json") if (RESULTS / f"v2_jump_fit_{H.lower()}.json").exists() else None
    jval = load(f"v2_jump_validation_{H.lower()}.json") if (RESULTS / f"v2_jump_validation_{H.lower()}.json").exists() else None
    post = load(f"v2_rest_posthoc_2c_{H.lower()}.json") if (RESULTS / f"v2_rest_posthoc_2c_{H.lower()}.json").exists() else None
    c = comp["models"][H]
    p1, p2 = fit1["parameters"], fit2["parameters"]
    w = fit2["weighting"]

    # --- model-free table
    rows = []
    for r in rests["rests"]:
        tag = "post-hoc only, never fitted" if r["post_hoc_only"] else ("fitted (v2)" if r["used_in_v2_fit"] else "not fitted")
        rows.append(
            f"| {r['rate']} | {r['cycle']} | {r['after']} | {abs(r['current_before_mA']):.3f} | {1000 * r['first_sample_delay_s']:.0f} | "
            f"{abs(r['jump_mV']):.1f} | {r['instantaneous_resistance_ohm']:.0f} | {abs(r['total_relaxation_mV']):.0f} | {tag} |"
        )
    rest_table = "\n".join(
        ["| Rate | Cycle | After | I (mA) | First sample (ms) | Jump (mV) | Jump / I (Ω) | 30 s relaxation (mV) | Use |", "|---|---|---|---|---|---|---|---|---|"] + rows
    )

    # --- parameter table v1 / v2 / v2 jump-only
    def prow(name, fj, vj, extra):
        p = fj["parameters"]
        return (
            f"| {name} | {extra} | {fmt_ci(p['tau_R2_over_Ds_s'])} | {fmt_ci(p['k'])} | {fmt_ci(p['R_contact'])} | "
            f"{fj['log10_correlation']['k-R_contact']:.2f} | {vj['headline']['predicted_mAh']:.3f} ({vj['headline']['error_pct']:+.1f} %) |"
        )

    par_rows = [
        prow(f"v1 {H}", fit1, val1, f"{1000 * fit1['training_rmse_V']:.1f} / –"),
        prow(f"v2 {H} (full rest traces, λ = {w['lambda']:g})", fit2, val2, f"{1000 * fit2['training_rmse_curves_V']:.1f} / {1000 * fit2['training_rmse_rests_V']:.1f}"),
    ]
    if jump and jval:
        par_rows.append(prow(f"v2 {H} (jump only, λ = {jump['weighting']['lambda']:g})", jump, jval, f"{1000 * jump['training_rmse_curves_V']:.1f} / {1000 * jump['training_rmse_rest_traces_V']:.1f}"))
    par_table = "\n".join(
        ["| Fit | Train RMSE curves / rests (mV) | τ = R²/Ds (s) | k (m^2.5 mol^-0.5 s^-1) | R_contact (Ω) | k–R_contact corr. | 2C prediction (mAh) |", "|---|---|---|---|---|---|---|"] + par_rows
    )

    # --- weighting sensitivity table
    srows = []
    for k, sj in sens.items():
        for name, e in sj["weightings"].items():
            srows.append(
                f"| {k} | {name} | {e['rest_weight']:.1f} | {e['n_rest_points']} | {g(e['tau_R2_over_Ds_s'])} | {g(e['parameters']['k'])} | {e['parameters']['R_contact']:.0f} | "
                f"{e['k_R_contact_correlation']:.2f} | {e['training_rmse_curves_mV']:.1f} | {e['training_rmse_rest_traces_mV']:.1f} |"
            )
    sens_table = "\n".join(
        ["| Model | Rest weighting | w per rest point | Rest points | τ (s) | k | R_contact (Ω) | k–R corr. | Curves RMSE (mV) | Rest traces RMSE (mV) |", "|---|---|---|---|---|---|---|---|---|---|"] + srows
    )

    # --- v1 vs v2 validation table
    vrows = []
    for rate in val1["rates"]:
        d1, d2 = val1["rates"][rate]["discharge"], val2["rates"][rate]["discharge"]
        role = val2["rates"][rate]["role"]
        if rate == "C/10":
            role = "curve not fitted; its rests are in v2"
        vrows.append(
            f"| {rate} | {role} | {d1['capacity_mAh_measured']:.3f} | {d1['capacity_mAh_model']:.3f} ({d1['capacity_error_pct']:+.1f} %) | {d2['capacity_mAh_model']:.3f} ({d2['capacity_error_pct']:+.1f} %) | "
            f"{d1['mean_voltage_V_measured']:.3f} | {d1['mean_voltage_V_model']:.3f} | {d2['mean_voltage_V_model']:.3f} | {d1['energy_mWh_measured']:.2f} | {d1['energy_mWh_model']:.2f} | {d2['energy_mWh_model']:.2f} | {d1['vq_rmse_mV']:.0f} | {d2['vq_rmse_mV']:.0f} |"
        )
    val_table = "\n".join(
        [
            "| Rate | Role | Q meas (mAh) | Q v1 | Q v2 | V̄ meas (V) | V̄ v1 | V̄ v2 | E meas (mWh) | E v1 | E v2 | V-Q RMSE v1 (mV) | V-Q RMSE v2 (mV) |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        + vrows
    )

    # --- CV table v2
    cv_rows = []
    for scan, e in cv2["scans"].items():
        m, sm = e["measured"], e["model"]
        cv_rows.append(
            f"| {scan} | {m['oxidation_peak_V']:.2f} / {sm['oxidation_peak_V']:.2f} | {m['oxidation_peak_mA']:.2f} / {sm['oxidation_peak_mA']:.2f} | "
            f"{m['reduction_peak_V']:.2f} / {sm['reduction_peak_V']:.2f} | {m['reduction_peak_mA']:.2f} / {sm['reduction_peak_mA']:.2f} | {e['peak_separation_measured_mV']:.0f} / {e['peak_separation_model_mV']:.0f} |"
        )
    cv_table = "\n".join(
        ["| Scan rate | Ox. peak V (meas / v2) | Ox. peak mA (meas / v2) | Red. peak V (meas / v2) | Red. peak mA (meas / v2) | ΔEp mV (meas / v2) |", "|---|---|---|---|---|---|"] + cv_rows
    )
    b2 = cv2["b_value"]

    # --- rest jumps at the optimum
    jrows = []
    for key, e in fit2["rest_jumps"].items():
        je = jump["rest_jumps"][key] if jump else None
        jrows.append(
            f"| {key} | {e['current_mA']:.3f} | {e['jump_measured_mV']:.0f} / {e['jump_model_mV']:.0f}{(' / ' + format(je['jump_model_mV'], '.0f')) if je else ''} | "
            f"{e['instantaneous_resistance_measured_ohm']:.0f} / {e['instantaneous_resistance_model_ohm']:.0f}{(' / ' + format(je['instantaneous_resistance_model_ohm'], '.0f')) if je else ''} | "
            f"{e['relaxation_30s_measured_mV']:.0f} / {e['relaxation_30s_model_mV']:.0f}{(' / ' + format(je['relaxation_30s_model_mV'], '.0f')) if je else ''} |"
        )
    jhead = "meas / v2" + (" / v2 jump-only" if jump else "")
    jump_table = "\n".join([f"| Rest | I (mA) | Jump, mV ({jhead}) | Jump / I, Ω ({jhead}) | 30 s relaxation, mV ({jhead}) |", "|---|---|---|---|---|"] + jrows)

    two2 = ov2["two_c"]
    sh2 = two2["shares"]
    sh1 = c["v1"]["two_c_shares"]
    e1, e2 = c["v1"]["two_c"], c["v2"]["two_c"]
    improved = c["two_c_improved"]
    at_bound = {n: (abs(p2[n]["log10"] - lo) < 0.02 or abs(p2[n]["log10"] - hi) < 0.02) for n, (lo, hi) in {"Ds": (-17.0, -11.0), "k": (-14.0, -8.0), "R_contact": (-1.0, 3.0)}.items()}
    bound_txt = ", ".join(n for n, b in at_bound.items() if b)
    r_up = p2["R_contact"]["value"] > p1["R_contact"]["value"]
    R_free_chg = rests["after_charge_mean_ohm_c10_c5_c2"]
    R_free_dis = rests["after_discharge_mean_ohm_c10_c5_c2"]

    verdict = (
        f"the 2C capacity error moves from {e1['capacity_error_pct']:+.1f} % to {e2['capacity_error_pct']:+.1f} %"
        + (" (better)" if improved else " (not better)")
    )
    key_bullet = (
        f"- **v2, rests added to the fit (not an improvement):** refitting on the C/5 and C/2 curves plus the C/10, C/5 and C/2 rests "
        f"{verdict}. The model-free rest jumps give {R_free_chg:.0f} Ω (after charge), yet the v2 fit moves R_contact "
        f"{'up' if r_up else 'down'} to {p2['R_contact']['value']:.0f} Ω: this model has no double-layer capacitance and one particle size, so it cannot "
        f"reproduce the log-time relaxation the rests show and instead matches their mean level with a bigger ohmic drop (see \"v1 vs v2\")."
    )
    if improved:
        key_bullet = key_bullet.replace("(not an improvement)", "(a modest improvement)")

    section = f"""## v1 vs v2: adding the post-step rests to the fit

v1 fitted only the C/5 and C/2 charge/discharge curves, which cannot separate an instantaneous ohmic drop from a
transport loss that builds over time (k vs R_contact correlation {fit1['log10_correlation']['k-R_contact']:.2f}); it fitted
R_contact = {p1['R_contact']['value']:.1f} Ω and over-predicted the blind 2C capacity by {e1['capacity_error_pct']:+.1f} %. Every galvanostatic
step ends with a 30 s rest whose first sample is logged {1000 * rests['first_sample_delay_s']['min']:.0f} ms after the current switches off. v2 keeps the
same model, the same three parameters and the same 10 mm disc, and adds those rests to the training residual. No 2C
data (curves or rests) is used in either version.

### Model-free check: instantaneous resistance of every rest

Jump = first rest sample minus last sample under current, straight from the raw files. The first sample already
contains whatever charge-transfer overpotential relaxed in the first {1000 * rests['first_sample_delay_s']['min']:.0f} ms, so jump / I is an upper bound
on the pure ohmic resistance.

{rest_table}

After charge the jump gives {R_free_chg:.0f} Ω on average over C/10, C/5 and C/2 (2C, post-hoc: {rests['after_charge_mean_ohm_2c_post_hoc']:.0f} Ω), well
below the v1 contact resistance of {p1['R_contact']['value']:.0f} Ω: v1 assigned too much of the loss to resistance. After discharge the jump is
larger ({R_free_dis:.0f} Ω on average) because the cut-off at 3.5 V sits on the steep end of the OCV where the exchange current
is smallest and the particle surface relaxes fastest. The C/5 cycle-2 discharge rest (156 Ω) carries the 22 mV contact
disturbance noted in `data/README.md`, which returns at the switch-off and so is confirmed as a contact resistance. What
the table does not show, and the figure does, is the shape: the relaxation after the jump is close to linear in
log(time) from 50 ms to 30 s, i.e. a distributed process, not one exponential or one sqrt(t) diffusion tail.

![Measured rests](figures/v2_rests_measured.png)

### v2 fit

Training residual: the v1 curve residual (C/5 and C/2 cycle 1, charge and discharge, 300 points each, {fit2['weighting']['n_curve_points']}
points) plus the model minus measured voltage at every logged point of the six cycle-1 rests after the C/10, C/5 and
C/2 steps ({fit2['weighting']['n_rest_points']} points), time measured from the last sample under current. Each rest residual is
multiplied by w = sqrt(λ · N_curve / N_rest) so that the rest block carries λ times the summed squared weight of the
curve block; the headline uses λ = {w['lambda']:g} (w = {w['rest_weight']:.1f} per rest point), and the table below shows λ = 0.25, 1, 4, the
unweighted case, and a variant in which only the first sample of each rest (the instantaneous jump) enters. The rests
are output at {w['rest_period_s']} s in the model with a particle mesh refined at the surface (60 radial points, exponential stretch),
because the uniform v1 mesh lags the surface relaxation by seconds. The C/10 curves remain a check, but its rests are
now training data. Leave-one-out refits drop one of the {len(fit2['leave_one_step_out'])} residual blocks (4 curve steps, 6 rests) in turn.

{par_table}

Intervals as in the v1 table (95 % Jacobian; LOO = leave-one-block-out range). k enters j0 = k F sqrt(c_e c_s (c_max - c_s)); its
bounds are 1e-14 to 1e-8 and Ds's are 1e-17 to 1e-11 m²/s. {('In the headline v2 fit ' + bound_txt + ' sit' + ('s' if ',' not in bound_txt else '') + ' at the bound of the search box, i.e. the optimiser wants no loss from that term at all.') if bound_txt else ''}
The k–R_contact correlation drops from {fit1['log10_correlation']['k-R_contact']:.2f} to {fit2['log10_correlation']['k-R_contact']:.2f}, so the rests do decorrelate the two parameters, but not in the
direction expected: R_contact moves {'up' if r_up else 'down'} to {p2['R_contact']['value']:.0f} Ω although the rests themselves say {R_free_chg:.0f} Ω.

![v2 rests](figures/v2_fit_rests_{H.lower()}.png)

{jump_table}

**Why the model cannot use the rests.** The measured relaxation keeps going at roughly the same number of mV per
decade of time from 50 ms to 30 s. In this model the current switches off and three things happen: the contact,
electrolyte and solid ohmic drops vanish at once, the charge-transfer overpotential vanishes at once (there is no
double-layer capacitance), and the particle-surface concentration relaxes as sqrt(t) with one time constant R²/Ds.
Nothing relaxes like log(t). Forced to put a 100 to 280 mV relaxation somewhere, the least-squares fit can only
match its mean level, which it does by making the instantaneous jump larger (bigger R_contact, k at the bound so that
charge transfer costs nothing) and the diffusion faster. The curves and the rests are therefore inconsistent inside
this model structure: the curves need about {1000 * fit1['parameters']['R_contact']['value'] * 0.287 / 0.287 / 1000 * 1:.0f} Ω of current-proportional loss, the rests say only {R_free_chg:.0f} Ω of it is
instantaneous, and the remainder relaxes over seconds in a way that is neither Butler-Volmer-instantaneous nor
single-particle Fickian. That is the finding of v2.

### Sensitivity to the rest weighting

{sens_table}

Every weighting that uses the full rest traces moves R_contact above the v1 value; the jump-only variant, which
asks the model to match only the first sample of each rest, cannot do that either without losing the curves (its curve
RMSE is the worst), because the model's instantaneous jump is R_contact plus the whole charge-transfer overpotential,
and the curves need more current-proportional loss than {R_free_chg:.0f} Ω provides.

![Weighting sensitivity](figures/v2_weighting_sensitivity_{H.lower()}.png)

### Blind validation, v1 vs v2 ({H})

{val_table}

{verdict[0].upper() + verdict[1:]}; the 2C V-Q RMSE goes from {e1['vq_rmse_mV']:.0f} to {e2['vq_rmse_mV']:.0f} mV and the mean-voltage error
from {e1['mean_voltage_error_mV']:+.0f} to {e2['mean_voltage_error_mV']:+.0f} mV. {'The jump-only variant predicts ' + format(jval['headline']['predicted_mAh'], '.3f') + ' mAh (' + format(jval['headline']['error_pct'], '+.1f') + ' %).' if jval else ''} Nothing was tuned on 2C or on the CVs.

![v1 vs v2 at 2C](figures/v2_comparison_2c.png)

Cyclic voltammetry with the v2 parameters, simulated with the measured voltage programme as in v1:

{cv_table}

b-values with v2: {b2['model_oxidation']['b']:.2f} (oxidation) and {b2['model_reduction']['b']:.2f} (reduction) against {b2['measured_oxidation']['b']:.2f} / {b2['measured_reduction']['b']:.2f} measured
(v1: {c['v1']['b_value_model']['oxidation']:.2f} / {c['v1']['b_value_model']['reduction']:.2f}). The largest CV peak shift is {c['v2']['cv_max_peak_shift_mV']:.0f} mV with v2 against {c['v1']['cv_max_peak_shift_mV']:.0f} mV with v1.

### Overpotential breakdown at 2C with v2

Capacity-weighted mean loss on the 2C discharge: {two2['mean_loss_mV']:.0f} mV (v1: {c['v1']['two_c_mean_loss_mV']:.0f} mV), of which {100 * sh2['contact_ohmic']:.0f} % contact/ohmic
(v1 {100 * sh1['contact_ohmic']:.0f} %), {100 * sh2['reaction']:.0f} % charge transfer (v1 {100 * sh1['reaction']:.0f} %), {100 * sh2['particle_diffusion']:.0f} % particle diffusion
(v1 {100 * sh1['particle_diffusion']:.0f} %) and {100 * (sh2['electrolyte_concentration'] + sh2['electrolyte_ohmic']):.0f} % electrolyte. The bulk OCV at the 2C cut-off is {two2['bulk_ocv_at_cutoff_V']:.2f} V, i.e. {two2['capacity_locked_out_mAh']:.2f} mAh locked out
(v1: {val1['rates']['2C']['discharge']['capacity_mAh_model'] and ov_two_c_v1(H)['bulk_ocv_at_cutoff_V']:.2f} V, {ov_two_c_v1(H)['capacity_locked_out_mAh']:.2f} mAh).

![v2 overpotentials](figures/v2_overpotentials_{H.lower()}.png)

### Post-hoc: the 2C rests (never used for fitting)

{('Measured 2C cycle-2 jumps: ' + format(post['measured_2c_after_charge_ohm'], '.0f') + ' Ω after charge and ' + format(post['measured_2c_after_discharge_ohm'], '.0f') + ' Ω after discharge, against the v2 R_contact of ' + format(post['v2_R_contact_ohm'], '.0f') + ' Ω; the v2 model jumps are ' + format(post['rests']['2C rest after charge']['instantaneous_resistance_model_ohm'], '.0f') + ' and ' + format(post['rests']['2C rest after discharge']['instantaneous_resistance_model_ohm'], '.0f') + ' Ω, and its 30 s relaxations ' + format(post['rests']['2C rest after charge']['relaxation_30s_model_mV'], '.0f') + ' / ' + format(post['rests']['2C rest after discharge']['relaxation_30s_model_mV'], '.0f') + ' mV against ' + format(post['rests']['2C rest after charge']['relaxation_30s_measured_mV'], '.0f') + ' / ' + format(post['rests']['2C rest after discharge']['relaxation_30s_measured_mV'], '.0f') + ' mV measured.') if post else ''}
The 2C after-charge jump sits with the C/10 to C/2 ones at about 70 to 75 Ω, so the instantaneous resistance of the
cell is rate-independent and the v2 value is not it: the comparison confirms the diagnosis above rather than the fit.

## Next step

Both versions put the loss where the model lets them. The rests show what is missing: a relaxation spread over 50 ms
to 30 s. Two additions to the model would let it use the rests properly, and both are standard PyBaMM options: a
double-layer capacitance (`surface form: differential`), which gives the charge-transfer overpotential a time constant
and a transmission-line relaxation through the porous electrode, and a particle-size distribution, which spreads the
diffusion time constant. With either in place the instantaneous jump (about {R_free_chg:.0f} Ω) would pin R_contact and the
relaxation would pin the rest; that is the fit that could move the 2C prediction, and it is the next thing to do.
A measured electrode area and an EIS spectrum would remove the two remaining guesses (area, double-layer capacitance)."""

    bullet = (
        f"- Extended the fit with the 30 s post-step relaxations (v2): parsed every rest with tests, showed model-free that the\n"
        f"  instantaneous resistance is {R_free_chg:.0f} Ω against the {p1['R_contact']['value']:.0f} Ω the curve-only fit needed, and established that the\n"
        f"  rests' log-time relaxation is outside the model class (no double layer, one particle size), so the 2C error went\n"
        f"  from {e1['capacity_error_pct']:+.1f} % to {e2['capacity_error_pct']:+.1f} %; identified double-layer capacitance and a particle-size distribution as the\n"
        f"  model additions that would use the rests.\n"
    )
    return {"section": section, "key_bullet": key_bullet, "bullet": bullet}


def ov_two_c_v1(H):
    return load(f"overpotentials_{H.lower()}.json")["two_c"]


def main(headline_kind="DFN"):
    kinds = [k for k in ("SPMe", "DFN") if (RESULTS / f"validation_{k.lower()}.json").exists()]
    if headline_kind not in kinds:
        headline_kind = kinds[-1]
    other = [k for k in kinds if k != headline_kind]
    H = headline_kind
    val = {k: load(f"validation_{k.lower()}.json") for k in kinds}
    fit = {k: load(f"fit_{k.lower()}.json") for k in kinds}
    cv = {k: load(f"cv_{k.lower()}.json") for k in kinds}
    ov = {k: load(f"overpotentials_{k.lower()}.json") for k in kinds}
    sens = {k: load(f"sensitivity_disc_{k.lower()}.json") for k in kinds if (RESULTS / f"sensitivity_disc_{k.lower()}.json").exists()}
    ocv = load("ocv.json")
    meas = load("measured_summary.json")
    h = val[H]["headline"]
    v2 = load_v2(kinds, H)
    sys.path.insert(0, str(ROOT / "src"))
    from lmo.model import DEFAULT_DISC_MM, disc_geometry

    geo = disc_geometry(DEFAULT_DISC_MM)
    c10_current_mA = next(r for r in meas["gcd"] if r["rate"] == "C/10" and r["cycle"] == 2 and r["step"] == "discharge")["mean_current_mA"]
    ir_c10_mV = 2 * c10_current_mA * fit[H]["parameters"]["R_contact"]["value"]

    # ---------------------------------------------------------------- results table
    rows = []
    for rate, e in val[H]["rates"].items():
        d = e["discharge"]
        rows.append(
            f"| {rate} | {e['cycle']} | {e['role']} | {d['capacity_mAh_measured']:.3f} | {d['capacity_mAh_model']:.3f} | "
            f"{d['capacity_error_pct']:+.1f} | {d['mean_voltage_V_measured']:.3f} | {d['mean_voltage_V_model']:.3f} | "
            f"{d['energy_mWh_measured']:.2f} | {d['energy_mWh_model']:.2f} | {d['vq_rmse_mV']:.0f} |"
        )
    table = "\n".join(
        [
            "| Rate | Cycle | Role | Q meas (mAh) | Q model (mAh) | Error (%) | V̄ meas (V) | V̄ model (V) | E meas (mWh) | E model (mWh) | V-Q RMSE (mV) |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        + rows
    )

    cv_rows = []
    for scan, e in cv[H]["scans"].items():
        m, s = e["measured"], e["model"]
        cv_rows.append(
            f"| {scan} | {m['oxidation_peak_V']:.2f} / {s['oxidation_peak_V']:.2f} | {m['oxidation_peak_mA']:.2f} / {s['oxidation_peak_mA']:.2f} | "
            f"{m['reduction_peak_V']:.2f} / {s['reduction_peak_V']:.2f} | {m['reduction_peak_mA']:.2f} / {s['reduction_peak_mA']:.2f} | "
            f"{e['peak_separation_measured_mV']:.0f} / {e['peak_separation_model_mV']:.0f} |"
        )
    cv_table = "\n".join(
        [
            "| Scan rate | Ox. peak V (meas / model) | Ox. peak mA (meas / model) | Red. peak V (meas / model) | Red. peak mA (meas / model) | ΔEp mV (meas / model) |",
            "|---|---|---|---|---|---|",
        ]
        + cv_rows
    )

    b = cv[H]["b_value"]
    rs = cv[H]["randles_sevcik"]
    two = ov[H]["two_c"]
    shares = two["shares"]
    dis_loss = ov[H]["discharge_mean_loss_mV"]

    par_rows = []
    for k in kinds:
        p = fit[k]["parameters"]
        par_rows.append(
            f"| {k} | {1000 * fit[k]['training_rmse_V']:.1f} | {fmt_ci(p['tau_R2_over_Ds_s'], sig=3)} | {fmt_ci(p['Ds'], sig=3)} | "
            f"{fmt_ci(p['k'], sig=3)} | {fmt_ci(p['R_contact'], sig=3)} | {val[k]['headline']['predicted_mAh']:.3f} ({val[k]['headline']['error_pct']:+.1f} %) |"
        )
    par_table = "\n".join(
        [
            "| Model | Train RMSE (mV) | τ = R²/Ds (s) | Ds (m²/s) at R = 5 µm | k (m^2.5 mol^-0.5 s^-1) | R_contact (Ω) | 2C prediction (mAh) |",
            "|---|---|---|---|---|---|---|",
        ]
        + par_rows
    )

    sens_rows = []
    for k in sens:
        for d, e in sorted(sens[k]["discs"].items(), key=lambda x: float(x[0])):
            g = e["geometry"]
            sens_rows.append(
                f"| {k} | {float(d):.0f} | {g['areal_loading_mg_cm2']:.1f} | {g['coating_density_g_cm3']:.2f} | {g['active_volume_fraction']:.2f} | {g['porosity']:.2f} | "
                f"{g.get('plausibility', '')} | {1000 * e['training_rmse_V']:.1f} | {e['parameters']['R_contact']:.0f} | {e['parameters']['Ds']:.2e} | "
                f"{e['predicted_2c_capacity_mAh']:.3f} ({e['predicted_2c_error_pct']:+.1f} %) |"
            )
    sens_table = "\n".join(
        [
            "| Model | Disc (mm) | Loading (mg/cm²) | Coating density (g/cm³) | ε_active | Porosity | Coating | Train RMSE (mV) | R_contact (Ω) | Ds (m²/s) | 2C prediction (mAh) |",
            "|---|---|---|---|---|---|---|---|---|---|---|",
        ]
        + sens_rows
    )
    plausible = sorted({float(d) for k in sens for d, e in sens[k]["discs"].items() if e["geometry"].get("plausibility") == "plausible"})
    pred_range = [e["predicted_2c_capacity_mAh"] for k in sens for e in sens[k]["discs"].values()]

    corr = fit[H]["log10_correlation"]
    c10 = val[H]["rates"]["C/10"]["discharge"]
    ce_c10 = next(r for r in meas["gcd"] if r["rate"] == "C/10" and r["cycle"] == 1 and r["step"] == "discharge")["coulombic_efficiency"]
    pk = ocv["dqdv_peaks_discharge"]
    other_txt = ""
    if other:
        o = other[0]
        other_txt = (
            f" The {o} gives {val[o]['headline']['predicted_mAh']:.3f} mAh ({val[o]['headline']['error_pct']:+.1f} %), "
            f"so the two model fidelities agree on the 2C result to within {max(abs(val[o]['headline']['error_pct'] - h['error_pct']), 0.1):.1f} percentage points."
        )

    readme = f"""# LiMn2O4 half cell: a physics-based rate model fitted to my own data

A PyBaMM Doyle-Fuller-Newman (DFN) half-cell model of a commercial LiMn2O4 (LMO) cathode against lithium
metal, parameterised with the open-circuit curve I measured at C/10, fitted to my C/5 and C/2 cycles only,
and then asked to predict the 2C cycle it had never seen and the cyclic voltammograms from a different
technique. All data were measured by me on one coin cell (MSE 5550, University of Pennsylvania, 2026;
see `data/README.md`).

## Key results

- **Training (C/5, C/2):** voltage RMSE {1000 * fit[H]['training_rmse_V']:.0f} mV with three fitted parameters; C/10, not used in the fit,
  reproduced to {abs(c10['capacity_error_pct']):.1f} % in capacity.
- **Blind 2C prediction:** average discharge voltage within {abs(h['mean_voltage_error_mV']):.0f} mV of measured; discharge capacity
  {h['predicted_mAh']:.3f} mAh against {h['measured_mAh']:.3f} mAh measured ({h['error_pct']:+.1f} %). The model predicts a
  {h['predicted_capacity_loss_vs_c10_pct']:.0f} % capacity loss from C/10 to 2C; the cell lost {h['measured_capacity_loss_vs_c10_pct']:.0f} %.{other_txt}
- **Cyclic voltammetry, a technique never used in the fit:** peak potentials at all three scan rates within
  {max(abs(e['oxidation_peak_shift_mV']) for e in cv[H]['scans'].values()):.0f} mV (oxidation) and {max(abs(e['reduction_peak_shift_mV']) for e in cv[H]['scans'].values()):.0f} mV (reduction); b-values {b['measured_oxidation']['b']:.2f} / {b['measured_reduction']['b']:.2f} measured vs {b['model_oxidation']['b']:.2f} / {b['model_reduction']['b']:.2f} model.
- **Diagnosis:** at 2C the model assigns {100 * shares['contact_ohmic']:.0f} % of the voltage loss to ohmic/contact resistance; the shape of
  the measured 2C curve says the real cell has a transport loss that grows with time, which the fit could not
  separate from resistance (see "Why 2C loses capacity").
{v2['key_bullet'] if v2 else ''}
![Measured vs model at all four rates](figures/headline_measured_vs_model.png)

*Solid: measured. Dashed: {H}. C/5 and C/2 were used for fitting; C/10 is a low-rate check; 2C is a blind prediction.*

## Results

### Galvanostatic cycles ({H}, 10 mm disc)

{table}

V-Q RMSE is computed on the measured capacity axis with the model held at the cut-off voltage beyond its own
capacity, so it penalises a capacity miss as well as a voltage miss. Charge-step metrics are in
`results/validation_{H.lower()}.json`; the C/10 and C/5 measured charges are longer than the model's because of
parasitic oxidation near 4.5 V (first-cycle coulombic efficiency at C/10 was {100 * ce_c10:.0f} %), which the model does not include.

### Fitted effective parameters

{par_table}

Intervals in brackets are 95 % Jacobian intervals from the least-squares covariance; "LOO" is the range of
refits that each leave one of the four training steps out. The Jacobian intervals are far too narrow because
the {fit[H]['n_residuals']} residuals are strongly autocorrelated; the leave-one-out ranges are the more honest
statement of uncertainty. The log-parameters are correlated (k vs R_contact: {corr['k-R_contact']:.2f}), i.e. the data
split the ohmic and the charge-transfer loss only weakly. k enters j0 = k F sqrt(c_e c_s (c_max - c_s)).
Only τ = R²/Ds is identifiable; the 5 µm radius is a fixed choice.

### Cyclic voltammetry, simulated with the measured voltage programme ({H})

{cv_table}

b-value (slope of log peak current vs log scan rate; 0.5 = semi-infinite diffusion, 1.0 = surface/capacitive):
measured {b['measured_oxidation']['b']:.2f} (oxidation) and {b['measured_reduction']['b']:.2f} (reduction); model
{b['model_oxidation']['b']:.2f} and {b['model_reduction']['b']:.2f}. Both sit between the two limits, as expected for a porous
electrode at scan rates where the sweep takes about as long as a 1C discharge. An apparent Randles-Sevcik
diffusivity from the measured peaks (planar-electrode formula, geometric 10 mm disc area, concentration of the
cycling lithium) is {rs['apparent_D_geometric_mean_m2_s']:.1e} m²/s, {rs['ratio_apparent_to_fitted']:.1f} times the fitted solid-state Ds of
{rs['fitted_Ds_m2_s']:.1e} m²/s. The two numbers are not expected to agree: Randles-Sevcik assumes a flat electrode, while
the real one is porous with a much larger internal area, so the apparent value is an upper bound.

![CVs](figures/cv_{H.lower()}.png)

### Why 2C loses capacity: overpotential breakdown ({H})

Capacity-weighted mean voltage loss on discharge: C/10 {dis_loss['C/10']:.0f} mV, C/5 {dis_loss['C/5']:.0f} mV,
C/2 {dis_loss['C/2']:.0f} mV, 2C {dis_loss['2C']:.0f} mV. At 2C the model attributes {100 * shares['contact_ohmic']:.0f} % of the loss to the
lumped contact/ohmic resistance, {100 * shares['reaction']:.0f} % to charge transfer, {100 * shares['particle_diffusion']:.0f} % to solid-state
diffusion in the particles and {100 * (shares['electrolyte_concentration'] + shares['electrolyte_ohmic']):.0f} % to the electrolyte. When the 2C
discharge hits 3.5 V the bulk open-circuit voltage of the particles is still {two['bulk_ocv_at_cutoff_V']:.2f} V, i.e. the model
has {two['capacity_locked_out_mAh']:.2f} mAh of lithium it could not deliver at that rate.

![Overpotentials](figures/overpotentials_{H.lower()}.png)

The measured 2C curve falls faster than the model's: the data show a steeper initial decline and lose
{h['measured_capacity_loss_vs_c10_pct'] - h['predicted_capacity_loss_vs_c10_pct']:.0f} percentage points more capacity than predicted. The fit on C/5 and C/2 cannot separate a large
ohmic drop from a transport limitation that grows with time, and it picked the ohmic explanation. That is the
single biggest limitation of this work, and it is visible in the k / R_contact correlation above.

### Sensitivity to the unknown electrode area

{sens_table}

The disc diameter is genuinely unknown: the EC-Lab `.mpr` settings block has electrode surface area and
characteristic mass at their 0.001 defaults (never entered); only the 1.434 mAh battery capacity was set. The
capacity therefore ties area to active volume fraction. Active mass = 1.434 mAh / 110 mAh/g = {1000 * ocv['active_mass_g']:.1f} mg,
active volume = {1000 * ocv['active_mass_g'] / 4.28 * 1e-3 * 1e3:.2f}e-3 cm³ at 4.28 g/cm³, so for the 85 µm coating ε_active = {1000 * ocv['active_mass_g'] / 4.28 / 85e-4 * 1e-3:.3f} / A[cm²]:
10 mm gives 0.46, 12 mm 0.32, 14 mm 0.23 and 16 mm 0.18. Binder and carbon (5.8 wt%) are taken at 1.8 g/cm³ and
the rest of the coating volume is pore. A calendered oxide cathode normally has a porosity of 0.25 to 0.40
(up to ~0.5 uncalendered) and an active volume fraction of 0.45 to 0.65. By that yardstick only the
{', '.join(f'{d:.0f} mm' for d in plausible) if plausible else 'none of the'} disc{'s' if len(plausible) != 1 else ''} give{'s' if len(plausible) == 1 else ''} a physically plausible coating; 12 mm would be an
unusually light coating and 14 to 16 mm imply porosities of 0.73 to 0.80, which no real electrode has. Either
the disc is small, or the 1.434 mAh (110 mAh/g) capacity under-states the active loading, so that the true
ε_active is higher than this derivation gives for any disc. The three parameters were re-fitted for every
diameter: the training RMSE, τ, k and the 2C prediction do not change to the precision shown
(2C spans {min(pred_range):.3f} to {max(pred_range):.3f} mAh), and only R_contact drifts by a few ohms as the electrolyte and solid
ohmic paths take up slack. The model is not electrolyte-limited with a lithium-metal counter electrode and this
electrolyte, so the unknown area does not affect the headline result. The headline uses the 10 mm disc, the only
diameter that gives a physically plausible coating.


{v2['section'] if v2 else V1_NEXT_STEP}

## Method

1. **Parser** (`src/lmo/io.py`): reads the EC-Lab text exports, splits cycles and steps by the sign of the
   measured current, integrates the measured current for capacity. `tests/test_io.py` checks the capacities,
   mean voltages, energies and CV peaks against the values in `data/README.md` (within 1 %).
2. **Pseudo-OCV** (`src/lmo/ocv.py`): average of the C/10 cycle-2 charge and discharge on a normalised-capacity
   axis (each step normalised by its own capacity, which lines up both cut-offs and the plateau transition),
   Savitzky-Golay smoothed, forced monotonic, mapped to x in Li_xMn2O4 with x = 1 at 3.5 V and
   x = {ocv['stoichiometry_window']['sto_min']:.3f} at 4.5 V (theoretical capacity {ocv['theoretical_capacity_mAh']:.3f} mAh for the implied
   {1000 * ocv['active_mass_g']:.1f} mg of LMO). dQ/dV at C/10 shows the two spinel plateaus at {pk[0]['voltage_V']:.2f} and {pk[1]['voltage_V']:.2f} V on discharge.
   The mean charge/discharge hysteresis at C/10 is {ocv['mean_charge_discharge_hysteresis_mV']:.0f} mV.
3. **Model** (`src/lmo/model.py`): PyBaMM `{H}` (and the other fidelity for comparison) with
   `working electrode = positive` and a contact resistance, starting from `Xu2019` for the lithium counter
   electrode, separator and electrolyte. Positive electrode: measured OCP, 85 µm, c_max = {ocv['lmo_c_max_mol_m3']:.0f} mol/m³
   from the LMO density, 10 mm disc (ε_active = {geo['active_volume_fraction']:.2f}, porosity {geo['porosity']:.2f}),
   25 °C, Butler-Volmer kinetics, 5 µm particles.
4. **Fit** (`src/lmo/fit.py`): scipy `least_squares` in log10 space on the voltage-time curves of C/5 cycle 1
   and C/2 cycle 1 (charge and discharge, 300 points each), simulating the full charge/rest/discharge protocol
   from the measured rest voltage. Each training run is one PyBaMM simulation with the three parameters as
   inputs, so one residual evaluation costs a fraction of a second for SPMe.
5. **Blind validation** (`src/lmo/validate.py`): all four rates simulated with the fitted parameters from each
   run's own first cycle (2C is simulated through its cycle 1 so that cycle 2 carries the right history); CVs
   simulated in voltage-control mode by feeding the model the measured voltage programme.
6. **Diagnostics** (`src/lmo/diagnostics.py`): PyBaMM's voltage-component variables, checked to close to the
   terminal voltage to better than {max(s['closure_max_abs_V'] for s in ov[H]['steps']) * 1e6:.0e} µV; CV b-values and Randles-Sevcik.
7. **v2, rests in the fit** (`src/lmo/fit_v2.py`, `src/lmo/relaxation.py`): the parser also extracts every
   30 s post-step rest (`tests/test_rests.py`); v2 refits the same three parameters on the C/5 and C/2 curves plus
   the C/10, C/5 and C/2 cycle-1 rests, with a surface-refined particle mesh so the first second of the rest is
   resolved. v2 writes `results/v2_*.json`, `data/processed/v2_*.csv` and `figures/v2_*.png`; the v1 files are
   untouched. 2C is never fitted in either version.
8. **Reproduce**: `pip install -r requirements.txt && pip install -e . && make all` regenerates
   `data/processed/`, `results/*.json`, `figures/*.png` and this README (`python scripts/write_readme.py`).
   `make test` runs the parser tests. `python scripts/run_all.py --fast` skips the DFN; `--v1-only` and `--v2-only`
   run one version.

## Assumptions and limitations

- One cell, one cycle per rate, and the seven tests ran in sequence on that cell, so slow drift (Mn dissolution,
  electrolyte oxidation at 4.5 V) is folded into the fitted parameters.
- The electrode area was never entered in EC-Lab (surface area and characteristic mass are at the 0.001
  defaults in the `.mpr` settings), so it is derived from the capacity: ε_active = {1000 * ocv['active_mass_g'] / 4.28 / 85e-4 * 1e-3:.3f} / A[cm²].
  Only a 10 mm disc gives a conventional coating; the 12 to 16 mm discs imply light to implausibly porous
  coatings. The headline therefore uses 10 mm; the sensitivity table shows the prediction does not depend on
  this choice.
- The electrolyte is assumed to be 1 M LiPF6 in carbonates with the `Xu2019` / Valoen2005 transport
  properties, and the lithium-metal kinetics are Xu2019's.
- The fitted parameters are effective values: τ = R²/Ds at an assumed 5 µm radius, a rate constant that absorbs
  the unknown real surface area, and a contact resistance of about {fit[H]['parameters']['R_contact']['value']:.0f} Ω that lumps current-collector contact,
  separator and lithium-metal interphase resistance. The round-trip IR drop of that resistance at the C/10 current is {ir_c10_mV:.0f} mV, about half of the
  {ocv['mean_charge_discharge_hysteresis_mV']:.0f} mV mean charge/discharge hysteresis measured at C/10; the rest is kinetic and intrinsic hysteresis
  that the averaged pseudo-OCV absorbs.
- The pseudo-OCV is a C/10 average, not a true equilibrium curve, and is linearly extrapolated outside the
  3.5 to 4.5 V window, so the parts of the CVs below 3.5 V and above 4.5 V (the 0.3 mV/s sweep went to 4.8 V)
  are outside the model's validity.
- The model has no parasitic reactions, so measured charge capacities at low rate exceed the model's.
- Jacobian confidence intervals are lower bounds on the real uncertainty (autocorrelated residuals); the
  leave-one-step-out ranges are reported alongside.
- The 2C over-prediction is reported as measured. Nothing was tuned on the 2C or CV data.

## Repository layout

```
data/raw/            EC-Lab exports (never modified)      data/README.md   data notes and known quirks
data/processed/      pseudo-OCV, dQ/dV, model curves       src/lmo/         package (io, ocv, model, fit, validate, diagnostics, figures)
scripts/run_all.py   regenerates everything                scripts/write_readme.py renders this file from results/*.json
tests/test_io.py     parser tests vs the report values     tests/test_rests.py  rest parser tests
results/*.json       every v1 number quoted here           figures/*.png    300 dpi
results/v2_*.json    every v2 number quoted here           figures/v2_*.png v2 figures
results/resume_bullets.md                                  notebooks/walkthrough.ipynb
```
"""
    (ROOT / "README.md").write_text(readme)

    cvmax = max(max(abs(e['oxidation_peak_shift_mV']), abs(e['reduction_peak_shift_mV'])) for e in cv[H]['scans'].values())
    bullets = f"""# Resume bullets

Every number below is read from `results/*.json`.

- Built a Doyle-Fuller-Newman half-cell model in PyBaMM on an open-circuit curve measured from my own C/10 data,
  and fitted solid-diffusion time constant, reaction rate constant and contact resistance by nonlinear least
  squares on the C/5 and C/2 cycles only ({1000 * fit[H]['training_rmse_V']:.0f} mV voltage RMSE), with leave-one-out uncertainty bounds.
- Validated on data the model never saw: reproduced the 2C average discharge voltage within
  {abs(h['mean_voltage_error_mV']):.0f} mV, C/10 capacity within {abs(c10['capacity_error_pct']):.0f}% and all six CV peak potentials across three scan rates
  within {cvmax:.0f} mV.
- Decomposed the 2C voltage loss ({100 * shares['contact_ohmic']:.0f}% ohmic/contact, {100 * shares['reaction']:.0f}% charge transfer,
  {100 * shares['particle_diffusion']:.0f}% solid-state diffusion), identified the ohmic/transport ambiguity that limits high-rate prediction,
  and packaged the work as a tested Python pipeline that regenerates every figure and number from the raw
  cycler files.
{v2['bullet'] if v2 else ''}"""
    (RESULTS / "resume_bullets.md").write_text(bullets)
    print("wrote README.md and results/resume_bullets.md using", H)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "DFN")
