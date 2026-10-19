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
    s = f"{g(v)} ({g(lo)} to {g(hi)})"
    if "leave_one_step_out_min" in p:
        s += f"; LOO {g(p['leave_one_step_out_min'] * scale)} to {g(p['leave_one_step_out_max'] * scale)}"
    return s + (f" {unit}" if unit else "")


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
    geo = sens[H]["discs"]["14.0"]["geometry"] if H in sens else None

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
                f"| {k} | {float(d):.0f} | {g['areal_loading_mg_cm2']:.1f} | {g['active_volume_fraction']:.2f} | {g['porosity']:.2f} | "
                f"{e['parameters']['R_contact']:.0f} | {e['parameters']['Ds']:.2e} | {e['predicted_2c_capacity_mAh']:.3f} ({e['predicted_2c_error_pct']:+.1f} %) |"
            )
    sens_table = "\n".join(
        [
            "| Model | Disc (mm) | Loading (mg/cm²) | ε_active | Porosity | R_contact (Ω) | Ds (m²/s) | 2C prediction (mAh) |",
            "|---|---|---|---|---|---|---|---|",
        ]
        + sens_rows
    )

    corr = fit[H]["log10_correlation"]
    c10 = val[H]["rates"]["C/10"]["discharge"]
    ce_c10 = next(r for r in meas["gcd"] if r["rate"] == "C/10" and r["cycle"] == 1 and r["step"] == "discharge")["coulombic_efficiency"]
    pk = ocv["dqdv_peaks_discharge"]
    other_txt = ""
    if other:
        o = other[0]
        other_txt = (
            f" The {o} gives {val[o]['headline']['predicted_mAh']:.3f} mAh ({val[o]['headline']['error_pct']:+.1f} %), "
            f"so the two model fidelities agree on the 2C result to within {abs(val[o]['headline']['error_pct'] - h['error_pct']):.1f} percentage points."
        )

    readme = f"""# LiMn2O4 half cell: a physics-based rate model fitted to my own data

A PyBaMM Doyle-Fuller-Newman (DFN) half-cell model of a commercial LiMn2O4 (LMO) cathode against lithium
metal, parameterised with the open-circuit curve I measured at C/10, fitted to my C/5 and C/2 cycles only,
and then asked to predict the 2C cycle it had never seen and the cyclic voltammograms from a different
technique. All data were measured by me on one coin cell (MSE 5550, University of Pennsylvania, 2026;
see `data/README.md`).

**Headline: the model fitted on C/5 and C/2 predicted the 2C discharge capacity as {h['predicted_mAh']:.3f} mAh against
{h['measured_mAh']:.3f} mAh measured, an over-prediction of {h['error_pct']:.1f} %.** It got the average 2C discharge voltage to
{abs(h['mean_voltage_error_mV']):.0f} mV and reproduced {h['predicted_capacity_loss_vs_c10_pct']:.0f} % of the measured
{h['measured_capacity_loss_vs_c10_pct']:.0f} % capacity loss between C/10 and 2C.{other_txt}
The same parameters put the three CV peak positions within {max(abs(e['oxidation_peak_shift_mV']) for e in cv[H]['scans'].values()):.0f} mV (oxidation)
and {max(abs(e['reduction_peak_shift_mV']) for e in cv[H]['scans'].values()):.0f} mV (reduction) of the measured ones.

![Measured vs model at all four rates](figures/headline_measured_vs_model.png)

*Solid: measured. Dashed: {H}. C/5 and C/2 were used for fitting; C/10 is a low-rate check; 2C is a blind prediction.*

## Results

### Galvanostatic cycles ({H}, 14 mm disc assumption)

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
diffusivity from the measured peaks (planar-electrode formula, geometric 14 mm disc area, concentration of the
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

The disc diameter is not known. For each diameter the active volume fraction and porosity were re-derived from
the capacity (active mass from the EC-Lab capacity and the 110 mAh/g rating, LMO density 4.28 g/cm³, binder and
carbon at 1.8 g/cm³) and the three parameters re-fitted. The implied porosities are high, which means either the
coating is unusually light or the nominal capacity under-states the loading; the 2C prediction barely moves,
because with a lithium-metal counter electrode and this electrolyte the half cell is not electrolyte-limited
in the model.

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
   from the LMO density, 14 mm disc (ε_active = {geo['active_volume_fraction'] if geo else float('nan'):.2f}, porosity {geo['porosity'] if geo else float('nan'):.2f}),
   25 °C, Butler-Volmer kinetics, 5 µm particles.
4. **Fit** (`src/lmo/fit.py`): scipy `least_squares` in log10 space on the voltage-time curves of C/5 cycle 1
   and C/2 cycle 1 (charge and discharge, 300 points each), simulating the full charge/rest/discharge protocol
   from the measured rest voltage. Each training run is one PyBaMM simulation with the three parameters as
   inputs, so one residual evaluation costs a fraction of a second for SPMe.
5. **Blind validation** (`src/lmo/validate.py`): all four rates simulated with the fitted parameters from each
   run's own first cycle (2C is simulated through its cycle 1 so that cycle 2 carries the right history); CVs
   simulated in voltage-control mode by feeding the model the measured voltage programme.
6. **Diagnostics** (`src/lmo/diagnostics.py`): PyBaMM's voltage-component variables, checked to close to the
   terminal voltage within {max(s['closure_max_abs_V'] for s in ov[H]['steps']) * 1e3:.1e} mV; CV b-values and Randles-Sevcik.
7. **Reproduce**: `pip install -r requirements.txt && pip install -e . && make all` regenerates
   `data/processed/`, `results/*.json`, `figures/*.png` and this README (`python scripts/write_readme.py`).
   `make test` runs the parser tests. `python scripts/run_all.py --fast` skips the DFN.

## Assumptions and limitations

- One cell, one cycle per rate, and the seven tests ran in sequence on that cell, so slow drift (Mn dissolution,
  electrolyte oxidation at 4.5 V) is folded into the fitted parameters.
- The electrode area is unknown; see the sensitivity table. The electrolyte is assumed to be 1 M LiPF6 in
  carbonates with the `Xu2019` / Valoen2005 transport properties, and the lithium-metal kinetics are Xu2019's.
- The fitted parameters are effective values: τ = R²/Ds at an assumed 5 µm radius, a rate constant that absorbs
  the unknown real surface area, and a contact resistance of about {fit[H]['parameters']['R_contact']['value']:.0f} Ω that lumps current-collector contact,
  separator and lithium-metal interphase resistance. The C/10 hysteresis of {ocv['mean_charge_discharge_hysteresis_mV']:.0f} mV is consistent with a
   resistance of that size at 0.14 mA.
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
tests/test_io.py     parser tests vs the report values     notebooks/walkthrough.ipynb
results/*.json       every number quoted here              figures/*.png    300 dpi
results/resume_bullets.md
```
"""
    (ROOT / "README.md").write_text(readme)

    bullets = f"""# Resume bullets

- Built a physics-based PyBaMM DFN/SPMe half-cell model of a LiMn2O4 cathode from my own galvanostatic and
  CV measurements, parameterised with a measured pseudo-OCV and three effective parameters (solid diffusion
  time constant, exchange-current rate constant, contact resistance) fitted by nonlinear least squares with
  uncertainty from Jacobian covariance and leave-one-out refits.
- Fitted on C/5 and C/2 only, the model predicted the unseen 2C discharge capacity to within
  {h['abs_error_pct']:.0f} % ({h['predicted_mAh']:.2f} vs {h['measured_mAh']:.2f} mAh) and the 2C average voltage to within
  {abs(h['mean_voltage_error_mV']):.0f} mV, and reproduced the oxidation/reduction peak potentials of three cyclic voltammograms
  at 0.1 to 0.3 mV/s to within {max(max(abs(e['oxidation_peak_shift_mV']), abs(e['reduction_peak_shift_mV'])) for e in cv[H]['scans'].values()):.0f} mV without any CV data in the fit.
- Used the model's overpotential decomposition to attribute the {h['measured_capacity_loss_vs_c10_pct']:.0f} % capacity loss between
  C/10 and 2C mainly to ohmic/contact resistance ({100 * shares['contact_ohmic']:.0f} % of the 2C voltage loss) and charge-transfer
  kinetics ({100 * shares['reaction']:.0f} %), with CV b-values of {b['measured_oxidation']['b']:.2f} to {b['measured_reduction']['b']:.2f} confirming mixed diffusion/surface control;
  packaged the pipeline as a tested Python package that regenerates every figure and number from raw data.
"""
    (RESULTS / "resume_bullets.md").write_text(bullets)
    print("wrote README.md and results/resume_bullets.md using", H)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "DFN")
