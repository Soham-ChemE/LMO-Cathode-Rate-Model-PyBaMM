# Raw data

Measured by Soham Kavathekar, MSE 5550 (Electrochemical Engineering of Materials), University of Pennsylvania,
Sep 28 to Oct 3, 2026. Full write-up: `MSE5550_assignment1_report.pdf`.

**Cell:** commercial LiMn2O4 electrode sheet (MTI Corp.): 15 µm Al foil, 85 µm single-side coating, 94.2 wt% active
material, 110 mAh/g rated. Potentials vs Li/Li+ (lithium-metal counter electrode, i.e. a half cell).
Capacity used for C-rates: 1.434 mAh (EC-Lab setting). Bio-Logic VMP-300, EC-Lab V11.33, channel 10.
All seven tests ran on the **same cell**, in this order:

| # | Test | Set current / scan rate | Window (V) | File |
|---|------|-------------------------|------------|------|
| 1 | Galvanostatic C/10 | 0.143 mA | 3.5-4.5 | `GCD_C-10.txt` (+ `.mpr`) |
| 2 | Galvanostatic C/5 | 0.287 mA | 3.5-4.5 | `GCD_C-5.txt` (+ `.mpr`) |
| 3 | Galvanostatic C/2 | 0.717 mA | 3.5-4.5 | `GCD_C-2.txt` |
| 4 | CV | 0.2 mV/s | 3.0-4.5 | `CV_0p2mVs.txt` |
| 5 | CV | 0.3 mV/s | 3.0-4.8 | `CV_0p3mVs.txt` |
| 6 | Galvanostatic 2C | 2.868 mA | 3.5-4.5 | `GCD_2C.txt` |
| 7 | CV | 0.1 mV/s | 3.0-4.5 | `CV_0p1mVs.txt` |

Galvanostatic protocol: charge to 4.5 V, discharge to 3.5 V, 30 s rest after each step, 2 cycles, record every
10 mV or 5 s. CV: 2 cycles. Text exports are tab-separated, latin-1, CRLF:
GCD columns `time/s, Ewe/V, <I>/mA, Capacity/mA.h, cycle number` (cycles 0 and 1);
CV columns `time/s, Ewe/V, <I>/mA, cycle number` (cycles 1 and 2).

## Known data issues (from the report; handle them, do not hide them)
- **C/10 current:** run on the 10 mA range, so current reads +0.147 mA charge / -0.140 mA discharge vs 0.143 set.
  Integrate the measured current; first-cycle charge capacity and coulombic efficiency (91%) are approximate.
- **C/5 contact disturbance:** at 19.2 h (second discharge) the voltage jumps ~1 s and settles ~22 mV lower for
  the rest of that discharge, current unchanged. Use C/5 cycle 1 for fitting.
- **2C:** followed a CV that ended at 3.0 V, so the first charge removed extra lithium (1.01 vs 0.85 mAh).
  Cycle 2 is the representative 2C cycle.
- **CV current-range switching (Auto range):** steps near ±1 mA are instrument artefacts. At 0.1 mV/s the anodic
  current sits near 1.03 mA from ~4.20-4.25 V then jumps to 1.6 mA; the real oxidation peak is at 4.36 V.
- **CV 0.1 mV/s** started from rest (3.96 V); it still contains one full 3.0 -> 4.5 -> 3.0 V cycle (cycle 2).
- The 0.2 mV/s file starts with a -8.8 mA inrush at 3.0 V (first seconds); ignore it.

## Summary values from the report (use as checks for the parser)
Discharge capacity (mAh): C/10 1.430, C/5 1.422, C/2 1.361, 2C 0.823 (cycle 2).
Average discharge voltage (V): 4.04, 4.02, 3.96, 3.70. Discharge energy (mWh): 5.77, 5.71, 5.39, 3.04.
CV peaks: oxidation 4.36 / 4.42 / 4.49 V (1.62 / 2.42 / 3.12 mA), reduction 3.81 / 3.72 / 3.64 V
(-1.28 / -2.06 / -2.73 mA) at 0.1 / 0.2 / 0.3 mV/s.
