# Resume bullets

Every number below is read from `results/*.json`.

- Built a Doyle-Fuller-Newman half-cell model in PyBaMM on an open-circuit curve measured from my own C/10 data,
  and fitted solid-diffusion time constant, reaction rate constant and contact resistance by nonlinear least
  squares on the C/5 and C/2 cycles only (26 mV voltage RMSE), with leave-one-out uncertainty bounds.
- Validated on data the model never saw: reproduced the 2C average discharge voltage within
  17 mV, C/10 capacity within 1% and all six CV peak potentials across three scan rates
  within 80 mV.
- Decomposed the 2C voltage loss (74% ohmic/contact, 17% charge transfer,
  6% solid-state diffusion), identified the ohmic/transport ambiguity that limits high-rate prediction,
  and packaged the work as a tested Python pipeline that regenerates every figure and number from the raw
  cycler files.
- Extended the fit with the 30 s post-step relaxations (v2): parsed every rest with tests, showed model-free that the
  instantaneous resistance is 72 Ω against the 97 Ω the curve-only fit needed, and established that the
  rests' log-time relaxation is outside the model class (no double layer, one particle size), so the 2C error went
  from +20.1 % to -16.3 %; identified double-layer capacitance and a particle-size distribution as the
  model additions that would use the rests.
- Diagnosed the model class from the rests and changed it (v2.1): adding double-layer capacitance let the same data
  fit the rests to 15 mV with a contact resistance of 89 Ω (leave-one-out 74 to 104 Ω, still above the
  model-free 72 Ω jump); the blind 2C capacity error went from +20.1 % (v1) to +27.1 % (v2.1), reported without any tuning on 2C.
