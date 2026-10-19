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
