# Resume bullets

- Built a physics-based PyBaMM DFN/SPMe half-cell model of a LiMn2O4 cathode from my own galvanostatic and
  CV measurements, parameterised with a measured pseudo-OCV and three effective parameters (solid diffusion
  time constant, exchange-current rate constant, contact resistance) fitted by nonlinear least squares with
  uncertainty from Jacobian covariance and leave-one-out refits.
- Fitted on C/5 and C/2 only, the model predicted the unseen 2C discharge capacity to within
  20 % (0.99 vs 0.82 mAh) and the 2C average voltage to within
  18 mV, and reproduced the oxidation/reduction peak potentials of three cyclic voltammograms
  at 0.1 to 0.3 mV/s to within 80 mV without any CV data in the fit.
- Used the model's overpotential decomposition to attribute the 42 % capacity loss between
  C/10 and 2C mainly to ohmic/contact resistance (77 % of the 2C voltage loss) and charge-transfer
  kinetics (15 %), with CV b-values of 0.60 to 0.69 confirming mixed diffusion/surface control;
  packaged the pipeline as a tested Python package that regenerates every figure and number from raw data.
