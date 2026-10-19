"""LMO half-cell rate model: parser, pseudo-OCV, PyBaMM model, fit, validation, diagnostics."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROCESSED_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

# Nominal capacity entered in EC-Lab, used to define the C-rates (mAh).
NOMINAL_CAPACITY_MAH = 1.434

# C-rate labels, set currents (mA) and file stems. Colours match the report.
RATES = {
    "C/10": {"file": "GCD_C-10", "current_mA": 0.143, "c_rate": 0.1, "color": "#1f77b4"},
    "C/5": {"file": "GCD_C-5", "current_mA": 0.287, "c_rate": 0.2, "color": "#d62728"},
    "C/2": {"file": "GCD_C-2", "current_mA": 0.717, "c_rate": 0.5, "color": "#2ca02c"},
    "2C": {"file": "GCD_2C", "current_mA": 2.868, "c_rate": 2.0, "color": "#9467bd"},
}

# Scan rate label -> file stem and scan rate in V/s.
CVS = {
    "0.1 mV/s": {"file": "CV_0p1mVs", "scan_rate_V_s": 1e-4, "v_min": 3.0, "v_max": 4.5},
    "0.2 mV/s": {"file": "CV_0p2mVs", "scan_rate_V_s": 2e-4, "v_min": 3.0, "v_max": 4.5},
    "0.3 mV/s": {"file": "CV_0p3mVs", "scan_rate_V_s": 3e-4, "v_min": 3.0, "v_max": 4.8},
}

# Which measured cycle the report treats as representative for each rate.
REPRESENTATIVE_CYCLE = {"C/10": 1, "C/5": 1, "C/2": 1, "2C": 2}
