"""Regenerate everything from data/raw: processed tables, results/*.json and figures/*.png.

Usage: python scripts/run_all.py [--fast]
  --fast  skip the DFN fit and the DFN disc sensitivity (SPMe only, a few minutes)
"""

from __future__ import annotations

import json
import sys
import time

from lmo import RESULTS_DIR, diagnostics, figures, fit, ocv, validate
from lmo.io import cv_peaks, gcd_summary
from lmo import CVS


def main(fast: bool = False):
    t0 = time.time()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("1. parser summary")
    summary = gcd_summary()
    summary.to_csv(RESULTS_DIR.parent / "data" / "processed" / "gcd_summary.csv", index=False)
    with open(RESULTS_DIR / "measured_summary.json", "w") as f:
        json.dump(
            {
                "gcd": summary.to_dict(orient="records"),
                "cv_peaks": {scan: cv_peaks(scan, 2) for scan in CVS},
            },
            f,
            indent=2,
        )

    print("2. pseudo-OCV and dQ/dV")
    ocv.run()

    kinds = ["SPMe"] if fast else ["SPMe", "DFN"]
    for kind in kinds:
        print(f"3-4. {kind} fit on C/5 and C/2 cycle 1")
        fit.run(kind)
        print(f"5. {kind} blind validation (all rates, CVs)")
        validate.predict_all(kind)
        validate.validate_cv(kind)
        print(f"6. {kind} overpotential breakdown")
        diagnostics.run(kind)
        print(f"3b. {kind} disc-diameter sensitivity")
        validate.disc_sensitivity(kind)

    print("7. figures")
    headline = "DFN" if not fast else "SPMe"
    for p in figures.make_all(headline, kinds):
        print("   ", p)

    write_summary(kinds)
    print(f"done in {(time.time() - t0) / 60:.1f} min")


def write_summary(kinds):
    """results/summary.json: the numbers the README quotes, gathered in one place."""
    out = {}
    for kind in kinds:
        val = json.load(open(RESULTS_DIR / f"validation_{kind.lower()}.json"))
        fitj = json.load(open(RESULTS_DIR / f"fit_{kind.lower()}.json"))
        cv = json.load(open(RESULTS_DIR / f"cv_{kind.lower()}.json"))
        ov = json.load(open(RESULTS_DIR / f"overpotentials_{kind.lower()}.json"))
        out[kind] = {
            "headline": val["headline"],
            "training_rmse_mV": 1000 * fitj["training_rmse_V"],
            "parameters": fitj["parameters"],
            "rates": {r: {"role": e["role"], "discharge": e["discharge"]} for r, e in val["rates"].items()},
            "cv": {s: {k: v for k, v in e.items() if k != "scan_rate_V_s"} for s, e in cv["scans"].items()},
            "b_value": cv["b_value"],
            "randles_sevcik": cv["randles_sevcik"],
            "two_c_overpotentials": ov["two_c"],
        }
    with open(RESULTS_DIR / "summary.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(fast="--fast" in sys.argv)
