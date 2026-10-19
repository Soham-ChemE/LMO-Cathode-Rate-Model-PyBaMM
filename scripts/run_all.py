"""Regenerate everything from data/raw: processed tables, results/*.json and figures/*.png.

Usage: python scripts/run_all.py [--fast] [--v1-only | --v2-only]
  --fast     skip the DFN (SPMe only, a few minutes)
  --v1-only  the v1 pipeline alone (fit on C/5 and C/2 curves; results/*.json, figures/*.png)
  --v2-only  the v2 pipeline alone (fit on curves plus rests; results/v2_*.json, figures/v2_*.png),
             which needs the v1 result files for the comparison
"""

from __future__ import annotations

import json
import sys
import time

from lmo import RESULTS_DIR, diagnostics, figures, fit, fit_v2, fit_v21, ocv, relaxation, validate
from lmo.io import cv_peaks, gcd_summary
from lmo import CVS


def main(fast: bool = False, v1: bool = True, v2: bool = True):
    t0 = time.time()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    kinds = ["SPMe"] if fast else ["SPMe", "DFN"]
    headline = "DFN" if not fast else "SPMe"
    if v1:
        run_v1(kinds, headline)
    if v2:
        run_v2(kinds, headline)
    write_summary(kinds)
    import importlib.util

    spec = importlib.util.spec_from_file_location("write_readme", RESULTS_DIR.parent / "scripts" / "write_readme.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main(headline)
    print(f"done in {(time.time() - t0) / 60:.1f} min")


def run_v1(kinds, headline):
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
    for p in figures.make_all(headline, kinds):
        print("   ", p)


def run_v2(kinds, headline):
    """v2: rests as data, fit on curves plus rests (never 2C), blind validation, post-hoc 2C rests."""
    print("v2-1. post-step rests: model-free resistance table")
    relaxation.measured_rests()
    for kind in kinds:
        print(f"v2-2. {kind} fit on C/5 and C/2 curves plus C/10, C/5, C/2 rests (with LOO; weighting sensitivity on the SPMe)")
        # The DFN costs ~8 min per fit: it runs the headline fit with LOO and the jump-only fit
        # without LOO; the weighting sensitivity and the jump-only LOO are SPMe only.
        fit_v2.run(kind, weighting_sensitivity=(kind == "SPMe"), jump_loo=(kind == "SPMe"))
        print(f"v2-3. {kind} blind validation (all rates, CVs)")
        validate.predict_all(kind, prefix="v2_", mesh="v2")
        validate.validate_cv(kind, prefix="v2_", mesh="v2")
        print(f"v2-4. {kind} overpotential breakdown")
        diagnostics.run(kind, prefix="v2_")
        print(f"v2-5. {kind} post-hoc 2C rests (not fitted)")
        relaxation.posthoc_2c_rests(kind)
    for kind in kinds:
        # The SPMe fits all three variants (the variant table); the DFN, whose size-distribution
        # variants cost an order of magnitude more, runs the headline variant from the SPMe optimum.
        variants = tuple(fit_v21.VARIANTS) if kind == "SPMe" else (fit_v21.HEADLINE_VARIANT,)
        print(f"v2.1-1. {kind}: variants {variants}, headline '{fit_v21.HEADLINE_VARIANT}' with LOO")
        fit_v21.run(kind, variants=variants)
        kw = fit_v21.sim_kwargs(fit_v21.HEADLINE_VARIANT)
        print(f"v2.1-2. {kind} blind validation, CVs, overpotentials, post-hoc 2C rests")
        validate.predict_all(kind, prefix="v2_1_", mesh="v2", sim_kwargs=kw)
        validate.validate_cv(kind, prefix="v2_1_", mesh="v2", sim_kwargs=kw)
        diagnostics.run(kind, prefix="v2_1_")
        relaxation.posthoc_2c_rests(kind, prefix="v2_1_", sim_kwargs=kw)
    print("v2-6. v1 vs v2 vs v2.1 comparison and figures")
    relaxation.comparison(kinds)
    for p in relaxation.make_figures(headline, kinds):
        print("   ", p)


def write_summary(kinds):
    """results/summary.json (v1) and results/v2_summary.json: the numbers the README quotes."""
    _write_summary(kinds, "")
    for prefix in ("v2_", "v2_1_"):
        if all((RESULTS_DIR / f"{prefix}validation_{k.lower()}.json").exists() for k in kinds):
            _write_summary(kinds, prefix)


def _write_summary(kinds, prefix):
    out = {}
    for kind in kinds:
        val = json.load(open(RESULTS_DIR / f"{prefix}validation_{kind.lower()}.json"))
        fitj = json.load(open(RESULTS_DIR / f"{prefix}fit_{kind.lower()}.json"))
        cv = json.load(open(RESULTS_DIR / f"{prefix}cv_{kind.lower()}.json"))
        ov = json.load(open(RESULTS_DIR / f"{prefix}overpotentials_{kind.lower()}.json"))
        out[kind] = {
            "headline": val["headline"],
            "training_rmse_mV": 1000 * (fitj.get("training_rmse_curves_V") or fitj["training_rmse_V"]),
            "parameters": fitj["parameters"],
            "rates": {r: {"role": e["role"], "discharge": e["discharge"]} for r, e in val["rates"].items()},
            "cv": {s: {k: v for k, v in e.items() if k != "scan_rate_V_s"} for s, e in cv["scans"].items()},
            "b_value": cv["b_value"],
            "randles_sevcik": cv["randles_sevcik"],
            "two_c_overpotentials": ov["two_c"],
        }
        if prefix:
            out[kind]["training_rmse_rests_mV"] = 1000 * fitj["training_rmse_rests_V"]
            out[kind]["log10_correlation"] = fitj["log10_correlation"]
            out[kind]["weighting"] = fitj["weighting"]
    with open(RESULTS_DIR / f"{prefix}summary.json", "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(fast="--fast" in sys.argv, v1="--v2-only" not in sys.argv, v2="--v1-only" not in sys.argv)
