"""v2.1: model variants that give the rest a relaxation the v2 model does not have.

v2 showed that the v1 model (single particle size, no double-layer capacitance, every
interfacial overpotential instantaneous) cannot use the 30 s rests: the measured voltage keeps
relaxing like log(t) from 50 ms to 30 s while the model only has an instantaneous jump and one
sqrt(t) diffusion tail. Two standard PyBaMM options add exactly the missing physics, and both
are supported by published work on relaxation after a current interruption:

* ``"surface form": "differential"`` adds the double-layer capacitance, so the charge-transfer
  overpotential (and, through the porous electrode, the electrolyte potential) relaxes with an
  RC time scale instead of vanishing at once. In the intermittent-current-interruption (ICI)
  literature this R_ct-C_dl process is what sets the voltage between ~2 ms and ~1 s after the
  interruption (Geng, Thiringer and Lacey, IEEE Trans. Transp. Electrif. 2021).
* ``"particle size": "distribution"`` (the MP-DFN / MPM of Kirk, Please and Chapman, J.
  Electrochem. Soc. 2021) spreads the solid-diffusion time constant over a lognormal
  particle-size distribution; Kirk et al. showed this is what the DFN lacks to reproduce the
  slow voltage relaxation after a discharge.

The variants add one fitted parameter each: an effective double-layer capacity C_dl (F/m2 of
particle surface, which also absorbs the unknown real surface area) and/or the relative
standard deviation sd_p of the area-weighted lognormal radius distribution (mean fixed at the
5 um of v1). Everything else is as in v2: same disc, same training data (C/5 and C/2 curves plus
the C/10, C/5, C/2 cycle-1 rests at lambda = 1), 2C never used.
"""

from __future__ import annotations

import json
import os
import time
from multiprocessing import Pool

import numpy as np
import pybamm
from scipy.optimize import least_squares

from . import PROCESSED_DIR, RESULTS_DIR
from .fit import BOUNDS_LOG10, PARAM_NAMES, covariance_report
from .fit_v2 import HEADLINE_LAMBDA, REST_KINDS, REST_TRAINING, TRAINING, V2Residual, block_label, rest_jumps
from .model import DEFAULT_DISC_MM, DEFAULT_FIT, PARTICLE_RADIUS_M

VARIANTS = {
    "dl": {"options": {"surface form": "differential"}, "extra": ("C_dl",)},
    "psd": {"options": {"particle size": "distribution"}, "extra": ("sd_p",)},
    "psd+dl": {"options": {"surface form": "differential", "particle size": "distribution"}, "extra": ("C_dl", "sd_p")},
}
EXTRA_BOUNDS_LOG10 = {"C_dl": (-2.0, 3.0), "sd_p": (-1.5, 0.3)}  # C_dl in F/m2 of particle surface; sd_p relative to the mean radius
EXTRA_DEFAULT = {"C_dl": 1.0, "sd_p": 0.5}
# The v2.1 cost surface has several local minima (the double layer trades against R_contact and
# k), so the SPMe fit is run from these starting points and the lowest cost is kept; the DFN
# starts from the SPMe optimum (v1 and v2 showed the two fidelities land within a few percent).
STARTS = [
    {"Ds": 1e-14, "k": 1e-11, "R_contact": 20.0, "C_dl": 1.0, "sd_p": 0.5},
    {"Ds": 1e-14, "k": 1e-11, "R_contact": 70.0, "C_dl": 5.0, "sd_p": 0.3},
    {"Ds": 1e-13, "k": 1e-9, "R_contact": 60.0, "C_dl": 0.2, "sd_p": 0.8},
    {"Ds": 3e-15, "k": 3e-11, "R_contact": 40.0, "C_dl": 20.0, "sd_p": 0.5},
]
R_MIN_REL, R_MAX_REL = 0.02, 4.0  # particle-size axis, relative to the mean radius
N_SIZES = 20
PREFIX = "v2_1_"
# The variant carried through validation. Chosen on the training data alone (curve and rest
# RMSE of the three variants, results/v2_1_variants_<model>.json), never on 2C.
HEADLINE_VARIANT = "dl"


def lognormal_area_weighted(R, R_av, sd_rel):
    """Area-weighted lognormal radius distribution with mean R_av and standard deviation
    sd_rel * R_av, normalised on [R_MIN_REL, R_MAX_REL] * R_av. Works with pybamm symbols so
    sd_rel can be an InputParameter."""
    s2 = pybamm.log(1 + sd_rel**2)
    sig = pybamm.sqrt(s2)
    mu = pybamm.log(R_av) - s2 / 2
    pdf = pybamm.exp(-((pybamm.log(R) - mu) ** 2) / (2 * s2)) / (R * sig * np.sqrt(2 * np.pi))
    cdf = lambda r: 0.5 * (1 + pybamm.erf((pybamm.log(r) - mu) / (sig * np.sqrt(2))))
    return pdf / (cdf(R_MAX_REL * R_av) - cdf(R_MIN_REL * R_av))


def param_hook(variant: str):
    extra = VARIANTS[variant]["extra"]

    def hook(params: pybamm.ParameterValues):
        if "sd_p" in extra:
            R_av = PARTICLE_RADIUS_M
            params.update(
                {
                    "Positive minimum particle radius [m]": R_MIN_REL * R_av,
                    "Positive maximum particle radius [m]": R_MAX_REL * R_av,
                    "Positive area-weighted particle-size distribution [m-1]": lambda R: lognormal_area_weighted(R, R_av, pybamm.InputParameter("sd_p")),
                },
                check_already_exists=False,
            )
        if "C_dl" in extra:
            params.update({"Positive electrode double-layer capacity [F.m-2]": pybamm.InputParameter("C_dl")}, check_already_exists=False)

    return hook


def sim_kwargs(variant: str) -> dict:
    v = VARIANTS[variant]
    kw = {"model_options": v["options"], "param_hook": param_hook(variant), "extra_inputs": v["extra"]}
    if "sd_p" in v["extra"]:
        kw["var_pts_extra"] = {"R_p": N_SIZES}
    return kw


def residual(variant: str, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, lam=HEADLINE_LAMBDA, mode="trace") -> V2Residual:
    names = list(PARAM_NAMES) + list(VARIANTS[variant]["extra"])
    return V2Residual(kind, disc_mm, lam, mode, sim_kwargs=sim_kwargs(variant), param_names=names)


def fit_variant(variant: str, kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, x0: dict | None = None, blocks=None, verbose=True, lam=HEADLINE_LAMBDA, starts=None):
    """Least squares from ``x0`` or, with ``starts`` (a list of dicts), from each start keeping the
    lowest cost. The result carries ``starts_tried`` (cost per start) for the record."""
    res = residual(variant, kind, disc_mm, lam)
    bounds = {**BOUNDS_LOG10, **EXTRA_BOUNDS_LOG10}
    lb = np.array([bounds[n][0] for n in res.param_names])
    ub = np.array([bounds[n][1] for n in res.param_names])
    t0 = time.time()
    best, tried = None, []
    for start in starts or [x0]:
        x = {**DEFAULT_FIT, **EXTRA_DEFAULT, **(start or {})}
        theta0 = np.array([np.log10(x[n]) for n in res.param_names])
        r = least_squares(lambda th: res.weighted(th, blocks), theta0, bounds=(lb, ub), diff_step=1e-3, xtol=1e-5, ftol=1e-6, verbose=2 if verbose else 0)
        tried.append({"start": {n: x[n] for n in res.param_names}, "cost": float(np.sum(r.fun**2)), "optimum": {n: float(v) for n, v in res.unpack(r.x).items()}})
        if best is None or tried[-1]["cost"] < best_cost:
            best, best_cost = r, tried[-1]["cost"]
    best.starts_tried = tried
    return best, res, time.time() - t0


def _loo_refit(args):
    variant, kind, disc_mm, best, keep, label = args
    r, res, _ = fit_variant(variant, kind, disc_mm, x0=best, blocks=keep, verbose=False)
    d = res.unpack(r.x)
    d["dropped"] = label
    return d


def summarise(result, res: V2Residual, variant: str, elapsed: float, loo=None) -> dict:
    fit = res.unpack(result.x)
    rep = covariance_report(result)
    rm = res.rmse_by_type(fit)
    names = res.param_names
    out = {
        "version": "v2.1",
        "variant": variant,
        "model_options": VARIANTS[variant]["options"],
        "model": res.kind,
        "disc_mm": res.disc_mm,
        "training_curves": [f"{r} cycle {c} charge+discharge" for r, c in TRAINING],
        "training_rests": [f"{r} cycle {c} rest after {a}" for r, c in REST_TRAINING for a in REST_KINDS],
        "two_c_used": False,
        "particle_radius_m": PARTICLE_RADIUS_M,
        "weighting": {"lambda": res.lam, "rest_weight": res.w_rest, "rest_points": "every logged sample of each rest", "n_curve_points": res.n_curve, "n_rest_points": res.n_rest},
        "training_rmse_curves_V": rm["curves_mV"] / 1000,
        "training_rmse_rests_V": rm["rests_mV"] / 1000,
        "training_rmse_rest_traces_V": rm["rest_traces_mV"] / 1000,
        "training_rmse_per_block_mV": rm["per_block_mV"],
        "weighted_rmse_V": rep["rmse_V"],
        "n_residuals": int(len(result.fun)),
        "fit_time_s": elapsed,
        "converged": bool(result.success),
        "starts_tried": getattr(result, "starts_tried", None),
        "parameters": {},
        "log10_correlation": {f"{a}-{b}": float(rep["corr"][i, j]) for i, a in enumerate(names) for j, b in enumerate(names) if i < j},
    }
    for i, n in enumerate(names):
        sd = rep["sd_log10"][i]
        entry = {
            "value": fit[n],
            "log10": float(result.x[i]),
            "sd_log10_jacobian": float(sd),
            "ci95_low_jacobian": float(10 ** (result.x[i] - 1.96 * sd)),
            "ci95_high_jacobian": float(10 ** (result.x[i] + 1.96 * sd)),
        }
        if loo:
            vals = np.array([d[n] for d in loo])
            entry["leave_one_step_out_min"] = float(vals.min())
            entry["leave_one_step_out_max"] = float(vals.max())
        out["parameters"][n] = entry
    tau = PARTICLE_RADIUS_M**2 / fit["Ds"]
    sd = rep["sd_log10"][0]
    out["parameters"]["tau_R2_over_Ds_s"] = {"value": tau, "ci95_low_jacobian": float(tau / 10 ** (1.96 * sd)), "ci95_high_jacobian": float(tau * 10 ** (1.96 * sd))}
    if loo:
        taus = np.array([PARTICLE_RADIUS_M**2 / d["Ds"] for d in loo])
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_min"] = float(taus.min())
        out["parameters"]["tau_R2_over_Ds_s"]["leave_one_step_out_max"] = float(taus.max())
    out["parameters"]["R_contact"]["unit"] = "Ohm"
    out["parameters"]["Ds"]["unit"] = "m2/s"
    out["parameters"]["k"]["unit"] = "m^2.5/(mol^0.5 s) in j0 = k F sqrt(c_e c_s (c_max - c_s))"
    if "C_dl" in out["parameters"]:
        out["parameters"]["C_dl"]["unit"] = "F/m2 of particle surface (effective; absorbs the unknown real surface area)"
    if "sd_p" in out["parameters"]:
        out["parameters"]["sd_p"]["unit"] = "area-weighted lognormal standard deviation relative to the 5 um mean radius"
    out["rest_jumps"] = rest_jumps(res, fit)
    return out


def fit_with_loo(variant: str, kind: str, disc_mm: float = DEFAULT_DISC_MM, x0: dict | None = None, starts=None):
    result, res, elapsed = fit_variant(variant, kind, disc_mm, x0=x0, starts=starts)
    best = res.unpack(result.x)
    jobs = [(variant, kind, disc_mm, best, [b for b in res.blocks if b != drop], block_label(drop)) for drop in res.blocks]
    workers = int(os.environ.get("LMO_WORKERS", "1"))
    if workers > 1:
        with Pool(workers) as pool:
            loo = pool.map(_loo_refit, jobs)
    else:
        loo = [_loo_refit(j) for j in jobs]
    summary = summarise(result, res, variant, elapsed, loo)
    summary["leave_one_step_out"] = loo
    return result, res, summary


def _starts_for(variant: str, kind: str, x0: dict | None):
    """SPMe: the multi-start list. DFN: the SPMe optimum of the same variant when available."""
    if x0:
        return [x0]
    if kind != "SPMe":
        tag = variant.replace("+", "_")
        for name in (f"{PREFIX}variants_spme.json",):
            p = RESULTS_DIR / name
            if p.exists():
                table = json.load(open(p))["variants"]
                if variant in table:
                    return [table[variant]["parameters"]]
    return STARTS


def run(kind: str = "SPMe", disc_mm: float = DEFAULT_DISC_MM, variants=tuple(VARIANTS), headline: str = HEADLINE_VARIANT, x0: dict | None = None) -> dict:
    """Fit every variant (training data only) for the variant table, then the headline variant
    with leave-one-block-out refits; writes v2_1_variants_<model>.json and v2_1_fit_<model>.json."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    table = {"model": kind, "headline_variant": headline, "chosen_on": "training curves and rests only (never 2C)", "variants": {}}
    headline_summary = None
    for variant in variants:
        starts = _starts_for(variant, kind, x0)
        if variant == headline:
            result, res, summary = fit_with_loo(variant, kind, disc_mm, starts=starts)
            headline_summary = summary
            with open(RESULTS_DIR / f"{PREFIX}fit_{kind.lower()}.json", "w") as f:
                json.dump(summary, f, indent=2)
            res.rest_frames(res.unpack(result.x)).to_csv(PROCESSED_DIR / f"{PREFIX}model_rests_{kind.lower()}.csv", index=False)
        else:
            result, res, elapsed = fit_variant(variant, kind, disc_mm, starts=starts, verbose=False)
            summary = summarise(result, res, variant, elapsed)
        table["variants"][variant] = {
            "model_options": VARIANTS[variant]["options"],
            "parameters": {n: summary["parameters"][n]["value"] for n in res.param_names},
            "tau_R2_over_Ds_s": summary["parameters"]["tau_R2_over_Ds_s"]["value"],
            "training_rmse_curves_mV": 1000 * summary["training_rmse_curves_V"],
            "training_rmse_rests_mV": 1000 * summary["training_rmse_rests_V"],
            "k_R_contact_correlation": summary["log10_correlation"]["k-R_contact"],
            "converged": summary["converged"],
            "fit_time_s": summary["fit_time_s"],
        }
        with open(RESULTS_DIR / f"{PREFIX}variants_{kind.lower()}.json", "w") as f:
            json.dump(table, f, indent=2)
    headline_summary["variants"] = table["variants"]
    return headline_summary


if __name__ == "__main__":
    import sys

    kind = sys.argv[1] if len(sys.argv) > 1 else "SPMe"
    out = run(kind)
    print(json.dumps({k: out[k] for k in ("parameters", "log10_correlation", "training_rmse_curves_V", "training_rmse_rests_V")}, indent=2))
