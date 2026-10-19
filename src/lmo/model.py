"""PyBaMM half-cell model of the LMO electrode against lithium metal.

Starts from the ``Xu2019`` half-cell parameter set (lithium-metal counter electrode, separator,
1 M LiPF6 carbonate electrolyte with Valoen2005 transport) and replaces the positive electrode
with the cell measured here: pseudo-OCV from the C/10 data, 85 um coating, LMO density
4.28 g/cm3, 94.2 wt% active, nominal capacity 1.434 mAh, 25 C.

The electrode area is unknown (coin-cell disc). For a chosen disc diameter the active-material
volume fraction follows from the active mass and the coating volume, and the porosity from what
is left after active material and binder/carbon. 14 mm is the default; 12 and 16 mm are run as
a sensitivity check.

Three effective parameters are fitted: particle diffusion time constant (R^2/Ds, implemented by
fitting Ds at a fixed radius), the reaction rate constant k of the Butler-Volmer exchange
current density, and a lumped contact resistance.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pybamm

from . import NOMINAL_CAPACITY_MAH, PROCESSED_DIR
from .ocv import LMO_C_MAX_MOL_M3, LMO_DENSITY_KG_M3, active_mass_g

COATING_THICKNESS_M = 85e-6
ACTIVE_WT_FRACTION = 0.942
INACTIVE_DENSITY_KG_M3 = 1800.0  # PVDF + carbon black blend, assumed
DEFAULT_DISC_MM = 14.0
PARTICLE_RADIUS_M = 5e-6  # fixed; only the time constant R^2/Ds is identifiable
T_K = 298.15
F = 96485.33212

# Default (starting) values of the fitted effective parameters.
DEFAULT_FIT = {"Ds": 1e-14, "k": 1e-11, "R_contact": 20.0}


def disc_geometry(disc_mm: float = DEFAULT_DISC_MM) -> dict:
    """Area, active volume fraction and porosity consistent with the capacity for a disc."""
    area = np.pi * (disc_mm * 1e-3 / 2) ** 2
    m_active = active_mass_g() * 1e-3
    m_inactive = m_active * (1 - ACTIVE_WT_FRACTION) / ACTIVE_WT_FRACTION
    v_coating = area * COATING_THICKNESS_M
    eps_s = m_active / LMO_DENSITY_KG_M3 / v_coating
    eps_inactive = m_inactive / INACTIVE_DENSITY_KG_M3 / v_coating
    porosity = 1.0 - eps_s - eps_inactive
    if porosity <= 0.05:
        raise ValueError(f"disc {disc_mm} mm leaves porosity {porosity:.2f}; capacity not consistent")
    return {
        "disc_mm": disc_mm,
        "area_m2": area,
        "areal_loading_mg_cm2": m_active * 1e6 / (area * 1e4),
        "active_volume_fraction": eps_s,
        "inactive_volume_fraction": eps_inactive,
        "porosity": porosity,
        "coating_density_g_cm3": (m_active + m_inactive) * 1e3 / (v_coating * 1e6),
    }


def load_ocp_table() -> np.ndarray:
    return np.loadtxt(PROCESSED_DIR / "lmo_ocp_vs_sto.csv", delimiter=",", skiprows=1)


_OCP_TABLE = None


def ocp_table() -> np.ndarray:
    global _OCP_TABLE
    if _OCP_TABLE is None:
        _OCP_TABLE = load_ocp_table()
    return _OCP_TABLE


def lmo_ocp(sto):
    table = ocp_table()
    return pybamm.Interpolant(table[:, 0], table[:, 1], sto, name="LMO measured OCP", interpolator="linear", extrapolate=True)


def ocv_to_sto(voltage: float) -> float:
    """Invert the measured pseudo-OCV (clipped to the measured window)."""
    table = ocp_table()
    # Voltage falls with stoichiometry, so reverse both columns for np.interp.
    v = table[::-1, 1]
    s = table[::-1, 0]
    # Keep strictly inside (0, 1): the exchange current contains sqrt(c_max - c_s).
    return float(np.clip(np.interp(voltage, v, s), 0.01, 0.995))


def ocv_of_sto(sto):
    table = ocp_table()
    return np.interp(sto, table[:, 0], table[:, 1])


def make_exchange_current(k):
    def j0(c_e, c_s_surf, c_s_max, T):
        # Guarded so that a sweep slightly past the measured OCV window (CV to 3.0 V) never
        # takes the square root of a negative number.
        floor = 1e-3 * c_s_max
        c_s = pybamm.maximum(pybamm.minimum(c_s_surf, c_s_max - floor), floor)
        return k * F * c_e**0.5 * c_s**0.5 * (c_s_max - c_s) ** 0.5

    return j0


def build_parameters(
    Ds: float = DEFAULT_FIT["Ds"],
    k: float = DEFAULT_FIT["k"],
    R_contact: float = DEFAULT_FIT["R_contact"],
    disc_mm: float = DEFAULT_DISC_MM,
    initial_sto: float = 1.0,
    v_min: float = 2.5,
    v_max: float = 5.0,
) -> pybamm.ParameterValues:
    """Parameter set for the cell. ``v_min``/``v_max`` are the model-level safety cut-offs;
    the 3.5 and 4.5 V protocol limits are applied by the experiment steps (C/10 started from
    3.33 V, below the protocol window)."""
    geo = disc_geometry(disc_mm)
    side = np.sqrt(geo["area_m2"])
    p = pybamm.ParameterValues("Xu2019")
    p.update(
        {
            "Nominal cell capacity [A.h]": NOMINAL_CAPACITY_MAH / 1000,
            "Current function [A]": NOMINAL_CAPACITY_MAH / 1000 / 5,
            "Electrode height [m]": side,
            "Electrode width [m]": side,
            "Positive electrode thickness [m]": COATING_THICKNESS_M,
            "Positive electrode active material volume fraction": geo["active_volume_fraction"],
            "Positive electrode porosity": geo["porosity"],
            "Maximum concentration in positive electrode [mol.m-3]": LMO_C_MAX_MOL_M3,
            "Initial concentration in positive electrode [mol.m-3]": initial_sto * LMO_C_MAX_MOL_M3,
            "Positive electrode OCP [V]": lmo_ocp,
            "Positive electrode OCP entropic change [V.K-1]": 0.0,
            "Positive particle radius [m]": PARTICLE_RADIUS_M,
            "Positive particle diffusivity [m2.s-1]": Ds,
            "Positive electrode exchange-current density [A.m-2]": make_exchange_current(k),
            "Positive electrode conductivity [S.m-1]": 10.0,
            "Positive electrode charge transfer coefficient": 0.5,
            "Contact resistance [Ohm]": R_contact,
            "Lower voltage cut-off [V]": v_min,
            "Upper voltage cut-off [V]": v_max,
            "Open-circuit voltage at 0% SOC [V]": 3.5,
            "Open-circuit voltage at 100% SOC [V]": 4.5,
            "Ambient temperature [K]": T_K,
            "Initial temperature [K]": T_K,
            "Reference temperature [K]": T_K,
        },
        check_already_exists=False,
    )
    return p


def make_model(kind: str = "SPMe", **options) -> pybamm.BaseModel:
    opts = {"working electrode": "positive", "contact resistance": "true"}
    opts.update(options)
    cls = {"SPMe": pybamm.lithium_ion.SPMe, "DFN": pybamm.lithium_ion.DFN, "SPM": pybamm.lithium_ion.SPM}[kind]
    return cls(opts, name=f"{kind} half cell")


def default_var_pts():
    return {"x_n": 10, "x_s": 10, "x_p": 30, "r_n": 10, "r_p": 30}


OUTPUT_VARIABLES = [
    "Time [s]",
    "Voltage [V]",
    "Current [A]",
    "Discharge capacity [A.h]",
    "Battery open-circuit voltage [V]",
    "Bulk open-circuit voltage [V]",
    "Surface open-circuit voltage [V]",
    "Battery particle concentration overpotential [V]",
    "X-averaged battery reaction overpotential [V]",
    "X-averaged battery concentration overpotential [V]",
    "X-averaged battery electrolyte ohmic losses [V]",
    "X-averaged battery solid phase ohmic losses [V]",
    "Contact overpotential [V]",
    "Lithium metal interface reaction overpotential [V]",
]


def gcd_experiment(
    current_mA: float | tuple[float, float],
    n_cycles: int = 1,
    rest_s: float = 30.0,
    v_min=3.5,
    v_max=4.5,
    start="charge",
    period_s: float | None = None,
):
    """Charge -> rest -> discharge -> rest, repeated. Currents in mA, EC-Lab style.

    ``current_mA`` may be a (charge, discharge) pair because the C/10 run read different
    currents in the two directions. ``period_s`` is the output interval; by default it is
    chosen so a full step has roughly 1500 points.
    """
    i_chg, i_dis = (current_mA, current_mA) if np.isscalar(current_mA) else current_mA
    if period_s is None:
        period_s = max(1.0, NOMINAL_CAPACITY_MAH / abs(i_dis) * 3600 / 1500)
    cyc = [
        f"Charge at {abs(i_chg)} mA until {v_max} V",
        f"Rest for {rest_s} seconds",
        f"Discharge at {abs(i_dis)} mA until {v_min} V",
        f"Rest for {rest_s} seconds",
    ]
    if start == "discharge":
        cyc = cyc[2:] + cyc[:2]
    return pybamm.Experiment(cyc * n_cycles, period=f"{period_s} seconds")


def make_simulation(experiment, params, kind: str = "SPMe", model=None) -> pybamm.Simulation:
    model = model or make_model(kind)
    solver = pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8)
    return pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=solver,
        var_pts=default_var_pts(),
    )


def solve(experiment, params, kind: str = "SPMe", model=None, inputs=None) -> pybamm.Solution:
    return make_simulation(experiment, params, kind, model).solve(inputs=inputs)


def step_frame(sol_step, kind: str, rate: str | None = None) -> pd.DataFrame:
    """Tabulate one experiment step with its voltage components."""
    t = sol_step["Time [s]"].entries
    out = {"time_s": t - t[0], "voltage_V": sol_step["Voltage [V]"].entries, "current_A": sol_step["Current [A]"].entries}
    out["capacity_mAh"] = np.concatenate(
        [[0.0], np.cumsum(0.5 * (np.abs(out["current_A"][1:]) + np.abs(out["current_A"][:-1])) * np.diff(t)) / 3.6]
    )
    for name, key in [
        ("ocv_V", "Battery open-circuit voltage [V]"),
        ("bulk_ocv_V", "Bulk open-circuit voltage [V]"),
        ("eta_particle_V", "Battery particle concentration overpotential [V]"),
        ("eta_reaction_V", "X-averaged battery reaction overpotential [V]"),
        ("eta_electrolyte_conc_V", "X-averaged battery concentration overpotential [V]"),
        ("eta_electrolyte_ohmic_V", "X-averaged battery electrolyte ohmic losses [V]"),
        ("eta_solid_ohmic_V", "X-averaged battery solid phase ohmic losses [V]"),
        ("eta_contact_V", "Contact overpotential [V]"),
        ("eta_li_metal_V", "Lithium metal interface reaction overpotential [V]"),
    ]:
        try:
            out[name] = sol_step[key].entries
        except KeyError:
            out[name] = np.full_like(t, np.nan)
    df = pd.DataFrame(out)
    df["step"] = kind
    if rate:
        df["rate"] = rate
    return df


def simulate_gcd(
    current_mA: float,
    start_voltage: float,
    fit: dict | None = None,
    kind: str = "SPMe",
    disc_mm: float = DEFAULT_DISC_MM,
    n_cycles: int = 1,
    rate: str | None = None,
    model=None,
):
    """Simulate the EC-Lab GCD protocol from a rest voltage. Returns (solution, list of step frames)."""
    fit = {**DEFAULT_FIT, **(fit or {})}
    sto0 = ocv_to_sto(start_voltage)
    params = build_parameters(disc_mm=disc_mm, initial_sto=sto0, **fit)
    sol = solve(gcd_experiment(current_mA, n_cycles), params, kind, model=model)
    return sol, solution_frames(sol, n_cycles, rate)


def solution_frames(sol, n_cycles: int = 1, rate: str | None = None) -> list[pd.DataFrame]:
    """Charge and discharge step frames of a charge/rest/discharge/rest experiment solution."""
    frames = []
    kinds = ["charge", "rest", "discharge", "rest"] * n_cycles
    steps = [s for c in sol.cycles for s in c.steps]
    for i, step in enumerate(steps):
        k = kinds[i % 4]
        if k == "rest":
            continue
        df = step_frame(step, k, rate)
        df["cycle"] = i // 4 + 1
        frames.append(df)
    return frames
