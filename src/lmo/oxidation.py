"""v3: a parasitic electrolyte-oxidation current at the positive electrode.

The cell loses charge to a side reaction near the 4.5 V cut-off (first-cycle coulombic
efficiency 91 % at C/10), and the rests the v2.1 model misses are the ones that start at that
cut-off. This module adds one anodic Tafel reaction on the LMO particle surface, in parallel
with the intercalation reaction:

    j_ox = j_ref * exp(alpha F (phi_s - phi_e - U_ref) / (R T))      [A per m2 of particle surface]

with U_ref fixed at 4.5 V, so that j_ref is the parasitic current density at the reference
potential and alpha sets the Tafel slope (59 mV per decade divided by alpha). Only those two
parameters are fitted; a free onset potential would be degenerate with j_ref. The reaction
consumes nothing the model tracks (the oxidised species is the electrolyte), so its only
effects are that part of the applied charging current never enters the particles, and that at
open circuit the oxidation is balanced by lithium re-intercalation (self-discharge), which pulls
the voltage down with a Tafel-limited, log(t)-like decay.

Implementation: PyBaMM's SEI machinery already provides a second interfacial reaction on the
positive electrode (the "SEI" option accepts a 2-tuple since v24.1), but all its built-in
kinetics are cathodic. ``OxidationCurrent`` subclasses ``pybamm.sei.SEIGrowth`` and replaces its
reaction rate with the anodic law above; the "SEI concentration" state is kept as a counter of
the charge passed by the side reaction, so ``"Loss of capacity to positive SEI [A.h]"`` is the
cumulative parasitic charge. The lithium-metal counter electrode gets no side reaction.
"""

from __future__ import annotations

import pybamm

OPTION = "positive electrode oxidation"  # lmo-level option handled by make_model: "none" or "Tafel"
U_REF_V = 4.5
PARAM_J_REF = "Positive electrode oxidation current density at reference potential [A.m-2]"
PARAM_ALPHA = "Positive electrode oxidation transfer coefficient"
PARAM_U_REF = "Positive electrode oxidation reference potential [V]"
PARASITIC_CHARGE_VARIABLE = "Loss of capacity to positive SEI [A.h]"
OXIDATION_CURRENT_VARIABLE = "Positive electrode SEI interfacial current density [A.m-2]"

# PyBaMM options that switch the second reaction on at the positive electrode only.
PYBAMM_OPTIONS = {"SEI": ("none", "reaction limited"), "SEI film resistance": "none"}
# Book-keeping parameters the SEI base class needs; they have no physical role here (the
# concentration state only counts charge) and are never fitted.
DUMMY_SEI_PARAMETERS = {
    "Initial SEI thickness [m]": 1.0e-8,  # with the molar volume below this scales the charge counter to ~300 mol/m3
    "SEI partial molar volume [m3.mol-1]": 1.0e-5,
    "Ratio of lithium moles to SEI moles": 1.0,
}


class OxidationCurrent(pybamm.sei.SEIGrowth):
    """Anodic Tafel side reaction on the positive particle surface (see module docstring)."""

    def get_coupled_variables(self, variables):
        domain, Domain = self.domain_Domain
        T = variables[f"{Domain} electrode temperature [K]"]
        delta_phi = variables[f"{Domain} electrode surface potential difference [V]"]
        j_ref = pybamm.Parameter(PARAM_J_REF)
        alpha = pybamm.Parameter(PARAM_ALPHA)
        U_ref = pybamm.Parameter(PARAM_U_REF)
        F_RT = self.param.F / (self.param.R * T)
        j_ox = j_ref * pybamm.exp(alpha * F_RT * (delta_phi - U_ref))
        variables.update(self._get_standard_reaction_variables(j_ox))
        variables.update(pybamm.sei.BaseModel.get_coupled_variables(self, variables))
        return variables

    def set_rhs(self, variables):
        # The parent integrates -a j / (F z) so that a cathodic SEI current builds up SEI; here
        # the state counts the (anodic) parasitic charge, so the sign is flipped and
        # "Loss of capacity to positive SEI [A.h]" is the charge passed by the oxidation.
        super().set_rhs(variables)
        (c_sei, rhs), = self.rhs.items()
        self.rhs = {c_sei: -rhs}


def oxidation_parameters(j_ref, alpha, u_ref: float = U_REF_V) -> dict:
    return {PARAM_J_REF: j_ref, PARAM_ALPHA: alpha, PARAM_U_REF: u_ref, **DUMMY_SEI_PARAMETERS}


def attach(model: pybamm.BaseModel) -> pybamm.BaseModel:
    """Replace the positive-electrode SEI submodel of an unbuilt model by the oxidation current
    and build it. The model must have been created with ``PYBAMM_OPTIONS`` and ``build=False``."""
    key = "positive primary sei"
    assert key in model.submodels, "model was not created with the positive-electrode SEI option"
    model.submodels[key] = OxidationCurrent(model.param, "positive", model.options, "primary")
    model.build_model()
    return model
