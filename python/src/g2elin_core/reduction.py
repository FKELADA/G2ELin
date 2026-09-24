"""The catalogue of model-order reductions: which states each element type
can give up, what the named levels are, and how a choice resolves to the
per-state modes :func:`g2elin_core.components.base.apply_reduction` wants.

This module holds *no* physics. It is a description of the choices, shared
by the symbolic model builders (``components/*.py``), the schema that
persists a choice (``network/schema.py``), the validation that rejects an
inconsistent one (``network/validation.py``), and the API that shows the
choices to the web UI. Keeping it in one place is what stops the level
names, the state groupings and the UI's checkboxes drifting apart.

**Granularity.** The unit of control is a *state group*, not an individual
state: the d- and q-axis halves of a stator flux, or the two integrators of
a current PI, are one physical approximation each and there is no
meaningful model in which one is kept and the other dropped. Groups are
also what make the named levels ("6th order", "droop only") expressible as
a handful of settings rather than a list of nineteen.

**Named levels vs. free choice.** A level is nothing more than a preset
mapping of groups to modes; anything a level can express, a set of explicit
group overrides can express too, and overrides are applied *on top of* a
level. So the UI can offer both the familiar names and the underlying
per-group ``dynamic | algebraic | frozen`` control without the two being
different mechanisms.

See :func:`g2elin_core.components.base.apply_reduction` for what the three
modes mean and why ``algebraic`` and ``frozen`` are not interchangeable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

from g2elin_core.components.base import ALGEBRAIC, DYNAMIC, FROZEN

__all__ = [
    "ALGEBRAIC", "DYNAMIC", "FROZEN", "MODES", "ElementModel", "StateGroup",
    "ELEMENTS", "element", "resolve_modes", "level_ids", "group_ids", "describe",
    "CATEGORIES", "SYNCHRONISATION", "CONTROL", "UNIT_ELECTRICAL", "NETWORK", "MIXED", "REFERENCE",
    "category_of_state", "group_of_state", "split_state_name",
    "sm_element", "gfm_element", "AVR_GROUPS", "PSS_GROUPS", "GFM_OUTER_GROUPS",
]

MODES = (DYNAMIC, ALGEBRAIC, FROZEN)

# What kind of phenomenon a state belongs to. Every state group carries one,
# so a mode can be classified by which category its participation falls in
# (see :func:`g2elin_core.modal.classify.classify_modes`) -- the eigenvalue
# map groups modes by it.
#
# The split is by *mechanism*, not by which box the state sits in. A
# converter has no rotor, but its droop law and a machine's swing equation
# do the same job -- they are how the unit stays in step with the grid -- so
# they share a category, and a question like "which modes are
# synchronisation modes?" has one answer across a mixed fleet.
SYNCHRONISATION = "synchronisation"
CONTROL = "control"
UNIT_ELECTRICAL = "unit_electrical"
NETWORK = "network"
MIXED = "mixed"
REFERENCE = "reference"

CATEGORIES: dict[str, dict[str, str]] = {
    SYNCHRONISATION: {
        "label": "Synchronisation",
        "note": "How a unit stays in step with the grid: a machine's swing equation, a "
                "grid-forming converter's droop angle and power filters, a grid-following "
                "converter's PLL. The slow modes -- local and inter-area oscillations live here.",
    },
    CONTROL: {
        "label": "Control",
        "note": "Regulators: AVR, exciter, PSS, governor, and a converter's cascaded voltage "
                "and current loops. Where a tuning change shows up first.",
    },
    UNIT_ELECTRICAL: {
        "label": "Unit electrical",
        "note": "The electromagnetic states inside a unit: stator and rotor flux, damper "
                "windings, filter currents and voltages, the DC link, the step-up transformer "
                "current. Fast, and the first thing model-order reduction removes.",
    },
    NETWORK: {
        "label": "Network",
        "note": "Bus voltages and branch currents -- the grid's own electromagnetic "
                "transients. Absent entirely from a quasi-stationary model.",
    },
    MIXED: {
        "label": "Mixed",
        "note": "No single category accounts for most of this mode: it couples across "
                "them. Often the interesting ones -- a control loop interacting with the "
                "network, for instance.",
    },
    REFERENCE: {
        "label": "Reference angle",
        "note": "Not a mode of the system at all: the model's dq frame has no absolute "
                "position, so turning every angle together changes nothing. Marginal by "
                "construction, and excluded from the stability verdict.",
    },
}


@dataclass(frozen=True)
class StateGroup:
    """One physical approximation a user can switch on or off."""

    id: str
    label: str
    symbols: tuple[str, ...]  # the DAE symbol names this group covers
    states: tuple[str, ...]  # their display names, for the UI
    allowed: tuple[str, ...] = (DYNAMIC, ALGEBRAIC)
    default: str = DYNAMIC
    note: str = ""
    category: str = UNIT_ELECTRICAL  # see CATEGORIES
    # Groups that must be algebraic too before this one can be. A control
    # loop's integrator has no term of its own in its own equation --
    # ``0 = Ki*(ref - x)`` says nothing about the integrator, only about
    # ``x`` -- so making it algebraic is only solvable once whatever it
    # regulates is an unknown as well. These are the minimum sets, found by
    # building every combination against CIGRE (see tests/test_model_order).
    requires: tuple[str, ...] = ()

    @property
    def locked(self) -> bool:
        return self.allowed == (DYNAMIC,)


@dataclass(frozen=True)
class ElementModel:
    kind: str
    label: str
    groups: tuple[StateGroup, ...]
    levels: dict[str, dict[str, str]]  # level id -> {group id: mode}, dynamic omitted
    level_labels: dict[str, str]
    level_notes: dict[str, str] = field(default_factory=dict)
    default_level: str = "full"

    def group(self, group_id: str) -> StateGroup:
        for g in self.groups:
            if g.id == group_id:
                return g
        # A ValueError, not a KeyError: this is reached from a pydantic
        # validator on user-supplied settings, and pydantic turns a
        # ValueError into a proper per-field validation error (a 422 with
        # the reason) while a KeyError escapes as a 500.
        raise ValueError(
            f"{self.kind!r} has no state group {group_id!r} -- have {[g.id for g in self.groups]}"
        )

    def modes_by_group(self, level: str | None = None, overrides: dict[str, str] | None = None) -> dict[str, str]:
        """Every group's mode, after applying ``level`` then ``overrides``."""
        level = level or self.default_level
        if level not in self.levels:
            raise ValueError(
                f"unknown {self.kind} model level {level!r} -- available: {sorted(self.levels)}"
            )
        modes = {g.id: g.default for g in self.groups}
        modes.update(self.levels[level])
        for group_id, mode in (overrides or {}).items():
            group = self.group(group_id)  # raises on an unknown group id
            if mode not in group.allowed:
                raise ValueError(
                    f"{self.kind}.{group_id} cannot be {mode!r} -- allowed: {list(group.allowed)}"
                    + (f" ({group.note})" if group.note else "")
                )
            modes[group_id] = mode
        return modes

    def unmet_requirements(self, by_group: dict[str, str]) -> list[tuple[str, tuple[str, ...]]]:
        """``(group, groups it still needs)`` for every group made algebraic
        without what it depends on -- see :attr:`StateGroup.requires`.

        Checked up front rather than left to surface as a singular Jacobian
        deep inside the linearisation, which says nothing about which
        setting caused it.
        """
        out = []
        for g in self.groups:
            if by_group.get(g.id) != ALGEBRAIC or not g.requires:
                continue
            missing = tuple(r for r in g.requires if by_group.get(r) != ALGEBRAIC)
            if missing:
                out.append((g.id, missing))
        return out

    def modes(self, level: str | None = None, overrides: dict[str, str] | None = None) -> dict[str, str]:
        """The per-*symbol* mode map ``apply_reduction`` takes."""
        by_group = self.modes_by_group(level, overrides)
        return {sym: by_group[g.id] for g in self.groups for sym in g.symbols}

    def matching_level(self, by_group: dict[str, str]) -> str | None:
        """The named level a group->mode map corresponds to, if any -- what
        the UI shows in its level picker after per-group edits."""
        for level_id in self.levels:
            if self.modes_by_group(level_id) == by_group:
                return level_id
        return None


# --- Synchronous machine ------------------------------------------------------
# The 19 states of components/sm.py. The nine controller states (governor,
# PSS, AVR) are groups of their own and stay dynamic in every named level:
# a reduced-order machine in a stability study keeps its full controls --
# that is the whole point of the study. They are still exposed, because
# switching one off is the cleanest way to ask what it contributes.
# The AVR and PSS groups depend on which regulator the machine carries, so
# they are built per variant (see _sm_controls) and spliced onto the fixed
# machine groups below.
_SM_MACHINE_GROUPS = (
        StateGroup(
            "trafo_current", "Step-up transformer current", ("igd", "igq"), ("i_gd", "i_gq"),
            note="must be algebraic when the network is quasi-stationary",
        ),
        StateGroup("stator_flux", "Stator flux", ("phi_d", "phi_q"), ("psi_d", "psi_q")),
        StateGroup(
            "field_flux", "Field winding flux", ("phi_fd",), ("psi_fd",),
            allowed=(DYNAMIC, ALGEBRAIC, FROZEN),
            note="a slow state (T'do ~ 5-10 s): freeze it for the classical model, "
                 "do not make it algebraic unless you mean an instantaneous field winding",
        ),
        StateGroup("damper_1d", "d-axis damper winding", ("phi_1d",), ("psi_1d",)),
        StateGroup("damper_1q", "q-axis damper winding", ("phi_1q",), ("psi_1q",)),
        StateGroup("damper_2q", "Second q-axis damper winding", ("phi_2q",), ("psi_2q",)),
        StateGroup(
            "swing", "Rotor swing", ("Dwr", "theta"), ("dw_r", "theta"), allowed=(DYNAMIC,),
            note="every machine model keeps the swing equation; theta defines the rotor frame",
            category=SYNCHRONISATION,
        ),
)

# --- the regulator groups, one per model --------------------------------------
# A model set to "none" contributes no group at all, rather than an empty one:
# there is nothing to reduce, and a row with no states in it would be a
# control the user cannot use.
GOVERNOR_GROUPS = {
    "g2elin": StateGroup(
        "governor", "Governor", ("Pm",), ("P_m",), allowed=(DYNAMIC, ALGEBRAIC, FROZEN),
        note="frozen = fixed mechanical power", category=CONTROL,
    ),
    "none": None,
}
PSS_GROUPS = {
    "g2elin": StateGroup(
        "pss", "Power system stabiliser", ("Dw1", "v1", "v2", "vpss"),
        ("dw_1", "v_1", "v_2", "v_pss"), allowed=(DYNAMIC, FROZEN),
        note="frozen = PSS held at its equilibrium output", category=CONTROL,
    ),
    "kundur": StateGroup(
        "pss", "Power system stabiliser (Kundur)", ("v1", "v2", "vpss"),
        ("v_1", "v_2", "v_pss"), allowed=(DYNAMIC, FROZEN),
        note="washout and two lead-lag stages, no input filter "
             "(Kundur Fig. E12.9); frozen = held at its equilibrium output",
        category=CONTROL,
    ),
}
AVR_GROUPS = {
    "g2elin": StateGroup(
        "avr", "AVR / exciter", ("e1", "e2", "efd", "e3"), ("e_1", "e_2", "e_fd", "e_3"),
        allowed=(DYNAMIC, ALGEBRAIC, FROZEN), note="frozen = constant field voltage",
        category=CONTROL,
    ),
    "kundur": StateGroup(
        "avr", "AVR / exciter (Kundur)", ("e1", "e_tgr"), ("e_1", "e_tgr"),
        allowed=(DYNAMIC, FROZEN),
        note="thyristor exciter with transient gain reduction (Kundur Fig. "
             "E12.9). It has no exciter lag, so the field voltage is algebraic "
             "rather than a state -- which is also why freezing this group "
             "holds the regulator but not the field voltage, since the gain "
             "path still follows V_ref and the stabiliser. Freeze the PSS too "
             "for a genuinely constant field voltage.",
        category=CONTROL,
    ),
}
PSS_GROUPS["none"] = None

_SM_LEVELS = {
    "full": {},
    "order8": {"trafo_current": ALGEBRAIC},
    "order6": {"trafo_current": ALGEBRAIC, "stator_flux": ALGEBRAIC},
    "order5": {"trafo_current": ALGEBRAIC, "stator_flux": ALGEBRAIC, "damper_2q": ALGEBRAIC},
    "order4": {
        "trafo_current": ALGEBRAIC, "stator_flux": ALGEBRAIC,
        "damper_1d": ALGEBRAIC, "damper_2q": ALGEBRAIC,
    },
    "order3": {
        "trafo_current": ALGEBRAIC, "stator_flux": ALGEBRAIC,
        "damper_1d": ALGEBRAIC, "damper_1q": ALGEBRAIC, "damper_2q": ALGEBRAIC,
    },
    "order2": {
        "trafo_current": ALGEBRAIC, "stator_flux": ALGEBRAIC, "field_flux": FROZEN,
        "damper_1d": ALGEBRAIC, "damper_1q": ALGEBRAIC, "damper_2q": ALGEBRAIC,
    },
}

_SM_LEVEL_LABELS = {
    "full": "Full order (EMT)",
    "order8": "8th order",
    "order6": "6th order",
    "order5": "5th order",
    "order4": "4th order (two-axis)",
    "order3": "3rd order (flux decay)",
    "order2": "2nd order (classical)",
}

_SM_LEVEL_NOTES = {
    "full": "Everything dynamic, including the step-up transformer current. Needs a dynamic network.",
    "order8": "Full Park machine. Transformer current algebraic, so it pairs with either network.",
    "order6": "Stator flux algebraic -- the standard machine model for stability studies "
              "(Sauer-Pai / Anderson-Fouad). Pairs with a quasi-stationary network.",
    "order5": "6th order without the second q-axis damper.",
    "order4": "Two-axis model: field flux and one q-axis damper (E'd, E'q).",
    "order3": "Flux-decay / one-axis model: field flux and swing.",
    "order2": "Classical model: swing only, with the field flux *frozen* at its "
              "equilibrium (constant E'), not made algebraic.",
}


@lru_cache(maxsize=None)
def sm_element(
    exciter: str = "g2elin", pss: str = "g2elin", governor: str = "g2elin"
) -> ElementModel:
    """The catalogue for a machine carrying one particular set of regulators.

    Only the ``governor``, ``pss`` and ``avr`` groups differ: the group *ids*
    are the same whichever model is chosen, so a level, an override or a
    saved network stays meaningful when a machine's exciter is swapped --
    only the symbols behind the group change. A regulator set to ``none``
    contributes no group.
    """
    for name, chosen, table in (
        ("exciter", exciter, AVR_GROUPS), ("PSS", pss, PSS_GROUPS),
        ("governor", governor, GOVERNOR_GROUPS),
    ):
        if chosen not in table:
            raise ValueError(f"unknown {name} {chosen!r} -- have {sorted(table)}")
    chosen_groups = tuple(
        g for g in (GOVERNOR_GROUPS[governor], PSS_GROUPS[pss], AVR_GROUPS[exciter])
        if g is not None
    )
    return ElementModel(
        kind="sm",
        label="Synchronous machine",
        groups=_SM_MACHINE_GROUPS + chosen_groups,
        levels=_SM_LEVELS,
        level_labels=_SM_LEVEL_LABELS,
        level_notes=_SM_LEVEL_NOTES,
    )


_SM = sm_element()

# --- Grid-forming converter ---------------------------------------------------
# The 15 states of components/gfm.py, ordered outward from the grid: the
# levels peel off one cascaded loop at a time, ending at the droop law plus
# its power measurement filters, which is the model RMS tools carry.
_GFM = ElementModel(
    kind="gfm",
    label="Grid-forming converter",
    groups=(
        StateGroup(
            "trafo_current", "Transformer current", ("igd", "igq"), ("i_gd", "i_gq"),
            note="must be algebraic when the network is quasi-stationary",
        ),
        StateGroup("filter", "LC filter", ("isd", "isq", "ved", "veq"), ("i_sd", "i_sq", "v_ed", "v_eq")),
        StateGroup("dc_link", "DC link", ("vdc", "idc"), ("v_dc", "i_dc")),
        StateGroup("current_loop", "Inner current loop", ("M_CLd", "M_CLq"), ("M_CLd", "M_CLq"),
                   requires=("filter",), category=CONTROL),
        StateGroup("voltage_loop", "Outer voltage loop", ("M_VLd", "M_VLq"), ("M_VLd", "M_VLq"),
                   requires=("filter",), category=CONTROL),
        StateGroup(
            "angle", "Control angle", ("theta",), ("theta",), allowed=(DYNAMIC,),
            note="theta defines the converter's own frame and can never be removed",
            category=SYNCHRONISATION,
        ),
    ),
    levels={
        "full": {},
        "no_trafo": {"trafo_current": ALGEBRAIC},
        "no_filter": {"trafo_current": ALGEBRAIC, "filter": ALGEBRAIC},
        "no_inner": {"trafo_current": ALGEBRAIC, "filter": ALGEBRAIC, "current_loop": ALGEBRAIC},
        "no_voltage": {
            "trafo_current": ALGEBRAIC, "filter": ALGEBRAIC,
            "current_loop": ALGEBRAIC, "voltage_loop": ALGEBRAIC,
        },
        "droop": {
            "trafo_current": ALGEBRAIC, "filter": ALGEBRAIC, "current_loop": ALGEBRAIC,
            "voltage_loop": ALGEBRAIC, "dc_link": ALGEBRAIC,
        },
    },
    level_labels={
        "full": "Full order (EMT)",
        "no_trafo": "No transformer current",
        "no_filter": "No LC filter",
        "no_inner": "No inner current loop",
        "no_voltage": "No cascaded voltage control",
        "droop": "Droop only (RMS)",
    },
    level_notes={
        "full": "Everything dynamic. Needs a dynamic network.",
        "no_trafo": "Transformer current algebraic, so it pairs with either network.",
        "no_filter": "Ideal LC filter -- the converter's terminal voltage follows its reference.",
        "no_inner": "The current loop tracks its reference instantly.",
        "no_voltage": "Both cascaded loops ideal: a controlled voltage source behind the filter impedance.",
        "droop": "The droop law and its power measurement filters only -- three states. "
                 "This is the grid-forming model RMS studies use.",
    },
)

# --- Grid-following converter -------------------------------------------------
# The 14 states of components/gfl.py. The PLL is the synchronising dynamic
# and survives every level: the smallest model is a current source behind
# its PLL, which is what phasor tools represent a grid-following unit as.
_GFL = ElementModel(
    kind="gfl",
    label="Grid-following converter",
    groups=(
        StateGroup(
            "trafo_current", "Transformer current", ("igd", "igq"), ("i_gd", "i_gq"),
            note="must be algebraic when the network is quasi-stationary",
        ),
        StateGroup("filter", "LCL filter", ("isd", "isq", "ved", "veq"), ("i_sd", "i_sq", "v_ed", "v_eq")),
        StateGroup("dc_link", "DC link", ("vdc", "idc"), ("v_dc", "i_dc")),
        StateGroup("current_loop", "Inner current loop", ("M_CLd", "M_CLq"), ("M_CLd", "M_CLq"),
                   requires=("filter",), category=CONTROL),
        StateGroup(
            "outer_loop", "Outer DC-voltage and reactive-power loops", ("M_d", "M_q"), ("M_d", "M_q"),
            requires=("filter", "dc_link"), category=CONTROL,
        ),
        StateGroup(
            "pll", "Phase-locked loop", ("M_pll", "theta_pll"), ("M_pll", "theta_pll"),
            allowed=(DYNAMIC,), category=SYNCHRONISATION,
            note="the PLL is how a grid-following unit synchronises; theta_pll defines its frame",
        ),
    ),
    levels={
        "full": {},
        "no_trafo": {"trafo_current": ALGEBRAIC},
        "no_filter": {"trafo_current": ALGEBRAIC, "filter": ALGEBRAIC},
        "no_inner": {"trafo_current": ALGEBRAIC, "filter": ALGEBRAIC, "current_loop": ALGEBRAIC},
        "no_dc": {
            "trafo_current": ALGEBRAIC, "filter": ALGEBRAIC,
            "current_loop": ALGEBRAIC, "dc_link": ALGEBRAIC,
        },
        "pll": {
            "trafo_current": ALGEBRAIC, "filter": ALGEBRAIC, "current_loop": ALGEBRAIC,
            "dc_link": ALGEBRAIC, "outer_loop": ALGEBRAIC,
        },
    },
    level_labels={
        "full": "Full order (EMT)",
        "no_trafo": "No transformer current",
        "no_filter": "No LCL filter",
        "no_inner": "No inner current loop",
        "no_dc": "No DC link",
        "pll": "PLL only (RMS)",
    },
    level_notes={
        "full": "Everything dynamic. Needs a dynamic network.",
        "no_trafo": "Transformer current algebraic, so it pairs with either network.",
        "no_filter": "Ideal LCL filter.",
        "no_inner": "The current loop tracks its reference instantly. Watch this one in weak grids: "
                    "at low SCR or low switching frequency the current loop couples with the PLL "
                    "and dropping it changes the answer.",
        "no_dc": "DC link held at its reference.",
        "pll": "A controlled current source behind its PLL -- two states. This is the "
               "grid-following model RMS studies use.",
    },
)

# --- Network ------------------------------------------------------------------
# One pseudo-element covering every passive element, because the choice is
# almost always made for the network as a whole. Its groups are the element
# kinds, so a user who does want a mixed network (dynamic lines, algebraic
# loads) can still say so.
_NETWORK = ElementModel(
    kind="network",
    label="Network",
    groups=(
        StateGroup("nodes", "Bus voltages (shunt capacitance)", ("vgd_g", "vgq_g"), ("v_{g_d}", "v_{g_q}"),
                   category=NETWORK),
        StateGroup("lines", "Line currents", ("ild_g", "ilq_g"), ("i_{l_d}", "i_{l_q}"), category=NETWORK),
        StateGroup("loads", "Load currents", ("icd_g", "icq_g"), ("i_{c_d}", "i_{c_q}"), category=NETWORK),
        StateGroup("shunts", "Shunt reactor currents", ("ild_g", "ilq_g"), ("i_{l_d}", "i_{l_q}"),
                   category=NETWORK),
        StateGroup(
            "transformers", "Branch transformer currents", ("ild_g", "ilq_g"), ("i_{l_d}", "i_{l_q}"),
            category=NETWORK,
        ),
    ),
    levels={
        "full": {},
        "quasi_stationary": {
            "nodes": ALGEBRAIC, "lines": ALGEBRAIC, "loads": ALGEBRAIC,
            "shunts": ALGEBRAIC, "transformers": ALGEBRAIC,
        },
    },
    level_labels={"full": "Dynamic (EMT)", "quasi_stationary": "Quasi-stationary (RMS)"},
    level_notes={
        "full": "Every passive element integrates its own L di/dt and C dv/dt. "
                "This is what makes a run electromagnetic-transient.",
        "quasi_stationary": "Every passive element is algebraic: the interconnection's own "
                            "elimination becomes the admittance-matrix solve a phasor tool does. "
                            "Removes the network's fast modes, which are the stiffest in the model.",
    },
)

# --- the grid-forming outer laws ----------------------------------------------
# What each of symGFM_types.m's control laws carries beyond the shared
# cascade. The group ids are shared where the physics is: every law that
# filters its power measurement calls that group "power_filter", so a saved
# setting survives swapping one law for another that also has one.
_POWER_FILTER = StateGroup(
    "power_filter", "Power measurement filters", ("pm", "qm"), ("p_m", "q_m"),
    note="algebraic = the control law sees instantaneous power",
    category=SYNCHRONISATION,
)
# One state, one group id, whichever law owns it: `Dw` is the frequency
# deviation the angle integrates in both the filtered droop and the VSM, so
# it keeps one id and changes only its label and what may be done to it.
# That is what keeps a saved reduction setting meaningful across a swap, and
# what keeps looking a state up by name single-valued.
_DROOP_FILTER = StateGroup(
    "frequency", "Droop filter", ("Dw",), ("dw",), allowed=(DYNAMIC, ALGEBRAIC),
    note="algebraic = the droop acts instantly again, i.e. plain droop",
    category=SYNCHRONISATION,
)
GFM_OUTER_GROUPS: dict[str, tuple[StateGroup, ...]] = {
    "droop": (_POWER_FILTER,),
    "droop_filtered": (_POWER_FILTER, _DROOP_FILTER),
    "dvoc": (
        _POWER_FILTER,
        StateGroup(
            "voltage_dynamics", "dVOC voltage amplitude", ("ved_ref",), ("v_ed_ref",),
            allowed=(DYNAMIC, ALGEBRAIC),
            note="the amplitude the oscillator carries; algebraic = it settles instantly",
            category=SYNCHRONISATION,
        ),
    ),
    # A virtual synchronous machine measures power directly -- its inertia is
    # what smooths the response -- so it has no power filters to reduce.
    "vsm": (
        StateGroup(
            "frequency", "Virtual swing equation", ("Dw",), ("dw",), allowed=(DYNAMIC,),
            note="the emulated rotor; removing it would leave no control law at all",
            category=SYNCHRONISATION,
        ),
        StateGroup(
            "flux", "Virtual flux", ("Phi",), ("Phi",), allowed=(DYNAMIC, ALGEBRAIC),
            note="the emulated field; its product with speed is the voltage reference",
            category=SYNCHRONISATION,
        ),
    ),
    # Matching control takes its frequency from the DC link, so the only
    # outer state is the measurement of that voltage.
    "matching": (
        StateGroup(
            "dc_measurement", "DC voltage measurement", ("vdc_m",), ("v_dc_m",),
            allowed=(DYNAMIC, ALGEBRAIC),
            note="the filtered DC-link voltage, which *is* this law's frequency signal",
            category=SYNCHRONISATION,
        ),
    ),
}


@lru_cache(maxsize=None)
def gfm_element(controller: str = "droop") -> ElementModel:
    """The catalogue for a converter running one particular control law.

    Only the outer-loop groups differ; the filter, DC link and the two
    cascaded loops are the same whichever law is running, which is why the
    named levels below mean the same thing for all of them.
    """
    if controller not in GFM_OUTER_GROUPS:
        raise ValueError(f"unknown GFM controller {controller!r} -- have {sorted(GFM_OUTER_GROUPS)}")
    groups = tuple(_GFM.groups[:-1]) + GFM_OUTER_GROUPS[controller] + (_GFM.groups[-1],)
    labels = dict(_GFM.level_labels)
    # The lowest level is "whatever the outer law is, on its own".
    labels["droop"] = {
        "droop": "Droop only (RMS)",
        "droop_filtered": "Filtered droop only (RMS)",
        "dvoc": "dVOC only (RMS)",
        "vsm": "Virtual machine only (RMS)",
        "matching": "Matching only (RMS)",
    }[controller]
    return ElementModel(
        kind="gfm", label=_GFM.label, groups=groups, levels=_GFM.levels,
        level_labels=labels, level_notes=_GFM.level_notes,
    )


ELEMENTS: dict[str, ElementModel] = {e.kind: e for e in (_NETWORK, _SM, gfm_element(), _GFL)}


# Block-name prefix -> (catalogue kind, the one group it can be). The passive
# elements all reuse the same two-state models, so "i_{l_d}" alone cannot say
# whether it is a line, a shunt or a branch transformer -- the block name is
# what distinguishes them. A frame or an infinite bus has no reducible group
# of its own; both carry the reference angle, so they count as synchronisation.
_BLOCK_PREFIX: dict[str, tuple[str, str | None]] = {
    "SM": ("sm", None), "GFM": ("gfm", None), "GFL": ("gfl", None),
    "Nd": ("network", "nodes"), "Ln": ("network", "lines"), "Ld": ("network", "loads"),
    "Sh": ("network", "shunts"), "Tr": ("network", "transformers"),
}
# A reference frame is nothing but its angle. An infinite bus carries an
# angle too, but also its own source currents, which are electrical.
_UNGROUPED_CATEGORY = {"Frame": SYNCHRONISATION, "IB": UNIT_ELECTRICAL}
_ANGLE_STATES = {"theta", "theta_pll"}


def split_state_name(state_name: str) -> tuple[str, str]:
    """``"psi_d_{SM_1}"`` -> ``("psi_d", "SM_1")``.

    The split is on the *last* ``_{``: a display name can itself contain
    braces, as ``"v_{g_d}_{Nd_3}"`` does.
    """
    idx = state_name.rfind("_{")
    if idx < 0:
        return state_name, ""
    return state_name[:idx], state_name[idx + 2:].rstrip("}")


def group_of_state(state_name: str) -> str | None:
    """The reduction group one assembled state belongs to, as
    ``"<kind>.<group>"``, or None when it belongs to none (a reference frame
    or an infinite bus)."""
    base, block = split_state_name(state_name)
    if not block:
        return None
    kind, fixed = _BLOCK_PREFIX.get(block.split("_")[0], (None, None))
    if kind is None:
        return None
    if fixed is not None:
        return f"{kind}.{fixed}"
    group_id = _STATE_TO_GROUP.get((kind, base))
    return f"{kind}.{group_id}" if group_id else None


def _group_category(kind: str, group_id: str) -> str:
    """A group's category, found across regulator variants as well."""
    for variant in (*PSS_GROUPS.values(), *AVR_GROUPS.values(), *GOVERNOR_GROUPS.values()):
        if kind == "sm" and variant is not None and variant.id == group_id:
            return variant.category
    if kind == "gfm":
        for groups in GFM_OUTER_GROUPS.values():
            for variant in groups:
                if variant.id == group_id:
                    return variant.category
    return element(kind).group(group_id).category


def category_of_state(state_name: str) -> str:
    """Which of :data:`CATEGORIES` one assembled state belongs to."""
    group = group_of_state(state_name)
    if group is None:
        base, block = split_state_name(state_name)
        if base in _ANGLE_STATES:
            return SYNCHRONISATION
        return _UNGROUPED_CATEGORY.get(block.split("_")[0], MIXED)
    kind, group_id = group.split(".", 1)
    return _group_category(kind, group_id)


# Which group a state display name belongs to, across *every* regulator
# variant -- a state name is looked up without knowing which unit produced
# it, and `e_tgr` only exists on machines with the Kundur exciter. The group
# ids are shared between variants, so this stays single-valued.
_STATE_TO_GROUP: dict[tuple[str, str], str] = {}
for _kind, _e in ELEMENTS.items():
    for _g in _e.groups:
        for _name in _g.states:
            _STATE_TO_GROUP[(_kind, _name)] = _g.id
for _variant in (*PSS_GROUPS.values(), *AVR_GROUPS.values(), *GOVERNOR_GROUPS.values()):
    if _variant is None:
        continue
    for _name in _variant.states:
        _STATE_TO_GROUP[("sm", _name)] = _variant.id
for _groups in GFM_OUTER_GROUPS.values():
    for _variant in _groups:
        for _name in _variant.states:
            # A state name has to mean one group. Two laws may each own a
            # state called `dw`, but if they filed it under different group
            # ids this lookup would answer with whichever was registered
            # last -- silently, and differently depending on import order.
            _seen = _STATE_TO_GROUP.get(("gfm", _name))
            if _seen is not None and _seen != _variant.id:
                raise AssertionError(
                    f"gfm state {_name!r} is in group {_seen!r} for one control law and "
                    f"{_variant.id!r} for another; give it one id in both"
                )
            _STATE_TO_GROUP[("gfm", _name)] = _variant.id


def element(
    kind: str, *, exciter: str | None = None, pss: str | None = None,
    governor: str | None = None, controller: str | None = None,
) -> ElementModel:
    """The catalogue for one element type.

    ``exciter``, ``pss`` and ``governor`` pick a synchronous machine's
    regulators, and ``controller`` a grid-forming converter's power-control
    law; their state groups depend on the models chosen. Callers that hold a
    unit pass them; callers asking about the element *type* leave them out
    and get the defaults.
    """
    if kind == "sm" and (exciter or pss or governor):
        return sm_element(exciter or "g2elin", pss or "g2elin", governor or "g2elin")
    if kind == "gfm" and controller:
        return gfm_element(controller)
    try:
        return ELEMENTS[kind]
    except KeyError:
        raise ValueError(f"no reduction catalogue for {kind!r} -- have {sorted(ELEMENTS)}") from None


def resolve_modes(kind: str, level: str | None = None, overrides: dict[str, str] | None = None) -> dict[str, str]:
    """``{symbol name: mode}`` for one element type, ready for
    :func:`g2elin_core.components.base.apply_reduction`."""
    return element(kind).modes(level, overrides)


def level_ids(kind: str) -> list[str]:
    return list(element(kind).levels)


def group_ids(kind: str) -> list[str]:
    return [g.id for g in element(kind).groups]


def _group_info(g: StateGroup) -> dict:
    return {
        "id": g.id,
        "label": g.label,
        "states": list(g.states),
        "allowed": list(g.allowed),
        "default": g.default,
        "locked": g.locked,
        "requires": list(g.requires),
        "note": g.note,
    }


#: Regulator models a synchronous machine can carry, and the state group each
#: contributes. The UI substitutes the chosen one into the machine's group
#: list, so switching an exciter changes what the per-group controls show
#: without the frontend knowing anything about either model.
SM_REGULATORS = {
    "governor": {
        "label": "Governor",
        "default": "g2elin",
        "options": [
            {"id": "g2elin", "label": "Droop into a first-order lag"},
            {"id": "none", "label": "None (constant mechanical power)"},
        ],
        "groups": GOVERNOR_GROUPS,
    },
    "exciter": {
        "label": "Exciter / AVR",
        "default": "g2elin",
        "options": [
            {"id": "g2elin", "label": "G2ELin original (transducer, amplifier, exciter, rate feedback)"},
            {"id": "kundur", "label": "Thyristor with TGR (Kundur Fig. E12.9)"},
        ],
        "groups": AVR_GROUPS,
    },
    "pss": {
        "label": "Power system stabiliser",
        "default": "g2elin",
        "options": [
            {"id": "g2elin", "label": "G2ELin original (input filter, washout, two lead-lags)"},
            {"id": "kundur", "label": "Washout and two lead-lags (Kundur Fig. E12.9)"},
            {"id": "none", "label": "None (no stabiliser fitted)"},
        ],
        "groups": PSS_GROUPS,
    },
}


#: A converter's outer power-control laws, and the state groups each brings.
#: Shaped like SM_REGULATORS so the UI can drive both from one code path.
GFM_CONTROLLERS = {
    "controller": {
        "label": "Power control law",
        "default": "droop",
        "options": [
            {"id": "droop", "label": "Droop"},
            {"id": "droop_filtered", "label": "Droop behind a filter (2nd-order response)"},
            {"id": "dvoc", "label": "dVOC (dispatchable virtual oscillator)"},
            {"id": "vsm", "label": "VSM (virtual synchronous machine)"},
            {"id": "matching", "label": "Matching (DC voltage sets frequency)"},
        ],
        "groups": GFM_OUTER_GROUPS,
    },
}


def describe(kind: str) -> dict:
    """The catalogue for one element type as plain JSON-able data -- what
    the API hands the web UI to build its pickers from, so the UI never
    hard-codes a level name or a state group.

    A machine also carries its regulator choices, each with the state group
    it contributes, so the UI can redraw the AVR and PSS rows when one is
    swapped.
    """
    e = element(kind)
    out = {
        "kind": e.kind,
        "label": e.label,
        "default_level": e.default_level,
        "levels": [
            {
                "id": level_id,
                "label": e.level_labels.get(level_id, level_id),
                "note": e.level_notes.get(level_id, ""),
                "modes": e.modes_by_group(level_id),
            }
            for level_id in e.levels
        ],
        "groups": [_group_info(g) for g in e.groups],
    }
    if kind == "gfm":
        out["regulators"] = [
            {
                "id": slot,
                "label": spec["label"],
                "default": spec["default"],
                # A law brings several groups where a regulator brings one,
                # so the UI is handed the whole set it should swap in.
                "options": [
                    {**opt, "group": None,
                     "groups": [_group_info(g) for g in spec["groups"][opt["id"]]]}
                    for opt in spec["options"]
                ],
            }
            for slot, spec in GFM_CONTROLLERS.items()
        ]
    if kind == "sm":
        out["regulators"] = [
            {
                "id": slot,
                "label": spec["label"],
                "default": spec["default"],
                "options": [
                    {**opt, "group": _group_info(spec["groups"][opt["id"]])}
                    for opt in spec["options"] if spec["groups"][opt["id"]] is not None
                ] + [
                    # A model with no states of its own still has to be
                    # offered; it just has no group to show under it.
                    {**opt, "group": None}
                    for opt in spec["options"] if spec["groups"][opt["id"]] is None
                ],
            }
            for slot, spec in SM_REGULATORS.items()
        ]
    return out
