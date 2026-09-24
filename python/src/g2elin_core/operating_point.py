"""Bridges a solved :class:`~g2elin_core.powerflow.PowerFlowResult` to the
per-component operating points ``sym*_subs``-equivalents need.

Ports ``script_generic.m``'s "SM Initializations" section (computing
``SM.Pt``, ``SM.Et``, flux linkages, etc. from ``bus_sol``) plus the
node/line/load operating-point loops in its "Substituting in the individual
state spaces" section.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from g2elin_core.components.gfl import GflOperatingPoint
from g2elin_core.components.gfm import GfmOperatingPoint
from g2elin_core.components.sm import SmOperatingPoint
from g2elin_core.network.breakers import (
    branch_transformer_indices, shunt_is_reactor, shunt_reactor_rx, transformer_ratio, transformer_rx,
    unit_transformers,
)
from g2elin_core.network.schema import Network
from g2elin_core.network.validation import validate_network
from g2elin_core.powerflow.pandapower_adapter import PowerFlowResult
from g2elin_core.pu_base import pu_base

# Default synchronous-machine electrical/control parameters, transcribed from
# script_generic.m Section V ("Defining Default Tuning") — the same values
# WSCC/script_WSCC.m and CIGRE/script_CIGRE_Islanded.m use (not overridden
# per case in the MATLAB source).
DEFAULT_SM_ELECTRICAL_PARAMS = dict(
    Ra=0.003, Ll=0.15, Lad=1.66, Laq=1.61, Lfd=0.165, Rfd=0.0006,
    L1d=0.1713, R1d=0.0284, L1q=0.7252, R1q=0.00619, L2q=0.125, R2q=0.02368,
)
DEFAULT_SM_MECHANICAL_PARAMS = dict(H=5.0, KD=0.0)
DEFAULT_SM_GOVERNOR_PARAMS = dict(mp=0.005, TG=0.2)  # mp already /100 (0.5%)
DEFAULT_SM_AVR_PARAMS = dict(Tr=20e-3, Ka=300.0, Ta=0.001, Ke=1.0, Te=0.0001, Kfd=0.001, Tfd=0.1)
DEFAULT_SM_PSS_PARAMS = dict(T_LP=0.03, K_PSS=2.0, T_HP=2.0, T1n=0.05, T1d=0.02, T2n=3.0, T2d=5.4)

# Kundur Fig. E12.9's regulators, each under the figure's own parameter names
# (components/sm.py). A machine carries one exciter and one stabiliser, so
# only the chosen models' keys are in its parameter dict -- nothing is shared
# with the originals and nothing is renamed to a common spelling.
#
# TA = TB makes the transient gain reduction transparent, leaving a plain
# high-gain thyristor exciter; TB > TA is what actually reduces the gain the
# regulator shows across the electromechanical band.
# Both sets are the book's own, from the data for Kundur Example 12.9:
#   (iii) thyristor exciter with TGR:  KA 200, TR 0.01, TA 1.0, TB 10.0
#   (iv)  thyristor exciter with high transient gain and PSS:
#         KA 200, TR 0.01, KSTAB 20.0, TW 10.0, T1 0.05, T2 0.02, T3 3.0, T4 5.4
# The book pairs the stabiliser with *high transient gain*, i.e. with TGR
# switched off (TA = TB): (iii) and (iv) are alternatives, not a combination.
# The TGR values here are case (iii), so a machine taking both defaults has
# TGR on; set TA = TB for case (iv).
KUNDUR_AVR_PARAMS = dict(TR=0.01, KA=200.0, TA=1.0, TB=10.0)
KUNDUR_PSS_PARAMS = dict(KSTAB=20.0, TW=10.0, T1=0.05, T2=0.02, T3=3.0, T4=5.4)

SM_EXCITER_PARAMS = {"g2elin": DEFAULT_SM_AVR_PARAMS, "kundur": KUNDUR_AVR_PARAMS}
SM_PSS_PARAMS = {"g2elin": DEFAULT_SM_PSS_PARAMS, "kundur": KUNDUR_PSS_PARAMS, "none": {}}
SM_GOVERNOR_PARAMS = {"g2elin": DEFAULT_SM_GOVERNOR_PARAMS, "none": {}}
DEFAULT_SM_AUX_LOAD_MVA_FRACTION = 1e-3  # SM.PL = 100 kW at Sb = 100 MVA -> RL_pu = Sb/PL = 1000


def sm_params(
    *, sn_mva: float, f_hz: float, rt_pu: float, lt_pu: float,
    exciter: str = "g2elin", pss: str = "g2elin", governor: str = "g2elin",
) -> dict:
    """Full parameter dict for :func:`g2elin_core.components.sm.linearize_sm`.

    The regulator half depends on which models the machine carries: each one
    brings its own parameters under its own names, so what a user edits is
    the model's own data sheet rather than a shared set of near-equivalents.
    A machine's key set therefore depends on ``exciter`` and ``pss``, which
    is why everything that asks what a unit may override has to be told them
    -- see :func:`overridable_param_keys`.
    """
    if exciter not in SM_EXCITER_PARAMS:
        raise ValueError(f"unknown exciter {exciter!r} -- have {sorted(SM_EXCITER_PARAMS)}")
    if pss not in SM_PSS_PARAMS:
        raise ValueError(f"unknown PSS {pss!r} -- have {sorted(SM_PSS_PARAMS)}")
    if governor not in SM_GOVERNOR_PARAMS:
        raise ValueError(f"unknown governor {governor!r} -- have {sorted(SM_GOVERNOR_PARAMS)}")
    rl_pu = 1.0 / DEFAULT_SM_AUX_LOAD_MVA_FRACTION  # = Sb / PL, dimensionless regardless of Sb
    return {
        "wb": 2 * math.pi * f_hz,
        "Rt": rt_pu,
        "Lt": lt_pu,
        "RL_pu": rl_pu,
        **DEFAULT_SM_ELECTRICAL_PARAMS,
        **DEFAULT_SM_MECHANICAL_PARAMS,
        **SM_GOVERNOR_PARAMS[governor],
        **SM_EXCITER_PARAMS[exciter],
        **SM_PSS_PARAMS[pss],
    }


# GFM (Droop) defaults, transcribed from script_generic.m Section V.
DEFAULT_GFM_FILTER_PARAMS = dict(R1_pu=0.01, L1_pu=0.1, C_pu=0.1)
_GFM_DC_LINK_R_OHM = 500e3
_GFM_DC_LINK_C_F = 50e-3
_GFM_DC_LINK_T_DC = 1e-3
_H_FIRST_ORDER = 3.0
_KD_OPT = 200.0
_DROOP_NQ = 0.0001
_GFM_INNER_CURRENT_LOOP_TR_S = 0.1e-3  # script_generic.m's Inner_current_loop.tr


def gfm_params(
    *, sn_mva: float, f_hz: float, un_kv: float, rt_pu: float, lt_pu: float,
    tr_cl: float = _GFM_INNER_CURRENT_LOOP_TR_S,
) -> dict:
    """Full parameter dict for :func:`g2elin_core.components.gfm.linearize_gfm`.

    ``tr_cl`` is the inner current loop's closed-loop rise-time tuning
    target (``script_generic.m``'s ``Inner_current_loop.tr``) -- the pole-
    placement formula below (``wn_cl = 3/(zita*tr_cl)``, then ``KpCL``/
    ``KiCL`` from it) is transcribed unchanged from that file. Exposed as a
    parameter (default: the same fixed value every preset used before this)
    so a caller can retune just this one loop, e.g. to trace how it moves
    the closed-loop eigenvalues -- see ``notebooks/tour.ipynb``.
    """
    wb_val = 2 * math.pi * f_hz
    base = pu_base(wb=wb_val, un_kv=un_kv, sn_mva=sn_mva)

    dc_r_pu = _GFM_DC_LINK_R_OHM / base.zb_dc
    dc_g_pu = 1.0 / dc_r_pu
    dc_c_pu = _GFM_DC_LINK_C_F / base.cb_dc

    mp = 1.0 / _KD_OPT
    wf = 1.0 / (2 * mp * _H_FIRST_ORDER)

    zita = 0.707
    tau_dc = dc_c_pu / (wb_val * dc_g_pu)
    tr_dc = 5e-3
    wn_dc = 3 / (zita * tr_dc)
    kpdc = ((3 * tau_dc / tr_dc) - 1) / (1.0 / dc_g_pu)

    kv = wb_val / DEFAULT_GFM_FILTER_PARAMS["C_pu"]
    tr_vl = 15e-3
    wn_vl = 3 / (tr_vl * zita)
    kp_vl = 2 * zita * wn_vl / kv
    ki_vl = wn_vl**2 / kv

    ki_cl_gain = 1.0 / DEFAULT_GFM_FILTER_PARAMS["R1_pu"]
    tau_i = DEFAULT_GFM_FILTER_PARAMS["L1_pu"] / (wb_val * DEFAULT_GFM_FILTER_PARAMS["R1_pu"])
    wn_cl = 3 / (tr_cl * zita)
    kp_cl = (2 * zita * wn_cl * tau_i - 1) / ki_cl_gain
    ki_cl = (wn_cl**2 * tau_i) / ki_cl_gain

    return {
        "wb": wb_val, "wff": 1.0,
        "Rf": DEFAULT_GFM_FILTER_PARAMS["R1_pu"], "Lf": DEFAULT_GFM_FILTER_PARAMS["L1_pu"],
        "Cf": DEFAULT_GFM_FILTER_PARAMS["C_pu"], "Rt": rt_pu, "Lt": lt_pu,
        "mp": mp, "nq": _DROOP_NQ, "wf": wf,
        "KpVL": kp_vl, "KiVL": ki_vl, "Kffi": 1.0,
        "KpCL": kp_cl, "KiCL": ki_cl, "Kffv": 1.0,
        "Cdc": dc_c_pu, "Gdc": dc_g_pu, "Kpdc": kpdc, "Tdc": _GFM_DC_LINK_T_DC,
    }


# GFL defaults, transcribed from script_generic.m Section V.
DEFAULT_GFL_FILTER_PARAMS = dict(R1_pu=0.015, L1_pu=0.1, C_pu=0.11)
_GFL_DC_LINK_R_OHM = 500e3
_GFL_DC_LINK_C_F = 50e-3
_GFL_DC_LINK_T_DC = 1e-3


def gfl_params(*, sn_mva: float, f_hz: float, un_kv: float, rt_pu: float, lt_pu: float) -> dict:
    """Full parameter dict for :func:`g2elin_core.components.gfl.linearize_gfl`."""
    wb_val = 2 * math.pi * f_hz
    base = pu_base(wb=wb_val, un_kv=un_kv, sn_mva=sn_mva)

    dc_r_pu = _GFL_DC_LINK_R_OHM / base.zb_dc
    dc_g_pu = 1.0 / dc_r_pu
    dc_c_pu = _GFL_DC_LINK_C_F / base.cb_dc

    zita = 0.707
    kd_gain = -1.0 / dc_g_pu
    tau_d = dc_c_pu / (wb_val * dc_g_pu)
    tr_d = 100e-3
    wn_d = 3 / (tr_d * zita)
    vdc_kp = (2 * zita * wn_d * tau_d - 1) / kd_gain
    vdc_ki = (wn_d**2 * tau_d) / kd_gain

    tr_q = 100e-3
    kq = -3 / tr_q

    ki_cl_gain = 1.0 / DEFAULT_GFL_FILTER_PARAMS["R1_pu"]
    tau_i = DEFAULT_GFL_FILTER_PARAMS["L1_pu"] / (wb_val * DEFAULT_GFL_FILTER_PARAMS["R1_pu"])
    tr_cl = 10e-3
    wn_cl = 3 / (tr_cl * zita)
    kp_cl = (2 * zita * wn_cl * tau_i - 1) / ki_cl_gain
    ki_cl = (wn_cl**2 * tau_i) / ki_cl_gain

    tr_pll = 50e-3
    wn_pll = 3 / (tr_pll * zita)
    kp_pll = 2 * zita * wn_pll / wb_val
    ki_pll = wn_pll**2 / wb_val

    return {
        "wb": wb_val, "wff": 1.0,
        "Rf": DEFAULT_GFL_FILTER_PARAMS["R1_pu"], "Lf": DEFAULT_GFL_FILTER_PARAMS["L1_pu"],
        "Cf": DEFAULT_GFL_FILTER_PARAMS["C_pu"], "Rt": rt_pu, "Lt": lt_pu,
        "Kpd": vdc_kp, "Kid": vdc_ki, "Kiq": kq,
        "KpCL": kp_cl, "KiCL": ki_cl, "Kffv": 1.0,
        "Cdc": dc_c_pu, "Gdc": dc_g_pu, "Tdc": _GFL_DC_LINK_T_DC,
        "Kppll": kp_pll, "Kipll": ki_pll,
    }


# Parameters set by the network itself (the base frequency), never per unit.
NON_OVERRIDABLE_PARAMS = frozenset({"wb"})


def overridable_param_keys(
    unit_type: str, *, exciter: str | None = None, pss: str | None = None,
    governor: str | None = None,
) -> frozenset[str]:
    """Names a ``DerUnit.params`` override may use for ``unit_type`` (empty for
    an infinite bus, which has no parameter set of its own).

    The key set depends on nothing but the type and, for a machine, its
    regulator models -- each of which brings its own parameters under its own
    names -- so any base values do. Callers holding a unit should pass
    ``exciter``/``pss``, or use :func:`unit_keys`.
    """
    probe = dict(sn_mva=100.0, f_hz=50.0, rt_pu=0.0, lt_pu=0.1)
    if unit_type == "sm":
        keys = sm_params(**probe, exciter=exciter or "g2elin", pss=pss or "g2elin",
                         governor=governor or "g2elin")
    elif unit_type == "gfm":
        keys = gfm_params(**probe, un_kv=20.0)
    elif unit_type == "gfl":
        keys = gfl_params(**probe, un_kv=20.0)
    else:
        return frozenset()
    return frozenset(keys) - NON_OVERRIDABLE_PARAMS


def unit_keys(der) -> frozenset[str]:
    """Names ``der`` may override, for the regulator models it carries."""
    if der.unit_type.value == "sm":
        return overridable_param_keys(
            "sm", exciter=der.exciter_model, pss=der.pss_model, governor=der.governor_model
        )
    return overridable_param_keys(der.unit_type.value)


def unit_transformer_rx(network: Network, der) -> tuple[float, float]:
    """``(Rt, Lt)`` in system per unit for ``der``'s dynamic model.

    By default the unit's *own* transformer (the one whose LV side is the
    unit's bus) -- the same element the power flow uses -- converted from
    per unit of the transformer's rating to per unit of the network base
    (``Z_sys = Z_tr * Sn_sys / Sn_tr``), so power flow and dynamics describe
    the same impedance. With ``network.units_use_first_transformer`` (the
    MATLAB tool's convention, kept to reproduce its results), every unit uses
    the network's *first* transformer's values unconverted. Falls back to the
    first transformer when the unit has none of its own (such a network
    fails validation for modal/EMT anyway), and to placeholders when the
    network has no transformer at all.
    """
    own = unit_transformers(network).get(der.bus)
    if network.units_use_first_transformer or own is None:
        first = network.transformers[0] if network.transformers else None
        return (first.r_pu, first.x_pu) if first else (0.0, 0.05)
    k = network.sn_mva / own.sn_mva
    return own.r_pu * k, own.x_pu * k


def unit_params(network: Network, der) -> dict:
    """The parameter dict ``der``'s model uses in ``network``: its type's
    defaults (with the models' shared-transformer convention: the network's
    first transformer) plus ``der.params`` overrides. Empty for an infinite
    bus."""
    rt, lt = unit_transformer_rx(network, der)
    base = dict(sn_mva=network.sn_mva, f_hz=network.f_hz, rt_pu=rt, lt_pu=lt)
    kind = der.unit_type.value
    if kind == "sm":
        defaults = sm_params(**base, exciter=der.exciter_model, pss=der.pss_model,
                             governor=der.governor_model)
    elif kind in ("gfm", "gfl"):
        fn = gfm_params if kind == "gfm" else gfl_params
        defaults = fn(**base, un_kv=network.bus(der.bus).vn_kv)
    else:
        return {}
    return apply_param_overrides(der, defaults, network.sn_mva)


# --- Per-unit rebasing (DerUnit.sn_mva) --------------------------------------
# Published machine data is per unit of the *machine's* rating, not the
# network's. With Z_base = V^2/S: an impedance in pu of a 900 MVA machine is
# 9x smaller once referred to a 100 MVA network (x_d = 1.8 -> 0.2), while
# inertia goes the other way (H = 6.5 s -> 58.5 s), because H is energy per
# unit of rating. Capacitance and admittance follow inertia (Y_base = S/V^2).
#
# Only quantities whose base is unambiguous are converted. Controller gains
# and time constants are *tuning*, not machine data: a rise time is a rise
# time on any base, and a gain's base depends on which signals it sits
# between, so rebasing them by rule would be guesswork. validate_network()
# warns when a rebased unit overrides one of those instead.
_REBASE_AS_IMPEDANCE = frozenset({  # x Sn_network / Sn_unit
    "Ra", "Ll", "Lad", "Laq", "Lfd", "Rfd", "L1d", "R1d", "L1q", "R1q", "L2q", "R2q",
    "RL_pu", "Rf", "Lf", "Rt", "Lt", "mp", "nq",
})
_REBASE_AS_INERTIA = frozenset({  # x Sn_unit / Sn_network
    "H", "KD", "Cf", "Cdc", "Gdc",
})
REBASED_PARAM_KEYS = _REBASE_AS_IMPEDANCE | _REBASE_AS_INERTIA


def rebase_params(params: dict, *, from_mva: float, to_mva: float) -> dict:
    """``params`` given per unit of ``from_mva``, expressed per unit of ``to_mva``."""
    if from_mva == to_mva:
        return dict(params)
    k = to_mva / from_mva
    out = dict(params)
    for name, value in params.items():
        if name in _REBASE_AS_IMPEDANCE:
            out[name] = value * k
        elif name in _REBASE_AS_INERTIA:
            out[name] = value / k
    return out


def apply_param_overrides(der, params: dict, network_sn_mva: float | None = None) -> dict:
    """``params`` with ``der.params`` applied on top; raises on an override
    name this unit type doesn't have (validate_network() reports the same).

    When the unit declares its own rating (``DerUnit.sn_mva``) its overrides
    are read per unit of *that* rating and converted to the network's base
    first -- see :func:`rebase_params`. ``params`` itself (the type's
    defaults) is already on the network base and is left alone.
    """
    unknown = set(der.params) - unit_keys(der)
    if unknown:
        raise ValueError(f"unit id={der.id} ({der.unit_type.value}): unknown parameter override(s) {sorted(unknown)}")
    overrides = der.params
    if der.sn_mva is not None and network_sn_mva is not None:
        overrides = rebase_params(overrides, from_mva=der.sn_mva, to_mva=network_sn_mva)
    return {**params, **overrides}


def _rotate_to_global(v_pu: float, angle_rad: float, theta_g_rad: float) -> tuple[float, float]:
    z = v_pu * complex(math.cos(angle_rad), math.sin(angle_rad)) * complex(
        math.cos(-theta_g_rad), math.sin(-theta_g_rad)
    )
    return z.real, z.imag


@dataclass
class NetworkOperatingPoint:
    """Everything needed to linearize every component in a solved network."""

    theta_g_rad: float
    sm_ops: dict[int, SmOperatingPoint]  # DER id -> operating point
    gfm_ops: dict[int, GfmOperatingPoint]
    gfl_ops: dict[int, GflOperatingPoint]
    ib_ops: dict[int, dict]  # DER id -> kwargs for linearize_ib()/ib_nonlinear_point()
    node_vg: dict[int, tuple[float, float]]  # bus id -> (vgd_g0, vgq_g0)
    line_i0: list[tuple[float, float]]  # per Network.lines index -> (ild_g0, ilq_g0)
    load_rx: dict[int, tuple[float, float]]  # Load index -> (r_pu, x_pu)
    # Shunt *reactors* only (a capacitor bank is part of its bus's own
    # capacitance, not a branch): Network.shunts index -> (r_pu, x_pu) and
    # the current it draws at the operating point.
    shunt_rx: dict[int, tuple[float, float]] = field(default_factory=dict)
    shunt_i0: dict[int, tuple[float, float]] = field(default_factory=dict)
    # Branch transformers only (a unit's step-up is inside its unit's model):
    # Network.transformers index -> (r_pu, x_pu) on the network base, and the
    # current through its series impedance at the operating point.
    transformer_rx: dict[int, tuple[float, float]] = field(default_factory=dict)
    transformer_i0: dict[int, tuple[float, float]] = field(default_factory=dict)


def compute_operating_point(network: Network, result: PowerFlowResult) -> NetworkOperatingPoint:
    bus_table = result.bus_table().set_index("bus")

    # network.validation.validate_network() is the single, comprehensive
    # structural pre-flight check (every DER behind its own transformer,
    # that transformer landing on a genuine grid bus, no Line/Load on a
    # DER's own bus, network connectivity, etc.) -- every hand-crafted
    # preset in this codebase satisfies it by construction, but a network
    # assembled by hand (the web UI's network editor / drag-and-drop
    # builder) can easily violate any part of it, previously crashing deep
    # inside with a bare KeyError/IndexError at whichever line hit it
    # first. Collecting every "error"-severity issue here instead of
    # stopping at the first one lets a caller fix everything in one pass.
    # "warning"-severity issues (e.g. an unsupported slack unit type) are
    # left for interconnect/network_assembly.py's own NotImplementedError
    # (-> HTTP 501, not 422 -- "not supported yet" is a different kind of
    # problem than "this network is structurally broken").
    errors = [issue.message for issue in validate_network(network) if issue.severity == "error"]
    if errors:
        raise ValueError(
            f"{len(errors)} problem(s) with this network's structure:\n"
            + "\n".join(f"- {m}" for m in errors)
        )

    def bus_vm_va(bus_id: int) -> tuple[float, float]:
        row = bus_table.loc[bus_id]
        return float(row["vm_pu"]), math.radians(float(row["va_degree"]))

    def bus_pq(bus_id: int) -> tuple[float, float]:
        row = bus_table.loc[bus_id]
        return float(row["p_net_gen_mw"]), float(row["q_net_gen_mvar"])

    transformer_by_lv_bus = {tr.lv_bus: tr for tr in network.transformers}

    # DER units. Each unit's transformer impedance comes from
    # unit_transformer_rx(): its own transformer by default, or DG#1's for
    # every unit in the MATLAB-compatible mode (see the g2elin_core.components.sm
    # module docstring). The slack unit is linearized first: its own theta0 defines the global
    # reference angle theta_g_0 that every other unit needs (theta_g_rad is
    # a required constructor argument but is provably unused when
    # is_slack=True, so a placeholder there is harmless — see SmOperatingPoint).
    ders_by_slack_first = sorted(network.der_units, key=lambda d: d.bus_type.value != "slack")

    def build_sm_op(der, theta_g_rad: float) -> SmOperatingPoint:
        tr = transformer_by_lv_bus[der.bus]
        rt, lt = unit_transformer_rx(network, der)
        params = apply_param_overrides(
            der,
            sm_params(sn_mva=network.sn_mva, f_hz=network.f_hz, rt_pu=rt, lt_pu=lt,
                      exciter=der.exciter_model, pss=der.pss_model,
                      governor=der.governor_model),
            network.sn_mva,
        )
        v_t, a_t = bus_vm_va(der.bus)
        p_net, q_net = bus_pq(der.bus)
        p_gross_pu = (p_net + der.p_cons_mw) / network.sn_mva
        q_gross_pu = (q_net + der.q_cons_mvar) / network.sn_mva
        v_g, a_g = bus_vm_va(tr.hv_bus)
        return SmOperatingPoint(
            params=params,
            v_terminal_pu=v_t,
            angle_terminal_rad=a_t,
            p_terminal_pu=p_gross_pu,
            q_terminal_pu=q_gross_pu,
            v_grid_pu=v_g,
            angle_grid_rad=a_g,
            p_ref_pu=p_gross_pu,
            theta_g_rad=theta_g_rad,
            # A unit only *owns* the frame in the MATLAB-compatible mode; with
            # a frame of its own (components/frame.py) every unit is ordinary.
            is_slack=network.frame_follows_slack and der.bus_type.value == "slack",
            exciter=der.exciter_model,
            pss=der.pss_model,
            governor=der.governor_model,
        )

    # The reference angle theta_g_0: the slack machine's own rotor angle, as
    # in the toolbox. With a frame of its own the frame simply starts there,
    # so every unit's theta - theta_g is the value computed here either way.
    sm_ops: dict[int, SmOperatingPoint] = {}
    theta_g_rad = 0.0
    slack_id = next((d.id for d in network.der_units if d.bus_type.value == "slack"), None)
    for der in ders_by_slack_first:
        if der.unit_type.value != "sm":
            continue
        op = build_sm_op(der, theta_g_rad)
        sm_ops[der.id] = op
        if der.id == slack_id:
            theta_g_rad = op.theta0
            if not network.frame_follows_slack:
                # Rebuilt now that the frame angle is known (its own theta_g).
                sm_ops[der.id] = build_sm_op(der, theta_g_rad)

    # GFM/GFL are never the slack (network_form.m never places them there),
    # so theta_g_rad is already final by the time these run.
    gfm_ops: dict[int, GfmOperatingPoint] = {}
    gfl_ops: dict[int, GflOperatingPoint] = {}
    for der in network.der_units:
        if der.unit_type.value not in ("gfm", "gfl"):
            continue
        tr = transformer_by_lv_bus[der.bus]
        un_kv = network.bus(der.bus).vn_kv
        v_t, a_t = bus_vm_va(der.bus)
        p_net, q_net = bus_pq(der.bus)
        p_pu = (p_net + der.p_cons_mw) / network.sn_mva
        q_pu = (q_net + der.q_cons_mvar) / network.sn_mva
        v_g, a_g = bus_vm_va(tr.hv_bus)
        if der.unit_type.value == "gfm":
            rt, lt = unit_transformer_rx(network, der)
            params = apply_param_overrides(der, gfm_params(
                sn_mva=network.sn_mva, f_hz=network.f_hz, un_kv=un_kv, rt_pu=rt, lt_pu=lt,
            ), network.sn_mva)
            gfm_ops[der.id] = GfmOperatingPoint(
                params=params, v_terminal_pu=v_t, angle_terminal_rad=a_t,
                p_terminal_pu=p_pu, q_terminal_pu=q_pu,
                v_grid_pu=v_g, angle_grid_rad=a_g, theta_g_rad=theta_g_rad,
            )
        else:
            rt, lt = unit_transformer_rx(network, der)
            params = apply_param_overrides(der, gfl_params(
                sn_mva=network.sn_mva, f_hz=network.f_hz, un_kv=un_kv, rt_pu=rt, lt_pu=lt,
            ), network.sn_mva)
            gfl_ops[der.id] = GflOperatingPoint(
                params=params, v_terminal_pu=v_t, angle_terminal_rad=a_t,
                p_terminal_pu=p_pu, q_terminal_pu=q_pu,
                v_grid_pu=v_g, angle_grid_rad=a_g, theta_g_rad=theta_g_rad,
            )

    # Nodes: every raw (non-DER) bus. Faithfully replicates the "every node
    # uses line #1's susceptance" quirk (see sm.py module docstring).
    der_buses = {d.bus for d in network.der_units}
    node_vg: dict[int, tuple[float, float]] = {}
    for bus in network.buses:
        if bus.id in der_buses:
            continue
        v, a = bus_vm_va(bus.id)
        node_vg[bus.id] = _rotate_to_global(v, a, theta_g_rad)

    wb_val = 2 * math.pi * network.f_hz

    # IB (always the slack, see components/ib.py) -- unlike SM/GFM/GFL it
    # has no internal control-loop parameters to derive, only its transformer
    # impedance (unit_transformer_rx, like every other unit).
    ib_ops: dict[int, dict] = {}
    for der in network.der_units:
        if der.unit_type.value != "infinite_bus":
            continue
        v_t, a_t = bus_vm_va(der.bus)
        p_net, q_net = bus_pq(der.bus)
        rt, lt = unit_transformer_rx(network, der)
        own_frame = network.frame_follows_slack and der.id == slack_id
        ib_ops[der.id] = dict(
            wb_val=wb_val, r_pu=rt, x_pu=lt, v_pu=v_t,
            p_mw=p_net + der.p_cons_mw, q_mvar=q_net + der.q_cons_mvar, sn_mva=network.sn_mva,
            # Where this source sits in the common frame (0 when it is the frame).
            delta_rad=0.0 if own_frame else math.radians(a_t) - theta_g_rad,
            is_slack=own_frame,
        )

    # Lines: solve the steady-state 2x2 for the line current.
    line_i0: list[tuple[float, float]] = []
    for ln in network.lines:
        vj, aj = bus_vm_va(ln.from_bus)
        vk, ak = bus_vm_va(ln.to_bus)
        vgdj, vgqj = _rotate_to_global(vj, aj, theta_g_rad)
        vgdk, vgqk = _rotate_to_global(vk, ak, theta_g_rad)
        coeff = np.array([[-ln.r_pu, ln.x_pu], [-ln.x_pu, -ln.r_pu]])
        rhs = np.array([-(vgdj - vgdk), -(vgqj - vgqk)])
        ild0, ilq0 = np.linalg.solve(coeff, rhs)
        line_i0.append((float(ild0), float(ilq0)))

    # Branch transformers (between two grid buses, as opposed to a unit's own
    # step-up, whose impedance lives inside that unit's model): the same RL
    # branch a line uses, fed through the ideal transformer's ratio. With a
    # complex ratio a*exp(j*phi), the voltage driving the series impedance is
    # v_hv/(a*exp(j*phi)) -- in dq, a scaling by 1/a and a rotation by -phi.
    transformer_i0: dict[int, tuple[float, float]] = {}
    transformer_rx_pu: dict[int, tuple[float, float]] = {}
    for idx in branch_transformer_indices(network):
        tr = network.transformers[idx]
        r_pu, x_pu = transformer_rx(tr, network.sn_mva)
        a, phi = transformer_ratio(tr)
        transformer_rx_pu[idx] = (r_pu, x_pu)
        vj, aj = bus_vm_va(tr.hv_bus)
        vk, ak = bus_vm_va(tr.lv_bus)
        vgdj, vgqj = _rotate_to_global(vj, aj, theta_g_rad)
        vgdk, vgqk = _rotate_to_global(vk, ak, theta_g_rad)
        cos_p, sin_p = math.cos(phi), math.sin(phi)
        vpd = (cos_p * vgdj + sin_p * vgqj) / a
        vpq = (-sin_p * vgdj + cos_p * vgqj) / a
        coeff = np.array([[-r_pu, x_pu], [-x_pu, -r_pu]])
        itr = np.linalg.solve(coeff, np.array([-(vpd - vgdk), -(vpq - vgqk)]))
        transformer_i0[idx] = (float(itr[0]), float(itr[1]))

    # Shunt reactors: an RL branch from the bus to zero volts. X comes from
    # the reactive power it absorbs at 1 pu (Q = V^2/X), R from its X/R --
    # exactly 0 would leave its resonance with the bus capacitance undamped.
    shunt_rx: dict[int, tuple[float, float]] = {}
    shunt_i0: dict[int, tuple[float, float]] = {}
    for idx, sh in enumerate(network.shunts):
        if not shunt_is_reactor(sh):
            continue
        r_pu, x_pu = shunt_reactor_rx(sh, network.sn_mva)
        shunt_rx[idx] = (r_pu, x_pu)
        vgd, vgq = node_vg[sh.bus]
        coeff = np.array([[-r_pu, x_pu], [-x_pu, -r_pu]])
        ish = np.linalg.solve(coeff, np.array([-vgd, -vgq]))
        shunt_i0[idx] = (float(ish[0]), float(ish[1]))

    # Loads: constant-impedance equivalent at the operating voltage.
    # theta = acos(PF) (matching Functions/*.m's script_generic.m Load
    # Parameters section) only gives the correct sign of X for inductive
    # (Q >= 0) loads, which is all that's used here.
    load_rx: dict[int, tuple[float, float]] = {}
    for idx, load in enumerate(network.loads):
        p_pu, q_pu = load.p_mw / network.sn_mva, load.q_mvar / network.sn_mva
        if q_pu == 0.0:
            # Not just the p=q=0 case (undefined apparent power, division by
            # zero right below): components/load.py's own dynamic model uses
            # this load's reactance x_pu = z_pu*sin(acos(p_pu/s_pu)) as Lc in
            # (wb/Lc)*(...) -- a purely resistive load (q_mvar=0, any p_mw)
            # makes sin(acos(+-1))=0, so x_pu=0 and that division is zero too.
            raise ValueError(
                f"load #{idx} (bus {load.bus}) has zero reactive power (q_mvar=0) -- its constant-"
                "impedance equivalent reactance would be exactly zero (division by zero downstream); "
                "give it a small nonzero q_mvar (positive for inductive, negative for capacitive)"
            )
        s_pu = math.hypot(p_pu, q_pu)
        v_pu, _ = bus_vm_va(load.bus)
        z_pu = v_pu**2 / s_pu
        pf = p_pu / s_pu
        theta = math.acos(pf)
        load_rx[idx] = (z_pu * math.cos(theta), z_pu * math.sin(theta))

    return NetworkOperatingPoint(
        theta_g_rad=theta_g_rad, sm_ops=sm_ops, gfm_ops=gfm_ops, gfl_ops=gfl_ops, ib_ops=ib_ops,
        node_vg=node_vg, line_i0=line_i0, load_rx=load_rx, shunt_rx=shunt_rx, shunt_i0=shunt_i0,
        transformer_rx=transformer_rx_pu, transformer_i0=transformer_i0,
    )
