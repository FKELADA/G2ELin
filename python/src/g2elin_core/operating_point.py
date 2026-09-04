"""Bridges a solved :class:`~g2elin_core.powerflow.PowerFlowResult` to the
per-component operating points ``sym*_subs``-equivalents need.

Ports ``script_generic.m``'s "SM Initializations" section (computing
``SM.Pt``, ``SM.Et``, flux linkages, etc. from ``bus_sol``) plus the
node/line/load operating-point loops in its "Substituting in the individual
state spaces" section.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from g2elin_core.components.gfl import GflOperatingPoint
from g2elin_core.components.gfm import GfmOperatingPoint
from g2elin_core.components.sm import SmOperatingPoint
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
DEFAULT_SM_MECHANICAL_PARAMS = dict(H=5.0, KD=0.0, mp=0.005, TG=0.2)  # mp already /100 (0.5%)
DEFAULT_SM_AVR_PARAMS = dict(Tr=20e-3, Ka=300.0, Ta=0.001, Ke=1.0, Te=0.0001, Kfd=0.001, Tfd=0.1)
DEFAULT_SM_PSS_PARAMS = dict(T_LP=0.03, K_PSS=2.0, T_HP=2.0, T1n=0.05, T1d=0.02, T2n=3.0, T2d=5.4)
DEFAULT_SM_AUX_LOAD_MVA_FRACTION = 1e-3  # SM.PL = 100 kW at Sb = 100 MVA -> RL_pu = Sb/PL = 1000


def sm_params(*, sn_mva: float, f_hz: float, rt_pu: float, lt_pu: float) -> dict:
    """Full parameter dict for :func:`g2elin_core.components.sm.linearize_sm`."""
    rl_pu = 1.0 / DEFAULT_SM_AUX_LOAD_MVA_FRACTION  # = Sb / PL, dimensionless regardless of Sb
    return {
        "wb": 2 * math.pi * f_hz,
        "Rt": rt_pu,
        "Lt": lt_pu,
        "RL_pu": rl_pu,
        **DEFAULT_SM_ELECTRICAL_PARAMS,
        **DEFAULT_SM_MECHANICAL_PARAMS,
        **DEFAULT_SM_AVR_PARAMS,
        **DEFAULT_SM_PSS_PARAMS,
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

    # DER units, faithfully replicating the "every SM uses DG#1's transformer
    # impedance" quirk (see g2elin_core.components.sm module docstring). The
    # slack unit is linearized first: its own theta0 defines the global
    # reference angle theta_g_0 that every other unit needs (theta_g_rad is
    # a required constructor argument but is provably unused when
    # is_slack=True, so a placeholder there is harmless — see SmOperatingPoint).
    first_tr = network.transformers[0]
    ders_by_slack_first = sorted(network.der_units, key=lambda d: d.bus_type.value != "slack")

    def build_sm_op(der, theta_g_rad: float) -> SmOperatingPoint:
        tr = transformer_by_lv_bus[der.bus]
        params = sm_params(sn_mva=network.sn_mva, f_hz=network.f_hz, rt_pu=first_tr.r_pu, lt_pu=first_tr.x_pu)
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
            is_slack=der.bus_type.value == "slack",
        )

    sm_ops: dict[int, SmOperatingPoint] = {}
    theta_g_rad = 0.0
    for der in ders_by_slack_first:
        if der.unit_type.value != "sm":
            continue
        op = build_sm_op(der, theta_g_rad)
        sm_ops[der.id] = op
        if op.is_slack:
            theta_g_rad = op.theta0

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
            params = gfm_params(
                sn_mva=network.sn_mva, f_hz=network.f_hz, un_kv=un_kv,
                rt_pu=first_tr.r_pu, lt_pu=first_tr.x_pu,
            )
            gfm_ops[der.id] = GfmOperatingPoint(
                params=params, v_terminal_pu=v_t, angle_terminal_rad=a_t,
                p_terminal_pu=p_pu, q_terminal_pu=q_pu,
                v_grid_pu=v_g, angle_grid_rad=a_g, theta_g_rad=theta_g_rad,
            )
        else:
            params = gfl_params(
                sn_mva=network.sn_mva, f_hz=network.f_hz, un_kv=un_kv,
                rt_pu=first_tr.r_pu, lt_pu=first_tr.x_pu,
            )
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
    # has no internal control-loop parameters to derive, just the same
    # transformer-impedance ("DG#1's transformer") quirk applied uniformly.
    ib_ops: dict[int, dict] = {}
    for der in network.der_units:
        if der.unit_type.value != "infinite_bus":
            continue
        v_t, _ = bus_vm_va(der.bus)
        p_net, q_net = bus_pq(der.bus)
        ib_ops[der.id] = dict(
            wb_val=wb_val, r_pu=first_tr.r_pu, x_pu=first_tr.x_pu, v_pu=v_t,
            p_mw=p_net + der.p_cons_mw, q_mvar=q_net + der.q_cons_mvar, sn_mva=network.sn_mva,
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

    # Loads: constant-impedance equivalent at the operating voltage.
    # theta = acos(PF) (matching Functions/*.m's script_generic.m Load
    # Parameters section) only gives the correct sign of X for inductive
    # (Q >= 0) loads, which is all that's used here.
    load_rx: dict[int, tuple[float, float]] = {}
    for idx, load in enumerate(network.loads):
        p_pu, q_pu = load.p_mw / network.sn_mva, load.q_mvar / network.sn_mva
        s_pu = math.hypot(p_pu, q_pu)
        v_pu, _ = bus_vm_va(load.bus)
        z_pu = v_pu**2 / s_pu
        pf = p_pu / s_pu
        theta = math.acos(pf)
        load_rx[idx] = (z_pu * math.cos(theta), z_pu * math.sin(theta))

    return NetworkOperatingPoint(
        theta_g_rad=theta_g_rad, sm_ops=sm_ops, gfm_ops=gfm_ops, gfl_ops=gfl_ops, ib_ops=ib_ops,
        node_vg=node_vg, line_i0=line_i0, load_rx=load_rx
    )
