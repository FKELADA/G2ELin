"""Preset network cases, ported from the MATLAB toolbox's ``preset_networks.m``.

Each preset is transcribed directly from the corresponding MATLAB source so
it can be diffed against the original if the values there ever change.

The WSCC and CIGRE-islanded families each share one raw topology across
several DER-mix variants (``WSCC/script_WSCC.m``'s ``switch model.name`` and
``CIGRE/script_CIGRE_Islanded.m``'s ``switch model.name`` both only vary the
``Y_DER`` table between cases, never the underlying buses/lines/loads) — so
each family is a shared ``_*_topology()`` builder plus one thin function per
variant supplying just its own DER rows, rather than duplicating the shared
part 8 (WSCC) or 4 (CIGRE) times over.
"""

from __future__ import annotations

from typing import NamedTuple

from .schema import Bus, BusType, DerUnit, GfmController, Line, Load, Network, Transformer, UnitType

# Base values, from WSCC/script_WSCC.m Section II ("General network parameters"):
#   f_hz = 60, Sn = 100 MVA, LV = 18 kV (generator terminals), HV = 230 kV (transmission)
_SN_MVA = 100.0
_HV_KV = 230.0
_LV_KV = 18.0


def _wscc9_topology() -> tuple[list[Bus], list[Line], list[Load]]:
    """The WSCC 9-bus raw transmission network + loads, shared by every
    DER-mix variant below. Ported from ``Functions/WSCC_raw.m``.

    Bus numbering matches the MATLAB tool: 1-3 are the HV terminals the
    generator transformers tap into, 4-6 are the HV load buses.
    """
    buses = [Bus(id=i, name=f"bus{i}", vn_kv=_HV_KV) for i in range(1, 7)]

    # Functions/WSCC_raw.m, Line.R_pu / Line.L_pu / Line.B_pu (modif = 1)
    lines = [
        Line(from_bus=1, to_bus=4, r_pu=0.0100, x_pu=0.0850, b_pu=0.176),
        Line(from_bus=1, to_bus=6, r_pu=0.0170, x_pu=0.0920, b_pu=0.158),
        Line(from_bus=2, to_bus=4, r_pu=0.0320, x_pu=0.1610, b_pu=0.306),
        Line(from_bus=2, to_bus=5, r_pu=0.0085, x_pu=0.0720, b_pu=0.149),
        Line(from_bus=3, to_bus=5, r_pu=0.0119, x_pu=0.1008, b_pu=0.209),
        Line(from_bus=3, to_bus=6, r_pu=0.0390, x_pu=0.1700, b_pu=0.358),
    ]

    # Functions/WSCC_raw.m, Load.P / Load.Q (MW / MVAr, nodes 4-6 only)
    loads = [
        Load(bus=4, p_mw=125.0, q_mvar=50.0),
        Load(bus=5, p_mw=100.0, q_mvar=35.0),
        Load(bus=6, p_mw=90.0, q_mvar=30.0),
    ]
    return buses, lines, loads


class _WsccDer(NamedTuple):
    """One ``Y_DER`` row from ``WSCC/script_WSCC.m``: ``[where DG-type load
    bustype V_init delta Pgen Qgen Pcons Qcons TR_R TR_XL Xd]`` -- ``load``
    and ``delta`` columns are always 0 for a DER's own bus and dropped here.
    """

    dispatch_bus: int  # "where?" -- the raw WSCC node (1, 2, or 3) this DER connects to
    unit_type: UnitType
    bus_type: BusType
    v_set_pu: float
    p_set_mw: float
    q_set_mvar: float
    p_cons_mw: float
    tr_r_pu: float
    tr_xl_pu: float
    xd_pu: float | None = None
    controller: GfmController | None = None


def _wscc9_network(name: str, ders: list[_WsccDer]) -> Network:
    """Assembles a WSCC preset from its DER rows the way
    ``Functions/network_form.m`` does: each DER gets its own new bus
    (7, 8, 9, ... in the order given here), connected to its
    ``dispatch_bus`` through a step-up transformer using that row's own
    ``TR_R``/``TR_XL``.
    """
    topo_buses, lines, loads = _wscc9_topology()
    der_buses = [Bus(id=7 + i, name=f"gen{i + 1}_terminal", vn_kv=_LV_KV) for i in range(len(ders))]
    der_units = [
        DerUnit(
            id=i + 1,
            bus=7 + i,
            unit_type=d.unit_type,
            bus_type=d.bus_type,
            v_set_pu=d.v_set_pu,
            p_set_mw=d.p_set_mw,
            q_set_mvar=d.q_set_mvar,
            p_cons_mw=d.p_cons_mw,
            controller=d.controller,
            xd_pu=d.xd_pu,
        )
        for i, d in enumerate(ders)
    ]
    transformers = [
        Transformer(
            hv_bus=d.dispatch_bus, lv_bus=7 + i, r_pu=d.tr_r_pu, x_pu=d.tr_xl_pu,
            sn_mva=_SN_MVA, name=f"gen{i + 1}_xfmr",
        )
        for i, d in enumerate(ders)
    ]
    return Network(
        name=name, f_hz=60.0, sn_mva=_SN_MVA,
        buses=topo_buses + der_buses, lines=lines, transformers=transformers,
        loads=loads, der_units=der_units,
    )


def wscc9_3sm() -> Network:
    """The classic WSCC 9-bus system, 3 synchronous machines.

    ``WSCC/script_WSCC.m``, case ``'WSCC_3SM'``. Bus numbering matches the
    MATLAB tool exactly: 7-9 are the LV generator terminal buses (behind
    their transformers), in DER declaration order.
    """
    return _wscc9_network("WSCC_3SM", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(2, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0625),
        _WsccDer(3, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0586),
    ])


def wscc9_2sm_1gfl() -> Network:
    """WSCC 9-bus, 2 synchronous machines + 1 grid-following converter.

    ``WSCC/script_WSCC.m``, case ``'WSCC_2SM_1GFL'``.
    """
    return _wscc9_network("WSCC_2SM_1GFL", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(2, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0625),
        _WsccDer(3, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586),
    ])


def wscc9_1sm_2gfl() -> Network:
    """WSCC 9-bus, 1 synchronous machine + 2 grid-following converters.

    ``WSCC/script_WSCC.m``, case ``'WSCC_1SM_2GFL'``.
    """
    return _wscc9_network("WSCC_1SM_2GFL", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(2, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0625),
        _WsccDer(3, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586),
    ])


def wscc9_1sm_1gfm_1gfl() -> Network:
    """WSCC 9-bus, 1 synchronous machine + 1 grid-forming (Droop) + 1
    grid-following converter.

    ``WSCC/script_WSCC.m``, case ``'WSCC_1SM_1GFM_1GFL'``.
    """
    return _wscc9_network("WSCC_1SM_1GFM_1GFL", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(2, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0625, controller=GfmController.DROOP),
        _WsccDer(3, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586),
    ])


def wscc9_1sm_2gfm() -> Network:
    """WSCC 9-bus, 1 synchronous machine + 2 grid-forming (Droop) converters.

    ``WSCC/script_WSCC.m``, case ``'WSCC_1SM_2GFM'``.
    """
    return _wscc9_network("WSCC_1SM_2GFM", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(2, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0625, controller=GfmController.DROOP),
        _WsccDer(3, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586, controller=GfmController.DROOP),
    ])


def wscc9_2sm_1gfm() -> Network:
    """WSCC 9-bus, 2 synchronous machines + 1 grid-forming (Droop) converter.

    ``WSCC/script_WSCC.m``, case ``'WSCC_2SM_1GFM'``. Declared in the same
    order as the MATLAB source (slack SM at node 1, GFM at node 3, second SM
    at node 2) -- note the second SM row there carries a real ``Xd`` (0.2),
    unlike the other variants' second SM row, transcribed as given rather
    than "corrected" to match.
    """
    return _wscc9_network("WSCC_2SM_1GFM", [
        _WsccDer(1, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, xd_pu=0.2),
        _WsccDer(3, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586, controller=GfmController.DROOP),
        _WsccDer(2, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0625, xd_pu=0.2),
    ])


def wscc9_1gfm_2gfl() -> Network:
    """WSCC 9-bus, 1 grid-forming (Droop) converter + 2 grid-following
    converters.

    ``WSCC/script_WSCC.m``, case ``'WSCC_1GFM_2GFL'``.

    **Power flow only** -- the slack unit here is a GFM, and this project's
    interconnection wiring (``interconnect/network_assembly.py``) only
    supports a synchronous-machine or infinite-bus slack so far (a GFM
    slack behaves fundamentally differently -- it doesn't have a P/Q
    dispatch to hold, it defines the voltage/frequency reference the rest
    of the network locks to -- and needs real modeling work, not just
    wiring). Power flow works fine (pandapower's slack handling doesn't
    care about unit type), but modal analysis / EMT raise
    ``NotImplementedError`` (surfaced as an HTTP 501 by the API), the same
    documented, controlled failure mode as any other unsupported DER
    combination in this codebase.
    """
    return _wscc9_network("WSCC_1GFM_2GFL", [
        _WsccDer(1, UnitType.GFM, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, controller=GfmController.DROOP),
        _WsccDer(2, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0625),
        _WsccDer(3, UnitType.GFL, BusType.PQ, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586),
    ])


def wscc9_3gfm() -> Network:
    """WSCC 9-bus, 3 grid-forming (Droop) converters.

    ``WSCC/script_WSCC.m``, case ``'WSCC_3GFM'``. **Power flow only** --
    same GFM-slack limitation as :func:`wscc9_1gfm_2gfl`, see its
    docstring.
    """
    return _wscc9_network("WSCC_3GFM", [
        _WsccDer(1, UnitType.GFM, BusType.SLACK, 1.025, 0.8 * _SN_MVA, 0.0, 0.001 * _SN_MVA, 0.0, 0.0576, controller=GfmController.DROOP),
        _WsccDer(2, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0625, controller=GfmController.DROOP),
        _WsccDer(3, UnitType.GFM, BusType.PV, 1.0, 1.05 * _SN_MVA, 0.0, 0.0, 0.0, 0.0586, controller=GfmController.DROOP),
    ])


# Base values, from CIGRE/script_CIGRE_Islanded.m Section II:
#   f_hz = 50, Sn = 2.5 MVA, LV = 10 kV (DER terminals), MV = 20 kV (feeder)
_CIGRE_SN_MVA = 2.5
_CIGRE_MV_KV = 20.0
_CIGRE_LV_KV = 10.0
# Every DER's own step-up transformer uses this same impedance in every
# CIGRE-islanded variant (confirmed across all 4 Y_DER tables in
# CIGRE/script_CIGRE_Islanded.m -- TR_R/TR_XL columns are 0.0002/0.0240 on
# every row, in every case).
_CIGRE_DER_TR_R_PU = 0.0002
_CIGRE_DER_TR_XL_PU = 0.0240


def _cigre_raw_topology() -> tuple[list[Bus], list[Line], list[Load]]:
    """The CIGRE MV feeder's raw nodes/lines/loads, shared by every
    islanded DER-mix variant below. Ported from ``Functions/CIGRE_raw.m``
    with switches ``S0=0`` (islanded), ``S1=1`` (closed), ``S2=0``,
    ``S3=0`` (open) -- i.e. the tie line between raw nodes 8 and 12 is
    present, the other two tie lines are not.
    """
    buses = [Bus(id=i, name=f"bus{i}", vn_kv=_CIGRE_MV_KV) for i in range(1, 15)]

    # Functions/CIGRE_raw.m, islanded branch (S0=0), lines 1-12 unconditional
    # + line 13 (tie 8-12, S1=1). Values computed from Line.R/XL/B/length at
    # Zb_MV = 160 ohm, Lb_MV = 0.509296 H (Sb=2.5 MVA, Un_MV=20 kV, wb=2*pi*50).
    lines = [
        Line(from_bus=2, to_bus=1, r_pu=0.008830, x_pu=0.012619, b_pu=0.021429, length_km=2.82),
        Line(from_bus=2, to_bus=3, r_pu=0.013840, x_pu=0.019779, b_pu=0.033587, length_km=4.42),
        Line(from_bus=4, to_bus=3, r_pu=0.001910, x_pu=0.002730, b_pu=0.004635, length_km=0.61),
        Line(from_bus=4, to_bus=5, r_pu=0.001754, x_pu=0.002506, b_pu=0.004255, length_km=0.56),
        Line(from_bus=6, to_bus=5, r_pu=0.004822, x_pu=0.006891, b_pu=0.011702, length_km=1.54),
        Line(from_bus=8, to_bus=3, r_pu=0.004071, x_pu=0.005817, b_pu=0.009879, length_km=1.30),
        Line(from_bus=8, to_bus=7, r_pu=0.005229, x_pu=0.007473, b_pu=0.012690, length_km=1.67),
        Line(from_bus=8, to_bus=9, r_pu=0.000751, x_pu=0.001074, b_pu=0.001824, length_km=0.24),
        Line(from_bus=10, to_bus=9, r_pu=0.002411, x_pu=0.003446, b_pu=0.005851, length_km=0.77),
        Line(from_bus=10, to_bus=11, r_pu=0.001002, x_pu=0.001432, b_pu=0.002432, length_km=0.32),
        Line(from_bus=13, to_bus=14, r_pu=0.015587, x_pu=0.011186, b_pu=0.002482, length_km=4.89),
        Line(from_bus=13, to_bus=12, r_pu=0.009531, x_pu=0.006840, b_pu=0.001517, length_km=2.99),
        Line(from_bus=8, to_bus=12, r_pu=0.006375, x_pu=0.004575, b_pu=0.001015, length_km=2.00, name="tie_S1"),
    ]

    # Functions/CIGRE_raw.m, Load.P / Load.Q (MW / MVAr), node 2 has none
    loads = [
        Load(bus=1, p_mw=1.200, q_mvar=0.35000),
        Load(bus=3, p_mw=0.550, q_mvar=0.26638),
        Load(bus=4, p_mw=0.445, q_mvar=0.11153),
        Load(bus=5, p_mw=0.340, q_mvar=0.08521),
        Load(bus=6, p_mw=0.570, q_mvar=0.27606),
        Load(bus=7, p_mw=0.675, q_mvar=0.41833),
        Load(bus=8, p_mw=0.605, q_mvar=0.15163),
        Load(bus=9, p_mw=0.090, q_mvar=0.05578),
        Load(bus=10, p_mw=0.565, q_mvar=0.14160),
        Load(bus=11, p_mw=0.750, q_mvar=0.18797),
        Load(bus=12, p_mw=0.605, q_mvar=0.29301),
        Load(bus=13, p_mw=0.040, q_mvar=0.02479),
        Load(bus=14, p_mw=0.800, q_mvar=0.20500),
    ]
    return buses, lines, loads


class _CigreDer(NamedTuple):
    """One ``Y_DER`` row from ``CIGRE/script_CIGRE_Islanded.m``: ``[where
    DG-type load bustype V_init delta Pgen Qgen Pcons Qcons TR_R TR_XL Xd]``
    -- ``TR_R``/``TR_XL`` are always ``_CIGRE_DER_TR_R_PU``/
    ``_CIGRE_DER_TR_XL_PU`` here (see that constant's own comment) so
    they're not repeated per row.
    """

    dispatch_bus: int  # "where?" -- the raw CIGRE feeder node this DER connects to
    unit_type: UnitType
    bus_type: BusType
    v_set_pu: float
    p_set_mw: float
    q_set_mvar: float
    p_cons_mw: float
    xd_pu: float | None = None
    controller: GfmController | None = None


def _cigre_network(name: str, ders: list[_CigreDer], terminal_names: list[str]) -> Network:
    """Assembles a CIGRE-islanded preset from its DER rows the way
    ``Functions/network_form.m`` does: each DER gets its own new terminal
    bus (15, 16, 17, ... in the order given here), connected to its
    ``dispatch_bus`` through a step-up transformer.
    """
    raw_buses, lines, loads = _cigre_raw_topology()
    der_buses = [Bus(id=15 + i, name=terminal_names[i], vn_kv=_CIGRE_LV_KV) for i in range(len(ders))]
    der_units = [
        DerUnit(
            id=i + 1,
            bus=15 + i,
            unit_type=d.unit_type,
            bus_type=d.bus_type,
            v_set_pu=d.v_set_pu,
            p_set_mw=d.p_set_mw,
            q_set_mvar=d.q_set_mvar,
            p_cons_mw=d.p_cons_mw,
            controller=d.controller,
            xd_pu=d.xd_pu,
        )
        for i, d in enumerate(ders)
    ]
    transformers = [
        Transformer(
            hv_bus=d.dispatch_bus, lv_bus=15 + i, r_pu=_CIGRE_DER_TR_R_PU, x_pu=_CIGRE_DER_TR_XL_PU,
            sn_mva=_CIGRE_SN_MVA, name=f"{terminal_names[i].removesuffix('_terminal')}_xfmr",
        )
        for i, d in enumerate(ders)
    ]
    return Network(
        name=name, f_hz=50.0, sn_mva=_CIGRE_SN_MVA,
        buses=raw_buses + der_buses, lines=lines, transformers=transformers,
        loads=loads, der_units=der_units,
    )


def cigre_islanded_1sm_1gfm_1gfl() -> Network:
    """CIGRE MV feeder, islanded, 1 SM + 1 GFM (Droop) + 1 GFL.

    ``CIGRE/script_CIGRE_Islanded.m``, case ``'CIGRE_Islanded_1SM_1GFM_1GFL'``.
    """
    return _cigre_network(
        "CIGRE_Islanded_1SM_1GFM_1GFL",
        [
            _CigreDer(3, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
            _CigreDer(12, UnitType.GFM, BusType.PV, 1.0, 1.0 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(9, UnitType.GFL, BusType.PQ, 1.0, 1.0 * _CIGRE_SN_MVA, 0.0, 0.0),
        ],
        ["sm_terminal", "gfm_terminal", "gfl_terminal"],
    )


def cigre_islanded_1sm_2gfm_1gfl() -> Network:
    """CIGRE LV/MV benchmark feeder, islanded, 1 SM + 2 GFM + 1 GFL.

    Ported from ``Functions/CIGRE_raw.m`` combined with the
    ``'CIGRE_Islanded_1SM_2GFM_1GFL'`` case of
    ``CIGRE/script_CIGRE_Islanded.m`` (the DER dispatch) and assembled the
    way ``Functions/network_form.m`` does. This is the case matching the
    Simulink model actually present in this repo,
    ``CIGRE/CIGRE_Islanded_1SM_2GFM_1GFL.slx``.

    Bus numbering matches the MATLAB tool: 1-14 are the CIGRE feeder's raw
    MV nodes, 15-18 are the LV DER terminal buses (behind their
    transformers), in DER declaration order (SM, GFM, GFM, GFL).
    """
    return _cigre_network(
        "CIGRE_Islanded_1SM_2GFM_1GFL",
        [
            _CigreDer(3, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.0, 0.0 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
            _CigreDer(12, UnitType.GFM, BusType.PV, 1.0, 0.7 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(6, UnitType.GFM, BusType.PV, 1.0, 0.7 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(9, UnitType.GFL, BusType.PQ, 1.0, 0.7 * _CIGRE_SN_MVA, 0.0, 0.0),
        ],
        ["sm_terminal", "gfm1_terminal", "gfm2_terminal", "gfl_terminal"],
    )


def cigre_islanded_1sm_3gfm_1gfl() -> Network:
    """CIGRE MV feeder, islanded, 1 SM + 3 GFM (Droop) + 1 GFL.

    ``CIGRE/script_CIGRE_Islanded.m``, case ``'CIGRE_Islanded_1SM_3GFM_1GFL'``.
    """
    return _cigre_network(
        "CIGRE_Islanded_1SM_3GFM_1GFL",
        [
            _CigreDer(3, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
            _CigreDer(12, UnitType.GFM, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(6, UnitType.GFM, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(11, UnitType.GFM, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(9, UnitType.GFL, BusType.PQ, 1.0, 0.7 * _CIGRE_SN_MVA, 0.0, 0.0),
        ],
        ["sm_terminal", "gfm1_terminal", "gfm2_terminal", "gfm3_terminal", "gfl_terminal"],
    )


def cigre_islanded_2sm_2gfm_2gfl() -> Network:
    """CIGRE MV feeder, islanded, 2 SM + 2 GFM (Droop) + 2 GFL.

    ``CIGRE/script_CIGRE_Islanded.m``, case ``'CIGRE_Islanded_2SM_2GFM_2GFL'``
    -- note the second SM row there also carries a real ``Xd`` (0.2), unlike
    the single-SM variants' non-slack SM rows elsewhere; transcribed as
    given.
    """
    return _cigre_network(
        "CIGRE_Islanded_2SM_2GFM_2GFL",
        [
            _CigreDer(7, UnitType.SYNCHRONOUS_MACHINE, BusType.SLACK, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
            _CigreDer(5, UnitType.GFM, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(11, UnitType.GFM, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
            _CigreDer(9, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
            _CigreDer(6, UnitType.GFL, BusType.PQ, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0),
            _CigreDer(10, UnitType.GFL, BusType.PQ, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0),
        ],
        ["sm1_terminal", "gfm1_terminal", "gfm2_terminal", "sm2_terminal", "gfl1_terminal", "gfl2_terminal"],
    )


# SMIB presets share CIGRE's own base values (Zb_lines = 160 ohm,
# Lb_lines = 0.509296 H at Sb = 2.5 MVA, Un_MV = 20 kV, wb = 2*pi*50) --
# confirmed numerically, not assumed: SMIB_raw.m and CIGRE_raw.m use the
# identical per-km line constants (R=0.501, XL=0.716 ohm/km, same B
# constant), and plugging CIGRE's own base values into SMIB_raw.m's
# formulas reproduces CIGRE's own line r_pu/x_pu/b_pu to full precision for
# a shared length -- i.e. SMIB_raw.m was written to plug into CIGRE's own
# base-value convention, not an independent one.
_SMIB_LV_KV = 10.0
# Functions/SMIB_raw.m: Line.R_pu/L_pu/C_pu at length=61 km, computed with
# the shared base values above (Zb=160, Lb=0.509296): R_pu=(0.501*61)/160,
# L_pu=(0.716/(2*pi*50)*61)/0.509296, C_pu=(47.493e-6*61)*160.
_SMIB_LINE_R_PU = 0.191006
_SMIB_LINE_X_PU = 0.272975
_SMIB_LINE_B_PU = 0.463532
# Functions/SMIB_raw.m: Load.P/Q (W), node 2 only, converted to MW/MVAr.
_SMIB_LOAD_P_MW = 1.0
_SMIB_LOAD_Q_MVAR = 0.25
# No CIGRE-style real TR_R/TR_XL exists for SMIB (see _smib()'s docstring) --
# reusing CIGRE's own DER-transformer value, a real number from this same
# codebase, rather than an invented one.
_SMIB_TR_R_PU = 0.0002
_SMIB_TR_X_PU = 0.0240


def _smib(
    *, unit_type: UnitType, controller: GfmController | None, bus_type: BusType,
    grid_unit: UnitType = UnitType.INFINITE_BUS,
) -> Network:
    """Shared builder for the SMIB presets and their SMSM counterparts below
    (``grid_unit`` is what stands in for "the grid" on bus 4: an infinite bus,
    or a synchronous machine acting as the slack)."""
    grid_tag = "ib" if grid_unit is UnitType.INFINITE_BUS else "grid_sm"
    buses = [
        Bus(id=1, name="gen_side", vn_kv=_CIGRE_MV_KV),
        Bus(id=2, name=f"{grid_tag}_side", vn_kv=_CIGRE_MV_KV),
        Bus(id=3, name="der_terminal", vn_kv=_SMIB_LV_KV),
        Bus(id=4, name=f"{grid_tag}_terminal", vn_kv=_SMIB_LV_KV),
    ]
    lines = [
        Line(from_bus=1, to_bus=2, r_pu=_SMIB_LINE_R_PU, x_pu=_SMIB_LINE_X_PU, b_pu=_SMIB_LINE_B_PU, length_km=61.0)
    ]
    loads = [Load(bus=2, p_mw=_SMIB_LOAD_P_MW, q_mvar=_SMIB_LOAD_Q_MVAR)]

    der_units = [
        DerUnit(
            id=1, bus=4, unit_type=grid_unit, bus_type=BusType.SLACK,
            v_set_pu=1.0, p_set_mw=0.0,
            xd_pu=0.2 if grid_unit is UnitType.SYNCHRONOUS_MACHINE else None,
        ),
        DerUnit(
            id=2, bus=3, unit_type=unit_type, bus_type=bus_type,
            v_set_pu=1.0, p_set_mw=1.0, q_set_mvar=0.0, controller=controller,
        ),
    ]
    transformers = [
        Transformer(hv_bus=2, lv_bus=4, r_pu=_SMIB_TR_R_PU, x_pu=_SMIB_TR_X_PU, sn_mva=_CIGRE_SN_MVA, name=f"{grid_tag}_xfmr"),
        Transformer(hv_bus=1, lv_bus=3, r_pu=_SMIB_TR_R_PU, x_pu=_SMIB_TR_X_PU, sn_mva=_CIGRE_SN_MVA, name="der_xfmr"),
    ]

    return Network(
        name=f"{'SMIB' if grid_unit is UnitType.INFINITE_BUS else 'SMSM'}_{unit_type.value}",
        f_hz=50.0,
        sn_mva=_CIGRE_SN_MVA,
        buses=buses,
        lines=lines,
        transformers=transformers,
        loads=loads,
        der_units=der_units,
    )


def sm_smib() -> Network:
    """Single-machine-infinite-bus: one synchronous machine against an
    infinite bus.

    **Not a straight MATLAB port, unlike every other preset in this
    module.** ``Functions/preset_networks.m`` names this case (``'SM_SMIB'``,
    alongside ``'GFM_SMIB'``/``'GFL_SMIB'``) but its ``switch`` body is
    empty — the driving script was never finished in the original tool.
    The topology, line, and load values below *are* real, transcribed from
    ``Functions/SMIB_raw.m`` (confirmed to share CIGRE's own base-value
    convention numerically, not just by inspection — see the module-level
    comment above). What's genuinely reconstructed rather than ported: the
    DER-unit attachment (which bus, which transformer impedance, dispatch)
    following this project's own ``network_form.m`` convention (a DER
    behind its own transformer, on a new bus) and reusing CIGRE's real
    DER-transformer impedance (0.0002/0.0240 pu) since no SMIB-specific
    value exists to transcribe.

    This is also the first preset to actually exercise
    :mod:`g2elin_core.components.ib` (the infinite-bus component was built
    early in this project but had no real test case to wire into the
    interconnection until this).
    """
    return _smib(unit_type=UnitType.SYNCHRONOUS_MACHINE, controller=None, bus_type=BusType.PV)


def gfm_smib() -> Network:
    """Single grid-forming converter (Droop control) against an infinite
    bus. See :func:`sm_smib` for what's ported vs. reconstructed here."""
    return _smib(unit_type=UnitType.GFM, controller=GfmController.DROOP, bus_type=BusType.PV)


def gfl_smib() -> Network:
    """Single grid-following converter against an infinite bus. See
    :func:`sm_smib` for what's ported vs. reconstructed here."""
    return _smib(unit_type=UnitType.GFL, controller=None, bus_type=BusType.PQ)


# --- Single machine vs. synchronous machine (SMSM) ---------------------------
# Same line, load and transformers as the SMIB presets above, but "the grid"
# on bus 4 is a synchronous machine (the slack, with the default SM
# parameters) instead of an infinite bus -- a unit against a finite-inertia,
# voltage-regulated source rather than an ideal one. Like SMIB, these are
# reconstructions on SMIB_raw.m's real topology, not ports of a MATLAB case.


def sm_smsm() -> Network:
    """One synchronous machine (PV, 1 MW) against a synchronous machine
    acting as the grid (slack)."""
    return _smib(unit_type=UnitType.SYNCHRONOUS_MACHINE, controller=None, bus_type=BusType.PV,
                 grid_unit=UnitType.SYNCHRONOUS_MACHINE)


def gfm_smsm() -> Network:
    """One grid-forming converter (Droop, PV, 1 MW) against a synchronous
    machine acting as the grid (slack)."""
    return _smib(unit_type=UnitType.GFM, controller=GfmController.DROOP, bus_type=BusType.PV,
                 grid_unit=UnitType.SYNCHRONOUS_MACHINE)


def gfl_smsm() -> Network:
    """One grid-following converter (PQ, 1 MW, absorbing 0.9 MVAr) against a
    synchronous machine acting as the grid (slack).

    Unlike :func:`gfl_smib`, the GFL absorbs reactive power: the 61 km SMIB
    line's charging (~1 MVAr) has to go somewhere, and a GFL at Q = 0 leaves
    all of it to the grid machine. An infinite bus absorbs that without
    complaint; a synchronous machine run that far under-excited is
    small-signal unstable (a real mode in its q-axis damper flux, +0.80 1/s,
    with the GFL terminal at 1.147 pu). At -0.9 MVAr the grid machine runs
    near unity power factor, voltages stay within 1.034 pu and the case is
    stable -- sweep the GFL's Q setpoint towards 0 on the root-locus page to
    watch that mode cross.
    """
    net = _smib(unit_type=UnitType.GFL, controller=None, bus_type=BusType.PQ,
                grid_unit=UnitType.SYNCHRONOUS_MACHINE)
    gfl = next(d for d in net.der_units if d.unit_type is UnitType.GFL)
    gfl.q_set_mvar = -0.9
    return net


# --- CIGRE MV feeder, interconnected ------------------------------------------
# The grid infeed's transformer uses the same impedance as every unit
# transformer here (_CIGRE_DER_TR_R_PU/_CIGRE_DER_TR_XL_PU), so all four
# transformers in the case are identical. (Functions/preset_networks.m's
# CIGRE case lists an HV/MV substation value, 0.016 + j1.92 ohm = 0.0001 +
# j0.012 pu on the MV base; it is deliberately not used, to keep the case
# uniform.)


def _cigre_interconnected_topology() -> tuple[list[Bus], list[Line], list[Load]]:
    """``Functions/CIGRE_raw.m`` with ``S0=1`` (interconnected), ``S1=1``,
    ``S2=S3=0`` (as the islanded presets): feeder 2 (line 11) starts at
    node 1 instead of node 14, node 14 is dropped, and its load is merged
    into node 1's (``Load.P_pu(1) + Load.P_pu(14)``)."""
    buses, lines, loads = _cigre_raw_topology()
    buses = [b for b in buses if b.id != 14]
    lines = [
        ln.model_copy(update={"from_bus": 13, "to_bus": 1}) if (ln.from_bus, ln.to_bus) == (13, 14) else ln
        for ln in lines
    ]
    load14 = next(ld for ld in loads if ld.bus == 14)
    loads = [
        ld.model_copy(update={"p_mw": ld.p_mw + load14.p_mw, "q_mvar": ld.q_mvar + load14.q_mvar}) if ld.bus == 1 else ld
        for ld in loads if ld.bus != 14
    ]
    return buses, lines, loads


def cigre_interconnected_1sm_1gfm_1gfl() -> Network:
    """CIGRE MV feeder, interconnected with the upstream grid, 1 SM + 1 GFM
    (Droop) + 1 GFL.

    ``Functions/preset_networks.m`` names this case
    (``'CIGRE_Interconnected_1SM_1GFM_1GFL'``) but, like the SMIB cases, its
    driving script was never written. **Ported:** the interconnected
    topology and loads (``CIGRE_raw.m``, ``S0=1``) and the base values.
    **Reconstructed:** the upstream grid, modelled as an infinite bus (the
    slack) behind a transformer at feeder node 1 -- given the same impedance
    as the unit transformers, so all four transformers are identical -- and
    the units placed and
    dispatched as in :func:`cigre_islanded_1sm_1gfm_1gfl` -- with the SM now
    a PV unit, since the grid takes the slack role.

    The grid terminal bus is modelled at MV (20 kV, the grid referred to the
    feeder side), so the substation transformer keeps the "unit on its own
    LV-side terminal bus" convention every unit here follows.
    """
    raw_buses, lines, loads = _cigre_interconnected_topology()
    ders = [
        _CigreDer(3, UnitType.SYNCHRONOUS_MACHINE, BusType.PV, 1.0, 0.5 * _CIGRE_SN_MVA, 0.0, 0.0, xd_pu=0.2),
        _CigreDer(12, UnitType.GFM, BusType.PV, 1.0, 1.0 * _CIGRE_SN_MVA, 0.0, 0.0, controller=GfmController.DROOP),
        _CigreDer(9, UnitType.GFL, BusType.PQ, 1.0, 1.0 * _CIGRE_SN_MVA, 0.0, 0.0),
    ]
    names = ["sm_terminal", "gfm_terminal", "gfl_terminal"]
    der_buses = [Bus(id=15 + i, name=names[i], vn_kv=_CIGRE_LV_KV) for i in range(len(ders))]
    grid_bus = Bus(id=18, name="grid_terminal", vn_kv=_CIGRE_MV_KV)
    der_units = [
        DerUnit(
            id=i + 1, bus=15 + i, unit_type=d.unit_type, bus_type=d.bus_type, v_set_pu=d.v_set_pu,
            p_set_mw=d.p_set_mw, q_set_mvar=d.q_set_mvar, p_cons_mw=d.p_cons_mw, controller=d.controller,
            xd_pu=d.xd_pu,
        )
        for i, d in enumerate(ders)
    ] + [DerUnit(id=4, bus=18, unit_type=UnitType.INFINITE_BUS, bus_type=BusType.SLACK, v_set_pu=1.0, p_set_mw=0.0)]
    transformers = [
        Transformer(
            hv_bus=d.dispatch_bus, lv_bus=15 + i, r_pu=_CIGRE_DER_TR_R_PU, x_pu=_CIGRE_DER_TR_XL_PU,
            sn_mva=_CIGRE_SN_MVA, name=f"{names[i].removesuffix('_terminal')}_xfmr",
        )
        for i, d in enumerate(ders)
    ] + [Transformer(hv_bus=1, lv_bus=18, r_pu=_CIGRE_DER_TR_R_PU, x_pu=_CIGRE_DER_TR_XL_PU,
                     sn_mva=_CIGRE_SN_MVA, name="hv_mv_xfmr")]
    return Network(
        name="CIGRE_Interconnected_1SM_1GFM_1GFL", f_hz=50.0, sn_mva=_CIGRE_SN_MVA,
        buses=raw_buses + der_buses + [grid_bus], lines=lines, transformers=transformers,
        loads=loads, der_units=der_units,
    )
