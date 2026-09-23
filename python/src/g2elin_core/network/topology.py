"""Computes a 2D layout for a :class:`~g2elin_core.network.schema.Network`'s
buses/lines/transformers, for the web UI's Network and Power Flow tabs.

A DER unit or load in this schema's own convention isn't a separate
geometric node — a DER sits on its own real :class:`Bus` (``network_form.m``
always gives it one, behind its own transformer), and a load just
references an existing bus id. So the graph to lay out is exactly
``network.buses`` (nodes) + ``network.lines``/``network.transformers``
(edges); DER/load presence is metadata decorating an existing bus node,
not a node of its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from g2elin_core.operating_point import (
    gfl_params, gfm_params, overridable_param_keys, rebase_params, sm_params, unit_transformer_rx,
)

from .schema import DerUnit, Network


@dataclass(frozen=True)
class BusNode:
    id: int
    name: str
    vn_kv: float
    x: float
    y: float
    unit_type: str | None  # DerUnit.unit_type.value if a DER sits here, else None
    load_p_mw: float | None  # sum of Load.p_mw at this bus, else None
    der_info: dict | None = None  # DerUnit's own fields + derived control params, if a DER sits here


@dataclass(frozen=True)
class TopologyEdge:
    kind: str  # "line" | "transformer"
    from_bus: int
    to_bus: int
    name: str
    params: dict = field(default_factory=dict)  # electrical parameters, keyed by name -- see below


@dataclass(frozen=True)
class NetworkTopology:
    nodes: list[BusNode]
    edges: list[TopologyEdge]


def _der_info(der: DerUnit, network: Network) -> dict:
    """The unit's own dispatch/rating fields plus, for SM/GFM/GFL, the
    derived electrical/control-parameter dict (``sm_params()`` and
    friends) actually used to build that unit's dynamic model -- exposed
    for the web UI's "click a unit to see its parameters" panel. IB has
    no separate control-parameter set (see ``components/ib.py``), so it
    only gets the dispatch fields.
    """
    info: dict = {
        "unit_type": der.unit_type.value,
        "bus_type": der.bus_type.value,
        "v_set_pu": der.v_set_pu,
        "p_set_mw": der.p_set_mw,
        "q_set_mvar": der.q_set_mvar,
        "p_cons_mw": der.p_cons_mw,
        "q_cons_mvar": der.q_cons_mvar,
        "controller": der.controller.value if der.controller else None,
        "xd_pu": der.xd_pu,
    }
    if not network.transformers:
        return info
    rt, lt = unit_transformer_rx(network, der)
    if der.unit_type.value == "sm":
        defaults = sm_params(sn_mva=network.sn_mva, f_hz=network.f_hz, rt_pu=rt, lt_pu=lt)
    elif der.unit_type.value in ("gfm", "gfl"):
        un_kv = network.bus(der.bus).vn_kv
        fn = gfm_params if der.unit_type.value == "gfm" else gfl_params
        defaults = fn(sn_mva=network.sn_mva, f_hz=network.f_hz, un_kv=un_kv, rt_pu=rt, lt_pu=lt)
    else:
        return info
    # control_params: what the model actually uses (defaults + this unit's
    # valid overrides, rebased when the unit states its own rating, exactly as
    # operating_point.apply_param_overrides does for the model itself);
    # control_params_default: the defaults alone, so the UI can show which
    # values were changed and restore them. Unknown names are dropped rather
    # than raised on here -- validate_network reports those.
    valid = overridable_param_keys(der.unit_type.value)
    overrides = {k: v for k, v in der.params.items() if k in valid}
    if der.sn_mva is not None:
        overrides = rebase_params(overrides, from_mva=der.sn_mva, to_mva=network.sn_mva)
    info["control_params"] = {**defaults, **overrides}
    info["control_params_default"] = defaults
    return info


# Above this many buses of its own, a network's unit terminals are hung off
# the bus they feed instead of being laid out as nodes. Below it the spring
# layout has room to place them properly, and does it better: a fixed offset
# pushes a machine onto the very branches its bus connects to.
_HANG_TERMINALS_ABOVE = 40
# How far a hung terminal sits from its bus, in the layout's own roughly
# [-1, 1] coordinates.
_TERMINAL_OFFSET = 0.16


def _unit_terminal_buses(network: Network) -> dict[int, int]:
    """``{terminal bus -> the grid bus it feeds}`` for every unit that sits on
    a bus of its own behind a single step-up, with no line or load on it."""
    line_buses = {b for ln in network.lines for b in (ln.from_bus, ln.to_bus)}
    load_buses = {ld.bus for ld in network.loads}
    out: dict[int, int] = {}
    for der in network.der_units:
        touching = [t for t in network.transformers if der.bus in (t.lv_bus, t.hv_bus)]
        if len(touching) == 1 and der.bus not in line_buses and der.bus not in load_buses:
            tr = touching[0]
            out[der.bus] = tr.hv_bus if tr.lv_bus == der.bus else tr.lv_bus
    return out


def _hung_terminal_positions(network, terminal_of, pos):
    """Place each unit terminal just off the bus it feeds, pointing away from
    that bus's own branches so the machine does not land on top of one, and
    fanned when a bus carries more than one unit."""
    neighbours: dict[int, list[int]] = {}
    for a, b in (
        [(ln.from_bus, ln.to_bus) for ln in network.lines]
        + [(t.hv_bus, t.lv_bus) for t in network.transformers]
    ):
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)

    per_grid_bus: dict[int, list[int]] = {}
    for terminal, grid_bus in terminal_of.items():
        per_grid_bus.setdefault(grid_bus, []).append(terminal)

    placed: dict[int, "np.ndarray"] = {}
    for grid_bus, terminals in per_grid_bus.items():
        anchor = pos.get(grid_bus)
        if anchor is None:
            continue
        # The direction the bus's branches least occupy: the reverse of their
        # summed unit vectors.
        away = np.zeros(2)
        for other in neighbours.get(grid_bus, []):
            q = pos.get(other)
            if q is None:
                continue
            d = anchor - q
            n = float(np.hypot(*d))
            if n > 1e-9:
                away += d / n
        base = math.atan2(away[1], away[0]) if float(np.hypot(*away)) > 1e-9 else 0.6 * grid_bus
        spread = math.pi / 3
        for k, terminal in enumerate(sorted(terminals)):
            offset = 0.0 if len(terminals) == 1 else (k / (len(terminals) - 1) - 0.5) * 2 * spread
            angle = base + offset
            placed[terminal] = anchor + _TERMINAL_OFFSET * np.array([math.cos(angle), math.sin(angle)])
    return placed


def compute_topology_layout(network: Network) -> NetworkTopology:
    """Kamada-Kawai layout (a graph-distance-preserving spring layout, reads
    well for the mostly-tree-like/lightly-meshed topologies these presets
    have) over the bus graph, normalized to roughly [-1, 1] per axis.
    Falls back to a spring layout if Kamada-Kawai can't be computed (e.g. a
    disconnected graph, which none of the current presets produce, but
    nothing here assumes it can't happen for a future one).
    """
    # A unit's own terminal bus -- the LV side of its step-up, carrying nothing
    # else -- is a modelling device rather than a place in the network. On a
    # large case laying each one out as a node of its own pulls the real
    # topology apart and spends space on it (54 extra nodes among the 118-bus
    # case's own), so there they are hung off the bus they feed instead. On a
    # small one the spring layout has room and places them better than any
    # fixed offset can, so it keeps them.
    all_terminals = _unit_terminal_buses(network)
    hang = len(network.buses) - len(all_terminals) > _HANG_TERMINALS_ABOVE
    terminal_of = all_terminals if hang else {}

    graph = nx.Graph()
    for bus in network.buses:
        if bus.id not in terminal_of:
            graph.add_node(bus.id)
    for line in network.lines:
        graph.add_edge(line.from_bus, line.to_bus)
    for tr in network.transformers:
        if tr.lv_bus not in terminal_of:
            graph.add_edge(tr.hv_bus, tr.lv_bus)

    try:
        pos = nx.kamada_kawai_layout(graph)
    except (nx.NetworkXException, ZeroDivisionError):
        pos = nx.spring_layout(graph, seed=0)

    if terminal_of:
        pos.update(_hung_terminal_positions(network, terminal_of, pos))

    der_by_bus = {der.bus: der for der in network.der_units}
    load_p_by_bus: dict[int, float] = {}
    for load in network.loads:
        load_p_by_bus[load.bus] = load_p_by_bus.get(load.bus, 0.0) + load.p_mw

    nodes = [
        BusNode(
            id=bus.id,
            name=bus.name or f"bus{bus.id}",
            vn_kv=bus.vn_kv,
            x=float(pos[bus.id][0]),
            y=float(pos[bus.id][1]),
            unit_type=der_by_bus[bus.id].unit_type.value if bus.id in der_by_bus else None,
            load_p_mw=load_p_by_bus.get(bus.id),
            der_info=_der_info(der_by_bus[bus.id], network) if bus.id in der_by_bus else None,
        )
        for bus in network.buses
    ]
    edges = [
        TopologyEdge(
            kind="line", from_bus=line.from_bus, to_bus=line.to_bus, name=line.name,
            params={"r_pu": line.r_pu, "x_pu": line.x_pu, "b_pu": line.b_pu, "length_km": line.length_km},
        )
        for line in network.lines
    ] + [
        TopologyEdge(
            kind="transformer", from_bus=tr.hv_bus, to_bus=tr.lv_bus, name=tr.name,
            params={"r_pu": tr.r_pu, "x_pu": tr.x_pu, "sn_mva": tr.sn_mva},
        )
        for tr in network.transformers
    ]
    return NetworkTopology(nodes=nodes, edges=edges)
