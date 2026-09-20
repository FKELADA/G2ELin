"""Builds the :class:`~g2elin_core.interconnect.assemble.Block` list and
wiring rules for a :class:`~g2elin_core.network.schema.Network`. An
infinite bus (IB) is always the slack (see ``components/ib.py``), so it
gets its own ``"ib_slack"`` kind alongside ``"sm_slack"`` — see
``_slack_kind()`` below.

``build_blocks_and_wiring`` is generic over anything satisfying
:class:`~g2elin_core.interconnect.assemble.PortSpec` — both
``LinearComponent`` (used by ``assemble_network``, for modal analysis) and
the nonlinear wrapper in ``timedomain/emt.py`` (for EMT simulation) build
the *same* topology from it, since the interconnection itself (which port
connects to which) doesn't care whether a block is linear or not.
"""

from __future__ import annotations

from g2elin_core.network.breakers import TYPE_LABEL, block_labels, service_state
from g2elin_core.network.schema import Network

from .assemble import AssembledSystem, Block, PortSpec, Wiring, assemble

_KIND_BY_UNIT_TYPE = {"sm": "sm", "gfm": "gfm", "gfl": "gfl", "infinite_bus": "ib"}
_SLACK_KIND_BY_UNIT_TYPE = {"sm": "sm_slack", "infinite_bus": "ib_slack"}
_TYPE_LABEL = TYPE_LABEL
# The output a unit's own frame follows (an infinite bus has none: its frame
# simply turns at the fixed speed of its source).
_SPEED_OUTPUT = {"sm": "w_r", "gfm": "w", "gfl": "w_pll"}


def build_blocks_and_wiring(
    network: Network,
    *,
    der_components: dict[int, PortSpec],
    node_components: dict[int, PortSpec],
    line_components: list[PortSpec],
    load_components: list[PortSpec],
    frame_components: dict[int, PortSpec] | None = None,
) -> tuple[list[Block], list[Wiring]]:
    """Wires one network's components together.

    ``frame_components`` (see ``components/frame.py``) are the dq reference
    frames everything is written in, one per island, keyed by the unit that
    island is referenced to. Each frame is a block of its own, following that
    unit's speed, so no unit is load-bearing for the model: any of them can
    be disconnected, and islands each keep their own frequency. Without them
    the frame is the slack unit itself (``Network.frame_follows_slack``, the
    MATLAB toolbox's convention), which then has to be a synchronous machine
    or an infinite bus and can never leave the model.
    """
    # Only the MATLAB-compatible frame needs a slack unit: it *is* the frame.
    # With frames of their own a network needn't name one at all -- which is
    # what a mid-run trip of the slack leaves behind (timedomain/events.py).
    slack_der = next((d for d in network.der_units if d.bus_type.value == "slack"), None)
    if not frame_components and slack_der is None:
        raise ValueError("a network whose frame follows its slack unit must have one")
    if not frame_components and slack_der.unit_type.value not in _SLACK_KIND_BY_UNIT_TYPE:
        raise NotImplementedError(
            "with the frame tied to the slack unit only a synchronous-machine or infinite-bus slack is "
            "wired up; switch Network.frame_follows_slack off to use a frame of its own"
        )
    slack_kind = _SLACK_KIND_BY_UNIT_TYPE.get(slack_der.unit_type.value) if slack_der else None

    missing = {d.id for d in network.der_units} - der_components.keys()
    if missing:
        raise ValueError(
            f"der_components is missing entries for DER id(s) {sorted(missing)} -- "
            f"every id in network.der_units needs one (a common cause: der_components was "
            f"built against a different Network than the one passed here)"
        )

    # Block names (and so every state/input/output name derived from them,
    # e.g. "dw_r_{SM_2}") are per-*type* counters -- "SM_1", "GFM_1",
    # "GFM_2", "GFL_1" -- not the DER's raw (type-agnostic) id, so a name
    # says what kind of unit it is without having to cross-reference the
    # network definition. Counted in `network.der_units` declaration order,
    # which every preset already lists in ascending id order. Line/load
    # numbers likewise come from network.breakers.block_labels(), so a
    # network with elements switched out keeps the full network's names.
    labels = block_labels(network)
    der_blocks: dict[int, Block] = {}
    for der in network.der_units:
        owns_frame = not frame_components and der.id == slack_der.id
        kind = slack_kind if owns_frame else _KIND_BY_UNIT_TYPE[der.unit_type.value]
        der_blocks[der.id] = Block(name=labels.der[der.id], kind=kind, comp=der_components[der.id])

    # One frame per island, named for its reference unit when there are
    # several ("Frame" stays "Frame" while the network is in one piece).
    frame_blocks: dict[int, Block] = {}
    if frame_components:
        several = len(frame_components) > 1
        for ref_id, comp in frame_components.items():
            name = f"Frame_{labels.der[ref_id]}" if several else "Frame"
            frame_blocks[ref_id] = Block(name=name, kind="frame", comp=comp)

    node_blocks: dict[int, Block] = {
        bus_id: Block(name=f"Nd_{bus_id}", kind="node", comp=comp)
        for bus_id, comp in node_components.items()
    }
    line_blocks = [
        Block(name=f"Ln_{labels.line[i] + 1}", kind="line", comp=comp) for i, comp in enumerate(line_components)
    ]
    load_blocks = [
        Block(name=f"Ld_{labels.load[i] + 1}", kind="load", comp=comp) for i, comp in enumerate(load_components)
    ]

    # Order matches script_generic.m's concatenation order: slack, other DGs,
    # nodes, lines, loads. Only cosmetic (state/output ordering), not required
    # for correctness, but kept for easy visual comparison with a MATLAB run.
    ordered_ders = ([slack_der] if slack_der else []) + [
        d for d in network.der_units if slack_der is None or d.id != slack_der.id
    ]
    blocks = (
        list(frame_blocks.values())
        + [der_blocks[d.id] for d in ordered_ders]
        + [node_blocks[b.id] for b in network.buses if b.id in node_blocks]
        + line_blocks
        + load_blocks
    )

    transformer_by_lv_bus = {tr.lv_bus: tr for tr in network.transformers}
    wiring: list[Wiring] = []

    # Which frame each bus belongs to: its island's (the slack unit's own
    # block in the MATLAB-compatible mode, where there is only ever one).
    if frame_blocks:
        island_of_bus: dict[int, Block] = {}
        for island in service_state(network).islands:
            frame = frame_blocks.get(island.reference)
            if frame is None:  # an island of a network built without one (shouldn't happen)
                frame = next(iter(frame_blocks.values()))
            for bus_id in island.buses:
                island_of_bus[bus_id] = frame
        frame_of = lambda bus_id: island_of_bus.get(bus_id, next(iter(frame_blocks.values())))  # noqa: E731
        # Each frame follows the speed of the unit its island is referenced to
        # (an infinite-bus reference turns at its own fixed speed: no input).
        for ref_id, frame in frame_blocks.items():
            speed = _SPEED_OUTPUT.get(next(d.unit_type.value for d in network.der_units if d.id == ref_id))
            if speed is not None:
                wiring.append(Wiring(frame, "w_in", [(der_blocks[ref_id], speed, 1.0)]))
    else:
        slack_block = der_blocks[slack_der.id]
        frame_of = lambda bus_id: slack_block  # noqa: E731

    # theta_g / wg: every unit that doesn't own the frame reads its angle, and
    # every node/line/load (and an infinite bus, for its own rotational terms)
    # its speed.
    for der in network.der_units:
        if not frame_blocks and slack_der is not None and der.id == slack_der.id:
            continue
        db = der_blocks[der.id]
        frame = frame_of(der.bus)
        wiring.append(Wiring(db, "theta_g", [(frame, "theta", 1.0)]))
        if db.kind == "ib":
            wiring.append(Wiring(db, "wg", [(frame, "wr", 1.0)]))
    for bus_id, b in node_blocks.items():
        wiring.append(Wiring(b, "wg", [(frame_of(bus_id), "wr", 1.0)]))
    for ln, b in zip(network.lines, line_blocks):
        wiring.append(Wiring(b, "wg", [(frame_of(ln.from_bus), "wr", 1.0)]))
    for ld, b in zip(network.loads, load_blocks):
        wiring.append(Wiring(b, "wg", [(frame_of(ld.bus), "wr", 1.0)]))

    # DG voltage input = its raw node's voltage output.
    for der in network.der_units:
        tr = transformer_by_lv_bus[der.bus]
        node_block = node_blocks[tr.hv_bus]
        db = der_blocks[der.id]
        wiring.append(Wiring(db, "vgd_g", [(node_block, "vgd_g", 1.0)]))
        wiring.append(Wiring(db, "vgq_g", [(node_block, "vgq_g", 1.0)]))

    # Line voltages = its from/to nodes' voltage outputs.
    for ln, block in zip(network.lines, line_blocks):
        from_node, to_node = node_blocks[ln.from_bus], node_blocks[ln.to_bus]
        wiring.append(Wiring(block, "vgdj_g", [(from_node, "vgd_g", 1.0)]))
        wiring.append(Wiring(block, "vgqj_g", [(from_node, "vgq_g", 1.0)]))
        wiring.append(Wiring(block, "vgdk_g", [(to_node, "vgd_g", 1.0)]))
        wiring.append(Wiring(block, "vgqk_g", [(to_node, "vgq_g", 1.0)]))

    # Load voltage = its node's voltage output.
    for ld, block in zip(network.loads, load_blocks):
        node_block = node_blocks[ld.bus]
        wiring.append(Wiring(block, "vgd_g", [(node_block, "vgd_g", 1.0)]))
        wiring.append(Wiring(block, "vgq_g", [(node_block, "vgq_g", 1.0)]))

    # Node current balance: +DG current in, +"to"-line current in,
    # -"from"-line current out, -load current out.
    der_by_grid_bus: dict[int, Block] = {
        transformer_by_lv_bus[d.bus].hv_bus: der_blocks[d.id] for d in network.der_units
    }
    for bus_id, node_block in node_blocks.items():
        terms_d: list[tuple[Block, str, float]] = []
        terms_q: list[tuple[Block, str, float]] = []
        if bus_id in der_by_grid_bus:
            terms_d.append((der_by_grid_bus[bus_id], "igd_g", 1.0))
            terms_q.append((der_by_grid_bus[bus_id], "igq_g", 1.0))
        for ln, block in zip(network.lines, line_blocks):
            if ln.to_bus == bus_id:
                terms_d.append((block, "ild_g", 1.0))
                terms_q.append((block, "ilq_g", 1.0))
            if ln.from_bus == bus_id:
                terms_d.append((block, "ild_g", -1.0))
                terms_q.append((block, "ilq_g", -1.0))
        for ld, block in zip(network.loads, load_blocks):
            if ld.bus == bus_id:
                terms_d.append((block, "icd_g", -1.0))
                terms_q.append((block, "icq_g", -1.0))
        wiring.append(Wiring(node_block, "ishd_g", terms_d))
        wiring.append(Wiring(node_block, "ishq_g", terms_q))

    return blocks, wiring


def assemble_network(
    network: Network,
    *,
    der_components: dict[int, PortSpec],
    node_components: dict[int, PortSpec],
    line_components: list[PortSpec],
    load_components: list[PortSpec],
    frame_components: dict[int, PortSpec] | None = None,
) -> AssembledSystem:
    blocks, wiring = build_blocks_and_wiring(
        network,
        der_components=der_components,
        node_components=node_components,
        line_components=line_components,
        load_components=load_components,
        frame_components=frame_components,
    )
    return assemble(blocks, wiring)
