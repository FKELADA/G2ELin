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

from g2elin_core.network.breakers import TYPE_LABEL, block_labels
from g2elin_core.network.schema import Network

from .assemble import AssembledSystem, Block, PortSpec, Wiring, assemble

_KIND_BY_UNIT_TYPE = {"sm": "sm", "gfm": "gfm", "gfl": "gfl"}
_SLACK_KIND_BY_UNIT_TYPE = {"sm": "sm_slack", "infinite_bus": "ib_slack"}
_TYPE_LABEL = TYPE_LABEL


def build_blocks_and_wiring(
    network: Network,
    *,
    der_components: dict[int, PortSpec],
    node_components: dict[int, PortSpec],
    line_components: list[PortSpec],
    load_components: list[PortSpec],
) -> tuple[list[Block], list[Wiring]]:
    slack_der = next(d for d in network.der_units if d.bus_type.value == "slack")
    if slack_der.unit_type.value not in _SLACK_KIND_BY_UNIT_TYPE:
        raise NotImplementedError("only a synchronous-machine or infinite-bus slack is wired up so far")
    slack_kind = _SLACK_KIND_BY_UNIT_TYPE[slack_der.unit_type.value]

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
        kind = slack_kind if der.id == slack_der.id else _KIND_BY_UNIT_TYPE[der.unit_type.value]
        der_blocks[der.id] = Block(name=labels.der[der.id], kind=kind, comp=der_components[der.id])

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
    ordered_ders = [slack_der] + [d for d in network.der_units if d.id != slack_der.id]
    blocks = (
        [der_blocks[d.id] for d in ordered_ders]
        + [node_blocks[b.id] for b in network.buses if b.id in node_blocks]
        + line_blocks
        + load_blocks
    )

    slack_block = der_blocks[slack_der.id]
    transformer_by_lv_bus = {tr.lv_bus: tr for tr in network.transformers}
    wiring: list[Wiring] = []

    # theta_g / wg: every non-slack DG's theta_g, and every node/line/load's
    # wg, equals the slack unit's own theta / wr.
    for der in network.der_units:
        if der.id == slack_der.id:
            continue
        wiring.append(Wiring(der_blocks[der.id], "theta_g", [(slack_block, "theta", 1.0)]))
    for b in list(node_blocks.values()) + line_blocks + load_blocks:
        wiring.append(Wiring(b, "wg", [(slack_block, "wr", 1.0)]))

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
) -> AssembledSystem:
    blocks, wiring = build_blocks_and_wiring(
        network,
        der_components=der_components,
        node_components=node_components,
        line_components=line_components,
        load_components=load_components,
    )
    return assemble(blocks, wiring)
