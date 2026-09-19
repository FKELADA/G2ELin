"""End-to-end: solved power flow -> linearized closed-loop small-signal model.

A synchronous-machine or infinite-bus slack is wired up (see
``interconnect/network_assembly.py``); GFM and GFL are supported as
non-slack units.
"""

from __future__ import annotations

from g2elin_core.components.gfl import linearize_gfl
from g2elin_core.components.gfm import linearize_gfm
from g2elin_core.components.ib import linearize_ib
from g2elin_core.components.line import linearize_line
from g2elin_core.components.load import linearize_load
from g2elin_core.components.node import linearize_node
from g2elin_core.components.sm import linearize_sm
from g2elin_core.interconnect import AssembledSystem, assemble_network
from g2elin_core.network.breakers import node_b_pu
from g2elin_core.network.schema import Network
from g2elin_core.operating_point import compute_operating_point
from g2elin_core.powerflow import PowerFlowResult


def linearize_network(network: Network, result: PowerFlowResult) -> AssembledSystem:
    if not result.converged:
        raise ValueError("power flow did not converge; can't linearize an unsolved operating point")

    op = compute_operating_point(network, result)
    wb_val = 2 * 3.141592653589793 * network.f_hz

    der_components = {}
    for der_id, sm_op in op.sm_ops.items():
        der_components[der_id] = linearize_sm(sm_op)
    for der_id, gfm_op in op.gfm_ops.items():
        der_components[der_id] = linearize_gfm(gfm_op)
    for der_id, gfl_op in op.gfl_ops.items():
        der_components[der_id] = linearize_gfl(gfl_op)
    for der_id, ib_kwargs in op.ib_ops.items():
        der_components[der_id] = linearize_ib(**ib_kwargs)

    missing = {d.id for d in network.der_units} - der_components.keys()
    if missing:
        raise NotImplementedError(f"unsupported DER unit type(s) for ids {sorted(missing)}")

    # Node quirk: every node uses the *first* line's susceptance (see
    # g2elin_core.components.sm module docstring) -- every bus's own dynamic
    # model (components/node.py) needs this as a nonzero denominator
    # (dv/dt = (wb/Cl)*(...)); a network with no lines at all (e.g.
    # hand-built with only transformers) has no "line #1" to borrow a value
    # from, and 0.0 isn't a safe fallback -- it's a real division by zero,
    # not just an unrealistic approximation. validate_network() already
    # catches this before it reaches here; this is the defensive backstop.
    b_pu_quirk = node_b_pu(network)
    if b_pu_quirk is None:
        raise ValueError(
            "this network has no Line elements -- every bus's own dynamic model needs a line-charging "
            "susceptance (b_pu) to linearize around, which this codebase always borrows from the first "
            "Line in the network; a network built entirely from transformers has no such source and "
            "can't run modal analysis or EMT. Add at least one Line (even a short one with a small b_pu)"
        )
    node_components = {
        bus_id: linearize_node(wb_val=wb_val, b_pu=b_pu_quirk, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        for bus_id, (vgd, vgq) in op.node_vg.items()
    }

    line_components = [
        linearize_line(wb_val=wb_val, r_pu=ln.r_pu, x_pu=ln.x_pu, wg0=1.0, ild_g0=i0[0], ilq_g0=i0[1])
        for ln, i0 in zip(network.lines, op.line_i0)
    ]

    load_components = []
    for idx, load in enumerate(network.loads):
        r_pu, x_pu = op.load_rx[idx]
        vgd, vgq = op.node_vg[load.bus]
        load_components.append(
            linearize_load(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        )

    return assemble_network(
        network,
        der_components=der_components,
        node_components=node_components,
        line_components=line_components,
        load_components=load_components,
    )
