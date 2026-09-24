"""End-to-end: solved power flow -> linearized closed-loop small-signal model.

A synchronous-machine or infinite-bus slack is wired up (see
``interconnect/network_assembly.py``); GFM and GFL are supported as
non-slack units.
"""

from __future__ import annotations

from g2elin_core.components.frame import linearize_frame
from g2elin_core.components.gfl import linearize_gfl
from g2elin_core.components.gfm import linearize_gfm
from g2elin_core.components.ib import linearize_ib
from g2elin_core.components.line import linearize_line
from g2elin_core.components.load import linearize_load
from g2elin_core.components.node import linearize_node
from g2elin_core.components.sm import linearize_sm
from g2elin_core.interconnect import AssembledSystem, assemble_network
from g2elin_core.network.breakers import frame_references, node_capacitances
from g2elin_core.network.schema import Network
from g2elin_core.operating_point import NetworkOperatingPoint, compute_operating_point
from g2elin_core.powerflow import PowerFlowResult


def linearize_network(network: Network, result: PowerFlowResult) -> AssembledSystem:
    """The closed-loop small-signal model of one network at one operating point."""
    return assemble_network(network, **linear_components(network, result))


def linear_components(
    network: Network, result: PowerFlowResult, *, op: NetworkOperatingPoint | None = None
) -> dict:
    """Every element's own linearised state-space, before they are wired
    together -- what :func:`linearize_network` assembles, and what a caller
    needs to know which block is which (see
    :mod:`g2elin_core.modal.parameters`).

    ``op`` is the network's operating point, computed here when not given. A
    caller that needs the operating point in its own right can pass the one
    it already has rather than have it built a second time.
    """
    if not result.converged:
        raise ValueError("power flow did not converge; can't linearize an unsolved operating point")

    if op is None:
        op = compute_operating_point(network, result)
    wb_val = 2 * 3.141592653589793 * network.f_hz
    # Which dynamics each element keeps (network/schema.ModelOptions). The
    # default is "all of them", i.e. exactly the model this pipeline built
    # before model-order reduction existed.
    net_modes = network.models.group_modes_for("network")
    fixed_f = network.models.fixed_network_frequency

    unit_modes = {d.id: network.unit_modes(d) for d in network.der_units if d.unit_type.value != "infinite_bus"}
    der_components = {}
    for der_id, sm_op in op.sm_ops.items():
        der_components[der_id] = linearize_sm(sm_op, modes=unit_modes[der_id])
    for der_id, gfm_op in op.gfm_ops.items():
        der_components[der_id] = linearize_gfm(gfm_op, modes=unit_modes[der_id])
    for der_id, gfl_op in op.gfl_ops.items():
        der_components[der_id] = linearize_gfl(gfl_op, modes=unit_modes[der_id])
    for der_id, ib_kwargs in op.ib_ops.items():
        der_components[der_id] = linearize_ib(**ib_kwargs)
    # One reference frame per island (components/frame.py), each following
    # the unit that island is referenced to -- an infinite bus turns at its
    # own fixed speed, a machine or grid-forming converter carries its frame
    # with it. All of them start at the operating point's reference angle.
    frame_components = None if network.frame_follows_slack else {
        ref: linearize_frame(wb_val=wb_val, theta0=op.theta_g_rad, driven=driven)
        for ref, driven in frame_references(network).items()
    }

    missing = {d.id for d in network.der_units} - der_components.keys()
    if missing:
        raise NotImplementedError(f"unsupported DER unit type(s) for ids {sorted(missing)}")

    # Each bus's own capacitance: half the charging of the lines on it plus
    # its capacitor banks, or -- in the MATLAB-compatible mode the ported
    # presets use -- the first line's charging shared by every bus. Raises
    # rather than divide by zero when a bus has none (node_capacitances).
    node_b = node_capacitances(network)
    if not node_b:
        raise ValueError(
            "this network has no Line elements -- with Network.nodes_share_first_line_b every bus's "
            "dynamic model borrows the first Line's charging susceptance (b_pu), and a network built "
            "entirely from transformers has no such source, so it can't run modal analysis or EMT. "
            "Switch that option off to give each bus its own capacitance, or add at least one Line."
        )
    node_components = {
        bus_id: linearize_node(
            wb_val=wb_val, b_pu=node_b[bus_id], wg0=1.0, vgd_g0=vgd, vgq_g0=vgq,
            mode=net_modes["nodes"], fixed_frequency=fixed_f,
        )
        for bus_id, (vgd, vgq) in op.node_vg.items()
    }

    line_components = [
        linearize_line(
            wb_val=wb_val, r_pu=ln.r_pu, x_pu=ln.x_pu, wg0=1.0, ild_g0=i0[0], ilq_g0=i0[1],
            mode=net_modes["lines"], fixed_frequency=fixed_f,
        )
        for ln, i0 in zip(network.lines, op.line_i0)
    ]

    load_components = []
    for idx, load in enumerate(network.loads):
        r_pu, x_pu = op.load_rx[idx]
        vgd, vgq = op.node_vg[load.bus]
        load_components.append(
            linearize_load(
                wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq,
                mode=net_modes["loads"], fixed_frequency=fixed_f,
            )
        )

    # Shunt reactors: the same RL branch a line uses, with its far end at
    # zero volts (the assembly wires it that way). Capacitor banks are
    # already in node_b above, and add no block at all.
    shunt_components = {
        idx: linearize_line(
            wb_val=wb_val, r_pu=rx[0], x_pu=rx[1], wg0=1.0,
            ild_g0=op.shunt_i0[idx][0], ilq_g0=op.shunt_i0[idx][1],
            mode=net_modes["shunts"], fixed_frequency=fixed_f,
        )
        for idx, rx in op.shunt_rx.items()
    }

    # Branch transformers: the same RL branch a line uses, with the ideal
    # transformer's ratio applied by the wiring (network_assembly).
    transformer_components = {
        idx: linearize_line(
            wb_val=wb_val, r_pu=rx[0], x_pu=rx[1], wg0=1.0,
            ild_g0=op.transformer_i0[idx][0], ilq_g0=op.transformer_i0[idx][1],
            mode=net_modes["transformers"], fixed_frequency=fixed_f,
        )
        for idx, rx in op.transformer_rx.items()
    }

    return dict(
        der_components=der_components,
        node_components=node_components,
        line_components=line_components,
        load_components=load_components,
        frame_components=frame_components,
        shunt_components=shunt_components,
        transformer_components=transformer_components,
    )
