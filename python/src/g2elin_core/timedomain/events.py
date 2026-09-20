"""Network events for EMT simulations: a breaker opening, a load step, or a
phase jump, applied at t = 0 to a model sitting at its operating point.

A breaker opening or a load step changes the model itself, so the run after
the event uses a **post-event model**: the same component instances as the
pre-event model (same parameters, same controller references and setpoints --
nothing is re-dispatched by a new power flow), rewired without the switched
element, or with one load's impedance changed. Its initial state is the
pre-event state, taken block by block by name (states are continuous across
the event; the algebraic variables re-solve instantly).

- **Breaker opening** on a line, a load, or a unit (a unit's transformer is
  inside the unit's model, so opening a unit transformer's breaker trips the
  unit). What the opening cuts off stays in the simulation: an islanded
  group of buses keeps evolving on its own (a load-only island decays, an
  island with a grid-forming unit keeps running) -- unlike the steady-state
  analyses, which drop de-energized parts (see ``network.breakers``). Any
  unit can be tripped, the power flow's slack included: the reference frame
  is a block of its own (``components/frame.py``) and simply moves to
  another unit. Only the MATLAB-compatible frame
  (``Network.frame_follows_slack``), which *is* the slack unit, keeps it.
- **Load step**: the load's active/reactive power change by ``dp_pct`` /
  ``dq_pct`` percent at the voltage of the operating point (the load is a
  constant impedance, so its conductance and susceptance scale by
  ``1 + dp_pct/100`` and ``1 + dq_pct/100``).
- **Phase jump** by ``angle_deg``: at a network bus, an instantaneous
  rotation of that bus's voltage phasor (the bus's shunt-capacitor state);
  at the infinite bus's own bus, a jump of the source's voltage angle --
  modelled, the frame being tied to the source, by rotating every other
  quantity by ``-angle`` (node voltages, line/load/source currents, and each
  other unit's angle). A phase jump doesn't change the model, only its
  initial state, so the linearized model can reproduce it too.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np

from g2elin_core.interconnect import build_blocks_and_wiring, compute_topology
from g2elin_core.network.breakers import BlockLabels, block_labels, frame_references
from g2elin_core.network.schema import Network

from .emt import NonlinearNetworkModel, nonlinear_frame_block, nonlinear_load_block

EVENT_KINDS = ("breaker", "load_step", "phase_jump")
BREAKER_ELEMENTS = ("line", "transformer", "load", "unit")
# Angle state of each non-slack unit kind (rotated by an infinite-bus phase jump).
_ANGLE_STATE = {"sm": "theta", "gfm": "theta", "gfl": "theta_pll"}
_COMMON_FRAME_KINDS = ("node", "line", "load", "ib_slack", "ib")


class EventError(ValueError):
    """The event can't be applied to this network (a 422 for the API)."""


@dataclass(frozen=True)
class NetworkEvent:
    kind: str                     # "breaker" | "load_step" | "phase_jump"
    element: str = ""             # breaker: "line" | "transformer" | "load" | "unit"
    index: int = 0                # line/transformer/load index in the full network, or a unit's id
    dp_pct: float = 0.0           # load_step
    dq_pct: float = 0.0
    bus: int | None = None        # phase_jump
    angle_deg: float = 0.0


@dataclass
class AppliedEvent:
    model: NonlinearNetworkModel  # the model to integrate from t = 0 (the pre-event one itself for a phase jump)
    x0: np.ndarray                # its initial state
    description: str
    # State offset w.r.t. the pre-event equilibrium when the event is a pure
    # initial-state change (phase jump) -- what the linear overlay needs;
    # None when the event changes the model (no linear equivalent).
    dx0: np.ndarray | None


def apply_event(model: NonlinearNetworkModel, event: NetworkEvent) -> AppliedEvent:
    if event.kind == "breaker":
        return _breaker(model, event)
    if event.kind == "load_step":
        return _load_step(model, event)
    if event.kind == "phase_jump":
        return _phase_jump(model, event)
    raise EventError(f"event kind must be one of {list(EVENT_KINDS)}")


# --- breaker ---------------------------------------------------------------
def _breaker(model: NonlinearNetworkModel, ev: NetworkEvent) -> AppliedEvent:
    net: Network = model.network
    lab = block_labels(net)
    lines, loads, trs, ders = list(net.lines), list(net.loads), list(net.transformers), list(net.der_units)
    line_l, load_l, tr_l = list(lab.line), list(lab.load), list(lab.transformer)
    buses = list(net.buses)

    def drop_unit(der_id: int) -> str:
        der = next((d for d in ders if d.id == der_id), None)
        if der is None:
            raise EventError(f"unit id {der_id} is already out of service")
        if der.bus_type.value == "slack" and getattr(net, "frame_follows_slack", False):
            raise EventError(
                f"unit id {der_id} is the slack, and this network's frame follows it "
                "(Network.frame_follows_slack) -- switch that off to be able to trip it"
            )
        ders.remove(der)
        for j, tr in enumerate(trs):
            if tr.lv_bus == der.bus:
                del trs[j], tr_l[j]
                break
        buses[:] = [b for b in buses if b.id != der.bus]
        return f"unit {der_id} ({lab.der[der_id]}) tripped"

    if ev.element == "line":
        if ev.index not in line_l:
            raise EventError(f"line #{ev.index} is already out of service (or doesn't exist)")
        k = line_l.index(ev.index)
        ln = lines[k]
        del lines[k], line_l[k]
        what = f"line #{ev.index} ({ln.from_bus} → {ln.to_bus}) opened"
    elif ev.element == "load":
        if ev.index not in load_l:
            raise EventError(f"load #{ev.index} is already out of service (or doesn't exist)")
        k = load_l.index(ev.index)
        what = f"load #{ev.index} (bus {loads[k].bus}) disconnected"
        del loads[k], load_l[k]
    elif ev.element == "unit":
        what = drop_unit(ev.index)
    elif ev.element == "transformer":
        if ev.index not in tr_l:
            raise EventError(f"transformer #{ev.index} is already out of service (or doesn't exist)")
        tr = trs[tr_l.index(ev.index)]
        der = next((d for d in ders if d.bus == tr.lv_bus), None)
        if der is None:
            raise EventError(
                f"transformer #{ev.index} doesn't connect a unit -- the dynamic model only has unit "
                "transformers (inside each unit's model)"
            )
        what = f"transformer #{ev.index} opened: " + drop_unit(der.id)
    else:
        raise EventError(f"breaker element must be one of {list(BREAKER_ELEMENTS)}")

    post = net.model_copy(update=dict(buses=buses, lines=lines, loads=loads, transformers=trs, der_units=ders))
    post._labels = BlockLabels(
        der={d.id: lab.der[d.id] for d in ders}, line=line_l, load=load_l, transformer=tr_l, node_b_pu=lab.node_b_pu,
    )
    new = rebuild(model, post)
    return AppliedEvent(new, new.initial_state(), what, None)


# --- load step ---------------------------------------------------------------
def _load_step(model: NonlinearNetworkModel, ev: NetworkEvent) -> AppliedEvent:
    net: Network = model.network
    lab = block_labels(net)
    if ev.index not in lab.load:
        raise EventError(f"load #{ev.index} is out of service (or doesn't exist)")
    k = lab.load.index(ev.index)
    fp, fq = 1 + ev.dp_pct / 100.0, 1 + ev.dq_pct / 100.0
    if fp < 0:
        raise EventError("the active-power step can't go below -100 %")
    if fq <= 0:
        raise EventError(
            "the reactive-power step must stay above -100 % (the load model needs a nonzero reactance; "
            "to disconnect the load, open its breaker instead)"
        )
    r, x = model.op.load_rx[k]
    y = 1 / complex(r, x)                 # G - jB
    y_new = complex(y.real * fp, y.imag * fq)
    z_new = 1 / y_new
    load = net.loads[k]
    vgd, vgq = model.op.node_vg[load.bus]
    comp = nonlinear_load_block(
        wb_val=2 * math.pi * net.f_hz, r_pu=z_new.real, x_pu=z_new.imag, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq,
    )
    new = rebuild(model, net, replace={f"Ld_{ev.index + 1}": comp})
    sign = lambda v: f"{v:+g} %"  # noqa: E731
    return AppliedEvent(
        new, new.initial_state(), f"load #{ev.index} (bus {load.bus}) stepped: P {sign(ev.dp_pct)}, Q {sign(ev.dq_pct)}",
        None,
    )


# --- phase jump --------------------------------------------------------------
def _phase_jump(model: NonlinearNetworkModel, ev: NetworkEvent) -> AppliedEvent:
    net: Network = model.network
    d = math.radians(ev.angle_deg)
    x_eq = model.initial_state()
    x0 = x_eq.copy()
    by_name = {b.name: b for b in model.blocks}

    def rotate(b, angle: float) -> None:
        s = b.state_off
        v = complex(x0[s], x0[s + 1]) * cmath.exp(1j * angle)
        x0[s], x0[s + 1] = v.real, v.imag

    ib = next((u for u in net.der_units if u.bus == ev.bus and u.unit_type.value == "infinite_bus"), None)
    node = by_name.get(f"Nd_{ev.bus}")
    if node is not None:
        rotate(node, d)
        what = f"phase jump of {ev.angle_deg:+g}° at bus {ev.bus}"
    elif ib is not None:
        ib_block = by_name[block_labels(net).der[ib.id]]
        if ib_block.kind == "ib":
            # The source carries its own angle: jump it, and nothing else.
            x0[ib_block.state_off + ib_block.comp.state_names.index("theta")] += d
        else:
            # MATLAB-compatible frame: this source *is* the frame, so the jump
            # is expressed the other way round -- everything else turns by -d.
            for b in model.blocks:
                if b.kind in _COMMON_FRAME_KINDS:
                    rotate(b, -d)
                elif b.kind in _ANGLE_STATE:
                    x0[b.state_off + b.comp.state_names.index(_ANGLE_STATE[b.kind])] -= d
        what = f"phase jump of {ev.angle_deg:+g}° of the infinite-bus source (bus {ev.bus})"
    else:
        raise EventError(
            f"bus {ev.bus} isn't a node of the dynamic model -- pick a network bus, or the infinite bus's own "
            "bus (a unit's terminal bus is inside the unit's model)"
        )
    # The model itself is unchanged: only its initial state moves.
    return AppliedEvent(model, x0, what, x0 - x_eq)


# --- rebuilding -------------------------------------------------------------
def rebuild(model: NonlinearNetworkModel, network: Network, replace: dict | None = None) -> NonlinearNetworkModel:
    """A model of ``network`` (``model.network`` or a sub-network of it with
    the same labels) from ``model``'s own component instances -- by block
    name, with ``replace`` overriding some -- starting from ``model``'s
    operating-point state and algebraic solution, block by block."""
    replace = replace or {}
    comps = {b.name: b.comp for b in model.blocks}
    comps.update(replace)
    lab = block_labels(network)
    # The reference frames of the network the event leaves behind: an island
    # split in two gets one each, and an island whose reference unit was just
    # tripped hands its frame to the unit that takes over. They start where
    # the frame they replace is now, so nothing jumps at the event.
    old_frames = [b for b in model.blocks if b.kind == "frame"]
    x_now = model.initial_state()
    theta_now = float(x_now[old_frames[0].state_off]) if old_frames else model.theta_g0
    frame_components = None
    if old_frames:
        wb_val = 2 * math.pi * network.f_hz
        frame_components = {
            ref: nonlinear_frame_block(wb_val=wb_val, theta0=theta_now, driven=driven)
            for ref, driven in frame_references(network).items()
        }
    blocks, wiring = build_blocks_and_wiring(
        network,
        der_components={d.id: comps[lab.der[d.id]] for d in network.der_units},
        node_components={b.id: comps[f"Nd_{b.id}"] for b in network.buses if f"Nd_{b.id}" in comps},
        line_components=[comps[f"Ln_{n + 1}"] for n in lab.line],
        load_components=[comps[f"Ld_{n + 1}"] for n in lab.load],
        frame_components=frame_components,
    )
    topology = compute_topology(blocks, wiring)
    z_offsets, off = [], 0
    for b in blocks:
        z_offsets.append(off)
        off += b.comp.n_z
    new = NonlinearNetworkModel(
        blocks=blocks, topology=topology, z_offsets=z_offsets, n_z=off, network=network,
        theta_g0=model.theta_g0, op=model.op,
    )

    # Carry the state and the algebraic solution over by block name. A block
    # the event creates -- a frame for an island that didn't exist before --
    # starts from its own initial values instead.
    x_pre = x_now
    z_pre, u_pre = model.solve_algebraic(x_pre, model.default_u_exo())
    old = {b.name: (b, zo) for b, zo in zip(model.blocks, model.z_offsets)}
    x0 = np.zeros(sum(b.comp.n_states for b in blocks))
    z0 = np.zeros(off)
    u0 = np.zeros(topology.n_u)
    for b, zo in zip(blocks, z_offsets):
        if b.name not in old:
            x0[b.state_off:b.state_off + b.comp.n_states] = b.comp.x0
            z0[zo:zo + b.comp.n_z] = b.comp.z0
            u0[b.input_off:b.input_off + b.comp.n_us + b.comp.n_ug] = b.comp.u0
            continue
        ob, ozo = old[b.name]
        x0[b.state_off:b.state_off + b.comp.n_states] = x_pre[ob.state_off:ob.state_off + ob.comp.n_states]
        z0[zo:zo + b.comp.n_z] = z_pre[ozo:ozo + ob.comp.n_z]
        n_u = b.comp.n_us + b.comp.n_ug
        u0[b.input_off:b.input_off + n_u] = u_pre[ob.input_off:ob.input_off + n_u]
    new.x_init = x0
    new.zu_guess = (z0, u0)
    return new
