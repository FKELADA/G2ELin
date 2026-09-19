"""Breakers: which elements are in service, and the network the analyses see.

Every line and transformer has a breaker at each end; every load and unit
has one between it and its bus (``Line.from_closed``/``to_closed``,
``Transformer.hv_closed``/``lv_closed``, ``Load.closed``,
``DerUnit.closed`` -- all closed by default).

Rules, shared by every analysis:

- A line or transformer with *either* breaker open is out of service. (A
  line open at one end only would, in reality, still charge from the other
  end; that charging current is neglected, as it is by the lumped model.)
- A load or unit with its breaker open is out of service.
- A bus is **energized** when it is connected to the slack unit's bus
  through in-service lines/transformers. Anything on a de-energized bus
  (an island cut off from the slack) is out of service too -- the tool has
  one reference (the slack) and doesn't simulate islanded operation from a
  steady state.
- The slack unit's own breaker can't be opened: it is the reference of the
  power flow and of the dynamic models' common frame.

Power flow keeps every element and marks these ones out of service
(pandapower ``in_service``), so its result tables keep the network's own
element numbering. The dynamic models (modal analysis, root locus, EMT) are
built from :func:`energized_network` -- the network with those elements
removed -- whose model names stay those of the full network (see
:class:`BlockLabels`), so e.g. line #3's states are always ``..._{Ln_4}``.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

from .schema import Network

TYPE_LABEL = {"sm": "SM", "gfm": "GFM", "gfl": "GFL", "infinite_bus": "IB"}


@dataclass(frozen=True)
class BlockLabels:
    """Names/numbers the dynamic model gives a network's elements.

    ``der``: unit id -> block name ("SM_2", "GFL_1", ...: per-type counters
    in declaration order). ``line``/``load``/``transformer``: for each
    element of *this* network, its 0-based index in the network it was
    taken from. ``node_b_pu``: the susceptance every node borrows from the
    first line (see ``pipeline.linearize_network``'s node quirk), taken from
    the full network so it doesn't change when line #1 is switched out.
    """

    der: dict[int, str]
    line: list[int]
    load: list[int]
    transformer: list[int]
    node_b_pu: float | None


def default_labels(network: Network) -> BlockLabels:
    counters: dict[str, int] = {}
    der: dict[int, str] = {}
    for d in network.der_units:
        label = TYPE_LABEL.get(d.unit_type.value, d.unit_type.value.upper())
        counters[label] = counters.get(label, 0) + 1
        der[d.id] = f"{label}_{counters[label]}"
    return BlockLabels(
        der=der,
        line=list(range(len(network.lines))),
        load=list(range(len(network.loads))),
        transformer=list(range(len(network.transformers))),
        node_b_pu=network.lines[0].b_pu if network.lines else None,
    )


def block_labels(network: Network) -> BlockLabels:
    """The labels a reduced network carries, else the network's own."""
    labels = getattr(network, "_labels", None)
    return labels if isinstance(labels, BlockLabels) else default_labels(network)


@dataclass(frozen=True)
class ServiceState:
    """Which elements of a network are in service (indices into its lists)."""

    slack_connected: bool
    energized_buses: frozenset[int]
    lines: tuple[bool, ...]
    transformers: tuple[bool, ...]
    loads: tuple[bool, ...]
    der_units: dict[int, bool] = field(default_factory=dict)

    @property
    def everything_in_service(self) -> bool:
        return (
            self.slack_connected and all(self.lines) and all(self.transformers)
            and all(self.loads) and all(self.der_units.values())
        )


def any_breaker_open(network: Network) -> bool:
    return (
        any(not (ln.from_closed and ln.to_closed) for ln in network.lines)
        or any(not (tr.hv_closed and tr.lv_closed) for tr in network.transformers)
        or any(not ld.closed for ld in network.loads)
        or any(not d.closed for d in network.der_units)
    )


def service_state(network: Network) -> ServiceState:
    slack = next(d for d in network.der_units if d.bus_type.value == "slack")
    adj: dict[int, list[int]] = {b.id: [] for b in network.buses}
    for ln in network.lines:
        if ln.from_closed and ln.to_closed:
            adj[ln.from_bus].append(ln.to_bus)
            adj[ln.to_bus].append(ln.from_bus)
    for tr in network.transformers:
        if tr.hv_closed and tr.lv_closed:
            adj[tr.hv_bus].append(tr.lv_bus)
            adj[tr.lv_bus].append(tr.hv_bus)

    energized: set[int] = set()
    if slack.closed:
        queue = deque([slack.bus])
        energized.add(slack.bus)
        while queue:
            for nxt in adj[queue.popleft()]:
                if nxt not in energized:
                    energized.add(nxt)
                    queue.append(nxt)

    return ServiceState(
        slack_connected=slack.closed,
        energized_buses=frozenset(energized),
        lines=tuple(ln.from_closed and ln.to_closed and ln.from_bus in energized for ln in network.lines),
        transformers=tuple(tr.hv_closed and tr.lv_closed and tr.hv_bus in energized for tr in network.transformers),
        loads=tuple(ld.closed and ld.bus in energized for ld in network.loads),
        der_units={d.id: d.closed and d.bus in energized for d in network.der_units},
    )


class SlackDisconnected(ValueError):
    pass


def energized_network(network: Network) -> Network:
    """``network`` without its out-of-service elements (see the module
    docstring), for the dynamic models. Returns ``network`` itself when
    every breaker is closed.

    A unit that is out of service takes its own transformer and terminal
    bus with it (the dynamic model has the transformer inside the unit).
    """
    if not any_breaker_open(network):
        return network
    st = service_state(network)
    if not st.slack_connected:
        raise SlackDisconnected(
            "the slack unit's breaker is open -- it is the reference of the power flow and of the dynamic "
            "models, so it can't be disconnected; make another unit the slack first"
        )
    full = block_labels(network)
    unit_bus_out = {d.bus for d in network.der_units if not st.der_units[d.id]}
    keep_tr = [
        j for j, tr in enumerate(network.transformers) if st.transformers[j] and tr.lv_bus not in unit_bus_out
    ]
    keep_ln = [i for i, ok in enumerate(st.lines) if ok]
    keep_ld = [i for i, ok in enumerate(st.loads) if ok]
    keep_der = [d for d in network.der_units if st.der_units[d.id]]
    buses = [b for b in network.buses if b.id in st.energized_buses and b.id not in unit_bus_out]

    reduced = network.model_copy(update=dict(
        buses=buses,
        lines=[network.lines[i] for i in keep_ln],
        transformers=[network.transformers[j] for j in keep_tr],
        loads=[network.loads[i] for i in keep_ld],
        der_units=keep_der,
    ))
    reduced._labels = BlockLabels(
        der={d.id: full.der[d.id] for d in keep_der},
        line=[full.line[i] for i in keep_ln],
        load=[full.load[i] for i in keep_ld],
        transformer=[full.transformer[j] for j in keep_tr],
        node_b_pu=full.node_b_pu,
    )
    return reduced


def node_b_pu(network: Network) -> float | None:
    """The susceptance every node uses (the full network's first line's)."""
    return block_labels(network).node_b_pu


def out_of_service_summary(network: Network) -> list[str]:
    """Human-readable list of what is out of service and why (empty when
    everything is in service)."""
    if not any_breaker_open(network):
        return []
    st = service_state(network)
    msgs: list[str] = []
    dead = sorted(b.id for b in network.buses if b.id not in st.energized_buses)
    if dead:
        msgs.append(f"de-energized bus(es) {dead} (cut off from the slack by open breakers)")
    lines = [i for i, ok in enumerate(st.lines) if not ok]
    if lines:
        msgs.append(f"line(s) #{', #'.join(map(str, lines))} out of service")
    trs = [j for j, ok in enumerate(st.transformers) if not ok]
    if trs:
        msgs.append(f"transformer(s) #{', #'.join(map(str, trs))} out of service")
    lds = [i for i, ok in enumerate(st.loads) if not ok]
    if lds:
        msgs.append(f"load(s) #{', #'.join(map(str, lds))} disconnected")
    ders = [i for i, ok in st.der_units.items() if not ok]
    if ders:
        msgs.append(f"unit(s) id {', '.join(map(str, ders))} disconnected")
    return msgs
