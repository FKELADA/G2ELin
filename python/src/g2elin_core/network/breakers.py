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
- Open breakers split the network into **islands**. An island is energized
  when it still holds a unit that can set its own voltage and frequency --
  a synchronous machine, a grid-forming converter or an infinite bus. That
  unit is the island's **reference** (the power flow's slack for it): the
  network's designated slack keeps the role in its own island, otherwise
  the largest grid former takes it, infinite bus first, then synchronous
  machine, then grid-forming converter.
- An island holding only grid-following converters and loads is **blacked
  out**: a grid-following converter needs a voltage to follow and cannot
  start one, and anti-islanding protection would trip it. Everything in
  such an island is out of service, as it is in every load-flow tool that
  requires a reference per island.
- Any breaker may be opened, the slack unit's included -- the reference
  then moves to another unit. Only a network left with no grid former at
  all has nothing to solve.

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

from .schema import BusType, Network

TYPE_LABEL = {"sm": "SM", "gfm": "GFM", "gfl": "GFL", "infinite_bus": "IB"}
# Units that can hold an island's voltage and frequency on their own, in the
# order they are preferred as its reference.
GRID_FORMING = ("infinite_bus", "sm", "gfm")


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
    shunt: list[int]
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
        shunt=list(range(len(network.shunts))),
        node_b_pu=network.lines[0].b_pu if network.lines else None,
    )


def block_labels(network: Network) -> BlockLabels:
    """The labels a reduced network carries, else the network's own."""
    labels = getattr(network, "_labels", None)
    return labels if isinstance(labels, BlockLabels) else default_labels(network)


@dataclass(frozen=True)
class Island:
    """One electrically connected group of buses, and the unit that is its
    power-flow reference (``None`` when it has no grid former: blacked out)."""

    buses: frozenset[int]
    reference: int | None


@dataclass(frozen=True)
class ServiceState:
    """Which elements of a network are in service (indices into its lists)."""

    islands: tuple[Island, ...]
    energized_buses: frozenset[int]
    lines: tuple[bool, ...]
    transformers: tuple[bool, ...]
    loads: tuple[bool, ...]
    shunts: tuple[bool, ...] = ()
    der_units: dict[int, bool] = field(default_factory=dict)

    @property
    def references(self) -> tuple[int, ...]:
        """The units acting as a power-flow reference, one per energized island."""
        return tuple(i.reference for i in self.islands if i.reference is not None)

    @property
    def primary_reference(self) -> int | None:
        """The reference of the largest energized island -- the one that takes
        the ``slack`` role when the designated slack unit is out of service."""
        live = [i for i in self.islands if i.reference is not None]
        return max(live, key=lambda i: (len(i.buses), -i.reference)).reference if live else None

    @property
    def everything_in_service(self) -> bool:
        return (
            all(self.lines) and all(self.transformers)
            and all(self.loads) and all(self.der_units.values())
        )


def any_breaker_open(network: Network) -> bool:
    return (
        any(not (ln.from_closed and ln.to_closed) for ln in network.lines)
        or any(not (tr.hv_closed and tr.lv_closed) for tr in network.transformers)
        or any(not ld.closed for ld in network.loads)
        or any(not sh.closed for sh in network.shunts)
        or any(not d.closed for d in network.der_units)
    )


def service_state(network: Network) -> ServiceState:
    adj: dict[int, list[int]] = {b.id: [] for b in network.buses}
    for ln in network.lines:
        if ln.from_closed and ln.to_closed:
            adj[ln.from_bus].append(ln.to_bus)
            adj[ln.to_bus].append(ln.from_bus)
    for tr in network.transformers:
        if tr.hv_closed and tr.lv_closed:
            adj[tr.hv_bus].append(tr.lv_bus)
            adj[tr.lv_bus].append(tr.hv_bus)

    # Islands: the connected groups of buses left by the closed breakers.
    seen: set[int] = set()
    groups: list[frozenset[int]] = []
    for bus in network.buses:
        if bus.id in seen:
            continue
        group, queue = {bus.id}, deque([bus.id])
        seen.add(bus.id)
        while queue:
            for nxt in adj[queue.popleft()]:
                if nxt not in group:
                    group.add(nxt)
                    seen.add(nxt)
                    queue.append(nxt)
        groups.append(frozenset(group))

    def rank(d) -> tuple:
        # The designated slack first, then by unit type, then the largest.
        return (d.bus_type is not BusType.SLACK, GRID_FORMING.index(d.unit_type.value), -abs(d.p_set_mw), d.id)

    islands: list[Island] = []
    for group in groups:
        formers = [
            d for d in network.der_units
            if d.closed and d.bus in group and d.unit_type.value in GRID_FORMING
        ]
        ref = min(formers, key=rank) if formers else None
        islands.append(Island(buses=group, reference=ref.id if ref else None))
    energized = {b for i in islands if i.reference is not None for b in i.buses}

    return ServiceState(
        islands=tuple(islands),
        energized_buses=frozenset(energized),
        lines=tuple(ln.from_closed and ln.to_closed and ln.from_bus in energized for ln in network.lines),
        transformers=tuple(tr.hv_closed and tr.lv_closed and tr.hv_bus in energized for tr in network.transformers),
        loads=tuple(ld.closed and ld.bus in energized for ld in network.loads),
        shunts=tuple(sh.closed and sh.bus in energized for sh in network.shunts),
        der_units={d.id: d.closed and d.bus in energized for d in network.der_units},
    )


class NoReferenceUnit(ValueError):
    """Nothing is left that could set a voltage and a frequency."""


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
    if not st.references:
        raise NoReferenceUnit(
            "no part of this network has a unit that could set its voltage and frequency -- open breakers "
            "have left only grid-following converters and loads, which can't energize an island on their own"
        )
    full = block_labels(network)
    unit_bus_out = {d.bus for d in network.der_units if not st.der_units[d.id]}
    keep_tr = [
        j for j, tr in enumerate(network.transformers) if st.transformers[j] and tr.lv_bus not in unit_bus_out
    ]
    keep_ln = [i for i, ok in enumerate(st.lines) if ok]
    keep_ld = [i for i, ok in enumerate(st.loads) if ok]
    keep_sh = [i for i, ok in enumerate(st.shunts) if ok]
    keep_der = [d for d in network.der_units if st.der_units[d.id]]
    # The reduced network needs exactly one unit marked as the slack (the
    # schema's own rule). It is the designated one while it is in service;
    # otherwise the largest energized island's reference takes the role.
    if not any(d.bus_type is BusType.SLACK for d in keep_der):
        primary = st.primary_reference
        keep_der = [
            d.model_copy(update={"bus_type": BusType.SLACK}) if d.id == primary else d for d in keep_der
        ]
    buses = [b for b in network.buses if b.id in st.energized_buses and b.id not in unit_bus_out]

    reduced = network.model_copy(update=dict(
        buses=buses,
        lines=[network.lines[i] for i in keep_ln],
        transformers=[network.transformers[j] for j in keep_tr],
        loads=[network.loads[i] for i in keep_ld],
        shunts=[network.shunts[i] for i in keep_sh],
        der_units=keep_der,
    ))
    reduced._labels = BlockLabels(
        der={d.id: full.der[d.id] for d in keep_der},
        line=[full.line[i] for i in keep_ln],
        load=[full.load[i] for i in keep_ld],
        transformer=[full.transformer[j] for j in keep_tr],
        shunt=[full.shunt[i] for i in keep_sh],
        node_b_pu=full.node_b_pu,
    )
    return reduced


def frame_references(network: Network) -> dict[int, bool]:
    """``{reference unit id: the frame follows its speed}`` -- one entry per
    energized island (see ``components/frame.py``). An infinite-bus reference
    drives nothing: its source already turns at a fixed speed, so its island's
    frame does too.
    """
    by_id = {d.id: d for d in network.der_units}
    return {
        i.reference: by_id[i.reference].unit_type.value != "infinite_bus"
        for i in service_state(network).islands if i.reference is not None
    }


def node_b_pu(network: Network) -> float | None:
    """The susceptance every node uses (the full network's first line's)."""
    return block_labels(network).node_b_pu


# A shunt reactor's resistance when it declares none: a typical X/R for one.
# Exactly 0 would leave its resonance with the bus capacitance undamped.
DEFAULT_SHUNT_REACTOR_XR = 50.0


def shunt_is_reactor(shunt) -> bool:
    """A reactor absorbs reactive power and needs a branch of its own; a
    capacitor bank generates it and folds into the bus's own capacitance."""
    return shunt.q_mvar > 0


def shunt_reactor_rx(shunt, sn_mva: float) -> tuple[float, float]:
    """``(r_pu, x_pu)`` of a shunt reactor's branch, on the network base.

    ``q_mvar`` is its nameplate, the reactive power it absorbs at nominal
    voltage through its reactance alone (Q = V^2/X). The resistance comes
    from its X/R, so the branch also draws a little active power -- which is
    why the power flow is given the P and Q of *this* R-X pair rather than
    the nameplate, so both sides of the tool see the same device.
    """
    x_pu = sn_mva / shunt.q_mvar
    r_pu = shunt.r_pu if shunt.r_pu is not None else x_pu / DEFAULT_SHUNT_REACTOR_XR
    return r_pu, x_pu


def shunt_reactor_pq_mw(shunt, sn_mva: float) -> tuple[float, float]:
    """``(p_mw, q_mvar)`` a reactor's R-X branch actually draws at 1 pu."""
    r_pu, x_pu = shunt_reactor_rx(shunt, sn_mva)
    denom = r_pu**2 + x_pu**2
    return sn_mva * r_pu / denom, sn_mva * x_pu / denom


def node_capacitances(network: Network) -> dict[int, float]:
    """Each bus's shunt capacitance ``Cl`` in per unit -- what its own dynamic
    model integrates (``components/node.py``: ``C dv/dt = i - jwC v``).

    Physically that is half the charging of every in-service line touching the
    bus (a pi-model puts half at each end) plus any capacitor bank on it. With
    ``Network.nodes_share_first_line_b`` every bus instead borrows the first
    line's charging, which is what the MATLAB toolbox does and what the ported
    presets keep so their numbers still match it.

    Covers the buses that get a node block of their own, which is every bus
    except a unit's own terminal bus -- that one lives inside the unit's model,
    behind its step-up transformer, and has no node equation to divide by.

    Counts every line and shunt the network still lists, open breakers
    included, exactly as the blocks do: dropping what is out of service is
    :func:`energized_network`'s job, and the reduced network it returns no
    longer lists them.

    Raises when a bus ends up with none and ``Network.min_node_b_pu`` says
    nothing, rather than inventing a value to divide by.
    """
    modelled = [b.id for b in network.buses if b.id not in {d.bus for d in network.der_units}]
    if network.nodes_share_first_line_b:
        shared = block_labels(network).node_b_pu
        if shared is None:
            return {}
        return {bus_id: shared for bus_id in modelled}

    cap = {bus_id: 0.0 for bus_id in modelled}
    for ln in network.lines:
        half = ln.b_pu / 2
        for end in (ln.from_bus, ln.to_bus):
            if end in cap:
                cap[end] += half
    for sh in network.shunts:
        if not shunt_is_reactor(sh) and sh.bus in cap:
            cap[sh.bus] += (-sh.q_mvar) / network.sn_mva

    empty = sorted(b for b, c in cap.items() if c <= 0.0)
    if empty:
        if network.min_node_b_pu is None:
            raise ValueError(
                f"bus(es) {empty} have no shunt capacitance: none of their in-service lines declare any "
                f"charging (b_pu) and they carry no capacitor bank. Every bus's dynamic model integrates "
                f"dv/dt = (wb/Cl)*(...), so Cl has to be greater than zero. Either give those lines their "
                f"charging, add a capacitor bank, set Network.min_node_b_pu to the value you want used "
                f"there (distribution-feeder data that omits charging is the usual reason), or set "
                f"Network.nodes_share_first_line_b to reproduce the MATLAB toolbox's single shared value."
            )
        for b in empty:
            cap[b] = network.min_node_b_pu
    return cap


def out_of_service_summary(network: Network) -> list[str]:
    """Human-readable list of what is out of service and why (empty when
    everything is in service)."""
    if not any_breaker_open(network):
        return []
    st = service_state(network)
    msgs: list[str] = []
    dead = sorted(b.id for b in network.buses if b.id not in st.energized_buses)
    if dead:
        blacked = [i for i in st.islands if i.reference is None and i.buses]
        why = ("split off with no unit able to set a voltage and a frequency"
               if any(len(i.buses) > 1 for i in blacked) else "cut off by open breakers")
        msgs.append(f"de-energized bus(es) {dead} ({why})")
    lines = [i for i, ok in enumerate(st.lines) if not ok]
    if lines:
        msgs.append(f"line(s) #{', #'.join(map(str, lines))} out of service")
    trs = [j for j, ok in enumerate(st.transformers) if not ok]
    if trs:
        msgs.append(f"transformer(s) #{', #'.join(map(str, trs))} out of service")
    lds = [i for i, ok in enumerate(st.loads) if not ok]
    if lds:
        msgs.append(f"load(s) #{', #'.join(map(str, lds))} disconnected")
    shs = [i for i, ok in enumerate(st.shunts) if not ok]
    if shs:
        msgs.append(f"shunt(s) #{', #'.join(map(str, shs))} disconnected")
    ders = [i for i, ok in st.der_units.items() if not ok]
    if ders:
        msgs.append(f"unit(s) id {', '.join(map(str, ders))} disconnected")
    return msgs
