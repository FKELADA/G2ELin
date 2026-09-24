"""Structural pre-flight checks for a :class:`~g2elin_core.network.schema.Network`
-- catches problems that would otherwise surface one at a time as a crash
(or a confusing wrong answer) deep inside power flow, modal analysis, or
EMT simulation, and reports *every* issue found in one pass instead of
just the first.

Every hand-crafted preset in this codebase satisfies all of these by
construction; a network assembled by hand (the web UI's network editor /
drag-and-drop builder) can easily violate any of them. pydantic's own
validators on ``Network`` already guarantee two structural properties
unconditionally (every bus reference resolves, exactly one slack DER) --
this module covers everything else: things that are *syntactically* valid
``Network`` JSON but still make power flow, modal analysis, or EMT
simulation fail, or silently misbehave, downstream.
"""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx

from g2elin_core import reduction

from .schema import DerUnit, Network, UnitType

# Only an SM or IB slack is wired up in interconnect/network_assembly.py --
# power flow doesn't care (pandapower's ext_grid is agnostic to unit
# type), so this is a *warning* for modal/EMT, not a power flow error.
_SUPPORTED_SLACK_UNIT_TYPES = {"sm", "infinite_bus"}
# Unit types operating_point.py actually derives a dynamics operating point
# for -- these are the ones that need their own transformer.
_DYNAMICS_UNIT_TYPES = {"sm", "gfm", "gfl", "infinite_bus"}


@dataclass(frozen=True)
class NetworkIssue:
    severity: str  # "error" (blocks the affected capabilities) | "warning" (a heads-up, not blocking)
    message: str
    affects: tuple[str, ...]  # which of "powerflow"/"modal"/"emt" this issue affects


def validate_network(network: Network) -> list[NetworkIssue]:
    """Every structural problem found, not just the first. An empty list
    means the network is safe to run power flow, modal analysis, and EMT
    simulation on -- modulo the actual numerics still converging for a
    given operating point, which no static check can guarantee.
    """
    issues: list[NetworkIssue] = []
    all_caps = ("powerflow", "modal", "emt")
    dynamics_caps = ("modal", "emt")

    bus_ids_seen: set[int] = set()
    for bus in network.buses:
        if bus.id in bus_ids_seen:
            issues.append(NetworkIssue("error", f"bus id {bus.id} is used more than once", all_caps))
        bus_ids_seen.add(bus.id)

    der_by_bus: dict[int, DerUnit] = {}
    der_ids_seen: set[int] = set()
    for der in network.der_units:
        if der.id in der_ids_seen:
            issues.append(NetworkIssue("error", f"DER unit id {der.id} is used more than once", all_caps))
        der_ids_seen.add(der.id)

        if der.bus in der_by_bus:
            other = der_by_bus[der.bus]
            issues.append(NetworkIssue(
                "error",
                f"DER units id={der.id} and id={other.id} are both on bus {der.bus} -- each DER unit "
                "needs its own private bus",
                all_caps,
            ))
        else:
            der_by_bus[der.bus] = der

        if der.unit_type == UnitType.NONE:
            issues.append(NetworkIssue(
                "error", f"DER unit id={der.id} has unit_type 'none' -- give it a real type or remove it",
                all_caps,
            ))

    for ln in network.lines:
        if ln.from_bus == ln.to_bus:
            issues.append(NetworkIssue(
                "error", f"a line has from_bus == to_bus == {ln.from_bus} (self-loop)", all_caps,
            ))
        for bus_id, end in ((ln.from_bus, "from_bus"), (ln.to_bus, "to_bus")):
            if bus_id in der_by_bus:
                d = der_by_bus[bus_id]
                issues.append(NetworkIssue(
                    "error",
                    f"a line's {end} (bus {bus_id}) is DER unit id={d.id} ({d.unit_type.value})'s own "
                    "bus -- a DER's own bus may only be reached through its own transformer, never a "
                    "line directly",
                    dynamics_caps,
                ))

    for tr in network.transformers:
        if tr.hv_bus == tr.lv_bus:
            issues.append(NetworkIssue(
                "error", f"a transformer has hv_bus == lv_bus == {tr.hv_bus} (self-loop)", all_caps,
            ))

    for idx, ld in enumerate(network.loads):
        if ld.bus in der_by_bus:
            d = der_by_bus[ld.bus]
            issues.append(NetworkIssue(
                "error",
                f"load #{idx} sits on bus {ld.bus}, DER unit id={d.id} ({d.unit_type.value})'s own bus "
                "-- a DER's own local consumption belongs on the unit itself (its p_cons_mw/q_cons_mvar "
                "fields), not as a separate Load on its private bus",
                dynamics_caps,
            ))
        if ld.q_mvar == 0:
            # Power flow treats this as a no-op or purely-resistive bus (both
            # perfectly valid), but operating_point.compute_operating_point()'s
            # constant-impedance equivalent reactance x_pu = z_pu*sin(acos(PF))
            # is exactly zero whenever q_mvar=0 -- regardless of p_mw -- and
            # components/load.py's own dynamic model divides by that reactance
            # (subsumes the old p_mw=q_mvar=0 "zero apparent power" case too,
            # which is undefined for the same underlying reason).
            issues.append(NetworkIssue(
                "error",
                f"load #{idx} (bus {ld.bus}) has zero reactive power (q_mvar=0) -- give it a small "
                "nonzero q_mvar (positive for inductive, negative for capacitive)",
                dynamics_caps,
            ))

    if not network.lines and network.der_units:
        # pipeline.linearize_network()'s "node quirk": every bus's own dynamic
        # model borrows its shunt (line-charging) capacitance from the *first*
        # Line's b_pu -- with zero Lines in the network there's nothing to
        # borrow, and the 0.0 fallback isn't actually safe: it makes every
        # node's own state equation (dv/dt = (wb/Cl)*(...)) divide by zero.
        # A transformer-only network (e.g. one DER + one transformer + one
        # load bus, no feeder) hits this even though power flow is fine with it.
        issues.append(NetworkIssue(
            "error",
            "this network has no Line elements -- every bus's own dynamic model needs a line-charging "
            "susceptance (b_pu) to linearize around, which this codebase always borrows from the first "
            "Line in the network; a network built entirely from transformers has no such source and "
            "can't run modal analysis or EMT. Add at least one Line (even a short one with a small b_pu)",
            dynamics_caps,
        ))

    if not network.transformers:
        if network.der_units:
            issues.append(NetworkIssue(
                "error",
                "this network has no transformers -- every SM/GFM/GFL/IB unit must sit behind its own "
                "transformer (a Transformer whose lv_bus is that unit's own bus)",
                dynamics_caps,
            ))
    else:
        transformer_by_lv_bus = {tr.lv_bus: tr for tr in network.transformers}
        transformer_by_hv_bus = {tr.hv_bus: tr for tr in network.transformers}
        for der in network.der_units:
            if der.unit_type.value not in _DYNAMICS_UNIT_TYPES:
                continue
            tr = transformer_by_lv_bus.get(der.bus)
            if tr is None:
                swapped = transformer_by_hv_bus.get(der.bus)
                if swapped is not None:
                    # The user did add a transformer between this DER and the rest of
                    # the network -- it's just wired backward. hv_bus/lv_bus direction
                    # matters here (see interconnect/network_assembly.py's own
                    # transformer_by_lv_bus lookup), not just a label: the DER's own
                    # bus must be lv_bus, the network-side bus must be hv_bus.
                    issues.append(NetworkIssue(
                        "error",
                        f"DER unit id={der.id} ({der.unit_type.value})'s transformer to bus {der.bus} has "
                        f"hv_bus/lv_bus swapped (currently hv_bus={swapped.hv_bus}, lv_bus={swapped.lv_bus}) "
                        f"-- a DER's own bus must be the transformer's lv_bus, not its hv_bus; swap the two",
                        dynamics_caps,
                    ))
                else:
                    issues.append(NetworkIssue(
                        "error",
                        f"DER unit id={der.id} ({der.unit_type.value}) at bus {der.bus} has no transformer "
                        f"connecting it to the rest of the network (no Transformer has lv_bus={der.bus})",
                        dynamics_caps,
                    ))
            elif tr.hv_bus in der_by_bus and der_by_bus[tr.hv_bus].id != der.id:
                other = der_by_bus[tr.hv_bus]
                issues.append(NetworkIssue(
                    "error",
                    f"DER unit id={der.id} ({der.unit_type.value})'s own transformer connects directly "
                    f"to DER unit id={other.id} ({other.unit_type.value})'s bus ({tr.hv_bus}) instead of "
                    "a plain grid bus",
                    dynamics_caps,
                ))

    # Parameter overrides must name a real parameter of that unit type -- a
    # typo would otherwise be silently ignored by the model.
    # Two transformers claiming the same LV bus: which grid bus a unit hangs
    # from would then be decided by list order, silently.
    lv_seen: dict[int, int] = {}
    unit_buses = {d.bus for d in network.der_units}
    for ti, tr in enumerate(network.transformers):
        if tr.lv_bus in unit_buses:
            if tr.lv_bus in lv_seen:
                issues.append(NetworkIssue(
                    "error",
                    f"transformers #{lv_seen[tr.lv_bus]} and #{ti} both have their LV side on bus "
                    f"{tr.lv_bus}, which carries a unit -- that unit's step-up has to be unambiguous, "
                    f"since its impedance goes inside the unit's own model and decides which grid bus "
                    f"the unit injects into",
                    all_caps,
                ))
            lv_seen[tr.lv_bus] = ti
        if tr.hv_bus == tr.lv_bus:
            issues.append(NetworkIssue("error", f"transformer #{ti} has both sides on bus {tr.hv_bus}", all_caps))

    for si, sh in enumerate(network.shunts):
        if sh.q_mvar == 0:
            issues.append(NetworkIssue(
                "warning", f"shunt #{si} (bus {sh.bus}) is rated 0 MVAr, so it does nothing", all_caps,
            ))
        if sh.bus in unit_buses:
            issues.append(NetworkIssue(
                "error",
                f"shunt #{si} sits on bus {sh.bus}, which is a unit's own terminal bus -- that bus is "
                f"inside the unit's model and has no node of its own to attach to. Put it on the grid "
                f"bus on the other side of the unit's transformer.",
                dynamics_caps,
            ))

    from g2elin_core.operating_point import REBASED_PARAM_KEYS, unit_keys  # lazy: it imports this module

    for der in network.der_units:
        if der.sn_mva is not None and not der.params:
            issues.append(NetworkIssue(
                "warning",
                f"unit id={der.id} declares sn_mva={der.sn_mva:g} MVA but overrides no parameters, so the "
                f"rating has no effect -- it states the base of `params`, and the built-in defaults are "
                f"already on the network base",
                dynamics_caps,
            ))
        if not der.params:
            continue
        unknown = sorted(set(der.params) - unit_keys(der))
        if unknown:
            issues.append(NetworkIssue(
                "error",
                f"unit id={der.id} ({der.unit_type.value}) overrides parameter(s) it doesn't have: {unknown}",
                dynamics_caps,
            ))
        # Gains and time constants are tuning, not machine data, so they are
        # taken as given on the network base (see operating_point.rebase_params).
        if der.sn_mva is not None and der.sn_mva != network.sn_mva:
            not_rebased = sorted(set(der.params) - set(unknown) - REBASED_PARAM_KEYS)
            if not_rebased:
                issues.append(NetworkIssue(
                    "warning",
                    f"unit id={der.id} is rated {der.sn_mva:g} MVA on a {network.sn_mva:g} MVA network, but "
                    f"{not_rebased} are gains or time constants, which are taken as given rather than "
                    f"converted between bases -- give them on the network base",
                    dynamics_caps,
                ))

    # Outside the MATLAB-compatible mode, a unit's Rt/Lt *are* its transformer
    # (see operating_point.unit_transformer_rx); an override would make the
    # dynamic model use a different impedance from the power flow's.
    if not network.units_use_first_transformer:
        for der in network.der_units:
            clash = sorted({"Rt", "Lt"} & set(der.params))
            if clash:
                issues.append(NetworkIssue(
                    "warning",
                    f"unit id={der.id} overrides {'/'.join(clash)}, which normally follow its own transformer -- "
                    "the power flow still uses the transformer's impedance, so the operating point and the "
                    "dynamic model disagree. Edit the transformer instead (or remove the override)",
                    dynamics_caps,
                ))

    # Only the MATLAB-compatible frame cares what the slack unit is: it *is*
    # the frame there. With a frame of its own any unit can be the slack.
    slack = next((d for d in network.der_units if d.bus_type.value == "slack"), None)
    if network.frame_follows_slack and slack is not None and slack.unit_type.value not in _SUPPORTED_SLACK_UNIT_TYPES:
        issues.append(NetworkIssue(
            "warning",
            f"the slack unit (id={slack.id}) is a {slack.unit_type.value}, not a synchronous machine or "
            f"infinite bus -- this network asks the dynamic models' frame to follow the slack "
            "(units_use_first_transformer's sibling, frame_follows_slack), which only those two can do; "
            "switch that off to use a frame of its own",
            dynamics_caps,
        ))

    # Connectivity: every bus must be reachable from the slack through
    # lines/transformers, or whichever part isn't has no reference bus for
    # power flow to solve against.
    graph = nx.Graph()
    graph.add_nodes_from(b.id for b in network.buses)
    graph.add_edges_from((ln.from_bus, ln.to_bus) for ln in network.lines)
    graph.add_edges_from((tr.hv_bus, tr.lv_bus) for tr in network.transformers)
    if network.buses and not nx.is_connected(graph):
        sizes = sorted((len(c) for c in nx.connected_components(graph)), reverse=True)
        issues.append(NetworkIssue(
            "error",
            f"the network isn't fully connected ({len(sizes)} separate groups of buses, sizes {sizes}) "
            "-- every bus must be reachable from the slack through lines/transformers",
            all_caps,
        ))

    # Breakers (network.breakers): what open breakers leave out, and which
    # unit is the reference of each island that is still energized.
    from .breakers import out_of_service_summary, service_state  # local: keeps this module's imports minimal

    st = service_state(network)
    if network.der_units and not st.references:
        issues.append(NetworkIssue(
            "error",
            "no unit is left that could set a voltage and a frequency: open breakers have taken every "
            "synchronous machine, grid-forming converter and infinite bus out of service (a grid-following "
            "converter can only follow a voltage, never start one)",
            all_caps,
        ))
    else:
        summary = out_of_service_summary(network)
        if summary:
            issues.append(NetworkIssue(
                "warning",
                "open breakers: " + "; ".join(summary) + " -- every analysis runs without them",
                all_caps,
            ))
        by_id = {d.id: d for d in network.der_units}
        if len(st.references) > 1:
            parts = [
                f"unit {i.reference} ({by_id[i.reference].unit_type.value}) for bus(es) {sorted(i.buses)}"
                for i in st.islands if i.reference is not None
            ]
            issues.append(NetworkIssue(
                "warning",
                f"open breakers have split the network into {len(st.references)} energized islands, each "
                "solved against its own reference: " + "; ".join(parts),
                all_caps,
            ))
        elif st.references and slack is not None and st.references[0] != slack.id:
            ref = st.references[0]
            issues.append(NetworkIssue(
                "warning",
                f"the designated slack (unit {slack.id}) is out of service -- unit {ref} "
                f"({by_id[ref].unit_type.value}) is the reference the power flow solves against instead",
                all_caps,
            ))

    issues.extend(_model_order_issues(network))
    return issues


def _model_order_issues(network: Network) -> list[NetworkIssue]:
    """Combinations of model levels that are mathematically buildable but
    physically inconsistent -- reported as warnings, never errors.

    They are warnings on purpose. A mixed-timescale model is a legitimate
    thing to ask for (studying what one element's fast dynamics contribute
    is exactly what this control gives you), and none of these combinations
    makes the model unsolvable. What they do is give an answer whose
    timescale is not the one the user probably thinks they asked for, which
    is precisely what a warning is for.
    """
    issues: list[NetworkIssue] = []
    dynamics_caps = ("modal", "emt")
    net_modes = network.models.group_modes_for("network")
    net_dynamic = [gid for gid, mode in net_modes.items() if mode == reduction.DYNAMIC]
    quasi_stationary = not net_dynamic

    if net_modes.get("nodes") == reduction.ALGEBRAIC and network.nodes_share_first_line_b:
        # A quasi-stationary bus equation says the current injected into the
        # bus equals its own shunt's, jB*V. That is the power flow's own bus
        # equation -- but only when B is the bus's real charging. With the
        # MATLAB-compatible option every bus borrows the *first line's*
        # susceptance instead, which no power flow ever saw, so the bus
        # voltages the model settles at are not the ones the power flow
        # solved for. (Measured on WSCC-9: bus voltages land 0.03 pu away
        # and the units start ~3000x further from equilibrium than with each
        # bus's own capacitance.) A dynamic network hides this as a fast
        # transient that decays in microseconds; a quasi-stationary one
        # cannot, so it becomes a visible start-up excursion.
        issues.append(NetworkIssue(
            "warning",
            "this network shares the first line's charging susceptance across every bus "
            "(nodes_share_first_line_b) while the network model is quasi-stationary. The two are "
            "inconsistent: the bus equation then demands a shunt current the power flow never "
            "solved for, so a time-domain run starts with a visible excursion. Switch that option "
            "off (each bus gets its own capacitance) for quasi-stationary runs, or keep the "
            "network dynamic.",
            dynamics_caps,
        ))

    if quasi_stationary:
        for der in network.der_units:
            if der.unit_type.value not in reduction.ELEMENTS:
                continue
            groups = network.unit_group_modes(der)
            kind = der.unit_type.value
            if groups.get("trafo_current") == reduction.DYNAMIC:
                issues.append(NetworkIssue(
                    "warning",
                    f"unit {der.id} ({kind}) still integrates its transformer current while the network "
                    "is quasi-stationary: the unit keeps an electromagnetic transient the grid it feeds "
                    "no longer has. Use a level of 'order8' or lower (machines) / 'no_trafo' or lower "
                    "(converters), or switch the network back to full dynamics.",
                    dynamics_caps,
                ))
            if kind == "sm" and groups.get("stator_flux") == reduction.DYNAMIC:
                issues.append(NetworkIssue(
                    "warning",
                    f"unit {der.id} keeps its stator flux transients against a quasi-stationary network. "
                    "The standard pairing is to drop both together (machine level 'order6' or lower); "
                    "keeping only the machine's produces modes with no physical counterpart.",
                    dynamics_caps,
                ))
            if kind in ("gfm", "gfl") and groups.get("filter") == reduction.DYNAMIC:
                issues.append(NetworkIssue(
                    "warning",
                    f"unit {der.id} ({kind}) keeps its filter dynamics against a quasi-stationary "
                    "network -- a mismatched pair, the same way stator flux is for a machine.",
                    dynamics_caps,
                ))
    elif len(net_dynamic) != len(net_modes):
        algebraic = sorted(set(net_modes) - set(net_dynamic))
        issues.append(NetworkIssue(
            "warning",
            "the network is split across timescales: "
            f"{', '.join(algebraic)} are quasi-stationary while {', '.join(sorted(net_dynamic))} still "
            "integrate. That is a valid model, but its results are neither EMT nor RMS.",
            dynamics_caps,
        ))

    # A unit whose per-group settings don't match any named level: worth
    # saying out loud, because the UI's level picker will show "custom" and
    # the combination may not be one anybody has validated.
    for der in network.der_units:
        if der.unit_type.value not in reduction.ELEMENTS:
            continue
        if network.unit_level(der) is None:
            issues.append(NetworkIssue(
                "warning",
                f"unit {der.id}'s state selection doesn't match any named model level -- a custom "
                "combination. Check it against the full-order model before trusting it (the model "
                "adequacy check does exactly that).",
                dynamics_caps,
            ))

    return issues
