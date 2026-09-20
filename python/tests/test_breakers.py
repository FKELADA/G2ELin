"""Breakers (g2elin_core.network.breakers) and EMT network events
(g2elin_core.timedomain.events): an open breaker takes its element out of
every analysis, with the full network's element numbering kept; a network
event applied at t=0 starts from the pre-event state and continues with the
post-event network."""

from __future__ import annotations

import numpy as np
import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core.modal import analyze
from g2elin_core.network.breakers import energized_network, service_state
from g2elin_core.network.presets import gfm_smib, wscc9_1sm_1gfm_1gfl, wscc9_3sm
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow

client = TestClient(app)


def _dump(net) -> dict:
    return net.model_dump(mode="json")


def test_everything_closed_is_unchanged():
    net = wscc9_3sm()
    assert energized_network(net) is net
    assert service_state(net).everything_in_service


def test_open_line_power_flow_keeps_numbering():
    net = wscc9_3sm()
    net.lines[3].to_closed = False
    r = client.post("/api/network/powerflow", json={"network": _dump(net)}).json()
    assert r["converged"]
    assert [row["in_service"] for row in r["lines"]] == [True, True, True, False, True, True]
    assert r["lines"][3]["p_from_mw"] == 0.0
    # The load that line fed still gets served through the rest of the mesh.
    assert all(b["vm_pu"] is not None for b in r["buses"])


def test_island_is_de_energized():
    net = wscc9_3sm()
    # Bus 5 hangs on two lines (2 -> 5 and 3 -> 5); opening both islands it with its load.
    ends = [i for i, ln in enumerate(net.lines) if 5 in (ln.from_bus, ln.to_bus)]
    for i in ends:
        net.lines[i].from_closed = False
    st = service_state(net)
    assert 5 not in st.energized_buses
    assert not any(st.loads[i] for i, ld in enumerate(net.loads) if ld.bus == 5)
    r = client.post("/api/network/powerflow", json={"network": _dump(net)}).json()
    assert r["converged"]
    bus5 = next(b for b in r["buses"] if b["bus"] == 5)
    assert bus5["vm_pu"] is None
    v = client.post("/api/network/validate", json={"network": _dump(net)}).json()
    assert v["ok"] and 5 not in v["service"]["energized_buses"]
    assert any("de-energized" in i["message"] for i in v["issues"])


def test_opening_the_slack_moves_the_reference():
    """The slack is a modelling device, not a physical requirement: with it
    out of service another grid former takes the role (see breakers)."""
    net = wscc9_3sm()
    slack = next(d for d in net.der_units if d.bus_type.value == "slack")
    slack.closed = False
    st = service_state(net)
    assert st.references and slack.id not in st.references
    assert st.energized_buses == {b.id for b in net.buses}   # the rest stays energized
    pf = client.post("/api/network/powerflow", json={"network": _dump(net)}).json()
    assert pf["converged"]
    assert [g["name"] for g in pf["external_grid"]] == [f"der{st.references[0]}"]
    v = client.post("/api/network/validate", json={"network": _dump(net)}).json()
    assert v["ok"] and any("designated slack" in i["message"] for i in v["issues"])
    assert client.post("/api/network/modal", json={"network": _dump(net)}).status_code == 200
    # The reduced network the dynamic models use carries the role over.
    assert [d.id for d in energized_network(net).der_units if d.bus_type.value == "slack"] == [st.references[0]]


def test_each_island_is_solved_against_its_own_reference():
    net = wscc9_3sm()
    for ln in net.lines:                                    # cut bus 1 (and its machine) loose
        if 1 in (ln.from_bus, ln.to_bus):
            ln.from_closed = False
    st = service_state(net)
    assert len(st.references) == 2
    assert sorted(len(i.buses) for i in st.islands) == [2, 7]
    pf = client.post("/api/network/powerflow", json={"network": _dump(net)}).json()
    assert pf["converged"] and len(pf["external_grid"]) == 2
    assert all(b["vm_pu"] is not None for b in pf["buses"])  # both islands are alive
    v = client.post("/api/network/validate", json={"network": _dump(net)}).json()
    assert v["ok"] and any("2 energized islands" in i["message"] for i in v["issues"])


def test_an_island_without_a_grid_former_is_blacked_out():
    """A grid-following converter can only follow a voltage, never start one."""
    net = wscc9_1sm_1gfm_1gfl()
    gfl = next(d for d in net.der_units if d.unit_type.value == "gfl")
    hv = next(t.hv_bus for t in net.transformers if t.lv_bus == gfl.bus)
    for ln in net.lines:
        if hv in (ln.from_bus, ln.to_bus):
            ln.from_closed = False
    st = service_state(net)
    dead = next(i for i in st.islands if i.reference is None)
    assert dead.buses == {hv, gfl.bus}
    assert not st.der_units[gfl.id]
    pf = client.post("/api/network/powerflow", json={"network": _dump(net)}).json()
    assert pf["converged"]
    assert next(b for b in pf["buses"] if b["bus"] == hv)["vm_pu"] is None


def test_a_network_with_no_grid_former_left_is_an_error():
    net = wscc9_1sm_1gfm_1gfl()
    for d in net.der_units:
        if d.unit_type.value in ("sm", "gfm", "infinite_bus"):
            d.closed = False
    v = client.post("/api/network/validate", json={"network": _dump(net)}).json()
    assert not v["ok"] and any("set a voltage and a frequency" in i["message"] for i in v["issues"])
    assert client.post("/api/network/powerflow", json={"network": _dump(net)}).status_code == 422
    assert client.post("/api/network/modal", json={"network": _dump(net)}).status_code == 422


def test_open_load_modal_equals_removed_load():
    """Opening a load's breaker is the same, for the dynamic model, as not
    having the load at all -- but the names keep the full network's numbers."""
    net = wscc9_3sm()
    net.loads[0].closed = False
    e = energized_network(net)
    sys_open = linearize_network(e, run_power_flow(e))

    removed = wscc9_3sm()
    del removed.loads[0]
    sys_removed = linearize_network(removed, run_power_flow(removed))
    ev_open = np.sort_complex(analyze(sys_open.A, sys_open.state_names).eigenvalues)
    ev_removed = np.sort_complex(analyze(sys_removed.A, sys_removed.state_names).eigenvalues)
    np.testing.assert_allclose(ev_open, ev_removed, rtol=1e-9, atol=1e-9)
    assert not any("{Ld_1}" in n for n in sys_open.state_names)
    assert any("{Ld_2}" in n for n in sys_open.state_names) and any("{Ld_3}" in n for n in sys_open.state_names)


def test_open_unit_drops_its_states():
    net = wscc9_3sm()
    net.der_units[2].closed = False
    r = client.post("/api/network/modal", json={"network": _dump(net)})
    assert r.status_code == 200
    s = client.post("/api/network/states", json={"network": _dump(net)}).json()
    assert not any("SM_3" in n for n in s["state_names"])
    assert "dw_r_{SM_2}" in s["state_names"]
    assert not any("unit3" in m["name"] for m in s["measurements"])


def _emt(net, event, **kw):
    body = dict(
        network=_dump(net), perturb_kind="event", event=event, t_final=0.3, t_pre=0.01,
        plot_states=["dw_r_{SM_2}", "dw_r_{SM_3}"], plot_measurements=["P_from_{line3}", "V_{bus5}"], **kw,
    )
    return client.post("/api/network/emt", json=body)


def test_event_line_trip():
    r = _emt(wscc9_3sm(), {"kind": "breaker", "element": "line", "index": 3})
    assert r.status_code == 200, r.text
    d = r.json()
    t = np.array(d["t"])
    p = np.array(d["measurements"]["P_from_{line3}"], dtype=float)
    assert abs(p[t < 0]).min() > 0.1        # carried power before the trip
    assert np.all(p[t > 0] == 0.0)          # none after
    assert "opened" in d["perturbed"]
    assert max(abs(v) for v in d["series"]["dw_r_{SM_3}"]) > 1e-5


def test_event_unit_trip_gaps_its_signals():
    d = _emt(wscc9_3sm(), {"kind": "breaker", "element": "unit", "index": 3}).json()
    t = np.array(d["t"])
    sm3 = d["series"]["dw_r_{SM_3}"]
    assert all(v is not None for v, tk in zip(sm3, t) if tk < 0)
    assert all(v is None for v, tk in zip(sm3, t) if tk > 0)
    # Losing a generator slows the others down.
    assert d["series"]["dw_r_{SM_2}"][-1] < 0


def test_event_load_step_and_no_linear_overlay():
    d = _emt(wscc9_3sm(), {"kind": "load_step", "index": 1, "dp_pct": 20, "dq_pct": 0}, linear_overlay=True).json()
    assert d["linear"] is None and "linear" in d["linear_note"]
    assert d["series"]["dw_r_{SM_2}"][-1] < 0  # more load -> slower


def test_event_phase_jump_has_linear_overlay():
    net = gfm_smib()
    ib = next(d for d in net.der_units if d.bus_type.value == "slack")
    body = dict(
        network=_dump(net), perturb_kind="event", event={"kind": "phase_jump", "bus": ib.bus, "angle_deg": 2},
        t_final=0.3, plot_states=[], linear_overlay=True,
    )
    r = client.post("/api/network/emt", json=body)
    assert r.status_code == 200, r.text
    d = r.json()
    assert d["linear"] is not None and "infinite-bus source" in d["perturbed"]


@pytest.mark.parametrize("event", [
    {"kind": "breaker", "element": "line", "index": 42},    # no such line
    {"kind": "load_step", "index": 0, "dq_pct": -100},      # zero reactance
    {"kind": "phase_jump", "bus": 7, "angle_deg": 5},       # a unit's own terminal bus
    {"kind": "surge"},
])
def test_bad_events_are_422(event):
    assert _emt(wscc9_3sm(), event).status_code == 422


def test_event_on_already_open_element_is_422():
    net = wscc9_3sm()
    net.lines[3].from_closed = False
    assert _emt(net, {"kind": "breaker", "element": "line", "index": 3}).status_code == 422


# --- The reference frame (g2elin_core.components.frame) ----------------------
def _eigs(net):
    e = energized_network(net)
    s = linearize_network(e, run_power_flow(e))
    return np.sort_complex(analyze(s.A, s.state_names).eigenvalues)


@pytest.mark.parametrize("build", [wscc9_3sm, gfm_smib])
def test_a_frame_of_its_own_reproduces_the_slack_tied_model(build):
    """The frame is only a choice of coordinates: standing on its own but
    following the same machine, it must give the same modes as the MATLAB
    convention it replaces (plus its own angle state)."""
    compat = build()
    compat.frame_follows_slack = True
    free = _eigs(build())
    fixed = _eigs(compat)
    assert len(free) == len(fixed) + 1                       # the frame angle
    for z in fixed:
        assert min(abs(free - z)) / (1 + abs(z)) < 1e-5


def test_the_slack_unit_can_be_tripped_in_an_emt_run():
    net = wscc9_3sm()
    slack = next(d for d in net.der_units if d.bus_type.value == "slack")
    body = dict(
        network=_dump(net), perturb_kind="event", t_final=0.3, t_pre=0.01,
        event={"kind": "breaker", "element": "unit", "index": slack.id},
        plot_states=["dw_r_{SM_1}", "dw_r_{SM_2}"],
    )
    r = client.post("/api/network/emt", json=body)
    assert r.status_code == 200, r.text
    d = r.json()
    t = np.array(d["t"])
    gone = d["series"]["dw_r_{SM_1}"]
    assert all(v is None for v, tk in zip(gone, t) if tk > 0)   # the slack's own states stop
    assert d["series"]["dw_r_{SM_2}"][-1] < 0                   # the machines left pick up its load
    assert "tripped" in d["perturbed"]


def test_tripping_the_slack_is_refused_when_the_frame_follows_it():
    net = wscc9_3sm()
    net.frame_follows_slack = True
    slack = next(d for d in net.der_units if d.bus_type.value == "slack")
    r = _emt(net, {"kind": "breaker", "element": "unit", "index": slack.id})
    assert r.status_code == 422 and "frame_follows_slack" in r.text


def test_an_islanded_group_keeps_running_at_its_own_frequency():
    """Opening the last line to a machine leaves it alone with no load: it
    speeds up, while the rest of the network, a generator short, slows down."""
    net = wscc9_3sm()
    first, second = [i for i, ln in enumerate(net.lines) if 1 in (ln.from_bus, ln.to_bus)][:2]
    net.lines[first].from_closed = False                     # already open before the run
    body = dict(
        network=_dump(net), perturb_kind="event", t_final=0.4, t_pre=0.01,
        event={"kind": "breaker", "element": "line", "index": second},
        plot_states=["dw_r_{SM_1}", "dw_r_{SM_2}"], plot_measurements=[],
    )
    r = client.post("/api/network/emt", json=body)
    assert r.status_code == 200, r.text
    d = r.json()
    islanded, rest = d["series"]["dw_r_{SM_1}"][-1], d["series"]["dw_r_{SM_2}"][-1]
    assert islanded > 1e-4 > 0 > rest                        # they drift apart, each on its own
