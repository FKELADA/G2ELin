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
from g2elin_core.network.presets import gfm_smib, wscc9_3sm
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


def test_slack_breaker_cannot_open():
    net = wscc9_3sm()
    next(d for d in net.der_units if d.bus_type.value == "slack").closed = False
    v = client.post("/api/network/validate", json={"network": _dump(net)}).json()
    assert not v["ok"]
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
    {"kind": "breaker", "element": "unit", "index": 1},     # the slack
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
