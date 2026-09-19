"""Root-locus parameter sweeps (/api/network/modal/sweep) and the unit
parameter defaults endpoint the editor uses (/api/units/defaults)."""

from __future__ import annotations

import json

import numpy as np
import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_api.schemas import SweepRequest, SweepTarget
from g2elin_api.sweep import match_to, sweep_values
from g2elin_core.network.presets import wscc9_1sm_1gfm_1gfl
from g2elin_core.operating_point import gfm_params

client = TestClient(app)


def _net() -> dict:
    return wscc9_1sm_1gfm_1gfl().model_dump()


def _gfm_id(net: dict) -> int:
    return next(d["id"] for d in net["der_units"] if d["unit_type"] == "gfm")


def _sweep(body: dict):
    r = client.post("/api/network/modal/sweep", json=body)
    return r, ([json.loads(x) for x in r.text.strip().split("\n")] if r.status_code == 200 else [])


def _req(start, stop, step):
    return SweepRequest(target=SweepTarget(element="network", field="f_hz"), start=start, stop=stop, step=step)


def test_sweep_values_include_both_ends():
    assert sweep_values(_req(1, 2, 0.25)) == pytest.approx([1, 1.25, 1.5, 1.75, 2])
    assert sweep_values(_req(2, 1, 0.5)) == pytest.approx([2, 1.5, 1])  # direction from the range
    assert sweep_values(_req(0, 1, 0.3)) == pytest.approx([0, 0.3, 0.6, 0.9, 1])  # end kept even off-grid


def test_match_to_follows_each_mode():
    prev = np.array([-1 + 10j, -1 - 10j, -100 + 0j])
    cur = np.array([-99 + 0j, -1.2 - 10.5j, -1.2 + 10.5j])
    assert np.allclose(match_to(prev, cur), [-1.2 + 10.5j, -1.2 - 10.5j, -99])


def test_sweep_gfm_gain_streams_one_line_per_value():
    net = _net()
    r, lines = _sweep({"network": net, "target": {"element": "unit", "key": _gfm_id(net), "field": "params.KpVL"},
                       "start": 0.05, "stop": 0.15, "step": 0.05})
    assert r.status_code == 200 and r.headers["content-type"].startswith("application/x-ndjson")
    head, *steps, done = lines
    assert head["values"] == pytest.approx([0.05, 0.1, 0.15])
    assert done == {"done": True, "n": 3}
    assert [s["i"] for s in steps] == [0, 1, 2] and all(s["ok"] for s in steps)
    n = {len(s["eig"]) for s in steps}
    assert len(n) == 1  # same model size at every step
    # The gain actually moves the eigenvalues.
    assert steps[0]["eig"] != steps[-1]["eig"]


def test_sweep_operating_point_parameter():
    net = _net()
    r, lines = _sweep({"network": net, "target": {"element": "load", "key": 0, "field": "p_mw"},
                       "start": 80, "stop": 120, "step": 40})
    steps = [x for x in lines if "i" in x]
    assert len(steps) == 2 and all(s["ok"] for s in steps)


def test_sweep_bad_value_is_reported_per_step_not_fatal():
    net = _net()
    # x_pu must be > 0: the 0 step fails validation, the others still run.
    r, lines = _sweep({"network": net, "target": {"element": "line", "key": 0, "field": "x_pu"},
                       "start": 0.0, "stop": 0.1, "step": 0.05})
    steps = [x for x in lines if "i" in x]
    assert steps[0]["ok"] is False and "error" in steps[0]
    assert all(s["ok"] for s in steps[1:])
    assert lines[-1]["done"]


@pytest.mark.parametrize("target,msg", [
    ({"element": "unit", "key": 999, "field": "p_set_mw"}, "no unit"),
    ({"element": "line", "key": 0, "field": "from_bus"}, "can't be swept"),
    ({"element": "unit", "key": 1, "field": "params.Nope"}, "has no parameter"),
    ({"element": "planet", "field": "mass"}, "unknown element"),
])
def test_sweep_bad_target_422(target, msg):
    r, _ = _sweep({"network": _net(), "target": target, "start": 0, "stop": 1, "step": 0.5})
    assert r.status_code == 422 and msg in r.json()["detail"]


def test_sweep_too_many_points_422():
    r, _ = _sweep({"network": _net(), "target": {"element": "network", "field": "f_hz"}, "start": 50, "stop": 60, "step": 0.001})
    assert r.status_code == 422


def test_unit_defaults_endpoint_matches_core():
    body = {"unit_type": "gfm", "sn_mva": 100, "f_hz": 60, "un_kv": 18, "rt_pu": 0.0, "lt_pu": 0.0576}
    r = client.post("/api/units/defaults", json=body)
    assert r.status_code == 200
    assert r.json()["params"] == pytest.approx(gfm_params(sn_mva=100, f_hz=60, un_kv=18, rt_pu=0.0, lt_pu=0.0576))
    assert client.post("/api/units/defaults", json={**body, "unit_type": "infinite_bus"}).json()["params"] == {}
    assert client.post("/api/units/defaults", json={**body, "unit_type": "windmill"}).status_code == 422


# --- Several parameters at once --------------------------------------------------

def test_extra_values_move_in_lockstep():
    from g2elin_api.schemas import SweepExtra
    from g2elin_api.sweep import extra_values

    req = SweepRequest(target=SweepTarget(element="network", field="f_hz"), start=50, stop=60, step=5,
                       extra=[SweepExtra(target=SweepTarget(element="load", key=0, field="p_mw"), start=1, stop=3)])
    values = sweep_values(req)
    assert extra_values(req, values) == [pytest.approx([1, 2, 3])]


def test_combined_sweep_equals_both_changes_applied_by_hand():
    net = _net()
    gfm = _gfm_id(net)
    r, lines = _sweep({
        "network": net, "target": {"element": "unit", "key": gfm, "field": "params.KpVL"},
        "start": 0.05, "stop": 0.15, "step": 0.1,
        "extra": [{"target": {"element": "load", "key": 0, "field": "p_mw"}, "start": 80, "stop": 120}],
    })
    head, *steps, done = lines
    assert head["extra_values"] == [pytest.approx([80, 120])]
    last = steps[-1]
    assert last["ok"]
    # The last step is the network with *both* parameters at their end values.
    by_hand = _net()
    next(d for d in by_hand["der_units"] if d["id"] == gfm)["params"] = {"KpVL": 0.15}
    by_hand["loads"][0]["p_mw"] = 120
    base = client.post("/api/network/modal", json={"network": by_hand}).json()
    assert sorted(round(re, 5) for re, _ in last["eig"]) == sorted(round(m["real"], 5) for m in base["modes"])


def test_same_parameter_twice_422():
    net = _net()
    t = {"element": "load", "key": 0, "field": "p_mw"}
    r, _ = _sweep({"network": net, "target": t, "start": 80, "stop": 90, "step": 10,
                   "extra": [{"target": t, "start": 1, "stop": 2}]})
    assert r.status_code == 422 and "more than once" in r.json()["detail"]


def test_bad_extra_target_422():
    r, _ = _sweep({"network": _net(), "target": {"element": "load", "key": 0, "field": "p_mw"}, "start": 80, "stop": 90,
                   "step": 10, "extra": [{"target": {"element": "line", "key": 0, "field": "from_bus"}, "start": 1, "stop": 2}]})
    assert r.status_code == 422


def test_length_sweep_scales_line_impedances():
    from g2elin_api.sweep import network_with
    from g2elin_core.network.schema import Network

    net = Network(**_net())
    line = net.lines[0]
    req = SweepRequest(target=SweepTarget(element="line", key=0, field="length_km"), start=1, stop=2, step=1)
    out = network_with(net, req, 2 * line.length_km).lines[0]
    assert out.length_km == pytest.approx(2 * line.length_km)
    for k in ("r_pu", "x_pu", "b_pu"):
        assert getattr(out, k) == pytest.approx(2 * getattr(line, k))
