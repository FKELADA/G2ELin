"""Loop tuning (g2elin_core.tuning): response time/damping <-> gains, checked
against the default tunings operating_point.py was derived with, and the
root-locus sweep of a loop's response time. Also the SMSM and CIGRE
interconnected presets."""

from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core import tuning
from g2elin_core.modal import analyze, reference_angle_modes
from g2elin_core.network.presets import (
    cigre_interconnected_1sm_1gfm_1gfl, gfl_smsm, gfm_smsm, sm_smsm, wscc9_1sm_1gfm_1gfl,
)
from g2elin_core.network.validation import validate_network
from g2elin_core.operating_point import gfl_params, gfm_params
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow

GFM = gfm_params(sn_mva=100, f_hz=60, un_kv=18, rt_pu=0.0, lt_pu=0.0576)
GFL = gfl_params(sn_mva=100, f_hz=60, un_kv=18, rt_pu=0.0, lt_pu=0.0576)

# The response times/damping operating_point.py derives the default gains from.
DEFAULT_TUNINGS = [
    ("gfm", GFM, "cl", 0.1e-3, 0.707), ("gfm", GFM, "vl", 15e-3, 0.707), ("gfm", GFM, "dc", 5e-3, None),
    ("gfl", GFL, "cl", 10e-3, 0.707), ("gfl", GFL, "pll", 50e-3, 0.707), ("gfl", GFL, "dcv", 100e-3, 0.707),
    ("gfl", GFL, "q", 100e-3, None),
]


@pytest.mark.parametrize("unit,params,loop_id,tr,zeta", DEFAULT_TUNINGS)
def test_default_gains_round_trip(unit, params, loop_id, tr, zeta):
    lp = tuning.loop(unit, loop_id)
    t = tuning.tuning_of(lp, params)
    assert t["tr"] == pytest.approx(tr, rel=1e-9)
    if zeta is not None:
        assert t["zeta"] == pytest.approx(zeta, rel=1e-9)
    # ... and back: the same tuning reproduces the default gains.
    for k, v in tuning.gains_for(lp, params, tr=tr, **({"zeta": zeta} if zeta else {})).items():
        assert v == pytest.approx(params[k], rel=1e-9)


def test_gains_for_keeps_the_other_quantity():
    lp = tuning.loop("gfm", "vl")
    new = {**GFM, **tuning.gains_for(lp, GFM, tr=30e-3)}
    t = tuning.tuning_of(lp, new)
    assert t["tr"] == pytest.approx(30e-3) and t["zeta"] == pytest.approx(0.707)


def test_droop_loop():
    lp = tuning.loop("gfm", "droop")
    assert tuning.tuning_of(lp, GFM)["H"] == pytest.approx(3.0)  # _H_FIRST_ORDER
    assert tuning.gains_for(lp, GFM, Tf=0.05) == {"wf": pytest.approx(20.0)}
    assert tuning.gains_for(lp, GFM, H=6.0)["wf"] == pytest.approx(GFM["wf"] / 2)


def test_unknown_loop():
    with pytest.raises(KeyError):
        tuning.loop("gfm", "pll")


client = TestClient(app)


def _sweep(net: dict, field: str, start, stop, step, unit="gfm"):
    key = next(d["id"] for d in net["der_units"] if d["unit_type"] == unit)
    r = client.post("/api/network/modal/sweep", json={"network": net, "target": {"element": "unit", "key": key, "field": field},
                                                      "start": start, "stop": stop, "step": step})
    return r, ([json.loads(x) for x in r.text.strip().split("\n")] if r.status_code == 200 else [])


def test_sweep_loop_response_time():
    net = wscc9_1sm_1gfm_1gfl().model_dump()
    r, lines = _sweep(net, "tune.vl.tr_ms", 10, 30, 10)
    steps = [x for x in lines if "i" in x]
    assert r.status_code == 200 and len(steps) == 3 and all(s["ok"] for s in steps)
    assert steps[0]["eig"] != steps[-1]["eig"]
    # t_r = 15 ms is the default tuning: that step must match the unmodified network.
    r2, lines2 = _sweep(net, "tune.vl.tr_ms", 15, 15, 1)
    base = client.post("/api/network/modal", json={"network": net}).json()
    got = sorted(round(re, 6) for re, _ in [x for x in lines2 if "i" in x][0]["eig"])
    want = sorted(round(m["real"], 6) for m in base["modes"])
    assert got == want


@pytest.mark.parametrize("field,msg", [("tune.pll.tr_ms", "no control loop"), ("tune.dc.zeta", "can be swept by"), ("tune.vl", "tune.<loop>")])
def test_sweep_bad_tuning_target_422(field, msg):
    r, _ = _sweep(wscc9_1sm_1gfm_1gfl().model_dump(), field, 1, 2, 1)
    assert r.status_code == 422 and msg in r.json()["detail"]


@pytest.mark.parametrize("build", [sm_smsm, gfm_smsm, gfl_smsm, cigre_interconnected_1sm_1gfm_1gfl])
def test_new_presets_are_valid_and_stable(build):
    net = build()
    assert not [i for i in validate_network(net) if i.severity == "error"]
    result = run_power_flow(net)
    assert result.converged
    system = linearize_network(net, result)
    modal = analyze(system.A, system.state_names)
    # Leaving out the reference-angle modes, which the formulation puts on
    # the imaginary axis (modal.reference_angle_modes).
    ref = set(reference_angle_modes(modal))
    assert max(z.real for j, z in enumerate(modal.eigenvalues) if j not in ref) < 1e-6


def test_cigre_interconnected_topology():
    net = cigre_interconnected_1sm_1gfm_1gfl()
    assert 14 not in {b.id for b in net.buses}
    assert any((ln.from_bus, ln.to_bus) == (13, 1) for ln in net.lines)  # feeder 2 starts at node 1
    load1 = next(ld for ld in net.loads if ld.bus == 1)
    assert load1.p_mw == pytest.approx(2.0) and load1.q_mvar == pytest.approx(0.555)  # loads 1 + 14 merged
    slack = next(d for d in net.der_units if d.bus_type.value == "slack")
    assert slack.unit_type.value == "infinite_bus"


def test_presets_exposed_by_api():
    ids = {p["id"] for p in client.get("/api/presets").json()}
    assert {"sm_smsm", "gfm_smsm", "gfl_smsm", "cigre_interconnected_1sm_1gfm_1gfl"} <= ids
