"""EMT measurement outputs (g2elin_core.timedomain.measurements): at the
operating point they must reproduce the power flow; over a simulation they
must be consistent (aligned with t, 3-phase balanced, nominal frequency
before the disturbance) and respond to the disturbance."""

from __future__ import annotations

import json

import numpy as np
import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core.network.presets import wscc9_1sm_1gfm_1gfl, wscc9_3sm
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network
from g2elin_core.timedomain.measurements import MeasurementSet

client = TestClient(app)


@pytest.fixture(scope="module")
def at_operating_point():
    net = wscc9_1sm_1gfm_1gfl()
    result = run_power_flow(net)
    model = build_nonlinear_network(net, result)
    ms = MeasurementSet(model)
    x0 = model.initial_state()
    z0, u0 = model.solve_algebraic(x0, model.default_u_exo())
    return net, result, ms, ms.evaluator(ms.names())(x0, z0, u0)


def test_catalog_covers_every_element(at_operating_point):
    net, _, ms, _ = at_operating_point
    names = set(ms.names())
    for b in (1, 4):
        assert {f"V_{{bus{b}}}", f"angle_{{bus{b}}}", f"f_{{bus{b}}}", f"f_inst_{{bus{b}}}", f"v_a_{{bus{b}}}"} <= names
    assert {"P_from_{line0}", "Q_to_{line5}", "P_{load2}", "P_grid_{unit2}", "f_{unit3}", "P_hv_{trafo1}"} <= names
    assert not any("bus7" in n for n in names)  # a unit's own terminal bus isn't a dynamic node


def test_bus_voltages_and_angles_match_power_flow(at_operating_point):
    net, result, _, m = at_operating_point
    bt = result.bus_table().set_index("bus")
    for bus in range(1, 7):
        assert m[f"V_{{bus{bus}}}"] == pytest.approx(bt.loc[bus, "vm_pu"], abs=1e-6)
        assert m[f"angle_{{bus{bus}}}"] == pytest.approx(bt.loc[bus, "va_degree"], abs=1e-4)


def test_branch_and_load_powers_match_power_flow(at_operating_point):
    net, result, _, m = at_operating_point
    sb, lt, ld, tt = net.sn_mva, result.line_table(), result.load_table(), result.trafo_table()
    for i in range(len(net.lines)):
        assert m[f"P_from_{{line{i}}}"] * sb == pytest.approx(lt.p_from_mw[i], abs=1e-3)
        assert m[f"Q_from_{{line{i}}}"] * sb == pytest.approx(lt.q_from_mvar[i], abs=1e-3)
        assert m[f"P_to_{{line{i}}}"] * sb == pytest.approx(-lt.p_to_mw[i], abs=1e-3)
        assert m[f"Q_to_{{line{i}}}"] * sb == pytest.approx(-lt.q_to_mvar[i], abs=1e-3)
    for i in range(len(net.loads)):
        assert m[f"P_{{load{i}}}"] * sb == pytest.approx(ld.p_mw[i], abs=1e-3)
    # The model's start point isn't an exact equilibrium (known), so the
    # slack's injection is only close; the dispatched units are exact.
    for j in range(len(net.transformers)):
        assert m[f"P_hv_{{trafo{j}}}"] * sb == pytest.approx(-tt.p_hv_mw[j], abs=0.05)


def test_frequencies_and_phases(at_operating_point):
    net, _, _, m = at_operating_point
    for u in (1, 2, 3):
        assert m[f"f_{{unit{u}}}"] == pytest.approx(net.f_hz, abs=1e-6)
    abc = [m[f"v_{p}_{{bus4}}"] for p in "abc"]
    assert sum(abc) == pytest.approx(0.0, abs=1e-9)
    assert max(abs(v) for v in abc) <= m["V_{bus4}"] + 1e-9


def test_states_endpoint_lists_measurements():
    r = client.post("/api/network/states", json={"network": wscc9_3sm().model_dump()}).json()
    ms = {m["name"]: m for m in r["measurements"]}
    assert ms["V_{bus4}"]["unit"] == "pu" and ms["f_{bus4}"]["unit"] == "Hz" and ms["angle_{bus4}"]["unit"] == "deg"
    assert ms["P_from_{line0}"]["group"].startswith("Line #0")


def _emt(**kw):
    body = {"network": wscc9_3sm().model_dump(), "perturb_kind": "input", "perturb_name": "P_ref_{SM_2}",
            "perturb_offset": 0.05, "t_final": 0.3, "plot_states": ["dw_r_{SM_2}"], **kw}
    r = client.post("/api/network/emt", json=body)
    assert r.status_code == 200, r.text
    return r.json()


def test_emt_measurements():
    names = ["V_{bus4}", "f_{bus4}", "f_inst_{bus4}", "v_a_{bus4}", "v_b_{bus4}", "v_c_{bus4}", "P_grid_{unit2}", "f_{unit2}"]
    b = _emt(plot_measurements=names, t_pre=0.01, dt=0.0005)
    t, meas = np.array(b["t"]), {k: np.array(v) for k, v in b["measurements"].items()}
    assert set(meas) == set(names) and all(len(v) == len(t) for v in meas.values())
    pre = t < 0
    assert pre.sum() >= 19  # dense pre-T0 samples (0.01 s at 0.5 ms) for the waveforms
    assert np.allclose(meas["f_{bus4}"][pre], 60.0) and np.allclose(meas["f_inst_{bus4}"][pre], 60.0)
    assert np.allclose(meas["v_a_{bus4}"] + meas["v_b_{bus4}"] + meas["v_c_{bus4}"], 0.0, atol=1e-9)
    # A sampled 60 Hz waveform: phase a actually oscillates before T0.
    assert np.ptp(meas["v_a_{bus4}"][pre]) > 1.0
    # The filter keeps the measured frequency far calmer than the raw one.
    post = t > 0
    assert np.ptp(meas["f_{bus4}"][post]) < 0.5 * np.ptp(meas["f_inst_{bus4}"][post])
    # The P_ref step of +0.05 pu raises SM 2's injection.
    p = meas["P_grid_{unit2}"]
    assert p[-1] > p[np.argmax(t >= 0)] + 0.01


def test_unknown_measurement_422():
    r = client.post("/api/network/emt", json={"network": wscc9_3sm().model_dump(), "perturb_name": "dw_r_{SM_2}",
                                               "t_final": 0.1, "plot_measurements": ["V_{bus99}"]})
    assert r.status_code == 422 and "unknown measurement" in r.json()["detail"]


def test_live_stream_carries_measurements():
    r = client.post("/api/network/emt/live", json={
        "network": wscc9_3sm().model_dump(), "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}",
        "perturb_offset": 0.002, "t_final": 0.02, "plot_states": ["dw_r_{SM_2}"], "t_pre": 0.005,
        "plot_measurements": ["V_{bus4}", "f_{bus4}"],
    })
    lines = [json.loads(x) for x in r.text.strip().split("\n")]
    steps = [x for x in lines if "t" in x]
    assert all(set(s["measurements"]) == {"V_{bus4}", "f_{bus4}"} for s in steps)
    assert steps[0]["t"] == pytest.approx(-0.005) and steps[0]["measurements"]["f_{bus4}"] == pytest.approx(60.0)
