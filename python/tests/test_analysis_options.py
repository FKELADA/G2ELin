"""Power-flow solver options, batch power flow, multi-output modal responses,
and the EMT pre-disturbance samples / linearised overlay -- the API additions
behind the web UI's Power Flow, Modal and EMT pages.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core.network.presets import wscc9_3sm

client = TestClient(app)


def _net() -> dict:
    return wscc9_3sm().model_dump()


def test_powerflow_algorithms_listed():
    algos = client.get("/api/powerflow/algorithms").json()
    assert {"nr", "iwamoto_nr", "fdbx", "fdxb", "gs", "bfsw"} <= set(algos)


@pytest.mark.parametrize("algorithm", ["nr", "iwamoto_nr", "fdbx", "fdxb"])
def test_powerflow_with_options_matches_default(algorithm):
    base = client.post("/api/network/powerflow", json={"network": _net()}).json()
    r = client.post("/api/network/powerflow", json={
        "network": _net(), "options": {"algorithm": algorithm, "max_iteration": 50, "tolerance_mva": 1e-8, "init": "flat"},
    })
    assert r.status_code == 200
    body = r.json()
    assert body["converged"] and body["algorithm"] == algorithm
    assert body["iterations"] is not None
    for a, b in zip(base["buses"], body["buses"]):
        assert a["vm_pu"] == pytest.approx(b["vm_pu"], abs=1e-6)


def test_preset_powerflow_still_accepts_a_bare_post():
    r = client.post("/api/presets/wscc9_3sm/powerflow")
    assert r.status_code == 200 and r.json()["converged"]


def test_powerflow_iteration_cap_reports_non_convergence():
    r = client.post("/api/network/powerflow", json={"network": _net(), "options": {"algorithm": "nr", "max_iteration": 1, "init": "flat"}})
    assert r.status_code == 200
    assert r.json()["converged"] is False


@pytest.mark.parametrize("options", [{"algorithm": "nope"}, {"init": "results"}, {"tolerance_mva": 0}, {"max_iteration": 0}])
def test_powerflow_bad_options_422(options):
    r = client.post("/api/network/powerflow", json={"network": _net(), "options": options})
    assert r.status_code == 422


def test_powerflow_solver_exception_is_422_not_500():
    # Backward/forward sweep needs a radial network; WSCC-9 is meshed.
    r = client.post("/api/network/powerflow", json={"network": _net(), "options": {"algorithm": "bfsw"}})
    assert r.status_code == 422
    assert "solver failed" in r.json()["detail"]


def test_batch_powerflow_ramps_from_base_to_targets():
    r = client.post("/api/network/powerflow/batch", json={
        "network": _net(), "load_p_scale": 1.5, "load_q_scale": 0.5, "der_scale": {"sm": 1.2}, "steps": 5,
    })
    assert r.status_code == 200
    snaps = r.json()["snapshots"]
    assert len(snaps) == 6
    assert snaps[0]["load_p_scale"] == pytest.approx(1.0) and snaps[-1]["load_p_scale"] == pytest.approx(1.5)
    assert snaps[-1]["load_q_scale"] == pytest.approx(0.5)
    assert snaps[-1]["der_scale"]["sm"] == pytest.approx(1.2)
    assert all(s["converged"] for s in snaps)
    # Network loads only -- units' own auxiliary loads (p_cons) are not load-scaled.
    load_p = [sum(ld["p_mw"] for ld in s["loads"] if not ld["name"].endswith("_aux_load")) for s in snaps]
    assert load_p[-1] == pytest.approx(1.5 * load_p[0], rel=1e-6)
    losses = [s["total_losses_mw"] for s in snaps]
    assert losses[-1] > losses[0]  # heavier loading -> more losses


def test_preset_batch_powerflow():
    r = client.post("/api/presets/wscc9_3sm/powerflow/batch", json={"load_p_scale": 0.8, "steps": 2})
    assert r.status_code == 200 and len(r.json()["snapshots"]) == 3


@pytest.mark.parametrize("body", [{"steps": 0}, {"load_p_scale": -0.1}, {"der_scale": {"windfarm": 1.0}}])
def test_batch_powerflow_bad_input_422(body):
    r = client.post("/api/network/powerflow/batch", json={"network": _net(), **body})
    assert r.status_code == 422


def test_step_response_several_outputs_match_single_output_calls():
    base = {"network": _net(), "input_name": "P_ref_{SM_2}", "amplitude": 0.1, "t_final": 1.0}
    multi = client.post("/api/network/modal/step_response", json={**base, "output_names": ["p_e_{SM_2}", "w_r_{SM_2}"]}).json()
    single = client.post("/api/network/modal/step_response", json={**base, "output_name": "w_r_{SM_2}"}).json()
    assert list(multi["series"]) == ["p_e_{SM_2}", "w_r_{SM_2}"]
    assert multi["y"] == multi["series"]["p_e_{SM_2}"]
    assert np.allclose(multi["series"]["w_r_{SM_2}"], single["y"])


def test_step_response_without_outputs_422():
    r = client.post("/api/network/modal/step_response", json={"network": _net(), "input_name": "P_ref_{SM_2}"})
    assert r.status_code == 422


def test_free_response_explicit_plot_states():
    r = client.post("/api/network/modal/free_response", json={
        "network": _net(), "perturb_state": "dw_r_{SM_2}", "plot_states": ["dw_r_{SM_1}", "dw_r_{SM_3}"],
    })
    assert r.status_code == 200
    assert list(r.json()["series"]) == ["dw_r_{SM_1}", "dw_r_{SM_3}"]
    bad = client.post("/api/network/modal/free_response", json={"network": _net(), "perturb_state": "dw_r_{SM_2}", "plot_states": ["nope"]})
    assert bad.status_code == 422


def _emt(**kw) -> dict:
    body = {"network": _net(), "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.002,
            "t_final": 0.3, "plot_states": ["dw_r_{SM_1}", "dw_r_{SM_2}"], "plot_outputs": ["p_e_{SM_2}"], **kw}
    r = client.post("/api/network/emt", json=body)
    assert r.status_code == 200, r.text
    return r.json()


def test_emt_t_pre_prepends_equilibrium_samples():
    body = _emt(t_pre=0.01)
    assert body["t"][:3] == [pytest.approx(-0.01), 0.0, 0.0]
    s = body["series"]["dw_r_{SM_2}"]
    assert s[0] == s[1] == pytest.approx(0.0, abs=1e-12)  # undisturbed before T0
    assert s[2] == pytest.approx(0.002, abs=1e-9)  # the offset applied at T0
    assert len(body["outputs"]["p_e_{SM_2}"]) == len(body["t"])


def test_emt_t_pre_out_of_range_422():
    r = client.post("/api/network/emt", json={"network": _net(), "perturb_name": "dw_r_{SM_2}", "t_final": 0.3, "t_pre": 5.0})
    assert r.status_code == 422


def test_emt_default_has_no_linear_overlay_and_no_pre_samples():
    body = _emt()
    assert body["linear"] is None
    assert body["t"][0] == 0.0


@pytest.mark.parametrize("kind,name,offset", [("state", "dw_r_{SM_2}", 0.002), ("input", "P_ref_{SM_2}", 0.01)])
def test_emt_linear_overlay_matches_nonlinear_deviation(kind, name, offset):
    # The nonlinear model starts slightly off-equilibrium (see
    # test_emt_simulation.test_network_operating_point_is_not_a_perfect_equilibrium),
    # so compare the *effect of the disturbance* -- EMT(disturbed) - EMT(undisturbed)
    # against the linear overlay's deviation -- for a small disturbance.
    kw = dict(perturb_kind=kind, perturb_name=name, linear_overlay=True, t_pre=0.01)
    a, z = _emt(perturb_offset=offset, **kw), _emt(perturb_offset=0.0, **kw)
    assert a["linear"]["t"] == a["t"]
    for grp in ("series", "outputs"):
        for n in a[grp]:
            d_nl = np.array(a[grp][n]) - np.array(z[grp][n])
            d_lin = np.array(a["linear"][grp][n]) - np.array(z["linear"][grp][n])
            assert np.max(np.abs(d_nl - d_lin)) < 0.1 * np.ptp(d_nl), n


def test_emt_live_streams_pre_samples_and_linear_overlay():
    r = client.post("/api/network/emt/live", json={
        "network": _net(), "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.002,
        "t_final": 0.05, "plot_states": ["dw_r_{SM_2}"], "t_pre": 0.01, "linear_overlay": True,
    })
    lines = [json.loads(x) for x in r.text.strip().split("\n")]
    assert lines[0]["t"] == pytest.approx(-0.01) and lines[1]["t"] == 0.0
    done = lines[-1]
    assert done["done"] and "linear" in done
    assert done["linear"]["t"][0] == pytest.approx(-0.01)
    assert done["n_steps"] == len(lines) - 3  # two pre-samples + the done line aren't solver steps
