"""API smoke tests — the HTTP layer over g2elin_core (P7). Uses FastAPI's
TestClient (in-process, no real server) so these run as part of the normal
suite; requires the `api` extra (`pip install -e ".[dev,api]"`).
"""

from __future__ import annotations

import json

import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from g2elin_api.main import app  # noqa: E402

client = TestClient(app)


def test_list_presets():
    r = client.get("/api/presets")
    assert r.status_code == 200
    ids = {p["id"] for p in r.json()}
    assert ids == {
        "wscc9_3sm", "wscc9_2sm_1gfl", "wscc9_1sm_2gfl", "wscc9_1sm_1gfm_1gfl",
        "wscc9_1sm_2gfm", "wscc9_2sm_1gfm", "wscc9_1gfm_2gfl", "wscc9_3gfm",
        "cigre_islanded_1sm_1gfm_1gfl", "cigre_islanded", "cigre_islanded_1sm_3gfm_1gfl",
        "cigre_islanded_2sm_2gfm_2gfl", "sm_smib", "gfm_smib", "gfl_smib",
    }


def test_unknown_preset_404s():
    r = client.post("/api/presets/nope/powerflow")
    assert r.status_code == 404


@pytest.mark.parametrize("preset_id", ["wscc9_3sm", "cigre_islanded", "sm_smib", "gfm_smib", "gfl_smib"])
def test_powerflow_endpoint(preset_id):
    r = client.post(f"/api/presets/{preset_id}/powerflow")
    assert r.status_code == 200
    body = r.json()
    assert body["converged"] is True
    assert len(body["buses"]) > 0
    assert body["total_losses_mw"] > 0


@pytest.mark.parametrize("preset_id", ["wscc9_3sm", "cigre_islanded", "sm_smib", "gfm_smib", "gfl_smib"])
def test_modal_endpoint(preset_id):
    r = client.post(f"/api/presets/{preset_id}/modal")
    assert r.status_code == 200
    body = r.json()
    assert body["stable"] is True
    assert body["n_states"] == len(body["modes"])
    assert len(body["state_names"]) == body["n_states"]
    assert len(body["participation"]) == body["n_states"]
    assert len(body["participation"][0]) == body["n_states"]


def test_modal_sensitivity_endpoint():
    r = client.post("/api/presets/wscc9_3sm/modal/sensitivity", json={"mode": 40})
    assert r.status_code == 200
    body = r.json()
    assert len(body["top"]) <= 8
    assert len(body["matrix"]) == len(body["state_names"])


def test_modal_sensitivity_bad_mode_422s():
    r = client.post("/api/presets/wscc9_3sm/modal/sensitivity", json={"mode": 9999})
    assert r.status_code == 422


def test_modal_mode_shape_endpoint():
    r = client.post("/api/presets/wscc9_3sm/modal/mode_shape", json={"mode": 40})
    assert r.status_code == 200
    body = r.json()
    assert len(body["states"]) == len(body["angles_deg"]) == 5


def test_modal_free_response_endpoint():
    r = client.post(
        "/api/presets/wscc9_3sm/modal/free_response",
        json={"perturb_state": "dw_r_{SM_2}", "offset": 0.02, "t_final": 0.5, "state_filter": "dw_r"},
    )
    assert r.status_code == 200
    body = r.json()
    assert set(body["series"].keys()) == {"dw_r_{SM_1}", "dw_r_{SM_2}", "dw_r_{SM_3}"}
    assert body["series"]["dw_r_{SM_2}"][0] == pytest.approx(0.02, abs=1e-6)


def test_modal_free_response_unknown_state_422s():
    r = client.post("/api/presets/wscc9_3sm/modal/free_response", json={"perturb_state": "not_a_real_state"})
    assert r.status_code == 422


def test_modal_step_response_endpoint():
    r = client.post(
        "/api/presets/wscc9_3sm/modal/step_response",
        json={"input_name": "P_ref_{SM_1}", "output_name": "p_e_{SM_1}", "amplitude": 0.1, "t_final": 0.5},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["t"]) == len(body["y"]) == 200
    assert body["y"][0] == pytest.approx(0.0, abs=1e-6)  # causal system, zero response before the step


def test_modal_step_response_unknown_io_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/modal/step_response",
        json={"input_name": "not_an_input", "output_name": "p_e_{SM_1}"},
    )
    assert r.status_code == 422


def test_timeseries_endpoint():
    r = client.post("/api/presets/wscc9_3sm/timeseries")
    assert r.status_code == 200
    body = r.json()
    labels = {s["label"] for s in body["snapshots"]}
    assert labels == {"70% load", "100% load", "130% load"}
    assert all(s["converged"] for s in body["snapshots"])


def test_static_index_served():
    r = client.get("/")
    assert r.status_code == 200
    assert "text/html" in r.headers["content-type"]
    assert "G2ELin" in r.text


def test_docs_tab_markup_present():
    r = client.get("/")
    assert 'data-panel="docs"' in r.text
    assert "/manual/index.html" in r.text


def test_manual_mount_serves_built_docs():
    # Skips rather than fails if the docs haven't been built locally (see
    # docs/sphinx/api_and_web.md) -- building them is a separate, explicit
    # `python tools/build_docs.py` step, not part of installing the package.
    r = client.get("/manual/index.html")
    if r.status_code == 404:
        pytest.skip("docs/sphinx/_build/html not present -- run `python tools/build_docs.py` first")
    assert r.status_code == 200
    assert "G2ELin documentation" in r.text


def test_states_endpoint():
    r = client.get("/api/presets/wscc9_3sm/states")
    assert r.status_code == 200
    body = r.json()
    names = body["state_names"]
    assert len(names) == 87
    assert "dw_r_{SM_2}" in names
    assert "theta_{SM_2}" in names
    assert "P_ref_{SM_2}" in body["input_names"]
    assert "p_e_{SM_2}" in body["output_names"]


def test_states_unknown_preset_404s():
    r = client.get("/api/presets/nope/states")
    assert r.status_code == 404


def test_emt_endpoint():
    # Default plot_states=[] -- server-side default is every dw_r_* state,
    # exercised here deliberately rather than passing an explicit list.
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.02, "t_final": 0.3},
    )
    assert r.status_code == 200
    body = r.json()
    assert body["perturbed"] == "dw_r_{SM_2}"
    assert body["perturb_kind"] == "state"
    assert set(body["series"].keys()) == {"dw_r_{SM_1}", "dw_r_{SM_2}", "dw_r_{SM_3}"}
    assert len(body["t"]) == 200
    assert all(len(vals) == 200 for vals in body["series"].values())
    # The perturbed state's own first sample should reflect the offset applied at t=0.
    assert body["series"]["dw_r_{SM_2}"][0] == pytest.approx(0.02, abs=1e-6)
    assert body["inputs"] == {}
    assert body["outputs"] == {}
    # default dt: t_final / (200 - 1)
    assert body["dt"] == pytest.approx(0.3 / 199, rel=1e-6)


def test_emt_input_perturbation():
    # A step in an exogenous reference held for the whole run, not an
    # initial-condition offset -- the system starts at its own equilibrium
    # and the equilibrium itself moves.
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={
            "perturb_kind": "input", "perturb_name": "P_ref_{SM_2}", "perturb_offset": 0.05, "t_final": 0.3,
            "plot_outputs": ["p_e_{SM_2}"],
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert body["perturbed"] == "P_ref_{SM_2}"
    assert body["perturb_kind"] == "input"
    # The state starts at equilibrium (no initial-condition jump) but the
    # commanded output should have visibly moved by the end of the run,
    # tracking the higher P_ref.
    p_e = body["outputs"]["p_e_{SM_2}"]
    assert p_e[-1] > p_e[0]


def _read_ndjson(response) -> list[dict]:
    return [json.loads(line) for line in response.text.strip().split("\n") if line.strip()]


def test_emt_live_endpoint():
    r = client.post(
        "/api/presets/wscc9_3sm/emt/live",
        json={
            "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.02, "t_final": 0.1,
            "plot_states": ["dw_r_{SM_2}"], "plot_inputs": ["V_ref_{SM_1}"], "plot_outputs": ["p_e_{SM_1}"],
        },
    )
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("application/x-ndjson")
    lines = _read_ndjson(r)
    assert len(lines) > 1  # a real step-by-step trace, not just the done line
    *steps, last = lines
    assert last == {"done": True, "perturbed": "dw_r_{SM_2}", "perturb_kind": "state", "n_steps": len(steps)}
    assert steps[0]["t"] == 0.0
    assert steps[0]["states"]["dw_r_{SM_2}"] == pytest.approx(0.02, abs=1e-6)
    assert "V_ref_{SM_1}" in steps[0]["inputs"]
    assert "p_e_{SM_1}" in steps[0]["outputs"]
    # t strictly increases step to step.
    assert all(b["t"] > a["t"] for a, b in zip(steps, steps[1:]))


def test_emt_live_default_plot_states_matches_one_shot_default():
    # Same "empty plot_states -> every dw_r_* state" default as the
    # one-shot endpoint (test_emt_endpoint above) -- the live route must
    # resolve it identically, not silently diverge.
    r = client.post(
        "/api/presets/wscc9_3sm/emt/live",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.02, "t_final": 0.1},
    )
    assert r.status_code == 200
    lines = _read_ndjson(r)
    assert set(lines[0]["states"].keys()) == {"dw_r_{SM_1}", "dw_r_{SM_2}", "dw_r_{SM_3}"}
    assert lines[0]["inputs"] == {}
    assert lines[0]["outputs"] == {}


def test_emt_live_t_final_out_of_bounds_422s():
    # Validated synchronously before any streaming starts -- a normal 422
    # JSON body, not a broken/empty stream.
    r = client.post(
        "/api/presets/wscc9_3sm/emt/live",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 999},
    )
    assert r.status_code == 422
    assert "t_final" in r.json()["detail"]


def test_emt_live_unknown_perturb_state_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt/live",
        json={"perturb_kind": "state", "perturb_name": "not_a_real_state", "t_final": 0.1},
    )
    assert r.status_code == 422


def test_emt_bad_perturb_kind_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "bogus", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3},
    )
    assert r.status_code == 422


def test_emt_unknown_perturb_input_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "input", "perturb_name": "not_a_real_input", "t_final": 0.3},
    )
    assert r.status_code == 422


def test_emt_custom_dt():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3, "dt": 0.01},
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["t"]) == 31  # 0.3/0.01 + 1
    assert body["dt"] == pytest.approx(0.01, rel=1e-6)


def test_emt_dt_too_small_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 2.0, "dt": 1e-5},
    )
    assert r.status_code == 422


def test_emt_explicit_plot_states():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={
            "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3,
            "plot_states": ["theta_{SM_1}", "theta_{SM_2}"],
        },
    )
    assert r.status_code == 200
    assert set(r.json()["series"].keys()) == {"theta_{SM_1}", "theta_{SM_2}"}


def test_emt_unknown_plot_state_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3, "plot_states": ["not_a_real_state"]},
    )
    assert r.status_code == 422


def test_emt_plot_inputs_and_outputs():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={
            "perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3,
            "plot_inputs": ["P_ref_{SM_1}"], "plot_outputs": ["p_e_{SM_1}", "V_t_{SM_1}"],
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert list(body["inputs"].keys()) == ["P_ref_{SM_1}"]
    assert set(body["outputs"].keys()) == {"p_e_{SM_1}", "V_t_{SM_1}"}
    assert len(body["outputs"]["p_e_{SM_1}"]) == 200


def test_emt_unknown_plot_output_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 0.3, "plot_outputs": ["not_a_real_output"]},
    )
    assert r.status_code == 422


def test_emt_unknown_perturb_state_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "not_a_real_state", "t_final": 0.3},
    )
    assert r.status_code == 422


def test_emt_ambiguous_perturb_state_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_", "t_final": 0.3},
    )
    assert r.status_code == 422


def test_emt_t_final_out_of_bounds_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "t_final": 10.0},
    )
    assert r.status_code == 422


def test_emt_solver_failure_is_a_clean_422_not_a_500(monkeypatch):
    # Originally found live: a 0.05 pu perturbation on dw_r_{SM_2} pushed the
    # coupled Newton solve past convergence, and simulate() raising a
    # RuntimeError wasn't caught by run_emt -- an unhandled 500. Fixed by
    # catching it there. Later, timedomain/emt.py's solve_algebraic gained an
    # hybr-then-lm fallback (see its docstring: needed once SMIB's tighter
    # coupling produced a ~0.23 pu initial-guess offset hybr alone couldn't
    # reach), which is robust enough that no perturbation on any preset
    # tried since reliably reproduces this failure through the real solver
    # any more -- a genuine improvement, but it means the *specific input*
    # this test used to rely on no longer proves anything. Testing the
    # error-handling branch itself (not "can we still find a failing input")
    # is what actually matters here, so this mocks simulate() directly.
    import g2elin_api.analysis as analysis_module

    def _boom(*args, **kwargs):
        raise RuntimeError("algebraic Newton solve failed: contrived for this test")

    monkeypatch.setattr(analysis_module, "simulate", _boom)
    r = client.post(
        "/api/presets/wscc9_3sm/emt",
        json={"perturb_kind": "state", "perturb_name": "dw_r_{SM_2}", "perturb_offset": 0.05, "t_final": 0.3},
    )
    assert r.status_code == 422
    assert "did not converge" in r.json()["detail"]


def test_roa_endpoint():
    r = client.post(
        "/api/presets/wscc9_3sm/roa",
        json={
            "axis_x_state": "theta_{SM_2}", "axis_x_range": 0.3,
            "axis_y_state": "dw_r_{SM_2}", "axis_y_range": 0.02,
            "grid_n": 3, "t_final": 0.4, "t_early": 0.1,
        },
    )
    assert r.status_code == 200
    body = r.json()
    assert len(body["in_roa"]) == 3 and all(len(row) == 3 for row in body["in_roa"])
    # The center grid point is the zero-perturbation baseline against itself
    # -- trivially "trending toward" (distance 0 at both checkpoints), same
    # sanity check tests/test_roa.py runs directly against the core module.
    assert body["in_roa"][1][1] is True
    assert body["failed"][1][1] is False
    assert body["early_distance"][1][1] == pytest.approx(0.0, abs=1e-9)


def test_roa_grid_n_out_of_bounds_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/roa",
        json={"axis_x_state": "theta_{SM_2}", "axis_y_state": "dw_r_{SM_2}", "grid_n": 20, "t_final": 0.4, "t_early": 0.1},
    )
    assert r.status_code == 422


def test_roa_t_early_not_less_than_t_final_422s():
    r = client.post(
        "/api/presets/wscc9_3sm/roa",
        json={"axis_x_state": "theta_{SM_2}", "axis_y_state": "dw_r_{SM_2}", "grid_n": 3, "t_final": 0.3, "t_early": 0.3},
    )
    assert r.status_code == 422
