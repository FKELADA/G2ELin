"""``/api/network/*`` -- the arbitrary-Network-JSON counterparts of
``test_api.py``'s preset-id-keyed endpoints. Mirrors that file's structure
and assertions, POSTing a preset's own ``Network.model_dump()`` as the body
instead of referencing it by id, plus dedicated cases for the pydantic
422-on-invalid-network path (not exercised by any preset-id endpoint, since
every preset is already known-valid).
"""

from __future__ import annotations

import json

from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core.network.presets import wscc9_1gfm_2gfl, wscc9_3sm

client = TestClient(app)


def _wscc_body() -> dict:
    return wscc9_3sm().model_dump()


def test_network_topology():
    r = client.post("/api/network/topology", json={"network": _wscc_body()})
    assert r.status_code == 200
    body = r.json()
    assert len(body["nodes"]) == 9
    assert len(body["edges"]) == 9  # 6 lines + 3 transformers


def test_network_powerflow_converges():
    r = client.post("/api/network/powerflow", json={"network": _wscc_body()})
    assert r.status_code == 200
    body = r.json()
    assert body["converged"] is True
    assert len(body["buses"]) == 9


def test_network_modal():
    r = client.post("/api/network/modal", json={"network": _wscc_body()})
    assert r.status_code == 200
    body = r.json()
    assert body["n_states"] > 0
    assert len(body["modes"]) == body["n_states"]


def test_network_modal_zero_q_load_is_422_not_500():
    # A load left at the web UI editor's own old "+ Add load" default
    # (p_mw=0, q_mvar=0) -- power flow treats it as a no-op bus and
    # converges, but its constant-impedance equivalent divides by its own
    # apparent power. Must come back as an actionable 422, not a bare
    # unhandled-exception 500.
    net = _wscc_body()
    net["loads"].append({"bus": net["buses"][0]["id"], "p_mw": 0.0, "q_mvar": 0.0, "name": ""})

    r = client.post("/api/network/powerflow", json={"network": net})
    assert r.status_code == 200
    assert r.json()["converged"] is True

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "zero reactive power" in r.json()["detail"]


def test_network_modal_purely_resistive_load_is_422_not_500():
    # A load with real, nonzero p_mw but q_mvar=0 exactly (purely resistive)
    # -- *not* the "zero apparent power" case. Power flow is fine with it,
    # but components/load.py's own dynamic model uses this load's reactance
    # x_pu = z_pu*sin(acos(p_pu/s_pu)) as a divisor, and sin(acos(+-1)) = 0
    # whenever q_mvar=0 regardless of how large p_mw is.
    net = _wscc_body()
    net["loads"].append({"bus": net["buses"][0]["id"], "p_mw": 5.0, "q_mvar": 0.0, "name": ""})

    r = client.post("/api/network/powerflow", json={"network": net})
    assert r.status_code == 200
    assert r.json()["converged"] is True

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "zero reactive power" in r.json()["detail"]


def test_network_validate_zero_q_load_is_flagged():
    net = _wscc_body()
    net["loads"].append({"bus": net["buses"][0]["id"], "p_mw": 0.0, "q_mvar": 0.0, "name": ""})

    r = client.post("/api/network/validate", json={"network": net})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is False
    assert any(i["severity"] == "error" and "zero reactive power" in i["message"] for i in body["issues"])


def test_network_modal_sensitivity():
    r = client.post("/api/network/modal/sensitivity", json={"network": _wscc_body(), "mode": 0})
    assert r.status_code == 200
    assert r.json()["mode"] == 0


def test_network_states():
    r = client.post("/api/network/states", json={"network": _wscc_body()})
    assert r.status_code == 200
    body = r.json()
    assert any("dw_r" in n for n in body["state_names"])


def test_network_timeseries():
    r = client.post("/api/network/timeseries", json={"network": _wscc_body()})
    assert r.status_code == 200
    assert len(r.json()["snapshots"]) == 3


def test_network_emt():
    states = client.post("/api/network/states", json={"network": _wscc_body()}).json()
    perturb = next(n for n in states["state_names"] if "dw_r" in n)
    r = client.post(
        "/api/network/emt",
        json={"network": _wscc_body(), "perturb_kind": "state", "perturb_name": perturb, "t_final": 0.2},
    )
    assert r.status_code == 200
    assert r.json()["perturbed"] == perturb


def test_network_emt_live():
    states = client.post("/api/network/states", json={"network": _wscc_body()}).json()
    perturb = next(n for n in states["state_names"] if "dw_r" in n)
    r = client.post(
        "/api/network/emt/live",
        json={"network": _wscc_body(), "perturb_kind": "state", "perturb_name": perturb,
              "perturb_offset": 0.02, "t_final": 0.1},
    )
    assert r.status_code == 200
    lines = [json.loads(line) for line in r.text.strip().split("\n") if line.strip()]
    assert len(lines) > 1
    assert lines[0]["t"] == 0.0
    assert lines[-1] == {"done": True, "perturbed": perturb, "perturb_kind": "state", "n_steps": len(lines) - 1}


def test_network_invalid_body_is_422_with_field_detail():
    # A negative x_pu -- Line.x_pu has Field(gt=0) -- should 422 with
    # pydantic's per-field detail list (not a hand-written HTTPException),
    # confirmed by shape here since that's the exact case
    # web/index.html's formatDetail() exists to render readably.
    net = wscc9_3sm().model_dump()
    net["lines"][0]["x_pu"] = -1.0
    r = client.post("/api/network/topology", json={"network": net})
    assert r.status_code == 422
    detail = r.json()["detail"]
    assert isinstance(detail, list)
    assert any("x_pu" in ".".join(str(p) for p in e["loc"]) for e in detail)


def test_network_missing_slack_is_422():
    net = wscc9_3sm().model_dump()
    net["der_units"][0]["bus_type"] = "pv"  # was slack -- now zero slack units
    r = client.post("/api/network/powerflow", json={"network": net})
    assert r.status_code == 422


def test_network_with_no_transformers_is_422_not_500():
    # Regression: a hand-built network (e.g. from the web UI's
    # drag-and-drop builder) can easily end up with zero transformers if
    # every DER was wired with "Wire: Line" instead of "Wire: Transformer"
    # -- operating_point.compute_operating_point() used to index
    # network.transformers[0] unconditionally, crashing with a bare
    # IndexError (-> unhandled 500) instead of a clear, actionable error.
    # Power flow doesn't care (pandapower has no such requirement), so this
    # specifically exercises modal/EMT, which do.
    net = wscc9_3sm().model_dump()
    net["transformers"] = []
    pf = client.post("/api/network/powerflow", json={"network": net})
    assert pf.status_code == 200 and pf.json()["converged"]

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "no transformers" in r.json()["detail"]

    r = client.post("/api/network/states", json={"network": net})
    assert r.status_code == 422
    assert "no transformers" in r.json()["detail"]


def test_network_der_not_behind_a_transformer_is_422_not_500():
    # Same class of bug as above, but a transformer exists somewhere in the
    # network -- it's just not the one connecting this DER (der.bus isn't
    # any transformer's lv_bus). Used to raise a bare KeyError.
    net = wscc9_3sm().model_dump()
    net["transformers"][0]["lv_bus"] = 99  # slack SM's own bus (7) is no longer any transformer's lv_bus
    net["buses"].append({"id": 99, "name": "stray", "vn_kv": 18.0})

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    detail = r.json()["detail"]
    assert "no transformer connecting" in detail
    assert "id=1" in detail  # the slack SM's own DER id


def test_network_der_transformer_hv_lv_swapped_is_422_with_actionable_message():
    # A transformer *does* connect this DER to the rest of the network, but
    # hv_bus/lv_bus are backward (the DER's own bus is hv_bus instead of
    # lv_bus) -- the natural mistake when drawing the connection in the
    # builder canvas starting *from* the DER. Must be diagnosed specifically
    # ("swapped") rather than reported as "no transformer at all".
    net = wscc9_3sm().model_dump()
    tr = net["transformers"][0]
    tr["hv_bus"], tr["lv_bus"] = tr["lv_bus"], tr["hv_bus"]

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    detail = r.json()["detail"]
    assert "swapped" in detail
    assert "id=1" in detail  # the slack SM's own DER id


def test_network_modal_no_lines_is_422_not_500():
    # A network built from a DER + its own transformer + a load bus, with no
    # Line elements at all (e.g. a single generator behind a transformer
    # serving one local load, no separate feeder). Power flow doesn't care,
    # but every bus's own dynamic model borrows its line-charging
    # susceptance from "line #1" -- with zero lines there's nothing to
    # borrow and the old 0.0 fallback was a real division by zero deep in
    # sympy (TypeError: Cannot convert complex to float), not a clean error.
    net = {
        "name": "no_lines", "f_hz": 60.0, "sn_mva": 100.0,
        "buses": [{"id": 1, "name": "der_bus", "vn_kv": 10.0}, {"id": 2, "name": "grid_bus", "vn_kv": 20.0}],
        "lines": [],
        "transformers": [{"hv_bus": 2, "lv_bus": 1, "r_pu": 0.0, "x_pu": 0.05, "sn_mva": 100.0, "name": ""}],
        "loads": [{"bus": 2, "p_mw": 5.0, "q_mvar": 1.0, "name": ""}],
        "der_units": [
            {"id": 1, "bus": 1, "unit_type": "infinite_bus", "bus_type": "slack", "v_set_pu": 1.0,
             "p_set_mw": 0.0, "q_set_mvar": 0.0, "p_cons_mw": 0.0, "q_cons_mvar": 0.0,
             "controller": None, "xd_pu": None},
        ],
    }
    r = client.post("/api/network/validate", json={"network": net})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is False
    assert any(i["severity"] == "error" and "no Line elements" in i["message"] for i in body["issues"])

    pf = client.post("/api/network/powerflow", json={"network": net})
    assert pf.status_code == 200 and pf.json()["converged"]

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "no Line elements" in r.json()["detail"]


def _minimal_2der_network(bad_transformer_hv_bus: int) -> dict:
    """A 3-bus network (1 grid bus + 2 DER terminal buses), for the two
    DER-topology regression tests below.
    """
    return {
        "name": "der_to_der", "f_hz": 60.0, "sn_mva": 100.0,
        "buses": [
            {"id": 1, "name": "grid", "vn_kv": 230.0},
            {"id": 2, "name": "sm_terminal", "vn_kv": 18.0},
            {"id": 3, "name": "gfm_terminal", "vn_kv": 18.0},
        ],
        "lines": [], "loads": [{"bus": 1, "p_mw": 5.0, "q_mvar": 1.0, "name": ""}],
        "der_units": [
            {"id": 1, "bus": 2, "unit_type": "sm", "bus_type": "slack", "v_set_pu": 1.0, "p_set_mw": 0.0,
             "q_set_mvar": 0.0, "p_cons_mw": 0.0, "q_cons_mvar": 0.0, "controller": None, "xd_pu": 0.2},
            {"id": 2, "bus": 3, "unit_type": "gfm", "bus_type": "pv", "v_set_pu": 1.0, "p_set_mw": 1.0,
             "q_set_mvar": 0.0, "p_cons_mw": 0.0, "q_cons_mvar": 0.0, "controller": "droop", "xd_pu": None},
        ],
        "transformers": [
            {"hv_bus": 1, "lv_bus": 2, "r_pu": 0.0, "x_pu": 0.05, "sn_mva": 100.0, "name": ""},
            {"hv_bus": bad_transformer_hv_bus, "lv_bus": 3, "r_pu": 0.0, "x_pu": 0.05, "sn_mva": 100.0, "name": ""},
        ],
    }


def test_network_der_transformer_to_another_der_bus_is_422_not_500():
    # Regression: connecting one DER's transformer directly to another
    # DER's own bus (instead of a plain grid bus) -- an easy mistake with
    # the drag-and-drop builder's "Wire: Transformer" mode -- used to crash
    # deep inside interconnect/network_assembly.py with a bare KeyError.
    net = _minimal_2der_network(bad_transformer_hv_bus=2)  # GFM's xfmr -> SM's own bus (2), not the grid bus (1)
    pf = client.post("/api/network/powerflow", json={"network": net})
    assert pf.status_code == 200 and pf.json()["converged"]

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "own transformer connects directly to" in r.json()["detail"]


def test_network_line_on_der_bus_is_422_not_500():
    # Regression: a line wired directly to/from a DER's own private bus
    # (instead of only through its own transformer) -- e.g. "Wire: Line"
    # used by mistake on a placed unit -- used to crash the same way.
    net = _minimal_2der_network(bad_transformer_hv_bus=1)  # both transformers correctly on the grid bus
    net["lines"].append(
        {"from_bus": 2, "to_bus": 3, "r_pu": 0.01, "x_pu": 0.1, "b_pu": 0.0, "length_km": 1.0, "name": ""}
    )

    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "own bus" in r.json()["detail"]


def test_network_validate_clean_network_is_ok():
    r = client.post("/api/network/validate", json={"network": _wscc_body()})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True
    assert body["issues"] == []


def test_network_validate_reports_every_issue_at_once():
    # The point of a dedicated validate endpoint: multiple independent
    # problems in one network all come back together, not one at a time.
    net = _minimal_2der_network(bad_transformer_hv_bus=2)  # GFM's xfmr -> SM's own bus
    net["lines"].append(
        {"from_bus": 2, "to_bus": 3, "r_pu": 0.01, "x_pu": 0.1, "b_pu": 0.0, "length_km": 1.0, "name": ""}
    )
    r = client.post("/api/network/validate", json={"network": net})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is False
    assert len(body["issues"]) >= 2
    assert all(i["severity"] == "error" for i in body["issues"])


def test_network_validate_unsupported_slack_is_a_warning_not_an_error():
    # A GFM/GFL slack: power flow works, modal/EMT don't support it
    # yet -- validate() should flag it as a warning, not block "ok".
    r = client.post("/api/network/validate", json={"network": wscc9_1gfm_2gfl().model_dump()})
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] is True  # no errors, just a warning
    assert any(i["severity"] == "warning" and "slack" in i["message"] for i in body["issues"])
