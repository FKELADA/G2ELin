"""Per-unit parameter overrides (``DerUnit.params``): applied on top of the
derived defaults by the operating point (so modal analysis and EMT both see
them), reported by /topology, and validated by name.
"""

from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from g2elin_api.main import app
from g2elin_core.network.presets import wscc9_1sm_1gfm_1gfl
from g2elin_core.operating_point import gfm_params, overridable_param_keys
from g2elin_core.powerflow import run_power_flow
from g2elin_core.operating_point import compute_operating_point

client = TestClient(app)


def _net(**gfm_overrides) -> dict:
    net = wscc9_1sm_1gfm_1gfl().model_dump()
    for der in net["der_units"]:
        if der["unit_type"] == "gfm":
            der["params"] = dict(gfm_overrides)
    return net


def test_overridable_keys():
    assert {"KpCL", "KiCL", "KpVL", "mp", "wf", "Rf"} <= overridable_param_keys("gfm")
    assert {"Kppll", "Kpd", "Kiq"} <= overridable_param_keys("gfl")
    assert {"H", "Ka", "K_PSS"} <= overridable_param_keys("sm")
    assert "wb" not in overridable_param_keys("sm")  # set by the network frequency
    assert overridable_param_keys("infinite_bus") == frozenset()


def test_override_reaches_the_operating_point():
    from g2elin_core.network.schema import Network

    net = Network(**_net(KpCL=0.123))
    op = compute_operating_point(net, run_power_flow(net))
    (gfm_op,) = op.gfm_ops.values()
    assert gfm_op.p["KpCL"] == 0.123
    # Everything not overridden keeps its default.
    bus = next(d.bus for d in net.der_units if d.unit_type.value == "gfm")
    defaults = gfm_params(sn_mva=net.sn_mva, f_hz=net.f_hz, un_kv=net.bus(bus).vn_kv,
                          rt_pu=net.transformers[0].r_pu, lt_pu=net.transformers[0].x_pu)
    assert gfm_op.p["KiCL"] == pytest.approx(defaults["KiCL"])


def test_override_changes_the_eigenvalues():
    base = client.post("/api/network/modal", json={"network": _net()}).json()
    # A much slower GFM voltage loop moves the closed-loop modes.
    tuned = client.post("/api/network/modal", json={"network": _net(KpVL=0.1 * 0.0035, KiVL=0.1 * 0.3)}).json()
    assert base["n_states"] == tuned["n_states"]
    a = sorted((m["real"], m["imag"]) for m in base["modes"])
    b = sorted((m["real"], m["imag"]) for m in tuned["modes"])
    assert a != b


def test_topology_reports_effective_and_default_params():
    r = client.post("/api/network/topology", json={"network": _net(KpCL=0.5)}).json()
    info = next(n["der_info"] for n in r["nodes"] if n["der_info"] and n["der_info"]["unit_type"] == "gfm")
    assert info["control_params"]["KpCL"] == 0.5
    assert info["control_params_default"]["KpCL"] != 0.5


def test_unknown_override_is_flagged_and_422():
    net = _net(NotAParam=1.0)
    v = client.post("/api/network/validate", json={"network": net}).json()
    assert not v["ok"]
    assert any("NotAParam" in i["message"] for i in v["issues"])
    r = client.post("/api/network/modal", json={"network": net})
    assert r.status_code == 422
    assert "NotAParam" in r.json()["detail"]
    # Power flow doesn't use dynamic parameters, so it still runs.
    assert client.post("/api/network/powerflow", json={"network": net}).json()["converged"]
