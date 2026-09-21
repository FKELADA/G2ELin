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




# --- DerUnit.sn_mva: machine data given on the machine's own base -------------
import numpy as np

from g2elin_core.network.presets import wscc9_3sm
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
def _rated(**params):
    net = wscc9_3sm()
    net.der_units[1].sn_mva = 900.0
    net.der_units[1].params = params
    return net


def test_rated_unit_equals_hand_converted_params():
    """A 900 MVA machine's own data must give exactly the model you get by
    converting it to the 100 MVA network base by hand."""
    machine = dict(H=6.5, Ra=0.0025, Ll=0.2, Lad=1.8, Laq=1.7, Lfd=0.165, KD=2.0)
    k = 100.0 / 900.0
    by_hand = {n: (v / k if n in ("H", "KD") else v * k) for n, v in machine.items()}

    rated = _rated(**machine)
    plain = wscc9_3sm()
    plain.der_units[1].params = by_hand

    A_rated = linearize_network(rated, run_power_flow(rated)).A
    A_plain = linearize_network(plain, run_power_flow(plain)).A
    assert np.allclose(A_rated, A_plain, rtol=0, atol=0)
    assert by_hand["Lad"] == pytest.approx(0.2)     # 1.8 pu on 900 MVA is 0.2 pu on 100 MVA
    assert by_hand["H"] == pytest.approx(58.5)      # 6.5 s of stored energy per 900 MVA


def test_rating_actually_changes_the_model():
    """Guards against the rebasing silently doing nothing."""
    machine = dict(H=6.5, Lad=1.8)
    rated = _rated(**machine)
    as_is = wscc9_3sm()
    as_is.der_units[1].params = dict(machine)
    A_rated = linearize_network(rated, run_power_flow(rated)).A
    A_as_is = linearize_network(as_is, run_power_flow(as_is)).A
    assert not np.allclose(A_rated, A_as_is, rtol=1e-6, atol=1e-6)


def test_rating_equal_to_the_network_base_is_a_no_op():
    net = wscc9_3sm()
    net.der_units[1].sn_mva = net.sn_mva
    net.der_units[1].params = dict(H=6.5, Lad=1.8)
    plain = wscc9_3sm()
    plain.der_units[1].params = dict(H=6.5, Lad=1.8)
    assert np.allclose(linearize_network(net, run_power_flow(net)).A,
                       linearize_network(plain, run_power_flow(plain)).A, rtol=0, atol=0)


def test_rebase_params_leaves_tuning_alone():
    from g2elin_core.operating_point import rebase_params

    out = rebase_params(dict(Lad=1.8, H=6.5, Ka=300.0, Tr=0.02), from_mva=900.0, to_mva=100.0)
    assert out["Lad"] == pytest.approx(0.2) and out["H"] == pytest.approx(58.5)
    assert out["Ka"] == 300.0 and out["Tr"] == 0.02       # gains/time constants are not machine data


def test_validation_warns_when_a_gain_is_given_on_a_rated_unit():
    net = _rated(H=6.5, Ka=250.0)
    warnings_ = [i for i in validate_network(net) if i.severity == "warning" and "Ka" in i.message]
    assert warnings_ and "taken as given" in warnings_[0].message
    assert not [i for i in validate_network(net) if i.severity == "error"]


def test_validation_warns_when_a_rating_has_no_effect():
    net = wscc9_3sm()
    net.der_units[1].sn_mva = 900.0
    msgs = [i.message for i in validate_network(net) if i.severity == "warning"]
    assert any("has no effect" in m for m in msgs)


def test_topology_reports_the_rebased_parameters_not_the_raw_ones():
    """/topology feeds the inspector's parameter panel. It built its own
    defaults+overrides merge and so skipped the rebasing, showing a 900 MVA
    machine's H as the 6.5 s typed in rather than the 58.5 s the model uses."""
    net = wscc9_3sm().model_dump()
    net["der_units"][1]["sn_mva"] = 900.0
    net["der_units"][1]["params"] = {"H": 6.5, "Lad": 1.8}
    nodes = client.post("/api/network/topology", json={"network": net}).json()["nodes"]
    rated = [n for n in nodes if n.get("der_info") and n["der_info"].get("control_params", {}).get("H") != 5.0]
    assert len(rated) == 1
    params = rated[0]["der_info"]["control_params"]
    assert params["H"] == pytest.approx(58.5)
    assert params["Lad"] == pytest.approx(0.2)
    # The defaults shown beside them stay on the network base.
    assert rated[0]["der_info"]["control_params_default"]["H"] == pytest.approx(5.0)
