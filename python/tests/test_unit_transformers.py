"""Which transformer impedance a unit's dynamic model uses: by default its
own transformer -- the same element the power flow uses, converted to the
system base -- or, in the MATLAB-compatible mode, the network's first
transformer for every unit (script_generic.m's Y_TR(1,:) convention)."""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.modal import analyze
from g2elin_core.network.presets import cigre_interconnected_1sm_1gfm_1gfl, wscc9_3sm
from g2elin_core.network.validation import validate_network
from g2elin_core.operating_point import compute_operating_point, unit_transformer_rx
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow


def test_each_unit_uses_its_own_transformer_by_default():
    net = wscc9_3sm()  # unit transformers 0.0576 / 0.0625 / 0.0586 pu
    own = {d.id: next(t for t in net.transformers if t.lv_bus == d.bus).x_pu for d in net.der_units}
    for der in net.der_units:
        assert unit_transformer_rx(net, der) == (pytest.approx(0.0), pytest.approx(own[der.id]))


def test_matlab_compatible_mode_uses_the_first_transformer_for_all():
    net = wscc9_3sm().model_copy(update={"units_use_first_transformer": True})
    for der in net.der_units:
        assert unit_transformer_rx(net, der) == (0.0, 0.0576)


def test_rating_is_converted_to_the_system_base():
    net = wscc9_3sm()
    tr = net.transformers[1].model_copy(update={"sn_mva": 50.0})  # half the 100 MVA system base
    net = net.model_copy(update={"transformers": [net.transformers[0], tr, net.transformers[2]]})
    der = next(d for d in net.der_units if d.bus == tr.lv_bus)
    assert unit_transformer_rx(net, der) == (pytest.approx(0.0), pytest.approx(2 * tr.x_pu))


def test_the_model_uses_it():
    net = wscc9_3sm()
    result = run_power_flow(net)
    op = compute_operating_point(net, result)
    lts = sorted(sm.p["Lt"] for sm in op.sm_ops.values())
    assert lts == pytest.approx([0.0576, 0.0586, 0.0625])


def test_both_modes_are_stable_and_close_on_wscc():
    base = wscc9_3sm()
    eig = {}
    for compat in (False, True):
        net = base.model_copy(update={"units_use_first_transformer": compat})
        r = run_power_flow(net)
        s = linearize_network(net, r)
        eig[compat] = analyze(s.A, s.state_names).eigenvalues
        assert eig[compat].real.max() < 1e-6
    # Same model size, and the electromechanical modes move only slightly.
    assert len(eig[False]) == len(eig[True])
    em = lambda e: np.sort(np.abs(e[(np.abs(e.imag) > 2 * np.pi * 0.5) & (np.abs(e.imag) < 2 * np.pi * 3)]))
    assert em(eig[False]) == pytest.approx(em(eig[True]), rel=0.05)


def test_rt_lt_override_is_flagged_outside_compat_mode():
    net = wscc9_3sm()
    net.der_units[1].params = {"Lt": 0.1}
    msgs = [i.message for i in validate_network(net) if i.severity == "warning"]
    assert any("Lt" in m and "transformer" in m for m in msgs)
    compat = net.model_copy(update={"units_use_first_transformer": True})
    assert not any("Lt" in i.message for i in validate_network(compat))


def test_cigre_interconnected_transformers_are_uniform():
    net = cigre_interconnected_1sm_1gfm_1gfl()
    assert len({(t.r_pu, t.x_pu, t.sn_mva) for t in net.transformers}) == 1
