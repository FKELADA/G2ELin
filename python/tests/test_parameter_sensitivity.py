"""Tracing a mode back to the physical parameter behind it.

The eigenvalue sensitivity says which entries of A a mode cares about; this
says which *parameters* those entries are made of, which is the thing anyone
can actually change. The worked example throughout is the one the maths
predicts: every entry of a machine's rotor-speed row carries 1/(2H), so an
electromechanical mode must come back to H.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

from g2elin_core import reduction as R
from g2elin_core.interconnect import assembly_parts
from g2elin_core.interconnect.network_assembly import build_blocks_and_wiring
from g2elin_core.modal import analyze, classify_modes, parameter_sensitivity
from g2elin_core.modal.parameters import _EntryReach, _block_index
from g2elin_core.components.sm import linearize_sm, sm_dae
from g2elin_core.network.presets import kundur_two_area, wscc9_3sm
from g2elin_core.operating_point import (
    REBASED_PARAM_KEYS, compute_operating_point, rebase_params,
)
from g2elin_core.pipeline import linear_components, linearize_network
from g2elin_core.powerflow import run_power_flow


@pytest.fixture(scope="module")
def wscc9():
    net = wscc9_3sm()
    result = run_power_flow(net)
    system = linearize_network(net, result)
    modal = analyze(system.A, system.state_names)
    return net, result, system, modal


def slowest_electromechanical(modal) -> int:
    kinds = classify_modes(modal)
    cand = [
        (abs(modal.eigenvalues[j]) / (2 * math.pi), j)
        for j in range(len(modal.eigenvalues))
        if modal.eigenvalues[j].imag > 0 and kinds[j].category == R.SYNCHRONISATION
    ]
    assert cand, "no synchronisation mode to test with"
    return sorted(cand)[-1][1]


# --- the claim the feature rests on -------------------------------------------


def test_a_parameter_is_detected_exactly_where_it_appears_symbolically():
    """The whole approach assumes "dA_ij/dp is non-zero" means "p appears in
    A_ij". The MATLAB toolbox printed its symbolic A entries, so that can be
    checked directly rather than argued."""
    from compare_matlab_symbolic import FUNCTIONS, IDENTIFIER, SYMBOLIC, parse_matlab_a

    if not SYMBOLIC.exists():
        pytest.skip("the MATLAB toolbox's symbolic output isn't in this checkout")

    dae = sm_dae(True)
    entries = parse_matlab_a(SYMBOLIC / "A_SG0.txt")
    params = list(dae.param_syms())
    rng = np.random.default_rng(7)
    base = {
        s: float(rng.uniform(0.4, 1.6))
        for s in set(params) | set(dae.state_syms) | set(dae.alg_syms) | set(dae.input_syms)
    }
    reference = dae.linearize(base).A
    floor = np.abs(reference).max() * 1e-7

    for sym in params:
        step = abs(base[sym]) * 1e-6
        up, down = dict(base), dict(base)
        up[sym] += step
        down[sym] -= step
        derivative = (dae.linearize(up).A - dae.linearize(down).A) / (2 * step)
        detected = np.abs(derivative) > floor
        for (i, j), text in entries.items():
            symbolic = sym.name in (set(IDENTIFIER.findall(text)) - FUNCTIONS)
            assert bool(detected[i - 1, j - 1]) == symbolic, (
                f"{sym.name} in A[{i}][{j}]: symbolic says {symbolic}, "
                f"derivative says {bool(detected[i - 1, j - 1])}"
            )


def test_inertia_lives_only_in_the_rotor_speed_row():
    """Why an electromechanical mode traces back to H: the swing equation
    divides by 2H, so H reaches A only through that one row."""
    from compare_matlab_symbolic import SYMBOLIC, parse_matlab_a

    if not SYMBOLIC.exists():
        pytest.skip("the MATLAB toolbox's symbolic output isn't in this checkout")
    entries = parse_matlab_a(SYMBOLIC / "A_SG0.txt")
    rows = {i for (i, _), text in entries.items() if re.search(r"\bH\b", text)}
    names = sm_dae(True).state_names
    assert rows, "H should appear somewhere"
    assert {names[i - 1] for i in rows} == {"dw_r"}


# --- the feature --------------------------------------------------------------


def test_inertia_is_the_top_parameter_for_an_electromechanical_mode(wscc9):
    net, result, _, modal = wscc9
    mode = slowest_electromechanical(modal)
    report = parameter_sensitivity(net, result, modal, mode)
    assert report.effects, "nothing was scanned"
    top = report.effects[0]
    assert top.parameter == "H", f"expected H to lead, got {top.parameter}"
    # ...and it should be the inertia of a machine that actually swings in it.
    speed_state = "dw_r_{" + top.unit_label + "}"
    index = modal.state_names.index(speed_state)
    assert abs(modal.participation[index, mode]) > 0.05


def test_the_most_sensitive_entry_names_the_parameter_inside_it(wscc9):
    """The bridge from "which entry of A" to "which parameter": the entry the
    mode is most sensitive to should report the parameters it is built from,
    and for a rotor-speed row that means H."""
    net, result, _, modal = wscc9
    report = parameter_sensitivity(net, result, modal, slowest_electromechanical(modal))
    assert report.entries
    top = report.entries[0]
    assert top.row_state.startswith("dw_r"), f"expected a rotor-speed row, got {top.row_state}"
    assert any(p.endswith(".H") for p in top.parameters), top.parameters


def test_an_exciter_mode_comes_back_to_exciter_parameters(wscc9):
    """The counterpart: a mode the AVR owns should not report H."""
    net, result, _, modal = wscc9
    kinds = classify_modes(modal)
    control = [
        j for j in range(len(modal.eigenvalues))
        if modal.eigenvalues[j].imag > 0 and kinds[j].category == R.CONTROL
        and 0.5 <= abs(modal.eigenvalues[j]) / (2 * math.pi) <= 1.0
    ]
    if not control:
        pytest.skip("this network has no exciter mode in the band")
    report = parameter_sensitivity(net, result, modal, control[0])
    leaders = {e.parameter for e in report.effects[:5]}
    assert leaders & {"Ka", "Ta", "Te", "Ke", "Tr", "Kfd", "Tfd", "Rfd", "Lfd", "Lad"}, leaders
    assert report.effects[0].parameter != "H"


def test_the_chain_rule_agrees_with_actually_moving_the_parameter(wscc9):
    """dlambda/dp is assembled from dlambda/dA and dA/dp rather than measured.
    Perturbing the parameter for real has to give the same answer.

    Two things this test has to get right, both of which it got wrong first
    time. The perturbation must be *small*: this is a first-order derivative,
    and on this mode a 1% change in H already disagrees by 10% while 0.01%
    agrees to 0.14%. And the moved eigenvalue must be tracked by eigenvector
    alignment, not by proximity -- WSCC-9's two local modes sit 0.045 apart,
    so the nearest eigenvalue to the old one is soon the *other* mode.
    """
    net, result, _, modal = wscc9
    mode = slowest_electromechanical(modal)
    report = parameter_sensitivity(net, result, modal, mode)
    effect = report.effects[0]

    bumped = effect.value * 1.0001
    units = [
        d.model_copy(update={"params": {**d.params, effect.parameter: bumped}})
        if d.id == effect.unit else d
        for d in net.der_units
    ]
    system = linearize_network(net.model_copy(update={"der_units": units}), result)
    after = analyze(system.A, system.state_names)

    before = modal.right_eigenvectors[:, mode]
    alignment = np.abs(before.conj() @ after.right_eigenvectors)
    moved = int(np.argmax(alignment))
    assert alignment[moved] > 0.99, "the mode could not be tracked across the perturbation"

    measured = (after.eigenvalues[moved] - modal.eigenvalues[mode]) / (bumped - effect.value)
    assert abs(measured - effect.d_lambda) < 0.01 * abs(effect.d_lambda)


def test_it_can_be_narrowed_to_one_unit_and_one_parameter(wscc9):
    net, result, _, modal = wscc9
    mode = slowest_electromechanical(modal)
    unit = net.der_units[0].id
    report = parameter_sensitivity(net, result, modal, mode, units=[unit], parameters=["H"])
    assert [(e.unit, e.parameter) for e in report.effects] == [(unit, "H")]


def test_a_mismatched_modal_result_is_refused(wscc9):
    """Passing a modal result from a different model would silently produce
    nonsense, since the sensitivity matrix would be indexed against the wrong
    states."""
    net, result, _, _ = wscc9
    other = kundur_two_area()
    other_system = linearize_network(other, run_power_flow(other))
    other_modal = analyze(other_system.A, other_system.state_names)
    with pytest.raises(ValueError, match="not the same model"):
        parameter_sensitivity(net, result, other_modal, 0)


# --- the two things the fast scan is allowed to assume -------------------------
#
# Scanning a thousand parameters is only affordable because perturbing one
# rebuilds one unit rather than the network (see the module docstring). That
# buys its speed with two assumptions, and both are checked here against the
# network actually being relinearised.


@pytest.fixture(scope="module")
def kundur():
    net = kundur_two_area()
    result = run_power_flow(net)
    system = linearize_network(net, result)
    return net, result, analyze(system.A, system.state_names)


def true_d_lambda(net, result, modal, mode, der_id, name, raw_value, from_value):
    """d(lambda)/dp measured by relinearising the whole network.

    Tracked by eigenvector alignment rather than by proximity: with several
    machines the nearest eigenvalue to the old one is soon a different mode.
    """
    units = [d.model_copy(update={"params": {**d.params, name: raw_value}})
             if d.id == der_id else d for d in net.der_units]
    system = linearize_network(net.model_copy(update={"der_units": units}), result)
    after = analyze(system.A, system.state_names)
    alignment = np.abs(modal.right_eigenvectors[:, mode].conj() @ after.right_eigenvectors)
    moved = int(np.argmax(alignment))
    assert alignment[moved] > 0.99, "the mode could not be tracked across the perturbation"
    return after.eigenvalues[moved] - modal.eigenvalues[mode]


def test_a_machine_on_its_own_rating_is_perturbed_on_the_right_base(kundur):
    """Kundur's machines are rated 900 MVA on a 100 MVA network, so their
    overrides are read per unit of 900 and converted (operating_point.py's
    rebasing). The scan reports the parameter as the *component* sees it,
    on the network base, and must perturb it on that same base -- reading
    the value from one base and moving it on the other would scale every
    impedance sensitivity by 9.
    """
    net, result, modal = kundur
    mode = slowest_electromechanical(modal)
    unit = net.der_units[0]
    assert unit.sn_mva is not None and unit.sn_mva != net.sn_mva, "this test needs a rebased unit"

    report = parameter_sensitivity(net, result, modal, mode, units=[unit.id],
                                   parameters=["H", "Ll", "Lad"])
    assert {e.parameter for e in report.effects} == {"H", "Ll", "Lad"}
    for effect in report.effects:
        assert effect.parameter in REBASED_PARAM_KEYS, "the point of the test is a rebased key"
        bumped = effect.value * 1.0001
        # Back to the unit's own base, which is what an override is read in.
        raw = rebase_params({effect.parameter: bumped},
                            from_mva=net.sn_mva, to_mva=unit.sn_mva)[effect.parameter]
        measured = true_d_lambda(net, result, modal, mode, unit.id, effect.parameter,
                                 raw, effect.value) / (bumped - effect.value)
        assert abs(measured - effect.d_lambda) < 0.01 * abs(effect.d_lambda), effect.parameter


def test_rebuilding_one_unit_gives_the_same_component_as_rebuilding_the_network(kundur):
    """The scan perturbs a parameter by rebuilding that unit's operating
    point alone (``RebuiltWithParams``). Everything an operating point is
    built from besides the parameters -- terminal voltage, dispatch -- is
    fixed by the power flow, so for an ordinary unit that has to be exact.
    """
    net, result, _ = kundur
    op = compute_operating_point(net, result)
    der = next(d for d in net.der_units if d.bus_type.value != "slack")
    base = op.sm_ops[der.id].p
    modes = net.unit_modes(der)

    for name in ("H", "Ra", "Ll", "Laq", "Ka"):
        moved = base[name] * 1.10
        fast = linearize_sm(op.sm_ops[der.id].with_params({**base, name: moved}), modes=modes)
        raw = rebase_params({name: moved}, from_mva=net.sn_mva, to_mva=der.sn_mva)[name]
        units = [d.model_copy(update={"params": {**d.params, name: raw}})
                 if d.id == der.id else d for d in net.der_units]
        slow = linear_components(net.model_copy(update={"der_units": units}),
                                 result)["der_components"][der.id]
        for label, a, b in (("A", fast.A, slow.A), ("B", fast.B, slow.B),
                            ("C", fast.C, slow.C), ("D", fast.D, slow.D)):
            assert np.allclose(a, b, rtol=0, atol=1e-12), f"{name}: {label} differs"


@pytest.mark.parametrize("frame_follows_slack", [False, True])
def test_the_slack_machine_moving_the_reference_angle_does_not_change_the_answer(
    frame_follows_slack,
):
    """The one unit the assumption above does not hold exactly for. The slack
    machine's rotor angle *is* the reference angle every other unit is built
    against, so Ra, Ll and Laq shift the whole frame -- and rebuilding that
    one unit does not propagate the shift to the others.

    It does not need to. The reference angle is a choice of frame, not a
    physical quantity, and an eigenvalue cannot depend on it. These are the
    three parameters where it would show up if that were wrong, under both
    frame conventions -- the slack carrying the frame (MATLAB-compatible) and
    an explicit frame component of its own.

    On WSCC-9 rather than Kundur because the reference here is a finite
    difference of the network relinearised, and Kundur's inter-area mode is
    too ill-conditioned to difference cleanly: shrinking the step there makes
    the *measured* value worse, not better.
    """
    net = wscc9_3sm().model_copy(update={"frame_follows_slack": frame_follows_slack})
    result = run_power_flow(net)
    system = linearize_network(net, result)
    modal = analyze(system.A, system.state_names)
    mode = slowest_electromechanical(modal)
    slack = next(d for d in net.der_units if d.bus_type.value == "slack")
    assert slack.sn_mva is None, "this network is on one base, so no rebasing to undo"

    report = parameter_sensitivity(net, result, modal, mode, units=[slack.id],
                                   parameters=["Ra", "Ll", "Laq"])
    assert len(report.effects) == 3
    for effect in report.effects:
        bumped = effect.value * 1.0001
        measured = true_d_lambda(net, result, modal, mode, slack.id, effect.parameter,
                                 bumped, effect.value) / (bumped - effect.value)
        assert abs(measured - effect.d_lambda) < 0.01 * abs(effect.d_lambda), effect.parameter


def test_entry_reach_agrees_with_the_matrix_it_avoids_building(wscc9):
    """Which parameters an entry of A is built from needs ``|dA_tot/dp|`` at
    a handful of entries -- not the whole matrix, which on a 118-bus model is
    27 MB per parameter and was, briefly, allocated that many times.

    ``_EntryReach`` reads those entries straight out of the component's own
    blocks. Against the matrix it is standing in for, it has to be exact.
    """
    net, result, _, _ = wscc9
    op = compute_operating_point(net, result)
    blocks, wiring = build_blocks_and_wiring(net, **linear_components(net, result, op=op))
    parts = assembly_parts(blocks, wiring)

    n = parts.A_tot.shape[0]
    rng = np.random.default_rng(0)
    cells = [(int(i), int(j)) for i, j in rng.integers(0, n, size=(40, 2))]
    reach = _EntryReach(parts, cells)

    GE = parts.topology.G @ parts.E_ol
    P, Q = GE @ parts.C_ol, parts.B_ol @ GE

    for der in net.der_units:
        index = _block_index(blocks, der.id, net)
        xs, us, ys = parts.block_slices(index)
        base, modes = op.sm_ops[der.id].p, net.unit_modes(der)
        for name in ("H", "Ra", "Ka"):
            step = base[name] * 1e-4
            up = linearize_sm(op.sm_ops[der.id].with_params({**base, name: base[name] + step}),
                              modes=modes)
            dn = linearize_sm(op.sm_ops[der.id].with_params({**base, name: base[name] - step}),
                              modes=modes)
            dA, dB = (up.A - dn.A) / (2 * step), (up.B - dn.B) / (2 * step)
            dC, dD = (up.C - dn.C) / (2 * step), (up.D - dn.D) / (2 * step)

            full = np.zeros_like(parts.A_tot)
            full[xs, xs] += dA
            full[xs, :] += dB @ P[us, :]
            full[:, xs] += Q[:, ys] @ dC
            full += Q[:, ys] @ dD @ P[us, :]

            expected = np.array([abs(full[i, j]) for i, j in cells])
            assert np.allclose(reach.at(index, dA, dB, dC, dD), expected, rtol=0, atol=1e-15)
