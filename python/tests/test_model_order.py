"""Model-order reduction: the mechanism, the named levels, and -- the part
that actually matters -- that a reduced model still answers the question the
full one does, in the band it claims to cover.

The acceptance criterion throughout is the one the literature uses for a
singular-perturbation reduction: the slow spectrum has to survive. A test
that only checked state counts would pass on a model that had quietly
stopped being the same system.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core import reduction
from g2elin_core.components.base import ALGEBRAIC, DYNAMIC, FROZEN
from g2elin_core.modal import check_adequacy
from g2elin_core.modal.adequacy import _group_of
from g2elin_core.network.presets import cigre_islanded_1sm_2gfm_1gfl, wscc9_3sm
from g2elin_core.network.schema import ModelOptions, Network
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network

ELECTROMECHANICAL_HZ = (0.1, 4.0)


def with_models(network: Network, **kwargs) -> Network:
    return network.model_copy(update={"models": ModelOptions(**kwargs)})


def modes_in_band(system, lo: float, hi: float) -> list[tuple[float, float]]:
    """(frequency, damping %) of each mode in a band, one entry per pair."""
    out = []
    for lam in np.linalg.eigvals(system.A):
        if lam.imag <= 0:
            continue
        hz = abs(lam) / (2 * np.pi)
        if lo <= hz <= hi:
            out.append((hz, -lam.real / abs(lam) * 100))
    return sorted(out)


def linearized(network: Network):
    result = run_power_flow(network)
    assert result.converged
    return linearize_network(network, result)


# --- the mechanism ------------------------------------------------------------


def test_a_fully_dynamic_reduction_is_the_full_model():
    """The default settings must reproduce the pre-reduction model exactly --
    not approximately. Everything else in this file is judged against it."""
    net = wscc9_3sm()
    a = linearized(net)
    b = linearized(with_models(net))
    assert a.state_names == b.state_names
    np.testing.assert_allclose(a.A, b.A, rtol=0, atol=0)


def test_algebraic_states_leave_the_state_vector_and_keep_their_names_out_of_it():
    net = with_models(wscc9_3sm(), network_level="quasi_stationary")
    full, reduced = linearized(wscc9_3sm()), linearized(net)
    assert reduced.A.shape[0] < full.A.shape[0]
    # Bus voltages, line and load currents are gone; machine states remain.
    assert not [n for n in reduced.state_names if n.startswith(("v_{g_d}", "i_{l_d}", "i_{c_d}"))]
    assert [n for n in reduced.state_names if n.startswith("psi_d")]


def test_a_frozen_state_becomes_a_parameter_not_an_algebraic_variable():
    """Truncation and residualization are different reductions -- freezing
    must remove the state from the model entirely, not relocate it."""
    from g2elin_core.components.sm import sm_dae

    frozen = sm_dae(False, (("phi_fd", FROZEN),))
    algebraic = sm_dae(False, (("phi_fd", ALGEBRAIC),))
    base = sm_dae(False)
    assert len(frozen.state_syms) == len(base.state_syms) - 1
    assert len(frozen.alg_syms) == len(base.alg_syms)  # not moved into z
    assert len(algebraic.state_syms) == len(base.state_syms) - 1
    assert len(algebraic.alg_syms) == len(base.alg_syms) + 1  # moved into z
    # The frozen state's equilibrium value is now one of the parameters.
    assert any(s.name == "phi_fd_0" for s in frozen.param_syms())


def test_freezing_and_making_algebraic_give_different_models():
    """If these two agreed, the distinction this whole design rests on would
    be decorative."""
    net = wscc9_3sm()
    a = linearized(with_models(net, network_level="quasi_stationary", sm_states={"field_flux": FROZEN}))
    b = linearized(with_models(net, network_level="quasi_stationary", sm_states={"field_flux": ALGEBRAIC}))
    assert a.A.shape == b.A.shape
    assert not np.allclose(a.A, b.A)


# --- the named levels ---------------------------------------------------------


@pytest.mark.parametrize("level", reduction.level_ids("sm"))
def test_every_machine_level_builds_and_has_the_expected_size(level):
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level=level)
    system = linearized(net)
    # 9 controller states always, plus the machine states the level keeps.
    e = reduction.element("sm")
    modes = e.modes_by_group(level)
    per_machine = sum(len(g.symbols) for g in e.groups if modes[g.id] == DYNAMIC)
    assert system.A.shape[0] == 3 * per_machine + 1  # +1 for the reference frame


def test_machine_levels_shrink_monotonically():
    sizes = [
        linearized(with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level=lv)).A.shape[0]
        for lv in ["full", "order8", "order6", "order5", "order4", "order3", "order2"]
    ]
    assert sizes == sorted(sizes, reverse=True)
    assert len(set(sizes)) == len(sizes)  # every level is genuinely different


@pytest.mark.parametrize("level", reduction.level_ids("gfm"))
def test_every_grid_forming_level_builds(level):
    net = with_models(cigre_islanded_1sm_2gfm_1gfl(), network_level="quasi_stationary",
                      sm_level="order6", gfm_level=level)
    assert linearized(net).A.shape[0] > 0


@pytest.mark.parametrize("level", reduction.level_ids("gfl"))
def test_every_grid_following_level_builds(level):
    net = with_models(cigre_islanded_1sm_2gfm_1gfl(), network_level="quasi_stationary",
                      sm_level="order6", gfl_level=level)
    assert linearized(net).A.shape[0] > 0


def test_the_smallest_converter_models_are_the_documented_sizes():
    """Droop-only is three states and PLL-only is two -- the models RMS
    tools carry. If these drift, the levels no longer mean what they say."""
    e_gfm, e_gfl = reduction.element("gfm"), reduction.element("gfl")
    gfm = e_gfm.modes_by_group("droop")
    gfl = e_gfl.modes_by_group("pll")
    assert sum(len(g.symbols) for g in e_gfm.groups if gfm[g.id] == DYNAMIC) == 3
    assert sum(len(g.symbols) for g in e_gfl.groups if gfl[g.id] == DYNAMIC) == 2


# --- the physics --------------------------------------------------------------


@pytest.mark.parametrize("preset", [wscc9_3sm, cigre_islanded_1sm_2gfm_1gfl])
def test_a_quasi_stationary_network_keeps_the_electromechanical_modes(preset):
    """The reduction's whole promise: drop the fast network dynamics, keep
    the answer to an electromechanical question."""
    net = preset()
    full = modes_in_band(linearized(net), *ELECTROMECHANICAL_HZ)
    reduced = modes_in_band(
        linearized(with_models(net, network_level="quasi_stationary")), *ELECTROMECHANICAL_HZ
    )
    assert len(reduced) == len(full)
    for (f_hz, f_d), (r_hz, r_d) in zip(full, reduced):
        assert abs(r_hz - f_hz) < 0.01, f"mode moved from {f_hz} to {r_hz} Hz"
        assert abs(r_d - f_d) < 1.0, f"damping moved from {f_d}% to {r_d}%"


def test_the_standard_stability_model_keeps_the_electromechanical_modes():
    """Quasi-stationary network + 6th-order machines is the pairing every
    stability tool uses. It has to reproduce the full model's slow modes."""
    net = wscc9_3sm()
    full = modes_in_band(linearized(net), *ELECTROMECHANICAL_HZ)
    reduced = modes_in_band(
        linearized(with_models(net, network_level="quasi_stationary", sm_level="order6")),
        *ELECTROMECHANICAL_HZ,
    )
    assert len(reduced) == len(full)
    for (f_hz, _), (r_hz, _) in zip(full, reduced):
        assert abs(r_hz - f_hz) < 0.06


def test_dropping_the_network_dynamics_removes_the_stiffness():
    """The point of the exercise is not fewer states, it is a spectrum an
    integrator can follow."""
    net = wscc9_3sm()
    fast_full = np.abs(np.linalg.eigvals(linearized(net).A)).max()
    fast_reduced = np.abs(np.linalg.eigvals(
        linearized(with_models(net, network_level="quasi_stationary", sm_level="order6")).A
    )).max()
    assert fast_reduced < fast_full / 100


# --- adequacy -----------------------------------------------------------------


def test_adequacy_reports_full_order_when_nothing_is_reduced():
    net = wscc9_3sm()
    report = check_adequacy(net, run_power_flow(net))
    assert report.verdict == "full_order"
    assert report.removed_states == []


def test_adequacy_accepts_the_standard_stability_model():
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order6")
    report = check_adequacy(net, run_power_flow(net))
    assert report.verdict in ("safe", "check")
    assert report.unmatched == 0
    assert report.n_states_reduced < report.n_states_full


def test_adequacy_rejects_a_reduction_this_network_cannot_take():
    """WSCC-9's damper windings carry real weight in its 4-5 Hz modes, so
    the classical machine model is not adequate here -- and the report has
    to say so, and name the states responsible."""
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order2")
    report = check_adequacy(net, run_power_flow(net))
    assert report.verdict == "unsafe"
    assert report.risks, "a rejected reduction should name the states that mattered"
    assert any("damper" in r.group or "flux" in r.group for r in report.risks)


def test_adequacy_works_on_a_network_with_an_open_breaker():
    """The comparison pairs states by name, and a network with elements
    switched out keeps the *full* network's numbering (network.breakers'
    labels). If the full-order twin lost those labels its state names would
    shift and every state would look removed -- a wrong answer that still
    reads like a real report."""
    from g2elin_core.network.breakers import energized_network

    net = wscc9_3sm().model_copy(update={"nodes_share_first_line_b": False, "min_node_b_pu": 1e-4})
    lines = [ln.model_copy(update={"from_closed": False}) if i == 0 else ln
             for i, ln in enumerate(net.lines)]
    net = net.model_copy(update={"lines": lines, "models": ModelOptions(
        network_level="quasi_stationary", sm_level="order6")})

    energized = energized_network(net)
    report = check_adequacy(energized, run_power_flow(energized))
    assert report.unmatched == 0
    # Only the states the reduction removes, not every state in the model.
    assert 0 < len(report.removed_states) < report.n_states_full
    assert all("?" not in _group_of(s, energized) for s in report.removed_states)


def test_a_mode_made_of_removed_states_is_not_counted_as_lost():
    """Dropping the damper windings deletes their own time constants. That
    is the reduction working; counting it as lost dynamics made every real
    reduction look unsafe, which would have made the verdict worthless."""
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order2")
    report = check_adequacy(net, run_power_flow(net))
    expected = [m for m in report.modes if m.expected_loss]
    assert expected, "order2 removes whole windings; their poles have to go with them"
    assert all(m.removed_share >= 0.5 for m in expected)
    # ...and they are excluded from the count that drives the verdict.
    assert report.unmatched == sum(1 for m in report.modes if not m.matched and not m.expected_loss)
    assert report.unmatched < sum(1 for m in report.modes if not m.matched)


def test_the_risk_list_has_no_duplicates():
    """A mode and its conjugate have identical participation factors, so
    scanning both listed every finding twice."""
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order2")
    report = check_adequacy(net, run_power_flow(net))
    keys = [(r.state, round(r.mode_hz, 6), round(r.participation, 6)) for r in report.risks]
    assert len(keys) == len(set(keys))


def test_adequacy_verdicts_get_worse_as_more_is_removed():
    """The verdict has to track the severity of the reduction, or it is not
    telling the user anything."""
    net = wscc9_3sm()
    result = run_power_flow(net)
    lost = {}
    for level in ("order6", "order4", "order3", "order2"):
        report = check_adequacy(
            with_models(net, network_level="quasi_stationary", sm_level=level), result)
        lost[level] = report.unmatched
    assert lost["order6"] == 0
    assert lost["order6"] <= lost["order4"] <= lost["order3"] <= lost["order2"]
    assert lost["order2"] > 0


def test_adequacy_narrows_its_verdict_as_the_band_narrows():
    """A reduction is adequate *for a question*. Asking only about the
    slowest modes has to be an easier test than asking about all of them."""
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order3")
    wide = check_adequacy(net, run_power_flow(net), band_hz=5.0)
    narrow = check_adequacy(net, run_power_flow(net), band_hz=1.0)
    assert narrow.unmatched <= wide.unmatched


# --- settings, validation and the catalogue -----------------------------------


def test_a_unit_can_override_the_network_default():
    net = wscc9_3sm()
    units = [d.model_copy(update={"level": "order2"}) if d.id == 3 else d for d in net.der_units]
    mixed = net.model_copy(update={"models": ModelOptions(
        network_level="quasi_stationary", sm_level="order6"), "der_units": units})
    assert mixed.unit_level(mixed.der_units[2]) == "order2"
    assert mixed.unit_level(mixed.der_units[0]) == "order6"
    uniform = with_models(net, network_level="quasi_stationary", sm_level="order6")
    assert linearized(mixed).A.shape[0] < linearized(uniform).A.shape[0]


def test_a_per_group_override_makes_the_level_read_custom():
    net = with_models(wscc9_3sm(), network_level="quasi_stationary",
                      sm_level="order6", sm_states={"avr": FROZEN})
    assert net.unit_level(net.der_units[0]) is None
    assert net.models.group_modes_for("sm")["avr"] == FROZEN


def test_an_unknown_level_or_a_forbidden_mode_is_refused():
    with pytest.raises(ValueError, match="unknown sm model level"):
        ModelOptions(sm_level="order7")
    with pytest.raises(ValueError, match="cannot be"):
        ModelOptions(sm_states={"swing": ALGEBRAIC})  # the swing equation always stays
    with pytest.raises(ValueError, match="no state group"):
        ModelOptions(sm_states={"not_a_group": ALGEBRAIC})


def test_a_mismatched_pairing_is_warned_about_but_not_blocked():
    """Keeping a machine's stator flux against a quasi-stationary network is
    a real modelling mistake, and also a legitimate thing to ask for on
    purpose. It warns; it does not refuse."""
    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order8")
    issues = validate_network(net)
    assert not [i for i in issues if i.severity == "error"]
    warnings = [i.message for i in issues if i.severity == "warning"]
    assert any("stator flux" in m for m in warnings)
    linearized(net)  # and it still builds


def test_the_standard_pairing_produces_no_model_order_warnings():
    net = with_models(
        wscc9_3sm().model_copy(update={"nodes_share_first_line_b": False, "min_node_b_pu": 1e-4}),
        network_level="quasi_stationary", sm_level="order6",
    )
    assert not [i for i in validate_network(net) if i.severity == "warning"]


@pytest.mark.parametrize(
    "kind,states,needs",
    [
        ("gfm", {"current_loop": ALGEBRAIC}, "filter"),
        ("gfm", {"voltage_loop": ALGEBRAIC}, "filter"),
        ("gfl", {"current_loop": ALGEBRAIC}, "filter"),
        ("gfl", {"outer_loop": ALGEBRAIC, "filter": ALGEBRAIC}, "dc_link"),
    ],
)
def test_an_unsolvable_control_cascade_is_refused_by_name(kind, states, needs):
    """A control loop's integrator doesn't appear in its own equation, so it
    can only go algebraic once what it regulates does. Left unchecked this
    builds and then dies as "singular algebraic Jacobian" -- true, and
    useless. The refusal has to name the group and the fix."""
    with pytest.raises(ValueError) as excinfo:
        ModelOptions(**{f"{kind}_states": states})
    assert needs in str(excinfo.value)


@pytest.mark.parametrize(
    "kind,states",
    [
        ("gfm", {"current_loop": ALGEBRAIC, "filter": ALGEBRAIC}),
        ("gfm", {"voltage_loop": ALGEBRAIC, "filter": ALGEBRAIC}),
        ("gfl", {"current_loop": ALGEBRAIC, "filter": ALGEBRAIC}),
        ("gfl", {"outer_loop": ALGEBRAIC, "filter": ALGEBRAIC, "dc_link": ALGEBRAIC}),
    ],
)
def test_the_satisfied_cascade_builds(kind, states):
    """The other half of the rule: once the requirement is met it must
    actually work, or the requirement is wrong."""
    net = cigre_islanded_1sm_2gfm_1gfl()
    net = net.model_copy(update={"models": ModelOptions(**{f"{kind}_states": states})})
    assert linearized(net).A.shape[0] > 0


def test_a_unit_override_that_breaks_a_cascade_is_refused():
    """A unit's settings and the network's are each valid alone; only their
    combination can be unsolvable, so it has to be checked where they meet."""
    net = cigre_islanded_1sm_2gfm_1gfl()
    gfm = next(d for d in net.der_units if d.unit_type.value == "gfm")
    units = [d.model_copy(update={"states": {"filter": DYNAMIC}}) if d.id == gfm.id else d
             for d in net.der_units]
    broken = net.model_copy(update={
        "models": ModelOptions(gfm_level="no_inner"), "der_units": units})
    with pytest.raises(ValueError, match="filter"):
        broken.unit_group_modes(broken.der_units[[d.id for d in broken.der_units].index(gfm.id)])


def test_every_group_with_requirements_actually_needs_them():
    """Guards the catalogue against a requirement nobody verified: every
    group that declares one must genuinely fail without it."""
    net = cigre_islanded_1sm_2gfm_1gfl()
    for kind in ("gfm", "gfl"):
        for group in reduction.element(kind).groups:
            if not group.requires:
                continue
            modes = {group.id: ALGEBRAIC, **{r: ALGEBRAIC for r in group.requires}}
            for dropped in group.requires:
                partial = {k: v for k, v in modes.items() if k != dropped}
                with pytest.raises(ValueError):
                    ModelOptions(**{f"{kind}_states": partial})
            # ...and pass with them.
            ok = net.model_copy(update={"models": ModelOptions(**{f"{kind}_states": modes})})
            linearized(ok)


def test_every_named_level_satisfies_its_own_cascade_requirements():
    for element in reduction.ELEMENTS.values():
        for level in element.levels:
            assert not element.unmet_requirements(element.modes_by_group(level)), \
                f"{element.kind}.{level} declares a cascade it does not satisfy"


def test_every_catalogued_level_resolves_to_allowed_modes():
    """The catalogue is what the UI and the schema both read; a level that
    named a mode its own group forbids would be a silent trap."""
    for kind, element in reduction.ELEMENTS.items():
        for level in element.levels:
            modes = element.modes_by_group(level)
            assert set(modes) == {g.id for g in element.groups}
            for group in element.groups:
                assert modes[group.id] in group.allowed, f"{kind}.{level}: {group.id}"


# --- the nonlinear side -------------------------------------------------------


def unit_state_rates(network: Network) -> float:
    """The largest rate of change of a *unit's* own states at t=0, ignoring
    the angles (which turn at wb by construction, equilibrium or not).

    How far the initial point is from a true equilibrium. The per-component
    operating points are built independently and are only approximately
    consistent with the coupled system, so this is never exactly zero --
    what matters is that reducing the model doesn't make it worse.
    """
    model = build_nonlinear_network(network, run_power_flow(network))
    xdot, _, _ = model.rhs(model.initial_state(), model.default_u_exo())
    unit = [
        abs(v) for v, n in zip(xdot, model.state_names)
        if "SM_" in n or "GFM_" in n or "GFL_" in n
        if not n.startswith("theta")
    ]
    return max(unit)


def test_reducing_the_model_does_not_move_the_initial_point():
    """A reduced model's operating point is the same one; only which
    variables carry it changes. So its units must start no further from
    equilibrium than the full model's do.

    This is the test that caught the bus-capacitance incompatibility below:
    a dynamic network absorbs an inconsistent initial point in fast node
    states that decay in microseconds, and a quasi-stationary one cannot, so
    anything inconsistent surfaces here instead of in a user's plot.
    """
    net = wscc9_3sm()
    full = unit_state_rates(net)
    reduced = unit_state_rates(with_models(net, network_level="quasi_stationary", sm_level="order6"))
    assert reduced <= max(full, 1e-3), f"reduced model starts at {reduced}, full at {full}"


def test_sharing_one_bus_capacitance_is_flagged_against_a_quasi_stationary_network():
    """The MATLAB-compatible bus capacitance gives every bus the *first*
    line's charging. A quasi-stationary bus equation is the power flow's own
    bus equation, so a susceptance the power flow never saw puts the two out
    of agreement -- measurably, which is why this is a warning and not a
    matter of taste."""
    own = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order6")
    assert own.nodes_share_first_line_b is False, "presets use their own capacitance now"
    shared = own.model_copy(update={"nodes_share_first_line_b": True})

    warnings = [i.message for i in validate_network(shared) if i.severity == "warning"]
    assert any("nodes_share_first_line_b" in m for m in warnings)
    assert not [i for i in validate_network(own) if "nodes_share_first_line_b" in i.message]
    assert unit_state_rates(own) < unit_state_rates(shared) / 10


def test_output_rates_match_the_states_that_still_have_derivatives():
    """output_rates() replaces a node's own dv/dt once the node is algebraic.
    Where both exist -- a dynamic network -- they must agree, or the
    quasi-stationary bus frequency is being computed by a different rule
    than the EMT one."""
    net = wscc9_3sm()
    model = build_nonlinear_network(net, run_power_flow(net))
    x = model.initial_state()
    x = x + 1e-4  # off equilibrium, so the derivatives aren't all ~0
    z, u = model.solve_algebraic(x, model.default_u_exo())
    ydot = model.output_rates(x, z, u)
    x_list, z_list, u_list = model._unpack_x(x), model._unpack_z(z), model._unpack_u(u)
    checked = 0
    for b, x_i, z_i, u_i in zip(model.blocks, x_list, z_list, u_list):
        if b.kind != "node":
            continue
        own = np.asarray(b.comp.f(x_i, z_i, u_i))  # the node's own C dv/dt
        via = ydot[b.output_off : b.output_off + 2]
        np.testing.assert_allclose(via, own, rtol=1e-6, atol=1e-6)
        checked += 1
    assert checked, "no dynamic nodes to compare against"


def test_output_rates_match_a_finite_difference_along_the_trajectory():
    """The whole point of output_rates() is that a quasi-stationary bus's
    dv/dt is derived from the algebraic system instead of differenced. That
    derivation is a page of linear algebra, so it gets checked against the
    difference it replaces."""
    net = with_models(
        wscc9_3sm().model_copy(update={"nodes_share_first_line_b": False, "min_node_b_pu": 1e-4}),
        network_level="quasi_stationary", sm_level="order6",
    )
    model = build_nonlinear_network(net, run_power_flow(net))
    u = model.default_u_exo()
    x = model.initial_state() + 1e-3  # off equilibrium, so the rates aren't ~0

    def outputs(xv):
        z, uu = model.solve_algebraic(xv, u)
        xl, zl, ul = model._unpack_x(xv), model._unpack_z(z), model._unpack_u(uu)
        return np.concatenate([b.comp.h(a, c, d) for b, a, c, d in zip(model.blocks, xl, zl, ul)])

    z, uu = model.solve_algebraic(x, u)
    analytic = model.output_rates(x, z, uu)
    xdot, _, _ = model.rhs(x, u, z, uu)
    h = 1e-7
    difference = (outputs(x + h * xdot) - outputs(x - h * xdot)) / (2 * h)

    big = np.abs(difference) > 1e-6
    assert big.sum() > 10, "not enough moving outputs to make this a real comparison"
    relative = np.abs(analytic[big] - difference[big]) / np.abs(difference[big])
    assert relative.max() < 1e-3


def test_bus_frequency_is_available_on_a_quasi_stationary_bus():
    from g2elin_core.timedomain.measurements import MeasurementSet

    net = with_models(wscc9_3sm(), network_level="quasi_stationary", sm_level="order6")
    model = build_nonlinear_network(net, run_power_flow(net))
    ms = MeasurementSet(model)
    name = next(n for n in ms.names() if n.startswith("f_{bus"))
    x = model.initial_state()
    z, u = model.solve_algebraic(x, model.default_u_exo())
    value = ms.evaluator([name])(x, z, u)[name]
    assert np.isfinite(value)
    assert abs(value - net.f_hz) < 0.5  # at rest, the bus runs at nominal


def test_a_mid_run_load_step_keeps_the_model_order():
    """An event rebuilds the element it changes. Rebuilding it at the
    default order would quietly put a dynamic load back into a
    quasi-stationary model, halfway through a run."""
    from g2elin_core.timedomain.events import NetworkEvent, apply_event

    net = with_models(
        wscc9_3sm().model_copy(update={"nodes_share_first_line_b": False, "min_node_b_pu": 1e-4}),
        network_level="quasi_stationary", sm_level="order6",
    )
    model = build_nonlinear_network(net, run_power_flow(net))
    before = model.initial_state().size
    applied = apply_event(model, NetworkEvent(kind="load_step", index=0, dp_pct=20.0, dq_pct=0.0))
    assert applied.model.initial_state().size == before
    assert all(b.comp.n_states == 0 for b in applied.model.blocks if b.kind == "load")


def test_a_phase_jump_at_a_quasi_stationary_bus_is_refused():
    """The jump is applied by rotating the bus's own voltage state. Without
    one there is nothing to rotate -- and the offset it would have written
    to lands in the next block's states, so this has to refuse rather than
    corrupt a different element."""
    from g2elin_core.timedomain.events import EventError, NetworkEvent, apply_event

    net = with_models(
        wscc9_3sm().model_copy(update={"nodes_share_first_line_b": False, "min_node_b_pu": 1e-4}),
        network_level="quasi_stationary", sm_level="order6",
    )
    model = build_nonlinear_network(net, run_power_flow(net))
    bus = next(int(b.name.split("_")[1]) for b in model.blocks if b.kind == "node")
    with pytest.raises(EventError, match="quasi-stationary"):
        apply_event(model, NetworkEvent(kind="phase_jump", bus=bus, angle_deg=10.0))


def test_a_phase_jump_still_works_with_a_dynamic_network():
    from g2elin_core.timedomain.events import NetworkEvent, apply_event

    net = wscc9_3sm()
    model = build_nonlinear_network(net, run_power_flow(net))
    bus = next(int(b.name.split("_")[1]) for b in model.blocks if b.kind == "node")
    applied = apply_event(model, NetworkEvent(kind="phase_jump", bus=bus, angle_deg=10.0))
    assert not np.allclose(applied.x0, model.initial_state())


def test_a_quasi_stationary_bus_labels_its_waveforms_as_reconstructions():
    from g2elin_core.timedomain.measurements import MeasurementSet

    net = with_models(wscc9_3sm(), network_level="quasi_stationary")
    ms = MeasurementSet(build_nonlinear_network(net, run_power_flow(net)))
    phase = next(m for m in ms.catalog() if m.name.startswith("v_a_{bus"))
    assert "phasor" in phase.label
