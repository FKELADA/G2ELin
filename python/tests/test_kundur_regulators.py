"""The Kundur exciter and stabiliser (Fig. E12.9), beside the originals.

The figure is a thyristor exciter -- a terminal-voltage transducer, a gain
and transient gain reduction, with no exciter lag -- and a stabiliser that is
a washout and two lead-lag stages. Both are selectable per machine
(``DerUnit.exciter`` / ``DerUnit.pss``); see
``docs/sphinx/design/controller-models.md``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from g2elin_core import reduction as R
from g2elin_core.components.sm import (
    ExciterKind, GovernorKind, PssKind, linearize_sm, sm_dae, sm_nonlinear_funcs,
    sm_nonlinear_point,
)
from g2elin_core.timedomain.emt import nonlinear_sm_block
from g2elin_core.modal import reference_angle_modes, unregulated_frequency_mode
from g2elin_core.modal.analysis import analyze
from g2elin_core.network.presets import kundur_two_area, kundur_two_area_classic, wscc9_3sm
from g2elin_core.network.schema import ExciterModel, GovernorModel, PssModel
from g2elin_api.analysis import emt_response, modal_response, modal_step_response_response
from g2elin_api.main import unit_defaults
from g2elin_api.schemas import EmtRequest, StepResponseRequest, UnitDefaultsRequest
from g2elin_core.components.sm import (
    EXCITER_PARAM_NAMES, GOVERNOR_PARAM_NAMES, PSS_PARAM_NAMES,
)
from g2elin_core.network.validation import validate_network
from g2elin_core.operating_point import compute_operating_point, sm_params, unit_params

REGULATOR_NAMES = {
    k: EXCITER_PARAM_NAMES[k] + PSS_PARAM_NAMES[k] for k in ("g2elin", "kundur")
}
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow


def with_regulators(net, exciter: str, pss: str, governor: str = "g2elin"):
    return net.model_copy(update={"der_units": [
        d.model_copy(update={"exciter": ExciterModel(exciter), "pss": PssModel(pss),
                             "governor": GovernorModel(governor)})
        for d in net.der_units
    ]})


def electromechanical(net, lo: float = 0.1, hi: float = 3.0) -> list[tuple[float, float]]:
    """(frequency Hz, damping %) of the modes the machines' speeds dominate."""
    system = linearize_network(net, run_power_flow(net))
    res = analyze(system.A, system.state_names)
    speeds = [i for i, n in enumerate(system.state_names) if n.startswith("dw_r")]
    out = []
    for i, ev in enumerate(res.eigenvalues):
        f = ev.imag / (2 * math.pi)
        if lo < f < hi and np.abs(res.participation[speeds, i].real).sum() > 0.25:
            out.append((f, -100 * ev.real / abs(ev)))
    return sorted(out)


def inter_area(net) -> tuple[float, float]:
    """The inter-area mode, picked by its frequency band rather than by being
    the slowest: with the governor responding, this network also has modes
    near 0.25 Hz that the machine speeds take part in."""
    band = electromechanical(net, 0.4, 0.9)
    assert band, "no inter-area mode found"
    return min(band)


# --- the structure of each block ------------------------------------------------
def test_the_kundur_exciter_has_no_field_voltage_state():
    """A thyristor bridge has no time constant of its own, so the field
    voltage follows the regulator instantaneously. That makes E_fd an
    algebraic variable rather than a state -- the one structural difference
    from the original exciter, and the reason the machine loses two states."""
    original = sm_dae(False, (), ExciterKind.G2ELIN, PssKind.G2ELIN)
    kundur = sm_dae(False, (), ExciterKind.KUNDUR, PssKind.G2ELIN)

    assert "e_fd" in original.state_names
    assert "e_fd" not in kundur.state_names
    assert "e_tgr" in kundur.state_names
    assert len(kundur.state_names) == len(original.state_names) - 2


def test_the_kundur_stabiliser_drops_the_input_filter():
    """Fig. E12.9's stabiliser takes the speed deviation straight into its
    washout; the original filters it first."""
    original = sm_dae(False, (), ExciterKind.G2ELIN, PssKind.G2ELIN)
    kundur = sm_dae(False, (), ExciterKind.G2ELIN, PssKind.KUNDUR)

    assert "dw_1" in original.state_names
    assert "dw_1" not in kundur.state_names
    assert len(kundur.state_names) == len(original.state_names) - 1


def test_the_default_machine_is_untouched():
    """Every saved network predating the choice must linearise to exactly the
    model it always did, so the default pair has to reproduce the original
    state vector -- names and order both."""
    assert sm_dae(False).state_names == [
        "i_gd", "i_gq", "psi_d", "psi_q", "psi_fd", "psi_1d", "psi_1q", "psi_2q",
        "dw_r", "theta", "P_m", "dw_1", "v_1", "v_2", "v_pss", "e_1", "e_2", "e_fd", "e_3",
    ]


@pytest.mark.parametrize("exciter", ["g2elin", "kundur"])
@pytest.mark.parametrize("pss", ["g2elin", "kundur"])
def test_the_regulator_states_start_at_rest(exciter, pss):
    """Every exciter state the operating point hands over has to be an
    equilibrium of its own equation, or the linearisation is taken about a
    point the model would immediately leave.

    Only the exciter: the stabiliser is driven by the speed signal, which at
    this operating point is not exactly at rest either -- see the test below,
    which pins that down rather than tolerating it.
    """
    net = with_regulators(wscc9_3sm(), exciter, pss)
    op = compute_operating_point(net, run_power_flow(net))
    avr_states = {"e_1", "e_2", "e_3", "e_fd", "e_tgr"}

    for der in net.der_units:
        unit_op, modes = op.sm_ops[der.id], net.unit_modes(der)
        key = tuple(sorted(modes.items()))
        dae = sm_dae(unit_op.is_slack, key, exciter, pss)
        funcs = sm_nonlinear_funcs(unit_op.is_slack, key, exciter, pss)
        rates = dict(zip(dae.state_names, np.ravel(funcs.f(*sm_nonlinear_point(unit_op, modes)))))
        for name in avr_states & rates.keys():
            assert abs(rates[name]) < 1e-8, f"SM_{der.id}.{name} is not at rest: {rates[name]:.3e}"


def test_the_kundur_stabiliser_inherits_the_operating_points_torque_imbalance():
    """Why the stabiliser above is not checked for being at rest.

    This operating point does not balance the swing equation exactly: the
    mechanical power is set from the *terminal* power while the electrical
    torque is the air-gap one, and the two differ by the stator loss. So
    dw_r/dt sits at a few times 1e-4 rather than zero. The Kundur stabiliser
    differentiates the speed signal directly, so its washout picks that up;
    the original's input filter starts from a speed *difference* and happens
    to hide it.

    The residual is not tolerated here but identified: the washout's output
    rate must be exactly the gain times the speed residual. That pins it to
    the known cause, so a genuine error in the block cannot hide inside it,
    and it will vanish on its own if the operating point is ever made
    torque-consistent.
    """
    net = with_regulators(wscc9_3sm(), "kundur", "kundur")
    op = compute_operating_point(net, run_power_flow(net))

    for der in net.der_units:
        unit_op, modes = op.sm_ops[der.id], net.unit_modes(der)
        key = tuple(sorted(modes.items()))
        dae = sm_dae(unit_op.is_slack, key, "kundur", "kundur")
        funcs = sm_nonlinear_funcs(unit_op.is_slack, key, "kundur", "kundur")
        rates = dict(zip(dae.state_names, np.ravel(funcs.f(*sm_nonlinear_point(unit_op, modes)))))

        # v_1 starts at zero, so its rate is K_PSS * d(dw_r)/dt exactly.
        assert rates["v_1"] == pytest.approx(unit_op.p["KSTAB"] * rates["dw_r"], rel=1e-9)
        assert abs(rates["dw_r"]) < 1e-3, "the torque imbalance has grown"


# --- what the catalogue says about them ------------------------------------------
def test_the_catalogue_follows_the_machines_regulators():
    """The group *ids* stay the same whichever model is chosen -- so a saved
    level or per-group override survives swapping an exciter -- while the
    states behind them change."""
    original = R.sm_element("g2elin", "g2elin")
    kundur = R.sm_element("kundur", "kundur")

    assert [g.id for g in original.groups] == [g.id for g in kundur.groups]
    assert original.group("avr").states == ("e_1", "e_2", "e_fd", "e_3")
    assert kundur.group("avr").states == ("e_1", "e_tgr")
    assert kundur.group("pss").states == ("v_1", "v_2", "v_pss")


def test_a_state_name_resolves_to_its_group_whichever_model_produced_it():
    """Mode classification and the eigenvalue map look a state up by name,
    without knowing which unit produced it, so a name that only exists on one
    variant still has to resolve."""
    assert R.group_of_state("e_tgr_{SM_1}") == "sm.avr"
    assert R.group_of_state("e_fd_{SM_1}") == "sm.avr"
    assert R.category_of_state("e_tgr_{SM_1}") == R.CONTROL


def test_the_ui_is_told_about_both_models():
    """The web UI builds its pickers from the catalogue and hard-codes no
    model name, so each option has to arrive with the group it contributes."""
    slots = {r["id"]: r for r in R.describe("sm")["regulators"]}
    assert set(slots) == {"governor", "exciter", "pss"}
    exciter = {o["id"]: o for o in slots["exciter"]["options"]}
    assert exciter["kundur"]["group"]["states"] == ["e_1", "e_tgr"]
    assert slots["exciter"]["default"] == "g2elin"
    # A model fitted as "none" is still offered, but carries no group: there
    # are no states under it to reduce.
    for slot in ("pss", "governor"):
        options = {o["id"]: o for o in slots[slot]["options"]}
        assert "none" in options and options["none"]["group"] is None
    assert "regulators" not in R.describe("network")


def test_only_a_machine_can_carry_them():
    """A converter has no exciter, and silently ignoring one set on it would
    hide a mistake rather than report it."""
    net = wscc9_3sm()
    machine = net.der_units[0]
    payload = machine.model_dump()
    payload.update(unit_type="gfm", exciter="kundur")
    with pytest.raises(ValueError, match="only meaningful for synchronous machines"):
        type(machine).model_validate(payload)


# --- what they do to the network the figure comes from ----------------------------
def test_the_kundur_presets_carry_the_books_own_regulators_and_values():
    """Both presets use the book's regulator models, its parameter values and
    its constant-torque assumption. The damped one is Example 12.9's case
    (iv) -- high transient gain (TA = TB, so no TGR) plus a stabiliser."""
    for net in (kundur_two_area(), kundur_two_area_classic()):
        for der in net.der_units:
            assert der.exciter_model == "kundur", net.name
            assert der.governor_model == "none", net.name
            p = unit_params(net, der)
            assert (p["KA"], p["TR"]) == (200.0, 0.01)
            assert p["TA"] == p["TB"], "case (iv) pairs the stabiliser with high transient gain"

    damped = kundur_two_area()
    for der in damped.der_units:
        assert der.pss_model == "kundur"
        p = unit_params(damped, der)
        assert (p["KSTAB"], p["TW"]) == (20.0, 10.0)
        assert (p["T1"], p["T2"], p["T3"], p["T4"]) == (0.05, 0.02, 3.0, 5.4)

    for der in kundur_two_area_classic().der_units:
        assert der.pss_model == "none"


# --- fitting a governor and a stabiliser, or not ----------------------------------
def test_a_machine_can_be_built_without_a_governor_or_a_stabiliser():
    """Both are optional equipment, not just parameters turned down: leaving
    one off removes its states rather than leaving them idle."""
    full = sm_dae(False, (), ExciterKind.KUNDUR, PssKind.KUNDUR, GovernorKind.G2ELIN)
    no_gov = sm_dae(False, (), ExciterKind.KUNDUR, PssKind.KUNDUR, GovernorKind.NONE)
    bare = sm_dae(False, (), ExciterKind.KUNDUR, PssKind.NONE, GovernorKind.NONE)

    assert "P_m" in full.state_names and "P_m" not in no_gov.state_names
    assert len(no_gov.state_names) == len(full.state_names) - 1
    assert not [n for n in bare.state_names if n.startswith("v_")]
    assert len(bare.state_names) == len(no_gov.state_names) - 3


def test_a_machine_without_a_governor_holds_its_mechanical_power():
    """Not a governor with its droop turned down: there is no droop and no
    state, and the mechanical power stays where the operating point put it.

    Checked against the equivalent that *is* a reduction -- freezing the
    governor state -- which has to give the same model.
    """
    net = with_regulators(wscc9_3sm(), "kundur", "kundur", governor="none")
    frozen = with_regulators(wscc9_3sm(), "kundur", "kundur")
    frozen = frozen.model_copy(update={"der_units": [
        d.model_copy(update={"states": {"governor": "frozen"}}) for d in frozen.der_units
    ]})

    a = linearize_network(net, run_power_flow(net))
    b = linearize_network(frozen, run_power_flow(frozen))
    assert a.state_names == b.state_names
    assert np.allclose(a.A, b.A, rtol=0, atol=1e-10)


def test_removing_a_regulator_removes_its_parameters_too():
    """A machine with no stabiliser has no stabiliser parameters to set, and
    an override for one is a mistake worth reporting."""
    base = dict(sn_mva=100.0, f_hz=60.0, rt_pu=0.0, lt_pu=0.15)
    bare = set(sm_params(**base, exciter="kundur", pss="none", governor="none"))
    assert not ({"KSTAB", "TW", "T1", "T2", "T3", "T4"} & bare)
    assert not ({"mp", "TG"} & bare)
    assert {"H", "KD"} <= bare, "inertia belongs to the machine, not the governor"

    net = with_regulators(wscc9_3sm(), "kundur", "none", governor="none")
    bad = net.model_copy(update={"der_units": [
        net.der_units[0].model_copy(update={"params": {"KSTAB": 20.0}}), *net.der_units[1:]
    ]})
    assert [p for p in validate_network(bad) if "KSTAB" in str(p)]


def test_the_catalogue_drops_the_group_of_a_regulator_that_is_not_fitted():
    """The state-group controls are built from the catalogue, so a machine
    without a governor must not be offered a governor row to reduce."""
    fitted = [g.id for g in R.sm_element("kundur", "kundur", "g2elin").groups]
    bare = [g.id for g in R.sm_element("kundur", "none", "none").groups]
    assert "governor" in fitted and "pss" in fitted
    assert "governor" not in bare and "pss" not in bare
    assert "avr" in bare and "swing" in bare


def test_a_fast_exciter_with_no_tgr_and_no_pss_destabilises_the_inter_area_mode():
    """Example 12.6's result, which is the whole point of the classic preset:
    a high-gain thyristor exciter with its transient gain reduction switched
    off (TC = TB) and no stabiliser drives the slow mode unstable, while the
    local modes stay damped."""
    modes = electromechanical(kundur_two_area_classic())
    inter = min(modes)
    assert 0.5 < inter[0] < 0.75, f"inter-area mode at {inter[0]:.3f} Hz"   # published ~0.55 Hz
    assert inter[1] < 0, "the inter-area mode should be negatively damped"
    assert all(z > 0 for f, z in modes if f > 0.9), "the local modes should stay damped"


def test_transient_gain_reduction_and_a_stabiliser_together_damp_it():
    """Fig. E12.9's own configuration: the same machines with TGR and the
    stabiliser switched on turn the unstable mode into a well-damped one."""
    damped = inter_area(kundur_two_area())
    assert damped[1] > 5.0, f"only {damped[1]:.1f}% damping at {damped[0]:.3f} Hz"


def tgr_damping(tc: float) -> float:
    """Inter-area damping with no stabiliser and the transient gain set by
    ``tc``: the regulator keeps gain KA below 1/TB and KA*TA/TB above it,
    so a smaller TA is a smaller gain across the electromechanical band."""
    net = kundur_two_area_classic()
    net = net.model_copy(update={"der_units": [
        d.model_copy(update={"params": {**d.params, "TA": tc, "TB": 10.0}})
        for d in net.der_units
    ]})
    return inter_area(net)[1]


def test_transient_gain_reduction_alone_is_not_a_reliable_fix():
    """The other half of Fig. E12.9, and the reason the stabiliser is there.

    Damping is *not* monotonic in the regulator's transient gain. Detuning
    the exciter from the preset's full gain first makes the inter-area mode
    worse, bottoming out around a fiftieth-scale band gain, and only helps
    once the regulator is turned down far enough to be barely acting. Which
    side of that trough a machine sits on depends on its loading -- Kundur
    Ch. 12's K5 changing sign -- which is exactly why exciter detuning is not
    the fix and a stabiliser is.
    """
    strong, trough, weak = tgr_damping(10.0), tgr_damping(2.0), tgr_damping(0.1)
    assert trough < strong, "detuning from full gain should first make it worse"
    assert trough < weak, "turning the regulator right down should recover damping"
    assert weak > 0, "a barely-acting exciter should leave the mode damped"


def pss_contribution(tc: float, k_pss: float = 20.0) -> float:
    """How many damping points the stabiliser adds at transient gain ``tc``."""
    net = kundur_two_area_classic()
    net = net.model_copy(update={"der_units": [
        d.model_copy(update={
            "pss": PssModel.KUNDUR,
            "params": {**d.params, "TA": tc, "TB": 10.0, "KSTAB": k_pss},
        })
        for d in net.der_units
    ]})
    return inter_area(net)[1] - tgr_damping(tc)


def test_the_stabiliser_acts_through_the_exciter():
    """The two halves of Fig. E12.9 are not independent, which is the
    practical reason the figure shows them together.

    A stabiliser has no actuator of its own: it adds a signal at the
    regulator's summing junction, and whatever the regulator does not pass
    never reaches the field. So its contribution scales with the gain the
    exciter shows across the electromechanical band -- turn the transient
    gain far enough down and the stabiliser stops working, which is the
    trap in treating TGR and a PSS as two interchangeable fixes.
    """
    contributions = [pss_contribution(tc) for tc in (10.0, 5.0, 2.0, 1.0, 0.5)]

    assert contributions[0] > 1.0, "at full transient gain the stabiliser should clearly help"
    assert contributions == sorted(contributions, reverse=True), (
        f"the stabiliser's contribution should fall with the exciter's gain: {contributions}"
    )
    assert pss_contribution(0.1) < 0.1, "a barely-acting exciter leaves the stabiliser nothing to do"


# --- each model's own parameters, not equivalents ---------------------------------
def test_each_regulator_model_brings_its_own_parameter_names():
    """A model's parameters are spelled as that model spells them. Nothing is
    shared between the two exciters or the two stabilisers, and nothing is
    renamed to a common set of near-equivalents: what a user types is what
    the model's own diagram shows.
    """
    base = dict(sn_mva=100.0, f_hz=60.0, rt_pu=0.0, lt_pu=0.15)
    original = set(sm_params(**base, exciter="g2elin", pss="g2elin"))
    kundur = set(sm_params(**base, exciter="kundur", pss="kundur"))

    assert {"TR", "KA", "TA", "TB"} <= kundur                 # Fig. E12.9's exciter
    assert {"KSTAB", "TW", "T1", "T2", "T3", "T4"} <= kundur  # Fig. E12.9's stabiliser
    # The original's names are gone entirely, not carried alongside.
    assert not ({"Tr", "Ka", "Ta", "Ke", "Te", "Kfd", "Tfd"} & kundur)
    assert not ({"K_PSS", "T_LP", "T_HP", "T1n", "T1d", "T2n", "T2d"} & kundur)
    assert not ({"TR", "KA", "TA", "TB", "KSTAB", "TW"} & original)
    # Only the regulators differ; the machine itself is the same either way.
    assert original - kundur <= set(REGULATOR_NAMES["g2elin"])
    assert kundur - original <= set(REGULATOR_NAMES["kundur"])


def test_an_override_for_the_other_model_is_rejected():
    """Setting Ka on a machine with the Kundur exciter is a mistake -- that
    model has KA and no Ka -- and has to be reported rather than silently
    ignored, which is what a shared parameter set would have done."""
    net = with_regulators(wscc9_3sm(), "kundur", "kundur")
    bad = net.model_copy(update={"der_units": [
        net.der_units[0].model_copy(update={"params": {"Ka": 300.0}}), *net.der_units[1:]
    ]})
    problems = [p for p in validate_network(bad) if "Ka" in str(p)]
    assert problems, "an override the chosen exciter has no parameter for should be reported"


def test_the_defaults_endpoint_serves_the_chosen_models_parameters():
    """The editor fetches a unit's defaults by type and base values; with the
    models in that request it gets the right names to draw fields for."""
    base = dict(unit_type="sm", sn_mva=100.0, f_hz=60.0, un_kv=20.0, rt_pu=0.0, lt_pu=0.15)
    original = unit_defaults(UnitDefaultsRequest(**base)).params
    kundur = unit_defaults(UnitDefaultsRequest(**base, exciter="kundur", pss="kundur")).params

    assert "Ka" in original and "KA" not in original
    assert "KA" in kundur and "Ka" not in kundur
    assert kundur["TA"] == 1.0 and kundur["TB"] == 10.0


# --- what removing every governor does to the verdict -----------------------------
def test_a_fleet_with_no_governor_has_a_free_frequency():
    """Nothing in these presets regulates frequency: every machine has
    constant mechanical power, and there is no grid-forming converter or
    infinite bus either. A uniform speed change then meets no restoring
    torque, so the common frequency is a free integrator."""
    assert not kundur_two_area().frequency_is_regulated()
    assert not kundur_two_area_classic().frequency_is_regulated()
    assert wscc9_3sm().frequency_is_regulated(), "this one keeps its governors"

    with_gov = with_regulators(kundur_two_area(), "kundur", "kundur", governor="g2elin")
    assert with_gov.frequency_is_regulated()


def test_the_free_frequency_is_not_reported_as_an_instability():
    """It is a marginal direction of the model, like the reference angle --
    the absence of a control loop, not a dynamic fault. Left in the verdict
    it reads as a very slow instability, because a *defective* zero lands a
    few times 1e-4 from the origin numerically rather than a few times 1e-11.

    The case that must still fail is the one next to it: the classic preset
    is genuinely unstable, and excluding the free frequency must not hide
    that.
    """
    damped = modal_response(kundur_two_area())
    assert damped.stable, f"max Re = {damped.max_real_part:+.3e}"
    assert damped.max_real_part < 0

    classic = modal_response(kundur_two_area_classic())
    assert not classic.stable, "the book's own case is unstable and must stay so"
    # The inter-area mode at ~0.61 Hz, not the free frequency at ~8e-3.
    assert classic.max_real_part > 1e-2


def test_the_free_frequency_mode_is_identified_not_guessed_at():
    """It is found as the aperiodic mode nearest the origin made only of
    angle and speed states, with the reference angles set aside -- no
    tolerance, because exactly one exists once frequency is unregulated.

    A network that *does* regulate frequency has no such mode, and the same
    call has to come back empty rather than seize on the slowest real mode
    it can find.
    """
    for preset in (kundur_two_area, kundur_two_area_classic):
        net = preset()
        system = linearize_network(net, run_power_flow(net))
        res = analyze(system.A, system.state_names)
        found = unregulated_frequency_mode(res, reference_angle_modes(res))
        assert found is not None, preset.__name__
        lam = res.eigenvalues[found]
        assert abs(lam.imag) < 1e-6, "a drifting frequency does not oscillate"
        assert abs(lam) < 1e-2, f"should sit at the origin, got {lam:+.3e}"

    regulated = wscc9_3sm()
    system = linearize_network(regulated, run_power_flow(regulated))
    res = analyze(system.A, system.state_names)
    assert unregulated_frequency_mode(res, reference_angle_modes(res)) is None


# --- the nonlinear path has to agree about which models a machine carries ---------
def test_the_time_domain_model_is_built_with_the_machines_own_regulators():
    """The linear and nonlinear builders both have to read the machine's
    regulators off its operating point.

    They did not: ``linearize_sm`` was updated and ``nonlinear_sm_block`` was
    not, so a Kundur machine's parameter vector was handed to the original
    machine's equations and every time-domain run on such a network died in
    lambdify with a length mismatch. The two builders are checked against
    each other here rather than only the one that happened to be used.
    """
    net = with_regulators(wscc9_3sm(), "kundur", "kundur", governor="none")
    op = compute_operating_point(net, run_power_flow(net))
    der = net.der_units[0]
    unit_op, modes = op.sm_ops[der.id], net.unit_modes(der)

    block = nonlinear_sm_block(unit_op, modes)
    linear = linearize_sm(unit_op, modes)
    assert block.n_states == linear.A.shape[0]

    # And it evaluates: the failure was a parameter vector of the wrong length.
    x0, z0, u0 = block.x0, block.z0, block.u0
    assert np.all(np.isfinite(block.f(x0, z0, u0)))
    assert np.all(np.isfinite(block.g(x0, z0, u0)))


@pytest.mark.parametrize("preset", [kundur_two_area, kundur_two_area_classic])
def test_a_time_domain_run_produces_a_trace_stable_or_not(preset):
    """Including the classic preset, which is genuinely unstable: an unstable
    system still has a trajectory, and refusing to plot one hides exactly the
    case a user most wants to look at."""
    req = EmtRequest(
        perturb_kind="state", perturb_name="dw_r_{SM_2}", perturb_offset=1e-3,
        t_final=3.0, plot_states=["dw_r_{SM_1}", "dw_r_{SM_2}"],
    )
    r = emt_response(preset(), req)
    assert len(r.t) > 50
    for name, values in r.series.items():
        v = np.asarray(values)
        assert np.all(np.isfinite(v)), name
        assert np.ptp(v) > 1e-6, f"{name} never moved"


def test_a_power_reference_with_no_governor_says_why_it_does_nothing():
    """P_ref is the governor's setpoint. With no governor it reaches nothing
    and the step response is exactly zero -- which is the right answer, but a
    flat line with no explanation reads as a broken plot."""
    dead = modal_step_response_response(kundur_two_area(), StepResponseRequest(
        input_name="P_ref_{SM_2}", output_names=["w_r_{SM_2}"], amplitude=0.01, t_final=3.0))
    assert max(abs(v) for v in dead.y) == 0.0
    assert "no governor" in dead.note

    live = modal_step_response_response(kundur_two_area(), StepResponseRequest(
        input_name="V_ref_{SM_2}", output_names=["w_r_{SM_2}"], amplitude=0.01, t_final=3.0))
    assert max(abs(v) for v in live.y) > 0
    assert live.note == "", "an input that works needs no excuse"


def test_a_run_longer_than_three_seconds_is_allowed():
    """How long a simulation is worth watching is the user's call. The slow
    things -- a governor, an inter-area mode, a frequency drift with no
    governor at all -- take tens of seconds to say anything, so a ceiling of
    three seconds decided that for them. What stays bounded is the genuine
    cost: the number of samples returned.
    """
    assert EmtRequest().t_final == 3.0, "three seconds is the default, not the limit"

    long_run = emt_response(kundur_two_area(), EmtRequest(
        perturb_kind="state", perturb_name="dw_r_{SM_2}", perturb_offset=1e-3,
        t_final=20.0, plot_states=["dw_r_{SM_1}"]))
    assert long_run.t[-1] == pytest.approx(20.0)
    assert np.all(np.isfinite(long_run.series["dw_r_{SM_1}"]))

    for bad in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(Exception, match="positive, finite"):
            emt_response(kundur_two_area(), EmtRequest(
                perturb_kind="state", perturb_name="dw_r_{SM_2}", t_final=bad))
