"""Classifying modes by what kind of phenomenon they are.

The classification drives how the eigenvalue map is drawn, so being wrong
here is worse than being absent: a mode labelled "control" that is really
electromechanical sends someone to tune the wrong thing.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core import reduction
from g2elin_core.modal import analyze, classify_modes, state_categories
from g2elin_core.modal.classify import DOMINANCE, SYNCHRONISATION_SHARE
from g2elin_core.network.presets import cigre_islanded_1sm_2gfm_1gfl, wscc9_3sm
from g2elin_core.network.schema import ModelOptions
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow


def classified(network):
    system = linearize_network(network, run_power_flow(network))
    modal = analyze(system.A, system.state_names)
    return system, modal, classify_modes(modal)


@pytest.fixture(scope="module")
def wscc9():
    return classified(wscc9_3sm())


# --- the state lookup ---------------------------------------------------------


def test_every_state_of_a_real_model_gets_a_category(wscc9):
    system, _, _ = wscc9
    cats = state_categories(system.state_names)
    assert len(cats) == len(system.state_names)
    assert set(cats) <= set(reduction.CATEGORIES)
    assert reduction.MIXED not in cats, "a state should always resolve to a real category"


@pytest.mark.parametrize(
    "state,expected",
    [
        ("dw_r_{SM_2}", reduction.SYNCHRONISATION),
        ("theta_{SM_1}", reduction.SYNCHRONISATION),
        ("theta_{Frame}", reduction.SYNCHRONISATION),
        ("theta_{GFM_1}", reduction.SYNCHRONISATION),
        ("p_m_{GFM_1}", reduction.SYNCHRONISATION),
        ("M_pll_{GFL_1}", reduction.SYNCHRONISATION),
        ("e_fd_{SM_1}", reduction.CONTROL),
        ("v_pss_{SM_3}", reduction.CONTROL),
        ("P_m_{SM_1}", reduction.CONTROL),
        ("M_CLd_{GFM_2}", reduction.CONTROL),
        ("psi_d_{SM_1}", reduction.UNIT_ELECTRICAL),
        ("psi_fd_{SM_1}", reduction.UNIT_ELECTRICAL),
        ("i_gd_{SM_1}", reduction.UNIT_ELECTRICAL),
        ("v_dc_{GFL_1}", reduction.UNIT_ELECTRICAL),
        ("i_d_{IB_1}", reduction.UNIT_ELECTRICAL),
        ("v_{g_d}_{Nd_3}", reduction.NETWORK),
        ("i_{l_d}_{Ln_2}", reduction.NETWORK),
        ("i_{c_d}_{Ld_1}", reduction.NETWORK),
    ],
)
def test_states_land_in_the_right_category(state, expected):
    assert reduction.category_of_state(state) == expected


def test_a_state_name_with_braces_in_it_splits_on_the_last_one():
    """``v_{g_d}_{Nd_3}`` has braces in the display name too, so splitting on
    the first ``_{`` would make the block ``g_d`` and lose the element."""
    assert reduction.split_state_name("v_{g_d}_{Nd_3}") == ("v_{g_d}", "Nd_3")
    assert reduction.split_state_name("psi_d_{SM_1}") == ("psi_d", "SM_1")


# --- the classification -------------------------------------------------------


def test_every_mode_is_classified(wscc9):
    _, modal, kinds = wscc9
    assert len(kinds) == len(modal.eigenvalues)
    assert [k.mode for k in kinds] == list(range(len(modal.eigenvalues)))
    assert all(k.category in reduction.CATEGORIES for k in kinds)


def test_shares_are_a_partition(wscc9):
    _, _, kinds = wscc9
    for k in kinds:
        assert sum(k.shares.values()) == pytest.approx(1.0, abs=1e-6)


def test_the_free_reference_angles_are_their_own_kind(wscc9):
    """They are coordinates, not dynamics. Calling them synchronisation modes
    because they are made of angle states would be true and useless."""
    _, modal, kinds = wscc9
    from g2elin_core.modal import reference_angle_modes

    reference = set(reference_angle_modes(modal))
    assert reference, "this model has free reference angles to find"
    assert {k.mode for k in kinds if k.category == reduction.REFERENCE} == reference


def test_the_electromechanical_modes_are_found(wscc9):
    """WSCC-9's local modes sit near 2.1 Hz with the rotor states dominating.
    If those come back as anything but synchronisation, the classification is
    not usable."""
    _, modal, kinds = wscc9
    found = [
        k for k in kinds
        if k.category == reduction.SYNCHRONISATION
        and 1.9 <= abs(modal.eigenvalues[k.mode]) / (2 * np.pi) <= 2.4
    ]
    assert found, "the ~2.1 Hz local modes should be synchronisation modes"
    assert max(k.share for k in found) > 0.5


def test_the_exciter_modes_are_not_called_electromechanical(wscc9):
    """WSCC-9's 0.66 and 0.76 Hz modes are field-flux and AVR: only 1-2% of
    them is rotor. They sit in the same frequency band as an inter-area mode,
    so frequency alone would misclassify them -- participation doesn't."""
    _, modal, kinds = wscc9
    for k in kinds:
        hz = abs(modal.eigenvalues[k.mode]) / (2 * np.pi)
        if 0.6 <= hz <= 0.8 and modal.eigenvalues[k.mode].imag > 0:
            assert k.category == reduction.CONTROL
            assert k.shares[reduction.SYNCHRONISATION] < 0.1


def test_a_synchronisation_mode_wins_on_participation_not_majority(wscc9):
    """A machine's controllers take part in its own local mode, and there are
    more of them than there are rotor states -- so a plain "largest share"
    rule would hand rotor modes to the controllers riding on them."""
    _, _, kinds = wscc9
    borderline = [
        k for k in kinds
        if k.category == reduction.SYNCHRONISATION and k.share < DOMINANCE
    ]
    assert borderline, "WSCC-9 has a local mode the controllers dominate by share"
    for k in borderline:
        assert k.share >= SYNCHRONISATION_SHARE


def test_a_quasi_stationary_model_has_no_network_modes():
    """Making the network algebraic deletes its states, so nothing can
    participate in a network mode any more."""
    net = wscc9_3sm()
    net = net.model_copy(update={"models": ModelOptions(
        network_level="quasi_stationary", sm_level="order6")})
    _, _, kinds = classified(net)
    assert not [k for k in kinds if k.category == reduction.NETWORK]
    assert [k for k in kinds if k.category == reduction.SYNCHRONISATION]


def test_converter_units_contribute_synchronisation_modes():
    """A converter has no rotor, but its droop angle and its PLL do the same
    job. A fleet question ("which are the synchronisation modes?") has to have
    one answer across machines and converters."""
    system, modal, kinds = classified(cigre_islanded_1sm_2gfm_1gfl())
    sync = [k for k in kinds if k.category == reduction.SYNCHRONISATION]
    assert sync
    owners = set()
    for k in sync:
        top = int(np.argmax(np.abs(modal.participation[:, k.mode])))
        _, block = reduction.split_state_name(system.state_names[top])
        owners.add(block.split("_")[0])
    assert owners & {"GFM", "GFL"}, f"only {owners} contributed synchronisation modes"
