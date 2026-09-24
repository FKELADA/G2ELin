"""Kundur's two-area system, against the book's own numbers.

Kundur, *Power System Stability and Control*, Example 12.6. This is the first
preset that is not a port of the MATLAB toolbox, and the first to need all four
of the things added for it: a unit rating (the machines are 900 MVA on a 100
MVA network), capacitor banks (the load buses do not hold up without them),
each bus's own capacitance, and published machine data.

What is checked against the book: the 400 MW tie flow, the three rotor-angle
differences, and the three electromechanical modes. What is *not* claimed is
an exact match of mode frequencies -- this tool's exciter, stabilizer and
governor are its own, not the book's, which is what
``kundur_two_area_classic`` exists to narrow (see the last test).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from g2elin_core import reduction
from g2elin_core.modal import analyze, classify_modes
from g2elin_core.network.presets import (
    kundur_machine_params, kundur_two_area, kundur_two_area_classic,
)
from g2elin_core.network.validation import validate_network
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow


def _solved(build=kundur_two_area):
    net = build()
    return net, run_power_flow(net)


def _electromechanical(net):
    """(frequency Hz, damping %) of the electromechanical modes.

    Uses the tool's own mode classifier rather than a fixed share of
    rotor-speed participation. An absolute cutoff is not a stable test: a
    strongly tuned stabiliser takes a real share of an electromechanical
    mode, so the rotor's share of it falls while the mode stays exactly what
    it was. That is the field's own definition -- rotor states participate
    *significantly*, not most -- and it is what classify_modes encodes.
    """
    sys_ = linearize_network(net, run_power_flow(net))
    res = analyze(sys_.A, sys_.state_names)
    kinds = classify_modes(res)
    return sorted(
        (ev.imag / (2 * math.pi), -ev.real / abs(ev) * 100)
        for j, ev in enumerate(res.eigenvalues)
        if ev.imag > 0 and kinds[j].category == reduction.SYNCHRONISATION
    )


# --- the machine data ----------------------------------------------------------
def test_machine_parameters_reproduce_the_published_reactances():
    """The book gives Xd/Xd'/Xd'' and the open-circuit time constants; the model
    wants mutual and leakage inductances with explicit windings. Inverting the
    conversion has to give the published numbers back."""
    p = kundur_machine_params(6.5)
    wb = 2 * math.pi * 60.0

    def par(*xs):
        return 1.0 / sum(1.0 / x for x in xs)

    xl = p["Ll"]
    assert xl + p["Lad"] == pytest.approx(1.8)                       # Xd
    assert xl + p["Laq"] == pytest.approx(1.7)                       # Xq
    assert xl + par(p["Lad"], p["Lfd"]) == pytest.approx(0.3)        # Xd'
    assert xl + par(p["Lad"], p["Lfd"], p["L1d"]) == pytest.approx(0.25)   # Xd''
    assert xl + par(p["Laq"], p["L1q"]) == pytest.approx(0.55)       # Xq'
    assert xl + par(p["Laq"], p["L1q"], p["L2q"]) == pytest.approx(0.25)   # Xq''
    assert (p["Lad"] + p["Lfd"]) / (wb * p["Rfd"]) == pytest.approx(8.0)   # Td0'
    assert (p["L1d"] + par(p["Lad"], p["Lfd"])) / (wb * p["R1d"]) == pytest.approx(0.03)
    assert (p["Laq"] + p["L1q"]) / (wb * p["R1q"]) == pytest.approx(0.4)   # Tq0'
    assert (p["L2q"] + par(p["Laq"], p["L1q"])) / (wb * p["R2q"]) == pytest.approx(0.05)


def test_the_machines_are_declared_on_their_own_900_mva_base():
    """The data above is per unit of the machine, not of the network -- which is
    what DerUnit.sn_mva is for."""
    net = kundur_two_area()
    assert {d.sn_mva for d in net.der_units} == {900.0}
    assert net.sn_mva == 100.0
    assert net.der_units[0].params["Lad"] == pytest.approx(1.6)      # as published
    # and the model sees it rebased onto the network
    from g2elin_core.operating_point import unit_params
    assert unit_params(net, net.der_units[0])["Lad"] == pytest.approx(1.6 * 100.0 / 900.0)


# --- the operating point, against the book -------------------------------------
def test_it_is_a_valid_network():
    assert [i for i in validate_network(kundur_two_area()) if i.severity == "error"] == []


def test_the_tie_carries_the_published_400_mw():
    net, result = _solved()
    lines = result.line_table()
    tie = lines[(lines.from_bus == 7) & (lines.to_bus == 8)]
    assert len(tie) == 2                                              # double circuit
    assert tie.p_from_mw.sum() == pytest.approx(400.0, abs=10.0)
    assert tie.p_from_mw.iloc[0] == pytest.approx(tie.p_from_mw.iloc[1], rel=1e-6)   # shared equally


def test_the_terminal_voltages_and_angles_match_the_book():
    """Kundur gives 1.03 at 20.2 deg, 1.01 at 10.5, 1.03 at -6.8 and 1.01 at
    -17.0 for G1..G4.

    The angles are checked absolutely, not only as differences: the preset
    holds G1 at its published 20.2 degrees (``DerUnit.angle_set_deg``), so
    the whole solution is on the book's own reference rather than 20.2
    degrees away from it.
    """
    net, result = _solved()
    t = result.bus_table().set_index("bus")
    for bus, vm, va in ((1, 1.03, 20.2), (2, 1.01, 10.5), (3, 1.03, -6.8), (4, 1.01, -17.0)):
        assert t.vm_pu[bus] == pytest.approx(vm, abs=0.005), f"G{bus} magnitude"
        assert t.va_degree[bus] == pytest.approx(va, abs=0.5), f"G{bus} angle"


def test_the_reference_angle_is_a_choice_of_origin_and_nothing_else():
    """Moving the slack unit's angle setpoint turns every angle together and
    changes nothing physical -- not a voltage, not a power flow, not an
    eigenvalue. That is what makes it safe to set it to a published value
    purely so the solved angles can be read against that source.
    """
    net = kundur_two_area()
    zeroed = net.model_copy(update={"der_units": [
        d.model_copy(update={"angle_set_deg": 0.0}) for d in net.der_units
    ]})
    a, b = run_power_flow(net), run_power_flow(zeroed)
    ta, tb = a.bus_table().set_index("bus"), b.bus_table().set_index("bus")

    assert np.allclose(ta.vm_pu, tb.vm_pu, atol=1e-10)
    shift = ta.va_degree - tb.va_degree
    assert np.allclose(shift, 20.2, atol=1e-6), "every angle should move by the same amount"

    # The slow modes -- the ones anyone reads. Compared to six significant
    # figures rather than exactly: the two runs are separate eigensolves of
    # the same model built in a different rotation, and this matrix has a
    # norm around 1e9, so they agree to a few times 1e-8 relative rather
    # than bit for bit. Anything the rotation actually changed would show up
    # far above that.
    def slow(network, result):
        sys_ = linearize_network(network, result)
        return np.array(sorted(
            (ev.imag / (2 * math.pi), -100 * ev.real / abs(ev))
            for ev in np.linalg.eigvals(sys_.A) if ev.imag > 0.1 and abs(ev) < 100
        ))
    before, after = slow(net, a), slow(zeroed, b)
    assert before.shape == after.shape
    assert np.allclose(before, after, rtol=1e-6, atol=1e-9)


def test_the_slack_machine_lands_on_its_scheduled_700_mw():
    net, result = _solved()
    assert float(result.ext_grid_table().p_mw.iloc[0]) == pytest.approx(700.0, abs=15.0)


def test_the_capacitor_banks_are_what_hold_the_load_buses_up():
    """They are not decoration: without them this operating point sags."""
    net, result = _solved()
    with_banks = result.bus_table().set_index("bus")["vm_pu"]
    bare = kundur_two_area()
    bare.shunts = []
    without = run_power_flow(bare).bus_table().set_index("bus")["vm_pu"]
    assert with_banks[7] > without[7] + 0.02
    assert with_banks[9] > without[9] + 0.02


# --- the modes the network is famous for ----------------------------------------
def test_it_has_one_inter_area_mode_and_two_local_ones():
    modes = _electromechanical(kundur_two_area())
    assert len(modes) >= 3
    inter, local1, local2 = modes[-3], modes[-2], modes[-1]
    assert 0.4 < inter[0] < 0.9                      # published ~0.55 Hz
    assert 1.0 < local1[0] < 1.5 and 1.0 < local2[0] < 1.5            # published ~1.1 Hz
    assert local1[0] < local2[0]
    assert inter[0] < local1[0]                      # the inter-area mode is the slow one


def test_the_inter_area_mode_is_the_two_areas_swinging_against_each_other():
    net = kundur_two_area()
    sys_ = linearize_network(net, run_power_flow(net))
    res = analyze(sys_.A, sys_.state_names)
    speeds = {n: i for i, n in enumerate(sys_.state_names) if n.startswith("dw_r")}
    best, best_f = None, None
    for i, ev in enumerate(res.eigenvalues):
        f = ev.imag / (2 * math.pi)
        if 0.4 < f < 0.9 and np.abs(res.participation[list(speeds.values()), i].real).sum() > 0.25:
            best, best_f = i, f
    assert best is not None, "no inter-area mode found"
    part = {n: abs(res.participation[i, best].real) for n, i in speeds.items()}
    # Every machine takes part, and area 2 (SM_3/SM_4) leads it.
    assert min(part.values()) > 0.01
    area2 = part["dw_r_{SM_3}"] + part["dw_r_{SM_4}"]
    assert area2 > 0.3 * sum(part.values())


def test_the_books_own_assumptions_give_a_negatively_damped_inter_area_mode():
    """What Example 12.6 is for: with a fast exciter, no stabilizer and constant
    mechanical torque, the inter-area mode is unstable. This tool's default
    controls damp it instead, which is the difference between the two presets.
    """
    classic = _electromechanical(kundur_two_area_classic())
    inter = min(classic)                                  # the slowest of them
    assert 0.5 < inter[0] < 0.75                          # published ~0.55 Hz
    assert inter[1] < 0                                   # negatively damped
    assert all(z > 0 for _, z in classic if _ > 0.9)      # the local modes stay damped

    with_controls = min(_electromechanical(kundur_two_area()))
    assert with_controls[1] > 0                           # the tool's PSS damps it
