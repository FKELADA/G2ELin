"""The five grid-forming power-control laws of ``symGFM_types.m``.

Droop is the one this component always had. The other four -- droop behind a
filter, dVOC, VSM and matching control -- differ from it and from each other
in exactly three things: which states the outer loop carries, how it forms
the frequency deviation, and how it forms the voltage reference. Everything
downstream is shared, which is what these tests lean on: a law is checked by
what it does differently, not by re-checking the cascade underneath it.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core import reduction as R
from g2elin_core.components.gfm import (
    GfmControlKind, OUTER_PARAM_NAMES, gfm_dae, gfm_nonlinear_funcs, gfm_nonlinear_point,
    linearize_gfm,
)
from g2elin_core.modal import analyze
from g2elin_core.network.presets import gfm_smib
from g2elin_core.network.schema import GfmController
from g2elin_core.network.validation import validate_network
from g2elin_core.operating_point import compute_operating_point, gfm_params, unit_params
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain.emt import nonlinear_gfm_block

LAWS = [c.value for c in GfmController]
PARAM_BASE = dict(sn_mva=100.0, f_hz=50.0, un_kv=20.0, rt_pu=0.0, lt_pu=0.05)


def with_law(net, law: str):
    return net.model_copy(update={"der_units": [
        d.model_copy(update={"controller": GfmController(law)}) if d.unit_type.value == "gfm" else d
        for d in net.der_units
    ]})


# --- the state vectors, against the file they come from ---------------------------
@pytest.mark.parametrize("law, outer", [
    ("droop", ["p_m", "theta", "q_m"]),
    ("droop_filtered", ["p_m", "dw", "theta", "q_m"]),
    ("dvoc", ["p_m", "theta", "q_m", "v_ed_ref"]),
    ("vsm", ["dw", "theta", "Phi"]),
    ("matching", ["v_dc_m", "theta"]),
])
def test_each_law_has_symgfm_types_own_state_vector(law, outer):
    """Order included, so a state vector here can be read straight against
    ``symGFM_types.m``'s ``stateVec`` for the same case. The plant states
    come before and the two cascaded loops after, identically for all five."""
    names = gfm_dae((), law).state_names
    assert names[:8] == ["i_sd", "i_sq", "i_gd", "i_gq", "v_ed", "v_eq", "v_dc", "i_dc"]
    assert names[-4:] == ["M_VLd", "M_VLq", "M_CLd", "M_CLq"]
    assert names[8:-4] == outer


def test_the_default_converter_is_untouched():
    """Droop is the default, and has to stay the 15-state model every saved
    network and every golden test was built against."""
    assert gfm_dae().state_names == gfm_dae((), GfmControlKind.DROOP).state_names
    assert len(gfm_dae().state_names) == 15


def test_a_virtual_machine_has_no_power_filters():
    """VSM measures power directly -- its inertia is what smooths the
    response -- so the filters droop needs are not merely retuned but
    absent. Matching measures no power at all."""
    for law in ("droop", "droop_filtered", "dvoc"):
        assert {"p_m", "q_m"} <= set(gfm_dae((), law).state_names), law
    for law in ("vsm", "matching"):
        assert not ({"p_m", "q_m"} & set(gfm_dae((), law).state_names)), law


# --- each law's own parameters ----------------------------------------------------
def test_each_law_brings_its_own_parameters():
    """As with a machine's regulators, nothing is shared beyond what is
    genuinely the same quantity: `wf` is the power-measurement filter and
    appears wherever there is one, while droop's `mp`/`nq`, dVOC's
    `eta`/`alfa` and VSM's `J`/`Dp`/`K`/`Dq` belong to their own law."""
    keys = {law: set(gfm_params(**PARAM_BASE, controller=law)) for law in LAWS}

    assert {"mp", "nq"} <= keys["droop"]
    assert "wc" in keys["droop_filtered"] and "wc" not in keys["droop"]
    assert {"eta", "alfa"} <= keys["dvoc"] and not ({"mp", "nq"} & keys["dvoc"])
    assert {"J", "Dp", "K", "Dq"} <= keys["vsm"] and "wf" not in keys["vsm"]
    assert "K_theta" in keys["matching"] and not ({"mp", "nq"} & keys["matching"])

    # Only the outer loop differs; the cascade underneath is the same set.
    shared = set.intersection(*keys.values())
    for law in LAWS:
        assert keys[law] - shared <= set(OUTER_PARAM_NAMES[law]), law


def test_the_laws_are_tuned_to_be_comparable():
    """script_generic.m writes every law's gains in terms of the droop
    tuning, so swapping one for another shows what the *law* changes rather
    than what a different tuning changes. These identities are that intent."""
    droop = gfm_params(**PARAM_BASE, controller="droop")
    mp, nq, wf = droop["mp"], droop["nq"], droop["wf"]

    assert gfm_params(**PARAM_BASE, controller="dvoc")["eta"] == pytest.approx(mp)
    assert gfm_params(**PARAM_BASE, controller="dvoc")["alfa"] == pytest.approx(1 / (2 * nq))
    vsm = gfm_params(**PARAM_BASE, controller="vsm")
    assert vsm["J"] == pytest.approx(1 / (mp * wf))
    assert vsm["Dp"] == pytest.approx(1 / mp)
    assert vsm["K"] == pytest.approx(1 / (nq * wf))
    assert vsm["Dq"] == pytest.approx(1 / nq)
    matching = gfm_params(**PARAM_BASE, controller="matching")
    assert matching["K_theta"] == pytest.approx(mp * droop["Kpdc"])


def test_an_override_for_another_law_is_rejected():
    """Setting mp on a VSM is a mistake -- that law has J and Dp -- and has
    to be reported rather than silently ignored."""
    bad = with_law(gfm_smib(), "vsm")
    bad = bad.model_copy(update={"der_units": [
        d.model_copy(update={"params": {"mp": 0.01}}) if d.unit_type.value == "gfm" else d
        for d in bad.der_units
    ]})
    assert [p for p in validate_network(bad) if "mp" in str(p)]


# --- every law starts where the power flow left it --------------------------------
@pytest.mark.parametrize("law", LAWS)
def test_the_outer_states_start_at_rest(law):
    """Each law's own states must be an equilibrium of its own equations, or
    the linearisation is taken about a point the model would leave at once.

    The plant states are excluded: the DC link carries a small residual at
    this operating point that predates all of this and is identical for
    every law, droop included.
    """
    net = with_law(gfm_smib(), law)
    op = compute_operating_point(net, run_power_flow(net))
    der = next(d for d in net.der_units if d.unit_type.value == "gfm")
    unit_op, modes = op.gfm_ops[der.id], net.unit_modes(der)
    key = tuple(sorted(modes.items()))

    dae = gfm_dae(key, law)
    rates = dict(zip(dae.state_names, np.ravel(
        gfm_nonlinear_funcs(key, law).f(*gfm_nonlinear_point(unit_op, modes)))))
    outer = set(dae.state_names[8:-4]) - {"theta"}
    assert outer, law
    for name in outer:
        assert abs(rates[name]) < 1e-8, f"{law}: {name} is not at rest ({rates[name]:.3e})"


@pytest.mark.parametrize("law", LAWS)
def test_the_whole_network_linearises_and_is_stable(law):
    """All five laws are tuned to the same equivalent inertia, so on a
    network droop holds up, they all should."""
    net = with_law(gfm_smib(), law)
    system = linearize_network(net, run_power_flow(net))
    res = analyze(system.A, system.state_names)
    assert np.all(np.isfinite(system.A))
    physical = [ev for ev in res.eigenvalues if abs(ev) > 1e-6]
    assert max(ev.real for ev in physical) < 1e-6, f"{law} is not stable"


@pytest.mark.parametrize("law", LAWS)
def test_the_time_domain_model_is_built_with_the_units_own_law(law):
    """The mistake the machine's regulators made once: the linear builder
    read the model off the operating point and the nonlinear one did not, so
    a converter's parameters met another law's equations. Checked here for
    every law rather than only the default."""
    net = with_law(gfm_smib(), law)
    op = compute_operating_point(net, run_power_flow(net))
    der = next(d for d in net.der_units if d.unit_type.value == "gfm")
    unit_op, modes = op.gfm_ops[der.id], net.unit_modes(der)

    block = nonlinear_gfm_block(unit_op, modes)
    assert block.n_states == linearize_gfm(unit_op, modes).A.shape[0]
    assert np.all(np.isfinite(block.f(block.x0, block.z0, block.u0)))
    assert np.all(np.isfinite(block.g(block.x0, block.z0, block.u0)))


# --- what the catalogue and the UI are told ---------------------------------------
def test_the_catalogue_follows_the_law():
    """The shared cascade keeps its group ids whichever law is running, so a
    saved level or per-group override survives a swap; only the outer groups
    change."""
    shared = {"trafo_current", "filter", "dc_link", "current_loop", "voltage_loop", "angle"}
    for law in LAWS:
        ids = {g.id for g in R.gfm_element(law).groups}
        assert shared <= ids, law
    assert "power_filter" in {g.id for g in R.gfm_element("droop").groups}
    assert "power_filter" not in {g.id for g in R.gfm_element("vsm").groups}
    assert {"frequency", "flux"} <= {g.id for g in R.gfm_element("vsm").groups}
    assert "dc_measurement" in {g.id for g in R.gfm_element("matching").groups}


def test_a_state_name_resolves_whichever_law_produced_it():
    """Mode classification looks a state up by name without knowing which
    unit produced it, so a name only one law has still has to resolve."""
    # `dw` is owned by two laws -- the filtered droop's filter state and the
    # VSM's swing state -- and both file it under one group id, so this stays
    # single-valued rather than answering with whichever registered last.
    for name, group in (("Phi", "flux"), ("v_dc_m", "dc_measurement"),
                        ("v_ed_ref", "voltage_dynamics"), ("dw", "frequency")):
        assert R.group_of_state(f"{name}_{{GFM_1}}") == f"gfm.{group}", name
        assert R.category_of_state(f"{name}_{{GFM_1}}") == R.SYNCHRONISATION


def test_the_ui_is_told_about_every_law():
    """A law brings a *set* of groups where a regulator brings one, so the
    catalogue hands the UI the whole set to swap in."""
    slots = {r["id"]: r for r in R.describe("gfm")["regulators"]}
    assert set(slots) == {"controller"}
    options = {o["id"]: o for o in slots["controller"]["options"]}
    assert set(options) == set(LAWS)
    assert [g["id"] for g in options["vsm"]["groups"]] == ["frequency", "flux"]
    assert [g["id"] for g in options["droop"]["groups"]] == ["power_filter"]
    assert slots["controller"]["default"] == "droop"


def test_a_unit_reports_the_parameters_of_the_law_it_runs():
    """unit_params is what the editor and the sensitivity scan read."""
    vsm = with_law(gfm_smib(), "vsm")
    der = next(d for d in vsm.der_units if d.unit_type.value == "gfm")
    p = unit_params(vsm, der)
    assert {"J", "Dp", "K", "Dq"} <= set(p)
    assert "mp" not in p and "nq" not in p
