"""SMIB preset tests -- exercises the infinite bus (`components.ib`) as a
slack for the first time in this project (see network/presets.py's
`sm_smib` docstring for what's transcribed vs. reconstructed).

No `tests/golden/`-style MATLAB reference exists for these (the driving
script was never finished in the original tool -- same honest caveat as
CIGRE/WSCC's own tests/golden/README.md), so these are structural/
convergence checks, not numeric-match checks.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.modal import analyze
from g2elin_core.network.presets import gfl_smib, gfm_smib, sm_smib
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow
from g2elin_core.timedomain import build_nonlinear_network, simulate

PRESETS = {"sm": sm_smib, "gfm": gfm_smib, "gfl": gfl_smib}


@pytest.mark.parametrize("name", PRESETS)
def test_power_flow_converges(name):
    net = PRESETS[name]()
    result = run_power_flow(net)
    assert result.converged
    # The DER-under-test is dispatched 1.0 MW; the IB (slack) picks up
    # whatever balances losses -- both should show up in the bus table.
    table = result.bus_table().set_index("bus")
    assert table.loc[3, "p_net_gen_mw"] == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("name", PRESETS)
def test_modal_analysis_is_stable(name):
    net = PRESETS[name]()
    result = run_power_flow(net)
    system = linearize_network(net, result)
    modal = analyze(system.A, system.state_names)
    # Real part == 0 (not < 0) for the largest eigenvalue is expected, not a
    # bug: the IB's own theta_up is a pure integrator (dtheta_up = wb*wup),
    # same structural zero mode a synchronous-machine slack's absolute
    # rotor angle produces (see components/ib.py, components/sm.py).
    assert modal.eigenvalues.real.max() < 1e-6


def test_emt_simulation_runs():
    # Only sm_smib -- cheapest of the three, and the one whose Newton solve
    # was the first to actually exercise the interconnected IB (see
    # timedomain/emt.py's solve_algebraic docstring: this preset's initial
    # guess is ~0.23 pu from the true coupled equilibrium, further than
    # WSCC-9/CIGRE's ~0.08 pu -- outside hybr's basin of convergence alone,
    # which is why solve_algebraic now falls back to lm on an hybr failure).
    net = sm_smib()
    result = run_power_flow(net)
    model = build_nonlinear_network(net, result)
    sim = simulate(model, (0.0, 0.2), rtol=1e-4, atol=1e-6, first_step=1e-8)
    assert np.all(np.isfinite(sim.x))
