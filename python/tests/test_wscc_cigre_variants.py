"""Structural/convergence tests for the 10 additional WSCC and CIGRE-islanded
DER-mix variants added alongside the original `wscc9_3sm`/
`cigre_islanded_1sm_2gfm_1gfl` (which keep their own `tests/golden/`
coverage, untouched by this file).

No MATLAB-run reference exists for these either (same honest caveat as
`test_smib_presets.py`), so these check convergence/structural sanity, not
numeric matches.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.modal import analyze
from g2elin_core.network.presets import (
    cigre_islanded_1sm_1gfm_1gfl,
    cigre_islanded_1sm_3gfm_1gfl,
    cigre_islanded_2sm_2gfm_2gfl,
    wscc9_1gfm_2gfl,
    wscc9_1sm_1gfm_1gfl,
    wscc9_1sm_2gfl,
    wscc9_1sm_2gfm,
    wscc9_2sm_1gfl,
    wscc9_2sm_1gfm,
    wscc9_3gfm,
)
from g2elin_core.pipeline import linearize_network
from g2elin_core.powerflow import run_power_flow

# Every new variant except the two GFM-slack WSCC ones (handled separately
# below -- power flow works, modal analysis doesn't, see their presets.py
# docstrings for why).
LINEARIZABLE_VARIANTS = {
    "wscc9_2sm_1gfl": wscc9_2sm_1gfl,
    "wscc9_1sm_2gfl": wscc9_1sm_2gfl,
    "wscc9_1sm_1gfm_1gfl": wscc9_1sm_1gfm_1gfl,
    "wscc9_1sm_2gfm": wscc9_1sm_2gfm,
    "wscc9_2sm_1gfm": wscc9_2sm_1gfm,
    "cigre_islanded_1sm_1gfm_1gfl": cigre_islanded_1sm_1gfm_1gfl,
    "cigre_islanded_1sm_3gfm_1gfl": cigre_islanded_1sm_3gfm_1gfl,
    "cigre_islanded_2sm_2gfm_2gfl": cigre_islanded_2sm_2gfm_2gfl,
}
GFM_SLACK_VARIANTS = {"wscc9_1gfm_2gfl": wscc9_1gfm_2gfl, "wscc9_3gfm": wscc9_3gfm}
ALL_VARIANTS = {**LINEARIZABLE_VARIANTS, **GFM_SLACK_VARIANTS}


@pytest.mark.parametrize("name", ALL_VARIANTS)
def test_power_flow_converges(name):
    net = ALL_VARIANTS[name]()
    result = run_power_flow(net)
    assert result.converged


@pytest.mark.parametrize("name", LINEARIZABLE_VARIANTS)
def test_modal_analysis_runs(name):
    net = LINEARIZABLE_VARIANTS[name]()
    result = run_power_flow(net)
    system = linearize_network(net, result)
    modal = analyze(system.A, system.state_names)
    assert np.all(np.isfinite(modal.eigenvalues))
    assert modal.participation.shape == (system.A.shape[0], system.A.shape[0])


@pytest.mark.parametrize("name", GFM_SLACK_VARIANTS)
def test_gfm_slack_variants_modal_analysis_not_implemented(name):
    # See wscc9_1gfm_2gfl's / wscc9_3gfm's own docstrings: a GFM slack isn't
    # wired up in interconnect/network_assembly.py yet, and should fail
    # this way (a controlled NotImplementedError -> HTTP 501 at the API
    # layer), not silently produce a wrong linearization.
    net = GFM_SLACK_VARIANTS[name]()
    result = run_power_flow(net)
    with pytest.raises(NotImplementedError):
        linearize_network(net, result)
