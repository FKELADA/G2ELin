"""Cross-validates every component's linearization against its own nonlinear
equations — the two are built from the same symbolic diffeq/algeq/output
vectors (see components/base.py) but through entirely different code paths
(symbolic Jacobian + substitution vs. lambdified nonlinear functions +
finite differences), so agreement here is real independent evidence of
correctness, not a tautology.

For each component instance, at its own linearization operating point:
1. The nonlinear differential/algebraic residuals should be ~0 (it's
   supposed to be an equilibrium).
2. A finite-difference Jacobian of the nonlinear f/g, put through the same
   Fx - Fz @ inv(Gz) @ Gx elimination formula, should reconstruct the
   symbolic A matrix from linearize_*().

This is also infrastructure for EMT (nonlinear) time-domain simulation
(feature 2.2) — the nonlinear f/g/h callables checked here are what a DAE
integrator would call at every step.

Finding this way (rather than by eyeballing eigenvalues): a real bug in
SmOperatingPoint's vdq0 rotation, which silently corrupted every non-slack
SM's operating point (masked for the slack because its bus angle is pinned
to exactly 0 by the power-flow reference). Fixed in components/sm.py; the
test that first caught it is test_sm_operating_point_self_consistency below.
"""

from __future__ import annotations

import numpy as np
import pytest

from g2elin_core.components.gfl import gfl_nonlinear_funcs, gfl_nonlinear_point, linearize_gfl
from g2elin_core.components.gfm import gfm_nonlinear_funcs, gfm_nonlinear_point, linearize_gfm
from g2elin_core.components.line import line_nonlinear_funcs, line_nonlinear_point, linearize_line
from g2elin_core.components.load import load_nonlinear_funcs, load_nonlinear_point, linearize_load
from g2elin_core.components.node import node_nonlinear_funcs, node_nonlinear_point, linearize_node
from g2elin_core.components.sm import sm_nonlinear_funcs, sm_nonlinear_point, linearize_sm
from g2elin_core.network.presets import cigre_islanded_1sm_2gfm_1gfl, wscc9_3sm
from g2elin_core.operating_point import compute_operating_point
from g2elin_core.powerflow import run_power_flow


def _fd_jacobian(fn, x0, z0, u0, p0, wrt: str, eps: float = 1e-6) -> np.ndarray:
    """Central-difference Jacobian of fn(x, z, u, p) wrt x or z."""
    vec = x0 if wrt == "x" else z0
    n_in = vec.shape[0]
    n_out = fn(x0, z0, u0, p0).shape[0]
    J = np.zeros((n_out, n_in))
    for j in range(n_in):
        d = np.zeros(n_in)
        d[j] = eps * max(1.0, abs(vec[j]))
        if wrt == "x":
            plus, minus = fn(x0 + d, z0, u0, p0), fn(x0 - d, z0, u0, p0)
        else:
            plus, minus = fn(x0, z0 + d, u0, p0), fn(x0, z0 - d, u0, p0)
        J[:, j] = (plus - minus) / (2 * d[j])
    return J


def assert_nonlinear_matches_linear(funcs, point, linear_comp, *, label: str) -> None:
    """Checks the nonlinear residuals and cross-validates the linearization.

    Three things make a naive "f(x0) == 0 everywhere, atol=1e-6" check wrong
    here, all expected rather than bugs — each was checked by hand (not
    just asserted away) before being encoded below:

    - A rotating machine's absolute angle state (``theta``/``theta_pll``)
      has derivative ``wb`` at "equilibrium", not 0 — it keeps advancing at
      synchronous speed; only its *difference* from the reference-frame
      angle is actually constant. Excluded entirely (not just loosened).
    - SM's ``i_gd``/``i_gq`` equations (``digd``/``digq``) carry a
      ``wb/Lt`` ~ 6500 coefficient, which amplifies even ordinary
      floating-point rotation/trig round-off (~1e-6, present for *every*
      unit, slack included — verified by hand: the pre-amplification
      residual ``ved0 - vgd0 - Rt*igd0 + Lt*igq0`` is ~1e-6) into an
      output that can look like a real ~0.01-0.03 error. On top of that
      noise floor, a *unit* whose own transformer impedance differs from
      DG#1's (see the next point) adds a real, larger (~10-30) contribution
      on the same two states. ``LARGE_COEFF_STATE_NAMES`` gets its own
      much looser bound (50) — comfortably above both effects, but two
      orders of magnitude below what the vdq0 rotation bug this test suite
      caught actually produced once amplified the same way (~500+).
    - Every unit's transformer impedance comes from DG#1's ``Y_TR`` row,
      not its own (a documented MATLAB-toolbox quirk preserved for
      fidelity — see components/sm.py) — a real mismatch between the
      power flow (which used each unit's own transformer) and the
      linearized model (which didn't), for every unit that isn't DG#1.
      This is what ``LARGE_COEFF_STATE_NAMES``'s loosened bound above also
      has to cover (CIGRE's preset happens to give every DG the same
      transformer impedance, so it's invisible there; WSCC's three
      different TR_XL values expose it).
    - GFM/GFL's DC-link voltage state (``v_dc``) carries a smaller
      (~1e-3) but real residual at every operating point: its equilibrium
      assumes ``idc0`` exactly balances ``md0*isd0 + mq0*isq0``
      (GFM_subs.m/GFL_subs.m's ``idc_0 = idc_ref_0 = p_ref_0/vdc_ref_0``),
      which ignores the filter resistance ``Rf``'s I^2R loss between the
      grid-side "p" the reference is set from and the converter-side power
      that expression actually represents — another originating-tool
      approximation. Covered by the default atol (5e-3), not excluded.

    The Jacobian comparison against the symbolic ``linearize_*()`` A matrix
    has none of these problems — a Jacobian is well-defined at any point,
    not just at a true equilibrium — so it stays at a strict tolerance and
    covers the excluded/loosened states too (a real bug there would still
    show up as a wrong A).
    """
    LARGE_COEFF_STATE_NAMES = {"i_gd", "i_gq"}
    x0, z0, u0, p0 = point

    f0 = funcs.f(x0, z0, u0, p0)
    rotating_idx = [i for i, n in enumerate(linear_comp.state_names) if n in ("theta", "theta_pll")]
    loose_idx = [i for i, n in enumerate(linear_comp.state_names) if n in LARGE_COEFF_STATE_NAMES]
    f0_checked = np.delete(f0, rotating_idx + loose_idx)
    assert np.allclose(f0_checked, 0, atol=5e-3), f"{label}: diffeq residual nonzero at operating point: {f0}"
    if loose_idx:
        assert np.allclose(f0[loose_idx], 0, atol=50), (
            f"{label}: i_gd/i_gq residual far beyond the expected noise+quirk range: {f0[loose_idx]}"
        )

    g0 = funcs.g(x0, z0, u0, p0)
    assert np.allclose(g0, 0, atol=1e-6), f"{label}: algeq residual nonzero at operating point: {g0}"

    if x0.shape[0] == 0:
        return
    n_x, n_z = x0.shape[0], z0.shape[0]
    Fx = _fd_jacobian(funcs.f, x0, z0, u0, p0, "x")
    Fz = _fd_jacobian(funcs.f, x0, z0, u0, p0, "z") if n_z else np.zeros((n_x, 0))
    Gx = _fd_jacobian(funcs.g, x0, z0, u0, p0, "x") if n_z else np.zeros((0, n_x))
    Gz = _fd_jacobian(funcs.g, x0, z0, u0, p0, "z") if n_z else np.zeros((0, 0))
    inv_gz = np.linalg.inv(Gz) if n_z else Gz
    A_fd = Fx - Fz @ inv_gz @ Gx

    assert np.allclose(A_fd, linear_comp.A, atol=1e-3, rtol=1e-3), (
        f"{label}: finite-difference A doesn't match the symbolic linearize_*() A\n"
        f"max abs diff = {np.max(np.abs(A_fd - linear_comp.A))}"
    )


@pytest.fixture(scope="module")
def wscc_op():
    network = wscc9_3sm()
    result = run_power_flow(network)
    return network, compute_operating_point(network, result)


@pytest.fixture(scope="module")
def cigre_op():
    network = cigre_islanded_1sm_2gfm_1gfl()
    result = run_power_flow(network)
    return network, compute_operating_point(network, result)


def test_sm_operating_point_self_consistency(wscc_op, cigre_op):
    """Every SM instance (slack and non-slack) across both presets."""
    for network, op in (wscc_op, cigre_op):
        for der_id, sm_op in op.sm_ops.items():
            funcs = sm_nonlinear_funcs(sm_op.is_slack)
            point = sm_nonlinear_point(sm_op)
            linear = linearize_sm(sm_op)
            assert_nonlinear_matches_linear(funcs, point, linear, label=f"{network.name} SM der {der_id}")


def test_gfm_operating_point_self_consistency(cigre_op):
    _, op = cigre_op
    assert op.gfm_ops, "expected at least one GFM in the CIGRE preset"
    for der_id, gfm_op in op.gfm_ops.items():
        funcs = gfm_nonlinear_funcs()
        point = gfm_nonlinear_point(gfm_op)
        linear = linearize_gfm(gfm_op)
        assert_nonlinear_matches_linear(funcs, point, linear, label=f"GFM der {der_id}")


def test_gfl_operating_point_self_consistency(cigre_op):
    _, op = cigre_op
    assert op.gfl_ops, "expected at least one GFL in the CIGRE preset"
    for der_id, gfl_op in op.gfl_ops.items():
        funcs = gfl_nonlinear_funcs()
        point = gfl_nonlinear_point(gfl_op)
        linear = linearize_gfl(gfl_op)
        assert_nonlinear_matches_linear(funcs, point, linear, label=f"GFL der {der_id}")


def test_line_node_load_self_consistency(wscc_op):
    network, op = wscc_op
    wb_val = 2 * np.pi * network.f_hz
    b_quirk = network.lines[0].b_pu

    for i, (ln, i0) in enumerate(zip(network.lines, op.line_i0)):
        vgdj, vgqj = op.node_vg[ln.from_bus]
        vgdk, vgqk = op.node_vg[ln.to_bus]
        funcs = line_nonlinear_funcs()
        point = line_nonlinear_point(
            wb_val=wb_val, r_pu=ln.r_pu, x_pu=ln.x_pu, wg0=1.0, ild_g0=i0[0], ilq_g0=i0[1],
            vgdj_g0=vgdj, vgqj_g0=vgqj, vgdk_g0=vgdk, vgqk_g0=vgqk,
        )
        linear = linearize_line(wb_val=wb_val, r_pu=ln.r_pu, x_pu=ln.x_pu, wg0=1.0, ild_g0=i0[0], ilq_g0=i0[1])
        assert_nonlinear_matches_linear(funcs, point, linear, label=f"Line {i}")

    for bus_id, (vgd, vgq) in op.node_vg.items():
        funcs = node_nonlinear_funcs()
        point = node_nonlinear_point(wb_val=wb_val, b_pu=b_quirk, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        linear = linearize_node(wb_val=wb_val, b_pu=b_quirk, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        assert_nonlinear_matches_linear(funcs, point, linear, label=f"Node {bus_id}")

    for idx, load in enumerate(network.loads):
        r_pu, x_pu = op.load_rx[idx]
        vgd, vgq = op.node_vg[load.bus]
        funcs = load_nonlinear_funcs()
        point = load_nonlinear_point(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        linear = linearize_load(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=1.0, vgd_g0=vgd, vgq_g0=vgq)
        assert_nonlinear_matches_linear(funcs, point, linear, label=f"Load {idx}")
