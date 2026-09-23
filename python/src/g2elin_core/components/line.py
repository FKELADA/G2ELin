"""RL-series transmission line, ported from ``Functions/symLine.m``.

Also serves shunt reactors and branch transformers -- both are the same RL
branch, differing only in how ``interconnect/network_assembly.py`` wires
their two ends (see that module).

With ``mode=ALGEBRAIC`` the branch's ``L di/dt`` is dropped and the current
becomes an algebraic variable of the series impedance: the quasi-stationary
branch a phasor tool uses. Nothing else about the model changes -- see
:func:`g2elin_core.components.base.apply_reduction`.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import (
    DYNAMIC, ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians,
    apply_reduction, build_dae, equilibrium_subs,
)

wb, Rl, Ll = sp.symbols("wb Rl Ll")
ild_g, ilq_g = sp.symbols("ild_g ilq_g")
wg, vgdj_g, vgqj_g, vgdk_g, vgqk_g = sp.symbols("wg vgdj_g vgqj_g vgdk_g vgqk_g")

_STATE_NAMES = ["i_{l_d}", "i_{l_q}"]


@lru_cache(maxsize=4)
def line_dae(mode: str = DYNAMIC, fixed_frequency: bool = False) -> ComponentDAE:
    """``mode`` is this branch's single state group -- see
    :mod:`g2elin_core.reduction`. ``fixed_frequency`` pins the ``w*L``
    speed-voltage term to nominal instead of following the frame, which is
    the classic phasor-network convention.
    """
    w = sp.Integer(1) if fixed_frequency else wg
    input_vec = [wg, vgdj_g, vgqj_g, vgdk_g, vgqk_g]

    dild_g = (wb / Ll) * (vgdj_g - vgdk_g - Rl * ild_g + w * Ll * ilq_g)
    dilq_g = (wb / Ll) * (vgqj_g - vgqk_g - Rl * ilq_g - w * Ll * ild_g)

    output_vec = [ild_g, ilq_g]  # outputeqVec_g only; outputeqVec_s is empty

    state_vec, diffeq_vec, state_names, alg_vec, algeq_vec, output_vec = apply_reduction(
        states=list(zip([ild_g, ilq_g], [dild_g, dilq_g], _STATE_NAMES)),
        modes={"ild_g": mode, "ilq_g": mode},
        alg_vec=[],
        algeq_vec=[],
        output_vec=output_vec,
    )

    return build_dae(
        state_vec=state_vec,
        alg_vec=alg_vec,
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=algeq_vec,
        output_vec=output_vec,
        n_us=0,
        n_ug=5,
        n_out_s=0,
        n_out_g=2,
        state_names=state_names,
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=4)
def line_nonlinear_funcs(mode: str = DYNAMIC, fixed_frequency: bool = False) -> NonlinearFuncs:
    return line_dae(mode, fixed_frequency).nonlinear_funcs()


@lru_cache(maxsize=4)
def line_nonlinear_jacobians(mode: str = DYNAMIC, fixed_frequency: bool = False) -> NonlinearJacobians:
    return line_dae(mode, fixed_frequency).nonlinear_jacobians()


def _line_operating_subs(
    *,
    wb_val: float,
    r_pu: float,
    x_pu: float,
    wg0: float,
    ild_g0: float,
    ilq_g0: float,
    vgdj_g0: float = 0.0,
    vgqj_g0: float = 0.0,
    vgdk_g0: float = 0.0,
    vgqk_g0: float = 0.0,
) -> dict:
    return equilibrium_subs({
        wb: wb_val, Rl: r_pu, Ll: x_pu, wg: wg0, ild_g: ild_g0, ilq_g: ilq_g0,
        vgdj_g: vgdj_g0, vgqj_g: vgqj_g0, vgdk_g: vgdk_g0, vgqk_g: vgqk_g0,
    })


def linearize_line(
    *, wb_val: float, r_pu: float, x_pu: float, wg0: float, ild_g0: float, ilq_g0: float,
    mode: str = DYNAMIC, fixed_frequency: bool = False,
) -> LinearComponent:
    """Linearize a line at an operating point.

    ``x_pu`` plays the role of ``Ll`` in ``(wb/Ll)*(...)``, matching the
    toolbox's convention of using per-unit reactance directly in place of an
    inductance. ``vgdj_g``/``vgqj_g``/``vgdk_g``/``vgqk_g`` (the two end
    voltages) don't appear in the linear system's coefficients (only ``wg``
    does, since it multiplies the states), so they're left at the
    :func:`_line_operating_subs` default of 0 here — harmless for
    linearization, but see :func:`line_nonlinear_point` for where the real
    values are needed (a nonlinear equilibrium check needs the true voltages).
    """
    subs = _line_operating_subs(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=wg0, ild_g0=ild_g0, ilq_g0=ilq_g0)
    return line_dae(mode, fixed_frequency).linearize(subs)


def line_nonlinear_point(
    *,
    wb_val: float,
    r_pu: float,
    x_pu: float,
    wg0: float,
    ild_g0: float,
    ilq_g0: float,
    vgdj_g0: float,
    vgqj_g0: float,
    vgdk_g0: float,
    vgqk_g0: float,
    mode: str = DYNAMIC,
    fixed_frequency: bool = False,
) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`line_nonlinear_funcs`, at a line's true
    operating point (unlike :func:`linearize_line`, the end voltages matter
    here).
    """
    subs = _line_operating_subs(
        wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=wg0, ild_g0=ild_g0, ilq_g0=ilq_g0,
        vgdj_g0=vgdj_g0, vgqj_g0=vgqj_g0, vgdk_g0=vgdk_g0, vgqk_g0=vgqk_g0,
    )
    return line_dae(mode, fixed_frequency).point_from_subs(subs)
