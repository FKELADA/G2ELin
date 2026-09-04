"""RL-series transmission line, ported from ``Functions/symLine.m``."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, Rl, Ll = sp.symbols("wb Rl Ll")
ild_g, ilq_g = sp.symbols("ild_g ilq_g")
wg, vgdj_g, vgqj_g, vgdk_g, vgqk_g = sp.symbols("wg vgdj_g vgqj_g vgdk_g vgqk_g")


@lru_cache(maxsize=1)
def line_dae() -> ComponentDAE:
    state_vec = [ild_g, ilq_g]
    input_vec = [wg, vgdj_g, vgqj_g, vgdk_g, vgqk_g]

    dild_g = (wb / Ll) * (vgdj_g - vgdk_g - Rl * ild_g + wg * Ll * ilq_g)
    dilq_g = (wb / Ll) * (vgqj_g - vgqk_g - Rl * ilq_g - wg * Ll * ild_g)
    diffeq_vec = [dild_g, dilq_g]

    output_vec = [ild_g, ilq_g]  # outputeqVec_g only; outputeqVec_s is empty

    return build_dae(
        state_vec=state_vec,
        alg_vec=[],
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=[],
        output_vec=output_vec,
        n_us=0,
        n_ug=5,
        n_out_s=0,
        n_out_g=2,
        state_names=["i_{l_d}", "i_{l_q}"],
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=1)
def line_nonlinear_funcs() -> NonlinearFuncs:
    return line_dae().nonlinear_funcs()


@lru_cache(maxsize=1)
def line_nonlinear_jacobians() -> NonlinearJacobians:
    return line_dae().nonlinear_jacobians()


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
    return {
        wb: wb_val, Rl: r_pu, Ll: x_pu, wg: wg0, ild_g: ild_g0, ilq_g: ilq_g0,
        vgdj_g: vgdj_g0, vgqj_g: vgqj_g0, vgdk_g: vgdk_g0, vgqk_g: vgqk_g0,
    }


def linearize_line(
    *, wb_val: float, r_pu: float, x_pu: float, wg0: float, ild_g0: float, ilq_g0: float
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
    return line_dae().linearize(subs)


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
) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`line_nonlinear_funcs`, at a line's true
    operating point (unlike :func:`linearize_line`, the end voltages matter
    here).
    """
    subs = _line_operating_subs(
        wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=wg0, ild_g0=ild_g0, ilq_g0=ilq_g0,
        vgdj_g0=vgdj_g0, vgqj_g0=vgqj_g0, vgdk_g0=vgdk_g0, vgqk_g0=vgqk_g0,
    )
    return line_dae().point_from_subs(subs)
