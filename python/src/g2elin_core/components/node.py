"""Network node (line-charging capacitance), ported from ``Functions/symNode.m``."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, Cl = sp.symbols("wb Cl")
vgd_g, vgq_g = sp.symbols("vgd_g vgq_g")
wg, ishd_g, ishq_g = sp.symbols("wg ishd_g ishq_g")


@lru_cache(maxsize=1)
def node_dae() -> ComponentDAE:
    state_vec = [vgd_g, vgq_g]
    input_vec = [wg, ishd_g, ishq_g]

    dvgd_g = (wb / Cl) * (ishd_g + wg * Cl * vgq_g)
    dvgq_g = (wb / Cl) * (ishq_g - wg * Cl * vgd_g)
    diffeq_vec = [dvgd_g, dvgq_g]

    output_vec = [vgd_g, vgq_g]

    return build_dae(
        state_vec=state_vec,
        alg_vec=[],
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=[],
        output_vec=output_vec,
        n_us=0,
        n_ug=3,
        n_out_s=0,
        n_out_g=2,
        state_names=["v_{g_d}", "v_{g_q}"],
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=1)
def node_nonlinear_funcs() -> NonlinearFuncs:
    return node_dae().nonlinear_funcs()


@lru_cache(maxsize=1)
def node_nonlinear_jacobians() -> NonlinearJacobians:
    return node_dae().nonlinear_jacobians()


def _node_operating_subs(*, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float) -> dict:
    ishd_g0 = -wg0 * b_pu * vgq_g0
    ishq_g0 = wg0 * b_pu * vgd_g0
    return {
        wb: wb_val, Cl: b_pu, wg: wg0, vgd_g: vgd_g0, vgq_g: vgq_g0, ishd_g: ishd_g0, ishq_g: ishq_g0,
    }


def linearize_node(
    *, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float
) -> LinearComponent:
    """Linearize a node's shunt (line-charging) capacitance at an operating point.

    ``b_pu`` (total line-charging susceptance seen at this node) plays the
    role of ``Cl`` in ``(wb/Cl)*(...)``.
    """
    subs = _node_operating_subs(wb_val=wb_val, b_pu=b_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return node_dae().linearize(subs)


def node_nonlinear_point(*, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`node_nonlinear_funcs`."""
    subs = _node_operating_subs(wb_val=wb_val, b_pu=b_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return node_dae().point_from_subs(subs)
