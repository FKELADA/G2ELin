"""Network node (line-charging capacitance), ported from ``Functions/symNode.m``.

With ``mode=ALGEBRAIC`` the bus's ``C dv/dt`` is dropped, leaving
``0 = i_sh + w*C*J*v``: the injected current equals the shunt's own
steady-state current, which is the bus row of a phasor tool's admittance
matrix. The capacitance does not disappear from the model -- it is still
the bus's shunt admittance -- it simply stops being integrated.

See :func:`g2elin_core.components.base.apply_reduction`.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import (
    DYNAMIC, ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians,
    apply_reduction, build_dae, equilibrium_subs,
)

wb, Cl = sp.symbols("wb Cl")
vgd_g, vgq_g = sp.symbols("vgd_g vgq_g")
wg, ishd_g, ishq_g = sp.symbols("wg ishd_g ishq_g")

_STATE_NAMES = ["v_{g_d}", "v_{g_q}"]


@lru_cache(maxsize=4)
def node_dae(mode: str = DYNAMIC, fixed_frequency: bool = False) -> ComponentDAE:
    """``mode`` is this bus's single state group -- see
    :mod:`g2elin_core.reduction`. ``fixed_frequency`` pins the ``w*C``
    term to nominal instead of following the frame.
    """
    w = sp.Integer(1) if fixed_frequency else wg
    input_vec = [wg, ishd_g, ishq_g]

    dvgd_g = (wb / Cl) * (ishd_g + w * Cl * vgq_g)
    dvgq_g = (wb / Cl) * (ishq_g - w * Cl * vgd_g)

    output_vec = [vgd_g, vgq_g]

    state_vec, diffeq_vec, state_names, alg_vec, algeq_vec, output_vec = apply_reduction(
        states=list(zip([vgd_g, vgq_g], [dvgd_g, dvgq_g], _STATE_NAMES)),
        modes={"vgd_g": mode, "vgq_g": mode},
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
        n_ug=3,
        n_out_s=0,
        n_out_g=2,
        state_names=state_names,
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=4)
def node_nonlinear_funcs(mode: str = DYNAMIC, fixed_frequency: bool = False) -> NonlinearFuncs:
    return node_dae(mode, fixed_frequency).nonlinear_funcs()


@lru_cache(maxsize=4)
def node_nonlinear_jacobians(mode: str = DYNAMIC, fixed_frequency: bool = False) -> NonlinearJacobians:
    return node_dae(mode, fixed_frequency).nonlinear_jacobians()


def _node_operating_subs(*, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float) -> dict:
    ishd_g0 = -wg0 * b_pu * vgq_g0
    ishq_g0 = wg0 * b_pu * vgd_g0
    return equilibrium_subs({
        wb: wb_val, Cl: b_pu, wg: wg0, vgd_g: vgd_g0, vgq_g: vgq_g0, ishd_g: ishd_g0, ishq_g: ishq_g0,
    })


def linearize_node(
    *, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float,
    mode: str = DYNAMIC, fixed_frequency: bool = False,
) -> LinearComponent:
    """Linearize a node's shunt (line-charging) capacitance at an operating point.

    ``b_pu`` (total line-charging susceptance seen at this node) plays the
    role of ``Cl`` in ``(wb/Cl)*(...)``.
    """
    subs = _node_operating_subs(wb_val=wb_val, b_pu=b_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return node_dae(mode, fixed_frequency).linearize(subs)


def node_nonlinear_point(
    *, wb_val: float, b_pu: float, wg0: float, vgd_g0: float, vgq_g0: float,
    mode: str = DYNAMIC, fixed_frequency: bool = False,
) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`node_nonlinear_funcs`."""
    subs = _node_operating_subs(wb_val=wb_val, b_pu=b_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return node_dae(mode, fixed_frequency).point_from_subs(subs)
