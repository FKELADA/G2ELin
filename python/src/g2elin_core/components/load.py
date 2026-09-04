"""Constant-impedance (RL) load branch, ported from ``Functions/symLoad.m``."""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, Rc, Lc = sp.symbols("wb Rc Lc")
icd_g, icq_g = sp.symbols("icd_g icq_g")
wg, vgd_g, vgq_g = sp.symbols("wg vgd_g vgq_g")


@lru_cache(maxsize=1)
def load_dae() -> ComponentDAE:
    state_vec = [icd_g, icq_g]
    input_vec = [wg, vgd_g, vgq_g]

    dicd_g = (wb / Lc) * (vgd_g - Rc * icd_g + wg * Lc * icq_g)
    dicq_g = (wb / Lc) * (vgq_g - Rc * icq_g - wg * Lc * icd_g)
    diffeq_vec = [dicd_g, dicq_g]

    output_vec = [icd_g, icq_g]

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
        state_names=["i_{c_d}", "i_{c_q}"],
        input_names=[],
        output_names=[],
    )


@lru_cache(maxsize=1)
def load_nonlinear_funcs() -> NonlinearFuncs:
    return load_dae().nonlinear_funcs()


@lru_cache(maxsize=1)
def load_nonlinear_jacobians() -> NonlinearJacobians:
    return load_dae().nonlinear_jacobians()


def _load_operating_subs(
    *, wb_val: float, r_pu: float, x_pu: float, wg0: float, vgd_g0: float, vgq_g0: float
) -> dict:
    # Steady state: 0 = vgd_g0 - Rc*icd_g0 + wg0*Lc*icq_g0 and
    #               0 = vgq_g0 - Rc*icq_g0 - wg0*Lc*icd_g0  -- solve the 2x2 linear system.
    coeff = np.array([[-r_pu, wg0 * x_pu], [-wg0 * x_pu, -r_pu]])
    rhs = np.array([-vgd_g0, -vgq_g0])
    icd_g0, icq_g0 = np.linalg.solve(coeff, rhs)

    return {
        wb: wb_val, Rc: r_pu, Lc: x_pu, wg: wg0, vgd_g: vgd_g0, vgq_g: vgq_g0,
        icd_g: icd_g0, icq_g: icq_g0,
    }


def linearize_load(
    *, wb_val: float, r_pu: float, x_pu: float, wg0: float, vgd_g0: float, vgq_g0: float
) -> LinearComponent:
    """Linearize a constant-impedance load at an operating point.

    ``x_pu`` (the load's equivalent reactance, ``V_pu^2/S_pu * sin(acos(PF))``)
    plays the role of ``Lc``, matching the toolbox's per-unit convention.
    """
    subs = _load_operating_subs(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return load_dae().linearize(subs)


def load_nonlinear_point(
    *, wb_val: float, r_pu: float, x_pu: float, wg0: float, vgd_g0: float, vgq_g0: float
) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`load_nonlinear_funcs`."""
    subs = _load_operating_subs(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, wg0=wg0, vgd_g0=vgd_g0, vgq_g0=vgq_g0)
    return load_dae().point_from_subs(subs)
