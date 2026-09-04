"""Upstream infinite bus, ported from ``Functions/symIB.m`` + ``IB_subs.m``.

An IB is always the slack DG (per ``Functions/network_form.m``'s comment:
"an Infinite bus can never be a normal DG, it is always assumed to be the
slack DG"), so unlike SM/GFM/GFL there's no non-slack variant to build.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, Rup, Lup = sp.symbols("wb Rup Lup")
igd_g, igq_g, theta_up = sp.symbols("igd_g igq_g theta_up")
wup, Vup, vgd_g, vgq_g = sp.symbols("wup Vup vgd_g vgq_g")


@lru_cache(maxsize=1)
def ib_dae() -> ComponentDAE:
    state_vec = [igd_g, igq_g, theta_up]
    input_vec = [wup, Vup, vgd_g, vgq_g]  # us = [wup, Vup], ug = [vgd_g, vgq_g]

    vgd, vgq = vgd_g, vgq_g

    digd_g = (wb / Lup) * (Vup - vgd - Rup * igd_g + wup * Lup * igq_g)
    digq_g = (wb / Lup) * (-vgq - Rup * igq_g - wup * Lup * igd_g)
    dtheta_up = wb * wup
    diffeq_vec = [digd_g, digq_g, dtheta_up]

    p_up = Vup * igd_g
    q_up = -Vup * igq_g
    output_vec = [p_up, q_up, igd_g, igq_g, theta_up, wup]  # outputs_s then outputs_g

    return build_dae(
        state_vec=state_vec,
        alg_vec=[],
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=[],
        output_vec=output_vec,
        n_us=2,
        n_ug=2,
        n_out_s=2,
        n_out_g=4,
        state_names=["i_d", "i_q", "theta"],
        input_names=["omega", "V_up"],
        output_names=["p_e", "q_e"],
    )


def _ib_operating_subs(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float
) -> dict:
    p_pu, q_pu = p_mw / sn_mva, q_mvar / sn_mva
    i = (complex(p_pu, -q_pu) / v_pu) if v_pu else 0j
    return {
        wb: wb_val,
        Rup: r_pu,
        Lup: x_pu,
        wup: 1.0,
        Vup: v_pu,
        vgd_g: v_pu,
        vgq_g: 0.0,
        igd_g: i.real,
        igq_g: i.imag,
        theta_up: 0.0,
    }


def linearize_ib(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float
) -> LinearComponent:
    """Linearize the infinite bus at its power-flow operating point.

    ``r_pu``/``x_pu`` are the DG's own interconnecting-transformer impedance
    (``Y_DER`` TR_R/TR_XL, in the slot the toolbox calls ``Rup``/``Lup``).
    ``p_mw``/``q_mvar`` are the power flowing out of the IB at its own bus,
    ``v_pu`` its voltage magnitude there (angle is the local reference, 0).
    """
    subs = _ib_operating_subs(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, v_pu=v_pu, p_mw=p_mw, q_mvar=q_mvar, sn_mva=sn_mva)
    return ib_dae().linearize(subs)


@lru_cache(maxsize=1)
def ib_nonlinear_funcs() -> NonlinearFuncs:
    return ib_dae().nonlinear_funcs()


@lru_cache(maxsize=1)
def ib_nonlinear_jacobians() -> NonlinearJacobians:
    return ib_dae().nonlinear_jacobians()


def ib_nonlinear_point(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float
) -> tuple:
    subs = _ib_operating_subs(wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, v_pu=v_pu, p_mw=p_mw, q_mvar=q_mvar, sn_mva=sn_mva)
    return ib_dae().point_from_subs(subs)
