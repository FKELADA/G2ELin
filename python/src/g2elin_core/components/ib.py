"""Upstream infinite bus, ported from ``Functions/symIB.m`` + ``IB_subs.m``.

The original toolbox made the IB the slack DG and the reference frame both
(per ``Functions/network_form.m``: "an Infinite bus can never be a normal
DG, it is always assumed to be the slack DG") -- that is ``is_slack=True``
here: the source sits at angle 0 of a frame turning with its own speed, and
the block hands that angle and speed out for the rest of the network to use.

``is_slack=False`` is the same source in a frame it does not own (see
``components/frame.py``): it reads the frame angle ``theta_g`` and speed
``wg`` like any other unit, and its source voltage enters at
``theta_up - theta_g``. The two coincide whenever the frame turns at the
source's own speed, which is the default -- what it buys is that the
infinite bus is no longer load-bearing for the model: it can be
disconnected, leaving the rest of the network islanded.
"""

from __future__ import annotations

import cmath
from functools import lru_cache

import sympy as sp

from .base import ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians, build_dae

wb, Rup, Lup = sp.symbols("wb Rup Lup")
igd_g, igq_g, theta_up = sp.symbols("igd_g igq_g theta_up")
wup, Vup, vgd_g, vgq_g = sp.symbols("wup Vup vgd_g vgq_g")
theta_g, wg = sp.symbols("theta_g wg")  # the frame this block sits in (is_slack=False)


@lru_cache(maxsize=2)
def ib_dae(is_slack: bool = True) -> ComponentDAE:
    state_vec = [igd_g, igq_g, theta_up]
    # us = [wup, Vup]; ug = [vgd_g, vgq_g], plus the frame ports when the
    # block doesn't define the frame itself.
    input_vec = [wup, Vup, vgd_g, vgq_g] if is_slack else [wup, Vup, theta_g, wg, vgd_g, vgq_g]

    # Source voltage in the common frame: on the d axis when this block *is*
    # the frame, at theta_up - theta_g otherwise. The frame's own speed drives
    # the rotational terms (they coincide when the frame turns with wup).
    if is_slack:
        ed, eq, w_frame = Vup, sp.Integer(0), wup
    else:
        delta = theta_up - theta_g
        ed, eq, w_frame = Vup * sp.cos(delta), Vup * sp.sin(delta), wg

    digd_g = (wb / Lup) * (ed - vgd_g - Rup * igd_g + w_frame * Lup * igq_g)
    digq_g = (wb / Lup) * (eq - vgq_g - Rup * igq_g - w_frame * Lup * igd_g)
    dtheta_up = wb * wup
    diffeq_vec = [digd_g, digq_g, dtheta_up]

    p_up = ed * igd_g + eq * igq_g
    q_up = eq * igd_g - ed * igq_g
    # outputs_s then outputs_g; only the frame-owning variant hands out theta/w.
    output_vec = [p_up, q_up, igd_g, igq_g] + ([theta_up, wup] if is_slack else [])

    return build_dae(
        state_vec=state_vec,
        alg_vec=[],
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=[],
        output_vec=output_vec,
        n_us=2,
        n_ug=2 if is_slack else 4,
        n_out_s=2,
        n_out_g=4 if is_slack else 2,
        state_names=["i_d", "i_q", "theta"],
        input_names=["omega", "V_up"],
        output_names=["p_e", "q_e"],
    )


def _ib_operating_subs(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float,
    delta_rad: float = 0.0, is_slack: bool = True,
) -> dict:
    """``delta_rad`` is the source's angle in the common frame (0 when this
    block *is* the frame, as it was in the original toolbox)."""
    p_pu, q_pu = p_mw / sn_mva, q_mvar / sn_mva
    rot = cmath.exp(1j * delta_rad)
    i = ((complex(p_pu, -q_pu) / v_pu) if v_pu else 0j) * rot
    v = complex(v_pu, 0.0) * rot
    subs = {
        wb: wb_val,
        Rup: r_pu,
        Lup: x_pu,
        wup: 1.0,
        Vup: v_pu,
        vgd_g: v.real,
        vgq_g: v.imag,
        igd_g: i.real,
        igq_g: i.imag,
        theta_up: delta_rad,
    }
    if not is_slack:
        subs[theta_g] = 0.0   # the frame starts at its own zero (see components/frame.py)
        subs[wg] = 1.0
    return subs


def linearize_ib(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float,
    delta_rad: float = 0.0, is_slack: bool = True,
) -> LinearComponent:
    """Linearize the infinite bus at its power-flow operating point.

    ``r_pu``/``x_pu`` are the DG's own interconnecting-transformer impedance
    (``Y_DER`` TR_R/TR_XL, in the slot the toolbox calls ``Rup``/``Lup``).
    ``p_mw``/``q_mvar`` are the power flowing out of the IB at its own bus,
    ``v_pu`` its voltage magnitude there (angle is the local reference, 0).
    """
    subs = _ib_operating_subs(
        wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, v_pu=v_pu, p_mw=p_mw, q_mvar=q_mvar, sn_mva=sn_mva,
        delta_rad=delta_rad, is_slack=is_slack,
    )
    return ib_dae(is_slack).linearize(subs)


@lru_cache(maxsize=2)
def ib_nonlinear_funcs(is_slack: bool = True) -> NonlinearFuncs:
    return ib_dae(is_slack).nonlinear_funcs()


@lru_cache(maxsize=2)
def ib_nonlinear_jacobians(is_slack: bool = True) -> NonlinearJacobians:
    return ib_dae(is_slack).nonlinear_jacobians()


def ib_nonlinear_point(
    *, wb_val: float, r_pu: float, x_pu: float, v_pu: float, p_mw: float, q_mvar: float, sn_mva: float,
    delta_rad: float = 0.0, is_slack: bool = True,
) -> tuple:
    subs = _ib_operating_subs(
        wb_val=wb_val, r_pu=r_pu, x_pu=x_pu, v_pu=v_pu, p_mw=p_mw, q_mvar=q_mvar, sn_mva=sn_mva,
        delta_rad=delta_rad, is_slack=is_slack,
    )
    return ib_dae(is_slack).point_from_subs(subs)
