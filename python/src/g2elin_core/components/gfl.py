"""Grid-following converter (LCL filter + PLL), ported from
``Functions/symGFL.m`` + ``Functions/GFL_subs.m``. A GFL is never the
slack DG, so unlike SM/GFM there's only one variant.
"""

from __future__ import annotations

import cmath
from collections.abc import Mapping
from functools import lru_cache

import sympy as sp

from .base import (
    ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians,
    apply_reduction, build_dae, equilibrium_subs, mode_key,
)

wb, wff, Rf, Lf, Cf, Rt, Lt = sp.symbols("wb wff Rf Lf Cf Rt Lt")
Kpd, Kid, Kiq, KpCL, KiCL, Kffv, Cdc, Gdc, Tdc, Kipll, Kppll = sp.symbols(
    "Kpd Kid Kiq KpCL KiCL Kffv Cdc Gdc Tdc Kipll Kppll"
)
M_d, M_q, M_CLd, M_CLq, M_pll, theta_pll = sp.symbols("M_d M_q M_CLd M_CLq M_pll theta_pll")
isd, isq, igd, igq, ved, veq, vdc, idc = sp.symbols("isd isq igd igq ved veq vdc idc")
md, mq, w_pll = sp.symbols("md mq w_pll")
idc_ref, vdc_ref, q_ref, theta_g, vgd_g, vgq_g = sp.symbols("idc_ref vdc_ref q_ref theta_g vgd_g vgq_g")

_STATE_NAMES = ["i_sd", "i_sq", "i_gd", "i_gq", "v_ed", "v_eq", "v_dc", "i_dc",
                "M_d", "M_q", "M_CLd", "M_CLq", "M_pll", "theta_pll"]


@lru_cache(maxsize=16)
def gfl_dae(modes: tuple[tuple[str, str], ...] = ()) -> ComponentDAE:
    """``modes`` is a :func:`~g2elin_core.components.base.mode_key` tuple
    naming the states this converter gives up -- see
    :mod:`g2elin_core.reduction`. The default ``()`` is the full 14-state
    model.
    """
    state_vec = [isd, isq, igd, igq, ved, veq, vdc, idc, M_d, M_q, M_CLd, M_CLq, M_pll, theta_pll]
    alg_vec = [md, mq, w_pll]
    us_vec = [vdc_ref, q_ref, idc_ref]
    ug_vec = [theta_g, vgd_g, vgq_g]
    input_vec = us_vec + ug_vec

    vgd = vgd_g * sp.cos(theta_g - theta_pll) - vgq_g * sp.sin(theta_g - theta_pll)
    vgq = vgd_g * sp.sin(theta_g - theta_pll) + vgq_g * sp.cos(theta_g - theta_pll)

    disd = (wb / Lf) * (md * vdc - ved - Rf * isd + w_pll * Lf * isq)
    disq = (wb / Lf) * (mq * vdc - veq - Rf * isq - w_pll * Lf * isd)
    digd = (wb / Lt) * (ved - vgd - Rt * igd + w_pll * Lt * igq)
    digq = (wb / Lt) * (veq - vgq - Rt * igq - w_pll * Lt * igd)
    dved = (wb / Cf) * (isd - igd + w_pll * Cf * veq)
    dveq = (wb / Cf) * (isq - igq - w_pll * Cf * ved)
    dvdc = (wb / Cdc) * (idc - Gdc * vdc - md * isd - mq * isq)
    didc = (1 / Tdc) * (idc_ref - idc)

    isd_ref = Kpd * (vdc_ref - vdc) + M_d
    dM_d = Kid * (vdc_ref - vdc)

    q = -ved * igq + veq * igd
    p = ved * igd + veq * igq
    dM_q = Kiq * (q_ref - q)

    dM_CLd = KiCL * (isd_ref - isd)
    dM_CLq = KiCL * (M_q - isq)
    dM_pll = Kipll * veq
    dtheta_pll = wb * w_pll

    diffeq_vec = [disd, disq, digd, digq, dved, dveq, dvdc, didc, dM_d, dM_q, dM_CLd, dM_CLq,
                  dM_pll, dtheta_pll]

    alg = [
        md - (1 / vdc) * (KpCL * (isd_ref - isd) + Kffv * ved - wff * Lf * isq + M_CLd),
        mq - (1 / vdc) * (KpCL * (M_q - isq) + Kffv * veq + wff * Lf * isd + M_CLq),
        w_pll - M_pll - Kppll * veq - wff,
    ]

    igd_g_out = igd * sp.cos(theta_pll - theta_g) - igq * sp.sin(theta_pll - theta_g)
    igq_g_out = igd * sp.sin(theta_pll - theta_g) + igq * sp.cos(theta_pll - theta_g)
    w_pll_out = M_pll + Kppll * veq + wff  # a *different* expression from the alg var w_pll
    Vt = sp.sqrt(ved**2 + veq**2)
    output_vec = [p, q, w_pll_out, Vt, igd_g_out, igq_g_out]

    state_vec, diffeq_vec, state_names, alg_vec, alg, output_vec = apply_reduction(
        states=list(zip(state_vec, diffeq_vec, _STATE_NAMES)),
        modes=dict(modes),
        alg_vec=alg_vec,
        algeq_vec=alg,
        output_vec=output_vec,
    )

    return build_dae(
        state_vec=state_vec,
        alg_vec=alg_vec,
        input_vec=input_vec,
        diffeq_vec=diffeq_vec,
        algeq_vec=alg,
        output_vec=output_vec,
        n_us=3,
        n_ug=3,
        n_out_s=4,
        n_out_g=2,
        state_names=state_names,
        input_names=["Vdc_ref", "Q_ref", "Idc_ref"],
        output_names=["p_e", "q_e", "w_pll", "V_t"],
    )


@lru_cache(maxsize=16)
def gfl_nonlinear_funcs(modes: tuple[tuple[str, str], ...] = ()) -> NonlinearFuncs:
    return gfl_dae(modes).nonlinear_funcs()


@lru_cache(maxsize=16)
def gfl_nonlinear_jacobians(modes: tuple[tuple[str, str], ...] = ()) -> NonlinearJacobians:
    return gfl_dae(modes).nonlinear_jacobians()


class GflOperatingPoint:
    """Mirrors ``GFL_subs.m``. Eg0/Ig0/Is0/Vm0 physics are identical to
    :class:`g2elin_core.components.gfm.GfmOperatingPoint` (same LCL-filter
    converter structure) — see there for the derivation. ``M_d0``/``M_q0``/
    ``M_CLd0``/``M_CLq0``/``M_pll0`` are derived from the inner-loop
    integrators' own equilibrium conditions (dM_x = 0 at the operating
    point), since ``GFL_subs.m`` reads them from data columns this port
    doesn't have an upstream initializer for; the resulting formulas match
    GFM's structurally-identical ones exactly, as a cross-check.
    """

    def __init__(
        self,
        *,
        params: dict,
        v_terminal_pu: float,
        angle_terminal_rad: float,
        p_terminal_pu: float,
        q_terminal_pu: float,
        v_grid_pu: float,
        angle_grid_rad: float,
        theta_g_rad: float,
    ):
        self.p = params
        Rf_, Lf_ = params["Rf"], params["Lf"]

        theta0 = angle_terminal_rad  # theta_pll_0: PLL locks to the bus's own angle
        Eg0 = v_terminal_pu * cmath.exp(1j * theta0)
        Ig0 = (p_terminal_pu - 1j * q_terminal_pu) / Eg0.conjugate()
        Ic0 = Eg0 * 1j * params["Cf"]
        Is0 = Ig0 + Ic0
        Vm0 = Eg0 + Is0 * (Rf_ + 1j * Lf_)

        Eg0dq = Eg0 * cmath.exp(-1j * theta0)
        Ig0dq = Ig0 * cmath.exp(-1j * theta0)
        Is0dq = Is0 * cmath.exp(-1j * theta0)
        Vm0dq = Vm0 * cmath.exp(-1j * theta0)

        self.theta_pll0 = theta0
        self.ved0, self.veq0 = Eg0dq.real, Eg0dq.imag  # veq0 == 0 exactly (PLL-locked frame)
        self.igd0, self.igq0 = Ig0dq.real, Ig0dq.imag
        self.isd0, self.isq0 = Is0dq.real, Is0dq.imag
        vmd0, vmq0 = Vm0dq.real, Vm0dq.imag

        grid_v = v_grid_pu * cmath.exp(1j * angle_grid_rad)
        vg_local = grid_v * cmath.exp(-1j * theta0)
        self.vgd0, self.vgq0 = vg_local.real, vg_local.imag
        vg_global = grid_v * cmath.exp(-1j * theta_g_rad)
        self.vgd_g0, self.vgq_g0 = vg_global.real, vg_global.imag
        self.theta_g0 = theta_g_rad

        self.p_ref0 = p_terminal_pu
        self.q_ref0 = q_terminal_pu
        self.vdc_ref0 = 1.0
        self.vdc0 = self.vdc_ref0
        self.idc_ref0 = self.p_ref0 / self.vdc_ref0
        self.idc0 = self.idc_ref0
        self.w_pll0 = 1.0

        self.md0, self.mq0 = vmd0 / self.vdc_ref0, vmq0 / self.vdc_ref0

        Kffv_, wff_, Kppll_ = params["Kffv"], params["wff"], params["Kppll"]
        self.M_pll0 = -Kppll_ * self.veq0
        self.M_d0 = self.isd0
        self.M_q0 = self.isq0
        self.M_CLd0 = vmd0 - self.ved0 * Kffv_ + self.isq0 * Lf_ * wff_
        self.M_CLq0 = vmq0 - self.veq0 * Kffv_ - self.isd0 * Lf_ * wff_


def _gfl_operating_subs(op: GflOperatingPoint) -> dict:
    p = op.p
    return equilibrium_subs({
        wb: p["wb"], wff: p["wff"], Rf: p["Rf"], Lf: p["Lf"], Cf: p["Cf"], Rt: p["Rt"], Lt: p["Lt"],
        Kpd: p["Kpd"], Kid: p["Kid"], Kiq: p["Kiq"], KpCL: p["KpCL"], KiCL: p["KiCL"],
        Kffv: p["Kffv"], Cdc: p["Cdc"], Gdc: p["Gdc"], Tdc: p["Tdc"], Kipll: p["Kipll"], Kppll: p["Kppll"],
        isd: op.isd0, isq: op.isq0, igd: op.igd0, igq: op.igq0, ved: op.ved0, veq: op.veq0,
        vdc: op.vdc0, idc: op.idc0, M_d: op.M_d0, M_q: op.M_q0, M_CLd: op.M_CLd0, M_CLq: op.M_CLq0,
        M_pll: op.M_pll0, theta_pll: op.theta_pll0, md: op.md0, mq: op.mq0, w_pll: op.w_pll0,
        idc_ref: op.idc_ref0, vdc_ref: op.vdc_ref0, q_ref: op.q_ref0,
        theta_g: op.theta_g0, vgd_g: op.vgd_g0, vgq_g: op.vgq_g0,
    })


def linearize_gfl(op: GflOperatingPoint, modes: Mapping[str, str] | None = None) -> LinearComponent:
    return gfl_dae(mode_key(modes)).linearize(_gfl_operating_subs(op))


def gfl_nonlinear_point(op: GflOperatingPoint, modes: Mapping[str, str] | None = None) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`gfl_nonlinear_funcs`."""
    return gfl_dae(mode_key(modes)).point_from_subs(_gfl_operating_subs(op))
