"""Grid-forming converter (Droop control only), ported from
``Functions/symGFM_types.m``'s ``'Droop'`` case + ``Functions/GFM_subs.m``.

Only the non-slack variant and the ``Droop`` outer-loop are implemented —
CIGRE's islanded preset uses ``GFM_P_control = 'Droop'`` and its GFM units
are never the slack. ``Droop+filter``/``Droop+PLL``/``dVOC``/``VSM``/
``Matching`` and the slack (``i_DG == 1``) branch follow the exact same
build_dae/linearize pattern (see ``symGFM_types.m``'s other ``switch``
cases) and would extend ``gfm_dae``/``GfmOperatingPoint`` the same way
``sm.py`` handles ``is_slack`` — not done here to keep scope to what
CIGRE's preset actually needs.
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

wb, wff, Rf, Lf, Cf, Rt, Lt, mp, nq, wf, KpVL, KiVL, Kffi, KpCL, KiCL, Kffv, Cdc, Gdc, Kpdc, Tdc = (
    sp.symbols("wb wff Rf Lf Cf Rt Lt mp nq wf KpVL KiVL Kffi KpCL KiCL Kffv Cdc Gdc Kpdc Tdc")
)
isd, isq, igd, igq, ved, veq, vdc, idc = sp.symbols("isd isq igd igq ved veq vdc idc")
pm, theta, qm, M_VLd, M_VLq, M_CLd, M_CLq = sp.symbols("pm theta qm M_VLd M_VLq M_CLd M_CLq")
md, mq, w = sp.symbols("md mq w")
p_ref, q_ref, ve_ref, w_ref, vdc_ref, theta_g, vgd_g, vgq_g = sp.symbols(
    "p_ref q_ref ve_ref w_ref vdc_ref theta_g vgd_g vgq_g"
)

_STATE_NAMES = ["i_sd", "i_sq", "i_gd", "i_gq", "v_ed", "v_eq", "v_dc", "i_dc",
                "p_m", "theta", "q_m", "M_VLd", "M_VLq", "M_CLd", "M_CLq"]


@lru_cache(maxsize=16)
def gfm_dae(modes: tuple[tuple[str, str], ...] = ()) -> ComponentDAE:
    """``modes`` is a :func:`~g2elin_core.components.base.mode_key` tuple
    naming the states this converter gives up -- see
    :mod:`g2elin_core.reduction`. The default ``()`` is the full 15-state
    model.
    """
    state_vec = [isd, isq, igd, igq, ved, veq, vdc, idc, pm, theta, qm, M_VLd, M_VLq, M_CLd, M_CLq]
    alg_vec = [md, mq, w]
    us_vec = [p_ref, q_ref, ve_ref, w_ref, vdc_ref]
    ug_vec = [theta_g, vgd_g, vgq_g]
    input_vec = us_vec + ug_vec

    vgd = vgd_g * sp.cos(theta_g - theta) - vgq_g * sp.sin(theta_g - theta)
    vgq = vgd_g * sp.sin(theta_g - theta) + vgq_g * sp.cos(theta_g - theta)

    disd = (wb / Lf) * (md * vdc - ved - Rf * isd + w * Lf * isq)
    disq = (wb / Lf) * (mq * vdc - veq - Rf * isq - w * Lf * isd)
    digd = (wb / Lt) * (ved - vgd - Rt * igd + w * Lt * igq)
    digq = (wb / Lt) * (veq - vgq - Rt * igq - w * Lt * igd)
    dved = (wb / Cf) * (isd - igd + w * Cf * veq)
    dveq = (wb / Cf) * (isq - igq - w * Cf * ved)
    dvdc = (wb / Cdc) * (idc - Gdc * vdc - md * isd - mq * isq)

    idc_ref = (p_ref / vdc_ref) + Kpdc * (vdc_ref - vdc)
    didc = (idc_ref - idc) / Tdc

    q = -ved * igq + veq * igd
    p = ved * igd + veq * igq
    Dw = mp * (p_ref - pm)
    dpm = wf * (p - pm)
    dqm = wf * (q - qm)
    ved_ref = ve_ref + (q_ref - qm) * nq
    veq_ref = 0
    dtheta = wb * w

    dM_VLd = KiVL * (ved_ref - ved)
    dM_VLq = KiVL * (veq_ref - veq)
    isd_ref = KpVL * (ved_ref - ved) + M_VLd + Kffi * igd - wff * Cf * veq
    isq_ref = KpVL * (veq_ref - veq) + M_VLq + Kffi * igq + wff * Cf * ved
    dM_CLd = KiCL * (isd_ref - isd)
    dM_CLq = KiCL * (isq_ref - isq)

    diffeq_vec = [disd, disq, digd, digq, dved, dveq, dvdc, didc, dpm, dtheta, dqm,
                  dM_VLd, dM_VLq, dM_CLd, dM_CLq]

    alg = [
        md - (1 / vdc) * (KpCL * (isd_ref - isd) + Kffv * ved - wff * Lf * isq + M_CLd),
        mq - (1 / vdc) * (KpCL * (isq_ref - isq) + Kffv * veq + wff * Lf * isd + M_CLq),
        w - Dw - w_ref,
    ]

    igd_g_out = igd * sp.cos(theta - theta_g) - igq * sp.sin(theta - theta_g)
    igq_g_out = igd * sp.sin(theta - theta_g) + igq * sp.cos(theta - theta_g)
    Vt = sp.sqrt(ved**2 + veq**2)
    output_vec = [p, q, w, Vt, igd_g_out, igq_g_out]

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
        n_us=5,
        n_ug=3,
        n_out_s=4,
        n_out_g=2,
        state_names=state_names,
        input_names=["P_ref", "Q_ref", "V_ref", "w_ref", "Vdc_ref"],
        output_names=["p_e", "q_e", "w", "V_t"],
    )


@lru_cache(maxsize=16)
def gfm_nonlinear_funcs(modes: tuple[tuple[str, str], ...] = ()) -> NonlinearFuncs:
    return gfm_dae(modes).nonlinear_funcs()


@lru_cache(maxsize=16)
def gfm_nonlinear_jacobians(modes: tuple[tuple[str, str], ...] = ()) -> NonlinearJacobians:
    return gfm_dae(modes).nonlinear_jacobians()


class GfmOperatingPoint:
    """Mirrors ``GFM_subs.m`` + the "GFM Initializations" block of ``script_generic.m``."""

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

        # Eg0 (= Vref*exp(i*Theta_Eg0)) rotated by its own angle is real by
        # construction, i.e. Eg0_d = Vref, Eg0_q = 0 always.
        theta0 = angle_terminal_rad
        Eg0 = v_terminal_pu * cmath.exp(1j * theta0)
        Ig0 = (p_terminal_pu - 1j * q_terminal_pu) / Eg0.conjugate()
        Ic0 = Eg0 * 1j * params["Cf"]
        Is0 = Ig0 + Ic0
        Vm0 = Eg0 + Is0 * (Rf_ + 1j * Lf_)

        Eg0dq = Eg0 * cmath.exp(-1j * theta0)
        Ig0dq = Ig0 * cmath.exp(-1j * theta0)
        Is0dq = Is0 * cmath.exp(-1j * theta0)
        Vm0dq = Vm0 * cmath.exp(-1j * theta0)

        self.theta0 = theta0
        self.ved0, self.veq0 = Eg0dq.real, Eg0dq.imag
        self.igd0, self.igq0 = Ig0dq.real, Ig0dq.imag
        self.isd0, self.isq0 = Is0dq.real, Is0dq.imag
        vmd0, vmq0 = Vm0dq.real, Vm0dq.imag

        # Grid-side voltage: raw node's voltage in this unit's own dq frame,
        # and in the global dq frame (see components/sm.py for the derivation).
        grid_v = v_grid_pu * cmath.exp(1j * angle_grid_rad)
        vg_local = grid_v * cmath.exp(-1j * theta0)
        self.vgd0, self.vgq0 = vg_local.real, vg_local.imag
        vg_global = grid_v * cmath.exp(-1j * theta_g_rad)
        self.vgd_g0, self.vgq_g0 = vg_global.real, vg_global.imag
        self.theta_g0 = theta_g_rad

        self.p_ref0 = p_terminal_pu
        self.q_ref0 = q_terminal_pu
        self.w_ref0 = 1.0
        self.ve_ref0 = v_terminal_pu  # ved_ref_0 = ve_ref_0 in GFM_subs.m
        self.vdc_ref0 = 1.0
        self.qm0 = self.q_ref0
        self.pm0 = self.p_ref0
        self.vdc0 = self.vdc_ref0
        self.w0 = 1.0
        self.Dw0 = 0.0

        vdc_ref0 = self.vdc_ref0
        self.md0, self.mq0 = vmd0 / vdc_ref0, vmq0 / vdc_ref0
        self.idc_ref0 = self.p_ref0 / vdc_ref0
        self.idc0 = self.idc_ref0

        Kffi_, Kffv_, Cf_, wff_ = params["Kffi"], params["Kffv"], params["Cf"], params["wff"]
        self.M_VLd0 = self.isd0 - self.igd0 * Kffi_ + self.veq0 * Cf_ * wff_
        self.M_VLq0 = self.isq0 - self.igq0 * Kffi_ - self.ved0 * Cf_ * wff_
        self.M_CLd0 = vmd0 - self.ved0 * Kffv_ + self.isq0 * Lf_ * wff_
        self.M_CLq0 = vmq0 - self.veq0 * Kffv_ - self.isd0 * Lf_ * wff_


def _gfm_operating_subs(op: GfmOperatingPoint) -> dict:
    p = op.p
    return equilibrium_subs({
        wb: p["wb"], wff: p["wff"], Rf: p["Rf"], Lf: p["Lf"], Cf: p["Cf"], Rt: p["Rt"], Lt: p["Lt"],
        mp: p["mp"], nq: p["nq"], wf: p["wf"], KpVL: p["KpVL"], KiVL: p["KiVL"], Kffi: p["Kffi"],
        KpCL: p["KpCL"], KiCL: p["KiCL"], Kffv: p["Kffv"], Cdc: p["Cdc"], Gdc: p["Gdc"],
        Kpdc: p["Kpdc"], Tdc: p["Tdc"],
        isd: op.isd0, isq: op.isq0, igd: op.igd0, igq: op.igq0, ved: op.ved0, veq: op.veq0,
        vdc: op.vdc0, idc: op.idc0, pm: op.pm0, theta: op.theta0, qm: op.qm0,
        M_VLd: op.M_VLd0, M_VLq: op.M_VLq0, M_CLd: op.M_CLd0, M_CLq: op.M_CLq0,
        md: op.md0, mq: op.mq0, w: op.w0,
        p_ref: op.p_ref0, q_ref: op.q_ref0, ve_ref: op.ve_ref0, w_ref: op.w_ref0, vdc_ref: op.vdc_ref0,
        theta_g: op.theta_g0, vgd_g: op.vgd_g0, vgq_g: op.vgq_g0,
    })


def linearize_gfm(op: GfmOperatingPoint, modes: Mapping[str, str] | None = None) -> LinearComponent:
    return gfm_dae(mode_key(modes)).linearize(_gfm_operating_subs(op))


def gfm_nonlinear_point(op: GfmOperatingPoint, modes: Mapping[str, str] | None = None) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`gfm_nonlinear_funcs`."""
    return gfm_dae(mode_key(modes)).point_from_subs(_gfm_operating_subs(op))
