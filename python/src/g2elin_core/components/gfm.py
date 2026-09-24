"""Grid-forming converter, ported from ``Functions/symGFM_types.m`` +
``Functions/GFM_subs.m``.

All five of that file's outer power-control laws are here -- droop, droop
behind a filter, dVOC, VSM and matching control -- selected per unit through
``DerUnit.controller``. They differ only in the outer loop (see
:func:`_outer_loop`); the cascaded voltage and current loops, the LC filter
and the DC link are shared.

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
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache

import sympy as sp

from .base import (
    ComponentDAE, LinearComponent, NonlinearFuncs, NonlinearJacobians,
    RebuiltWithParams,
    apply_reduction, build_dae, equilibrium_subs, mode_key,
)

wb, wff, Rf, Lf, Cf, Rt, Lt, mp, nq, wf, KpVL, KiVL, Kffi, KpCL, KiCL, Kffv, Cdc, Gdc, Kpdc, Tdc = (
    sp.symbols("wb wff Rf Lf Cf Rt Lt mp nq wf KpVL KiVL Kffi KpCL KiCL Kffv Cdc Gdc Kpdc Tdc")
)
isd, isq, igd, igq, ved, veq, vdc, idc = sp.symbols("isd isq igd igq ved veq vdc idc")
pm, theta, qm, M_VLd, M_VLq, M_CLd, M_CLq = sp.symbols("pm theta qm M_VLd M_VLq M_CLd M_CLq")
md, mq, w = sp.symbols("md mq w")

# The outer power-control laws other than plain droop (symGFM_types.m). Each
# brings states and parameters of its own; the names are that file's.
Dw, Phi, vdc_m, ved_ref_s = sp.symbols("Dw Phi vdc_m ved_ref")
wc, eta, alfa, J, Dp, K, Dq, K_theta = sp.symbols("wc eta alfa J Dp K Dq K_theta")
p_ref, q_ref, ve_ref, w_ref, vdc_ref, theta_g, vgd_g, vgq_g = sp.symbols(
    "p_ref q_ref ve_ref w_ref vdc_ref theta_g vgd_g vgq_g"
)

#: The states every control law shares, in order. The outer law's own states
#: are spliced in between the DC link and the voltage loop, which is where
#: symGFM_types.m puts them.
_PLANT_STATE_NAMES = ["i_sd", "i_sq", "i_gd", "i_gq", "v_ed", "v_eq", "v_dc", "i_dc"]
_INNER_STATE_NAMES = ["M_VLd", "M_VLq", "M_CLd", "M_CLq"]


class GfmControlKind(str, Enum):
    """Which outer power-control law a grid-forming converter runs.

    These are ``symGFM_types.m``'s ``P_type`` cases. They differ in three
    things and nothing else: which states the outer loop carries, how it
    forms the frequency deviation ``Dw``, and how it forms the voltage
    reference the cascaded loops then track. Everything downstream -- the
    voltage loop, the current loop, the filter, the DC link -- is identical.
    """

    DROOP = "droop"
    DROOP_FILTERED = "droop_filtered"
    DVOC = "dvoc"
    VSM = "vsm"
    MATCHING = "matching"


@dataclass(frozen=True)
class _OuterLoop:
    """One power-control law, as the pieces :func:`gfm_dae` splices in."""

    states: tuple[tuple[sp.Symbol, sp.Expr, str], ...]
    dw: sp.Expr        # the frequency deviation the angle integrates
    ved_ref: sp.Expr   # the d-axis voltage reference the voltage loop tracks
    veq_ref: sp.Expr = sp.Integer(0)


def _outer_loop(kind: GfmControlKind, *, p: sp.Expr, q: sp.Expr) -> _OuterLoop:
    """The outer power-control law, transcribed from ``symGFM_types.m``.

    ``theta`` belongs to every one of them -- integrating the frequency is
    what makes a converter grid-*forming* -- so each law lists it in the
    place that file's own ``stateVec`` puts it, rather than the caller
    guessing a position that happens to be right for four of the five.
    """
    dpm, dqm = wf * (p - pm), wf * (q - qm)
    p_filter, q_filter = (pm, dpm, "p_m"), (qm, dqm, "q_m")
    angle = (theta, wb * w, "theta")
    droop_v_ref = ve_ref + (q_ref - qm) * nq

    if kind is GfmControlKind.DROOP:
        return _OuterLoop(
            states=(p_filter, angle, q_filter),
            dw=mp * (p_ref - pm), ved_ref=droop_v_ref,
        )

    if kind is GfmControlKind.DROOP_FILTERED:
        # The droop law itself behind a first-order filter, which is what
        # gives the converter a second-order (inertia-like) power response
        # rather than the first-order one plain droop has.
        return _OuterLoop(
            states=(p_filter, (Dw, mp * wc * (p_ref - pm) - wc * Dw, "dw"), angle, q_filter),
            dw=Dw, ved_ref=droop_v_ref,
        )

    if kind is GfmControlKind.DVOC:
        # Dispatchable virtual oscillator control: the voltage *amplitude* is
        # a state with dynamics of its own, driven by the reactive error and
        # pulled back toward ve_ref by the alfa term, and both channels are
        # normalised by the square of the amplitude.
        dved_ref = (
            eta * ((q_ref / ve_ref**2) - (qm / ved_ref_s**2))
            + (eta * alfa / ve_ref**2) * (ve_ref**2 - ved_ref_s**2)
        ) * ved_ref_s
        return _OuterLoop(
            states=(p_filter, angle, q_filter, (ved_ref_s, dved_ref, "v_ed_ref")),
            dw=eta * ((p_ref / ve_ref**2) - (pm / ved_ref_s**2)),
            ved_ref=ved_ref_s,
        )

    if kind is GfmControlKind.VSM:
        # A virtual synchronous machine: a swing equation on Dw with inertia
        # J and damping Dp, and a flux state Phi whose product with speed is
        # the voltage reference. It measures p and q *directly*, so it has no
        # power filters at all -- the inertia is what smooths the response.
        return _OuterLoop(
            states=(
                (Dw, (1 / (J * wff)) * (p_ref - p) - (Dp / J) * Dw, "dw"),
                angle,
                (Phi, (1 / K) * (q_ref - q) + (Dq / K) * (ve_ref - ved), "Phi"),
            ),
            dw=Dw, ved_ref=w * Phi,
        )

    # Matching control: the DC-link voltage *is* the frequency signal, which
    # is what a machine's speed does physically. No power measurement at all.
    return _OuterLoop(
        states=((vdc_m, wf * (vdc - vdc_m), "v_dc_m"), angle),
        dw=K_theta * (vdc_m - vdc_ref), ved_ref=ve_ref,
    )


#: The parameters each control law owns, under symGFM_types.m's own names.
#: A converter runs one law, so only its parameters are in its dict.
OUTER_SYMBOLS: dict[str, tuple[sp.Symbol, ...]] = {
    GfmControlKind.DROOP: (mp, nq, wf),
    GfmControlKind.DROOP_FILTERED: (mp, nq, wf, wc),
    GfmControlKind.DVOC: (eta, alfa, wf),
    GfmControlKind.VSM: (J, Dp, K, Dq),
    GfmControlKind.MATCHING: (K_theta, wf),
}
OUTER_PARAM_NAMES = {k: tuple(x.name for x in v) for k, v in OUTER_SYMBOLS.items()}


@lru_cache(maxsize=64)
def gfm_dae(
    modes: tuple[tuple[str, str], ...] = (),
    controller: GfmControlKind = GfmControlKind.DROOP,
) -> ComponentDAE:
    """``modes`` is a :func:`~g2elin_core.components.base.mode_key` tuple
    naming the states this converter gives up -- see
    :mod:`g2elin_core.reduction`.

    ``controller`` chooses the outer power-control law. The default is plain
    droop, the 15-state model this component has always built. The others
    carry different outer states: droop-behind-a-filter and dVOC add one to
    droop's three, VSM and matching drop the power filters entirely for two
    and three states respectively.
    """
    controller = GfmControlKind(controller)
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
    outer = _outer_loop(controller, p=p, q=q)

    dM_VLd = KiVL * (outer.ved_ref - ved)
    dM_VLq = KiVL * (outer.veq_ref - veq)
    isd_ref = KpVL * (outer.ved_ref - ved) + M_VLd + Kffi * igd - wff * Cf * veq
    isq_ref = KpVL * (outer.veq_ref - veq) + M_VLq + Kffi * igq + wff * Cf * ved
    dM_CLd = KiCL * (isd_ref - isd)
    dM_CLq = KiCL * (isq_ref - isq)

    plant_states = list(zip(
        [isd, isq, igd, igq, ved, veq, vdc, idc],
        [disd, disq, digd, digq, dved, dveq, dvdc, didc],
        _PLANT_STATE_NAMES,
    ))
    inner_states = list(zip(
        [M_VLd, M_VLq, M_CLd, M_CLq], [dM_VLd, dM_VLq, dM_CLd, dM_CLq], _INNER_STATE_NAMES,
    ))
    all_states = plant_states + list(outer.states) + inner_states

    alg = [
        md - (1 / vdc) * (KpCL * (isd_ref - isd) + Kffv * ved - wff * Lf * isq + M_CLd),
        mq - (1 / vdc) * (KpCL * (isq_ref - isq) + Kffv * veq + wff * Lf * isd + M_CLq),
        w - outer.dw - w_ref,
    ]

    igd_g_out = igd * sp.cos(theta - theta_g) - igq * sp.sin(theta - theta_g)
    igq_g_out = igd * sp.sin(theta - theta_g) + igq * sp.cos(theta - theta_g)
    Vt = sp.sqrt(ved**2 + veq**2)
    output_vec = [p, q, w, Vt, igd_g_out, igq_g_out]

    state_vec, diffeq_vec, state_names, alg_vec, alg, output_vec = apply_reduction(
        states=all_states,
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
def gfm_nonlinear_funcs(
    modes: tuple[tuple[str, str], ...] = (),
    controller: GfmControlKind = GfmControlKind.DROOP,
) -> NonlinearFuncs:
    return gfm_dae(modes, controller).nonlinear_funcs()


@lru_cache(maxsize=16)
def gfm_nonlinear_jacobians(
    modes: tuple[tuple[str, str], ...] = (),
    controller: GfmControlKind = GfmControlKind.DROOP,
) -> NonlinearJacobians:
    return gfm_dae(modes, controller).nonlinear_jacobians()


class GfmOperatingPoint(RebuiltWithParams):
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
        controller: GfmControlKind = GfmControlKind.DROOP,
    ):
        # locals() here, before anything else runs, is exactly the
        # arguments -- see RebuiltWithParams.
        self._built_from = {k: v for k, v in locals().items() if k != "self"}
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
        self.controller = GfmControlKind(controller)
        # The outer laws other than droop carry a state of their own. Each
        # sits where its own equation is at rest with p = p_ref, q = q_ref
        # and the terminal voltage where the power flow put it:
        #   dVOC   the amplitude state settles at the reference amplitude,
        #          which zeroes both its reactive error and its alfa pull;
        #   VSM    ved_ref = w*Phi has to equal that same amplitude, and
        #          w is 1 at rest;
        #   Match. the measured DC voltage is the DC voltage, which the DC
        #          loop holds at its reference.
        self.ved_ref0 = self.ve_ref0
        self.Phi0 = self.ve_ref0 / self.w0
        self.vdc_m0 = self.vdc0

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
        KpVL: p["KpVL"], KiVL: p["KiVL"], Kffi: p["Kffi"],
        KpCL: p["KpCL"], KiCL: p["KiCL"], Kffv: p["Kffv"], Cdc: p["Cdc"], Gdc: p["Gdc"],
        Kpdc: p["Kpdc"], Tdc: p["Tdc"],
        # The outer law's own parameters, under its own names. Only the
        # chosen law's are in `p`, and only its symbols reach the equations.
        **{sym: p[sym.name] for sym in OUTER_SYMBOLS[op.controller]},
        isd: op.isd0, isq: op.isq0, igd: op.igd0, igq: op.igq0, ved: op.ved0, veq: op.veq0,
        vdc: op.vdc0, idc: op.idc0, pm: op.pm0, theta: op.theta0, qm: op.qm0,
        M_VLd: op.M_VLd0, M_VLq: op.M_VLq0, M_CLd: op.M_CLd0, M_CLq: op.M_CLq0,
        # The outer laws' own states. Only the chosen law's appear in its
        # equations, so the rest are simply unused entries here.
        Dw: op.Dw0, Phi: op.Phi0, vdc_m: op.vdc_m0, ved_ref_s: op.ved_ref0,
        md: op.md0, mq: op.mq0, w: op.w0,
        p_ref: op.p_ref0, q_ref: op.q_ref0, ve_ref: op.ve_ref0, w_ref: op.w_ref0, vdc_ref: op.vdc_ref0,
        theta_g: op.theta_g0, vgd_g: op.vgd_g0, vgq_g: op.vgq_g0,
    })


def linearize_gfm(op: GfmOperatingPoint, modes: Mapping[str, str] | None = None) -> LinearComponent:
    return gfm_dae(mode_key(modes), op.controller).linearize(_gfm_operating_subs(op))


def gfm_nonlinear_point(op: GfmOperatingPoint, modes: Mapping[str, str] | None = None) -> tuple:
    """``(x0, z0, u0, p0)`` for :func:`gfm_nonlinear_funcs`."""
    return gfm_dae(mode_key(modes), op.controller).point_from_subs(_gfm_operating_subs(op))
