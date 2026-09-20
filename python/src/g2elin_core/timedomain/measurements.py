"""Measurement outputs for EMT simulations: quantities an engineer would read
off a meter -- power flows, frequency, voltage magnitude and angle, and
instantaneous 3-phase voltages -- computed from each element's own model
variables at every sample.

All model variables live in a common dq frame -- the frame block a network's
island is referenced to (``components/frame.py``), or the slack unit itself
in the MATLAB-compatible mode. Its ``theta``/``wr`` are read here. From them:

- **Bus voltage** ``v = vd + j vq`` (a network node's two states):
  magnitude ``|v|`` [pu]; angle ``atan2(vq, vd) + theta_g0`` [deg], where
  ``theta_g0`` is the frame rotation the operating point applied, so angles
  equal the power flow's at t = 0; instantaneous frequency
  ``f = f_n (wg + dphi/dt / wb)`` [Hz] with ``wg`` the frame speed and
  ``dphi/dt = (vd dvq/dt - vq dvd/dt)/|v|^2`` taken analytically from the
  node equations -- it carries the network's fast electromagnetic
  transients (several Hz swings in the first ~100 ms), so the default
  "measured" frequency passes it through a one-cycle first-order filter
  (time constant ``1/f_n``), as a frequency meter or PLL would;
  instantaneous phase voltages ``v_a,b,c = |v| cos(Theta - k 2pi/3)`` [pu,
  peak phase value, amplitude-invariant pu] with ``Theta = (theta_frame(t) -
  theta_frame(0)) + phi + theta_g0``, the voltage phasor's absolute angle.
- **Power** ``S = v i* = (vd id + vq iq) + j (vq id - vd iq)`` [pu on the
  network base]:
  lines at both ends (the line current flows from-bus -> to-bus; "from" is
  what enters the line, "to" what it delivers; each end includes half the
  line's own charging ``|v|^2 b/2``, which the model keeps in the bus
  nodes, so these are the pi-line's terminal flows -- the power flow's
  convention), loads (consumed), units at
  their grid-side bus (injected -- this is also their transformer's HV-side
  flow, since a unit's transformer is part of the unit's own model).
- **Unit frequency** [Hz]: the unit's own speed/frequency output (rotor speed
  for an SM, droop frequency for a GFM, PLL frequency for a GFL, the grid
  frequency for an infinite bus) times ``f_n``.

A unit's LV terminal bus is not a node of the dynamic model (the unit's
model includes its transformer and connects to the HV bus), so it has no bus
measurements of its own; the unit's own ``V_t`` output gives its terminal
voltage magnitude.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np

from g2elin_core.interconnect.assemble import _OUTG_PORTS
from g2elin_core.interconnect.network_assembly import _TYPE_LABEL
from g2elin_core.network.breakers import block_labels

TWO_PI_3 = 2 * math.pi / 3
_UNIT_SPEED_OUTPUT = {"sm": "w_r", "gfm": "w", "gfl": "w_pll"}


@dataclass(frozen=True)
class Measurement:
    name: str
    group: str  # e.g. "Bus 4", "Line #0 (1 → 4)"
    label: str  # e.g. "Voltage magnitude"
    unit: str


class MeasurementSet:
    """The measurements available for one EMT model, and how to evaluate
    them from a solved ``(x, z, u)`` sample."""

    def __init__(self, model) -> None:
        self.model = model
        net = model.network
        self.f_n = net.f_hz
        self.wb = 2 * math.pi * net.f_hz
        by_name = {b.name: b for b in model.blocks}
        # Element numbers are those of the full network (a network with
        # elements switched out keeps them, see network.breakers).
        lab = block_labels(net)
        self.nodes = {bus_id: by_name[f"Nd_{bus_id}"] for bus_id in (b.id for b in net.buses) if f"Nd_{bus_id}" in by_name}
        self.lines = [(lab.line[i], ln, by_name[f"Ln_{lab.line[i] + 1}"]) for i, ln in enumerate(net.lines)]
        self.loads = [(lab.load[i], ld, by_name[f"Ld_{lab.load[i] + 1}"]) for i, ld in enumerate(net.loads)]
        self.units = []
        hv_of = {t.lv_bus: (lab.transformer[j], t) for j, t in enumerate(net.transformers)}
        for der in net.der_units:
            block = by_name[lab.der[der.id]]
            j, tr = hv_of.get(der.bus, (None, None))
            self.units.append((der, block, j, tr))
        # The frame everything is written in: its own block, or the slack
        # unit when the frame follows it (Network.frame_follows_slack). With
        # several islands the first frame is used for the shared quantities
        # (the absolute phase of the 3-phase waveforms); every other
        # measurement is local to its element.
        self.frame = next(
            (b for b in model.blocks if b.kind == "frame"),
            next((b for d, b, _, _ in self.units if b.kind.endswith("_slack")), None),
        )
        self.filtered: set[str] = set()      # names smoothed along time (see smooth())
        self.frequencies: set[str] = set()   # bus frequencies (ideal value before the disturbance)
        self._catalog = self._build_catalog()
        self._theta_frame0: float | None = None

    # --- catalogue -----------------------------------------------------------
    def _build_catalog(self) -> dict[str, tuple[Measurement, Callable]]:
        cat: dict[str, tuple[Measurement, Callable]] = {}

        def add(name, group, label, unit, fn):
            cat[name] = (Measurement(name, group, label, unit), fn)

        net = self.model.network
        for bus in net.buses:
            nb = self.nodes.get(bus.id)
            if nb is None:
                continue  # a unit's own terminal bus: not a node of the dynamic model
            g = f"Bus {bus.id}" + (f" ({bus.name})" if bus.name else "")
            add(f"V_{{bus{bus.id}}}", g, "Voltage magnitude", "pu", lambda c, nb=nb: abs(c.v(nb)))
            add(f"angle_{{bus{bus.id}}}", g, "Voltage angle", "deg", lambda c, nb=nb: c.angle_deg(nb))
            add(f"f_{{bus{bus.id}}}", g, "Frequency (measured, 1-cycle filter)", "Hz", lambda c, nb=nb: c.freq(nb))
            self.filtered.add(f"f_{{bus{bus.id}}}")
            add(f"f_inst_{{bus{bus.id}}}", g, "Frequency (instantaneous)", "Hz", lambda c, nb=nb: c.freq(nb))
            self.frequencies.update({f"f_{{bus{bus.id}}}", f"f_inst_{{bus{bus.id}}}"})
            for k, ph in enumerate("abc"):
                add(f"v_{ph}_{{bus{bus.id}}}", g, f"Phase-{ph} voltage (instantaneous)", "pu", lambda c, nb=nb, k=k: c.phase(nb, k))
        for i, ln, lb in self.lines:
            g = f"Line #{i} ({ln.from_bus} → {ln.to_bus})" + (f" {ln.name}" if ln.name else "")
            nj, nk = self.nodes.get(ln.from_bus), self.nodes.get(ln.to_bus)
            if nj is None or nk is None:
                continue
            add(f"P_from_{{line{i}}}", g, "Active power entering at the from-bus", "pu", lambda c, lb=lb, nj=nj: c.s(nj, c.x2(lb)).real)
            add(f"Q_from_{{line{i}}}", g, "Reactive power entering at the from-bus", "pu",
                lambda c, lb=lb, nj=nj, b=ln.b_pu: c.s(nj, c.x2(lb)).imag - abs(c.v(nj)) ** 2 * b / 2)
            add(f"P_to_{{line{i}}}", g, "Active power delivered at the to-bus", "pu", lambda c, lb=lb, nk=nk: c.s(nk, c.x2(lb)).real)
            add(f"Q_to_{{line{i}}}", g, "Reactive power delivered at the to-bus", "pu",
                lambda c, lb=lb, nk=nk, b=ln.b_pu: c.s(nk, c.x2(lb)).imag + abs(c.v(nk)) ** 2 * b / 2)
        for i, ld, lb in self.loads:
            nb = self.nodes.get(ld.bus)
            if nb is None:
                continue
            g = f"Load #{i} (bus {ld.bus})" + (f" {ld.name}" if ld.name else "")
            add(f"P_{{load{i}}}", g, "Active power consumed", "pu", lambda c, lb=lb, nb=nb: c.s(nb, c.x2(lb)).real)
            add(f"Q_{{load{i}}}", g, "Reactive power consumed", "pu", lambda c, lb=lb, nb=nb: c.s(nb, c.x2(lb)).imag)
        for der, ub, j, tr in self.units:
            nb = self.nodes.get(tr.hv_bus) if tr is not None else None
            g = f"Unit {der.id} ({_TYPE_LABEL[der.unit_type.value]}, {ub.name})"
            if nb is not None:
                add(f"P_grid_{{unit{der.id}}}", g, "Active power injected into the grid bus", "pu", lambda c, ub=ub, nb=nb: c.s(nb, c.igrid(ub)).real)
                add(f"Q_grid_{{unit{der.id}}}", g, "Reactive power injected into the grid bus", "pu", lambda c, ub=ub, nb=nb: c.s(nb, c.igrid(ub)).imag)
                tg = f"Transformer #{j} ({tr.hv_bus} → {tr.lv_bus})" + (f" {tr.name}" if tr.name else "")
                add(f"P_hv_{{trafo{j}}}", tg, "Active power delivered to the HV bus", "pu", lambda c, ub=ub, nb=nb: c.s(nb, c.igrid(ub)).real)
                add(f"Q_hv_{{trafo{j}}}", tg, "Reactive power delivered to the HV bus", "pu", lambda c, ub=ub, nb=nb: c.s(nb, c.igrid(ub)).imag)
            add(f"f_{{unit{der.id}}}", g, "Unit frequency (rotor / droop / PLL)", "Hz", lambda c, ub=ub: c.unit_freq(ub))
        return cat

    def catalog(self) -> list[Measurement]:
        return [m for m, _ in self._catalog.values()]

    def names(self) -> list[str]:
        return list(self._catalog)

    # --- along time ------------------------------------------------------------
    def smoother(self, names: list[str]) -> "Callable[[float, dict[str, float]], dict[str, float]]":
        """Stateful ``step(t, values) -> values`` applying the measurement
        filter (first order, time constant one cycle) to the filtered
        names, sample by sample -- exact for any, even irregular, spacing.
        Starts from the nominal frequency."""
        tau = 1.0 / self.f_n
        names = [n for n in names if n in self.filtered]
        state: dict[str, float] = {n: self.f_n for n in names}
        last = {"t": None}

        def step(t: float, values: dict[str, float]) -> dict[str, float]:
            if last["t"] is not None and names:
                a = 1.0 - math.exp(-max(t - last["t"], 0.0) / tau)
                for n in names:
                    state[n] += a * (values[n] - state[n])
            last["t"] = t
            out = dict(values)
            out.update({n: state[n] for n in names})
            return out

        return step

    def smooth(self, t: np.ndarray, series: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """The filtered names of ``series`` (sampled at ``t``) filtered along time."""
        step = self.smoother(list(series))
        rows = [step(float(t[k]), {n: float(v[k]) for n, v in series.items()}) for k in range(len(t))]
        return {n: np.array([r[n] for r in rows]) for n in series}

    def before_disturbance(self, names: list[str], t: np.ndarray) -> dict[str, np.ndarray]:
        """Each measurement over ``t`` (t <= 0) at the undisturbed operating
        point: constants, except the 3-phase voltages, which keep rotating at
        nominal frequency, and bus frequencies, which read the nominal value
        (the model's initial point isn't an exact equilibrium, so its
        instantaneous derivative would show a start-up artefact here)."""
        m = self.model
        x0 = m.initial_state()
        z0, u0 = m.solve_algebraic(x0, m.default_u_exo())
        measure = self.evaluator(names)
        out = {n: np.empty(len(t)) for n in names}
        for k, tk in enumerate(t):
            vals = measure(x0, z0, u0, frame_offset=self.wb * float(tk))
            for n in names:
                out[n][k] = self.f_n if n in self.frequencies else vals[n]
        return out

    # --- evaluation ------------------------------------------------------------
    def evaluator(self, names: list[str]) -> Callable[[np.ndarray, np.ndarray, np.ndarray], dict[str, float]]:
        """``measure(x, z, u) -> {name: value}`` for ``names`` (raises KeyError
        on an unknown name), for ``NonlinearNetworkModel.recover_signals``."""
        fns = [(n, self._catalog[n][1]) for n in names]
        if self._theta_frame0 is None:
            m = self.model
            x0 = m.initial_state()
            z0, u0 = m.solve_algebraic(x0, m.default_u_exo())
            self._theta_frame0 = _Sample(self, x0, z0, u0).frame_theta()

        def measure(x, z, u, frame_offset: float = 0.0) -> dict[str, float]:
            c = _Sample(self, x, z, u, frame_offset)
            return {n: float(fn(c)) for n, fn in fns}

        return measure


class _Sample:
    """One solved sample, with the per-block lookups the measurements use."""

    def __init__(self, ms: MeasurementSet, x, z, u, frame_offset: float = 0.0) -> None:
        self.ms, self.m = ms, ms.model
        self.x, self.x_list = x, ms.model._unpack_x(x)
        self.z_list, self.u_list = ms.model._unpack_z(z), ms.model._unpack_u(u)
        self.frame_offset = frame_offset
        self._y: dict[int, np.ndarray] = {}

    def _idx(self, b) -> int:
        return self.m.blocks.index(b)

    def y(self, b) -> np.ndarray:
        i = self._idx(b)
        if i not in self._y:
            self._y[i] = np.asarray(b.comp.h(self.x_list[i], self.z_list[i], self.u_list[i]), dtype=float)
        return self._y[i]

    def x2(self, b) -> complex:
        xs = self.x_list[self._idx(b)]
        return complex(xs[0], xs[1])

    def v(self, node) -> complex:
        return self.x2(node)

    def s(self, node, i: complex) -> complex:
        v = self.v(node)
        return complex(v.real * i.real + v.imag * i.imag, v.imag * i.real - v.real * i.imag)

    def igrid(self, b) -> complex:
        y, ports = self.y(b), _OUTG_PORTS[b.kind]
        n = b.comp.n_out_s
        return complex(y[n + ports["igd_g"]], y[n + ports["igq_g"]])

    def frame_speed(self) -> float:
        b = self.ms.frame
        return float(self.y(b)[b.comp.n_out_s + _OUTG_PORTS[b.kind]["wr"]])

    def frame_theta(self) -> float:
        b = self.ms.frame
        return float(self.y(b)[b.comp.n_out_s + _OUTG_PORTS[b.kind]["theta"]])

    def phi(self, node) -> float:
        v = self.v(node)
        return math.atan2(v.imag, v.real)

    def angle_deg(self, node) -> float:
        a = math.degrees(self.phi(node) + self.m.theta_g0)
        return (a + 180.0) % 360.0 - 180.0

    def freq(self, node) -> float:
        i = self._idx(node)
        dv = np.asarray(node.comp.f(self.x_list[i], self.z_list[i], self.u_list[i]), dtype=float)
        v = self.v(node)
        mag2 = v.real ** 2 + v.imag ** 2
        dphi = (v.real * dv[1] - v.imag * dv[0]) / mag2 if mag2 > 0 else 0.0
        return self.ms.f_n * (self.frame_speed() + dphi / self.ms.wb)

    def phase(self, node, k: int) -> float:
        theta = (self.frame_theta() - self.ms._theta_frame0) + self.frame_offset + self.phi(node) + self.m.theta_g0
        return abs(self.v(node)) * math.cos(theta - k * TWO_PI_3)

    def unit_freq(self, b) -> float:
        if b.kind.startswith("ib"):
            # An infinite bus runs at its own (exogenous) speed input.
            return self.ms.f_n * float(self.u_list[self.ms.model.blocks.index(b)][0])
        names = b.comp.output_names
        key = _UNIT_SPEED_OUTPUT.get(b.kind.removesuffix("_slack"))
        if key in names:
            return self.ms.f_n * float(self.y(b)[names.index(key)])
        return float("nan")
