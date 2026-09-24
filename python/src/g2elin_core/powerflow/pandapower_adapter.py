"""Static power flow via pandapower.

Replaces ``Functions/Power_Fl.m`` / ``Functions/Load_flow.m`` (a hand-rolled
sparse Newton-Raphson solver). pandapower's ``gen``/``ext_grid`` elements
give native PV-bus (voltage-controlled) and slack-bus behavior, so unlike
power-grid-model this needs no outer voltage-control loop — see the
"Static power flow" section of the migration plan for why PGM was dropped.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import pandapower as pp
import pandas as pd

from g2elin_core.network.breakers import (
    NoReferenceUnit, service_state, shunt_is_reactor, shunt_reactor_pq_mw,
)
from g2elin_core.network.schema import BusType, Network, UnitType

_LARGE_Q_LIMIT_MVAR = 9999.0  # matches Power_Fl.m's default qg_max/qg_min when unspecified
_PLACEHOLDER_MAX_I_KA = 999.0  # thermal rating isn't modeled upstream; not used by the LF solver


def build_pandapower_net(network: Network) -> tuple[pp.pandapowerNet, dict[int, int]]:
    """Build a pandapower network from a :class:`Network`.

    Returns the net and a mapping from ``Network`` bus ids to pandapower's
    internal bus indices (pandapower reindexes on its own, so this mapping
    is needed to translate results back).

    Elements switched out by open breakers (see ``network.breakers``) are
    created but marked out of service, and de-energized buses too, so every
    result table keeps the network's own element order.

    Each energized island gets its own reference (``ext_grid``), as every
    load-flow tool requires -- the designated slack in its own island, the
    largest grid former in any other (``breakers.service_state``).
    """
    st = service_state(network)
    if not st.references:
        raise NoReferenceUnit(
            "the power flow has no reference: open breakers have left no synchronous machine, grid-forming "
            "converter or infinite bus in service to set a voltage and a frequency"
        )
    references = set(st.references)
    net = pp.create_empty_network(name=network.name, f_hz=network.f_hz, sn_mva=network.sn_mva)

    bus_index: dict[int, int] = {
        bus.id: pp.create_bus(
            net, vn_kv=bus.vn_kv, name=bus.name or str(bus.id), in_service=bus.id in st.energized_buses,
        )
        for bus in network.buses
    }

    for li, line in enumerate(network.lines):
        vn_kv = network.bus(line.from_bus).vn_kv
        if vn_kv != network.bus(line.to_bus).vn_kv:
            raise ValueError(
                f"line {line.from_bus}->{line.to_bus} spans two voltage levels; "
                "model the voltage change as a Transformer instead"
            )
        z_base_ohm = vn_kv**2 / network.sn_mva
        omega = 2 * math.pi * network.f_hz
        b_siemens_total = line.b_pu / z_base_ohm
        c_farad_total = b_siemens_total / omega

        pp.create_line_from_parameters(
            net,
            from_bus=bus_index[line.from_bus],
            to_bus=bus_index[line.to_bus],
            length_km=line.length_km,
            r_ohm_per_km=line.r_pu * z_base_ohm / line.length_km,
            x_ohm_per_km=line.x_pu * z_base_ohm / line.length_km,
            c_nf_per_km=(c_farad_total / line.length_km) * 1e9,
            max_i_ka=_PLACEHOLDER_MAX_I_KA,
            name=line.name,
            in_service=st.lines[li],
        )

    for ti, tr in enumerate(network.transformers):
        vk_percent = math.hypot(tr.r_pu, tr.x_pu) * 100.0
        vkr_percent = tr.r_pu * 100.0
        pp.create_transformer_from_parameters(
            net,
            hv_bus=bus_index[tr.hv_bus],
            lv_bus=bus_index[tr.lv_bus],
            sn_mva=tr.sn_mva,
            # An off-nominal ratio is the HV winding rated above its bus:
            # pandapower's ratio is vn_hv_kv/vn_lv_kv against the two buses'
            # own nominal voltages, so scaling vn_hv_kv by tap_ratio is exactly
            # the a of the dynamic model (network/breakers.transformer_ratio).
            vn_hv_kv=network.bus(tr.hv_bus).vn_kv * tr.tap_ratio,
            vn_lv_kv=network.bus(tr.lv_bus).vn_kv,
            vkr_percent=vkr_percent,
            vk_percent=max(vk_percent, 1e-6),
            shift_degree=tr.shift_degree,
            pfe_kw=0.0,
            i0_percent=0.0,
            name=tr.name,
            in_service=st.transformers[ti],
        )

    for li, load in enumerate(network.loads):
        pp.create_load(
            net, bus=bus_index[load.bus], p_mw=load.p_mw, q_mvar=load.q_mvar, name=load.name,
            in_service=st.loads[li],
        )

    # Shunts: pandapower's own sign convention is this schema's (q_mvar > 0
    # absorbs), and its shunt is rated at the bus's nominal voltage, which is
    # the voltage Shunt.q_mvar is quoted at.
    for si, shunt in enumerate(network.shunts):
        # A capacitor bank is a pure susceptance. A reactor has an X/R, and
        # gets the P and Q of the very R-X branch its dynamic model uses
        # (network/breakers.shunt_reactor_pq_mw), so power flow and dynamics
        # agree on what it draws.
        if shunt_is_reactor(shunt):
            p_mw, q_mvar = shunt_reactor_pq_mw(shunt, network.sn_mva)
        else:
            p_mw, q_mvar = 0.0, shunt.q_mvar
        pp.create_shunt(
            net, bus=bus_index[shunt.bus], q_mvar=q_mvar, p_mw=p_mw, name=shunt.name,
            in_service=st.shunts[si],
        )

    for der in network.der_units:
        b = bus_index[der.bus]
        name = f"der{der.id}"
        on = st.der_units[der.id]
        if der.id in references:
            # This island's reference: a slack bus for it, whatever the unit's
            # own declared bus type (it may be a PV machine picking up the role).
            pp.create_ext_grid(
                net, bus=b, vm_pu=der.v_set_pu, va_degree=der.angle_set_deg, name=name
            )
        elif der.bus_type is BusType.PV:
            pp.create_gen(
                net,
                bus=b,
                p_mw=der.p_set_mw,
                vm_pu=der.v_set_pu,
                min_q_mvar=-_LARGE_Q_LIMIT_MVAR,
                max_q_mvar=_LARGE_Q_LIMIT_MVAR,
                name=name,
                in_service=on,
            )
        else:  # PQ-dispatched DER (e.g. a GFL unit not under voltage control)
            pp.create_sgen(net, bus=b, p_mw=der.p_set_mw, q_mvar=der.q_set_mvar, name=name, in_service=on)

        if der.p_cons_mw or der.q_cons_mvar:
            pp.create_load(
                net, bus=b, p_mw=der.p_cons_mw, q_mvar=der.q_cons_mvar, name=f"{name}_aux_load", in_service=on,
            )

    return net, bus_index


@dataclass
class PowerFlowResult:
    """Solved power flow, with results indexed by the original ``Network`` bus id."""

    network: Network
    net: pp.pandapowerNet
    bus_index: dict[int, int]
    converged: bool

    def bus_table(self) -> pd.DataFrame:
        """Per-bus results, laid out to mirror ``Power_Fl.m``'s ``bus_sol`` columns."""
        rows = []
        inv_index = {v: k for k, v in self.bus_index.items()}
        for pp_idx, row in self.net.res_bus.iterrows():
            bus_id = inv_index[pp_idx]
            p_gen = -row.p_mw  # pandapower's res_bus convention is load-positive
            q_gen = -row.q_mvar
            rows.append(
                {
                    "bus": bus_id,
                    "vm_pu": row.vm_pu,
                    "va_degree": row.va_degree,
                    "p_net_gen_mw": p_gen,
                    "q_net_gen_mvar": q_gen,
                }
            )
        return pd.DataFrame(rows).sort_values("bus").reset_index(drop=True)

    def total_losses_mw(self) -> float:
        return float(self.net.res_line.pl_mw.fillna(0).sum() + self.net.res_trafo.pl_mw.fillna(0).sum())

    def _element_rows(self, res_df: pd.DataFrame, element_df: pd.DataFrame, bus_cols: list[str]) -> list[dict]:
        """Joins a pandapower ``res_*`` table with its element table's
        ``name``/bus columns, remapping bus columns from pandapower's own
        internal indices back to this project's own ``Network`` bus ids
        (same as :meth:`bus_table`, generalized to every other element kind
        ``build_pandapower_net`` creates).
        """
        inv_index = {v: k for k, v in self.bus_index.items()}
        rows = []
        for idx in res_df.index:
            row = {"name": element_df.at[idx, "name"], "in_service": bool(element_df.at[idx, "in_service"])}
            for col in bus_cols:
                row[col] = inv_index[int(element_df.at[idx, col])]
            row.update(res_df.loc[idx].to_dict())
            rows.append(row)
        return rows

    def line_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_line, self.net.line, ["from_bus", "to_bus"]))

    def trafo_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_trafo, self.net.trafo, ["hv_bus", "lv_bus"]))

    def load_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_load, self.net.load, ["bus"]))

    def gen_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_gen, self.net.gen, ["bus"]))

    def sgen_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_sgen, self.net.sgen, ["bus"]))

    def ext_grid_table(self) -> pd.DataFrame:
        return pd.DataFrame(self._element_rows(self.net.res_ext_grid, self.net.ext_grid, ["bus"]))


def run_power_flow(network: Network, **runpp_kwargs) -> PowerFlowResult:
    """Build and solve the power flow for ``network``."""
    net, bus_index = build_pandapower_net(network)
    runpp_kwargs.setdefault("calculate_voltage_angles", True)
    try:
        pp.runpp(net, **runpp_kwargs)
        converged = True
    except pp.powerflow.LoadflowNotConverged:
        converged = False
    return PowerFlowResult(network=network, net=net, bus_index=bus_index, converged=converged)
