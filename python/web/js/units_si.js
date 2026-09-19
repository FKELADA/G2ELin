// Per-unit <-> SI conversions for the editor. Per unit is what the network
// stores; the SI value is always derived from it with the current base
// values, so changing a bus voltage or the base power updates every SI figure.
//
// Bases (g2elin_core/pu_base.py, amplitude-invariant, but Zb reduces to the
// usual form): Zb = Un^2/Sn [ohm], Lb = Zb/wb [H], Cb = 1/(Zb wb) [F];
// DC side: Zb_dc = 8/3 Zb, Cb_dc = 3/8 Cb. Lines use the same Zb as the
// power-flow adapter (vn_kv^2 / sn_mva at the line's voltage level).

function puBases(unKv) {
  const net = state.network;
  const wb = 2 * Math.PI * net.f_hz;
  const zb = (unKv * unKv) / net.sn_mva;
  return { wb, zb, lb: zb / wb, cb: 1 / (zb * wb), zbDc: (8 / 3) * zb, cbDc: (3 / 8) * (1 / (zb * wb)) };
}
function busKv(id) { return state.network.buses.find(b => b.id === id)?.vn_kv; }

// [{unit, factor, title}] for an editor field, with si = pu * factor.
function siSpecs(kind, obj, field) {
  if (!state.network || !obj) return [];
  if (kind === "line" && ["r_pu", "x_pu", "b_pu"].includes(field)) {
    const kv = busKv(obj.from_bus);
    if (!kv) return [];
    const { zb } = puBases(kv), len = obj.length_km;
    const [u, f] = field === "b_pu" ? ["µS", 1e6 / zb] : ["Ω", zb];
    return [{ unit: u, factor: f, title: "whole line" }, ...(len > 0 ? [{ unit: `${u}/km`, factor: f / len, title: "per km" }] : [])];
  }
  if (kind === "transformer" && ["r_pu", "x_pu"].includes(field)) {
    const kv = busKv(obj.hv_bus);
    if (!kv || !(obj.sn_mva > 0)) return [];
    return [{ unit: "Ω (HV side)", factor: (kv * kv) / obj.sn_mva, title: "referred to the HV side, on the transformer rating" }];
  }
  if (kind === "der" && field === "v_set_pu") {
    const kv = busKv(obj.bus);
    return kv ? [{ unit: "kV", factor: kv, title: "line-to-line" }] : [];
  }
  if (kind === "der" && field === "xd_pu") {
    const kv = busKv(obj.bus);
    return kv ? [{ unit: "Ω", factor: puBases(kv).zb, title: "at the unit's terminal voltage" }] : [];
  }
  return [];
}

// A unit's electrical parameters (unit_params.js), on the unit's own base
// (its terminal bus voltage, the network base power). Rotor quantities of the
// synchronous machine are referred to the stator.
const UNIT_PARAM_SI = {
  R: { keys: ["Rf", "Rt", "Ra", "Rfd", "R1d", "R1q", "R2q", "RL_pu"], unit: "Ω", f: b => b.zb },
  L: { keys: ["Lf", "Lt", "Ll", "Lad", "Laq", "Lfd", "L1d", "L1q", "L2q"], unit: "mH", f: b => b.lb * 1e3 },
  C: { keys: ["Cf"], unit: "µF", f: b => b.cb * 1e6 },
  Cdc: { keys: ["Cdc"], unit: "mF", f: b => b.cbDc * 1e3 },
  Gdc: { keys: ["Gdc"], unit: "µS", f: b => 1e6 / b.zbDc },
};
function unitParamSi(der, key) {
  const kv = busKv(der.bus);
  if (!kv) return null;
  for (const g of Object.values(UNIT_PARAM_SI)) if (g.keys.includes(key)) return { unit: g.unit, factor: g.f(puBases(kv)) };
  return null;
}

const siFmt = v => (Number.isFinite(v) ? +v.toPrecision(6) : "");
