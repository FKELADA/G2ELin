// Editable control & electrical parameters of one unit (SM / GFM / GFL), and
// a control-loop tuner for GFM/GFL. Edits are stored as overrides in
// der.params (only values that differ from the derived default); the
// backend applies them on top of sm_params()/gfm_params()/gfl_params().
//
// The tuner uses the same pole-placement formulas as operating_point.py,
// in both directions: response time t_r and damping ζ -> Kp, Ki, and back.

// The regulator models a machine can carry, and the parameters each one
// brings. Every model uses its *own* names -- nothing is shared between them
// and nothing is renamed to a common spelling -- so a machine's parameter set
// depends on which pair it carries. The backend's key set agrees
// (operating_point.SM_EXCITER_PARAMS / SM_PSS_PARAMS); an override under a
// name the chosen model does not have is rejected there.
const SM_REGULATORS = [
  ["governor", "Governor", {
    g2elin: ["Droop into a first-order lag", ["mp", "TG"]],
    none: ["None (constant mechanical power)", []],
  }],
  ["exciter", "Exciter / AVR", {
    g2elin: ["G2ELin original", ["Tr", "Ka", "Ta", "Ke", "Te", "Kfd", "Tfd"]],
    kundur: ["Thyristor with TGR (Kundur Fig. E12.9)", ["TR", "KA", "TA", "TB"]],
  }],
  ["pss", "Power system stabiliser", {
    g2elin: ["G2ELin original", ["K_PSS", "T_LP", "T_HP", "T1n", "T1d", "T2n", "T2d"]],
    kundur: ["Washout and two lead-lags (Kundur Fig. E12.9)", ["KSTAB", "TW", "T1", "T2", "T3", "T4"]],
    none: ["None (no stabiliser fitted)", []],
  }],
];

// A converter's outer power-control law, and the parameters each one brings.
// Same shape as SM_REGULATORS: the law is chosen on the unit, and only its
// own parameters are shown (and accepted by the backend).
const GFM_CONTROLLERS = [
  ["controller", "Power control law", {
    droop: ["Droop", ["mp", "nq", "wf"]],
    droop_filtered: ["Droop behind a filter (2nd-order response)", ["mp", "nq", "wf", "wc"]],
    dvoc: ["dVOC (dispatchable virtual oscillator)", ["eta", "alfa", "wf"]],
    vsm: ["VSM (virtual synchronous machine)", ["J", "Dp", "K", "Dq"]],
    matching: ["Matching (DC voltage sets frequency)", ["K_theta", "wf"]],
  }],
];

const PARAM_GROUPS = {
  sm: [
    ["Stator & rotor windings (pu)", ["Ra", "Ll", "Lad", "Laq", "Lfd", "Rfd", "L1d", "R1d", "L1q", "R1q", "L2q", "R2q"]],
    ["Mechanical", ["H", "KD"]],
    ["Transformer & auxiliary load (pu)", ["Rt", "Lt", "RL_pu"]],
  ],
  gfm: [
    ["LC filter & transformer (pu)", ["Rf", "Lf", "Cf", "Rt", "Lt"]],
    ["Voltage loop", ["KpVL", "KiVL", "Kffi"]],
    ["Current loop (inner)", ["KpCL", "KiCL", "Kffv"]],
    ["DC link", ["Cdc", "Gdc", "Kpdc", "Tdc"]],
  ],
  gfl: [
    ["LC filter & transformer (pu)", ["Rf", "Lf", "Cf", "Rt", "Lt"]],
    ["Current loop (inner)", ["KpCL", "KiCL", "Kffv"]],
    ["Phase-locked loop", ["Kppll", "Kipll"]],
    ["DC-voltage loop (outer)", ["Kpd", "Kid"]],
    ["Reactive-power loop (outer)", ["Kiq"]],
    ["DC link", ["Cdc", "Gdc", "Tdc"]],
  ],
};

// The two regulator groups a machine shows, for the models it carries.
function unitModelSlots(der) {
  if (der.unit_type === "sm") return SM_REGULATORS;
  if (der.unit_type === "gfm") return GFM_CONTROLLERS;
  return [];
}

function regulatorGroups(der) {
  if (!unitModelSlots(der).length) return [];
  // A model fitted as "none" brings no parameters, so it contributes no
  // group -- the picker for it still shows, above.
  return unitModelSlots(der).map(([slot, title, models]) => {
    const fallback = models.g2elin ? "g2elin" : "droop";
    const entry = models[der[slot] || fallback] || models[fallback];
    return [title, entry[1]];
  }).filter(([, keys]) => keys.length);
}

// The droop-equivalent tuning a set of parameters stands for. Display only:
// anything that *writes* goes through /api/units/retune, so the five laws'
// formulas live in one place. Kept in step with
// operating_point.gfm_outer_tuning.
function gfmTuning(p, law) {
  let mp, nq, wf;
  if (law === "droop" || law === "droop_filtered") { mp = p.mp; nq = p.nq; wf = p.wf; }
  else if (law === "dvoc") { mp = p.eta; nq = 1 / (2 * p.alfa); wf = p.wf; }
  else if (law === "vsm") { mp = 1 / p.Dp; nq = 1 / p.Dq; wf = p.Dp / p.J; }
  else if (law === "matching") { mp = p.Kpdc ? p.K_theta / p.Kpdc : NaN; nq = p.nq ?? 1e-4; wf = p.wf; }
  else return null;
  if (![mp, nq, wf].every(x => Number.isFinite(x) && x > 0)) return null;
  return { mp, nq, wf, H: 1 / (2 * mp * wf) };
}

const PARAM_HELP = {
  Ra: "stator resistance", Ll: "stator leakage inductance", Lad: "d-axis mutual inductance", Laq: "q-axis mutual inductance",
  Lfd: "field leakage inductance", Rfd: "field resistance", L1d: "d damper inductance", R1d: "d damper resistance",
  L1q: "1st q damper inductance", R1q: "1st q damper resistance", L2q: "2nd q damper inductance", R2q: "2nd q damper resistance",
  H: "inertia constant (s)", KD: "damping coefficient", mp: "frequency droop (pu)", TG: "governor time constant (s)",
  Tr: "voltage transducer time constant (s)", Ka: "exciter gain", Ta: "exciter time constant (s)", Ke: "exciter constant",
  Te: "exciter time constant (s)", Kfd: "stabilising feedback gain", Tfd: "stabilising feedback time constant (s)",
  K_PSS: "PSS gain", T_LP: "PSS low-pass time constant (s)", T_HP: "PSS washout time constant (s)",
  T1n: "lead-lag 1 numerator (s)", T1d: "lead-lag 1 denominator (s)", T2n: "lead-lag 2 numerator (s)", T2d: "lead-lag 2 denominator (s)",
  TR: "voltage transducer time constant (s)", KA: "exciter gain",
  TA: "transient gain reduction, lead (s)", TB: "transient gain reduction, lag (s) — TA = TB switches TGR off",
  KSTAB: "stabiliser gain", TW: "washout time constant (s)",
  T1: "lead-lag 1 lead (s)", T2: "lead-lag 1 lag (s)", T3: "lead-lag 2 lead (s)", T4: "lead-lag 2 lag (s)",
  Rt: "transformer resistance", Lt: "transformer reactance", RL_pu: "auxiliary load resistance",
  Rf: "filter resistance", Lf: "filter inductance", Cf: "filter capacitance",
  nq: "voltage/reactive-power droop", wf: "power-measurement filter cut-off (rad/s)",
  wc: "droop filter cut-off (rad/s)", eta: "dVOC synchronisation gain", alfa: "dVOC amplitude-regulation gain",
  J: "virtual inertia constant (s)", Dp: "virtual damping (pu)",
  K: "virtual excitation gain", Dq: "virtual voltage droop (pu)",
  K_theta: "DC-voltage to frequency gain (pu)",
  KpVL: "voltage loop proportional gain", KiVL: "voltage loop integral gain", Kffi: "current feed-forward gain",
  KpCL: "current loop proportional gain", KiCL: "current loop integral gain", Kffv: "voltage feed-forward gain",
  Cdc: "DC-link capacitance (pu)", Gdc: "DC-link conductance (pu)", Kpdc: "DC-voltage control gain", Tdc: "DC source time constant (s)",
  Kppll: "PLL proportional gain", Kipll: "PLL integral gain", Kpd: "DC-voltage loop proportional gain", Kid: "DC-voltage loop integral gain",
  Kiq: "reactive-power loop integral gain", wff: "feed-forward frequency (pu)", wb: "base angular frequency (rad/s) — set by the network frequency",
};

// --- Loop tuning maths (mirrors operating_point.gfm_params / gfl_params) ----
// "pi2": PI loop on a first-order plant K/(1 + tau s):   Kp = (2ζωn·tau − 1)/K, Ki = ωn²·tau/K
// "pi0": PI loop on an integrator (plant gain k):       Kp = 2ζωn/k,          Ki = ωn²/k
// "p1":  proportional loop on a first-order plant:       Kp·(1/G) = 3·tau/t_r − 1
// "i1":  integral-only loop:                             Ki = −3/t_r
// "droop": wf = 1/(2·mp·H)   (droop with emulated inertia H)
// with ωn = 3/(ζ·t_r) throughout, as in script_generic.m.
function loopDefs(unitType, p) {
  const wb = p.wb;
  if (unitType === "gfm") return [
    { id: "cl", deps: ["Rf", "Lf"], name: "Current loop (inner)", kind: "pi2", kp: "KpCL", ki: "KiCL", K: 1 / p.Rf, tau: p.Lf / (wb * p.Rf), plant: `K = 1/Rf = ${fmtSmart(1 / p.Rf)}, τ = Lf/(ωb·Rf) = ${fmtSmart(p.Lf / (wb * p.Rf))} s` },
    { id: "vl", deps: ["Cf"], name: "Voltage loop (outer)", kind: "pi0", kp: "KpVL", ki: "KiVL", k: wb / p.Cf, plant: `k = ωb/Cf = ${fmtSmart(wb / p.Cf)}` },
    { id: "dc", deps: ["Cdc", "Gdc"], name: "DC-link voltage loop", kind: "p1", kp: "Kpdc", G: p.Gdc, tau: p.Cdc / (wb * p.Gdc), plant: `τ = Cdc/(ωb·Gdc) = ${fmtSmart(p.Cdc / (wb * p.Gdc))} s` },
    { id: "droop", deps: [], name: "Droop & emulated inertia (power loop)", kind: "droop", mp: "mp", wf: "wf", plant: "wf = 1/(2·mp·H)" },
  ];
  if (unitType === "gfl") return [
    { id: "cl", deps: ["Rf", "Lf"], name: "Current loop (inner)", kind: "pi2", kp: "KpCL", ki: "KiCL", K: 1 / p.Rf, tau: p.Lf / (wb * p.Rf), plant: `K = 1/Rf = ${fmtSmart(1 / p.Rf)}, τ = Lf/(ωb·Rf) = ${fmtSmart(p.Lf / (wb * p.Rf))} s` },
    { id: "pll", deps: [], name: "Phase-locked loop", kind: "pi0", kp: "Kppll", ki: "Kipll", k: wb, plant: `k = ωb = ${fmtSmart(wb)}` },
    { id: "dcv", deps: ["Cdc", "Gdc"], name: "DC-voltage loop (outer)", kind: "pi2", kp: "Kpd", ki: "Kid", K: -1 / p.Gdc, tau: p.Cdc / (wb * p.Gdc), plant: `K = −1/Gdc = ${fmtSmart(-1 / p.Gdc)}, τ = Cdc/(ωb·Gdc) = ${fmtSmart(p.Cdc / (wb * p.Gdc))} s` },
    { id: "q", deps: [], name: "Reactive-power loop (outer)", kind: "i1", ki: "Kiq", plant: "Kiq = −3/t_r" },
  ];
  return [];
}

// t_r [s], ζ -> gains
function gainsFrom(loop, tr, z) {
  const wn = 3 / (z * tr);
  if (loop.kind === "pi2") return { kp: (2 * z * wn * loop.tau - 1) / loop.K, ki: (wn * wn * loop.tau) / loop.K, wn };
  if (loop.kind === "pi0") return { kp: (2 * z * wn) / loop.k, ki: (wn * wn) / loop.k, wn };
  if (loop.kind === "p1") return { kp: ((3 * loop.tau) / tr - 1) * loop.G };
  if (loop.kind === "i1") return { ki: -3 / tr };
  return {};
}
// gains -> t_r [s], ζ (null when the gains don't correspond to a real tuning)
function tuningFrom(loop, kp, ki) {
  if (loop.kind === "pi2") {
    const wn2 = (ki * loop.K) / loop.tau;
    if (!(wn2 > 0)) return null;
    const wn = Math.sqrt(wn2), z = (kp * loop.K + 1) / (2 * wn * loop.tau);
    return z > 0 ? { tr: 3 / (z * wn), z, wn } : null;
  }
  if (loop.kind === "pi0") {
    const wn2 = ki * loop.k;
    if (!(wn2 > 0)) return null;
    const wn = Math.sqrt(wn2), z = (kp * loop.k) / (2 * wn);
    return z > 0 ? { tr: 3 / (z * wn), z, wn } : null;
  }
  if (loop.kind === "p1") { const d = kp / loop.G + 1; return d > 0 ? { tr: (3 * loop.tau) / d } : null; }
  if (loop.kind === "i1") return ki < 0 ? { tr: -3 / ki } : null;
  return null;
}

// A unit's Rt/Lt: its own transformer (the one the power flow uses),
// converted from the transformer's rating to the system base -- or, in the
// MATLAB-compatible mode, the network's first transformer for every unit.
// Mirrors g2elin_core.operating_point.unit_transformer_rx.
function ownTransformer(der) { return state.network.transformers.find(t => t.lv_bus === der.bus) || null; }
function unitTransformerRx(der) {
  const net = state.network, own = ownTransformer(der);
  if (net.units_use_first_transformer || !own) {
    const first = net.transformers[0];
    return first ? { rt: first.r_pu, lt: first.x_pu, tr: first, linked: false } : { rt: 0.0, lt: 0.05, tr: null, linked: false };
  }
  const k = net.sn_mva / own.sn_mva;
  return { rt: own.r_pu * k, lt: own.x_pu * k, tr: own, linked: true };
}
const LINKED_TR_FIELD = { Rt: "r_pu", Lt: "x_pu" };

// Plant values the tuning depends on (see loopDefs().deps).
const PLANT_PARAMS = ["Rf", "Lf", "Cf", "Cdc", "Gdc"];

const UnitParams = {
  _defaultsCache: new Map(),

  // A unit's default parameter set depends only on its type and base values
  // (network base power and frequency, its bus voltage, and -- by the models'
  // shared-transformer convention -- the network's first transformer), so it
  // is fetched from /api/units/defaults rather than derived from the whole
  // network: that works for a unit that isn't connected yet, and while the
  // network as a whole is still invalid.
  async defaultsFor(der) {
    const net = state.network;
    const bus = net.buses.find(b => b.id === der.bus);
    const { rt, lt, tr } = unitTransformerRx(der);
    const body = { unit_type: der.unit_type, sn_mva: net.sn_mva, f_hz: net.f_hz, un_kv: bus ? bus.vn_kv : 20, rt_pu: rt, lt_pu: lt };
    // A machine's parameter set depends on its regulators, so they are part
    // of what the defaults are fetched (and cached) for.
    if (der.unit_type === "sm") {
      body.exciter = der.exciter || "g2elin";
      body.pss = der.pss || "g2elin";
      body.governor = der.governor || "g2elin";
    } else if (der.unit_type === "gfm") {
      body.controller = der.controller || "droop";
    }
    const key = JSON.stringify(body);
    if (!this._defaultsCache.has(key)) this._defaultsCache.set(key, api("/api/units/defaults", jsonPost(body)).then(r => r.params));
    try { return { defaults: await this._defaultsCache.get(key), noTransformer: !tr }; }
    catch (e) { this._defaultsCache.delete(key); throw e; }
  },

  // The droop-equivalent tuning a converter's outer loop stands for,
  // whichever law it runs. Every law's gains are written in terms of these
  // (script_generic.m), which is what makes the laws comparable -- so they
  // can be set on a VSM or a matching converter that carries no parameter
  // by these names at all.
  tuningHtml(der, defaults) {
    const t = gfmTuning({ ...defaults, ...(der.params || {}) }, der.controller || "droop");
    if (!t) return "";
    const row = (key, label, help, value) => `<div class="urow" data-key="${key}">
        <label for="tune-${der.id}-${key}" title="${esc(help)}"><code>${esc(label)}</code><span>${esc(help)}</span></label>
        <input id="tune-${der.id}-${key}" type="number" step="any" data-tuning="${key}" value="${+value.toPrecision(8)}">
        <span></span><span></span></div>`;
    return `<div class="ugroup"><div class="ugroup-title">Equivalent tuning</div>
      <p class="muted" style="font-size:0.74rem;margin:0 0 0.3rem">Every control law is tuned from these, so setting one here rewrites whichever law this unit runs \u2014 and swapping laws keeps them. H and \u03c9f say the same thing twice (H = 1/(2\u00b7mp\u00b7\u03c9f)); setting H moves \u03c9f.</p>
      ${row("mp", "mp", "active power / frequency droop (pu)", t.mp)}
      ${row("nq", "nq", "reactive power / voltage droop (pu)", t.nq)}
      ${row("wf", "\u03c9f", "power-measurement filter cut-off (rad/s)", t.wf)}
      ${row("H", "H", "equivalent inertia constant (s)", t.H)}
      <div class="notice warn-bg" data-role="tune-err" style="display:none;font-size:0.76rem"></div></div>`;
  },

  // Re-express a converter's outer loop: another law, another tuning, or
  // both. The formulas are the backend's (operating_point.gfm_retuned) and
  // are deliberately not mirrored here -- there are five laws and they all
  // have to agree. Only what differs from the new defaults is stored, so
  // der.params stays a set of overrides rather than a full copy.
  async retune(box, der, { to_controller = null, tuning = null } = {}) {
    const from = der.controller || "droop";
    const current = { ...box._defaults, ...(der.params || {}) };
    const r = await api("/api/units/retune", jsonPost({
      controller: from, params: current, to_controller, tuning,
    }));
    if (to_controller) der.controller = to_controller;
    const { defaults } = await this.defaultsFor(der);
    der.params = Object.fromEntries(Object.entries(r.params).filter(([k, v]) =>
      !(k in defaults) || Math.abs(v - defaults[k]) > 1e-9 * Math.max(1, Math.abs(v))));
    box.dataset.key = "";
    networkChanged();
    await this.render(box, der);
  },

  // box: container; der: the DerUnit (mutated in place)
  async render(box, der) {
    let got;
    try { got = await this.defaultsFor(der); }
    catch (e) { box.innerHTML = errorHtml(e); box.dataset.key = ""; return; }
    // The inspector may have moved to another element while this loaded.
    if (!document.body.contains(box)) return;
    const { defaults, noTransformer } = got;
    if (!Object.keys(defaults).length) { box.innerHTML = ""; box.dataset.key = ""; return; }
    const link = unitTransformerRx(der);
    // Don't rebuild (and lose focus/typing) when nothing structural changed.
    const { Rt: _rt, Lt: _lt, ...rest } = defaults;
    const key = `${der.id}|${der.unit_type}|${link.linked}|${der.exciter || ""}|${der.pss || ""}|${der.governor || ""}|${der.controller || ""}|${JSON.stringify(rest)}`;
    if (box.dataset.key === key && box.querySelector(".uparams")) {
      Object.assign(box._defaults, defaults);  // Rt/Lt follow the transformer
      this.refreshMarks(box, der, box._defaults);
      return;
    }
    box.dataset.key = key;
    box._defaults = defaults;
    der.params ||= {};
    const groups = [...(PARAM_GROUPS[der.unit_type] || []), ...regulatorGroups(der)];
    const listed = new Set(groups.flatMap(g => g[1]));
    const extra = Object.keys(defaults).filter(k => !listed.has(k) && k !== "wb");
    const allGroups = [...groups, ...(extra.length ? [["Other", extra]] : [])];
    const n = Object.keys(der.params).length;
    box.innerHTML = `${loopDefs(der.unit_type, defaults).length ? `<details class="uparams-wrap" open style="margin-top:0.9rem"><summary class="uparams-summary">Control loop tuner</summary><div class="tuner"></div></details>` : ""}
      <details class="uparams-wrap" open style="margin-top:0.7rem">
        <summary class="uparams-summary">Control & electrical parameters <span class="badge${n ? " warn" : ""}" data-role="count">${n ? `${n} changed` : "defaults"}</span></summary>
        <div class="uparams">
          <p class="muted" style="font-size:0.74rem;margin:0.5rem 0 0.2rem">Values are per unit unless stated. Changed values are highlighted; ↺ restores the default. Base frequency ωb = ${fmtSmart(defaults.wb)} rad/s (set by the network).${link.linked ? ` Rt/Lt are this unit's transformer (${esc(link.tr.name || `${link.tr.hv_bus} → ${link.tr.lv_bus}`)}), converted to the system base: editing them edits the transformer, so power flow and dynamics stay consistent.`
     : state.network.units_use_first_transformer ? " MATLAB-compatible mode: Rt/Lt follow the network's first transformer for every unit."
     : noTransformer ? " Rt/Lt will follow this unit's transformer; there is none yet, so they show placeholder defaults."
     : " This unit has no transformer of its own yet; Rt/Lt show the network's first transformer until it has one."}</p>
          <div class="notice warn-bg" data-role="plant-note" style="display:none;font-size:0.76rem;margin:0.4rem 0"></div>
          ${der.unit_type === "gfm" ? this.tuningHtml(der, defaults) : ""}
          ${unitModelSlots(der).length ? `<div class="ugroup uregulators"><div class="ugroup-title">${der.unit_type === "sm" ? "Regulator models" : "Control law"}</div>
            ${unitModelSlots(der).map(([slot, title, models]) => {
              const chosen = der[slot] || (models.g2elin ? "g2elin" : "droop");
              const help = "Each model brings its own parameters, listed below. Choosing \u2018none\u2019 removes the equipment altogether \u2014 its states and its parameters go with it, rather than being left idle.";
              return `<div class="ureg-row">
                <label for="ureg-${der.id}-${slot}" title="${esc(help)}"><code>${esc(title)}</code></label>
                <select id="ureg-${der.id}-${slot}" data-regulator="${slot}" title="${esc(help)}">
                  ${Object.entries(models).map(([id, [label]]) =>
                    `<option value="${id}"${id === chosen ? " selected" : ""}>${esc(label)}</option>`).join("")}
                </select></div>`;
            }).join("")}</div>` : ""}
          ${allGroups.map(([title, keys]) => `<div class="ugroup"><div class="ugroup-title">${esc(title)}</div>
            ${keys.filter(k => k in defaults).map(k => `
              <div class="urow" data-key="${k}">
                <label for="up-${der.id}-${k}" title="${esc(PARAM_HELP[k] || k)}"><code>${esc(k)}</code><span>${esc(PARAM_HELP[k] || "")}</span></label>
                <input id="up-${der.id}-${k}" type="number" step="any" data-param="${k}"${link.linked && LINKED_TR_FIELD[k] ? ` data-linked="${LINKED_TR_FIELD[k]}" title="Linked to the unit's transformer"` : ""} value="${+(der.params[k] ?? defaults[k]).toPrecision(8)}">
                ${(() => { const si = unitParamSi(der, k); return si
                  ? `<span class="usi"><input type="number" step="any" data-param-si="${k}" title="SI equivalent" value="${siFmt((der.params[k] ?? defaults[k]) * si.factor)}"><span class="u">${si.unit}</span></span>`
                  : "<span></span>"; })()}
                <button type="button" class="ghost small ureset" data-reset="${k}" title="Restore default (${fmtSmart(defaults[k])})">↺</button>
              </div>`).join("")}</div>`).join("")}
          <div class="controls" style="margin-top:0.6rem"><button type="button" class="ghost small" data-role="reset-all">Restore all defaults</button></div>
        </div>
      </details>`;

    // Swapping a regulator changes which parameters exist, so anything
    // overridden under a name the new model does not have is dropped rather
    // than left behind to be rejected by the backend. The panel is rebuilt
    // from scratch: its key includes the models.
    box.querySelectorAll("select[data-regulator]").forEach(sel => sel.addEventListener("change", async () => {
      const slot = sel.dataset.regulator;
      // A converter's laws are all written in terms of the same droop
      // tuning, so swapping one carries that tuning across rather than
      // dropping the unit back to the defaults. A machine's regulators have
      // no such correspondence -- a rate-feedback gain is not a transient
      // gain reduction -- so theirs are simply discarded.
      if (der.unit_type === "gfm" && slot === "controller") {
        try { await this.retune(box, der, { to_controller: sel.value }); return; }
        catch (e) { /* fall through to the plain swap below */ }
      }
      der[slot] = sel.value;
      const kept = new Set(regulatorGroups(der).flatMap(g => g[1]));
      const ownedByAny = new Set(unitModelSlots(der).flatMap(([, , models]) =>
        Object.values(models).flatMap(([, keys]) => keys)));
      Object.keys(der.params || {}).forEach(k => {
        if (ownedByAny.has(k) && !kept.has(k)) delete der.params[k];
      });
      box.dataset.key = "";
      networkChanged();
      this.render(box, der);
    }));

    box.querySelectorAll("input[data-tuning]").forEach(inp => inp.addEventListener("change", async () => {
      const value = parseFloat(inp.value);
      if (!Number.isFinite(value) || value <= 0) return this.render(box, der);
      try { await this.retune(box, der, { tuning: { [inp.dataset.tuning]: value } }); }
      catch (e) {
        const slot = box.querySelector('[data-role="tune-err"]');
        if (slot) { slot.textContent = String(e.message || e); slot.style.display = ""; }
      }
    }));

    // SI value typed -> per-unit value (the stored one) -> the usual path.
    box.querySelectorAll("input[data-param-si]").forEach(inp => inp.addEventListener("input", () => {
      const k = inp.dataset.paramSi, si = unitParamSi(der, k), v = parseFloat(inp.value);
      if (!si || !Number.isFinite(v)) return;
      const pu = box.querySelector(`input[data-param="${k}"]`);
      pu.value = +(v / si.factor).toPrecision(10);
      pu.dispatchEvent(new Event("input"));
    }));
    box.querySelectorAll("input[data-param]").forEach(inp => inp.addEventListener("input", () => {
      const v = parseFloat(inp.value);
      if (!Number.isFinite(v)) return;
      if (inp.dataset.linked) {
        // Rt/Lt linked to the transformer: edit the transformer itself (in
        // per unit of its own rating), never a separate override.
        const tr = ownTransformer(der);
        if (!tr || !(v > 0 || (inp.dataset.linked === "r_pu" && v >= 0))) return;
        tr[inp.dataset.linked] = v * tr.sn_mva / state.network.sn_mva;
        delete der.params[inp.dataset.param];
        defaults[inp.dataset.param] = v;
        networkChanged();
        this.refreshMarks(box, der, defaults);
        return;
      }
      this.setParam(der, defaults, inp.dataset.param, v);
      this.refreshMarks(box, der, defaults);
      this.renderTuner(box, der, defaults);
    }));
    box.querySelectorAll("[data-reset]").forEach(btn => btn.addEventListener("click", () => {
      const k = btn.dataset.reset;
      this.setParam(der, defaults, k, defaults[k]);
      box.querySelector(`input[data-param="${k}"]`).value = +defaults[k].toPrecision(8);
      this.refreshMarks(box, der, defaults);
      this.renderTuner(box, der, defaults);
    }));
    box.querySelector('[data-role="reset-all"]').addEventListener("click", () => {
      if (!Object.keys(der.params).length) return;
      der.params = {};
      networkChanged();
      box.querySelectorAll("input[data-param]").forEach(inp => { inp.value = +defaults[inp.dataset.param].toPrecision(8); });
      this.refreshMarks(box, der, defaults);
      this.renderTuner(box, der, defaults);
    });
    this.refreshMarks(box, der, defaults);
    this.renderTuner(box, der, defaults);
  },

  effective(der, defaults) { return { ...defaults, ...(der.params || {}) }; },

  // Store only real changes: a value equal to its default removes the override.
  setParam(der, defaults, k, v) {
    der.params ||= {};
    const same = Math.abs(v - defaults[k]) <= 1e-12 * Math.max(1, Math.abs(defaults[k]));
    if (same) delete der.params[k]; else der.params[k] = v;
    networkChanged();
  },

  refreshMarks(box, der, defaults) {
    box.querySelectorAll("input[data-linked]").forEach(inp => {
      if (inp !== document.activeElement) inp.value = +(der.params?.[inp.dataset.param] ?? defaults[inp.dataset.param]).toPrecision(8);
    });
    // SI values follow the per-unit ones (and the base values: bus voltage,
    // base power, frequency).
    box.querySelectorAll("input[data-param-si]").forEach(inp => {
      if (inp === document.activeElement) return;
      const k = inp.dataset.paramSi, si = unitParamSi(der, k);
      inp.value = si ? siFmt((der.params?.[k] ?? defaults[k]) * si.factor) : "";
    });
    box.querySelectorAll(".urow").forEach(row => {
      const changed = row.dataset.key in (der.params || {});
      row.classList.toggle("changed", changed);
    });
    const changedPlant = PLANT_PARAMS.filter(k => k in (der.params || {}));
    const note = box.querySelector('[data-role="plant-note"]');
    if (note) {
      note.style.display = changedPlant.length ? "" : "none";
      note.innerHTML = changedPlant.length ? `<span><b>${changedPlant.join(", ")} changed.</b> The controller gains are not re-derived automatically — they're still the ones tuned for the default ${changedPlant.length > 1 ? "values" : "value"}. The loop tuner above flags the affected loops and can re-apply the default tuning to the new plant.</span>` : "";
    }
    const n = Object.keys(der.params || {}).length;
    const c = box.querySelector('[data-role="count"]');
    if (c) { c.textContent = n ? `${n} changed` : "defaults"; c.classList.toggle("warn", !!n); }
  },

  renderTuner(box, der, defaults) {
    const host = box.querySelector(".tuner");
    if (!host) return;
    const p = this.effective(der, defaults);
    const loops = loopDefs(der.unit_type, p);
    // Keep a loop's half-typed state if the user is inside it.
    const active = document.activeElement && host.contains(document.activeElement) ? document.activeElement.closest(".loop")?.dataset.loop : null;
    host.innerHTML = `<p class="muted" style="font-size:0.74rem;margin:0.5rem 0 0.6rem">Type a response time t_r and damping ζ to get Kp/Ki — or type Kp/Ki to see the resulting t_r and ζ. ωn = 3/(ζ·t_r), as in the default tuning. “Apply” writes the gains to this unit.</p>`
      + loops.map(l => this.loopCard(l, p, this.staleTuning(l, der, defaults))).join("");
    loops.forEach(l => this.bindLoop(host.querySelector(`.loop[data-loop="${l.id}"]`), l, der, defaults, box));
    host.querySelectorAll("[data-retune]").forEach(btn => btn.addEventListener("click", () => {
      const l = loops.find(x => x.id === btn.dataset.retune);
      const st = this.staleTuning(l, der, defaults);
      if (!st) return;
      Object.entries(st.gains).forEach(([k, v]) => {
        this.setParam(der, defaults, k, v);
        const inp = box.querySelector(`input[data-param="${k}"]`);
        if (inp) inp.value = +v.toPrecision(8);
      });
      this.refreshMarks(box, der, defaults);
      this.renderTuner(box, der, defaults);
    }));
    if (active) host.querySelector(`.loop[data-loop="${active}"] input`)?.focus();
  },

  // A loop whose plant values (l.deps) were changed while its gains were
  // left at their defaults: those gains were tuned for the default plant.
  // Returns the changed plant values and the gains that would give the
  // default tuning (same t_r and ζ) on the new plant, or null.
  staleTuning(l, der, defaults) {
    const params = der.params || {};
    const changed = (l.deps || []).filter(k => k in params);
    const gainKeys = [l.kp, l.ki].filter(Boolean);
    if (!changed.length || gainKeys.some(k => k in params)) return null;
    const defLoop = loopDefs(der.unit_type, defaults).find(x => x.id === l.id);
    const t = tuningFrom(defLoop, defaults[l.kp], defaults[l.ki]);
    if (!t) return { changed, gains: {} };
    const g = gainsFrom(l, t.tr, t.z || 1);
    const gains = {};
    if (l.kp && "kp" in g) gains[l.kp] = g.kp;
    if (l.ki && "ki" in g) gains[l.ki] = g.ki;
    return { changed, gains, tr: t.tr, z: t.z };
  },

  loopCard(l, p, stale = null) {
    const num = (role, label, value, unit = "", title = "") => `<div class="field"><label title="${esc(title)}">${label}${unit ? ` <span class="unit">(${unit})</span>` : ""}</label><input type="number" step="any" data-role="${role}" value="${value}"></div>`;
    const r = v => (Number.isFinite(v) ? +v.toPrecision(6) : "");
    let fields = "";
    if (l.kind === "pi2" || l.kind === "pi0") {
      const t = tuningFrom(l, p[l.kp], p[l.ki]);
      fields = `<div class="loop-grid">
        ${num("tr", "t_r", t ? r(t.tr * 1e3) : "", "ms", "response time")}${num("z", "ζ", t ? r(t.z) : "", "", "damping ratio")}
        <span class="loop-arrow">⇄</span>
        ${num("kp", l.kp, r(p[l.kp]))}${num("ki", l.ki, r(p[l.ki]))}</div>
        <p class="loop-info" data-role="info"></p>`;
    } else if (l.kind === "p1") {
      const t = tuningFrom(l, p[l.kp]);
      fields = `<div class="loop-grid">${num("tr", "t_r", t ? r(t.tr * 1e3) : "", "ms", "response time")}<span class="loop-arrow">⇄</span>${num("kp", l.kp, r(p[l.kp]))}</div><p class="loop-info" data-role="info"></p>`;
    } else if (l.kind === "i1") {
      const t = tuningFrom(l, 0, p[l.ki]);
      fields = `<div class="loop-grid">${num("tr", "t_r", t ? r(t.tr * 1e3) : "", "ms", "response time")}<span class="loop-arrow">⇄</span>${num("ki", l.ki, r(p[l.ki]))}</div><p class="loop-info" data-role="info"></p>`;
    } else if (l.kind === "droop") {
      const H = 1 / (2 * p.mp * p.wf);
      fields = `<div class="loop-grid">${num("mp", "mp (droop)", r(p.mp), "pu")}${num("H", "H (inertia)", r(H), "s")}<span class="loop-arrow">⇄</span>${num("wf", "wf", r(p.wf), "rad/s")}</div><p class="loop-info" data-role="info"></p>`;
    }
    return `<div class="loop" data-loop="${l.id}"><div class="loop-head"><b>${esc(l.name)}</b><span class="spacer"></span>
      <button type="button" class="small" data-role="apply" disabled>Apply</button></div>
      <p class="loop-plant">${esc(l.plant)}</p>${stale ? `<div class="loop-stale"><b>${stale.changed.join(", ")} changed</b> — these gains are still the ones tuned for the default plant, so the loop no longer has its intended response (below: what it gives now).
        ${Object.keys(stale.gains).length ? `<button type="button" class="secondary small" data-retune="${l.id}">Re-tune to the default t_r = ${fmtSmart(stale.tr * 1e3)} ms${stale.z ? `, ζ = ${fmtSmart(stale.z)}` : ""}</button>` : ""}</div>` : ""}${fields}</div>`;
  },

  bindLoop(card, l, der, defaults, box) {
    const q = role => card.querySelector(`[data-role="${role}"]`);
    const val = role => parseFloat(q(role)?.value);
    const set = (role, v) => { if (q(role)) q(role).value = Number.isFinite(v) ? +v.toPrecision(6) : ""; };
    const current = () => this.effective(der, defaults);
    const pending = {};
    const info = msg => { q("info").innerHTML = msg; };
    const refreshApply = () => {
      const p = current();
      const differs = Object.entries(pending).some(([k, v]) => Number.isFinite(v) && Math.abs(v - p[k]) > 1e-12 * Math.max(1, Math.abs(p[k])));
      q("apply").disabled = !differs;
    };
    const describe = t => {
      if (!t) return `<span class="warn">These gains don't correspond to a stable, positively damped tuning of this loop.</span>`;
      return [t.wn ? `ωn = ${fmtSmart(t.wn)} rad/s (${fmtSmart(t.wn / (2 * Math.PI))} Hz)` : null, `t_r = ${fmtSmart(t.tr * 1e3)} ms`, t.z ? `ζ = ${fmtSmart(t.z)}` : null].filter(Boolean).join(" · ");
    };
    const fromTuning = () => {
      const tr = val("tr") / 1e3, z = l.kind === "pi2" || l.kind === "pi0" ? val("z") : 1;
      if (!(tr > 0) || !(z > 0)) { info(`<span class="warn">t_r and ζ must be positive.</span>`); return; }
      const g = gainsFrom(l, tr, z);
      if ("kp" in g) { set("kp", g.kp); pending[l.kp] = g.kp; }
      if ("ki" in g) { set("ki", g.ki); pending[l.ki] = g.ki; }
      const plantSign = l.kind === "pi2" ? Math.sign(l.K) : 1;
      const warn = ("kp" in g && g.kp * plantSign < 0) ? ` · <span class="warn">${l.kp} has the wrong sign: this t_r is slower than the plant itself</span>` : "";
      info(describe(tuningFrom(l, g.kp, g.ki) || { tr, z, wn: g.wn }) + warn);
      refreshApply();
    };
    const fromGains = () => {
      const kp = val("kp"), ki = val("ki");
      if (q("kp")) pending[l.kp] = kp;
      if (q("ki")) pending[l.ki] = ki;
      const t = tuningFrom(l, kp, ki);
      if (t) { set("tr", t.tr * 1e3); if (t.z) set("z", t.z); }
      info(describe(t));
      refreshApply();
    };
    if (l.kind === "droop") {
      const fromDroop = () => {
        const mp = val("mp"), H = val("H");
        if (!(mp > 0) || !(H > 0)) { info(`<span class="warn">mp and H must be positive.</span>`); return; }
        const wf = 1 / (2 * mp * H);
        set("wf", wf); pending.mp = mp; pending.wf = wf;
        info(`Power filter time constant 1/wf = ${fmtSmart(1 / wf)} s · droop gain 1/mp = ${fmtSmart(1 / mp)}`);
        refreshApply();
      };
      const fromWf = () => {
        const mp = val("mp"), wf = val("wf");
        if (!(mp > 0) || !(wf > 0)) { info(`<span class="warn">mp and wf must be positive.</span>`); return; }
        set("H", 1 / (2 * mp * wf)); pending.mp = mp; pending.wf = wf;
        info(`Emulated inertia H = ${fmtSmart(1 / (2 * mp * wf))} s · power filter time constant ${fmtSmart(1 / wf)} s`);
        refreshApply();
      };
      q("mp").addEventListener("input", fromDroop);
      q("H").addEventListener("input", fromDroop);
      q("wf").addEventListener("input", fromWf);
      const p = current();
      info(`Emulated inertia H = ${fmtSmart(1 / (2 * p.mp * p.wf))} s · power filter time constant ${fmtSmart(1 / p.wf)} s`);
    } else {
      ["tr", "z"].forEach(r => q(r)?.addEventListener("input", fromTuning));
      ["kp", "ki"].forEach(r => q(r)?.addEventListener("input", fromGains));
      const p = current();
      info(describe(tuningFrom(l, p[l.kp], p[l.ki])));
    }
    q("apply").addEventListener("click", () => {
      Object.entries(pending).forEach(([k, v]) => {
        if (!Number.isFinite(v)) return;
        this.setParam(der, defaults, k, v);
        const inp = box.querySelector(`input[data-param="${k}"]`);
        if (inp) inp.value = +v.toPrecision(8);
      });
      this.refreshMarks(box, der, defaults);
      this.renderTuner(box, der, defaults);
    });
  },
};
