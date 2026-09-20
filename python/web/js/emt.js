// EMT simulation page: integrates the nonlinear network model after a state
// offset, an input step, or a network event (breaker opening, load step,
// phase jump) at T0 = 0. Trajectories start T_PRE before T0 so
// the undisturbed x0 is visible; optionally the linearised model's response
// to the same disturbance is overlaid as dotted lines in the same colours.
//
// What is plotted is organised in *scopes*, like an oscilloscope/Simulink
// scope: each scope is one plot panel with its own mix of signals (states,
// inputs, outputs, measurements), on one set of axes or stacked subplots.
// The run requests every signal any scope uses.

const EMT_T_PRE = 0.01;

// Signal keys carry their kind, since names can repeat across kinds (e.g. an
// SM's "theta" is both a state and an output).
const SIGNAL_KINDS = {
  s: { title: "States", req: "plot_states", res: "series", lin: "series" },
  i: { title: "Inputs", req: "plot_inputs", res: "inputs", lin: "inputs" },
  o: { title: "Outputs", req: "plot_outputs", res: "outputs", lin: "outputs" },
  m: { title: "Measurements", req: "plot_measurements", res: "measurements", lin: null },
};
const sigKind = key => key.slice(0, key.indexOf(":"));
const sigName = key => key.slice(key.indexOf(":") + 1);

const EmtPage = {
  inited: false,
  names: null,
  scopes: [],
  measInfo: {},
  _scopeId: 0,

  init() {
    if (this.inited) return;
    this.inited = true;
    $("#page-emt").innerHTML = `
      <div class="page-head"><div class="crumb">Analysis</div><h1>EMT Simulation</h1>
        <p>Integrates the network's full nonlinear differential-algebraic equations (a coupled Newton solve at every step) after a disturbance at T0 = 0 — a state offset, an input step, or a network event (a breaker opening, a load step, a phase jump) — the ground truth that modal analysis linearises. Plots start ${EMT_T_PRE} s before T0 to show the initial operating point x0.</p></div>
      <div id="emt-context"></div>
      <div class="card">
        <div class="card-title">Disturbance at T0</div>
        <div class="btn-group" role="group" aria-label="Disturbance type" style="margin-bottom:0.9rem">
          <button class="secondary toggle small" data-dist="signal" aria-pressed="true">State offset / input step</button>
          <button class="secondary toggle small" data-dist="event" aria-pressed="false">Network event</button>
        </div>
        <div class="controls" id="emt-signal-box">
          <div class="field"><label for="emt-el-kind">Element type</label><select id="emt-el-kind"></select></div>
          <div class="field"><label for="emt-el">Element</label><select id="emt-el" style="min-width:230px"></select></div>
          <div class="field"><label for="emt-perturb">Perturbed state / stepped input</label><select id="emt-perturb" style="min-width:250px"></select></div>
          <div class="field"><label for="emt-amp">Amplitude</label><input type="number" id="emt-amp" step="any" value="0.02"></div>
        </div>
        <div id="emt-event-box" style="display:none">
          <div class="controls">
            <div class="field"><label for="emt-ev-kind">Event</label><select id="emt-ev-kind">
              <option value="breaker">Breaker opening (line / transformer / load / unit)</option>
              <option value="load_step">Load increase / decrease</option>
              <option value="phase_jump">Phase jump at a bus</option></select></div>
            <div class="field" data-ev="breaker"><label for="emt-ev-breaker">Element to disconnect</label><select id="emt-ev-breaker" style="min-width:300px"></select></div>
            <div class="field" data-ev="load_step"><label for="emt-ev-load">Load</label><select id="emt-ev-load" style="min-width:220px"></select></div>
            <div class="field" data-ev="load_step"><label for="emt-ev-dp">ΔP (% of the load)</label><input type="number" id="emt-ev-dp" step="any" value="20"></div>
            <div class="field" data-ev="load_step"><label for="emt-ev-dq">ΔQ (% of the load)</label><input type="number" id="emt-ev-dq" step="any" value="0"></div>
            <div class="field" data-ev="phase_jump"><label for="emt-ev-bus">Bus</label><select id="emt-ev-bus" style="min-width:240px"></select></div>
            <div class="field" data-ev="phase_jump"><label for="emt-ev-angle">Jump (degrees)</label><input type="number" id="emt-ev-angle" step="any" value="10"></div>
          </div>
          <p class="muted" id="emt-ev-note" style="font-size:0.76rem;margin:0.5rem 0 0"></p>
        </div>
        <div class="controls" style="margin-top:0.9rem">
          <div class="field"><label for="emt-tf">Duration after T0 (s, max 3)</label><input type="number" id="emt-tf" step="0.1" min="0.05" max="3" value="1.0"></div>
          <div class="field"><label for="emt-dt">Output timestep (s)</label><input type="number" id="emt-dt" step="any" placeholder="auto"></div>
        </div>
        <div class="card-title" style="margin-top:1.2rem">Scopes <span class="card-sub">— each scope is one plot panel; put any states, inputs, outputs and measurements in it</span></div>
        <div id="emt-scopes"></div>
        <div class="controls" style="margin-top:0.6rem;align-items:center">
          <button class="secondary small" id="emt-add-scope">+ Add scope</button>
          <label class="check" title="Buses, lines and loads have a lot of states (voltages and currents) that are rarely the point of a study. Their measurements stay available either way."><input type="checkbox" id="emt-hide-net" checked> Hide bus, line and load states</label>
          <span class="muted" id="emt-meas-note" style="font-size:0.76rem">Measurements: power flows (lines at both ends, loads, units / transformers), bus voltage magnitude, angle and frequency, instantaneous 3-phase bus voltages, and each unit's own frequency — in pu on the network base unless stated.</span>
        </div>
        <div class="controls" style="margin-top:1.1rem;align-items:center">
          <label class="check"><input type="checkbox" id="emt-linear"> Overlay the equivalent linearised response (dotted, same colours; states, inputs and outputs)</label>
          <label class="check"><input type="checkbox" id="emt-live" checked> Trace live <span class="muted" style="font-size:0.76rem">(streams solver steps; slower)</span></label>
        </div>
        <div class="controls" style="margin-top:1rem;align-items:center">
          <button id="emt-run"><svg viewBox="0 0 24 24" fill="currentColor"><path d="M8 5v14l11-7z"/></svg>Run EMT simulation</button>
          <button class="secondary" id="emt-stop" style="display:none">Stop</button>
          <span class="spinner" id="emt-spin"></span>
        </div>
      </div>
      <div id="emt-results" class="stack" style="margin-top:1rem"></div>`;
    $("#emt-run").addEventListener("click", () => this.run());
    $$("#page-emt [data-dist]").forEach(b => b.addEventListener("click", () => this.setDist(b.dataset.dist)));
    $("#emt-hide-net").addEventListener("change", () => this.refreshSignalOptions());
    $("#emt-el-kind").addEventListener("change", () => { this.elKey = ""; this.fillElementSelects(); this.fillPerturbSelect(); });
    $("#emt-el").addEventListener("change", e => { this.elKey = e.target.value; this.fillPerturbSelect(); });
    $("#emt-ev-kind").addEventListener("change", () => this.showEventFields());
    $("#emt-add-scope").addEventListener("click", () => { this.scopes.push(this.newScope()); this.renderScopes(); });
    on("network:loaded", () => {
      this.names = null;
      $("#emt-results").innerHTML = "";
      // Opened directly on this page, the network may arrive after it was shown.
      if ($("#page-emt").classList.contains("active")) this.onShow();
    });
    on("network:changed", () => { this.names = null; if (this.dist === "event") this.fillEventOptions(); });
  },

  dist: "signal",
  setDist(kind) {
    this.dist = kind;
    $$("#page-emt [data-dist]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.dist === kind)));
    $("#emt-signal-box").style.display = kind === "signal" ? "" : "none";
    $("#emt-event-box").style.display = kind === "event" ? "" : "none";
    if (kind === "event") this.fillEventOptions();
  },

  // Event targets: only elements in service now (open breakers already took
  // the others out), never the slack unit.
  fillEventOptions() {
    const net = state.network;
    if (!net) return;
    const sv = serviceState(net);
    const keep = (sel, html) => { const prev = sel.value; sel.innerHTML = html; if ([...sel.options].some(o => o.value === prev)) sel.value = prev; };
    const unitByBus = Object.fromEntries(net.der_units.map(d => [d.bus, d]));
    const lines = net.lines.map((l, i) => [l, i]).filter(([, i]) => sv.lines[i]);
    const trs = net.transformers.map((t, i) => [t, i]).filter(([t, i]) => sv.transformers[i] && unitByBus[t.lv_bus] && unitByBus[t.lv_bus].bus_type !== "slack" && sv.units[unitByBus[t.lv_bus].id]);
    const loads = net.loads.map((l, i) => [l, i]).filter(([, i]) => sv.loads[i]);
    const units = net.der_units.filter(d => d.bus_type !== "slack" && sv.units[d.id]);
    keep($("#emt-ev-breaker"),
      (lines.length ? `<optgroup label="Lines">${lines.map(([l, i]) => `<option value="line:${i}">Line #${i} (${l.from_bus} → ${l.to_bus})${l.name ? ` ${esc(l.name)}` : ""}</option>`).join("")}</optgroup>` : "")
      + (trs.length ? `<optgroup label="Unit transformers">${trs.map(([t, i]) => `<option value="transformer:${i}">Transformer #${i} (${t.hv_bus} → ${t.lv_bus}) — trips unit ${unitByBus[t.lv_bus].id}</option>`).join("")}</optgroup>` : "")
      + (loads.length ? `<optgroup label="Loads">${loads.map(([l, i]) => `<option value="load:${i}">Load #${i} at bus ${l.bus} (${fmtSmart(l.p_mw)} MW)</option>`).join("")}</optgroup>` : "")
      + (units.length ? `<optgroup label="Units">${units.map(d => `<option value="unit:${d.id}">${esc(UNIT_NAME[d.unit_type] || d.unit_type)} ${d.id} at bus ${d.bus} (${fmtSmart(d.p_set_mw)} MW)</option>`).join("")}</optgroup>` : ""));
    keep($("#emt-ev-load"), loads.map(([l, i]) => `<option value="${i}">Load #${i} at bus ${l.bus} (${fmtSmart(l.p_mw)} MW, ${fmtSmart(l.q_mvar)} MVAr)</option>`).join(""));
    // Phase-jump buses: the network buses (a unit's own terminal bus is inside
    // its model), plus the infinite bus's source.
    const slack = net.der_units.find(d => d.bus_type === "slack");
    const nodes = net.buses.filter(b => !unitByBus[b.id] && sv.energized.has(b.id));
    keep($("#emt-ev-bus"),
      (slack && slack.unit_type === "infinite_bus" ? `<option value="${slack.bus}">Infinite-bus source (bus ${slack.bus})</option>` : "")
      + nodes.map(b => `<option value="${b.id}">Bus ${b.id}${b.name ? ` (${esc(b.name)})` : ""}</option>`).join(""));
    this.showEventFields();
  },

  showEventFields() {
    const kind = $("#emt-ev-kind").value;
    $$("#emt-event-box [data-ev]").forEach(f => { f.style.display = f.dataset.ev === kind ? "" : "none"; });
    $("#emt-ev-note").textContent = {
      breaker: "The element is disconnected at T0 and the simulation continues with the rest of the network, from the pre-event state and with unchanged setpoints. Whatever the opening islands keeps running on its own (a load-only island decays). Opening a unit transformer trips its unit. The slack can't be tripped (it is the reference frame).",
      load_step: "The load's P and Q change by these percentages at T0 (as a constant impedance: at the operating-point voltage). To disconnect a load, use a breaker opening.",
      phase_jump: "At a network bus: the bus voltage phasor is rotated instantaneously (its shunt capacitor's state). At the infinite-bus source: the source's voltage angle jumps — the usual grid-code phase-jump test. A phase jump only moves the initial state, so the linearised model can reproduce it.",
    }[kind];
  },

  event() {
    const kind = $("#emt-ev-kind").value;
    if (kind === "breaker") {
      const v = $("#emt-ev-breaker").value;
      if (!v) return null;
      const [element, idx] = v.split(":");
      return { kind, element, index: +idx };
    }
    if (kind === "load_step") {
      const v = $("#emt-ev-load").value;
      return v === "" ? null : { kind, index: +v, dp_pct: parseFloat($("#emt-ev-dp").value) || 0, dq_pct: parseFloat($("#emt-ev-dq").value) || 0 };
    }
    const b = $("#emt-ev-bus").value;
    return b === "" ? null : { kind, bus: +b, angle_deg: parseFloat($("#emt-ev-angle").value) || 0 };
  },

  newScope(signals = []) {
    this._scopeId += 1;
    return { id: this._scopeId, title: `Scope ${this.scopes.length + 1}`, signals, stacked: false };
  },

  async onShow() {
    $("#emt-context").innerHTML = networkContextHtml({ plot: true });
    bindContextPlot($("#emt-context"));
    if (!state.network) return;
    if (this.names && this.names.version === state.version) return;
    const sel = $("#emt-perturb");
    sel.innerHTML = `<option>loading…</option>`; sel.disabled = true;
    const version = state.version;
    setSpinner($("#emt-spin"), "Building the nonlinear model…");
    try {
      const r = await netPost("states");
      if (version !== state.version) return;
      this.names = { version, data: r };
      if (this.dist === "event") this.fillEventOptions();
      this.measInfo = Object.fromEntries((r.measurements || []).map(m => [m.name, m]));
      // Every signal this model has, keyed "<kind>:<name>"; "raw" is the model
      // name, which says (see core.js) which element it belongs to.
      this.allOptions = [
        ...r.state_names.map(n => ({ name: `s:${n}`, raw: n, group: "States" })),
        ...r.input_names.map(n => ({ name: `i:${n}`, raw: n, group: "Inputs" })),
        ...r.output_names.map(n => ({ name: `o:${n}`, raw: n, group: "Outputs" })),
        ...(r.measurements || []).map(m => ({ name: `m:${m.name}`, raw: m.name, group: `Measurements · ${m.group}`, label: `${m.label} (${m.unit})` })),
      ];
      const valid = new Set(this.allOptions.map(o => o.name));
      this.scopes.forEach(sc => { sc.signals = sc.signals.filter(k => valid.has(k)); });
      if (!this.scopes.length || this.scopes.every(sc => !sc.signals.length)) {
        const dw = r.state_names.filter(n => n.includes("dw_r")).slice(0, 8).map(n => `s:${n}`);
        this.scopes = [{ ...this.newScope(dw), title: "Rotor speed deviations" }];
      }
      this.refreshSignalOptions();
    } catch (e) {
      sel.innerHTML = `<option>unavailable</option>`;
      $("#emt-results").innerHTML = errorHtml(e);
    } finally { setSpinner($("#emt-spin"), ""); }
  },

  // The states of buses, lines and loads are hidden unless asked for: there
  // are many of them (two per element) and a study rarely starts there. What
  // is already plotted is never hidden away.
  hideNetStates() { return $("#emt-hide-net")?.checked !== false; },

  refreshSignalOptions() {
    const all = this.allOptions || [];
    const kept = new Set(this.scopes.flatMap(sc => sc.signals));
    this.options = this.hideNetStates()
      ? all.filter(o => !(o.name.startsWith("s:") && isNetworkElementSignal(o.raw)) || kept.has(o.name))
      : all;
    this.fillElementSelects();
    this.fillPerturbSelect();
    this.renderScopes();
  },

  // Only elements that actually have a perturbable signal are offered.
  perturbOptions() {
    const net = this.hideNetStates();
    return (this.options || []).filter(o => (o.name.startsWith("s:") && !(net && isNetworkElementSignal(o.raw))) || o.name.startsWith("i:"));
  },

  fillElementSelects() {
    const keys = new Set();
    this.perturbOptions().forEach(o => { const e = signalElement(o.raw); if (e) keys.add(e.key); });
    const kindSel = $("#emt-el-kind"), elSel = $("#emt-el");
    if (!kindSel) return;
    this.elKind = kindSel.value || "";
    if (this.elKey && !keys.has(this.elKey)) this.elKey = "";
    kindSel.innerHTML = elementKindOptionsHtml(this.elKind, keys);
    this.elKind = kindSel.value;
    elSel.innerHTML = elementOptionsHtml(this.elKey, { kind: this.elKind, keys });
    this.elKey = elSel.value;
  },

  fillPerturbSelect() {
    const sel = $("#emt-perturb");
    if (!sel) return;
    const prev = sel.value;
    const match = o => {
      const e = signalElement(o.raw);
      if (this.elKey) return e && e.key === this.elKey;
      if (this.elKind) return e && e.kind === this.elKind;
      return true;
    };
    const opts = this.perturbOptions().filter(match);
    const group = (kind, label) => {
      const list = opts.filter(o => o.name.startsWith(kind));
      return list.length ? `<optgroup label="${label}">${list.map(o => `<option value="${esc(o.raw)}" data-kind="${kind === "s:" ? "state" : "input"}">${esc(o.raw)}</option>`).join("")}</optgroup>` : "";
    };
    sel.innerHTML = group("s:", "States (initial-condition offset)") + group("i:", "Inputs (step held from T0)");
    if (!sel.options.length) { sel.innerHTML = `<option value="">(no state or input here)</option>`; sel.disabled = true; return; }
    sel.disabled = false;
    const wanted = [...sel.options].some(o => o.value === prev) ? prev
      : ([...sel.options].find(o => o.value.includes("dw_r")) || sel.options[0]).value;
    sel.value = wanted;
  },

  // How a signal reads to a person: measurements by their description.
  label(key) {
    const k = sigKind(key), n = sigName(key);
    if (k === "m") { const m = this.measInfo[n]; return m ? `${n} — ${m.label} (${m.unit})` : n; }
    return k === "s" ? n : `${n} (${SIGNAL_KINDS[k].title.toLowerCase().slice(0, -1)})`;
  },

  renderScopes() {
    const box = $("#emt-scopes");
    box.innerHTML = "";
    this.scopes.forEach((sc, idx) => {
      const node = el(`<div class="scope-edit">
        <input type="text" class="scope-title" value="${esc(sc.title)}" aria-label="Scope title">
        <div class="scope-picker"></div>
        <label class="check" title="One subplot per signal instead of one shared axis"><input type="checkbox"${sc.stacked ? " checked" : ""}> Stacked</label>
        ${this.scopes.length > 1 ? `<button class="ghost small" title="Remove this scope">✕</button>` : ""}
      </div>`);
      node.querySelector(".scope-title").addEventListener("input", e => { sc.title = e.target.value; });
      node.querySelector("input[type=checkbox]").addEventListener("change", e => { sc.stacked = e.target.checked; });
      node.querySelector("button.ghost")?.addEventListener("click", () => { this.scopes.splice(idx, 1); this.renderScopes(); });
      new SignalPicker(node.querySelector(".scope-picker"), {
        options: this.options || [], selected: sc.signals, colors: true, placeholder: "add signal", byElement: true,
        display: key => sigName(key) + (sigKind(key) === "i" ? " (input)" : sigKind(key) === "o" ? " (output)" : ""),
        onChange: list => { sc.signals = list; this.checkTimestep(); },
      });
      box.appendChild(node);
    });
    this.checkTimestep();
  },

  // Every signal any scope uses, split by kind for the request.
  wanted() {
    const out = { plot_states: [], plot_inputs: [], plot_outputs: [], plot_measurements: [] };
    const seen = new Set();
    this.scopes.forEach(sc => sc.signals.forEach(key => {
      if (seen.has(key)) return;
      seen.add(key);
      out[SIGNAL_KINDS[sigKind(key)].req].push(sigName(key));
    }));
    return out;
  },

  request() {
    const sel = $("#emt-perturb");
    const opt = sel.selectedOptions[0];
    const dtRaw = $("#emt-dt").value.trim();
    const ev = this.dist === "event";
    return {
      perturb_kind: ev ? "event" : (opt ? opt.dataset.kind : "state"),
      perturb_name: ev ? "" : sel.value,
      perturb_offset: parseFloat($("#emt-amp").value) || 0,
      event: ev ? this.event() : null,
      t_final: parseFloat($("#emt-tf").value) || 1.0,
      dt: dtRaw ? parseFloat(dtRaw) : null,
      ...this.wanted(),
      t_pre: EMT_T_PRE,
      linear_overlay: $("#emt-linear").checked,
    };
  },

  // 3-phase waveforms need a fine timestep (a 50/60 Hz sine sampled every
  // few ms aliases): pick ~40 samples per cycle unless the user chose finer.
  checkTimestep() {
    const note = $("#emt-meas-note");
    const wave = this.scopes.some(sc => sc.signals.some(k => /^m:v_[abc]_/.test(k)));
    if (!wave || !state.network || !note) return;
    const f = state.network.f_hz, tf = parseFloat($("#emt-tf").value) || 1;
    const want = 1 / (40 * f), cur = parseFloat($("#emt-dt").value);
    if (!(cur > 0) || cur > 1 / (20 * f)) {
      const dt = Math.max(want, tf / 1999);
      $("#emt-dt").value = +dt.toPrecision(3);
      note.innerHTML = `<span class="warn">Timestep set to ${fmtSmart(dt * 1000)} ms</span> so the 3-phase waveforms are resolved (~${Math.round(1 / (f * dt))} samples per cycle)${dt > want * 1.01 ? " — limited by the 2000-sample maximum; shorten the duration for a finer one" : ""}.`;
    }
  },

  async run() {
    if (!this.names) { await this.onShow(); if (!this.names) return; }
    this.checkTimestep();
    const req = this.request();
    const out = $("#emt-results"), btn = $("#emt-run"), spin = $("#emt-spin");
    if (!this.scopes.some(sc => sc.signals.length)) { out.innerHTML = `<div class="notice warn-bg">Add at least one signal to a scope.</div>`; return; }
    if (req.perturb_kind === "event" && !req.event) { out.innerHTML = `<div class="notice warn-bg">Pick the element the event applies to.</div>`; return; }
    btn.disabled = true;
    out.innerHTML = "";
    try {
      if ($("#emt-live").checked) await this.runLive(req);
      else {
        setSpinner(spin, `Integrating${req.plot_inputs.length || req.plot_outputs.length || req.plot_measurements.length ? " (recovering inputs/outputs/measurements adds a Newton solve per sample)" : ""}${req.linear_overlay ? " + linearised response" : ""}…`);
        const r = await netPost("emt", req);
        this.renderResult(r, req);
      }
    } catch (e) {
      out.innerHTML = errorHtml(e);
    } finally { btn.disabled = false; setSpinner(spin, ""); }
  },

  renderResult(r, req, extraNote = "") {
    const out = $("#emt-results");
    const kind = r.perturb_kind === "input" ? "a step held from T0 in input" : "an initial-condition offset at T0 in state";
    const what = r.perturb_kind === "event"
      ? `Network event at T0: <b>${esc(r.perturbed || "")}</b>`
      : `Perturbed with ${kind} <b>${esc(r.perturbed)}</b> of ${fmt(req.perturb_offset, 4)}`;
    out.innerHTML = `<div class="card"><p class="status-line">${what} · ${fmt(req.t_final, 2)} s after T0${r.dt ? ` · dt = ${fmt(r.dt, 5)} s` : ""}${extraNote}</p>
      ${r.perturb_kind === "event" ? `<p class="muted" style="font-size:0.76rem;margin:0.4rem 0 0">Signals of an element the event removes stop at T0 (its power flows read 0).</p>` : ""}
      ${req.linear_overlay && r.linear_note ? `<div class="notice warn-bg" style="margin-top:0.6rem">${esc(r.linear_note)}</div>` : ""}
      ${r.linear ? `<div class="legend"><span><span class="line-swatch" style="border-color:var(--text-secondary)"></span>Nonlinear EMT</span><span><span class="line-swatch dotted" style="border-color:var(--text-secondary)"></span>Linearised model (same disturbance, same operating point)</span></div>
        <p class="muted" style="font-size:0.76rem;margin:0.5rem 0 0">Note: the nonlinear model's initial point is not an exact equilibrium (every node borrows line #1's susceptance), so the EMT traces include a small initial transient even without a disturbance. The linearised response has no such transient — part of the gap between the curves comes from that, not from nonlinearity.</p>` : ""}</div>`;
    this.scopes.forEach(sc => {
      // This scope's traces, in its own colour order, with the linear overlay where there is one.
      const traces = sc.signals.map((key, i) => {
        const kd = SIGNAL_KINDS[sigKind(key)], n = sigName(key);
        const y = (r[kd.res] || {})[n];
        const lin = kd.lin && r.linear ? (r.linear[kd.lin] || {})[n] : null;
        return y ? { key, label: this.label(key), color: seriesColor(i), y, lin } : null;
      }).filter(Boolean);
      if (!traces.length) return;
      const card = el(`<div class="card"><div class="card-title">${esc(sc.title || "Scope")} <span class="card-sub">(${traces.length} signal${traces.length > 1 ? "s" : ""})</span></div></div>`);
      const toSeries = tr => [{ name: tr.label, t: r.t, y: tr.y, color: tr.color }, ...(tr.lin ? [{ name: tr.label, t: r.linear.t, y: tr.lin, color: tr.color, dash: true }] : [])];
      if (sc.stacked && traces.length > 1) {
        const group = {}, wrap = el(`<div class="subplots"></div>`);
        traces.forEach(tr => wrap.appendChild(lineChart(toSeries(tr), { title: tr.label, height: 180, group, markers: [{ x: 0, label: "T0" }] })));
        card.appendChild(wrap);
      } else {
        card.appendChild(lineChart(traces.flatMap(toSeries), { height: 300, markers: [{ x: 0, label: "T0" }] }));
      }
      card.insertAdjacentHTML("beforeend", legendHtml(traces.map(tr => ({ name: tr.label, color: tr.color }))));
      out.appendChild(card);
    });
  },

  // Live tracing: NDJSON stream, charts redrawn at most ~8x per second.
  async runLive(req) {
    const out = $("#emt-results"), stop = $("#emt-stop");
    const res = await fetch("/api/network/emt/live", jsonPost({ network: state.network, ...req }));
    if (!res.ok) {
      const b = await res.json().catch(() => ({ detail: res.statusText }));
      throw new Error(formatDetail(b.detail) || `HTTP ${res.status}`);
    }
    const reader = res.body.getReader(), dec = new TextDecoder();
    stop.style.display = ""; stop.onclick = () => reader.cancel();
    const data = { t: [], series: {}, inputs: {}, outputs: {}, measurements: {} };
    let buf = "", done = null, error = null, lastDraw = 0;
    const draw = force => {
      const now = performance.now();
      if (!force && now - lastDraw < 120) return;
      lastDraw = now;
      this.renderResult({ ...data, perturb_kind: req.perturb_kind, perturbed: req.perturb_kind === "event" ? "(network event)" : req.perturb_name, dt: null }, req, ` · <span class="muted">tracing live — ${data.t.length} points</span>`);
    };
    setSpinner($("#emt-spin"), "Tracing live…");
    try {
      outer: while (true) {
        const { value, done: end } = await reader.read();
        if (end) break;
        buf += dec.decode(value, { stream: true });
        let i;
        while ((i = buf.indexOf("\n")) >= 0) {
          const line = buf.slice(0, i); buf = buf.slice(i + 1);
          if (!line.trim()) continue;
          const msg = JSON.parse(line);
          if (msg.error) { error = msg.error; break outer; }
          if (msg.done) { done = msg; break outer; }
          data.t.push(msg.t);
          for (const [grp, key] of [["series", "states"], ["inputs", "inputs"], ["outputs", "outputs"], ["measurements", "measurements"]]) {
            for (const [n, v] of Object.entries(msg[key] || {})) (data[grp][n] ||= []).push(v);
          }
          draw(false);
        }
      }
    } catch { if (!error) done = done || { stopped: true }; }
    finally { stop.style.display = "none"; }
    if (error) throw new Error(error);
    const note = done && !done.stopped ? ` · traced live, ${done.n_steps} solver steps` : ` · <span class="warn">stopped</span> — partial trace`;
    this.renderResult({ ...data, perturb_kind: req.perturb_kind, perturbed: (done && done.perturbed) || req.perturb_name, dt: null, linear: done && done.linear }, req, note);
    if (done && done.linear_error) out.insertAdjacentHTML("afterbegin", `<div class="notice warn-bg">Linearised overlay unavailable: ${esc(done.linear_error)}</div>`);
  },
};
