// Root locus (Modal Analysis sub-page): sweep any one network parameter over
// a range, re-solving the power flow and the closed-loop eigenvalues at every
// value (POST /api/network/modal/sweep, streamed), and draw each mode's locus
// coloured by the parameter value.

// Perceptually ordered low -> high (viridis-like), for the parameter value.
const LOCUS_STOPS = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
const locusColor = t => lerpStops(LOCUS_STOPS, t);

const SWEEP_FIELDS = {
  network: [["f_hz", "Frequency (Hz)"], ["sn_mva", "Base power (MVA)"]],
  bus: [["vn_kv", "Nominal voltage (kV)"]],
  line: [["r_pu", "Resistance R (pu)"], ["x_pu", "Reactance X (pu)"], ["b_pu", "Shunt susceptance B (pu)"], ["length_km", "Length (km) — no effect: impedances are per unit"]],
  transformer: [["r_pu", "Resistance R (pu)"], ["x_pu", "Reactance X (pu)"], ["sn_mva", "Rating (MVA)"]],
  load: [["p_mw", "Active power P (MW)"], ["q_mvar", "Reactive power Q (MVAr)"]],
  unit: [["p_set_mw", "Active power setpoint (MW)"], ["q_set_mvar", "Reactive power setpoint (MVAr)"], ["v_set_pu", "Voltage setpoint (pu)"],
    ["p_cons_mw", "Auxiliary load P (MW)"], ["q_cons_mvar", "Auxiliary load Q (MVAr)"], ["xd_pu", "Transient reactance Xd (pu) — SCR only, no effect on dynamics", ["sm"]]],
};
const ELEMENT_LABELS = { network: "Network", bus: "Bus", line: "Line", transformer: "Transformer", load: "Load", unit: "Unit" };

const RootLocus = {
  // Parameter rows: rows[0] is the main one (from/to/step -> the number of
  // values); every further row moves from its own "from" to "to" in lockstep
  // with it, so each step shows the combined effect of all of them.
  rows: [],
  result: null,      // {version, labels, values, extraValues, steps: [{i, value, ok, eig|error}], done}
  view: { axes: "symlog", onlyOsc: false, connect: true, baseline: true, track: null },
  controller: null,

  newRow(prev) {
    return { element: prev?.element ?? "unit", key: null, field: null, from: null, to: null, step: null, defaults: null };
  },

  render(body) {
    const net = state.network;
    if (this.result && this.result.version !== state.version) this.result = null;
    if (!this.rows.length) this.rows = [this.newRow()];
    body.innerHTML = `
      <div class="card"><div class="card-title">Parameters to sweep <span class="card-sub">— add several to see their combined effect: at every step each one is at the same fraction of its range</span></div>
        <div id="rl-rows"></div>
        <div class="controls" style="margin-top:0.8rem;align-items:center">
          <button class="secondary small" id="rl-add">+ Add parameter</button>
          <span class="status-line" id="rl-count"></span>
          <span style="flex:1"></span>
          <button id="rl-run">Run sweep</button><button class="secondary" id="rl-stop" style="display:none">Stop</button>
        </div>
        <div id="rl-progress" style="margin-top:0.8rem"></div>
      </div>
      <div class="card"><div class="card-title">Root locus <span class="card-sub" id="rl-sub"></span><span class="spacer"></span><button class="ghost small" id="rl-video" title="Record the sweep as a video, one frame per value">Video</button><button class="ghost small" id="rl-reset">Reset zoom</button></div>
        <div class="controls" style="margin-bottom:0.6rem;align-items:center">
          <div class="field"><label for="rl-axes">Axes</label><select id="rl-axes"><option value="symlog">Symlog (all modes)</option><option value="linear">Linear (zoomed region)</option></select></div>
          <div class="field rl-lin"><label for="rl-re0">Real from</label><input type="number" step="any" id="rl-re0"></div>
          <div class="field rl-lin"><label for="rl-re1">Real to</label><input type="number" step="any" id="rl-re1"></div>
          <div class="field rl-lin"><label for="rl-f1">Max |f| (Hz)</label><input type="number" step="any" id="rl-f1"></div>
          <label class="check"><input type="checkbox" id="rl-osc"> Oscillatory modes only</label>
          <label class="check"><input type="checkbox" id="rl-connect" checked> Connect loci</label>
          <label class="check"><input type="checkbox" id="rl-base" checked> Current eigenvalues</label>
        </div>
        <div class="rl-main" id="rl-main">
          <div id="rl-plot"><p class="empty">Choose a parameter and a range, then run the sweep.</p></div>
          <aside class="rl-side" id="rl-part" hidden></aside>
        </div>
        <div class="plot-legends" id="rl-legend"></div>
      </div>
      <div class="card"><div class="card-title">Most affected modes <span class="card-sub">— click a row to highlight its locus</span></div><div id="rl-table"><p class="empty">No sweep yet.</p></div></div>`;

    $("#rl-video").addEventListener("click", () => this.recordSweep());
    $("#rl-add").addEventListener("click", () => { this.rows.push(this.newRow(this.rows[this.rows.length - 1])); this.renderRows(); });
    $("#rl-run").addEventListener("click", () => this.run());
    $("#rl-stop").addEventListener("click", () => this.controller && this.controller.abort());
    $("#rl-axes").value = this.view.axes;
    $("#rl-osc").checked = this.view.onlyOsc; $("#rl-connect").checked = this.view.connect; $("#rl-base").checked = this.view.baseline;
    $("#rl-axes").addEventListener("change", () => { this.view.axes = $("#rl-axes").value; this.zoom = null; this.draw(); });
    $("#rl-osc").addEventListener("change", e => { this.view.onlyOsc = e.target.checked; this.draw(); this.drawTable(); });
    $("#rl-connect").addEventListener("change", e => { this.view.connect = e.target.checked; this.draw(); });
    $("#rl-base").addEventListener("change", e => { this.view.baseline = e.target.checked; this.draw(); });
    ["rl-re0", "rl-re1", "rl-f1"].forEach(id => $(`#${id}`).addEventListener("input", () => { this.zoom = null; this.draw(); }));
    $("#rl-reset").addEventListener("click", () => { this.zoom = null; this.draw(); });
    if (!net) return;
    this.renderRows();
    if (this.result) { this.draw(); this.drawTable(); this.renderProgress(); }
  },

  // --- Parameter rows ---
  renderRows() {
    const box = $("#rl-rows");
    box.innerHTML = "";
    this.rows.forEach((row, i) => {
      const main = i === 0;
      const node = el(`<div class="rl-row${main ? " main" : ""}">
        <span class="rl-tag">${main ? "Main" : `+ ${i}`}</span>
        <div class="field"><label>Element type</label><select data-r="el">${Object.entries(ELEMENT_LABELS).map(([k, v]) => `<option value="${k}">${v}</option>`).join("")}</select></div>
        <div class="field" data-r="key-field"><label>Element</label><select data-r="key" style="width:170px"></select></div>
        <div class="field"><label>Parameter</label><select data-r="field" style="width:280px"></select></div>
        <div class="field"><label>Current</label><span class="status-line" data-r="current" style="padding:0.5rem 0;min-width:3.5rem"></span></div>
        <div class="field"><label>From</label><input type="number" step="any" data-r="from"></div>
        <div class="field"><label>To</label><input type="number" step="any" data-r="to"></div>
        <div class="field"><label>${main ? "Step" : "Step (implied)"}</label>${main ? `<input type="number" step="any" data-r="step">` : `<span class="status-line" data-r="istep" style="padding:0.5rem 0;min-width:5.6rem;display:inline-block"></span>`}</div>
        ${main ? "" : `<button class="ghost small" data-r="remove" title="Remove this parameter" style="margin-bottom:0.3rem">✕</button>`}
      </div>`);
      const q = r => node.querySelector(`[data-r="${r}"]`);
      q("el").value = row.element;
      q("el").addEventListener("change", () => { row.element = q("el").value; row.key = null; row.field = null; this.fillKeys(row, node); });
      q("key").addEventListener("change", () => { row.key = +q("key").value; row.field = null; this.fillFields(row, node); });
      q("field").addEventListener("change", () => { row.field = q("field").value; this.showCurrent(row, node, true); });
      ["from", "to", "step"].forEach(k => q(k)?.addEventListener("input", () => { row[k] = parseFloat(q(k).value); this.updateCount(); }));
      q("remove")?.addEventListener("click", () => { this.rows.splice(i, 1); this.renderRows(); });
      box.appendChild(node);
      this.fillKeys(row, node, true);
    });
    this.updateCount();
  },

  keysFor(element) {
    const net = state.network;
    if (element === "network") return [];
    if (element === "bus") return net.buses.map(b => [b.id, `${b.id} — ${b.name || "bus"}`]);
    if (element === "line") return net.lines.map((l, i) => [i, `#${i}: ${l.from_bus} → ${l.to_bus}${l.name ? ` (${l.name})` : ""}`]);
    if (element === "transformer") return net.transformers.map((t, i) => [i, `#${i}: ${t.hv_bus} → ${t.lv_bus}${t.name ? ` (${t.name})` : ""}`]);
    if (element === "load") return net.loads.map((l, i) => [i, `#${i} at bus ${l.bus}${l.name ? ` (${l.name})` : ""}`]);
    if (element === "unit") return net.der_units.map(d => [d.id, `${UNIT_LABEL[d.unit_type] || d.unit_type} id ${d.id} (bus ${d.bus})`]);
    return [];
  },

  // keepRange: re-rendering existing rows keeps the ranges already typed.
  fillKeys(row, node, keepRange = false) {
    const q = r => node.querySelector(`[data-r="${r}"]`);
    const keys = this.keysFor(row.element);
    q("key-field").style.display = row.element === "network" ? "none" : "";
    q("key").innerHTML = keys.length ? keys.map(([k, lab]) => `<option value="${k}">${esc(lab)}</option>`).join("") : `<option value="">none in this network</option>`;
    if (row.key === null || !keys.some(([k]) => k === row.key)) row.key = keys.length ? keys[0][0] : null;
    if (row.key !== null) q("key").value = row.key;
    this.fillFields(row, node, keepRange);
  },

  element(row) {
    const net = state.network, { element, key } = row;
    if (element === "network") return net;
    if (element === "bus") return net.buses.find(b => b.id === key);
    if (element === "unit") return net.der_units.find(d => d.id === key);
    return ({ line: net.lines, transformer: net.transformers, load: net.loads }[element] || [])[key];
  },

  async fillFields(row, node, keepRange = false) {
    const sel = node.querySelector('[data-r="field"]');
    const elm = this.element(row);
    if (!elm) { sel.innerHTML = ""; this.showCurrent(row, node, !keepRange); return; }
    // A third entry restricts a field to certain unit types -- a transient
    // reactance is a machine's, and offering it on a converter would sweep
    // something that unit does not have.
    const offered = SWEEP_FIELDS[row.element]
      .filter(([, , only]) => !only || (row.element === "unit" && only.includes(elm.unit_type)));
    let html = `<optgroup label="${esc(ELEMENT_LABELS[row.element])}">${offered.map(([f, lab]) => `<option value="${f}">${esc(lab)}</option>`).join("")}</optgroup>`;
    row.defaults = null;
    if (row.element === "unit" && elm.unit_type !== "infinite_bus") {
      try {
        const { defaults } = await UnitParams.defaultsFor(elm);
        row.defaults = defaults;
        const groups = PARAM_GROUPS[elm.unit_type] || [];
        const loops = loopDefs(elm.unit_type, { ...defaults, ...(elm.params || {}) });
        if (loops.length) {
          const q = l => l.kind === "droop"
            ? [["H", "emulated inertia H (s)", "H (s)"], ["Tf_ms", "power filter time constant 1/wf (ms)", "T_f (ms)"]]
            : [["tr_ms", "response time t_r (ms)", "t_r (ms)"], ...(l.kind === "pi2" || l.kind === "pi0" ? [["zeta", "damping ζ", "ζ"]] : [])];
          html = `<optgroup label="Control loop tuning (Kp/Ki re-derived at each value)">${loops.flatMap(l => q(l).map(([id, what, short]) =>
            `<option value="tune.${l.id}.${id}" data-short="${esc(l.name)} ${esc(short)}">${esc(l.name)} — ${esc(what)}</option>`)).join("")}</optgroup>` + html;
        }
        html += groups.map(([title, keys]) => `<optgroup label="${esc(title)}">${keys.filter(k => k in defaults).map(k => `<option value="params.${k}">${esc(k)} — ${esc(PARAM_HELP[k] || "")}</option>`).join("")}</optgroup>`).join("");
      } catch { /* unit fields only */ }
    }
    sel.innerHTML = html;
    if (row.field && [...sel.options].some(o => o.value === row.field)) sel.value = row.field;
    else {
      const pref = { gfm: "tune.vl.tr_ms", gfl: "tune.pll.tr_ms", sm: "params.H" }[elm.unit_type];
      sel.value = pref && [...sel.options].some(o => o.value === pref) ? pref : sel.options[0]?.value;
      keepRange = false;
    }
    row.field = sel.value;
    this.showCurrent(row, node, !keepRange || row.from === null);
  },

  currentValue(row) {
    const elm = this.element(row), f = row.field;
    if (!elm || !f) return null;
    if (f.startsWith("tune.")) {
      // Same tuning maths as the loop tuner (unit_params.js / g2elin_core/tuning.py).
      const [, loopId, qty] = f.split(".");
      const p = { ...(row.defaults || {}), ...(elm.params || {}) };
      const l = loopDefs(elm.unit_type, p).find(x => x.id === loopId);
      if (!l) return null;
      if (l.kind === "droop") return qty === "H" ? 1 / (2 * p.mp * p.wf) : 1000 / p.wf;
      const t = tuningFrom(l, p[l.kp], p[l.ki]);
      if (!t) return null;
      return qty === "zeta" ? t.z : t.tr * 1000;
    }
    if (f.startsWith("params.")) { const k = f.slice(7); return elm.params && k in elm.params ? elm.params[k] : row.defaults?.[k] ?? null; }
    return elm[f] ?? null;
  },

  showCurrent(row, node, resetRange) {
    const q = r => node.querySelector(`[data-r="${r}"]`);
    const v = this.currentValue(row);
    q("current").textContent = v === null ? "–" : fmtSmart(v);
    if (resetRange && v !== null) {
      const lo = v === 0 ? 0 : v * 0.5, hi = v === 0 ? 1 : v * 1.5;
      const r = x => +x.toPrecision(4);
      row.from = r(Math.min(lo, hi)); row.to = r(Math.max(lo, hi));
      if (row === this.rows[0]) row.step = r(Math.abs(hi - lo) / 10);
    }
    if (row.from !== null) q("from").value = row.from;
    if (row.to !== null) q("to").value = row.to;
    if (q("step") && row.step !== null) q("step").value = row.step;
    this.updateCount();
  },

  range() {
    const m = this.rows[0];
    const a = m?.from, b = m?.to, s = Math.abs(m?.step);
    if (![a, b, s].every(Number.isFinite) || s === 0) return null;
    // Same count as the server's sweep_values(): whole steps, plus the end
    // point when the range isn't a whole number of steps.
    const whole = Math.floor(Math.abs(b - a) / s + 1e-9);
    const offGrid = Math.abs(whole * s - Math.abs(b - a)) > 1e-9 * Math.max(1, Math.abs(b));
    return { start: a, stop: b, step: s, n: whole + 1 + (offGrid ? 1 : 0) };
  },

  updateCount() {
    const r = this.range();
    const out = $("#rl-count");
    if (!out) return;
    // Implied steps of the extra rows (they share the main row's number of values).
    $$("#rl-rows .rl-row").forEach((node, i) => {
      const is = node.querySelector('[data-r="istep"]'), row = this.rows[i];
      if (is) is.textContent = r && r.n > 1 && Number.isFinite(row.from) && Number.isFinite(row.to) ? `≈ ${fmtSmart((row.to - row.from) / (r.n - 1))}` : "–";
    });
    if (!r) { out.innerHTML = `<span class="warn">main parameter: enter from, to and a non-zero step</span>`; return; }
    const n = state.network ? state.network.buses.length : 10;
    const extra = this.rows.length > 1 ? ` · ${this.rows.length} parameters together` : "";
    out.innerHTML = r.n > 201 ? `<span class="bad">${r.n} values — at most 201</span>` : `${r.n} values${extra} · roughly ${Math.max(1, Math.round(r.n * (0.3 + n * 0.08)))} s`;
  },

  rowLabel(row, node) {
    const opt = node.querySelector('[data-r="field"]').selectedOptions[0];
    const fieldLabel = opt?.dataset.short || opt?.textContent.split(" — ")[0] || row.field;
    const elLabel = row.element === "network" ? "Network" : node.querySelector('[data-r="key"]').selectedOptions[0]?.textContent;
    return `${elLabel} · ${fieldLabel}`;
  },

  // --- Run ---
  async run() {
    const r = this.range();
    if (!r || r.n > 201 || !this.rows[0].field) { this.updateCount(); return; }
    const nodes = $$("#rl-rows .rl-row");
    const target = row => ({ element: row.element, key: row.element === "network" ? null : row.key, field: row.field });
    const extras = this.rows.slice(1);
    if (extras.some(x => !x.field || !Number.isFinite(x.from) || !Number.isFinite(x.to))) {
      $("#rl-count").innerHTML = `<span class="bad">every added parameter needs a from and a to value</span>`;
      return;
    }
    const keys = this.rows.map(row => JSON.stringify(target(row)));
    if (new Set(keys).size !== keys.length) { $("#rl-count").innerHTML = `<span class="bad">the same parameter is listed twice</span>`; return; }
    const labels = this.rows.map((row, i) => this.rowLabel(row, nodes[i]));
    this.result = {
      version: state.version, labels, label: labels.length > 1 ? `${labels[0]} + ${labels.length - 1} more` : labels[0],
      values: [], extraValues: [], steps: [], done: false, stopped: false,
    };
    this.zoom = null;
    this.view.track = null;
    this.stopParticipation();
    const btn = $("#rl-run"), stop = $("#rl-stop");
    btn.disabled = true; stop.style.display = "";
    this.controller = new AbortController();
    this.renderProgress();
    try {
      const body = {
        network: state.network, target: target(this.rows[0]), start: r.start, stop: r.stop, step: r.step,
        extra: extras.map(x => ({ target: target(x), start: x.from, stop: x.to })),
      };
      const res = await fetch("/api/network/modal/sweep", { ...jsonPost(body), signal: this.controller.signal });
      if (!res.ok) {
        const b = await res.json().catch(() => ({ detail: res.statusText }));
        throw new Error(formatDetail(b.detail) || `HTTP ${res.status}`);
      }
      const reader = res.body.getReader(), dec = new TextDecoder();
      let buf = "", last = 0;
      while (true) {
        const { value, done } = await reader.read();
        if (done) break;
        buf += dec.decode(value, { stream: true });
        let i;
        while ((i = buf.indexOf("\n")) >= 0) {
          const line = buf.slice(0, i); buf = buf.slice(i + 1);
          if (!line.trim()) continue;
          const msg = JSON.parse(line);
          if (msg.start) { this.result.values = msg.values; this.result.extraValues = msg.extra_values || []; }
          else if (msg.state_names) { this.result.stateNames = msg.state_names; this.result.steps.push(msg); }
          else if (msg.done) this.result.done = true;
          else this.result.steps.push(msg);
        }
        if (performance.now() - last > 250) { last = performance.now(); this.redrawIfShown(); }
      }
    } catch (e) {
      if (e.name === "AbortError") this.result.stopped = true;
      else { this.result.error = e.message; }
    } finally {
      btn.disabled = false; stop.style.display = "none"; this.controller = null;
      this.redrawIfShown();
    }
  },

  // Every swept parameter's value at step i, as [label, value] pairs.
  stepValues(i) {
    const r = this.result;
    return r.labels.map((lab, k) => [lab, k === 0 ? r.values[i] : r.extraValues[k - 1]?.[i]]);
  },

  redrawIfShown() {
    if (!$("#rl-plot")) return;
    this.renderProgress(); this.draw(); this.drawTable(); this.drawParticipation();
  },

  renderProgress() {
    const box = $("#rl-progress"), r = this.result;
    if (!box || !r) return;
    const n = r.values.length || 1, k = r.steps.length, bad = r.steps.filter(s => !s.ok);
    const pct = Math.round((100 * k) / n);
    box.innerHTML = r.error ? errorHtml(r.error)
      : `<div class="rl-bar"><div style="width:${pct}%"></div></div>
        <p class="status-line" style="margin-top:0.35rem">${r.done ? '<span class="ok">done</span>' : r.stopped ? '<span class="warn">stopped</span>' : "sweeping…"} · ${k}/${r.values.length || "?"} values
        ${bad.length ? ` · <span class="warn">${bad.length} did not solve</span> <span class="muted">(${esc(bad[0].error.slice(0, 140))}${bad.length > 1 ? " …" : ""})</span>` : ""}</p>`;
    $("#rl-sub").textContent = r.label ? `— ${r.label}` : "";
  },

  // tracks[k] = [{value, re, im}] -- the k-th mode across the solved steps
  tracks() {
    // view.upTo limits the loci to the values solved so far -- how the video
    // is drawn, frame by frame (recordSweep).
    const upTo = this.view.upTo;
    const ok = (this.result?.steps || []).filter(s => s.ok && (upTo === undefined || s.i <= upTo));
    if (!ok.length) return [];
    const n = ok[0].eig.length;
    const tracks = Array.from({ length: n }, () => []);
    ok.forEach(s => { if (s.eig.length === n) s.eig.forEach(([re, im], k) => tracks[k].push({ value: s.value, i: s.i, re, im })); });
    return tracks.map((pts, k) => ({ k, pts })).filter(t => !this.view.onlyOsc || t.pts.some(p => Math.abs(p.im) > 1e-6));
  },

  // --- Plot ---
  draw() {
    const box = $("#rl-plot");
    if (!box || !this.result) return;
    const tracks = this.tracks();
    if (!tracks.length) { box.innerHTML = `<p class="empty">${this.result.done || this.result.stopped ? "No value solved — nothing to plot." : "Waiting for the first value…"}</p>`; return; }
    const vals = this.result.values.length ? this.result.values : this.result.steps.map(s => s.value);
    const vmin = Math.min(...vals), vmax = Math.max(...vals);
    const tOf = v => (vmax > vmin ? (v - vmin) / (vmax - vmin) : 0.5);
    const W = 900, H = 520, PAD = 56;
    const lin = this.view.axes === "linear";
    $$(".rl-lin").forEach(f => { f.style.display = lin ? "" : "none"; });
    const all = tracks.flatMap(t => t.pts);
    let sx, sy, grid = "";
    if (!lin) {
      const xMax = Math.max(1, ...all.map(p => Math.abs(symlog(p.re)))), yMax = Math.max(1, ...all.map(p => Math.abs(symlog(p.im / (2 * Math.PI)))));
      sx = re => PAD + ((symlog(re) + xMax) / (2 * xMax)) * (W - 2 * PAD);
      sy = im => H - PAD - ((symlog(im / (2 * Math.PI)) + yMax) / (2 * yMax)) * (H - 2 * PAD);
      grid = eigenGridLines(v => PAD + ((v + xMax) / (2 * xMax)) * (W - 2 * PAD), v => H - PAD - ((v + yMax) / (2 * yMax)) * (H - 2 * PAD), xMax, yMax, W, H, PAD)
        + dampingGuideLines(v => PAD + ((v + xMax) / (2 * xMax)) * (W - 2 * PAD), v => H - PAD - ((v + yMax) / (2 * yMax)) * (H - 2 * PAD), 0.05, "var(--critical)");
    } else {
      // Default region: the slow modes (|λ| < 300 rad/s), where root loci are usually read.
      if (!$("#rl-re0").value || !$("#rl-f1").value) {
        const slow = all.filter(p => Math.hypot(p.re, p.im) < 300);
        const src = slow.length ? slow : all;
        $("#rl-re0").value = +(Math.min(...src.map(p => p.re)) * 1.1 - 0.5).toPrecision(3);
        $("#rl-re1").value = +Math.max(1, Math.max(...src.map(p => p.re)) * 1.1 + 0.5).toPrecision(3);
        $("#rl-f1").value = +Math.max(1, Math.max(...src.map(p => Math.abs(p.im) / (2 * Math.PI))) * 1.15).toPrecision(3);
      }
      const re0 = +$("#rl-re0").value, re1 = +$("#rl-re1").value, f1 = Math.abs(+$("#rl-f1").value) || 1;
      sx = re => PAD + ((re - re0) / ((re1 - re0) || 1)) * (W - 2 * PAD);
      sy = im => H - PAD - ((im / (2 * Math.PI) + f1) / (2 * f1)) * (H - 2 * PAD);
      const xt = niceTicks(re0, re1, 8), yt = niceTicks(-f1, f1, 6);
      grid = xt.ticks.map(v => `<line class="eigen-gridline" x1="${sx(v)}" y1="${PAD}" x2="${sx(v)}" y2="${H - PAD}"/><text class="eigen-ticklabel" x="${sx(v)}" y="${H - PAD + 14}" text-anchor="middle">${tickLabel(v, xt.step)}</text>`).join("")
        + yt.ticks.map(v => `<line class="eigen-gridline" x1="${PAD}" y1="${sy(v * 2 * Math.PI)}" x2="${W - PAD}" y2="${sy(v * 2 * Math.PI)}"/><text class="eigen-ticklabel" x="${PAD - 6}" y="${sy(v * 2 * Math.PI) + 3}" text-anchor="end">${tickLabel(v, yt.step)}</text>`).join("");
      // 5 % damping line: |f| = slope·|re|/(2π)
      const slope = Math.sqrt(1 - 0.05 ** 2) / 0.05;
      const reEnd = Math.min(re0, -1e-9);
      grid += `<polyline points="${sx(0)},${sy(0)} ${sx(reEnd)},${sy(-slope * reEnd)}" fill="none" stroke="var(--critical)" stroke-dasharray="4,3"/>
        <polyline points="${sx(0)},${sy(0)} ${sx(reEnd)},${sy(slope * reEnd)}" fill="none" stroke="var(--critical)" stroke-dasharray="4,3"/>`;
    }
    const zero = `<line x1="${PAD}" y1="${sy(0)}" x2="${W - PAD}" y2="${sy(0)}" stroke="var(--baseline)"/><line x1="${sx(0)}" y1="${PAD}" x2="${sx(0)}" y2="${H - PAD}" stroke="var(--baseline)" stroke-width="1.5"/>`;
    const hl = this.view.track;
    let lines = "", dots = "";
    tracks.forEach(t => {
      const faded = hl !== null && t.k !== hl;
      if (this.view.connect && t.pts.length > 1) {
        lines += `<polyline points="${t.pts.map(p => `${sx(p.re).toFixed(1)},${sy(p.im).toFixed(1)}`).join(" ")}" fill="none" stroke="${t.k === hl ? "var(--text-primary)" : "#9a9890"}" stroke-width="${t.k === hl ? 2 : 1}" opacity="${faded ? 0.15 : 0.7}"/>`;
      }
      t.pts.forEach((p, j) => {
        dots += `<circle cx="${sx(p.re).toFixed(1)}" cy="${sy(p.im).toFixed(1)}" r="${t.k === hl ? 4.5 : 3}" fill="${locusColor(tOf(p.value))}" opacity="${faded ? 0.2 : 1}" data-k="${t.k}" data-j="${j}" class="rl-pt"/>`;
      });
    });
    let base = "";
    if (this.view.baseline && state.modal && state.modal.version === state.version) {
      base = state.modal.data.modes.filter(m => !this.view.onlyOsc || Math.abs(m.imag) > 1e-6)
        .map(m => `<circle cx="${sx(m.real).toFixed(1)}" cy="${sy(m.imag).toFixed(1)}" r="5" fill="none" stroke="var(--text-primary)" stroke-width="1.2" pointer-events="none"/>`).join("");
    }
    box.innerHTML = `<svg class="eigenmap-svg" viewBox="0 0 ${W} ${H}" style="width:100%;height:auto;display:block">
      <defs><clipPath id="rl-clip"><rect x="${PAD}" y="${PAD}" width="${W - 2 * PAD}" height="${H - 2 * PAD}"/></clipPath></defs>
      ${grid}${zero}<g clip-path="${lin ? "url(#rl-clip)" : ""}">${lines}${base}${dots}
        <circle id="rl-cursor" r="8" fill="none" stroke="var(--text-primary)" stroke-width="2.5" opacity="0" pointer-events="none"/></g>
      <text x="${W - PAD}" y="${H - 8}" text-anchor="end" class="eigen-ticklabel" style="font-size:10.5px">Real part (1/s${lin ? "" : ", symlog"})</text>
      <text x="${PAD}" y="${PAD - 10}" class="eigen-ticklabel" style="font-size:10.5px">Frequency, Hz${lin ? "" : " (symlog)"}</text></svg>`;
    // Keep the projection, so the participation panel can move its marker
    // along the locus without redrawing the whole plot at every frame.
    this._proj = { sx, sy };
    const svg = box.querySelector("svg");
    const prev = this.zoom ? this.zoom.get() : null;
    this.zoom = attachSvgZoomPan(svg);
    if (prev) this.zoom.set(prev);
    svg.addEventListener("mouseover", e => {
      const c = e.target.closest(".rl-pt");
      if (!c) return;
      const t = tracks.find(x => x.k === +c.dataset.k), p = t.pts[+c.dataset.j];
      const wn = Math.hypot(p.re, p.im), z = wn ? -p.re / wn : 0;
      showTooltip(`<div class="tt-title">Mode ${t.k}</div>` + ttTable([
        ...this.stepValues(p.i).map(([lab, v]) => [lab, fmtSmart(v)]), { sep: `step ${p.i + 1} of ${this.result.values.length}` }, ["λ", `${p.re.toExponential(3)} ${p.im >= 0 ? "+" : "−"} j${Math.abs(p.im).toExponential(3)}`],
        ["Frequency", `${fmt(Math.abs(p.im) / (2 * Math.PI), 3)} Hz`], ["Damping", `${fmt(100 * z, 2)} %`]]), e);
    });
    svg.addEventListener("mousemove", e => { if (e.target.closest(".rl-pt")) moveTooltip(e); });
    svg.addEventListener("mouseover", e => {
      // Hovering the tracked mode's own locus scrubs its participation bars.
      const c = e.target.closest(".rl-pt");
      if (!c || +c.dataset.k !== this.view.track || !this.part) return;
      this.pauseParticipation();
      this.setParticipationStep(+c.dataset.j);
    });
    svg.addEventListener("mouseout", e => { if (e.target.closest(".rl-pt")) hideTooltip(); });
    svg.addEventListener("click", e => {
      const c = e.target.closest(".rl-pt");
      if (svg.__justPanned) return;
      this.view.track = c ? (+c.dataset.k === this.view.track ? null : +c.dataset.k) : null;
      this.draw(); this.drawTable(); this.drawParticipation();
    });
    $("#rl-legend").innerHTML = gradientLegendHtml(this.result.labels[0], LOCUS_STOPS, vmin, vmax, 4)
      + (this.result.labels.length > 1 ? `<span class="muted" style="font-size:0.74rem">colour = main parameter; ${this.result.labels.slice(1).map(esc).join(", ")} move${this.result.labels.length > 2 ? "" : "s"} with it</span>` : "")
      + (this.view.baseline ? `<span class="legend" style="margin:0"><span><svg width="12" height="12" style="vertical-align:-2px;margin-right:0.3em"><circle cx="6" cy="6" r="4.5" fill="none" stroke="#0b0b0b"/></svg>Eigenvalues at the current value</span></span>` : "")
      + `<span class="muted" style="font-size:0.74rem">Scroll to zoom · drag to pan · click a point to follow its mode</span>`;
  },

  // --- Participation of the selected mode, along the sweep ---
  // Clicking a mode opens a bar plot of the states that make it up -- the
  // same reading as the single-mode participation page, but moving: it plays
  // through the sweep, so the composition of the mode can be watched changing
  // with the parameter (a mode handing over from one machine to another, an
  // inner loop taking over as a gain rises). The bars keep a fixed order (by
  // each state's largest participation over the whole sweep) so only their
  // lengths move.
  PART_TOP: 15,
  PART_MS: 320,

  stopParticipation() {
    if (this._partTimer) clearInterval(this._partTimer);
    this._partTimer = null;
    this.part = null;
  },
  pauseParticipation() {
    if (this._partTimer) clearInterval(this._partTimer);
    this._partTimer = null;
    if (this.part) this.part.playing = false;
    const btn = $("#rl-part [data-role=\"play\"]");
    if (btn) btn.textContent = "▶ Play";
  },
  playParticipation() {
    if (!this.part || this.part.pts.length < 2) return;
    this.pauseParticipation();
    this.part.playing = true;
    const btn = $("#rl-part [data-role=\"play\"]");
    if (btn) btn.textContent = "❚❚ Pause";
    this._partTimer = setInterval(() => {
      if (!this.part || !$("#rl-part") || !$("#page-modal")?.classList.contains("active")) return this.pauseParticipation();
      this.setParticipationStep((this.part.step + 1) % this.part.pts.length);
    }, this.PART_MS);
  },

  // Everything the panel shows for one mode, from the steps solved so far.
  participationData(k) {
    const r = this.result;
    if (!r || !r.stateNames) return null;
    const pts = r.steps.filter(s => s.ok && s.part && s.part.length > k)
      .map(s => ({ i: s.i, value: s.value, eig: s.eig[k], rows: s.part[k] }));
    if (!pts.length) return null;
    const best = new Map();
    pts.forEach(p => p.rows.forEach(([i, v]) => best.set(i, Math.max(best.get(i) || 0, v))));
    const states = [...best.entries()].sort((a, b) => b[1] - a[1]).slice(0, this.PART_TOP).map(([i]) => i);
    const at = pts.map(p => { const m = new Map(p.rows); return states.map(i => m.get(i) || 0); });
    return { k, pts, states, at, max: Math.max(0.05, ...best.values()) };
  },

  drawParticipation() {
    const box = $("#rl-part");
    if (!box) return;
    const k = this.view.track;
    if (k === null || k === undefined) {
      this.stopParticipation();
      $("#rl-cursor")?.setAttribute("opacity", "0");
      box.hidden = true;
      $("#rl-main")?.classList.remove("with-side");
      return;
    }
    const d = this.participationData(k);
    box.hidden = false;
    $("#rl-main")?.classList.add("with-side");
    if (!d) {
      box.innerHTML = `<div class="rl-side-head"><h4>Mode ${k}</h4></div><p class="empty">No participation for this mode yet.</p>`;
      return;
    }
    const same = this.part && this.part.k === k && this.part.pts.length === d.pts.length && $("#rl-part .rl-bars");
    const playing = this.part ? this.part.playing : true;
    const step = same ? Math.min(this.part.step, d.pts.length - 1) : 0;
    this.part = { ...d, step, playing };
    if (!same) {
      box.innerHTML = `<div class="rl-side-head"><h4>Mode ${k} — participation</h4>
          <button class="ghost small" data-role="close" aria-label="Close">✕</button></div>
        <p class="status-line" data-role="at"></p>
        <div class="rl-play"><button class="secondary small" data-role="play">▶ Play</button>
          <input type="range" min="0" max="${d.pts.length - 1}" value="${step}" data-role="slider" aria-label="Sweep step"></div>
        <div class="rl-bars">${d.states.map(i => `<div class="rl-bar-row"><span class="rl-bar-label" title="${esc(this.result.stateNames[i])}">${esc(this.result.stateNames[i])}</span>
          <span class="rl-bar-track"><span class="rl-bar-fill"></span></span><span class="rl-bar-val"></span></div>`).join("")}</div>
        <p class="muted" style="font-size:0.72rem;margin:0.6rem 0 0">Top ${d.states.length} states of this mode, in the order of their largest participation over the sweep. Hover a point of its locus to jump to that value.</p>`;
      // Bars are HTML, not a figure the generic export can serialize, so the
      // numbers behind them get their own CSV (exports.js).
      exportToolbar(box, [
        ["CSV", "Download this mode's participation across the sweep",
          () => exportCsv(this.participationRows(), exportName(box))],
        ["Video", "Record this mode's participation playing through the sweep",
          () => this.recordParticipation()],
      ]);
      box.querySelector('[data-role="close"]').addEventListener("click", () => {
        this.view.track = null; this.draw(); this.drawTable(); this.drawParticipation();
      });
      box.querySelector('[data-role="play"]').addEventListener("click", () => {
        if (this.part.playing) this.pauseParticipation(); else this.playParticipation();
      });
      box.querySelector('[data-role="slider"]').addEventListener("input", e => {
        this.pauseParticipation();
        this.setParticipationStep(+e.target.value);
      });
      if (playing) this.playParticipation(); else this.pauseParticipation();
    } else {
      const slider = box.querySelector('[data-role="slider"]');
      slider.max = String(d.pts.length - 1);
    }
    this.setParticipationStep(step);
  },

  // [[parameter value, participation of each of the top states], ...]
  participationRows() {
    const d = this.part, names = this.result.stateNames;
    const rows = [[this.result.labels[0] || "value", ...d.states.map(i => names[i])]];
    d.pts.forEach((pt, k) => rows.push([pt.value, ...d.at[k]]));
    return rows;
  },

  setParticipationStep(step) {
    const box = $("#rl-part");
    if (!box || !this.part) return;
    const d = this.part;
    d.step = Math.max(0, Math.min(step, d.pts.length - 1));
    const p = d.pts[d.step], vals = d.at[d.step];
    const [re, im] = p.eig;
    const wn = Math.hypot(re, im), z = wn ? -re / wn : 0;
    box.querySelector('[data-role="at"]').innerHTML =
      `${this.stepValues(p.i).map(([lab, v]) => `${esc(lab.split(" · ").pop())} = <b>${fmtSmart(v)}</b>`).join(" · ")}
       <br>${fmt(Math.abs(im) / (2 * Math.PI), 3)} Hz · damping ${fmt(100 * z, 2)} % · step ${p.i + 1}/${this.result.values.length}`;
    const slider = box.querySelector('[data-role="slider"]');
    if (slider && +slider.value !== d.step) slider.value = String(d.step);
    this.moveLocusCursor();
    box.querySelectorAll(".rl-bar-row").forEach((row, j) => {
      const v = vals[j] || 0;
      row.querySelector(".rl-bar-fill").style.width = `${(100 * v / d.max).toFixed(1)}%`;
      row.querySelector(".rl-bar-val").textContent = v ? v.toFixed(3) : "—";
    });
  },

  // The ring on the locus at the value the bars are showing.
  moveLocusCursor() {
    const c = $("#rl-cursor");
    if (!c) return;
    const d = this.part;
    if (!d || this.view.track !== d.k || !this._proj) { c.setAttribute("opacity", "0"); return; }
    const [re, im] = d.pts[d.step].eig;
    c.setAttribute("cx", this._proj.sx(re).toFixed(1));
    c.setAttribute("cy", this._proj.sy(im).toFixed(1));
    c.setAttribute("opacity", "1");
  },

  // --- The sweep as a video ---
  // One frame per solved value, the loci growing as the parameter moves. The
  // frames are the page's own plot: each one is drawn with view.upTo set,
  // serialized exactly as the PNG export does (exports.js), and painted onto
  // a canvas that MediaRecorder is recording. MP4 when the browser can write
  // it, WebM otherwise -- both play everywhere that matters, and neither
  // needs anything installed.
  VIDEO_FPS: 25,          // what the recorder samples the canvas at
  VIDEO_BUDGET_MS: 9000,  // roughly how long the sweep should take to play

  videoMimeType() {
    const wanted = [
      ["video/mp4;codecs=avc1.42E01E", "mp4"], ["video/mp4", "mp4"],
      ["video/webm;codecs=vp9", "webm"], ["video/webm", "webm"],
    ];
    for (const [mime, ext] of wanted) {
      if (window.MediaRecorder && MediaRecorder.isTypeSupported(mime)) return { mime, ext };
    }
    return null;
  },

  // Records a playback: `paint(ctx, k)` draws value k onto the canvas and
  // returns a repaint function for it. A fixed capture rate, and the canvas
  // repainted while each value is held -- pushing frames by hand
  // (captureStream(0) + requestFrame) dropped half of them, and a canvas left
  // untouched emits no frames at all, which came out as a video a fraction of
  // its intended length.
  async recordVideo({ btn, name, width, height, count, paint, done }) {
    const fmt_ = this.videoMimeType();
    if (!fmt_) { alert("This browser can't record video (MediaRecorder is unavailable)."); return; }
    const canvas = document.createElement("canvas");
    canvas.width = Math.round(width);
    canvas.height = Math.round(height);
    const ctx = canvas.getContext("2d");
    const rec = new MediaRecorder(canvas.captureStream(this.VIDEO_FPS), {
      mimeType: fmt_.mime, videoBitsPerSecond: 6e6,
    });
    const chunks = [];
    rec.ondataavailable = e => { if (e.data.size) chunks.push(e.data); };
    const holdMs = Math.max(120, Math.min(400, this.VIDEO_BUDGET_MS / count));
    let frame = () => {};
    const repaint = setInterval(() => frame(), Math.round(1000 / this.VIDEO_FPS));
    btn.disabled = true;
    const was = btn.textContent;
    try {
      rec.start();
      for (let k = 0; k < count; k++) {
        btn.textContent = `${Math.round((100 * (k + 1)) / count)} %`;
        frame = (await paint(ctx, k)) || frame;
        frame();
        await new Promise(r => setTimeout(r, k === count - 1 ? 1600 : holdMs));
      }
      await new Promise(res => { rec.onstop = res; rec.stop(); });
      downloadBlob(`${exportSlug(name)}.${fmt_.ext}`, new Blob(chunks, { type: fmt_.mime }));
    } catch (e) {
      alert(`Could not record the video: ${e.message}`);
    } finally {
      clearInterval(repaint);
      btn.disabled = false;
      btn.textContent = was;
      done?.();
    }
  },

  // The sweep itself: the loci growing, value by value, with a caption strip.
  async recordSweep() {
    const steps = (this.result?.steps || []).filter(s => s.ok);
    if (steps.length < 2) { alert("Run a sweep first — a video needs at least two solved values."); return; }
    const svg = $("#rl-plot svg");
    if (!svg) return;
    const vb = svg.viewBox.baseVal, scale = 1.5;
    const keepUpTo = this.view.upTo, keepTrack = this.view.track;
    const net = state.networkLabel ? `${state.networkLabel}-` : "";
    await this.recordVideo({
      btn: $("#rl-video"),
      name: `${net}root-locus-${this.result.labels[0] || "sweep"}`,
      width: vb.width * scale,
      height: (vb.height + 34) * scale,   // room for the caption strip
      count: steps.length,
      done: () => {
        this.view.upTo = keepUpTo;
        this.view.track = keepTrack;
        this.draw();
        this.drawParticipation();
      },
      paint: async (ctx, k) => {
        const step = steps[k];
        this.view.upTo = step.i;
        this.draw();
        const url = URL.createObjectURL(new Blob([serializeSvg($("#rl-plot svg"))], { type: "image/svg+xml;charset=utf-8" }));
        try {
          const img = new Image();
          await new Promise((res, rej) => { img.onload = res; img.onerror = rej; img.src = url; });
          const W = Math.round(vb.width * scale), H = Math.round((vb.height + 34) * scale);
          return () => {
            ctx.fillStyle = "#ffffff";
            ctx.fillRect(0, 0, W, H);
            ctx.drawImage(img, 0, 0, W, Math.round(vb.height * scale));
            ctx.fillStyle = "#0b0b0b";
            ctx.font = `${Math.round(13 * scale)}px "IBM Plex Mono", monospace`;
            const vals = this.stepValues(step.i).map(([lab, v]) => `${lab.split(" · ").pop()} = ${fmtSmart(v)}`).join("   ");
            ctx.fillText(`${vals}   ·   value ${step.i + 1}/${this.result.values.length}`, 10 * scale, H - 12 * scale);
            const tracked = this.view.track;
            const eig = tracked === null || tracked === undefined ? null : step.eig[tracked];
            if (eig) {
              const [re, im] = eig, wn = Math.hypot(re, im), z = wn ? -100 * re / wn : 0;
              ctx.fillStyle = z < 0 ? "#d03b3b" : "#52514e";
              ctx.textAlign = "right";
              ctx.fillText(`mode ${tracked}: ${fmt(Math.abs(im) / (2 * Math.PI), 2)} Hz · damping ${fmt(z, 2)} %`,
                           W - 10 * scale, H - 12 * scale);
              ctx.textAlign = "left";
            }
          };
        } finally { URL.revokeObjectURL(url); }
      },
    });
  },

  // The participation panel playing through the sweep. Its bars are HTML, so
  // the frames are drawn straight onto the canvas rather than serialized.
  async recordParticipation() {
    const d = this.part;
    if (!d || d.pts.length < 2) { alert("Select a mode first — a video needs at least two solved values."); return; }
    const names = this.result.stateNames;
    const S = 2;                                   // drawn at twice the size
    const rowH = 26, top = 92, labelW = 250, valueW = 72, pad = 16;
    const W = 900, H = top + d.states.length * rowH + pad;
    const net = state.networkLabel ? `${state.networkLabel}-` : "";
    const keepStep = d.step;
    this.pauseParticipation();
    await this.recordVideo({
      btn: $('#rl-part .fig-tools button:last-child'),
      name: `${net}mode-${d.k}-participation`,
      width: W * S, height: H * S,
      count: d.pts.length,
      done: () => this.setParticipationStep(keepStep),
      paint: (ctx, k) => {
        const pt = d.pts[k], vals = d.at[k];
        const [re, im] = pt.eig, wn = Math.hypot(re, im), z = wn ? -100 * re / wn : 0;
        return () => {
          ctx.save();
          ctx.scale(S, S);
          ctx.fillStyle = "#ffffff";
          ctx.fillRect(0, 0, W, H);
          ctx.fillStyle = "#0b0b0b";
          ctx.font = '600 19px "IBM Plex Sans", system-ui, sans-serif';
          ctx.fillText(`Mode ${d.k} — participation`, pad, 30);
          ctx.font = '13px "IBM Plex Mono", monospace';
          ctx.fillStyle = "#52514e";
          const vs = this.stepValues(pt.i).map(([lab, v]) => `${lab.split(" · ").pop()} = ${fmtSmart(v)}`).join("   ");
          ctx.fillText(`${vs}   ·   value ${pt.i + 1}/${this.result.values.length}`, pad, 54);
          ctx.fillStyle = z < 0 ? "#d03b3b" : "#52514e";
          ctx.fillText(`${fmt(Math.abs(im) / (2 * Math.PI), 3)} Hz · damping ${fmt(z, 2)} %`, pad, 76);
          d.states.forEach((si, i) => {
            const y = top + i * rowH, v = vals[i] || 0;
            ctx.font = '12px "IBM Plex Mono", monospace';
            ctx.fillStyle = "#52514e";
            ctx.textAlign = "right";
            ctx.fillText(names[si], labelW, y + 13);
            ctx.textAlign = "left";
            const x0 = labelW + 14, trackW = W - x0 - valueW - pad;
            ctx.fillStyle = "#f7f6f3";
            ctx.fillRect(x0, y + 2, trackW, 15);
            ctx.fillStyle = "#2a78d6";
            ctx.fillRect(x0, y + 2, Math.max(0, (v / d.max) * trackW), 15);
            ctx.fillStyle = "#898781";
            ctx.textAlign = "right";
            ctx.fillText(v ? v.toFixed(3) : "—", W - pad, y + 13);
            ctx.textAlign = "left";
          });
          ctx.restore();
        };
      },
    });
  },

  // Modes ranked by how far they move -- the ones this parameter really affects.
  drawTable() {
    const box = $("#rl-table");
    if (!box || !this.result) return;
    const rows = this.tracks().filter(t => t.pts.length > 1).map(t => {
      const a = t.pts[0], b = t.pts[t.pts.length - 1];
      const damp = p => { const w = Math.hypot(p.re, p.im); return w ? -100 * p.re / w : 100; };
      const minD = Math.min(...t.pts.map(damp));
      const moved = Math.max(...t.pts.map(p => Math.hypot(p.re - a.re, p.im - a.im) / (1 + Math.hypot(a.re, a.im))));
      const unstable = t.pts.find(p => p.re > 1e-6);
      return { t, a, b, minD, moved, unstable, damp };
    }).sort((x, y) => y.moved - x.moved).slice(0, 25);
    if (!rows.length) { box.innerHTML = `<p class="empty">Need at least two solved values.</p>`; return; }
    const lam = p => `${p.re.toExponential(2)} ${p.im >= 0 ? "+" : "−"} j${Math.abs(p.im).toExponential(2)}`;
    box.innerHTML = `<div class="tablewrap" style="max-height:420px"><table><thead><tr><th>Mode</th><th>λ at start</th><th>λ at end</th><th>Freq start → end (Hz)</th><th>Damping start → end (%)</th><th>Min damping (%)</th><th>Relative shift</th><th>Stability</th></tr></thead><tbody>
      ${rows.map(r => `<tr class="clickable${r.t.k === this.view.track ? " hl" : ""}" data-k="${r.t.k}"><td>${r.t.k}</td><td>${lam(r.a)}</td><td>${lam(r.b)}</td>
        <td>${fmt(Math.abs(r.a.im) / (2 * Math.PI), 3)} → ${fmt(Math.abs(r.b.im) / (2 * Math.PI), 3)}</td><td>${fmt(r.damp(r.a), 2)} → ${fmt(r.damp(r.b), 2)}</td>
        <td>${fmt(r.minD, 2)}</td><td>${fmtSmart(r.moved)}</td>
        <td class="name">${r.unstable ? `<span class="bad">unstable from step ${r.unstable.i + 1} (${esc(this.result.labels[0].split(" · ").pop())} = ${fmtSmart(r.unstable.value)})</span>` : '<span class="ok">stable</span>'}</td></tr>`).join("")}</tbody></table></div>`;
    $$("#rl-table tr[data-k]").forEach(tr => tr.addEventListener("click", () => {
      this.view.track = +tr.dataset.k === this.view.track ? null : +tr.dataset.k;
      this.draw(); this.drawTable(); this.drawParticipation();
    }));
  },
};
