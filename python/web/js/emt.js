// EMT simulation page: integrates the nonlinear network model after a state
// offset or an input step at T0 = 0. Trajectories start T_PRE before T0 so
// the undisturbed x0 is visible; optionally the linearised model's response
// to the same disturbance is overlaid as dotted lines in the same colours.

const EMT_T_PRE = 0.01;

const EmtPage = {
  inited: false,
  pickers: {},
  names: null,
  running: null,

  init() {
    if (this.inited) return;
    this.inited = true;
    $("#page-emt").innerHTML = `
      <div class="page-head"><div class="crumb">Analysis</div><h1>EMT Simulation</h1>
        <p>Integrates the network's full nonlinear differential-algebraic equations (a coupled Newton solve at every step) after perturbing one state or stepping one input at T0 = 0 — the ground truth that modal analysis linearises. Plots start ${EMT_T_PRE} s before T0 to show the initial operating point x0.</p></div>
      <div id="emt-context"></div>
      <div class="card">
        <div class="card-title">Disturbance</div>
        <div class="controls">
          <div class="field"><label for="emt-perturb">Perturb (state offset or input step)</label><select id="emt-perturb" style="min-width:260px"></select></div>
          <div class="field"><label for="emt-amp">Amplitude</label><input type="number" id="emt-amp" step="any" value="0.02"></div>
          <div class="field"><label for="emt-tf">Duration after T0 (s, max 3)</label><input type="number" id="emt-tf" step="0.1" min="0.05" max="3" value="1.0"></div>
          <div class="field"><label for="emt-dt">Output timestep (s)</label><input type="number" id="emt-dt" step="any" placeholder="auto"></div>
        </div>
        <div class="card-title" style="margin-top:1.2rem">Signals to plot</div>
        <div class="controls" style="align-items:flex-start">
          <div class="field"><label>States</label><div id="emt-p-states"></div></div>
          <div class="field"><label>Inputs</label><div id="emt-p-inputs"></div></div>
          <div class="field"><label>Outputs</label><div id="emt-p-outputs"></div></div>
        </div>
        <div class="controls" style="margin-top:1.1rem;align-items:center">
          <label class="check"><input type="checkbox" id="emt-linear"> Overlay the equivalent linearised response (dotted, same colours)</label>
          <label class="check"><input type="checkbox" id="emt-live"> Trace live <span class="muted" style="font-size:0.76rem">(streams solver steps; slower)</span></label>
        </div>
        <div class="controls" style="margin-top:1rem;align-items:center">
          <button id="emt-run"><svg viewBox="0 0 24 24" fill="currentColor"><path d="M8 5v14l11-7z"/></svg>Run EMT simulation</button>
          <button class="secondary" id="emt-stop" style="display:none">Stop</button>
          <span class="spinner" id="emt-spin"></span>
        </div>
      </div>
      <div id="emt-results" class="stack" style="margin-top:1rem"></div>`;
    ["states", "inputs", "outputs"].forEach(k => {
      this.pickers[k] = new SignalPicker($(`#emt-p-${k}`), { options: [], colors: true, placeholder: `add ${k.slice(0, -1)}` });
    });
    $("#emt-run").addEventListener("click", () => this.run());
    on("network:loaded", () => { this.names = null; $("#emt-results").innerHTML = ""; });
    on("network:changed", () => { this.names = null; });
  },

  async onShow() {
    $("#emt-context").innerHTML = networkContextHtml();
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
      const prevSel = sel.value;
      sel.innerHTML = `<optgroup label="States (initial-condition offset)">${r.state_names.map(n => `<option value="${esc(n)}" data-kind="state">${esc(n)}</option>`).join("")}</optgroup>
        <optgroup label="Inputs (step held from T0)">${r.input_names.map(n => `<option value="${esc(n)}" data-kind="input">${esc(n)}</option>`).join("")}</optgroup>`;
      const firstDw = r.state_names.find(n => n.includes("dw_r"));
      sel.value = [...sel.options].some(o => o.value === prevSel) ? prevSel : (firstDw || r.state_names[0]);
      sel.disabled = false;
      const keepOrDefault = (picker, names, def) => {
        const kept = picker.get().filter(n => names.includes(n));
        picker.setOptions(names);
        picker.set(kept.length ? kept : def);
      };
      keepOrDefault(this.pickers.states, r.state_names, r.state_names.filter(n => n.includes("dw_r")).slice(0, 8));
      keepOrDefault(this.pickers.inputs, r.input_names, []);
      keepOrDefault(this.pickers.outputs, r.output_names, []);
    } catch (e) {
      sel.innerHTML = `<option>unavailable</option>`;
      $("#emt-results").innerHTML = errorHtml(e);
    } finally { setSpinner($("#emt-spin"), ""); }
  },

  request() {
    const sel = $("#emt-perturb");
    const opt = sel.selectedOptions[0];
    const dtRaw = $("#emt-dt").value.trim();
    return {
      perturb_kind: opt ? opt.dataset.kind : "state",
      perturb_name: sel.value,
      perturb_offset: parseFloat($("#emt-amp").value) || 0,
      t_final: parseFloat($("#emt-tf").value) || 1.0,
      dt: dtRaw ? parseFloat(dtRaw) : null,
      plot_states: this.pickers.states.get(),
      plot_inputs: this.pickers.inputs.get(),
      plot_outputs: this.pickers.outputs.get(),
      t_pre: EMT_T_PRE,
      linear_overlay: $("#emt-linear").checked,
    };
  },

  async run() {
    if (!this.names) { await this.onShow(); if (!this.names) return; }
    const req = this.request();
    const out = $("#emt-results"), btn = $("#emt-run"), spin = $("#emt-spin");
    btn.disabled = true;
    out.innerHTML = "";
    try {
      if ($("#emt-live").checked) await this.runLive(req);
      else {
        setSpinner(spin, `Integrating${req.plot_inputs.length || req.plot_outputs.length ? " (recovering inputs/outputs adds a Newton solve per sample)" : ""}${req.linear_overlay ? " + linearised response" : ""}…`);
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
    out.innerHTML = `<div class="card"><p class="status-line">Perturbed with ${kind} <b>${esc(r.perturbed)}</b> of ${fmt(req.perturb_offset, 4)} · ${fmt(req.t_final, 2)} s after T0${r.dt ? ` · dt = ${fmt(r.dt, 5)} s` : ""}${extraNote}</p>
      ${r.linear ? `<div class="legend"><span><span class="line-swatch" style="border-color:var(--text-secondary)"></span>Nonlinear EMT</span><span><span class="line-swatch dotted" style="border-color:var(--text-secondary)"></span>Linearised model (same disturbance, same operating point)</span></div>
        <p class="muted" style="font-size:0.76rem;margin:0.5rem 0 0">Note: the nonlinear model's initial point is not an exact equilibrium (every node borrows line #1's susceptance), so the EMT traces include a small initial transient even without a disturbance. The linearised response has no such transient — part of the gap between the curves comes from that, not from nonlinearity.</p>` : ""}</div>`;
    [["States", r.series, r.linear && r.linear.series], ["Inputs", r.inputs, r.linear && r.linear.inputs], ["Outputs", r.outputs, r.linear && r.linear.outputs]].forEach(([title, series, lin]) => {
      const names = Object.keys(series || {});
      if (!names.length) return;
      const s = [];
      names.forEach((n, i) => {
        s.push({ name: n, t: r.t, y: series[n], color: seriesColor(i) });
        if (lin && lin[n]) s.push({ name: n, t: r.linear.t, y: lin[n], color: seriesColor(i), dash: true });
      });
      const card = el(`<div class="card"><div class="card-title">${title} <span class="card-sub">(${names.length})</span></div></div>`);
      card.appendChild(lineChart(s, { height: 300, markers: [{ x: 0, label: "T0" }] }));
      card.insertAdjacentHTML("beforeend", legendHtml(names.map((n, i) => ({ name: n, color: seriesColor(i) }))));
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
    const data = { t: [], series: {}, inputs: {}, outputs: {} };
    let buf = "", done = null, error = null, lastDraw = 0;
    const draw = force => {
      const now = performance.now();
      if (!force && now - lastDraw < 120) return;
      lastDraw = now;
      this.renderResult({ ...data, perturb_kind: req.perturb_kind, perturbed: req.perturb_name, dt: null }, req, ` · <span class="muted">tracing live — ${data.t.length} points</span>`);
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
          for (const [grp, key] of [["series", "states"], ["inputs", "inputs"], ["outputs", "outputs"]]) {
            for (const [n, v] of Object.entries(msg[key])) (data[grp][n] ||= []).push(v);
          }
          draw(false);
        }
      }
    } catch { if (!error) done = done || { stopped: true }; }
    finally { stop.style.display = "none"; }
    if (error) throw new Error(error);
    const note = done && !done.stopped ? ` · traced live, ${done.n_steps} solver steps` : ` · <span class="warn">stopped</span> — partial trace`;
    this.renderResult({ ...data, perturb_kind: req.perturb_kind, perturbed: req.perturb_name, dt: null, linear: done && done.linear }, req, note);
    if (done && done.linear_error) out.insertAdjacentHTML("afterbegin", `<div class="notice warn-bg">Linearised overlay unavailable: ${esc(done.linear_error)}</div>`);
  },
};
