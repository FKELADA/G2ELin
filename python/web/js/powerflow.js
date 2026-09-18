// Power Flow page: solver settings, single power flow on the initial operating
// point, and a batch (time-series-like) ramp of power flows with load and
// unit-setpoint scaling. Results are drawn on the network diagram (heatmap
// colours, per-element hover) and listed in a table below.

const PF_TABLES = {
  buses: "Bus results", lines: "Line results", transformers: "Transformer results", loads: "Load results",
  generators: "Generator (PV) results", static_generators: "Static generator (PQ) results", external_grid: "Slack (external grid) results",
};
const BUS_METRICS = {
  vm_pu: { label: "Voltage magnitude", unit: " pu", digits: 3, get: b => b.vm_pu, center: 1.0 },
  va_degree: { label: "Voltage angle", unit: "°", digits: 2, get: b => b.va_degree, center: 0.0 },
};
const EDGE_METRICS = {
  p: { label: "Active power flow |P|", unit: " MW", line: r => Math.abs(r.p_from_mw), trafo: r => Math.abs(r.p_hv_mw) },
  q: { label: "Reactive power flow |Q|", unit: " MVAr", line: r => Math.abs(r.q_from_mvar), trafo: r => Math.abs(r.q_hv_mvar) },
  loss: { label: "Active losses", unit: " MW", line: r => r.pl_mw, trafo: r => r.pl_mw },
  i: { label: "Current", unit: " kA", line: r => r.i_ka, trafo: r => r.i_hv_ka },
};

const PowerFlowPage = {
  inited: false,
  view: null,
  result: null,        // {version, kind: "single"|"batch", data, idx}
  batchOpen: false,
  selectedRow: null,

  init() {
    if (this.inited) return;
    this.inited = true;
    const page = $("#page-powerflow");
    page.innerHTML = `
      <div class="page-head"><div class="crumb">Analysis</div><h1>Power Flow</h1>
        <p>Steady-state AC solution: every bus voltage and angle, and every branch flow, at the current dispatch — the operating point modal analysis and EMT simulation start from. Run it once on the initial operating point, or as a batch while scaling loads and unit setpoints.</p></div>
      <div id="pf-context"></div>
      <div class="card">
        <div class="card-title">Solver settings</div>
        <div class="controls">
          <div class="field"><label for="pf-algo">Solver</label><select id="pf-algo" style="min-width:280px"></select></div>
          <div class="field"><label for="pf-maxit">Max iterations</label><input type="number" id="pf-maxit" min="1" step="1" placeholder="auto"></div>
          <div class="field"><label for="pf-tol">Tolerance (MVA)</label><input type="number" id="pf-tol" step="any" value="1e-8"></div>
          <div class="field"><label for="pf-init">Initialisation</label><select id="pf-init"><option value="auto">auto</option><option value="flat">flat start</option><option value="dc">DC power flow</option></select></div>
        </div>
        <div class="pf-actions" style="margin-top:1rem">
          <button id="pf-run"><svg viewBox="0 0 24 24" fill="currentColor"><path d="M8 5v14l11-7z"/></svg>Run power flow</button>
          <button class="secondary toggle" id="pf-batch-toggle" aria-pressed="false"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M4 18h4V9H4zM10 18h4V5h-4zM16 18h4v-6h-4z"/></svg>Batch power flow</button>
          <span class="spinner" id="pf-spinner"></span>
        </div>
        <div class="batch-panel" id="pf-batch">
          <div class="card-title" style="margin-bottom:0.3rem">Batch power flow <span class="card-sub">— a ramp of independent power flows from the initial operating point (100 %) to the targets below</span></div>
          <p class="muted" style="font-size:0.78rem;margin:0 0 0.9rem">Unit sliders scale each unit's active <i>and</i> reactive power setpoints. The slack unit's P, and PV units' Q, are outcomes of the power flow, so scaling those setpoints has no effect.</p>
          <div class="sliders" id="pf-sliders"></div>
          <div class="controls" style="margin-top:1rem;align-items:center">
            <div class="field"><label for="pf-steps">Steps</label><input type="number" id="pf-steps" min="1" max="100" value="10"></div>
            <button id="pf-run-batch">Run batch power flow</button>
            <button class="ghost" id="pf-reset-sliders">Reset to 100 %</button>
          </div>
        </div>
      </div>
      <div class="card">
        <div class="card-title">Network <span class="card-sub" id="pf-plot-sub">— hover an element for its results</span></div>
        <div class="plot-opts">
          <div class="field"><label for="pf-bus-metric">Colour buses by</label><select id="pf-bus-metric">${Object.entries(BUS_METRICS).map(([k, m]) => `<option value="${k}">${m.label}</option>`).join("")}</select></div>
          <div class="field"><label for="pf-edge-metric">Colour lines by</label><select id="pf-edge-metric">${Object.entries(EDGE_METRICS).map(([k, m]) => `<option value="${k}">${m.label}</option>`).join("")}</select></div>
          <span style="flex:1"></span>
          <p class="status-line" id="pf-status"></p>
        </div>
        <div id="pf-scrubber"></div>
        <div id="pf-canvas"></div>
        <div class="plot-legends" id="pf-legends"></div>
      </div>
      <div id="pf-batch-charts"></div>
      <div class="card" id="pf-results-card">
        <div class="card-title">Results</div>
        <div id="pf-table-area"><p class="empty">Run a power flow to see results.</p></div>
      </div>`;

    this.view = new NetworkView($("#pf-canvas"), {
      editable: false,
      tooltip: sel => this.tooltip(sel),
      onSelect: sel => this.onSelect(sel),
    });
    this.fillAlgorithms();
    $("#pf-run").addEventListener("click", () => this.runSingle());
    $("#pf-batch-toggle").addEventListener("click", () => this.toggleBatch());
    $("#pf-run-batch").addEventListener("click", () => this.runBatch());
    $("#pf-reset-sliders").addEventListener("click", () => { $$("#pf-sliders input[type=range]").forEach(r => { r.value = 100; r.dispatchEvent(new Event("input")); }); });
    $("#pf-bus-metric").addEventListener("change", () => this.redraw());
    $("#pf-edge-metric").addEventListener("change", () => this.redraw());
    on("network:loaded", () => this.onNetwork(true));
    on("network:changed", () => this.onNetwork(false));
    on("network:layout", () => { this.view.render(); this.view.fit(); });
    this.onNetwork(true);
  },

  onShow() {
    $("#pf-context").innerHTML = networkContextHtml();
    this.view.render();
    if (!this._fitted && state.network) { this.view.fit(); this._fitted = true; }
  },

  async fillAlgorithms() {
    let algos = { nr: "Newton-Raphson" };
    try { algos = await api("/api/powerflow/algorithms"); } catch { /* keep the default */ }
    $("#pf-algo").innerHTML = Object.entries(algos).map(([k, v]) => `<option value="${k}">${esc(v)}</option>`).join("");
  },

  onNetwork(loaded) {
    this.result = null;
    this.buildSliders();
    if (loaded) this._fitted = false;
    if ($("#page-powerflow").classList.contains("active")) this.onShow();
    this.redraw();
  },

  options() {
    const maxit = $("#pf-maxit").value.trim();
    const tol = parseFloat($("#pf-tol").value);
    return {
      algorithm: $("#pf-algo").value || "nr",
      max_iteration: maxit ? parseInt(maxit, 10) : "auto",
      tolerance_mva: Number.isFinite(tol) && tol > 0 ? tol : 1e-8,
      init: $("#pf-init").value,
    };
  },

  buildSliders() {
    const box = $("#pf-sliders");
    if (!box) return;
    const types = state.network ? [...new Set(state.network.der_units.map(d => d.unit_type))].filter(t => t !== "infinite_bus") : [];
    const defs = [
      { key: "load_p", label: "Loads — active power P" },
      { key: "load_q", label: "Loads — reactive power Q" },
      ...types.map(t => ({ key: `der:${t}`, label: `${UNIT_LABEL[t]} setpoints — P & Q`, color: UNIT_COLOR[t] })),
    ];
    const prev = Object.fromEntries($$("#pf-sliders input[type=range]").map(r => [r.dataset.key, r.value]));
    box.innerHTML = defs.map(d => `<div class="slider-field">
        <div class="slider-head"><label for="pf-s-${d.key}">${d.color ? `<span class="swatch" style="display:inline-block;width:9px;height:9px;border-radius:50%;background:${d.color};margin-right:0.35rem"></span>` : ""}${esc(d.label)}</label><span class="val" id="pf-v-${d.key}"></span></div>
        <input type="range" id="pf-s-${d.key}" data-key="${d.key}" min="0" max="200" step="1" value="${prev[d.key] ?? 100}">
        <div class="ticks"><span>0%</span><span>50%</span><span>100%</span><span>150%</span><span>200%</span></div></div>`).join("");
    $$("#pf-sliders input[type=range]").forEach(r => {
      const show = () => { $(`#pf-v-${CSS.escape(r.dataset.key)}`).textContent = `${r.value} %`; };
      r.addEventListener("input", show);
      show();
    });
  },

  toggleBatch() {
    this.batchOpen = !this.batchOpen;
    $("#pf-batch").classList.toggle("open", this.batchOpen);
    $("#pf-batch-toggle").setAttribute("aria-pressed", String(this.batchOpen));
  },

  async runSingle() {
    const btn = $("#pf-run"), spin = $("#pf-spinner");
    btn.disabled = true; setSpinner(spin, "Solving…");
    const version = state.version;
    try {
      const r = await netPost("powerflow", { options: this.options() });
      if (version !== state.version) return;
      this.result = { version, kind: "single", data: r, idx: 0 };
      this.afterResult();
    } catch (e) {
      this.result = null; this.afterResult();
      $("#pf-table-area").innerHTML = errorHtml(e);
    } finally { btn.disabled = false; setSpinner(spin, ""); }
  },

  async runBatch() {
    const btn = $("#pf-run-batch"), spin = $("#pf-spinner");
    const sliders = Object.fromEntries($$("#pf-sliders input[type=range]").map(r => [r.dataset.key, parseInt(r.value, 10) / 100]));
    const der_scale = {};
    Object.entries(sliders).forEach(([k, v]) => { if (k.startsWith("der:")) der_scale[k.slice(4)] = v; });
    const steps = Math.max(1, Math.min(100, parseInt($("#pf-steps").value, 10) || 10));
    btn.disabled = true; setSpinner(spin, `Solving ${steps + 1} power flows…`);
    const version = state.version;
    try {
      const r = await netPost("powerflow/batch", {
        options: this.options(), load_p_scale: sliders.load_p ?? 1, load_q_scale: sliders.load_q ?? 1, der_scale, steps,
      });
      if (version !== state.version) return;
      this.result = { version, kind: "batch", data: r, idx: r.snapshots.length - 1 };
      this.afterResult();
    } catch (e) {
      $("#pf-table-area").innerHTML = errorHtml(e);
    } finally { btn.disabled = false; setSpinner(spin, ""); }
  },

  current() {
    const r = this.result;
    if (!r || r.version !== state.version) return null;
    return r.kind === "single" ? r.data : r.data.snapshots[r.idx];
  },
  allResults() {
    const r = this.result;
    if (!r || r.version !== state.version) return [];
    return (r.kind === "single" ? [r.data] : r.data.snapshots).filter(s => s.converged);
  },

  afterResult() {
    this.renderScrubber();
    this.renderBatchCharts();
    this.redraw();
    this.renderTable();
  },

  // Colour domains are fixed across a whole batch so colours stay comparable
  // while scrubbing through the snapshots.
  domains() {
    const all = this.allResults();
    const bm = BUS_METRICS[$("#pf-bus-metric").value], em = EDGE_METRICS[$("#pf-edge-metric").value];
    let dev = 0, emax = 0;
    all.forEach(s => {
      s.buses.forEach(b => { dev = Math.max(dev, Math.abs(bm.get(b) - bm.center)); });
      s.lines.forEach(l => { emax = Math.max(emax, em.line(l) || 0); });
      s.transformers.forEach(t => { emax = Math.max(emax, em.trafo(t) || 0); });
    });
    dev = Math.max(dev, bm.center === 1 ? 0.005 : 0.5);
    return { bm, em, lo: bm.center - dev, hi: bm.center + dev, emax: emax || 1 };
  },

  redraw() {
    if (!this.view) return;
    const cur = this.current();
    const legends = $("#pf-legends");
    if (!cur || !cur.converged) {
      this.view.setStyle(null);
      legends.innerHTML = unitLegendHtml();
      $("#pf-status").innerHTML = cur && !cur.converged ? `<span class="bad">did not converge</span>` : (state.network ? "not solved yet" : "");
      return;
    }
    const d = this.domains();
    const byBus = Object.fromEntries(cur.buses.map(b => [b.bus, b]));
    this.view.setStyle({
      bus: b => byBus[b.id] ? { fill: divergingColor((d.bm.get(byBus[b.id]) - d.lo) / (d.hi - d.lo)) } : {},
      edge: (kind, i) => {
        const row = kind === "line" ? cur.lines[i] : cur.transformers[i];
        if (!row) return {};
        const v = kind === "line" ? d.em.line(row) : d.em.trafo(row);
        const t = (v || 0) / d.emax;
        return { stroke: heatColor(t), width: 2 + 5 * t };
      },
    });
    legends.innerHTML = gradientLegendHtml(`Buses: ${d.bm.label}`, DIVERGING, d.lo, d.hi, d.bm.digits, d.bm.unit)
      + gradientLegendHtml(`Lines: ${d.em.label}`, SEQ_HEAT, 0, d.emax, 3, d.em.unit)
      + `<span class="muted" style="font-size:0.74rem">Unit buses keep their type colour as a ring.</span>`;
    const r = this.result;
    const diag = r.kind === "single" && r.data.iterations != null ? ` · ${r.data.iterations} iterations` : "";
    const time = r.kind === "single" && r.data.solve_time_s != null ? ` · ${fmt(r.data.solve_time_s * 1000, 1)} ms` : "";
    $("#pf-status").innerHTML = `<span class="ok">converged</span>${diag}${time} · losses ${fmt(cur.total_losses_mw, 3)} MW`;
  },

  renderScrubber() {
    const box = $("#pf-scrubber");
    const r = this.result;
    if (!r || r.kind !== "batch") { box.innerHTML = ""; return; }
    const n = r.data.snapshots.length;
    box.innerHTML = `<div class="scrubber">
      <button class="secondary small" id="pf-play" title="Play through the snapshots">▶ Play</button>
      <span class="field-label">Snapshot</span>
      <input type="range" id="pf-snap" min="0" max="${n - 1}" step="1" value="${r.idx}">
      <span class="snap-label" id="pf-snap-label"></span></div>`;
    const input = $("#pf-snap");
    const update = () => {
      r.idx = parseInt(input.value, 10);
      const s = r.data.snapshots[r.idx];
      const ders = Object.entries(s.der_scale).map(([t, f]) => `${UNIT_LABEL[t] || t} ${Math.round(f * 100)}%`).join(" · ");
      $("#pf-snap-label").innerHTML = `${r.idx}/${n - 1} · load P ${Math.round(s.load_p_scale * 100)}% · Q ${Math.round(s.load_q_scale * 100)}%${ders ? " · " + ders : ""}${s.converged ? "" : ' · <span class="bad">not converged</span>'}`;
      this.redraw();
      this.renderTable();
    };
    input.addEventListener("input", update);
    $("#pf-play").addEventListener("click", () => {
      clearInterval(this._play);
      let k = 0;
      input.value = 0; update();
      this._play = setInterval(() => {
        k += 1;
        if (k >= n || !document.body.contains(input)) { clearInterval(this._play); return; }
        input.value = k; update();
      }, 450);
    });
    update();
  },

  renderBatchCharts() {
    const box = $("#pf-batch-charts");
    const r = this.result;
    if (!r || r.kind !== "batch") { box.innerHTML = ""; return; }
    const snaps = r.data.snapshots;
    const t = snaps.map((_, i) => i);
    const xFormat = v => {
      const s = snaps[Math.max(0, Math.min(snaps.length - 1, Math.round(v)))];
      return `#${Math.round(v)} (P ${Math.round(s.load_p_scale * 100)}%)`;
    };
    const buses = state.network.buses.map(b => b.id);
    const vSeries = buses.map((id, i) => ({
      name: `bus ${id}`, t, color: seriesColor(i),
      y: snaps.map(s => { const b = s.buses.find(x => x.bus === id); return b ? b.vm_pu : NaN; }),
    }));
    const lossSeries = [{ name: "Total losses (MW)", t, color: seriesColor(0), y: snaps.map(s => s.converged ? s.total_losses_mw : NaN) }];
    const slack = [{ name: "Slack P (MW)", t, color: seriesColor(1), y: snaps.map(s => s.converged && s.external_grid[0] ? s.external_grid[0].p_mw : NaN) }];
    box.innerHTML = `<div class="card"><div class="card-title">Batch overview <span class="card-sub">— across the ${snaps.length} snapshots</span></div><div class="two-col"><div id="pf-ch-v"></div><div id="pf-ch-l"></div></div></div>`;
    const opts = { xUnit: "", xLabel: "snapshot", xFormat, xTickFormat: v => `#${Math.round(v)}`, height: 250 };
    $("#pf-ch-v").appendChild(lineChart(vSeries, { ...opts, title: "Bus voltage magnitudes (pu)" }));
    $("#pf-ch-v").insertAdjacentHTML("beforeend", legendHtml(vSeries));
    $("#pf-ch-l").appendChild(lineChart([...lossSeries, ...slack], { ...opts, title: "Losses and slack injection (MW)" }));
    $("#pf-ch-l").insertAdjacentHTML("beforeend", legendHtml([...lossSeries, ...slack]));
  },

  // "Show" selector sits directly above the table it controls.
  renderTable() {
    const area = $("#pf-table-area");
    const cur = this.current();
    if (!cur) { area.innerHTML = `<p class="empty">Run a power flow to see results.</p>`; return; }
    if (!cur.converged) { area.innerHTML = `<p class="status-line"><span class="bad">This power flow did not converge</span> — try another solver, more iterations or a different initialisation.</p>`; return; }
    const prev = $("#pf-show") ? $("#pf-show").value : "buses";
    area.innerHTML = `<div class="controls" style="margin-bottom:0.6rem">
        <div class="field"><label for="pf-show">Show</label><select id="pf-show">${Object.entries(PF_TABLES).map(([k, v]) => `<option value="${k}"${k === prev ? " selected" : ""}>${v} (${cur[k].length})</option>`).join("")}</select></div></div>
      <div id="pf-table"></div>`;
    $("#pf-show").addEventListener("change", () => { this.selectedRow = null; this.renderTableBody(); });
    this.renderTableBody();
  },

  renderTableBody() {
    const cur = this.current();
    const kind = $("#pf-show").value;
    const rows = cur[kind];
    const box = $("#pf-table");
    if (!rows.length) { box.innerHTML = `<p class="empty">No ${PF_TABLES[kind].toLowerCase()} for this network.</p>`; return; }
    const cols = kind === "buses" ? ["bus", "vm_pu", "va_degree", "p_net_gen_mw", "q_net_gen_mvar"] : Object.keys(rows[0]);
    const label = { bus: "Bus", vm_pu: "V (pu)", va_degree: "Angle (deg)", p_net_gen_mw: "P net gen (MW)", q_net_gen_mvar: "Q net gen (MVAr)" };
    const hl = this.selectedRow && this.selectedRow.kind === kind ? this.selectedRow.index : -1;
    box.innerHTML = `<div class="tablewrap"><table><thead><tr>${cols.map(c => `<th>${esc(label[c] || c)}</th>`).join("")}</tr></thead><tbody>${rows.map((row, i) =>
      `<tr class="${i === hl ? "hl" : ""}" data-i="${i}">${cols.map(c => `<td${typeof row[c] === "string" ? ' class="name"' : ""}>${typeof row[c] === "number" ? fmtSmart(row[c]) : esc(row[c] ?? "")}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`;
    const hlRow = box.querySelector("tr.hl");
    if (hlRow) hlRow.scrollIntoView({ block: "nearest" });
  },

  // Clicking an element on the plot jumps the table to that element's row.
  onSelect(sel) {
    this.view.selected = sel;
    this.view.render();
    const cur = this.current();
    if (!sel || !cur || !cur.converged) return;
    let kind, index;
    if (sel.kind === "bus") { kind = "buses"; index = cur.buses.findIndex(b => b.bus === sel.id); }
    else { kind = sel.kind === "line" ? "lines" : "transformers"; index = sel.index; }
    this.selectedRow = { kind, index };
    $("#pf-show").value = kind;
    this.renderTableBody();
  },

  tooltip(sel) {
    const cur = this.current();
    if (!cur || !cur.converged) return elementTooltip(sel);
    const net = state.network;
    if (sel.kind === "bus") {
      const b = net.buses.find(x => x.id === sel.id);
      const r = cur.buses.find(x => x.bus === sel.id);
      if (!b || !r) return elementTooltip(sel);
      const rows = [["Voltage", `${fmt(r.vm_pu, 4)} pu`], ["Angle", `${fmt(r.va_degree, 3)}°`], ["P net generation", `${fmt(r.p_net_gen_mw, 3)} MW`], ["Q net generation", `${fmt(r.q_net_gen_mvar, 3)} MVAr`]];
      const der = net.der_units.find(d => d.bus === sel.id);
      if (der) {
        const name = `der${der.id}`;
        const gen = [...cur.external_grid, ...cur.generators, ...cur.static_generators].find(g => g.name === name);
        rows.push({ sep: `${UNIT_NAME[der.unit_type] || der.unit_type} (${der.bus_type})` });
        if (gen) rows.push(["P output", `${fmt(gen.p_mw, 3)} MW`], ["Q output", `${fmt(gen.q_mvar, 3)} MVAr`]);
      }
      const loads = cur.loads.filter(l => l.bus === sel.id);
      if (loads.length) { rows.push({ sep: "Loads" }); loads.forEach(l => rows.push([l.name || "load", `${fmt(l.p_mw, 3)} MW, ${fmt(l.q_mvar, 3)} MVAr`])); }
      return `<div class="tt-title">${esc(b.name || `bus${b.id}`)} <span class="muted">· bus ${b.id}</span></div>` + ttTable(rows);
    }
    const isLine = sel.kind === "line";
    const r = isLine ? cur.lines[sel.index] : cur.transformers[sel.index];
    if (!r) return elementTooltip(sel);
    const rows = isLine
      ? [["Buses", `${r.from_bus} → ${r.to_bus}`], ["P from / to", `${fmt(r.p_from_mw, 3)} / ${fmt(r.p_to_mw, 3)} MW`], ["Q from / to", `${fmt(r.q_from_mvar, 3)} / ${fmt(r.q_to_mvar, 3)} MVAr`],
        ["Losses", `${fmtSmart(r.pl_mw)} MW, ${fmtSmart(r.ql_mvar)} MVAr`], ["Current", `${fmtSmart(r.i_ka)} kA`], ["V from / to", `${fmt(r.vm_from_pu, 4)} / ${fmt(r.vm_to_pu, 4)} pu`]]
      : [["HV → LV bus", `${r.hv_bus} → ${r.lv_bus}`], ["P HV / LV", `${fmt(r.p_hv_mw, 3)} / ${fmt(r.p_lv_mw, 3)} MW`], ["Q HV / LV", `${fmt(r.q_hv_mvar, 3)} / ${fmt(r.q_lv_mvar, 3)} MVAr`],
        ["Losses", `${fmtSmart(r.pl_mw)} MW, ${fmtSmart(r.ql_mvar)} MVAr`], ["Current HV", `${fmtSmart(r.i_hv_ka)} kA`], ["Loading", `${fmt(r.loading_percent, 2)} %`]];
    return `<div class="tt-title">${isLine ? "Line" : "Transformer"} ${esc(r.name || `#${sel.index}`)}</div>` + ttTable(rows);
  },
};
