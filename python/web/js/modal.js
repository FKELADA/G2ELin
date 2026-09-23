// Modal analysis: one linearisation (cached per network version), seven
// sub-pages over it. Free-motion and step-response pages hold "channels"
// (one disturbance each) with any number of subplots (a set of signals each).

const MODAL_VIEWS = {
  eigenmap: { title: "Eigenvalue map", desc: "Every eigenvalue of the closed-loop state matrix on a symlog map, with the 5 % and 70.7 % damping-ratio guides. Click a mode (on the map or in the table) to select it for the other views." },
  participation: { title: "Participation heatmap", desc: "Participation factors of states (rows) in modes (columns) — which physical states drive which oscillation. Each column sums to 1." },
  single: { title: "Single-mode participation", desc: "How much every state participates in one chosen mode." },
  sensitivity: { title: "Sensitivity heatmap", desc: "Sensitivity of one eigenvalue to each entry of the state matrix (masked to structurally non-zero entries) — which couplings move the mode most." },
  shape: { title: "Mode shape", desc: "Relative phase of the most-participating states in one mode: who swings against whom." },
  free: { title: "Free-motion response", desc: "Closed-form (modal-expansion) response of the linearised model to an initial-condition offset in one state. Each channel is one perturbation; split the signals over as many subplots as you like." },
  root: { title: "Root locus", desc: "Sweep any one network parameter across a range: at every value the power flow and the eigenvalues are recomputed, and each mode's path is drawn coloured by the parameter value. The table ranks the modes the parameter moves most." },
  step: { title: "Step response", desc: "Response of chosen outputs of the linearised model to a step in one input. Each channel is one input step; split its outputs over as many subplots as you like." },
};

const ModalPage = {
  inited: false,
  view: "eigenmap",
  mode: null,
  pending: null,
  channels: { free: [], step: [] },
  channelNamesKey: null,
  // Buses, lines and loads have two states each and are rarely what a study
  // is about, so they are kept out of the pickers unless asked for.
  hideNetStates: true,

  init() {
    if (this.inited) return;
    this.inited = true;
    $("#page-modal").innerHTML = `
      <div class="page-head"><div class="crumb">Analysis · Modal analysis</div><h1 id="modal-title"></h1><p id="modal-desc"></p></div>
      <div id="modal-context"></div>
      <div class="card"><div class="controls" style="align-items:center">
        <button class="secondary" id="modal-run">Re-run modal analysis</button>
        <span class="spinner" id="modal-spinner"></span>
        <p class="status-line" id="modal-status"></p></div></div>
      <div id="modal-body" class="stack" style="margin-top:1rem"></div>`;
    $("#modal-run").addEventListener("click", () => this.load(true));
    on("network:loaded", () => {
      this.mode = null;
      if ($("#page-modal").classList.contains("active")) this.show(this.view);
    });
  },

  show(view) {
    this.view = MODAL_VIEWS[view] ? view : "eigenmap";
    $("#modal-title").textContent = MODAL_VIEWS[this.view].title;
    $("#modal-desc").textContent = MODAL_VIEWS[this.view].desc;
    $("#modal-context").innerHTML = networkContextHtml({ plot: true });
    bindContextPlot($("#modal-context"));
    this.load(false);
  },

  // Runs (or reuses) the modal analysis for the current network version.
  async load(force) {
    const body = $("#modal-body"), status = $("#modal-status"), spin = $("#modal-spinner");
    if (!state.network) { body.innerHTML = ""; status.textContent = ""; return; }
    if (!force && state.modal && state.modal.version === state.version) { this.renderStatus(); this.renderView(); return; }
    const version = state.version;
    body.innerHTML = "";
    status.textContent = "";
    setSpinner(spin, "Linearising the network (the first run derives the symbolic models and can take a while)…");
    $("#modal-run").disabled = true;
    try {
      if (!this.pending || this.pending.version !== version || force) {
        this.pending = { version, promise: netPost("modal") };
      }
      const data = await this.pending.promise;
      if (version !== state.version) return;
      state.modal = { version, data };
      this.renderStatus();
      this.renderView();
    } catch (e) {
      if (version === state.version) body.innerHTML = errorHtml(e);
      this.pending = null;
    } finally { setSpinner(spin, ""); $("#modal-run").disabled = false; }
  },

  renderStatus() {
    const m = state.modal.data;
    const ref = (m.reference_modes || []).length;
    $("#modal-status").innerHTML = `${m.n_states} states · ${m.stable ? '<span class="ok">small-signal stable</span>' : '<span class="bad">unstable</span>'} · max Re(λ) = ${m.max_real_part.toExponential(3)}`
      + (ref ? ` · <span class="muted" title="Nothing pins the absolute position of the dq frame, so the model always has a marginal direction: turn every angle by the same amount and nothing physical changes. These modes sit at the origin by construction and are left out of the verdict above.">${ref} reference-angle mode${ref > 1 ? "s" : ""} (λ ≈ 0, not counted)</span>` : "");
  },

  sortedModes() { return [...state.modal.data.modes].sort((a, b) => b.real - a.real); },
  defaultMode() {
    const modes = this.sortedModes();
    return (modes.find(m => m.imag > 1e-6) || modes[0]).mode;
  },
  modeLabel(m) { return `Mode ${m.mode}: ${fmt(m.undamped_hz, 3)} Hz, ζ = ${fmt(m.damping_pct, 2)} % (λ = ${m.real.toExponential(2)} ${m.imag >= 0 ? "+" : "−"} j${Math.abs(m.imag).toExponential(2)})`; },
  modeSelectHtml() {
    if (this.mode === null) this.mode = this.defaultMode();
    return `<div class="field mode-picker"><label for="modal-mode">Mode</label><select id="modal-mode">${this.sortedModes().map(m => `<option value="${m.mode}"${m.mode === this.mode ? " selected" : ""}>${esc(this.modeLabel(m))}</option>`).join("")}</select></div>`;
  },
  bindModeSelect(render) {
    $("#modal-mode").addEventListener("change", e => { this.mode = parseInt(e.target.value, 10); render(); });
  },

  renderView() {
    const fn = {
      eigenmap: this.viewEigenmap, participation: this.viewParticipation, single: this.viewSingle,
      sensitivity: this.viewSensitivity, shape: this.viewShape, free: () => this.viewChannels("free"), step: () => this.viewChannels("step"),
      root: () => RootLocus.render($("#modal-body")),
    }[this.view];
    fn.call(this);
  },

  // --- Eigenvalue map ---
  viewEigenmap() {
    const m = state.modal.data;
    if (this.mode === null) this.mode = this.defaultMode();
    // Categories the user has switched off. Kept across redraws of this view
    // but not across networks -- a different model has different modes.
    if (!this._hiddenCats) this._hiddenCats = new Set();
    const cats = m.categories || [];
    const counts = {};
    m.modes.forEach(md => { const c = md.category || "mixed"; counts[c] = (counts[c] || 0) + 1; });
    const present = cats.filter(c => counts[c.id]);
    const body = $("#modal-body");
    body.innerHTML = `<div class="card"><div class="card-title">Eigenvalues<span class="card-sub">— shape is what kind of mode it is, colour is how stable</span><span class="spacer"></span><button class="ghost small" id="eig-reset">Reset zoom</button></div>
        <div class="legend kind-legend" id="eig-kinds">${present.map(c => `
          <button class="kind-toggle" data-cat="${c.id}" aria-pressed="true" title="${esc(c.note)}">
            ${modeMarkerSwatch(c.id)}${esc(c.label)} <span class="muted">${counts[c.id]}</span></button>`).join("")}
          <span class="muted" style="margin-left:0.3rem">click a kind to hide it</span></div>
        <div id="eig-map"></div>
        <div class="legend"><span><span class="swatch" style="background:var(--good)"></span>Stable</span><span><span class="swatch" style="background:var(--warning)"></span>Marginal</span><span><span class="swatch" style="background:var(--critical)"></span>Unstable</span>
          <span><span class="line-swatch" style="border-color:var(--critical);border-top-style:dashed"></span>5 % damping</span><span><span class="line-swatch" style="border-color:var(--good);border-top-style:dashed"></span>70.7 % damping</span>
          <span class="muted">Scroll to zoom · drag to pan · click a mode</span></div>
        <p class="muted" id="eig-kind-note" style="font-size:0.78rem;margin:0.6rem 0 0"></p></div>
      <div class="card" id="eig-detail"></div>
      <div class="card"><div class="card-title">Modes <span class="card-sub">— sorted by real part, least stable first</span></div><div id="eig-table"></div></div>`;
    const drawMap = () => {
      const view = this._eigZoom ? this._eigZoom.get() : null;
      $("#eig-map").innerHTML = eigenvalueMapSvg(m.modes, this.mode, this._hiddenCats);
      const svg = $("#eig-map svg");
      this._eigZoom = attachSvgZoomPan(svg);
      if (view) this._eigZoom.set(view);
      $$("#eig-map .pt").forEach(pt => {
        const mode = m.modes.find(x => x.mode === +pt.dataset.mode);
        pt.addEventListener("mouseenter", evt => showTooltip(`<div class="tt-title">Mode ${mode.mode}</div>` + ttTable([["λ", `${mode.real.toExponential(3)} ${mode.imag >= 0 ? "+" : "−"} j${Math.abs(mode.imag).toExponential(3)}`], ["Frequency", `${fmt(mode.undamped_hz, 3)} Hz`], ["Damping", `${fmt(mode.damping_pct, 2)} %`], ["Top state", `${esc(mode.state1)} (${fmt(mode.part1_pct, 1)} %)`],
          ["Kind", `${esc(kindLabel(m, mode.category))} — ${fmt((mode.category_share || 0) * 100, 0)} % of the mode`]]), evt));
        pt.addEventListener("mousemove", moveTooltip);
        pt.addEventListener("mouseleave", hideTooltip);
        pt.addEventListener("click", () => { if (!svg.__justPanned) { this.mode = mode.mode; drawMap(); drawDetail(); drawTable(); } });
      });
    };
    const drawDetail = () => {
      const mode = m.modes.find(x => x.mode === this.mode);
      const st = statusOf(mode.real);
      $("#eig-detail").innerHTML = `<div class="card-title">Selected: mode ${mode.mode} <span class="badge" style="color:var(--${st.cls})">${st.label}</span></div>
        <div class="stats"><div class="stat"><span class="n">${mode.real.toExponential(3)}</span><span class="l">Real part (1/s)</span></div>
          <div class="stat"><span class="n">${fmt(mode.undamped_hz, 3)} Hz</span><span class="l">Frequency</span></div>
          <div class="stat"><span class="n">${fmt(mode.damping_pct, 2)} %</span><span class="l">Damping ratio</span></div>
          <div class="stat"><span class="n" style="font-size:0.9rem">${esc(mode.state1)} · ${esc(mode.state2)} · ${esc(mode.state3)}</span><span class="l">Top participating states</span></div>
          <div class="stat"><span class="n" style="font-size:0.9rem">${esc(kindLabel(m, mode.category))}</span><span class="l">Kind — ${fmt((mode.category_share || 0) * 100, 0)} % of the mode</span></div></div>
        <div class="controls" style="margin-top:0.9rem"><a href="#/modal/single"><button class="secondary small">Participation →</button></a><a href="#/modal/shape"><button class="secondary small">Mode shape →</button></a><a href="#/modal/sensitivity"><button class="secondary small">Sensitivity →</button></a></div>`;
    };
    const drawTable = () => {
      const rows = this.sortedModes().filter(md => !this._hiddenCats.has(md.category || "mixed")).map(md => {
        const st = statusOf(md.real);
        return `<tr class="clickable${md.mode === this.mode ? " hl" : ""}" data-mode="${md.mode}"><td>${md.mode}</td><td class="name"><span class="swatch" style="display:inline-block;width:9px;height:9px;border-radius:50%;background:var(--${st.cls});margin-right:0.4em"></span>${st.label}</td>
          <td class="name" title="${esc(kindShares(md))}">${modeMarkerSwatch(md.category || "mixed")}${esc(kindLabel(m, md.category))}</td>
          <td>${md.real.toExponential(3)}</td><td>${md.imag.toExponential(3)}</td><td>${fmt(md.undamped_hz, 3)}</td><td>${fmt(md.damping_pct, 2)}</td>
          <td class="name">${esc(md.state1)} (${fmt(md.part1_pct, 1)}%)</td><td class="name">${esc(md.state2)} (${fmt(md.part2_pct, 1)}%)</td><td class="name">${esc(md.state3)} (${fmt(md.part3_pct, 1)}%)</td></tr>`;
      }).join("");
      $("#eig-table").innerHTML = `<div class="tablewrap"><table><thead><tr><th>Mode</th><th>Status</th><th>Kind</th><th>Real</th><th>Imag</th><th>Freq (Hz)</th><th>Damping (%)</th><th>Top state</th><th>2nd</th><th>3rd</th></tr></thead><tbody>${rows}</tbody></table></div>`;
      $$("#eig-table tr[data-mode]").forEach(tr => tr.addEventListener("click", () => { this.mode = +tr.dataset.mode; drawMap(); drawDetail(); drawTable(); }));
    };
    const syncLegend = () => {
      $$("#eig-kinds .kind-toggle").forEach(b => {
        const off = this._hiddenCats.has(b.dataset.cat);
        b.setAttribute("aria-pressed", String(!off));
        b.classList.toggle("off", off);
      });
      // With one kind left showing, say what it is -- that is the moment the
      // explanation is wanted, and the moment there is room for it.
      const shown = present.filter(c => !this._hiddenCats.has(c.id));
      $("#eig-kind-note").textContent = shown.length === 1 ? shown[0].note : "";
    };
    $$("#eig-kinds .kind-toggle").forEach(btn => btn.addEventListener("click", () => {
      const c = btn.dataset.cat;
      if (this._hiddenCats.has(c)) this._hiddenCats.delete(c); else this._hiddenCats.add(c);
      syncLegend(); drawMap(); drawTable();
    }));
    this._eigZoom = null;
    drawMap(); drawDetail(); drawTable(); syncLegend();
    $("#eig-reset").addEventListener("click", () => this._eigZoom && this._eigZoom.reset());
  },

  // --- Participation heatmap ---
  viewParticipation() {
    const m = state.modal.data;
    const body = $("#modal-body");
    const nS = m.state_names.length, nM = m.modes.length;
    body.innerHTML = `<div class="card"><div class="controls">
        <div class="field"><label for="pp-filter">States containing</label><input type="text" id="pp-filter" placeholder="e.g. SM_1 or dw_r"></div>
        <div class="field"><label for="pp-s0">States from</label><input type="number" id="pp-s0" value="0" min="0" max="${nS - 1}"></div>
        <div class="field"><label for="pp-s1">to</label><input type="number" id="pp-s1" value="${Math.min(30, nS) - 1}" min="0" max="${nS - 1}"></div>
        <div class="field"><label for="pp-m0">Modes from</label><input type="number" id="pp-m0" value="0" min="0" max="${nM - 1}"></div>
        <div class="field"><label for="pp-m1">to</label><input type="number" id="pp-m1" value="${Math.min(20, nM) - 1}" min="0" max="${nM - 1}"></div>
        <label class="check" style="margin-bottom:0.5rem"><input type="checkbox" id="pp-sort"> Modes least-stable first</label>
      </div></div><div class="card" id="pp-out"></div>`;
    const draw = () => {
      const q = $("#pp-filter").value.trim().toLowerCase();
      let rows = m.state_names.map((n, i) => i);
      if (q) rows = rows.filter(i => m.state_names[i].toLowerCase().includes(q));
      const s0 = +$("#pp-s0").value || 0, s1 = +$("#pp-s1").value;
      if (!q) rows = rows.filter(i => i >= s0 && i <= s1);
      let cols = m.modes.map(md => md.mode).sort((a, b) => a - b);
      if ($("#pp-sort").checked) cols = this.sortedModes().map(md => md.mode);
      const m0 = +$("#pp-m0").value || 0, m1 = +$("#pp-m1").value;
      cols = cols.slice(m0, m1 + 1);
      if (!rows.length || !cols.length) { $("#pp-out").innerHTML = `<p class="empty">Nothing in range.</p>`; return; }
      $("#pp-out").innerHTML = `<p class="status-line" style="margin-bottom:0.6rem">${rows.length} states × ${cols.length} modes</p>`
        + heatmapTable(rows.map(i => cols.map(c => m.participation[i][c])), rows.map(i => m.state_names[i]), cols.map(c => `M${c}`));
    };
    $$("#modal-body input").forEach(i => i.addEventListener(i.type === "checkbox" ? "change" : "input", draw));
    draw();
  },

  // --- Single-mode participation ---
  viewSingle() {
    const m = state.modal.data;
    $("#modal-body").innerHTML = `<div class="card"><div class="controls">${this.modeSelectHtml()}
      <div class="field"><label for="sp-top">Show</label><select id="sp-top"><option value="15">top 15 states</option><option value="30">top 30 states</option><option value="all">all states</option></select></div></div></div>
      <div class="card" id="sp-out"></div>`;
    const draw = () => {
      const idx = m.state_names.map((_, i) => i).sort((a, b) => m.participation[b][this.mode] - m.participation[a][this.mode]);
      const top = $("#sp-top").value;
      const shown = top === "all" ? m.state_names.map((_, i) => i) : idx.slice(0, +top);
      $("#sp-out").innerHTML = `<p class="status-line" style="margin-bottom:0.6rem">Participation of ${top === "all" ? "every state" : `the top ${shown.length} states`} in mode ${this.mode}</p>`
        + barChartHtml(shown.map(i => m.state_names[i]), shown.map(i => m.participation[i][this.mode]));
    };
    this.bindModeSelect(draw);
    $("#sp-top").addEventListener("change", draw);
    draw();
  },

  // --- Sensitivity ---
  viewSensitivity() {
    $("#modal-body").innerHTML = `<div class="card"><div class="controls">${this.modeSelectHtml()}<span class="spinner" id="sens-spin"></span></div></div><div id="sens-out" class="stack"></div>`;
    const draw = async () => {
      const out = $("#sens-out"), mode = this.mode;
      out.innerHTML = "";
      setSpinner($("#sens-spin"), "Computing…");
      try {
        const r = await netPost("modal/sensitivity", { mode });
        if (mode !== this.mode) return;
        const n = r.state_names.length;
        // Only rows/columns that hold any sensitivity at all -- the full n x n
        // matrix is mostly structural zeros.
        const keep = r.state_names.map((_, i) => i).filter(i => r.matrix[i].some(v => v > 0) || r.matrix.some(row => row[i] > 0));
        const top = r.top.map(e => `<tr><td class="name">${esc(e.row_state)}</td><td class="name">${esc(e.col_state)}</td><td>${fmtSmart(e.value)}</td></tr>`).join("");
        out.innerHTML = `<div class="card"><div class="card-title">Largest sensitivities, mode ${mode}</div><div class="tablewrap"><table><thead><tr><th>Row state (∂f of)</th><th>Column state (w.r.t.)</th><th>|∂λ/∂A|</th></tr></thead><tbody>${top}</tbody></table></div></div>
          <div class="card"><div class="card-title">Sensitivity matrix <span class="card-sub">— ${keep.length <= 40 ? `${keep.length} of ${n} states with non-zero entries` : `${keep.length} states have non-zero entries; showing the 40 with the largest`}</span></div><div id="sens-heat"></div></div>`;
        let shown = keep;
        if (keep.length > 40) {
          const score = i => Math.max(...r.matrix[i], ...r.matrix.map(row => row[i]));
          shown = [...keep].sort((a, b) => score(b) - score(a)).slice(0, 40).sort((a, b) => a - b);
        }
        $("#sens-heat").innerHTML = heatmapTable(shown.map(i => shown.map(j => r.matrix[i][j])), shown.map(i => r.state_names[i]), shown.map(i => r.state_names[i]));
      } catch (e) { out.innerHTML = errorHtml(e); }
      finally { setSpinner($("#sens-spin"), ""); }
    };
    this.bindModeSelect(draw);
    draw();
  },

  // --- Mode shape ---
  viewShape() {
    $("#modal-body").innerHTML = `<div class="card"><div class="controls">${this.modeSelectHtml()}<span class="spinner" id="shape-spin"></span></div></div><div class="card" id="shape-out"></div>`;
    const draw = async () => {
      const out = $("#shape-out"), mode = this.mode;
      setSpinner($("#shape-spin"), "Computing…");
      try {
        const r = await netPost("modal/mode_shape", { mode });
        if (mode !== this.mode) return;
        out.innerHTML = `<div class="card-title">Mode ${mode} — relative phase of the 5 most-participating states</div>${polarChartHtml(r.states, r.angles_deg)}
          <div class="legend" style="justify-content:center">${r.states.map((s, i) => `<span><span class="swatch" style="background:${seriesColor(i)}"></span>${esc(s)} (${fmt(r.angles_deg[i], 1)}°)</span>`).join("")}</div>`;
      } catch (e) { out.innerHTML = errorHtml(e); }
      finally { setSpinner($("#shape-spin"), ""); }
    };
    this.bindModeSelect(draw);
    draw();
  },

  // --- Channels (free motion / step response) ---
  resetChannelsIfNeeded() {
    const m = state.modal.data;
    const key = JSON.stringify([m.state_names, m.input_names, m.output_names]);
    if (key === this.channelNamesKey) return;
    this.channelNamesKey = key;
    this.channels = { free: [this.newChannel("free")], step: [this.newChannel("step")] };
  },

  newChannel(kind) {
    const m = state.modal.data;
    this._chId = (this._chId || 0) + 1;
    if (kind === "free") {
      const dw = m.state_names.filter(n => n.includes("dw_r"));
      const rotor = m.state_names.filter(n => n.includes("dw_r") || n.startsWith("theta"));
      return { id: this._chId, perturb: (dw[0] || rotor[0] || m.state_names[0]), offset: 0.05, tFinal: 2.0, subplots: [dw.length ? dw : m.state_names.slice(0, 3)], result: null };
    }
    const input = m.input_names.find(n => n.startsWith("P_ref")) || m.input_names[0];
    const suffix = input && input.includes("_{") ? input.slice(input.indexOf("_{")) : "";
    const outs = m.output_names.filter(n => suffix && n.endsWith(suffix));
    const pick = (re) => outs.filter(n => re.test(n));
    const sp1 = pick(/^p_e|^P_|^p_/).slice(0, 2), sp2 = pick(/^w|^dw|^omega/).slice(0, 2);
    const subplots = [sp1.length ? sp1 : (outs.slice(0, 2).length ? outs.slice(0, 2) : m.output_names.slice(0, 1))];
    if (sp2.length) subplots.push(sp2);
    return { id: this._chId, input, amplitude: 0.1, tFinal: 2.0, subplots, result: null };
  },

  // The states a channel offers: the units' own unless the network's are asked
  // for. Whatever is already plotted or perturbed is always kept.
  channelStates(keep = []) {
    const names = state.modal.data.state_names;
    if (!this.hideNetStates) return names;
    const k = new Set(keep);
    return names.filter(n => k.has(n) || !isNetworkElementSignal(n));
  },

  viewChannels(kind) {
    this.resetChannelsIfNeeded();
    const body = $("#modal-body");
    body.innerHTML = `<div class="card"><div class="controls" style="align-items:center">
        <button id="ch-add">+ Add channel</button><button class="secondary" id="ch-plot-all">Plot all channels</button>
        ${kind === "free" ? `<label class="check" title="Buses, lines and loads have a lot of states (voltages and currents) that are rarely the point of a study."><input type="checkbox" id="ch-hide-net"${this.hideNetStates ? " checked" : ""}> Hide bus, line and load states</label>` : ""}
        <span class="muted" style="font-size:0.8rem">A channel = one ${kind === "free" ? "initial-condition perturbation" : "input step"}; each subplot shows its own set of ${kind === "free" ? "states" : "outputs"} (crosshairs are synced within a channel).</span></div></div>
      <div id="ch-list"></div>`;
    $("#ch-hide-net")?.addEventListener("change", e => { this.hideNetStates = e.target.checked; this.renderChannels(kind); });
    $("#ch-add").addEventListener("click", () => {
      const last = this.channels[kind][this.channels[kind].length - 1];
      const ch = this.newChannel(kind);
      if (last) {
        ch.subplots = last.subplots.map(s => [...s]);
        // Default the new channel to the next state/input after the last one's.
        const m = state.modal.data, names = kind === "free" ? m.state_names : m.input_names;
        const key = kind === "free" ? "perturb" : "input";
        ch[key] = names[(names.indexOf(last[key]) + 1) % names.length];
      }
      this.channels[kind].push(ch);
      this.renderChannels(kind);
      this.plotChannel(kind, ch);
    });
    $("#ch-plot-all").addEventListener("click", () => this.channels[kind].forEach(ch => this.plotChannel(kind, ch)));
    this.renderChannels(kind);
    this.channels[kind].forEach(ch => { if (!ch.result) this.plotChannel(kind, ch); });
  },

  renderChannels(kind) {
    const list = $("#ch-list");
    list.innerHTML = "";
    this.channels[kind].forEach((ch, ci) => list.appendChild(this.channelEl(kind, ch, ci)));
  },

  channelEl(kind, ch, ci) {
    const m = state.modal.data;
    const color = seriesColor(ci);
    const node = el(`<div class="channel" data-ch="${ch.id}">
      <div class="channel-head"><span class="ch-dot" style="background:${color}"></span><h3>Channel ${ci + 1}</h3><span class="muted" style="font-size:0.8rem" data-role="summary"></span><span class="spacer"></span>
        <span class="spinner" data-role="spin"></span><button class="small" data-role="plot">Plot</button>${this.channels[kind].length > 1 ? `<button class="ghost small" data-role="remove">Remove</button>` : ""}</div>
      <div class="channel-body">
        <div class="controls">
          <div class="field"><label>Element type</label><select data-role="el-kind"></select></div>
          <div class="field"><label>Element</label><select data-role="el" style="min-width:210px"></select></div>
          ${kind === "free"
          ? `<div class="field"><label>Perturbed state</label><select data-role="perturb" style="min-width:210px"></select></div>
             <div class="field"><label>Offset</label><input type="number" step="any" data-role="amp" value="${ch.offset}"></div>`
          : `<div class="field"><label>Stepped input</label><select data-role="perturb" style="min-width:210px"></select></div>
             <div class="field"><label>Step amplitude</label><input type="number" step="any" data-role="amp" value="${ch.amplitude}"></div>`}
          <div class="field"><label>Duration (s, max 20)</label><input type="number" step="0.1" min="0.01" max="20" data-role="tf" value="${ch.tFinal}"></div>
        </div>
        <div style="margin-top:0.9rem" data-role="subplots"></div>
        <button class="secondary small" style="margin-top:0.6rem" data-role="add-sp">+ Add subplot</button>
        <div class="channel-plots" data-role="plots"></div>
      </div></div>`);
    const q = r => node.querySelector(`[data-role="${r}"]`);
    const key = kind === "free" ? "perturb" : "input";
    // What can be perturbed/stepped, and the elements those signals belong to.
    const sigNames = kind === "free" ? this.channelStates([ch.perturb]) : m.input_names;
    const names = kind === "free" ? this.channelStates(ch.subplots.flat()) : m.output_names;
    const elKeys = new Set();
    sigNames.forEach(n => { const e = signalElement(n); if (e) elKeys.add(e.key); });
    // Element type -> element -> signal, so the list stays short and readable.
    const fillPerturb = () => {
      const list = sigNames.filter(n => {
        const e = signalElement(n);
        if (ch.elKey) return e && e.key === ch.elKey;
        if (ch.elKind) return e && e.kind === ch.elKind;
        return true;
      });
      q("perturb").innerHTML = list.length
        ? list.map(n => `<option${n === ch[key] ? " selected" : ""}>${esc(n)}</option>`).join("")
        : `<option value="">(none here)</option>`;
      q("perturb").disabled = !list.length;
      if (list.length && !list.includes(ch[key])) { ch[key] = list[0]; q("perturb").value = list[0]; }
    };
    const fillElements = () => {
      q("el-kind").innerHTML = elementKindOptionsHtml(ch.elKind || "", elKeys);
      ch.elKind = q("el-kind").value;
      q("el").innerHTML = elementOptionsHtml(ch.elKey || "", { kind: ch.elKind, keys: elKeys });
      ch.elKey = q("el").value;
      fillPerturb();
    };
    // Start on the element of whatever the channel already perturbs.
    if (ch.elKey === undefined) { const e = signalElement(ch[key]); ch.elKind = e ? e.kind : ""; ch.elKey = e ? e.key : ""; }
    fillElements();
    const updateSummary = () => { q("summary").textContent = kind === "free" ? `— ${ch.perturb} + ${ch.offset}` : `— step ${ch.amplitude} in ${ch.input}`; };
    updateSummary();
    q("el-kind").addEventListener("change", e => { ch.elKind = e.target.value; ch.elKey = ""; fillElements(); updateSummary(); });
    q("el").addEventListener("change", e => { ch.elKey = e.target.value; fillPerturb(); updateSummary(); });
    q("perturb").addEventListener("change", e => { if (kind === "free") ch.perturb = e.target.value; else ch.input = e.target.value; updateSummary(); });
    q("amp").addEventListener("input", e => { const v = parseFloat(e.target.value); if (Number.isFinite(v)) { if (kind === "free") ch.offset = v; else ch.amplitude = v; updateSummary(); } });
    q("tf").addEventListener("input", e => { const v = parseFloat(e.target.value); if (Number.isFinite(v)) ch.tFinal = v; });
    q("plot").addEventListener("click", () => this.plotChannel(kind, ch));
    if (q("remove")) q("remove").addEventListener("click", () => { this.channels[kind] = this.channels[kind].filter(c => c !== ch); this.renderChannels(kind); });
    const drawSubplots = () => {
      const box = q("subplots");
      box.innerHTML = "";
      ch.subplots.forEach((sp, si) => {
        const row = el(`<div class="subplot-row"><span class="sp-tag">Subplot ${si + 1}</span><div class="sp-picker"></div>${ch.subplots.length > 1 ? `<button class="ghost small">Remove</button>` : ""}</div>`);
        new SignalPicker(row.querySelector(".sp-picker"), {
          options: names, selected: sp, colors: true, byElement: true,
          placeholder: kind === "free" ? "add state" : "add output",
          onChange: list => { ch.subplots[si] = list; },
        });
        const rm = row.querySelector("button.ghost");
        if (rm) rm.addEventListener("click", () => { ch.subplots.splice(si, 1); drawSubplots(); if (ch.result) this.drawChannelPlots(kind, ch, node); });
        box.appendChild(row);
      });
    };
    drawSubplots();
    q("add-sp").addEventListener("click", () => { ch.subplots.push([]); drawSubplots(); });
    if (ch.result) this.drawChannelPlots(kind, ch, node);
    return node;
  },

  async plotChannel(kind, ch) {
    const node = $(`#ch-list [data-ch="${ch.id}"]`);
    if (!node) return;
    const spin = node.querySelector('[data-role="spin"]'), plots = node.querySelector('[data-role="plots"]');
    const union = [...new Set(ch.subplots.flat())];
    if (!union.length) { plots.innerHTML = `<p class="empty">Add at least one signal to a subplot.</p>`; return; }
    setSpinner(spin, "Computing…");
    try {
      const r = kind === "free"
        ? await netPost("modal/free_response", { perturb_state: ch.perturb, offset: ch.offset, t_final: ch.tFinal, plot_states: union })
        : await netPost("modal/step_response", { input_name: ch.input, output_names: union, amplitude: ch.amplitude, t_final: ch.tFinal });
      ch.result = { t: r.t, series: r.series, label: kind === "free" ? `${ch.perturb} + ${ch.offset}` : `step ${ch.amplitude} in ${ch.input}` };
      this.drawChannelPlots(kind, ch, node);
    } catch (e) { plots.innerHTML = errorHtml(e); }
    finally { setSpinner(spin, ""); }
  },

  drawChannelPlots(kind, ch, node) {
    const plots = node.querySelector('[data-role="plots"]');
    plots.innerHTML = `<p class="status-line" style="margin-bottom:0.3rem">Linearised response to ${esc(ch.result.label)}${kind === "free" ? " (deviation from equilibrium)" : " (deviation from the operating point)"}</p>`;
    const group = {};
    const wrap = el(`<div class="subplots"></div>`);
    ch.subplots.forEach((sp, si) => {
      const names = sp.filter(n => ch.result.series[n]);
      if (!names.length) return;
      const series = names.map((n, i) => ({ name: n, t: ch.result.t, y: ch.result.series[n], color: seriesColor(i) }));
      const chart = lineChart(series, { title: `Subplot ${si + 1}`, group, height: 220, xLabel: "t" });
      chart.insertAdjacentHTML("beforeend", legendHtml(series));
      wrap.appendChild(chart);
    });
    if (!wrap.children.length) wrap.innerHTML = `<p class="empty">The selected signals changed since this was computed — press Plot.</p>`;
    plots.appendChild(wrap);
  },
};


// --- Mode kinds (see g2elin_core.modal.classify) ------------------------------
// The labels and explanations come from the modal response, so nothing here
// hard-codes a category; these are only the lookup and the tooltip text.
function kindLabel(modalData, categoryId) {
  const c = (modalData.categories || []).find(x => x.id === categoryId);
  return c ? c.label : (categoryId || "Mixed");
}

function kindShares(mode) {
  const shares = mode.category_shares || {};
  return Object.entries(shares)
    .filter(([, v]) => v > 0.005)
    .sort((a, b) => b[1] - a[1])
    .map(([k, v]) => `${k}: ${(v * 100).toFixed(0)}%`)
    .join(" · ");
}
