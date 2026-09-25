// Network page: choose a preset (plotted immediately) or build from scratch,
// edit on the canvas, and edit every element parameter in the inspector that
// opens on the right when an element is clicked.

const FIELD_DEFS = {
  network: [
    { key: "name", label: "Network name", type: "text", full: true },
    { key: "f_hz", label: "Frequency", unit: "Hz", type: "number" },
    { key: "sn_mva", label: "Base power", unit: "MVA", type: "number", help: "All per-unit values are on this base." },
    { key: "units_use_first_transformer", label: "MATLAB-compatible transformers", type: "bool", full: true,
      help: "Off (default): each unit's dynamic model uses its own transformer — the one the power flow uses. On: every unit uses the network's first transformer, as the MATLAB tool does (to reproduce its results)." },
    { key: "nodes_share_first_line_b", label: "MATLAB-compatible bus capacitance", type: "bool", full: true,
      help: "Off (default): each bus's dynamic model uses its own capacitance — half the charging of the lines meeting at it, plus its capacitor banks. On: every bus uses the first line's charging, whatever is connected to it, as the MATLAB tool does. The ported presets have it on so they still reproduce that tool's numbers." },
    { key: "min_node_b_pu", label: "Minimum bus capacitance", unit: "pu", type: "number", nullable: true, full: true,
      help: "Only used for a bus whose lines declare no charging at all, as distribution-feeder data often does. Without it such a network is refused rather than run on an invented number, since every bus's dynamic model divides by its capacitance." },
  ],
  bus: [
    { key: "id", label: "Bus id", type: "readonly" },
    { key: "name", label: "Name", type: "text" },
    { key: "vn_kv", label: "Nominal voltage", unit: "kV", type: "number" },
  ],
  line: [
    { key: "from_bus", label: "From bus", type: "bus" },
    { key: "to_bus", label: "To bus", type: "bus" },
    { key: "r_pu", label: "Resistance R", unit: "pu", type: "number" },
    { key: "x_pu", label: "Reactance X", unit: "pu", type: "number" },
    { key: "b_pu", label: "Shunt susceptance B", unit: "pu", type: "number", help: "Total line charging." },
    { key: "length_km", label: "Length", unit: "km", type: "number", help: "R, X and B are the whole line's; changing the length rescales them at constant Ω/km. It has no other effect on the results." },
    { key: "name", label: "Name", type: "text", full: true },
    { key: "from_closed", label: "Breaker at from-bus closed", type: "bool", breaker: true },
    { key: "to_closed", label: "Breaker at to-bus closed", type: "bool", breaker: true, help: "Either breaker open takes the line out of service." },
  ],
  transformer: [
    { key: "hv_bus", label: "HV bus", type: "bus" },
    { key: "lv_bus", label: "LV bus", type: "bus", help: "A unit's own terminal bus must be the LV side." },
    { key: "r_pu", label: "Resistance R", unit: "pu", type: "number", help: "Per unit of the transformer's rating. Also the Rt of the unit on its LV side (modal analysis & EMT)." },
    { key: "x_pu", label: "Reactance X", unit: "pu", type: "number", help: "Also the Lt of the unit on its LV side (modal analysis & EMT)." },
    { key: "sn_mva", label: "Rating", unit: "MVA", type: "number" },
    { key: "tap_ratio", label: "Tap ratio", unit: "pu", type: "number", branchOnly: true,
      help: "Off-nominal turns ratio on the HV side: v_HV = tap × v_LV at no load. 1.0 is the ratio the two buses' nominal voltages already imply. Only for a transformer between two grid buses — a unit's step-up ratio is part of that unit's own model." },
    { key: "shift_degree", label: "Phase shift", unit: "°", type: "number", branchOnly: true,
      help: "Phase-shifting transformer: the angle from HV to LV. A few degrees moves a lot of active power. As for the tap, only for a transformer between two grid buses." },
    { key: "name", label: "Name", type: "text" },
    { key: "hv_closed", label: "Breaker at HV side closed", type: "bool", breaker: true },
    { key: "lv_closed", label: "Breaker at LV side closed", type: "bool", breaker: true, help: "Either breaker open takes the transformer out of service — for a unit's transformer, the unit too." },
  ],
  load: [
    { key: "p_mw", label: "Active power P", unit: "MW", type: "number" },
    { key: "q_mvar", label: "Reactive power Q", unit: "MVAr", type: "number" },
    { key: "name", label: "Name", type: "text", full: true },
    { key: "closed", label: "Breaker closed", type: "bool", breaker: true, full: true, help: "Open = the load is disconnected from its bus." },
  ],
  shunt: [
    { key: "q_mvar", label: "Reactive power Q", unit: "MVAr", type: "number",
      help: "At nominal voltage. Negative = a capacitor bank (generates reactive power, raises the voltage); positive = a reactor (absorbs it, lowers the voltage). Same sign convention as a load." },
    { key: "r_pu", label: "Reactor resistance R", unit: "pu", type: "number", nullable: true,
      help: "A reactor's series resistance. Empty derives one from X/R = 50, a typical shunt reactor; exactly 0 leaves its resonance with the bus undamped. A capacitor bank has no branch of its own, so this does nothing for one." },
    { key: "name", label: "Name", type: "text", full: true },
    { key: "closed", label: "Breaker closed", type: "bool", breaker: true, full: true, help: "Open = the device is disconnected from its bus." },
  ],
  der: [
    { key: "id", label: "Unit id", type: "readonly" },
    { key: "unit_type", label: "Unit type", type: "select", options: ["sm", "gfm", "gfl", "infinite_bus"] },
    { key: "bus_type", label: "Load-flow bus type", type: "select", options: ["slack", "pv", "pq"], help: "Exactly one unit must be the slack." },
    { key: "v_set_pu", label: "Voltage setpoint", unit: "pu", type: "number" },
    { key: "angle_set_deg", label: "Voltage angle setpoint", unit: "\u00b0", type: "number",
      help: "The reference angle, on the slack unit only \u2014 every other bus angle is solved against it. Moving it rotates every angle together and changes nothing physical: no power flow, no current, no eigenvalue. Set it to read bus angles against a published reference (Kundur's two-area case puts G1 at 20.2\u00b0)." },
    { key: "p_set_mw", label: "Active power setpoint", unit: "MW", type: "number" },
    { key: "q_set_mvar", label: "Reactive power setpoint", unit: "MVAr", type: "number" },
    { key: "p_cons_mw", label: "Auxiliary load P", unit: "MW", type: "number" },
    { key: "q_cons_mvar", label: "Auxiliary load Q", unit: "MVAr", type: "number" },
    { key: "xd_pu", label: "Transient reactance Xd", unit: "pu", type: "number", nullable: true, only: ["sm"], help: "Used for SCR calculations only." },
    { key: "sn_mva", label: "Unit rating", unit: "MVA", type: "number", nullable: true,
      help: "The base this unit's parameter overrides are given on. Published machine data is per unit of the machine's own rating: set this to 900 to type a 900 MVA machine's x_d = 1.8 and H = 6.5 s in as published, and they become 0.2 pu and 58.5 s on a 100 MVA network. Leave empty when the overrides are already on the network base. Gains and time constants are taken as given either way." },
    { key: "closed", label: "Unit breaker closed", type: "bool", breaker: true, full: true, help: "Open = the unit is disconnected (with its transformer). The slack unit's breaker can't be opened." },
  ],
};

// Table view: same schema, one table per element kind.
const TABLE_KINDS = {
  buses: { label: "Buses", defs: FIELD_DEFS.bus },
  lines: { label: "Lines", defs: FIELD_DEFS.line },
  transformers: { label: "Transformers", defs: FIELD_DEFS.transformer },
  loads: { label: "Loads", defs: [{ key: "bus", label: "Bus", type: "bus" }, ...FIELD_DEFS.load] },
  shunts: { label: "Shunts", defs: [{ key: "bus", label: "Bus", type: "bus" }, ...FIELD_DEFS.shunt] },
  der_units: { label: "DER units", defs: [FIELD_DEFS.der[0], { key: "bus", label: "Bus", type: "bus" }, ...FIELD_DEFS.der.slice(1)] },
};

const TARGET_KIND = { net: "network", bus: "bus", der: "der", load: "load", shunt: "shunt", line: "line", transformer: "transformer" };

// A capacitor bank or a reactor, from its sign (schema.Shunt).
const shuntKindLabel = s => (Number(s.q_mvar) > 0 ? "Reactor" : "Capacitor bank");

const NetworkPage = {
  view: null,
  selected: null,
  inited: false,

  init() {
    if (this.inited) return;
    this.inited = true;
    const page = $("#page-network");
    page.innerHTML = `
      <div class="page-head"><div class="crumb">Workspace</div><h1>Network</h1>
        <p>Pick a preset or start from scratch. Drag buses to arrange the diagram, drop new elements from the palette, and click any element to edit all of its parameters on the right.</p></div>
      <div class="card">
        <div class="net-toolbar">
          <div class="field"><label for="net-preset">Preset network</label><select id="net-preset" style="min-width:300px"></select></div>
          <button class="secondary" id="net-new">New empty network</button>
          <span class="sep"></span>
          <button class="secondary" id="net-import">Import JSON</button>
          <button class="secondary" id="net-export">Export JSON</button>
          <input type="file" id="net-import-file" accept="application/json,.json" hidden>
          <span style="flex:1"></span>
          <span id="net-modified"></span>
          <button class="ghost" id="net-reset" style="display:none">Reset to preset</button>
        </div>
        <div id="net-summary" style="margin-top:1rem"></div>
      </div>
      <div class="card" id="net-saved"></div>
      <div class="workspace" id="net-workspace">
        <div class="card">
          <div class="canvas-tools">
            <div class="palette"><span class="palette-label">Drag to add</span>
              <span class="chip" data-chip="bus"><span class="dot" style="background:#77766f"></span>Bus</span>
              <span class="chip" data-chip="load"><span class="dot" style="background:var(--text-secondary);border-radius:1px"></span>Load</span>
              ${Object.keys(UNIT_LABEL).map(u => `<span class="chip" data-chip="${u}" title="${UNIT_NAME[u]}">${unitGlyphHtml(u, UNIT_COLOR[u], 15)}${UNIT_LABEL[u]}</span>`).join("")}
            </div>
            <div class="btn-group" role="group" aria-label="Canvas mode">
              <button class="secondary toggle small" data-mode="select" aria-pressed="true">Select / move</button>
              <button class="secondary toggle small" data-mode="wire-line" aria-pressed="false">Draw line</button>
              <button class="secondary toggle small" data-mode="wire-transformer" aria-pressed="false">Draw transformer</button>
            </div>
            <span style="flex:1"></span>
            <label class="check"><input type="checkbox" id="net-labels" checked> Labels</label>
            <button class="ghost small" id="net-autolayout">Auto-layout</button>
            <button class="ghost small" id="net-settings">Network settings</button>
          </div>
          <div id="net-canvas"></div>
          <p class="status-line" id="net-mode-hint" style="margin-top:0.5rem"></p>
          ${unitLegendHtml()}
        </div>
        <aside class="card inspector" id="net-inspector" style="display:none"></aside>
      </div>
      <div class="card"><div class="card-title">Model order <span class="card-sub">which dynamics the models keep — this is what makes a time-domain run EMT or RMS</span></div><div id="net-model-order"><p class="empty">No network loaded.</p></div></div>
      <div class="card"><div class="card-title">Validation <span class="card-sub" id="net-valid-sub"></span></div><div id="net-validation"><p class="empty">No network loaded.</p></div></div>
      <details class="card" id="net-tables"><summary>Table view — bulk edit every element</summary><div id="net-tables-body"></div></details>`;

    this.view = new NetworkView($("#net-canvas"), {
      editable: true,
      tooltip: elementTooltip,
      onSelect: sel => this.openInspector(sel),
      onWire: (mode, a, b) => this.commitWire(mode, a, b),
      onBreaker: br => {
        if (!toggleBreaker(br)) return;
        this.view.render();
        if (this.selected) this.openInspector(this.selected);
        if ($("#net-tables").open) this.renderTables();
      },
    });
    $$("#page-network [data-mode]").forEach(btn => btn.addEventListener("click", () => this.setMode(btn.dataset.mode)));
    this.setMode("select");
    $$("#page-network .chip").forEach(chip => this.attachChip(chip));
    $("#net-labels").addEventListener("change", e => { this.view.opts.showLabels = e.target.checked; this.view.render(); });
    $("#net-autolayout").addEventListener("click", async () => { try { await autoLayout(); this.view.fit(); } catch (e) { alert(e.message); } });
    $("#net-settings").addEventListener("click", () => this.openInspector({ kind: "network" }));
    $("#net-new").addEventListener("click", () => {
      if (state.network && !confirm("Start a new empty network? The current one is not saved.")) return;
      setNetwork({ name: "custom_network", f_hz: 50.0, sn_mva: 100.0, buses: [], lines: [], transformers: [], loads: [], shunts: [], der_units: [] }, { label: "New network" });
      $("#net-preset").value = "";
    });
    $("#net-reset").addEventListener("click", () => { if (state.presetId) this.loadPreset(state.presetId); });
    $("#net-export").addEventListener("click", () => this.exportJson());
    $("#net-import").addEventListener("click", () => $("#net-import-file").click());
    $("#net-import-file").addEventListener("change", e => this.importJson(e.target.files[0]));
    $("#net-preset").addEventListener("change", e => { if (e.target.value) this.loadPreset(e.target.value); });
    $("#net-tables").addEventListener("toggle", () => { if ($("#net-tables").open) this.renderTables(); });
    this.bindTableEvents();
    this.fillPresetSelect();
    SavedNetworks.mount();

    on("network:loaded", () => { this.closeInspector(); this.refreshAll(); });
    on("network:changed", () => ModelOrder.refreshSummary());
    on("network:layout", () => { this.view.render(); this.view.fit(); });
    on("network:changed", () => { this.renderSummary(); this.scheduleValidation(); });
    on("network:derived", () => this.refreshControlParams());
    on("network:changed", () => this.refreshControlParams());
    on("network:changed", () => this.refreshUnitModelOrder());
    this.refreshAll();
  },

  onShow() { this.view.render(); if (!this._fitted && state.network) { this.view.fit(); this._fitted = true; } },

  fillPresetSelect() {
    const sel = $("#net-preset");
    if (!sel) return;
    const groups = new Map();
    state.presets.forEach(p => { const f = presetFamily(p.id); if (!groups.has(f)) groups.set(f, []); groups.get(f).push(p); });
    sel.innerHTML = `<option value="">— choose a preset —</option>` + [...groups].map(([f, ps]) =>
      `<optgroup label="${esc(f)}">${ps.map(p => `<option value="${esc(p.id)}">${esc(p.name)}</option>`).join("")}</optgroup>`).join("");
    if (state.presetId) sel.value = state.presetId;
  },

  async loadPreset(id) {
    try {
      const net = await api(`/api/presets/${encodeURIComponent(id)}/network`);
      const p = state.presets.find(x => x.id === id);
      setNetwork(net, { label: p ? p.name : id, presetId: id });
      if ($("#net-preset")) $("#net-preset").value = id;
    } catch (e) { alert(`Could not load preset: ${e.message}`); }
  },

  refreshAll() {
    this.renderSummary();
    this.view.render();
    this.scheduleValidation(true);
    // Re-mounted rather than refreshed: which unit types a network has
    // decides which model-order sections exist at all.
    ModelOrder.mount($("#net-model-order")).catch(e => { $("#net-model-order").innerHTML = errorHtml(e); });
    if ($("#net-tables")?.open) this.renderTables();
  },

  renderSummary() {
    const box = $("#net-summary");
    if (!box) return;
    const n = state.network;
    $("#net-reset").style.display = state.presetId && isModified() ? "" : "none";
    $("#net-modified").innerHTML = !n ? "" : (state.presetId ? (isModified() ? '<span class="badge warn">modified from preset</span>' : '<span class="badge accent">preset</span>') : '<span class="badge">custom network</span>');
    if (!n) { box.innerHTML = `<p class="empty">No network loaded — choose a preset above or start a new one.</p>`; return; }
    const p = state.presets.find(x => x.id === state.presetId);
    const types = [...new Set(n.der_units.map(d => UNIT_LABEL[d.unit_type] || d.unit_type))].join(", ") || "none";
    box.innerHTML = `<div class="stats">
        <div class="stat"><span class="n">${n.buses.length}</span><span class="l">Buses</span></div>
        <div class="stat"><span class="n">${n.lines.length}</span><span class="l">Lines</span></div>
        <div class="stat"><span class="n">${n.transformers.length}</span><span class="l">Transformers</span></div>
        <div class="stat"><span class="n">${n.loads.length}</span><span class="l">Loads</span></div>
        <div class="stat"><span class="n">${n.der_units.length}</span><span class="l">Units (${esc(types)})</span></div>
        <div class="stat"><span class="n">${fmt(n.f_hz, 0)} Hz</span><span class="l">Frequency</span></div>
        <div class="stat"><span class="n">${fmtSmart(n.sn_mva)} MVA</span><span class="l">Base power</span></div>
      </div>${p ? `<p class="muted" style="margin:0.7rem 0 0;font-size:0.85rem;max-width:90ch">${esc(p.description)}</p>` : ""}`;
    updateSidebarNetwork();
  },

  setMode(mode) {
    this.view.setMode(mode);
    $$("#page-network [data-mode]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.mode === mode)));
    $("#net-mode-hint").textContent = mode === "select"
      ? "Click an element to edit it · click a breaker (■) to open/close it · drag a bus to move it · drag the background to pan · scroll to zoom."
      : `Drag from one bus to another to add a ${mode === "wire-line" ? "line" : "transformer"}.`;
  },

  // --- Palette drag & drop ---
  attachChip(chip) {
    chip.addEventListener("pointerdown", evt => {
      evt.preventDefault();
      if (!state.network) { alert("Load a preset or start a new network first."); return; }
      const ghost = el(`<div class="drag-ghost">${esc(chip.textContent)}</div>`);
      document.body.appendChild(ghost);
      const move = e => { ghost.style.left = `${e.clientX}px`; ghost.style.top = `${e.clientY}px`; };
      move(evt);
      const up = e => {
        window.removeEventListener("pointermove", move);
        window.removeEventListener("pointerup", up);
        ghost.remove();
        if (!this.view.containsClient(e.clientX, e.clientY)) return;
        const p = this.view.clientToDiagram(e.clientX, e.clientY);
        this.placeElement(chip.dataset.chip, p.x, p.y);
      };
      window.addEventListener("pointermove", move);
      window.addEventListener("pointerup", up);
    });
  },

  placeElement(kind, x, y) {
    const net = state.network;
    const defaultVn = () => {
      const counts = {};
      net.buses.forEach(b => { counts[b.vn_kv] = (counts[b.vn_kv] || 0) + 1; });
      const top = Object.keys(counts).sort((a, b) => counts[b] - counts[a])[0];
      return top ? parseFloat(top) : 20.0;
    };
    if (kind === "bus") {
      const id = nextId(net.buses);
      net.buses.push({ id, name: `bus${id}`, vn_kv: defaultVn() });
      state.positions[id] = { x, y };
      this.afterEdit({ kind: "bus", id });
      return;
    }
    if (kind === "load") {
      const b = this.view.busAt(x, y, 30);
      if (!b) { alert("Drop a load onto an existing bus."); return; }
      net.loads.push({ bus: b.id, p_mw: 0.0, q_mvar: 0.001 * net.sn_mva, name: "" });
      this.afterEdit({ kind: "bus", id: b.id });
      return;
    }
    // A unit sits on its own terminal bus behind a transformer (network_form.m's
    // convention). Dropped onto/near an existing bus, it's connected to it
    // right away; dropped in free space it stays unconnected until wired.
    const host = this.view.busAt(x, y, 40);
    const busId = nextId(net.buses);
    net.buses.push({ id: busId, name: `${UNIT_LABEL[kind] || kind}${busId}`, vn_kv: host ? host.vn_kv : defaultVn() });
    state.positions[busId] = host ? { x: state.positions[host.id].x + 70, y: state.positions[host.id].y - 50 } : { x, y };
    if (host) net.transformers.push({ hv_bus: host.id, lv_bus: busId, r_pu: 0.0, x_pu: 0.05, sn_mva: net.sn_mva, name: "" });
    net.der_units.push({
      id: nextId(net.der_units), bus: busId, unit_type: kind, bus_type: net.der_units.some(d => d.bus_type === "slack") ? "pq" : "slack",
      v_set_pu: 1.0, angle_set_deg: 0.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0,
      controller: kind === "gfm" ? "droop" : null, xd_pu: kind === "sm" ? 0.2 : null,
    });
    this.afterEdit({ kind: "bus", id: busId });
  },

  commitWire(mode, a, b) {
    const net = state.network;
    if (mode === "wire-line") {
      net.lines.push({ from_bus: a, to_bus: b, r_pu: 0.01, x_pu: 0.1, b_pu: 0.001, length_km: 1.0, name: "" });
      this.afterEdit({ kind: "line", index: net.lines.length - 1 });
    } else {
      // A unit's own bus must be the LV side, whichever way the wire was drawn.
      const aDer = net.der_units.some(d => d.bus === a);
      const [hv, lv] = aDer ? [b, a] : [a, b];
      net.transformers.push({ hv_bus: hv, lv_bus: lv, r_pu: 0.0, x_pu: 0.05, sn_mva: net.sn_mva, name: "" });
      this.afterEdit({ kind: "transformer", index: net.transformers.length - 1 });
    }
  },

  // After a structural edit: bump version, redraw, select the new element.
  afterEdit(select) {
    networkChanged();
    if (select !== undefined) { this.view.selected = select; this.openInspector(select); }
    this.view.render();
    if ($("#net-tables").open) this.renderTables();
  },

  // --- Inspector ---
  closeInspector() {
    this.selected = null;
    this.view.selected = null;
    $("#net-inspector").style.display = "none";
    $("#net-workspace").classList.remove("with-inspector");
  },

  openInspector(sel) {
    if (!sel || !state.network) { this.closeInspector(); this.view.render(); return; }
    this.selected = sel;
    if (sel.kind !== "network") this.view.selected = sel;
    const box = $("#net-inspector");
    box.style.display = "";
    $("#net-workspace").classList.add("with-inspector");
    const net = state.network;
    let html = "";
    const head = (kind, title, deletable = true) => `<div class="insp-head"><div><div class="insp-kind">${kind}</div><h3>${esc(title)}</h3></div><span class="spacer"></span>
      ${deletable ? `<button class="danger small" data-act="delete">Delete</button>` : ""}<button class="ghost small" data-act="close" aria-label="Close">✕</button></div>`;
    if (sel.kind === "network") {
      html = head("Network", "Network settings", false) + `<div class="insp-section" style="border:0;padding-top:0.3rem">${this.formHtml("network", net, "net")}</div>`;
    } else if (sel.kind === "bus") {
      const b = net.buses.find(x => x.id === sel.id);
      if (!b) { this.closeInspector(); return; }
      const der = net.der_units.find(d => d.bus === b.id);
      const loads = net.loads.map((l, i) => [l, i]).filter(([l]) => l.bus === b.id);
      html = head(der ? `Bus · ${UNIT_NAME[der.unit_type] || der.unit_type}` : "Bus", b.name || `bus${b.id}`)
        + `<div class="insp-section" style="border:0;padding-top:0.3rem">${this.formHtml("bus", b, "bus")}</div>`;
      if (der) {
        const di = net.der_units.indexOf(der);
        html += `<div class="insp-section"><h4><span class="swatch" style="display:inline-block;width:10px;height:10px;border-radius:50%;background:${UNIT_COLOR[der.unit_type]}"></span>${esc(UNIT_NAME[der.unit_type] || "Unit")}<span class="spacer"></span><button class="ghost small" data-act="remove-der" data-index="${di}">Remove unit</button></h4>
          ${this.formHtml("der", der, `der:${di}`)}
          ${this.retypeNoteHtml(der)}
          <div id="insp-model-order"></div>
          <div id="insp-ctrl"></div></div>`;
      } else {
        html += `<div class="insp-section"><h4>Generation unit</h4><div class="controls" style="align-items:center"><select id="insp-attach-type">${Object.keys(UNIT_LABEL).map(u => `<option value="${u}">${UNIT_NAME[u]}</option>`).join("")}</select><button class="secondary small" data-act="attach-der">Attach to this bus</button></div>
          <p class="muted" style="font-size:0.74rem;margin:0.4rem 0 0">Units normally sit on their own terminal bus behind a transformer — dropping one from the palette onto a bus does that for you.</p></div>`;
      }
      html += `<div class="insp-section"><h4>Loads at this bus<span class="spacer"></span><button class="secondary small" data-act="add-load">+ Add load</button></h4>
        ${loads.length ? loads.map(([l, i]) => `<div class="sub-card"><div class="controls" style="margin-bottom:0.4rem"><b style="font-size:0.8rem">Load #${i}</b><span style="flex:1"></span><button class="ghost small" data-act="delete-load" data-index="${i}">Remove</button></div>${this.formHtml("load", l, `load:${i}`)}</div>`).join("") : `<p class="empty">No load.</p>`}</div>`;
      const shunts = (net.shunts || []).map((s, i) => [s, i]).filter(([s]) => s.bus === b.id);
      html += `<div class="insp-section"><h4>Shunt compensation<span class="spacer"></span><button class="secondary small" data-act="add-cap">+ Capacitor bank</button><button class="secondary small" data-act="add-reactor">+ Reactor</button></h4>
        ${shunts.length ? shunts.map(([s, i]) => `<div class="sub-card"><div class="controls" style="margin-bottom:0.4rem"><b style="font-size:0.8rem">${esc(shuntKindLabel(s))} #${i}</b><span style="flex:1"></span><button class="ghost small" data-act="delete-shunt" data-index="${i}">Remove</button></div>${this.formHtml("shunt", s, `shunt:${i}`)}</div>`).join("") : `<p class="empty">None. A capacitor bank raises this bus's voltage and adds to its own capacitance; a reactor lowers it and gets a branch of its own.</p>`}</div>`;
      const conns = [...net.lines.map((l, i) => ["line", i, l.from_bus === b.id ? l.to_bus : l.to_bus === b.id ? l.from_bus : null]),
        ...net.transformers.map((t, i) => ["transformer", i, t.hv_bus === b.id ? t.lv_bus : t.lv_bus === b.id ? t.hv_bus : null])].filter(c => c[2] !== null);
      html += `<div class="insp-section"><h4>Connections</h4>${conns.length ? `<div class="controls">${conns.map(([k, i, o]) => `<button class="secondary small" data-act="goto" data-kind="${k}" data-index="${i}">${k === "line" ? "Line" : "Trafo"} → bus ${o}</button>`).join("")}</div>` : `<p class="empty">Not connected — use “Draw line” / “Draw transformer”.</p>`}</div>`;
    } else {
      const arr = sel.kind === "line" ? net.lines : net.transformers;
      const e = arr[sel.index];
      if (!e) { this.closeInspector(); return; }
      const [a, b] = sel.kind === "line" ? [e.from_bus, e.to_bus] : [e.hv_bus, e.lv_bus];
      html = head(sel.kind === "line" ? "Line" : "Transformer", e.name || `${a} → ${b}`)
        + `<div class="insp-section" style="border:0;padding-top:0.3rem">${this.formHtml(sel.kind, e, `${sel.kind}:${sel.index}`)}</div>`;
    }
    box.innerHTML = html;
    this.bindInspector(box);
    this.refreshControlParams();
    this.refreshUnitModelOrder();
    this.view.render();
  },

  formHtml(defsKey, obj, target) {
    const net = state.network;
    const busOpts = cur => net.buses.map(b => `<option value="${b.id}"${b.id === cur ? " selected" : ""}>${b.id} — ${esc(b.name || "")}</option>`).join("");
    return `<div class="form-grid">${FIELD_DEFS[defsKey].map(f => {
      // A field belonging to one kind of unit is left out of the others
      // entirely rather than shown greyed: an empty box with a name on it
      // invites someone to wonder what it would do.
      if (f.only && !f.only.includes(obj.unit_type)) return "";
      const v = obj[f.key];
      const id = `f-${target.replace(":", "-")}-${f.key}`;
      const lab = `<label for="${id}">${esc(f.label)}${f.unit ? ` <span class="unit">(${f.unit})</span>` : ""}</label>`;
      const help = f.help ? `<span class="hint">${esc(f.help)}</span>` : "";
      const cls = `field${f.full ? " full" : ""}`;
      const data = `data-target="${target}" data-field="${f.key}"`;
      if (f.type === "readonly") return `<div class="${cls}">${lab}<input id="${id}" type="text" value="${esc(v)}" disabled></div>`;
      if (f.type === "bool") {
        const on = f.breaker ? v !== false : !!v;
        return `<div class="${cls}${f.breaker ? " brk-field" : ""}"><label class="check" style="font-size:0.84rem"><input id="${id}" type="checkbox" ${data} data-type="bool"${on ? " checked" : ""}> ${esc(f.label)}${f.breaker && !on ? ' <span class="badge crit">open</span>' : ""}</label>${help}</div>`;
      }
      if (f.type === "bus") return `<div class="${cls}">${lab}<select id="${id}" ${data} data-type="int">${busOpts(v)}</select>${help}</div>`;
      if (f.type === "select") {
        return `<div class="${cls}">${lab}<select id="${id}" ${data}>${f.options.map(o => `<option value="${o}"${o === (v ?? "") ? " selected" : ""}>${o === "" ? "(none)" : o}</option>`).join("")}</select>${help}</div>`;
      }
      const type = f.type === "number" ? "number" : "text";
      const si = f.type === "number" ? siSpecs(defsKey, obj, f.key) : [];
      if (si.length) {
        // Per unit (stored) next to its SI equivalent(s); editing one updates the other.
        return `<div class="field full dual">${lab}<div class="dual-row">
          <span class="dual-in"><input id="${id}" type="number" step="any" ${data}${f.nullable ? ' data-nullable="1"' : ""} value="${esc(v ?? "")}"${f.nullable ? ' placeholder="(none)"' : ""}><span class="u">pu</span></span>
          ${si.map((sp, i) => `<span class="dual-in" title="${esc(sp.title)}"><input type="number" step="any" data-si-target="${target}" data-si-field="${f.key}" data-si-i="${i}" value="${v == null ? "" : siFmt(v * sp.factor)}"><span class="u">${esc(sp.unit)}</span></span>`).join("")}
        </div>${help}</div>`;
      }
      return `<div class="${cls}">${lab}<input id="${id}" type="${type}"${type === "number" ? ' step="any"' : ""} ${data}${f.nullable ? ' data-nullable="1"' : ""} value="${esc(v ?? "")}"${f.nullable ? ' placeholder="(none)"' : ""}>${help}</div>`;
    }).join("")}</div>`;
  },

  targetObj(target) {
    const net = state.network;
    const [kind, idx] = target.split(":");
    if (kind === "net") return net;
    if (kind === "bus") return net.buses.find(b => b.id === this.selected.id);
    if (kind === "der") return net.der_units[+idx];
    if (kind === "load") return net.loads[+idx];
    if (kind === "shunt") return (net.shunts || [])[+idx];
    if (kind === "line") return net.lines[+idx];
    if (kind === "transformer") return net.transformers[+idx];
    return null;
  },

  // Applies one form edit to the network; returns true when the inspector
  // needs re-rendering (fields that change what else is editable).
  // Shown once, right after a type change dropped overrides (retypeUnitParams).
  retypeNoteHtml(der) {
    const note = this.retypeNote;
    if (!note || note.id !== der.id || note.unit_type !== der.unit_type) return "";
    this.retypeNote = null;
    return `<div class="notice warn-bg" style="margin-top:0.7rem;font-size:0.78rem">
      <b>${note.dropped.length} parameter override${note.dropped.length === 1 ? "" : "s"} dropped</b> —
      a ${esc(UNIT_NAME[der.unit_type] || der.unit_type)} has no ${esc(note.dropped.slice(0, 6).join(", "))}${note.dropped.length > 6 ? ", …" : ""}.
      This unit now uses its own type's defaults; set them again below if you need to.</div>`;
  },

  applyField(obj, field, input) {
    let v = input.value;
    if (input.dataset.type === "bool") v = input.checked;
    if (field === "closed" && v === false && obj.bus_type === "slack") {
      alert("The slack unit's breaker can't be opened: it is the reference of the power flow and of the dynamic models. Make another unit the slack first.");
      input.checked = true;
      return false;
    }
    else if (input.dataset.type === "int") v = parseInt(v, 10);
    else if (input.type === "number") {
      if (v === "") { if (!input.dataset.nullable) return false; v = null; }
      else { v = parseFloat(v); if (!Number.isFinite(v)) return false; }
    } else if (input.tagName === "SELECT" && v === "") v = null;
    if (field === "length_km" && "r_pu" in obj) {
      if (!(v > 0)) return false;
      const k = v / obj.length_km;
      if (Number.isFinite(k) && k > 0) { obj.r_pu *= k; obj.x_pu *= k; obj.b_pu *= k; }
    }
    obj[field] = v;
    let rerender = input.tagName === "SELECT" || field.endsWith("closed");
    if (field === "unit_type") {
      // The model a unit carries belongs to its type: a converter has a
      // control law, a machine has an exciter, a stabiliser and a governor.
      // Leaving the old type's behind fails validation ("only meaningful
      // for ..."), so they go with the parameters they named.
      if (v === "gfm" && !obj.controller) obj.controller = "droop";
      if (v !== "gfm") obj.controller = null;
      if (v !== "sm") { obj.exciter = null; obj.pss = null; obj.governor = null; }
      if (v !== "sm") obj.xd_pu = null;
      this.retypeUnitParams(obj);
    }
    if (field === "bus_type" && v === "slack") {
      state.network.der_units.forEach(d => { if (d !== obj && d.bus_type === "slack") d.bus_type = "pv"; });
    }
    return rerender;
  },

  // Changing a unit's type makes parameter overrides written for the old one
  // meaningless -- a machine's flux linkages say nothing to a converter -- and
  // the network then fails validation ("overrides parameter(s) it doesn't
  // have"). So drop the ones the new type hasn't got, and keep those it shares
  // (the two converter types share their filter and inner current loop). What
  // was dropped is reported in the inspector rather than thrown away quietly.
  async retypeUnitParams(der) {
    if (!der.params || !Object.keys(der.params).length) return;
    let valid;
    try {
      valid = new Set(Object.keys((await UnitParams.defaultsFor(der)).defaults));
    } catch {
      return;   // server unreachable, or the network is too incomplete to say
    }
    const kept = Object.fromEntries(Object.entries(der.params).filter(([k]) => valid.has(k)));
    const dropped = Object.keys(der.params).filter(k => !(k in kept));
    if (!dropped.length) return;
    der.params = kept;
    this.retypeNote = { id: der.id, unit_type: der.unit_type, dropped };
    networkChanged();
    if (this.selected) this.openInspector(this.selected);
  },

  // Refresh pu inputs (e.g. impedances rescaled by a length change) and the
  // SI values next to them, except the input being typed in.
  syncInputs(box) {
    box.querySelectorAll("input[data-target][data-field]").forEach(inp => {
      if (inp === document.activeElement || inp.type !== "number") return;
      const obj = this.targetObj(inp.dataset.target), v = obj?.[inp.dataset.field];
      if (typeof v === "number" && parseFloat(inp.value) !== v) inp.value = +v.toPrecision(10);
    });
    box.querySelectorAll("input[data-si-target]").forEach(inp => {
      if (inp === document.activeElement) return;
      const obj = this.targetObj(inp.dataset.siTarget);
      const sp = siSpecs(TARGET_KIND[inp.dataset.siTarget.split(":")[0]], obj, inp.dataset.siField)[+inp.dataset.siI];
      const v = obj?.[inp.dataset.siField];
      inp.value = sp && typeof v === "number" ? siFmt(v * sp.factor) : "";
    });
  },

  bindInspector(box) {
    box.querySelectorAll("input[data-si-target]").forEach(inp => inp.addEventListener("input", () => {
      const obj = this.targetObj(inp.dataset.siTarget);
      const sp = siSpecs(TARGET_KIND[inp.dataset.siTarget.split(":")[0]], obj, inp.dataset.siField)[+inp.dataset.siI];
      const v = parseFloat(inp.value);
      if (!sp || !Number.isFinite(v)) return;
      const pu = box.querySelector(`input[data-target="${inp.dataset.siTarget}"][data-field="${inp.dataset.siField}"]`);
      pu.value = +(v / sp.factor).toPrecision(10);
      pu.dispatchEvent(new Event("input"));
    }));
    box.querySelectorAll("[data-target]").forEach(input => {
      const evtName = input.tagName === "SELECT" || input.type === "checkbox" ? "change" : "input";
      input.addEventListener(evtName, () => {
        const obj = this.targetObj(input.dataset.target);
        if (!obj) return;
        const rerender = this.applyField(obj, input.dataset.field, input);
        networkChanged();
        this.view.render();
        if (rerender) this.openInspector(this.selected);
        else if (input.dataset.field === "name") {
          const h = box.querySelector(".insp-head h3");
          if (h && this.selected.kind !== "network" && input.dataset.target.startsWith(this.selected.kind)) h.textContent = input.value || h.textContent;
        }
        if ($("#net-tables").open) this.renderTables();
        this.syncInputs(box);
      });
    });
    box.querySelectorAll("[data-act]").forEach(btn => btn.addEventListener("click", () => this.inspectorAction(btn)));
  },

  inspectorAction(btn) {
    const net = state.network, sel = this.selected, act = btn.dataset.act;
    if (act === "close") { this.closeInspector(); this.view.render(); return; }
    if (act === "goto") { this.openInspector({ kind: btn.dataset.kind, index: +btn.dataset.index }); return; }
    if (act === "add-load") { net.loads.push({ bus: sel.id, p_mw: 0.0, q_mvar: 0.001 * net.sn_mva, name: "" }); this.afterEdit(sel); return; }
    if (act === "delete-load") { net.loads.splice(+btn.dataset.index, 1); this.afterEdit(sel); return; }
    if (act === "add-cap" || act === "add-reactor") {
      const q = 0.05 * net.sn_mva;
      (net.shunts ||= []).push({ bus: sel.id, q_mvar: act === "add-cap" ? -q : q, r_pu: null, closed: true, name: "" });
      this.afterEdit(sel); return;
    }
    if (act === "delete-shunt") { net.shunts.splice(+btn.dataset.index, 1); this.afterEdit(sel); return; }
    if (act === "remove-der") { net.der_units.splice(+btn.dataset.index, 1); this.afterEdit(sel); return; }
    if (act === "attach-der") {
      const kind = $("#insp-attach-type").value;
      net.der_units.push({
        id: nextId(net.der_units), bus: sel.id, unit_type: kind, bus_type: net.der_units.some(d => d.bus_type === "slack") ? "pq" : "slack",
        v_set_pu: 1.0, angle_set_deg: 0.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0,
        controller: kind === "gfm" ? "droop" : null, xd_pu: kind === "sm" ? 0.2 : null,
      });
      this.afterEdit(sel);
      return;
    }
    if (act === "delete") {
      if (sel.kind === "bus") {
        const id = sel.id;
        const nL = net.lines.filter(l => l.from_bus === id || l.to_bus === id).length;
        const nT = net.transformers.filter(t => t.hv_bus === id || t.lv_bus === id).length;
        const nLd = net.loads.filter(l => l.bus === id).length, nD = net.der_units.filter(d => d.bus === id).length;
        const extra = [nL && `${nL} line(s)`, nT && `${nT} transformer(s)`, nLd && `${nLd} load(s)`, nD && `${nD} unit(s)`].filter(Boolean);
        if (extra.length && !confirm(`Delete bus ${id} together with its ${extra.join(", ")}?`)) return;
        net.buses = net.buses.filter(b => b.id !== id);
        net.lines = net.lines.filter(l => l.from_bus !== id && l.to_bus !== id);
        net.transformers = net.transformers.filter(t => t.hv_bus !== id && t.lv_bus !== id);
        net.loads = net.loads.filter(l => l.bus !== id);
        net.der_units = net.der_units.filter(d => d.bus !== id);
        delete state.positions[id];
      } else if (sel.kind === "line") net.lines.splice(sel.index, 1);
      else if (sel.kind === "transformer") net.transformers.splice(sel.index, 1);
      this.closeInspector();
      this.afterEdit(undefined);
      this.closeInspector();
      this.view.render();
    }
  },

    // This unit's own model level, overriding the network default for its
    // type. Rendered in the inspector so the choice sits next to the
    // parameters it decides the fate of.
    refreshUnitModelOrder() {
    const box = $("#insp-model-order");
    if (!box || !this.selected || this.selected.kind !== "bus") return;
    const der = state.network.der_units.find(d => d.bus === this.selected.id);
    if (!der || !ModelOrder.catalogue || !ModelOrder.catalogue[der.unit_type]) { if (box) box.innerHTML = ""; return; }
    const net = ModelOrder.matchingLevel(der.unit_type, ModelOrder.groupModes(der.unit_type));
    const own = der.level || der.states && Object.keys(der.states).length;
    box.innerHTML = `<details class="mo-section"${own ? " open" : ""}>
      <summary>Model order${own ? ` <span class="badge warn">overridden</span>` : ` <span class="muted">— following the network default</span>`}</summary>
      ${ModelOrder.controlsHtml(der.unit_type, der)}
      ${own ? `<button class="ghost small" data-act="mo-reset">Follow the network default (${esc(ModelOrder.levelLabel(der.unit_type, net))})</button>` : ""}
    </details>`;
    ModelOrder.bind(box);
    const reset = $("[data-act='mo-reset']", box);
    if (reset) reset.addEventListener("click", () => { der.level = null; der.states = {}; networkChanged(); this.refreshUnitModelOrder(); });
  },

  refreshControlParams() {
    const box = $("#insp-ctrl");
    if (!box || !this.selected || this.selected.kind !== "bus") return;
    const der = state.network.der_units.find(d => d.bus === this.selected.id);
    if (!der) return;
    if (der.unit_type === "infinite_bus") { box.innerHTML = ""; return; }
    UnitParams.render(box, der);
  },

  // --- Validation ---
  scheduleValidation(now = false) {
    clearTimeout(this._valTimer);
    this._valTimer = setTimeout(() => this.validate(), now ? 0 : 400);
  },
  async validate() {
    const out = $("#net-validation");
    if (!out) return;
    if (!state.network) { out.innerHTML = `<p class="empty">No network loaded.</p>`; return; }
    const version = state.version;
    try {
      const r = await netPost("validate");
      if (version !== state.version) return;
      $("#net-valid-sub").textContent = r.ok ? "" : "— fix the errors before running analyses";
      out.innerHTML = r.issues.length
        ? `<ul class="issues">${r.issues.map(i => `<li><b class="${i.severity === "error" ? "bad" : "warn"}">${i.severity}</b> ${esc(i.message)} <span class="muted">(affects: ${esc(i.affects.join(", "))})</span></li>`).join("")}</ul>`
        : `<p class="status-line"><span class="ok">No issues found</span> — ready for power flow, modal analysis and EMT.</p>`;
    } catch (e) {
      if (version !== state.version) return;
      out.innerHTML = `<ul class="issues"><li><b class="bad">error</b> ${esc(e.message)}</li></ul>`;
    }
  },

  // --- Table view ---
  renderTables() {
    const body = $("#net-tables-body");
    if (!state.network) { body.innerHTML = `<p class="empty">No network loaded.</p>`; return; }
    body.innerHTML = Object.entries(TABLE_KINDS).map(([kind, { label, defs }]) => {
      const rows = state.network[kind] || [];   // a network JSON saved before this kind existed
      const head = defs.map(f => `<th>${esc(f.label)}${f.unit ? ` (${f.unit})` : ""}</th>`).join("") + "<th></th>";
      const bodyRows = rows.map((row, idx) => `<tr>${defs.map(f => {
        const v = row[f.key];
        const data = `data-kind="${kind}" data-index="${idx}" data-field="${f.key}"`;
        if (f.type === "readonly") return `<td class="name">${esc(v)}</td>`;
        if (f.type === "bool") return `<td style="text-align:center"><input type="checkbox" ${data} data-type="bool"${(f.breaker ? v !== false : v) ? " checked" : ""}></td>`;
        if (f.type === "bus") return `<td><select ${data} data-type="int">${state.network.buses.map(b => `<option value="${b.id}"${b.id === v ? " selected" : ""}>${b.id}</option>`).join("")}</select></td>`;
        if (f.type === "select") return `<td><select ${data}>${f.options.map(o => `<option value="${o}"${o === (v ?? "") ? " selected" : ""}>${o || "(none)"}</option>`).join("")}</select></td>`;
        const type = f.type === "number" ? "number" : "text";
        return `<td><input type="${type}"${type === "number" ? ' step="any"' : ""} ${data}${f.nullable ? ' data-nullable="1"' : ""} value="${esc(v ?? "")}"></td>`;
      }).join("")}<td><button class="ghost small" data-del-kind="${kind}" data-del-index="${idx}">Delete</button></td></tr>`).join("");
      return `<div class="kind-block"><h4>${label} <span class="muted" style="font-weight:400">(${rows.length})</span></h4>
        <div class="tablewrap" style="max-height:340px"><table class="edit"><thead><tr>${head}</tr></thead><tbody>${bodyRows}</tbody></table></div>
        <button class="secondary small" style="margin-top:0.4rem" data-add-kind="${kind}">+ Add ${label.toLowerCase().replace(/s$/, "").replace("der unit", "unit")}</button></div>`;
    }).join("");
  },

  bindTableEvents() {
    const body = $("#net-tables-body");
    const onEdit = evt => {
      const t = evt.target;
      if (!t.dataset.kind) return;
      const obj = state.network[t.dataset.kind][+t.dataset.index];
      this.applyField(obj, t.dataset.field, t);
      networkChanged();
      this.view.render();
      if (t.tagName === "SELECT") this.renderTables();
      else {
        // e.g. a length change rescales the same row's R/X/B
        body.querySelectorAll(`input[data-kind="${t.dataset.kind}"][data-index="${t.dataset.index}"]`).forEach(inp => {
          const v = obj[inp.dataset.field];
          if (inp !== t && typeof v === "number" && parseFloat(inp.value) !== v) inp.value = +v.toPrecision(10);
        });
      }
      if (this.selected) this.openInspector(this.selected);
    };
    body.addEventListener("input", evt => { if (evt.target.tagName !== "SELECT") onEdit(evt); });
    body.addEventListener("change", evt => { if (evt.target.tagName === "SELECT") onEdit(evt); });
    body.addEventListener("click", evt => {
      const t = evt.target, net = state.network;
      if (t.dataset.addKind) {
        const kind = t.dataset.addKind, b0 = net.buses[0]?.id ?? 1, b1 = net.buses[1]?.id ?? b0;
        if (kind === "buses") { const id = nextId(net.buses); net.buses.push({ id, name: `bus${id}`, vn_kv: net.buses[0]?.vn_kv ?? 20 }); }
        if (kind === "lines") net.lines.push({ from_bus: b0, to_bus: b1, r_pu: 0.01, x_pu: 0.1, b_pu: 0.001, length_km: 1.0, name: "" });
        if (kind === "transformers") net.transformers.push({ hv_bus: b0, lv_bus: b1, r_pu: 0.0, x_pu: 0.05, sn_mva: net.sn_mva, name: "" });
        if (kind === "loads") net.loads.push({ bus: b0, p_mw: 0.0, q_mvar: 0.001 * net.sn_mva, name: "" });
        if (kind === "der_units") net.der_units.push({ id: nextId(net.der_units), bus: b0, unit_type: "gfl", bus_type: net.der_units.some(d => d.bus_type === "slack") ? "pq" : "slack", v_set_pu: 1.0, angle_set_deg: 0.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0, controller: null, xd_pu: null });
        this.afterEdit(undefined);
      } else if (t.dataset.delKind) {
        net[t.dataset.delKind].splice(+t.dataset.delIndex, 1);
        this.closeInspector();
        this.afterEdit(undefined);
        this.closeInspector();
      }
    });
  },

  // --- Import / export ---
  exportJson() {
    if (!state.network) return;
    const blob = new Blob([JSON.stringify({ network: state.network, positions: state.positions }, null, 2)], { type: "application/json" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${(state.network.name || "network").replace(/[^\w.-]+/g, "_")}.json`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  },
  async importJson(file) {
    if (!file) return;
    try {
      const data = JSON.parse(await file.text());
      const net = data.network || data;
      if (!Array.isArray(net.buses)) throw new Error("not a G2ELin network (no 'buses' list)");
      ["lines", "transformers", "loads", "der_units"].forEach(k => { net[k] ||= []; });
      setNetwork(net, { label: net.name || file.name, positions: data.positions || null });
      $("#net-preset").value = "";
    } catch (e) { alert(`Could not import: ${e.message}`); }
    $("#net-import-file").value = "";
  },
};

function nextId(rows) { return rows.length ? Math.max(...rows.map(r => r.id)) + 1 : 1; }

function presetFamily(id) {
  if (id.startsWith("wscc9_")) return "WSCC 9-bus";
  if (id.startsWith("cigre_islanded")) return "CIGRE MV islanded";
  if (id.startsWith("cigre_interconnected")) return "CIGRE MV interconnected";
  if (id.endsWith("_smib")) return "Single machine – infinite bus";
  if (id.endsWith("_smsm")) return "Single machine – synchronous machine";
  return "Other";
}
