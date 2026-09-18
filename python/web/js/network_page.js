// Network page: choose a preset (plotted immediately) or build from scratch,
// edit on the canvas, and edit every element parameter in the inspector that
// opens on the right when an element is clicked.

const FIELD_DEFS = {
  network: [
    { key: "name", label: "Network name", type: "text", full: true },
    { key: "f_hz", label: "Frequency", unit: "Hz", type: "number" },
    { key: "sn_mva", label: "Base power", unit: "MVA", type: "number", help: "All per-unit values are on this base." },
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
    { key: "length_km", label: "Length", unit: "km", type: "number" },
    { key: "name", label: "Name", type: "text", full: true },
  ],
  transformer: [
    { key: "hv_bus", label: "HV bus", type: "bus" },
    { key: "lv_bus", label: "LV bus", type: "bus", help: "A unit's own terminal bus must be the LV side." },
    { key: "r_pu", label: "Resistance R", unit: "pu", type: "number" },
    { key: "x_pu", label: "Reactance X", unit: "pu", type: "number" },
    { key: "sn_mva", label: "Rating", unit: "MVA", type: "number" },
    { key: "name", label: "Name", type: "text" },
  ],
  load: [
    { key: "p_mw", label: "Active power P", unit: "MW", type: "number" },
    { key: "q_mvar", label: "Reactive power Q", unit: "MVAr", type: "number" },
    { key: "name", label: "Name", type: "text", full: true },
  ],
  der: [
    { key: "id", label: "Unit id", type: "readonly" },
    { key: "unit_type", label: "Unit type", type: "select", options: ["sm", "gfm", "gfl", "infinite_bus"] },
    { key: "bus_type", label: "Load-flow bus type", type: "select", options: ["slack", "pv", "pq"], help: "Exactly one unit must be the slack." },
    { key: "v_set_pu", label: "Voltage setpoint", unit: "pu", type: "number" },
    { key: "p_set_mw", label: "Active power setpoint", unit: "MW", type: "number" },
    { key: "q_set_mvar", label: "Reactive power setpoint", unit: "MVAr", type: "number" },
    { key: "p_cons_mw", label: "Auxiliary load P", unit: "MW", type: "number" },
    { key: "q_cons_mvar", label: "Auxiliary load Q", unit: "MVAr", type: "number" },
    { key: "controller", label: "GFM outer control", type: "select", options: ["", "droop", "droop_filtered", "dvoc", "vsm", "matching"], gfmOnly: true },
    { key: "xd_pu", label: "Transient reactance Xd", unit: "pu", type: "number", nullable: true, help: "Used for SCR calculations only." },
  ],
};

// Table view: same schema, one table per element kind.
const TABLE_KINDS = {
  buses: { label: "Buses", defs: FIELD_DEFS.bus },
  lines: { label: "Lines", defs: FIELD_DEFS.line },
  transformers: { label: "Transformers", defs: FIELD_DEFS.transformer },
  loads: { label: "Loads", defs: [{ key: "bus", label: "Bus", type: "bus" }, ...FIELD_DEFS.load] },
  der_units: { label: "DER units", defs: [FIELD_DEFS.der[0], { key: "bus", label: "Bus", type: "bus" }, ...FIELD_DEFS.der.slice(1)] },
};

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
      <div class="workspace" id="net-workspace">
        <div class="card">
          <div class="canvas-tools">
            <div class="palette"><span class="palette-label">Drag to add</span>
              <span class="chip" data-chip="bus"><span class="dot" style="background:#77766f"></span>Bus</span>
              <span class="chip" data-chip="load"><span class="dot" style="background:var(--text-secondary);border-radius:1px"></span>Load</span>
              ${Object.keys(UNIT_LABEL).map(u => `<span class="chip" data-chip="${u}" title="${UNIT_NAME[u]}"><span class="dot" style="background:${UNIT_COLOR[u]}"></span>${UNIT_LABEL[u]}</span>`).join("")}
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
      <div class="card"><div class="card-title">Validation <span class="card-sub" id="net-valid-sub"></span></div><div id="net-validation"><p class="empty">No network loaded.</p></div></div>
      <details class="card" id="net-tables"><summary>Table view — bulk edit every element</summary><div id="net-tables-body"></div></details>`;

    this.view = new NetworkView($("#net-canvas"), {
      editable: true,
      tooltip: elementTooltip,
      onSelect: sel => this.openInspector(sel),
      onWire: (mode, a, b) => this.commitWire(mode, a, b),
    });
    $$("#page-network [data-mode]").forEach(btn => btn.addEventListener("click", () => this.setMode(btn.dataset.mode)));
    this.setMode("select");
    $$("#page-network .chip").forEach(chip => this.attachChip(chip));
    $("#net-labels").addEventListener("change", e => { this.view.opts.showLabels = e.target.checked; this.view.render(); });
    $("#net-autolayout").addEventListener("click", async () => { try { await autoLayout(); this.view.fit(); } catch (e) { alert(e.message); } });
    $("#net-settings").addEventListener("click", () => this.openInspector({ kind: "network" }));
    $("#net-new").addEventListener("click", () => {
      if (state.network && !confirm("Start a new empty network? The current one is not saved.")) return;
      setNetwork({ name: "custom_network", f_hz: 50.0, sn_mva: 100.0, buses: [], lines: [], transformers: [], loads: [], der_units: [] }, { label: "New network" });
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

    on("network:loaded", () => { this.closeInspector(); this.refreshAll(); });
    on("network:layout", () => { this.view.render(); this.view.fit(); });
    on("network:changed", () => { this.renderSummary(); this.scheduleValidation(); });
    on("network:derived", () => this.refreshControlParams());
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
        <div class="stat"><span class="n">${fmt(n.sn_mva, 0)} MVA</span><span class="l">Base power</span></div>
      </div>${p ? `<p class="muted" style="margin:0.7rem 0 0;font-size:0.85rem;max-width:90ch">${esc(p.description)}</p>` : ""}`;
    updateSidebarNetwork();
  },

  setMode(mode) {
    this.view.setMode(mode);
    $$("#page-network [data-mode]").forEach(b => b.setAttribute("aria-pressed", String(b.dataset.mode === mode)));
    $("#net-mode-hint").textContent = mode === "select"
      ? "Click an element to edit it · drag a bus to move it · drag the background to pan · scroll to zoom."
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
      v_set_pu: 1.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0,
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
          <div id="insp-ctrl"></div></div>`;
      } else {
        html += `<div class="insp-section"><h4>Generation unit</h4><div class="controls" style="align-items:center"><select id="insp-attach-type">${Object.keys(UNIT_LABEL).map(u => `<option value="${u}">${UNIT_NAME[u]}</option>`).join("")}</select><button class="secondary small" data-act="attach-der">Attach to this bus</button></div>
          <p class="muted" style="font-size:0.74rem;margin:0.4rem 0 0">Units normally sit on their own terminal bus behind a transformer — dropping one from the palette onto a bus does that for you.</p></div>`;
      }
      html += `<div class="insp-section"><h4>Loads at this bus<span class="spacer"></span><button class="secondary small" data-act="add-load">+ Add load</button></h4>
        ${loads.length ? loads.map(([l, i]) => `<div class="sub-card"><div class="controls" style="margin-bottom:0.4rem"><b style="font-size:0.8rem">Load #${i}</b><span style="flex:1"></span><button class="ghost small" data-act="delete-load" data-index="${i}">Remove</button></div>${this.formHtml("load", l, `load:${i}`)}</div>`).join("") : `<p class="empty">No load.</p>`}</div>`;
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
    this.view.render();
  },

  formHtml(defsKey, obj, target) {
    const net = state.network;
    const busOpts = cur => net.buses.map(b => `<option value="${b.id}"${b.id === cur ? " selected" : ""}>${b.id} — ${esc(b.name || "")}</option>`).join("");
    return `<div class="form-grid">${FIELD_DEFS[defsKey].map(f => {
      const v = obj[f.key];
      const id = `f-${target.replace(":", "-")}-${f.key}`;
      const lab = `<label for="${id}">${esc(f.label)}${f.unit ? ` <span class="unit">(${f.unit})</span>` : ""}</label>`;
      const help = f.help ? `<span class="hint">${esc(f.help)}</span>` : "";
      const cls = `field${f.full ? " full" : ""}`;
      const data = `data-target="${target}" data-field="${f.key}"`;
      if (f.type === "readonly") return `<div class="${cls}">${lab}<input id="${id}" type="text" value="${esc(v)}" disabled></div>`;
      if (f.type === "bus") return `<div class="${cls}">${lab}<select id="${id}" ${data} data-type="int">${busOpts(v)}</select>${help}</div>`;
      if (f.type === "select") {
        const disabled = f.gfmOnly && obj.unit_type !== "gfm" ? " disabled" : "";
        return `<div class="${cls}">${lab}<select id="${id}" ${data}${disabled}>${f.options.map(o => `<option value="${o}"${o === (v ?? "") ? " selected" : ""}>${o === "" ? "(none)" : o}</option>`).join("")}</select>${help}</div>`;
      }
      const type = f.type === "number" ? "number" : "text";
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
    if (kind === "line") return net.lines[+idx];
    if (kind === "transformer") return net.transformers[+idx];
    return null;
  },

  // Applies one form edit to the network; returns true when the inspector
  // needs re-rendering (fields that change what else is editable).
  applyField(obj, field, input) {
    let v = input.value;
    if (input.dataset.type === "int") v = parseInt(v, 10);
    else if (input.type === "number") {
      if (v === "") { if (!input.dataset.nullable) return false; v = null; }
      else { v = parseFloat(v); if (!Number.isFinite(v)) return false; }
    } else if (input.tagName === "SELECT" && v === "") v = null;
    obj[field] = v;
    let rerender = input.tagName === "SELECT";
    if (field === "unit_type") {
      if (v === "gfm" && !obj.controller) obj.controller = "droop";
      if (v !== "gfm") obj.controller = null;
    }
    if (field === "bus_type" && v === "slack") {
      state.network.der_units.forEach(d => { if (d !== obj && d.bus_type === "slack") d.bus_type = "pv"; });
    }
    return rerender;
  },

  bindInspector(box) {
    box.querySelectorAll("[data-target]").forEach(input => {
      const evtName = input.tagName === "SELECT" ? "change" : "input";
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
    if (act === "remove-der") { net.der_units.splice(+btn.dataset.index, 1); this.afterEdit(sel); return; }
    if (act === "attach-der") {
      const kind = $("#insp-attach-type").value;
      net.der_units.push({
        id: nextId(net.der_units), bus: sel.id, unit_type: kind, bus_type: net.der_units.some(d => d.bus_type === "slack") ? "pq" : "slack",
        v_set_pu: 1.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0,
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

  // Derived (read-only) control parameters for the selected unit -- updated
  // in place when /topology comes back, so typing isn't interrupted.
  refreshControlParams() {
    const box = $("#insp-ctrl");
    if (!box || !this.selected || this.selected.kind !== "bus") return;
    const info = state.derInfo[this.selected.id];
    const cp = info && info.control_params;
    box.innerHTML = cp
      ? `<details style="margin-top:0.8rem"><summary style="cursor:pointer;font-size:0.8rem;font-weight:600">Control & electrical parameters <span class="muted" style="font-weight:400">(derived, read-only)</span></summary>
          <div class="tablewrap" style="margin-top:0.5rem;max-height:300px"><table class="ro-table"><tbody>${Object.entries(cp).map(([k, v]) => `<tr><td class="name">${esc(k)}</td><td>${fmtSmart(v)}</td></tr>`).join("")}</tbody></table></div></details>`
      : "";
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
      const rows = state.network[kind];
      const head = defs.map(f => `<th>${esc(f.label)}${f.unit ? ` (${f.unit})` : ""}</th>`).join("") + "<th></th>";
      const bodyRows = rows.map((row, idx) => `<tr>${defs.map(f => {
        const v = row[f.key];
        const data = `data-kind="${kind}" data-index="${idx}" data-field="${f.key}"`;
        if (f.type === "readonly") return `<td class="name">${esc(v)}</td>`;
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
        if (kind === "der_units") net.der_units.push({ id: nextId(net.der_units), bus: b0, unit_type: "gfl", bus_type: net.der_units.some(d => d.bus_type === "slack") ? "pq" : "slack", v_set_pu: 1.0, p_set_mw: 0.0, q_set_mvar: 0.0, p_cons_mw: 0.0, q_cons_mvar: 0.0, controller: null, xd_pu: null });
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
  if (id.endsWith("_smib")) return "Single machine – infinite bus";
  return "Other";
}
