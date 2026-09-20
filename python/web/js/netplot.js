// One-line network diagram, shared by the Network page (editable: drag buses,
// wire lines/transformers, drop elements from a palette) and the analysis
// pages (read-only, coloured by results). Bus positions live in
// state.positions, so every page shows the same layout.
//
// A selection is {kind: "bus", id} | {kind: "line", index} | {kind: "transformer", index}.
//
// Breakers (small squares: filled = closed, hollow red = open) sit at both ends
// of every line/transformer, on every load's connection and on every unit;
// with opts.onBreaker they toggle on click. Whatever open breakers take out of
// service is drawn dashed / faded (serviceState in core.js).

class NetworkView {
  constructor(host, opts = {}) {
    this.host = host;
    this.opts = opts;               // {editable, onSelect, onBreaker(br), tooltip(sel), style: {bus(bus, der), edge(kind, index)}, showLabels}
    this.mode = "select";           // editable only: "select" | "wire-line" | "wire-transformer"
    this.selected = null;
    host.classList.add("net-host");
    this.svg = svgEl("svg", { viewBox: `0 0 ${DIAGRAM_W} ${DIAGRAM_H}`, role: "img", "aria-label": "Network diagram" });
    this.layer = svgEl("g");
    this.svg.appendChild(this.layer);
    host.appendChild(this.svg);
    this.emptyNote = el(`<div class="net-empty" style="display:none"></div>`);
    host.appendChild(this.emptyNote);
    const zoomBox = el(`<div class="net-zoom"><button type="button" title="Zoom in">+</button><button type="button" title="Zoom out">−</button><button type="button" title="Reset view" style="font-size:0.75rem">⟲</button></div>`);
    host.appendChild(zoomBox);
    this.zoom = attachSvgZoomPan(this.svg, { canStartPan: evt => !evt.target.closest(".nd-bus") || !this.opts.editable });
    const [zin, zout, zreset] = zoomBox.querySelectorAll("button");
    zin.addEventListener("click", () => this.zoom.zoomIn());
    zout.addEventListener("click", () => this.zoom.zoomOut());
    zreset.addEventListener("click", () => this.fit());
    this.svg.addEventListener("click", evt => {
      if (this.svg.__justPanned) return;
      if (evt.target === this.svg || evt.target === this.layer || evt.target.classList.contains("nd-bg")) this.select(null, true);
    });
    this.wire = null;
    if (opts.editable) this.attachWireHandlers();
  }

  setStyle(style) { this.opts.style = style; this.render(); }
  setMode(mode) { this.mode = mode; this.cancelWire(); }

  select(sel, fromUser = false) {
    this.selected = sel;
    this.render();
    if (fromUser && this.opts.onSelect) this.opts.onSelect(sel);
  }

  // Fit the viewBox to the buses' bounding box.
  fit() {
    const pts = Object.values(state.positions);
    if (!pts.length) { this.zoom.reset(); return; }
    const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    const pad = 60;
    let x0 = Math.min(...xs) - pad, x1 = Math.max(...xs) + pad, y0 = Math.min(...ys) - pad, y1 = Math.max(...ys) + pad;
    const ar = DIAGRAM_W / DIAGRAM_H;
    let w = x1 - x0, h = y1 - y0;
    if (w / h > ar) { const nh = w / ar; y0 -= (nh - h) / 2; h = nh; } else { const nw = h * ar; x0 -= (nw - w) / 2; w = nw; }
    this.zoom.set({ x: x0, y: y0, width: w, height: h });
  }

  clientToDiagram(cx, cy) {
    const r = this.svg.getBoundingClientRect(), vb = this.svg.viewBox.baseVal;
    return { x: vb.x + (cx - r.left) * vb.width / r.width, y: vb.y + (cy - r.top) * vb.height / r.height };
  }
  containsClient(cx, cy) {
    const r = this.svg.getBoundingClientRect();
    return cx >= r.left && cx <= r.right && cy >= r.top && cy <= r.bottom;
  }
  busAt(x, y, threshold = 22) {
    let best = null, bestD = threshold;
    for (const b of state.network?.buses || []) {
      const p = state.positions[b.id];
      if (!p) continue;
      const d = Math.hypot(p.x - x, p.y - y);
      if (d < bestD) { bestD = d; best = b; }
    }
    return best;
  }

  render() {
    const net = state.network;
    this.layer.innerHTML = "";
    this.cancelWire();
    if (!net || !net.buses.length) {
      this.emptyNote.style.display = "";
      this.emptyNote.textContent = !net ? "No network loaded." : (this.opts.editable ? "Empty network — drag a Bus or a unit from the palette onto the canvas." : "This network has no buses yet.");
      return;
    }
    this.emptyNote.style.display = "none";
    ensurePositions();
    const style = this.opts.style || {};
    const pos = id => state.positions[id] || { x: DIAGRAM_W / 2, y: DIAGRAM_H / 2 };
    const updaters = {};
    const addUpd = (id, fn) => { (updaters[id] ||= []).push(fn); };
    this.layer.appendChild(svgEl("rect", { class: "nd-bg", x: -5000, y: -5000, width: 10000, height: 10000, fill: "transparent" }));

    const derByBus = {};
    net.der_units.forEach(d => { derByBus[d.bus] = d; });
    const loadsByBus = {};
    net.loads.forEach((l, i) => { (loadsByBus[l.bus] ||= []).push([l, i]); });
    const service = serviceState(net);
    const busR = id => (derByBus[id] ? 16 : 8);
    // Breaker offset from a bus centre along an edge (further out at a unit,
    // whose own breaker badge sits right against its circle).
    const brkGap = id => busR(id) + (derByBus[id] ? 19 : 9);
    const sel = this.selected;
    const isSel = (kind, key) => sel && sel.kind === kind && (kind === "bus" ? sel.id === key : sel.index === key);

    const drawEdge = (kind, index, a, b) => {
      const st = (style.edge && style.edge(kind, index)) || {};
      const e = kind === "line" ? net.lines[index] : net.transformers[index];
      const inService = kind === "line" ? service.lines[index] : service.transformers[index];
      const g = svgEl("g", { class: inService ? "" : "nd-out" });
      const main = svgEl("line", { class: "nd-edge", stroke: st.stroke || (kind === "transformer" ? "#8f8e88" : "#9d9c95"), "stroke-width": st.width || 2.2 });
      const hit = svgEl("line", { class: "nd-hit" });
      const hov = svgEl("line", { class: `nd-hover${isSel(kind, index) ? " nd-sel" : ""}` });
      g.append(main, hit, hov);
      let c1 = null, c2 = null;
      if (kind === "transformer") {
        c1 = svgEl("circle", { class: "nd-trafo", r: 6 });
        c2 = svgEl("circle", { class: "nd-trafo", r: 6 });
        g.append(c1, c2);
      }
      const [endA, endB] = kind === "line" ? ["from", "to"] : ["hv", "lv"];
      const bkA = this.breakerNode({ kind, index, end: endA }, e[`${endA}_closed`] !== false);
      const bkB = this.breakerNode({ kind, index, end: endB }, e[`${endB}_closed`] !== false);
      g.append(bkA, bkB);
      const upd = () => {
        const p = pos(a), q = pos(b);
        [main, hit, hov].forEach(l => { l.setAttribute("x1", p.x); l.setAttribute("y1", p.y); l.setAttribute("x2", q.x); l.setAttribute("y2", q.y); });
        const L = Math.hypot(q.x - p.x, q.y - p.y) || 1;
        const ux = (q.x - p.x) / L, uy = (q.y - p.y) / L;
        if (c1) {
          const mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2;
          c1.setAttribute("cx", mx - ux * 4); c1.setAttribute("cy", my - uy * 4);
          c2.setAttribute("cx", mx + ux * 4); c2.setAttribute("cy", my + uy * 4);
        }
        // Breakers near each end, unless the edge is too short to fit them.
        const ga = Math.min(brkGap(a), L * 0.3), gb = Math.min(brkGap(b), L * 0.3);
        bkA.setAttribute("transform", `translate(${p.x + ux * ga},${p.y + uy * ga})`);
        bkB.setAttribute("transform", `translate(${q.x - ux * gb},${q.y - uy * gb})`);
      };
      upd();
      addUpd(a, upd); addUpd(b, upd);
      const selObj = { kind, index };
      hit.addEventListener("click", evt => { evt.stopPropagation(); this.select(selObj, true); });
      this.attachTooltip(hit, selObj);
      this.layer.appendChild(g);
    };
    net.lines.forEach((l, i) => drawEdge("line", i, l.from_bus, l.to_bus));
    net.transformers.forEach((t, i) => drawEdge("transformer", i, t.hv_bus, t.lv_bus));

    net.buses.forEach(b => {
      const der = derByBus[b.id];
      const st = (style.bus && style.bus(b, der)) || {};
      const r = der ? 16 : 8;
      const dead = !service.energized.has(b.id) || (der && !service.units[der.id]);
      const g = svgEl("g", { class: `nd-bus${isSel("bus", b.id) ? " sel" : ""}${dead ? " nd-out" : ""}`, "data-bus-id": b.id });
      g.appendChild(svgEl("circle", { class: "ring", r: r + 5, fill: "none", stroke: "none" }));
      const fill = st.fill || (der ? UNIT_COLOR[der.unit_type] || "#666" : "#77766f");
      const core = svgEl("circle", { class: "core", r, fill, "data-bus-id": b.id });
      if (st.fill && der) { core.setAttribute("stroke", UNIT_COLOR[der.unit_type] || "#666"); core.setAttribute("stroke-width", 3.5); }
      g.appendChild(core);
      if (der) {
        const t = svgEl("text", { class: "nd-unit-label", fill: st.fill ? (isDark(st.fill) ? "#fff" : "#111") : "#fff" });
        t.textContent = UNIT_LABEL[der.unit_type] || "?";
        g.appendChild(t);
      }
      // Loads: one stub each (fanned out when a bus has several), with its breaker.
      const loads = loadsByBus[b.id] || [];
      loads.forEach(([ld, li], k) => {
        const dx = (k - (loads.length - 1) / 2) * 15;
        const lg = svgEl("g", { class: service.loads[li] ? "" : "nd-out" });
        lg.appendChild(svgEl("line", { x1: dx * 0.4, y1: r, x2: dx, y2: r + 17, stroke: "var(--text-secondary)", "stroke-width": 1.5 }));
        lg.appendChild(svgEl("polygon", { class: "nd-load", points: `${dx - 5.5},${r + 16} ${dx + 5.5},${r + 16} ${dx},${r + 25}` }));
        const bk = this.breakerNode({ kind: "load", index: li }, ld.closed !== false);
        bk.setAttribute("transform", `translate(${dx * 0.7},${r + 9})`);
        lg.appendChild(bk);
        g.appendChild(lg);
      });
      // A unit's own breaker: a badge on its circle, towards its transformer.
      if (der) {
        const tr = net.transformers.find(t => t.lv_bus === b.id);
        const bk = this.breakerNode({ kind: "unit", id: der.id }, der.closed !== false);
        const place = () => {
          const p = pos(b.id), q = tr ? pos(tr.hv_bus) : { x: p.x - 1, y: p.y };
          const L = Math.hypot(q.x - p.x, q.y - p.y) || 1;
          bk.setAttribute("transform", `translate(${(q.x - p.x) / L * (r + 7)},${(q.y - p.y) / L * (r + 7)})`);
        };
        place();
        addUpd(b.id, place);
        if (tr) addUpd(tr.hv_bus, place);
        g.appendChild(bk);
      }
      if (this.opts.showLabels !== false) {
        const lab = svgEl("text", { class: "nd-bus-label", x: r + 4, y: -r + 1 });
        lab.textContent = b.name || `bus${b.id}`;
        g.appendChild(lab);
      }
      const upd = () => { const p = pos(b.id); g.setAttribute("transform", `translate(${p.x},${p.y})`); };
      upd();
      this.attachTooltip(core, { kind: "bus", id: b.id });
      this.attachBusInteractions(g, b, upd, updaters);
      this.layer.appendChild(g);
    });
  }

  // One breaker symbol (positioned by the caller).
  breakerNode(br, closed) {
    const g = svgEl("g", {
      class: `nd-brk${closed ? "" : " open"}${this.opts.onBreaker ? " live" : ""}`,
      "data-brk": br.kind === "unit" ? `unit:${br.id}` : `${br.kind}:${br.index}${br.end ? `:${br.end}` : ""}`,
    });
    g.appendChild(svgEl("rect", { class: "hit", x: -8, y: -8, width: 16, height: 16 }));
    g.appendChild(svgEl("rect", { class: "sq", x: -4, y: -4, width: 8, height: 8, rx: 1 }));
    if (!closed) g.appendChild(svgEl("line", { class: "slash", x1: -6, y1: 6, x2: 6, y2: -6 }));
    const tip = () => {
      const label = breakerLabel(state.network, br);
      return `<div class="tt-title">Breaker — ${closed ? "closed" : '<span style="color:var(--critical)">open</span>'}</div>${ttTable([["Element", label]])}${this.opts.onBreaker ? `<div class="muted" style="font-size:0.72rem;margin-top:0.3rem">Click to ${closed ? "open" : "close"} it.</div>` : ""}`;
    };
    g.addEventListener("mouseenter", evt => { if (!this.dragging) showTooltip(tip(), evt); });
    g.addEventListener("mousemove", evt => moveTooltip(evt));
    g.addEventListener("mouseleave", hideTooltip);
    // Own the press, so it neither drags the bus nor pans the view.
    g.addEventListener("pointerdown", evt => evt.stopPropagation());
    g.addEventListener("click", evt => {
      evt.stopPropagation();
      if (!this.opts.onBreaker) return;
      hideTooltip();
      this.opts.onBreaker(br);
    });
    return g;
  }

  attachTooltip(node, sel) {
    if (!this.opts.tooltip) return;
    node.addEventListener("mouseenter", evt => { if (!this.dragging) { const h = this.opts.tooltip(sel); if (h) showTooltip(h, evt); } });
    node.addEventListener("mousemove", evt => { if (!this.dragging) moveTooltip(evt); });
    node.addEventListener("mouseleave", hideTooltip);
  }

  attachBusInteractions(g, bus, upd, updaters) {
    let moved = false;
    g.addEventListener("pointerdown", evt => {
      if (!this.opts.editable || evt.button !== 0) return;
      evt.stopPropagation();
      if (this.mode !== "select") { this.startWire(bus.id, evt); return; }
      g.setPointerCapture(evt.pointerId);
      const sx = evt.clientX, sy = evt.clientY;
      moved = false;
      const onMove = e => {
        if (!moved && Math.hypot(e.clientX - sx, e.clientY - sy) < 3) return;
        moved = true; this.dragging = true; hideTooltip();
        state.positions[bus.id] = this.clientToDiagram(e.clientX, e.clientY);
        upd();
        (updaters[bus.id] || []).forEach(fn => fn());
      };
      const onUp = () => {
        g.removeEventListener("pointermove", onMove);
        this.dragging = false;
        if (moved) emit("network:layout-moved");
      };
      g.addEventListener("pointermove", onMove);
      g.addEventListener("pointerup", onUp, { once: true });
      g.addEventListener("lostpointercapture", onUp, { once: true });
    });
    g.addEventListener("click", evt => {
      evt.stopPropagation();
      if (moved) { moved = false; return; }
      // In the drawing modes the svg's pointerup handles a click (it holds
      // the pointer capture there); here only select mode needs handling.
      if (this.mode === "select") this.select({ kind: "bus", id: bus.id }, true);
    });
  }

  // --- Wiring (editable only) ---
  startWire(fromId, evt) {
    const p = state.positions[fromId];
    const line = svgEl("line", { class: "nd-wire-preview", x1: p.x, y1: p.y, x2: p.x, y2: p.y });
    this.layer.appendChild(line);
    this.wire = { fromId, line, x: evt.clientX, y: evt.clientY };
    this.svg.setPointerCapture(evt.pointerId);
  }
  cancelWire() {
    if (this.wire) { this.wire.line.remove(); this.wire = null; }
  }
  attachWireHandlers() {
    this.svg.addEventListener("pointermove", evt => {
      if (!this.wire) return;
      const p = this.clientToDiagram(evt.clientX, evt.clientY);
      this.wire.line.setAttribute("x2", p.x); this.wire.line.setAttribute("y2", p.y);
    });
    this.svg.addEventListener("pointerup", evt => {
      if (!this.wire) return;
      const { fromId: from, x, y } = this.wire;
      this.cancelWire();
      // A press-and-release without dragging is a click: select the bus
      // (so the inspector opens in the drawing modes too), don't draw.
      if (Math.hypot(evt.clientX - x, evt.clientY - y) < 5) { this.select({ kind: "bus", id: from }, true); return; }
      const p = this.clientToDiagram(evt.clientX, evt.clientY);
      const target = this.busAt(p.x, p.y, 24);
      if (target && target.id !== from && this.opts.onWire) this.opts.onWire(this.mode, from, target.id);
    });
  }
}

// Default hover content: an element's own parameters (no results).
function elementTooltip(sel) {
  const net = state.network;
  if (!net || !sel) return null;
  if (sel.kind === "bus") {
    const b = net.buses.find(x => x.id === sel.id);
    if (!b) return null;
    const sv = serviceState(net);
    const rows = [["Bus id", b.id], ["Nominal voltage", `${fmt(b.vn_kv, 2)} kV`]];
    if (!sv.energized.has(b.id)) rows.push(["State", "de-energized (open breakers)"]);
    const der = net.der_units.find(d => d.bus === b.id);
    if (der) {
      rows.push({ sep: UNIT_NAME[der.unit_type] || der.unit_type }, ["Bus type", der.bus_type]);
      // Which unit each island is solved against can differ from the
      // designated slack once breakers are open (see serviceState).
      if (sv.references.includes(der.id)) rows.push(["Role", sv.references.length > 1 ? "reference of its island" : "reference (slack)"]);
      rows.push(["P set", `${fmt(der.p_set_mw, 3)} MW`], ["Q set", `${fmt(der.q_set_mvar, 3)} MVAr`], ["V set", `${fmt(der.v_set_pu, 3)} pu`], ["Breaker", der.closed === false ? "open" : "closed"]);
    }
    const loads = net.loads.map((l, i) => [l, i]).filter(([l]) => l.bus === b.id);
    if (loads.length) {
      rows.push({ sep: `Load${loads.length > 1 ? "s" : ""}` });
      loads.forEach(([l, i]) => rows.push([l.name || `load #${i}`, `${fmt(l.p_mw, 3)} MW, ${fmt(l.q_mvar, 3)} MVAr${l.closed === false ? " (open)" : ""}`]));
    }
    return `<div class="tt-title">${esc(b.name || `bus${b.id}`)}</div>` + ttTable(rows);
  }
  const e = sel.kind === "line" ? net.lines[sel.index] : net.transformers[sel.index];
  if (!e) return null;
  const [a, b] = sel.kind === "line" ? [e.from_bus, e.to_bus] : [e.hv_bus, e.lv_bus];
  const rows = [["Buses", `${a} → ${b}`], ["R", `${fmtSmart(e.r_pu)} pu`], ["X", `${fmtSmart(e.x_pu)} pu`]];
  const [ea, eb] = sel.kind === "line" ? ["from", "to"] : ["hv", "lv"];
  if (e[`${ea}_closed`] === false || e[`${eb}_closed`] === false) rows.push(["Breakers", `${e[`${ea}_closed`] === false ? "open" : "closed"} / ${e[`${eb}_closed`] === false ? "open" : "closed"} — out of service`]);
  if (sel.kind === "line") rows.push(["B", `${fmtSmart(e.b_pu)} pu`], ["Length", `${fmt(e.length_km, 3)} km`]);
  else rows.push(["Rating", `${fmt(e.sn_mva, 2)} MVA`]);
  return `<div class="tt-title">${sel.kind === "line" ? "Line" : "Transformer"} ${esc(e.name || `#${sel.index}`)}</div>` + ttTable(rows);
}

function unitLegendHtml() {
  return `<div class="legend">${Object.keys(UNIT_LABEL).map(u => `<span><span class="swatch" style="background:${UNIT_COLOR[u]}"></span>${UNIT_LABEL[u]} — ${UNIT_NAME[u]}</span>`).join("")}
    <span><span class="swatch" style="background:#77766f"></span>Network bus</span>
    <span><svg width="12" height="12" style="vertical-align:-2px;margin-right:0.3em"><polygon points="1,2 11,2 6,11" fill="var(--text-secondary)"/></svg>Load</span>
    <span><svg width="20" height="12" style="vertical-align:-2px;margin-right:0.3em"><circle cx="7" cy="6" r="5" fill="#fff" stroke="#52514e"/><circle cx="13" cy="6" r="5" fill="#fff" stroke="#52514e" fill-opacity="0.6"/></svg>Transformer</span>
    <span><svg width="12" height="12" style="vertical-align:-2px;margin-right:0.3em"><rect x="2" y="2" width="8" height="8" rx="1" fill="var(--text-primary)"/></svg>Breaker closed</span>
    <span><svg width="14" height="14" style="vertical-align:-3px;margin-right:0.3em"><rect x="3" y="3" width="8" height="8" rx="1" fill="var(--surface-1)" stroke="var(--critical)" stroke-width="1.6"/><line x1="1" y1="13" x2="13" y2="1" stroke="var(--critical)" stroke-width="1.4"/></svg>Breaker open</span>
    <span><svg width="22" height="12" style="vertical-align:-2px;margin-right:0.3em"><line x1="1" y1="6" x2="21" y2="6" stroke="#9d9c95" stroke-width="2" stroke-dasharray="4 3" opacity="0.6"/></svg>Out of service</span></div>`;
}
