// Shared state, API access, events and small UI helpers used by every page.
//
// The whole UI works on one client-held `Network` JSON (state.network): a
// preset is loaded by fetching its Network definition once, after which every
// analysis goes through /api/network/* with that JSON. Editing is therefore
// just mutating state.network and calling networkChanged() -- nothing is ever
// persisted server-side.

const state = {
  presets: [],
  presetId: null,         // preset the current network came from (null for a hand-built one)
  presetBaseline: null,   // JSON of that preset as loaded, to tell "modified" apart
  network: null,
  networkLabel: "",
  version: 0,             // bumped on every edit; results remember the version they were computed for
  derInfo: {},            // bus id -> der_info (incl. derived control params) from /topology
  positions: {},          // bus id -> {x, y} in diagram units, shared by every network plot
  pf: null,               // {version, single|batch results}
  modal: null,            // {version, data}
  emtNames: null,         // {version, data} -- /states response
};

// --- Events ------------------------------------------------------------------
const bus = new EventTarget();
function on(evt, fn) { bus.addEventListener(evt, e => fn(e.detail)); }
function emit(evt, detail) { bus.dispatchEvent(new CustomEvent(evt, { detail })); }

// --- API -----------------------------------------------------------------------
// FastAPI's automatic 422 for an invalid body has `detail` as a list of
// {loc, msg} objects; every hand-written HTTPException uses a plain string.
function formatDetail(detail) {
  if (Array.isArray(detail)) return detail.map(e => `${(e.loc || []).join(".")}: ${e.msg}`).join("; ");
  return detail;
}

async function api(path, opts) {
  const res = await fetch(path, opts);
  if (!res.ok) {
    const body = await res.json().catch(() => ({ detail: res.statusText }));
    throw new Error(formatDetail(body.detail) || `HTTP ${res.status}`);
  }
  return res.json();
}

function jsonPost(body) {
  return { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
}

// Every analysis call: POST /api/network/<kind> with the current network.
function netPost(kind, extra = {}) {
  if (!state.network) return Promise.reject(new Error("no network loaded -- pick or build one on the Network page"));
  return api(`/api/network/${kind}`, jsonPost({ network: state.network, ...extra }));
}

// --- Network lifecycle ---------------------------------------------------------
function setNetwork(network, { label, presetId = null, positions = null } = {}) {
  state.network = network;
  state.presetId = presetId;
  state.presetBaseline = presetId ? JSON.stringify(network) : null;
  state.networkLabel = label || network.name;
  state.positions = positions || {};
  state.needsLayout = !positions;
  state.derInfo = {};
  bumpVersion();
  emit("network:loaded");
  refreshDerivedNetworkInfo(true);
}

function bumpVersion() {
  state.version += 1;
  state.pf = null;
  state.modal = null;
  state.emtNames = null;
}

// Call after any in-place edit of state.network.
function networkChanged() {
  bumpVersion();
  emit("network:changed");
  refreshDerivedNetworkInfo(false);
}

// Which elements open breakers leave in service -- the same rules as the
// server's network.breakers.service_state: a line/transformer with either
// breaker open, a load/unit with its breaker open, and anything cut off from
// the slack's bus are out of service.
function serviceState(net) {
  const out = { slackConnected: true, energized: new Set(), lines: [], transformers: [], loads: [], units: {} };
  if (!net) return out;
  const slack = net.der_units.find(d => d.bus_type === "slack");
  const adj = new Map(net.buses.map(b => [b.id, []]));
  const link = (a, b) => { adj.get(a)?.push(b); adj.get(b)?.push(a); };
  net.lines.forEach(l => { if (l.from_closed !== false && l.to_closed !== false) link(l.from_bus, l.to_bus); });
  net.transformers.forEach(t => { if (t.hv_closed !== false && t.lv_closed !== false) link(t.hv_bus, t.lv_bus); });
  out.slackConnected = !slack || slack.closed !== false;
  if (slack && out.slackConnected && adj.has(slack.bus)) {
    const queue = [slack.bus];
    out.energized.add(slack.bus);
    while (queue.length) for (const n of adj.get(queue.shift()) || []) if (!out.energized.has(n)) { out.energized.add(n); queue.push(n); }
  } else if (!slack) net.buses.forEach(b => out.energized.add(b.id));  // no slack yet: nothing to grey out
  const on = b => out.energized.has(b);
  out.lines = net.lines.map(l => l.from_closed !== false && l.to_closed !== false && on(l.from_bus));
  out.transformers = net.transformers.map(t => t.hv_closed !== false && t.lv_closed !== false && on(t.hv_bus));
  out.loads = net.loads.map(l => l.closed !== false && on(l.bus));
  net.der_units.forEach(d => { out.units[d.id] = d.closed !== false && on(d.bus); });
  return out;
}

// Breakers: {kind: "line"|"transformer", index, end: "from"|"to"|"hv"|"lv"} |
// {kind: "load", index} | {kind: "unit", id}. Field holding its state:
function breakerField(br) {
  return br.kind === "load" || br.kind === "unit" ? "closed" : `${br.end}_closed`;
}
function breakerTarget(net, br) {
  if (br.kind === "line") return net.lines[br.index];
  if (br.kind === "transformer") return net.transformers[br.index];
  if (br.kind === "load") return net.loads[br.index];
  return net.der_units.find(d => d.id === br.id);
}
function breakerLabel(net, br) {
  const o = breakerTarget(net, br);
  if (!o) return "";
  if (br.kind === "line") return `Line #${br.index} (${o.from_bus} → ${o.to_bus}), ${br.end === "from" ? `from-bus end (bus ${o.from_bus})` : `to-bus end (bus ${o.to_bus})`}`;
  if (br.kind === "transformer") return `Transformer #${br.index} (${o.hv_bus} → ${o.lv_bus}), ${br.end === "hv" ? `HV end (bus ${o.hv_bus})` : `LV end (bus ${o.lv_bus})`}`;
  if (br.kind === "load") return `Load #${br.index} at bus ${o.bus}`;
  return `${UNIT_NAME[o.unit_type] || "Unit"} ${o.id} at bus ${o.bus}`;
}
// Toggles a breaker (refusing to open the slack's); returns false if refused.
function toggleBreaker(br) {
  const net = state.network, o = breakerTarget(net, br);
  if (!o) return false;
  const f = breakerField(br), closed = o[f] !== false;
  if (closed && br.kind === "unit" && o.bus_type === "slack") {
    alert("The slack unit's breaker can't be opened: it is the reference of the power flow and of the dynamic models. Make another unit the slack first.");
    return false;
  }
  o[f] = !closed;
  networkChanged();
  return true;
}

function isModified() {
  return state.presetBaseline !== null && JSON.stringify(state.network) !== state.presetBaseline;
}

// der_info (dispatch + derived control parameters) comes from /topology;
// refreshed after edits, debounced. With seedPositions, buses without a
// position also get one from the server's Kamada-Kawai layout.
let derivedTimer = null;
function refreshDerivedNetworkInfo(seedPositions) {
  clearTimeout(derivedTimer);
  const version = state.version;
  derivedTimer = setTimeout(async () => {
    try {
      const topo = await netPost("topology");
      if (version !== state.version) return;
      state.derInfo = Object.fromEntries(topo.nodes.filter(n => n.der_info).map(n => [n.id, n.der_info]));
      if (seedPositions) { applyLayout(topo, state.needsLayout); state.needsLayout = false; }
      emit("network:derived");
    } catch { /* an invalid mid-edit network has no topology yet -- keep the last one */ }
  }, seedPositions ? 0 : 450);
}

const DIAGRAM_W = 1000, DIAGRAM_H = 600, DIAGRAM_PAD = 70;
function applyLayout(topo, overwriteAll) {
  const xs = topo.nodes.map(n => n.x), ys = topo.nodes.map(n => n.y);
  if (!xs.length) return;
  const xMin = Math.min(...xs), xMax = Math.max(...xs), yMin = Math.min(...ys), yMax = Math.max(...ys);
  const sx = v => DIAGRAM_PAD + ((v - xMin) / ((xMax - xMin) || 1)) * (DIAGRAM_W - 2 * DIAGRAM_PAD);
  const sy = v => DIAGRAM_PAD + ((v - yMin) / ((yMax - yMin) || 1)) * (DIAGRAM_H - 2 * DIAGRAM_PAD);
  topo.nodes.forEach(n => {
    if (overwriteAll || !state.positions[n.id]) state.positions[n.id] = { x: sx(n.x), y: sy(n.y) };
  });
  emit("network:layout");
}

async function autoLayout() {
  const topo = await netPost("topology");
  applyLayout(topo, true);
}

// A bus added outside the canvas (e.g. from the table view) gets a position
// next to a connected bus, or near the centre.
function ensurePositions() {
  const net = state.network;
  if (!net) return;
  for (const b of net.buses) {
    if (state.positions[b.id]) continue;
    const nb = neighbours(b.id).find(id => state.positions[id]);
    const base = nb !== undefined ? state.positions[nb] : { x: DIAGRAM_W / 2, y: DIAGRAM_H / 2 };
    state.positions[b.id] = { x: base.x + 40 + Math.random() * 30, y: base.y + (Math.random() - 0.5) * 60 };
  }
}

function neighbours(busId) {
  const net = state.network, out = [];
  net.lines.forEach(l => { if (l.from_bus === busId) out.push(l.to_bus); else if (l.to_bus === busId) out.push(l.from_bus); });
  net.transformers.forEach(t => { if (t.hv_bus === busId) out.push(t.lv_bus); else if (t.lv_bus === busId) out.push(t.hv_bus); });
  return out;
}

// --- Formatting ------------------------------------------------------------------
function fmt(n, digits = 4) {
  if (typeof n !== "number" || Number.isNaN(n)) return "-";
  return n.toFixed(digits);
}
function fmtSmart(v) {
  if (v === null || v === undefined) return "-";
  if (typeof v !== "number") return esc(String(v));
  if (v === 0) return "0";
  if (Number.isInteger(v)) return String(v);
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(3);
  return v.toFixed(a >= 100 ? 2 : 4);
}
function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}
function seriesColor(i) { return `var(--series-${(i % 8) + 1})`; }
// Resolved hex values, for places that need a real colour (canvas-free SVG
// attributes work with var(), but luminance maths doesn't).
const SERIES_HEX = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"];

const UNIT_LABEL = { sm: "SM", gfm: "GFM", gfl: "GFL", infinite_bus: "IB" };
const UNIT_COLOR = { sm: "#2a78d6", gfm: "#eb6834", gfl: "#1baf7a", infinite_bus: "#4a3aa7" };
const UNIT_NAME = { sm: "Synchronous machine", gfm: "Grid-forming converter", gfl: "Grid-following converter", infinite_bus: "Infinite bus" };

function el(html) {
  const t = document.createElement("template");
  t.innerHTML = html.trim();
  return t.content.firstElementChild;
}
function $(sel, root = document) { return root.querySelector(sel); }
function $$(sel, root = document) { return Array.from(root.querySelectorAll(sel)); }

function errorHtml(e) { return `<div class="notice err-bg"><span><b class="bad">Error:</b> ${esc(e.message || e)}</span></div>`; }

function setSpinner(node, text) {
  node.textContent = text || "";
  node.classList.toggle("on", !!text);
}

// --- Floating tooltip ------------------------------------------------------------
const tooltipEl = document.createElement("div");
tooltipEl.id = "tooltip";
document.addEventListener("DOMContentLoaded", () => document.body.appendChild(tooltipEl));
function showTooltip(html, evt) {
  tooltipEl.innerHTML = html;
  tooltipEl.style.display = "block";
  moveTooltip(evt);
}
function moveTooltip(evt) {
  if (!evt) return;
  const pad = 14, r = tooltipEl.getBoundingClientRect();
  let x = evt.clientX + pad, y = evt.clientY + pad;
  if (x + r.width > window.innerWidth - 8) x = evt.clientX - pad - r.width;
  if (y + r.height > window.innerHeight - 8) y = evt.clientY - pad - r.height;
  tooltipEl.style.left = `${Math.max(4, x)}px`;
  tooltipEl.style.top = `${Math.max(4, y)}px`;
}
function hideTooltip() { tooltipEl.style.display = "none"; }
function ttTable(rows) {
  return `<table>${rows.map(r => r.sep
    ? `<tr class="tt-sep"><td colspan="2">${esc(r.sep)}</td></tr>`
    : `<tr><td>${r.sw ? `<span class="tt-swatch" style="background:${r.sw}"></span>` : ""}${esc(r[0])}</td><td>${r[1]}</td></tr>`).join("")}</table>`;
}

// --- Signal picker ---------------------------------------------------------------
// A chip list + searchable checkbox popover over a (possibly grouped) list of
// names -- replaces native <select multiple>, which is unusable for the
// 50-200 state/output names these models have.
class SignalPicker {
  constructor(host, { options = [], selected = [], placeholder = "add signal", colors = false, onChange = null, display = null } = {}) {
    this.host = host;
    this.options = options;      // [{name, group?}] or [name]
    this.selected = [...selected];
    this.placeholder = placeholder;
    this.colors = colors;
    this.display = display || (n => n);  // how a selected name reads on its chip
    this.onChange = onChange;
    this.render();
  }
  setOptions(options, keepSelection = true) {
    this.options = options;
    const names = new Set(this.optionNames());
    this.selected = keepSelection ? this.selected.filter(n => names.has(n)) : [];
    this.render();
  }
  optionNames() { return this.options.map(o => typeof o === "string" ? o : o.name); }
  get() { return [...this.selected]; }
  set(list) { this.selected = [...list]; this.render(); }
  changed() { this.render(); if (this.onChange) this.onChange(this.get()); }
  render() {
    this.host.innerHTML = "";
    const box = el(`<div class="picker"></div>`);
    this.selected.forEach((n, i) => {
      const chip = el(`<span class="pchip">${this.colors ? `<span class="sw" style="background:${seriesColor(i)}"></span>` : ""}${esc(this.display(n))}<button type="button" aria-label="Remove ${esc(n)}">&times;</button></span>`);
      chip.querySelector("button").addEventListener("click", ev => {
        ev.stopPropagation();
        this.selected = this.selected.filter(x => x !== n);
        this.changed();
      });
      box.appendChild(chip);
    });
    const add = el(`<button type="button" class="padd">+ ${esc(this.placeholder)}</button>`);
    add.addEventListener("click", ev => { ev.stopPropagation(); this.openPopover(box); });
    box.appendChild(add);
    this.host.appendChild(box);
  }
  openPopover(box) {
    closePopovers();
    const pop = el(`<div class="picker-pop" role="dialog">
      <input type="search" placeholder="Filter…" aria-label="Filter signals">
      <div class="pactions"><button type="button" class="ghost small" data-a="all">Select shown</button><button type="button" class="ghost small" data-a="none">Clear shown</button></div>
      <div class="plist"></div></div>`);
    pop.addEventListener("click", ev => ev.stopPropagation());
    const list = pop.querySelector(".plist");
    const input = pop.querySelector("input");
    const shown = () => {
      const q = input.value.trim().toLowerCase();
      return this.options.map(o => typeof o === "string" ? { name: o } : o).filter(o => !q || o.name.toLowerCase().includes(q) || (o.label || "").toLowerCase().includes(q) || (o.group || "").toLowerCase().includes(q));
    };
    const draw = () => {
      let lastGroup = null, html = "";
      const items = shown();
      items.slice(0, 400).forEach(o => {
        if (o.group && o.group !== lastGroup) { html += `<div class="pgroup">${esc(o.group)}</div>`; lastGroup = o.group; }
        html += `<label class="popt"><input type="checkbox" value="${esc(o.name)}"${this.selected.includes(o.name) ? " checked" : ""}>${o.label ? `<span>${esc(o.label)} <span class="muted">${esc(this.display(o.name))}</span></span>` : esc(this.display(o.name))}</label>`;
      });
      if (items.length > 400) html += `<div class="pgroup">${items.length - 400} more — refine the filter</div>`;
      if (!items.length) html = `<p class="empty">No match.</p>`;
      list.innerHTML = html;
    };
    list.addEventListener("change", ev => {
      const n = ev.target.value;
      if (ev.target.checked) { if (!this.selected.includes(n)) this.selected.push(n); }
      else this.selected = this.selected.filter(x => x !== n);
      this.renderChipsOnly();
      if (this.onChange) this.onChange(this.get());
    });
    pop.querySelector('[data-a="all"]').addEventListener("click", () => {
      shown().slice(0, 400).forEach(o => { if (!this.selected.includes(o.name)) this.selected.push(o.name); });
      draw(); this.renderChipsOnly(); if (this.onChange) this.onChange(this.get());
    });
    pop.querySelector('[data-a="none"]').addEventListener("click", () => {
      const s = new Set(shown().map(o => o.name));
      this.selected = this.selected.filter(n => !s.has(n));
      draw(); this.renderChipsOnly(); if (this.onChange) this.onChange(this.get());
    });
    input.addEventListener("input", draw);
    draw();
    box.appendChild(pop);
    input.focus();
  }
  // Re-render chips while keeping the open popover in place.
  renderChipsOnly() {
    const pop = this.host.querySelector(".picker-pop");
    const scroll = pop ? pop.querySelector(".plist").scrollTop : 0;
    this.render();
    if (pop) {
      const nb = this.host.querySelector(".picker");
      nb.appendChild(pop);
      pop.querySelector(".plist").scrollTop = scroll;
    }
  }
}
function closePopovers() { $$(".picker-pop").forEach(p => p.remove()); }
document.addEventListener("click", closePopovers);
document.addEventListener("keydown", e => { if (e.key === "Escape") closePopovers(); });

// Network-context strip shown at the top of analysis pages.
function networkContextHtml() {
  if (!state.network) {
    return `<div class="notice warn-bg"><span>No network loaded yet — <a href="#/network">choose a preset or build one</a> on the Network page.</span></div>`;
  }
  const n = state.network;
  const sv = serviceState(n);
  const out = sv.lines.filter(x => !x).length + sv.transformers.filter(x => !x).length + sv.loads.filter(x => !x).length
    + Object.values(sv.units).filter(x => !x).length;
  return `<div class="card" style="padding:0.7rem 1rem"><div class="controls" style="align-items:center">
    <span class="field-label">Network</span><b style="font-size:0.9rem">${esc(state.networkLabel)}</b>
    ${isModified() ? '<span class="badge warn">modified</span>' : ""}
    ${out ? `<span class="badge crit" title="Open breakers: these elements (and anything they cut off from the slack) are left out of every analysis">${out} element${out > 1 ? "s" : ""} out of service</span>` : ""}
    <span class="muted" style="font-size:0.8rem;font-family:var(--font-mono)">${n.buses.length} buses · ${n.lines.length} lines · ${n.transformers.length} trafos · ${n.loads.length} loads · ${n.der_units.length} DER</span>
    <span style="flex:1"></span><a href="#/network" style="font-size:0.8rem">Change / edit →</a></div></div>`;
}
