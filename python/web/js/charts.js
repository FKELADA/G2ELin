// Charts: a general time-series line chart (shared by modal free/step
// response, EMT and batch power flow), colour scales, and the modal-analysis
// specific charts (eigenvalue map, heatmap table, bar and polar charts).

const SVG_NS = "http://www.w3.org/2000/svg";
function svgEl(tag, attrs) {
  const n = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs || {})) n.setAttribute(k, v);
  return n;
}

// --- Colour scales -------------------------------------------------------------
function lerpStops(stops, t) {
  t = Math.min(1, Math.max(0, Number.isFinite(t) ? t : 0.5));
  const seg = (stops.length - 1) * t;
  const i = Math.min(stops.length - 2, Math.floor(seg));
  const f = seg - i;
  const a = stops[i], b = stops[i + 1];
  return `rgb(${a.map((c, k) => Math.round(c + (b[k] - c) * f)).join(",")})`;
}
// Low = blue, neutral in the middle, high = red (bus voltages).
const DIVERGING = [[33, 84, 160], [104, 158, 222], [226, 228, 230], [240, 140, 106], [192, 40, 45]];
// Light -> dark red (line loading/flows).
const SEQ_HEAT = [[250, 233, 190], [246, 178, 94], [226, 108, 50], [178, 34, 38]];
// Light -> dark blue (heatmap tables).
const SEQ_BLUE = [[235, 242, 252], [120, 170, 230], [13, 54, 107]];
const divergingColor = t => lerpStops(DIVERGING, t);
const heatColor = t => lerpStops(SEQ_HEAT, t);
const sequentialBlue = t => lerpStops(SEQ_BLUE, t);
function gradientCss(stops) { return `linear-gradient(90deg, ${stops.map(s => `rgb(${s.join(",")})`).join(", ")})`; }
function isDark(rgbStr) {
  const m = rgbStr.match(/\d+/g);
  if (!m) return false;
  const [r, g, b] = m.map(Number);
  return 0.299 * r + 0.587 * g + 0.114 * b < 140;
}
function gradientLegendHtml(title, stops, lo, hi, digits = 3, unit = "") {
  return `<span class="gradient-legend"><span class="gl-title">${esc(title)}</span>${fmt(lo, digits)}${unit}
    <span class="bar" style="background:${gradientCss(stops)}"></span>${fmt(hi, digits)}${unit}</span>`;
}

// --- Line chart ------------------------------------------------------------------
// series: [{name, t, y, color, dash?: bool, width?}]. Each series carries its
// own time grid (EMT trajectories, their linear overlay and live-traced
// points are all sampled differently). opts.group: an object shared by
// several charts to sync their crosshairs (subplots of one channel).
function niceTicks(lo, hi, n = 5) {
  const span = hi - lo || Math.abs(hi) || 1;
  const step0 = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return { ticks: out, step };
}
function tickLabel(v, step) {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(1);
  const d = Math.max(0, Math.min(6, -Math.floor(Math.log10(step)) + 1));
  return v.toFixed(d);
}
function nearestIndex(arr, x) {
  let lo = 0, hi = arr.length - 1;
  if (hi < 0) return -1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (arr[mid] < x) lo = mid; else hi = mid;
  }
  return Math.abs(arr[lo] - x) <= Math.abs(arr[hi] - x) ? lo : hi;
}

function lineChart(series, opts = {}) {
  const W = opts.width || 780, H = opts.height || 260;
  const PADL = 64, PADR = 18, PADT = 14, PADB = opts.xTitle ? 42 : 30;
  const xUnit = opts.xUnit ?? "s";
  const wrap = document.createElement("div");
  wrap.className = "chart";
  if (opts.title) wrap.appendChild(el(`<p class="chart-title">${esc(opts.title)}</p>`));
  const live = series.filter(s => s.t.length);
  if (!live.length) { wrap.appendChild(el(`<p class="empty">No data.</p>`)); return wrap; }

  // The full extent of the data, measured once. A zoom is a window over
  // this, never a re-measure, so zooming out always lands back exactly here.
  const full = (() => {
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    live.forEach(s => {
      xMin = Math.min(xMin, s.t[0]); xMax = Math.max(xMax, s.t[s.t.length - 1]);
      s.y.forEach(v => { if (Number.isFinite(v)) { yMin = Math.min(yMin, v); yMax = Math.max(yMax, v); } });
    });
    if (opts.xMax !== undefined) xMax = Math.max(xMax, opts.xMax);
    if (!Number.isFinite(yMin)) { yMin = -1; yMax = 1; }
    const pad = (yMax - yMin) * 0.08 || Math.max(Math.abs(yMax), 1e-6) * 0.05 || 0.01;
    return { xMin, xMax, y0: yMin - pad, y1: yMax + pad };
  })();

  // State that has to outlive a redraw: which signals are switched off, and
  // the window a zoom has set.
  const hidden = new Set();
  let win = null;
  let impl = null;                                   // the current drawing
  const handle = { showAt: t => impl && impl.showAt(t), hide: () => impl && impl.hide() };
  if (opts.group) (opts.group.charts ||= []).push(handle);

  const zoomable = opts.zoomable !== false;
  let mode = opts.zoomMode || "x";
  let bar = null;

  const draw = () => {
    const { xMin, xMax, y0, y1 } = win || full;
    const sx = v => PADL + ((v - xMin) / ((xMax - xMin) || 1)) * (W - PADL - PADR);
    const sy = v => (H - PADB) - ((v - y0) / ((y1 - y0) || 1)) * (H - PADT - PADB);
    const svg = svgEl("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": opts.title || "Line chart" });

    const yt = niceTicks(y0, y1, 5);
    yt.ticks.forEach(v => {
      if (v < y0 || v > y1) return;
      const y = sy(v).toFixed(1);
      svg.appendChild(svgEl("line", { x1: PADL, x2: W - PADR, y1: y, y2: y, stroke: "var(--gridline)" }));
      const tx = svgEl("text", { x: PADL - 8, y: +y + 3, "text-anchor": "end", class: "axis-text" });
      tx.textContent = tickLabel(v, yt.step);
      svg.appendChild(tx);
    });
    const xt = niceTicks(xMin, xMax, 6);
    xt.ticks.forEach(v => {
      if (v < xMin || v > xMax) return;
      const x = sx(v).toFixed(1);
      svg.appendChild(svgEl("line", { x1: x, x2: x, y1: H - PADB, y2: H - PADB + 4, stroke: "var(--baseline)" }));
      const tx = svgEl("text", { x, y: H - PADB + 16, "text-anchor": "middle", class: "axis-text" });
      const tf = opts.xTickFormat || opts.xFormat;
      tx.textContent = tf ? tf(v) : tickLabel(v, xt.step) + xUnit;
      svg.appendChild(tx);
    });
    svg.appendChild(svgEl("line", { x1: PADL, x2: W - PADR, y1: H - PADB, y2: H - PADB, stroke: "var(--baseline)" }));
    if (opts.xTitle) {
      const tx = svgEl("text", { x: (PADL + W - PADR) / 2, y: H - 6, "text-anchor": "middle", class: "axis-title" });
      tx.textContent = opts.xTitle;
      svg.appendChild(tx);
    }
    if (opts.yTitle) {
      const tx = svgEl("text", { x: 12, y: (PADT + H - PADB) / 2, "text-anchor": "middle", class: "axis-title", transform: `rotate(-90 12 ${(PADT + H - PADB) / 2})` });
      tx.textContent = opts.yTitle;
      svg.appendChild(tx);
    }
    (opts.markers || []).forEach(m => {
      if (m.x < xMin || m.x > xMax) return;
      const x = sx(m.x).toFixed(1);
      svg.appendChild(svgEl("line", { x1: x, x2: x, y1: PADT, y2: H - PADB, stroke: "var(--text-muted)", "stroke-dasharray": "2 3" }));
      const tx = svgEl("text", { x: +x + 4, y: PADT + 9, class: "axis-text" });
      tx.textContent = m.label;
      svg.appendChild(tx);
    });

    // Everything is clipped to the plot area, so a zoomed-in trace does not
    // draw over the axes.
    const clipId = `clip-${Math.random().toString(36).slice(2)}`;
    const defs = svgEl("defs", {});
    const cp = svgEl("clipPath", { id: clipId });
    cp.appendChild(svgEl("rect", { x: PADL, y: PADT, width: W - PADL - PADR, height: H - PADT - PADB }));
    defs.appendChild(cp);
    svg.appendChild(defs);
    const plot = svgEl("g", { "clip-path": `url(#${clipId})` });
    svg.appendChild(plot);

    // Downsample very long series for drawing (hover still reads full data).
    const pathEls = live.map(s => {
      const n = s.t.length, stride = Math.max(1, Math.floor(n / 2500));
      let d = "", gap = true;  // a missing value (null/NaN) breaks the line
      for (let j = 0; j < n; j += stride) {
        if (!Number.isFinite(s.y[j])) { gap = true; continue; }
        d += `${gap ? "M" : "L"}${sx(s.t[j]).toFixed(1)},${sy(s.y[j]).toFixed(1)}`;
        gap = false;
      }
      if (stride > 1 && n && Number.isFinite(s.y[n - 1]) && !gap) d += `L${sx(s.t[n - 1]).toFixed(1)},${sy(s.y[n - 1]).toFixed(1)}`;
      const node = svgEl("path", {
        d, fill: "none", stroke: s.color, "stroke-width": s.width || (s.dash ? 2.2 : 1.8), "stroke-linejoin": "round",
        "stroke-linecap": "round", "data-series": s.name, ...(s.dash ? { "stroke-dasharray": "1.5 4" } : {}),
      });
      plot.appendChild(node);
      return node;
    });

    const hoverLine = svgEl("line", { class: "hover-line", y1: PADT, y2: H - PADB });
    hoverLine.style.display = "none";
    svg.appendChild(hoverLine);
    const dots = live.map(s => {
      const c = svgEl("circle", { r: 3.5, class: "hover-dot", fill: s.color });
      c.style.display = "none";
      plot.appendChild(c);
      return c;
    });
    const band = svgEl("rect", { class: "zoom-band", y: PADT, height: H - PADT - PADB, width: 0 });
    band.style.display = "none";
    svg.appendChild(band);
    const capture = svgEl("rect", { class: "capture", x: PADL, y: PADT, width: W - PADL - PADR, height: H - PADT - PADB });
    svg.appendChild(capture);

    const applyHidden = () => live.forEach((s, i) => {
      const off = hidden.has(s.name);
      pathEls[i].style.display = off ? "none" : "";
      if (off) dots[i].style.display = "none";
    });
    applyHidden();

    const showAt = tx => {
      const xPix = sx(tx);
      hoverLine.setAttribute("x1", xPix); hoverLine.setAttribute("x2", xPix);
      hoverLine.style.display = (tx >= xMin && tx <= xMax) ? "" : "none";
      return live.map((s, i) => {
        if (hidden.has(s.name)) { dots[i].style.display = "none"; return null; }
        const k = nearestIndex(s.t, tx);
        const v = s.y[k];
        dots[i].setAttribute("cx", sx(s.t[k])); dots[i].setAttribute("cy", sy(v));
        dots[i].style.display = Number.isFinite(v) ? "" : "none";
        return { s, v };
      }).filter(Boolean);
    };
    const hide = () => { hoverLine.style.display = "none"; dots.forEach(d => { d.style.display = "none"; }); };
    impl = { showAt, hide, applyHidden };

    // Pixel -> data, through the SVG's own transform so it survives any
    // CSS scaling of the chart.
    const at = evt => {
      const ctm = svg.getScreenCTM();
      if (!ctm) return null;
      const pt = svg.createSVGPoint(); pt.x = evt.clientX; pt.y = evt.clientY;
      const q = pt.matrixTransform(ctm.inverse());
      return {
        px: q.x, py: q.y,
        x: xMin + ((q.x - PADL) / (W - PADL - PADR)) * (xMax - xMin),
        y: y0 + ((H - PADB - q.y) / (H - PADT - PADB)) * (y1 - y0),
      };
    };

    let drag = null;
    capture.addEventListener("mousemove", evt => {
      const q = at(evt);
      if (!q) return;
      if (drag) {
        // The rubber band shows exactly what the release will keep, so a
        // horizontal zoom draws a full-height band and not a box the user
        // then finds was interpreted differently.
        const x0 = Math.min(drag.px, q.px), x1p = Math.max(drag.px, q.px);
        const yA = Math.min(drag.py, q.py), yB = Math.max(drag.py, q.py);
        const horiz = mode !== "y", vert = mode !== "x";
        band.setAttribute("x", horiz ? x0 : PADL);
        band.setAttribute("width", horiz ? Math.max(0, x1p - x0) : W - PADL - PADR);
        band.setAttribute("y", vert ? yA : PADT);
        band.setAttribute("height", vert ? Math.max(0, yB - yA) : H - PADT - PADB);
        band.style.display = "";
        hideTooltip();
        return;
      }
      const vals = showAt(q.x);
      if (opts.group) opts.group.charts.forEach(c => { if (c !== handle) c.showAt(q.x); });
      const xLab = opts.xFormat ? opts.xFormat(q.x) : `${fmt(q.x, 4)} ${xUnit}`;
      const rows = vals.slice(0, 16).map(({ s, v }) => [s.name + (s.dash ? " (lin.)" : ""), fmtSmart(v)]);
      showTooltip(`<div class="tt-title">${esc(opts.xLabel || "t")} = ${esc(xLab)}</div>` +
        ttTable(rows.map((r, i) => Object.assign(r, { sw: vals[i].s.color }))) +
        (vals.length > 16 ? `<div class="muted">+${vals.length - 16} more</div>` : ""), evt);
    });
    capture.addEventListener("mouseleave", () => {
      hide(); hideTooltip();
      if (opts.group) opts.group.charts.forEach(c => c.hide());
    });

    if (zoomable) {
      capture.addEventListener("mousedown", evt => {
        if (evt.button !== 0) return;
        evt.preventDefault();
        drag = at(evt);
        hide(); hideTooltip();
      });
      capture.addEventListener("mouseup", evt => {
        const from = drag, to = at(evt);
        drag = null;
        band.style.display = "none";
        if (!from || !to) return;
        // A click, or a band too thin to mean anything, is not a zoom: a
        // stray click would otherwise blow the axes up to a single point.
        const next = { ...(win || full) };
        const wideX = Math.abs(to.x - from.x) > (next.xMax - next.xMin) * 0.01;
        const wideY = Math.abs(to.y - from.y) > (next.y1 - next.y0) * 0.01;
        const takeX = mode !== "y" && wideX, takeY = mode !== "x" && wideY;
        if (!takeX && !takeY) return;
        if (takeX) { next.xMin = Math.min(from.x, to.x); next.xMax = Math.max(from.x, to.x); }
        if (takeY) { next.y0 = Math.min(from.y, to.y); next.y1 = Math.max(from.y, to.y); }
        win = next;
        redraw();
      });
      capture.addEventListener("mouseleave", () => { drag = null; band.style.display = "none"; });
      capture.addEventListener("dblclick", () => { if (win) { win = null; redraw(); } });
    }

    return svg;
  };

  let svgNode = null;
  const redraw = () => {
    const next = draw();
    if (svgNode) wrap.replaceChild(next, svgNode); else wrap.appendChild(next);
    svgNode = next;
    if (bar) bar.refresh();
  };

  if (zoomable) {
    const modes = [["x", "X", "Drag to window the time axis"],
                   ["y", "Y", "Drag to window the value axis"],
                   ["box", "Box", "Drag a box to window both"]];
    bar = el(`<div class="chart-zoombar">${modes.map(([id, label, tip]) =>
      `<button type="button" class="ghost small${id === mode ? " on" : ""}" data-zmode="${id}" title="${tip}">${label}</button>`).join("")}
      <button type="button" class="ghost small" data-zreset title="Back to the full extent (or double-click the plot)">Reset</button>
      <span class="muted zoom-note"></span></div>`);
    bar.refresh = () => {
      bar.querySelectorAll("[data-zmode]").forEach(b => b.classList.toggle("on", b.dataset.zmode === mode));
      bar.querySelector("[data-zreset]").disabled = !win;
      bar.querySelector(".zoom-note").textContent = win ? "zoomed" : "";
    };
    bar.querySelectorAll("[data-zmode]").forEach(b => b.addEventListener("click", () => {
      mode = b.dataset.zmode; bar.refresh();
    }));
    bar.querySelector("[data-zreset]").addEventListener("click", () => { win = null; redraw(); });
    wrap.appendChild(bar);
  }

  wrap.__chart = {
    series: live,
    names: [...new Set(live.map(s => s.name))],
    isHidden: name => hidden.has(name),
    toggle(name) { hidden.has(name) ? hidden.delete(name) : hidden.add(name); impl.applyHidden(); },
    isolate(name) {
      const others = this.names.filter(n => n !== name);
      const already = !hidden.has(name) && others.every(n => hidden.has(n));
      hidden.clear();
      if (!already) others.forEach(n => hidden.add(n));   // a second double-click brings them back
      impl.applyHidden();
    },
  };

  redraw();
  if (bar) bar.refresh();
  return wrap;
}

// A legend whose entries switch their signal off and on: click to hide one,
// double-click to show only it (and again to bring the rest back).
function legendHtml(items) {
  return `<div class="legend legend-toggle">${items.map(it =>
    `<span data-series="${esc(it.name)}" role="button" tabindex="0" title="Click to hide this signal · double-click to show only it"><span class="line-swatch${it.dash ? " dotted" : ""}" style="border-color:${it.color}"></span>${esc(it.name)}</span>`).join("")}</div>`;
}

// The charts a legend entry speaks for: the ones in its own card (so a scope
// split over stacked subplots switches all of them together).
function legendCharts(item) {
  const scope = item.closest(".card, .channel, aside") || document;
  return [...scope.querySelectorAll(".chart")].filter(c => c.__chart);
}
function legendApply(item, fn) {
  const charts = legendCharts(item);
  charts.forEach(c => fn(c.__chart, item.dataset.series));
  const first = charts[0];
  if (!first) return;
  item.parentElement.querySelectorAll("[data-series]").forEach(sp => {
    sp.classList.toggle("off", first.__chart.isHidden(sp.dataset.series));
  });
}
let _legendClick = null;
document.addEventListener("click", ev => {
  const item = ev.target.closest(".legend-toggle [data-series]");
  if (!item || ev.detail > 1) return;
  clearTimeout(_legendClick);
  _legendClick = setTimeout(() => { _legendClick = null; legendApply(item, (c, n) => c.toggle(n)); }, 220);
});
document.addEventListener("dblclick", ev => {
  const item = ev.target.closest(".legend-toggle [data-series]");
  if (!item) return;
  clearTimeout(_legendClick);   // the click that came with it doesn't count
  _legendClick = null;
  legendApply(item, (c, n) => c.isolate(n));
});

// --- Pan/zoom by viewBox ---------------------------------------------------------
// Wheel zooms around the cursor, background drag pans. `canStartPan(evt)` lets
// a caller (the network editor) keep drags on elements for itself.
function attachSvgZoomPan(svg, { canStartPan = () => true, onChange = null } = {}) {
  const vb = svg.viewBox.baseVal;
  const initial = { x: vb.x, y: vb.y, width: vb.width, height: vb.height };
  const cur = { ...initial };
  const apply = () => { svg.setAttribute("viewBox", `${cur.x} ${cur.y} ${cur.width} ${cur.height}`); if (onChange) onChange(cur); };
  const zoomBy = (factor, cx, cy) => {
    const newW = Math.min(initial.width * 6, Math.max(initial.width * 0.05, cur.width * factor));
    const f = newW / cur.width;
    cx ??= cur.x + cur.width / 2; cy ??= cur.y + cur.height / 2;
    cur.x = cx - (cx - cur.x) * f; cur.y = cy - (cy - cur.y) * f;
    cur.width = newW; cur.height *= f;
    apply();
  };
  svg.addEventListener("wheel", evt => {
    evt.preventDefault();
    const r = svg.getBoundingClientRect();
    const mx = cur.x + (evt.clientX - r.left) * cur.width / r.width;
    const my = cur.y + (evt.clientY - r.top) * cur.height / r.height;
    zoomBy(evt.deltaY > 0 ? 1.15 : 1 / 1.15, mx, my);
  }, { passive: false });
  let start = null;
  svg.addEventListener("pointerdown", evt => {
    if (evt.button !== 0 || !canStartPan(evt)) return;
    start = { x: evt.clientX, y: evt.clientY, vb: { ...cur }, moved: false, id: evt.pointerId };
  });
  svg.addEventListener("pointermove", evt => {
    if (!start) return;
    if (!start.moved && Math.hypot(evt.clientX - start.x, evt.clientY - start.y) < 4) return;
    if (!start.moved) { start.moved = true; svg.setPointerCapture(start.id); svg.classList.add("panning", "dragging"); }
    const r = svg.getBoundingClientRect();
    cur.x = start.vb.x - (evt.clientX - start.x) * start.vb.width / r.width;
    cur.y = start.vb.y - (evt.clientY - start.y) * start.vb.height / r.height;
    apply();
  });
  const end = () => { if (start && start.moved) svg.__justPanned = true; start = null; svg.classList.remove("panning", "dragging"); setTimeout(() => { svg.__justPanned = false; }, 0); };
  svg.addEventListener("pointerup", end);
  svg.addEventListener("pointercancel", end);
  return {
    reset: () => { Object.assign(cur, initial); apply(); },
    zoomIn: () => zoomBy(1 / 1.3), zoomOut: () => zoomBy(1.3),
    get: () => ({ ...cur }), set: v => { Object.assign(cur, v); apply(); },
  };
}

// --- Modal-analysis charts -------------------------------------------------------
function statusOf(real) {
  if (real > 1e-6) return { cls: "critical", label: "Unstable" };
  if (real > -1e-3) return { cls: "warning", label: "Marginal" };
  return { cls: "good", label: "Stable" };
}
function symlog(v) { return Math.sign(v) * Math.log10(1 + Math.abs(v)); }
function invSymlog(v) { return Math.sign(v) * (Math.pow(10, Math.abs(v)) - 1); }

// Constant-damping-ratio guide lines: straight through the origin in linear
// (real, Hz) space, sampled log-spaced in |real| and symlog-transformed like
// every point here, so they read as curves.
function dampingGuideLines(sx, sy, zeta, color) {
  const slope = Math.sqrt((1 - zeta * zeta) / (zeta * zeta));
  const reals = [];
  for (let t = -3; t <= 9; t += 0.25) reals.push(-Math.pow(10, t));
  const branch = sign => reals.map(re => `${sx(symlog(re)).toFixed(1)},${sy(symlog(sign * slope * Math.abs(re) / (2 * Math.PI))).toFixed(1)}`).join(" ");
  return `<polyline points="${branch(1)}" fill="none" stroke="${color}" stroke-width="1" stroke-dasharray="4,3"/>
    <polyline points="${branch(-1)}" fill="none" stroke="${color}" stroke-width="1" stroke-dasharray="4,3"/>`;
}

function eigenGridLines(sx, sy, xMax, yMax, W, H, PAD) {
  const N = 4;
  const ticks = n => Array.from({ length: 2 * n + 1 }, (_, i) => (i - n) / n);
  const xTicks = ticks(N).map(t => t * xMax), yTicks = ticks(N).map(t => t * yMax);
  return xTicks.filter(v => v !== 0).map(v => `<line class="eigen-gridline" x1="${sx(v).toFixed(1)}" y1="${PAD}" x2="${sx(v).toFixed(1)}" y2="${H - PAD}"/>`).join("")
    + yTicks.filter(v => v !== 0).map(v => `<line class="eigen-gridline" x1="${PAD}" y1="${sy(v).toFixed(1)}" x2="${W - PAD}" y2="${sy(v).toFixed(1)}"/>`).join("")
    + xTicks.map(v => `<text class="eigen-ticklabel" x="${sx(v).toFixed(1)}" y="${H - PAD + 13}" text-anchor="middle">${invSymlog(v).toExponential(1)}</text>`).join("")
    + yTicks.map(v => `<text class="eigen-ticklabel" x="${PAD - 6}" y="${sy(v).toFixed(1)}" text-anchor="end" dominant-baseline="middle">${invSymlog(v).toExponential(1)}</text>`).join("");
}

// One marker shape per mode category. Colour on this map already means
// stability, so shape is the channel that is free -- and it is the one that
// survives greyscale printing and colour-blindness anyway. Drawn as a path
// around (0,0) so the same code can place and scale every shape.
const MODE_MARKERS = {
  synchronisation: r => `M0,${-r * 1.25}L${r * 1.15},${r * 0.75}L${-r * 1.15},${r * 0.75}Z`,  // triangle
  control:         r => `M${-r},${-r}h${2 * r}v${2 * r}h${-2 * r}Z`,                          // square
  unit_electrical: r => `M0,${-r * 1.3}L${r * 1.3},0L0,${r * 1.3}L${-r * 1.3},0Z`,            // diamond
  network:         r => `M${-r},0a${r},${r} 0 1,0 ${2 * r},0a${r},${r} 0 1,0 ${-2 * r},0`,    // circle
  mixed:           r => `M0,${-r * 1.35}L${r * 0.42},${-r * 0.42}L${r * 1.35},0L${r * 0.42},${r * 0.42}`
                        + `L0,${r * 1.35}L${-r * 0.42},${r * 0.42}L${-r * 1.35},0L${-r * 0.42},${-r * 0.42}Z`,
  reference:       r => `M${-r},${-r}L${r},${r}M${-r},${r}L${r},${-r}`,                       // cross
};
const MODE_MARKER_FALLBACK = MODE_MARKERS.network;

function modeMarkerPath(category, r) {
  return (MODE_MARKERS[category] || MODE_MARKER_FALLBACK)(r);
}

// The legend entry for one category, as a small standalone SVG.
function modeMarkerSwatch(category) {
  // The cross has no interior, so it needs the hollow class to be stroked in
  // the series colour -- .pt's default stroke is the card background, which
  // would draw it invisibly.
  const hollow = category === "reference" ? " hollow" : "";
  return `<svg width="14" height="14" viewBox="-7 -7 14 14" style="vertical-align:-3px;margin-right:0.35em" aria-hidden="true">
    <path d="${modeMarkerPath(category, 4.4)}" class="pt good${hollow}"/></svg>`;
}

function eigenvalueMapSvg(modes, selectedMode, hidden = null) {
  const W = 820, H = 440, PAD = 50;
  const xs = modes.map(m => symlog(m.real)), ys = modes.map(m => symlog(m.imag / (2 * Math.PI)));
  const xMax = Math.max(1, ...xs.map(Math.abs)), yMax = Math.max(1, ...ys.map(Math.abs));
  const sx = v => PAD + ((v + xMax) / (2 * xMax)) * (W - 2 * PAD);
  const sy = v => H - PAD - ((v + yMax) / (2 * yMax)) * (H - 2 * PAD);
  const shown = hidden ? modes.filter(m => !hidden.has(m.category || "mixed")) : modes;
  const points = shown.map(m => {
    const st = statusOf(m.real);
    const sel = m.mode === selectedMode;
    const cat = m.category || "mixed";
    // The selected mode is drawn last so it is never hidden under a neighbour.
    return `<path d="${modeMarkerPath(cat, sel ? 6.2 : 4.2)}" transform="translate(${sx(symlog(m.real)).toFixed(1)},${sy(symlog(m.imag / (2 * Math.PI))).toFixed(1)})" class="pt ${st.cls}${sel ? " sel" : ""}${cat === "reference" ? " hollow" : ""}" data-mode="${m.mode}" data-cat="${cat}"></path>`;
  }).sort((a, b) => (a.includes(" sel") ? 1 : 0) - (b.includes(" sel") ? 1 : 0)).join("");
  return `<svg class="eigenmap-svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="Eigenvalue map (scroll to zoom, drag to pan)" style="width:100%;height:auto;display:block">
    ${eigenGridLines(sx, sy, xMax, yMax, W, H, PAD)}
    <line x1="${PAD}" y1="${sy(0)}" x2="${W - PAD}" y2="${sy(0)}" stroke="var(--baseline)" stroke-width="1"/>
    <line x1="${sx(0)}" y1="${PAD}" x2="${sx(0)}" y2="${H - PAD}" stroke="var(--baseline)" stroke-width="1.5"/>
    ${dampingGuideLines(sx, sy, 0.05, "var(--critical)")}
    ${dampingGuideLines(sx, sy, 0.707, "var(--good)")}
    ${points}
    <text x="${W - PAD}" y="${sy(0) - 8}" text-anchor="end" class="eigen-ticklabel" style="font-size:10.5px">Real part (symlog) — stable ← | → unstable</text>
    <text x="${sx(0) + 8}" y="${PAD + 4}" class="eigen-ticklabel" style="font-size:10.5px">Frequency, Hz (symlog)</text>
  </svg>`;
}

function heatmapTable(matrix, rowLabels, colLabels) {
  const vmax = Math.max(...matrix.flat().map(Math.abs), 1e-12);
  const header = `<tr><th></th>${colLabels.map(c => `<th>${esc(c)}</th>`).join("")}</tr>`;
  const rows = matrix.map((row, i) => `<tr><th style="text-align:right">${esc(rowLabels[i])}</th>${row.map((v, j) => {
    const t = Math.abs(v) / vmax, bg = sequentialBlue(t);
    return `<td style="background:${bg};color:${t > 0.55 ? "#fff" : "var(--text-primary)"}" title="${esc(rowLabels[i])} / ${esc(colLabels[j])}: ${fmt(v, 4)}">${fmt(v, 2)}</td>`;
  }).join("")}</tr>`).join("");
  return `<div class="tablewrap" style="max-height:640px"><table class="heat-table"><thead>${header}</thead><tbody>${rows}</tbody></table></div>`;
}

function barChartHtml(labels, values, color = "var(--series-1)") {
  const barH = 15, gap = 4, PADL = 200, PADR = 70, PADT = 8, W = 820;
  const H = PADT + labels.length * (barH + gap) + 4;
  const vmax = Math.max(...values.map(Math.abs), 1e-12);
  return `<svg viewBox="0 0 ${W} ${H}" style="width:100%;height:auto;display:block">${labels.map((lab, i) => {
    const y = PADT + i * (barH + gap), w = (Math.abs(values[i]) / vmax) * (W - PADL - PADR);
    return `<text x="${PADL - 8}" y="${y + barH * 0.72}" text-anchor="end" class="net-label">${esc(lab)}</text>
      <rect x="${PADL}" y="${y}" width="${w.toFixed(1)}" height="${barH}" fill="${color}" rx="2"><title>${esc(lab)}: ${fmt(values[i], 4)}</title></rect>
      <text x="${PADL + w + 5}" y="${y + barH * 0.72}" class="net-label">${fmt(values[i], 3)}</text>`;
  }).join("")}</svg>`;
}

function polarChartHtml(states, anglesDeg) {
  const W = 460, H = 460, cx = W / 2, cy = H / 2, R = W / 2 - 60;
  let grid = `<circle cx="${cx}" cy="${cy}" r="${R}" fill="none" stroke="var(--gridline)"/><circle cx="${cx}" cy="${cy}" r="${R / 2}" fill="none" stroke="var(--gridline)"/>`;
  for (let a = 0; a < 360; a += 30) {
    const r = a * Math.PI / 180, x = cx + R * Math.cos(r), y = cy - R * Math.sin(r);
    grid += `<line x1="${cx}" y1="${cy}" x2="${x.toFixed(1)}" y2="${y.toFixed(1)}" stroke="var(--gridline)"/>
      <text x="${(cx + (R + 16) * Math.cos(r)).toFixed(1)}" y="${(cy - (R + 16) * Math.sin(r) + 3).toFixed(1)}" text-anchor="middle" class="eigen-ticklabel">${a}°</text>`;
  }
  const arrows = states.map((s, i) => {
    const r = anglesDeg[i] * Math.PI / 180, x = cx + R * Math.cos(r), y = cy - R * Math.sin(r);
    return `<line x1="${cx}" y1="${cy}" x2="${x.toFixed(1)}" y2="${y.toFixed(1)}" stroke="${seriesColor(i)}" stroke-width="3" stroke-linecap="round"><title>${esc(s)}: ${fmt(anglesDeg[i], 1)}°</title></line>
      <circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="4" fill="${seriesColor(i)}"/>`;
  }).join("");
  return `<svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:${W}px;height:auto;display:block;margin:0 auto">${grid}${arrows}</svg>`;
}
