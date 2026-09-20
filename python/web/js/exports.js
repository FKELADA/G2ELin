// Exporting what is on screen: any figure as PNG or SVG, any table as CSV,
// and any chart's own data as CSV.
//
// Nothing on the pages has to opt in. Every figure in this app is an <svg>
// (line charts, the eigenvalue map, root loci, mode shapes, participation
// bars, the network diagram) and every result table is a <table>, so a
// MutationObserver decorates them as they appear: a small toolbar in the
// corner, shown on hover. Charts carry their own series (lineChart puts them
// on the element), which is what "CSV" exports there -- the numbers behind
// the picture, not the picture.
//
// A figure gets its looks from the page's stylesheet, and most of those rules
// are descendant selectors (".chart svg text", ".eigenmap-svg .pt"): none of
// them match once the <svg> stands on its own, which would export a black
// rectangle. So the styles are read off the live elements with
// getComputedStyle and written onto the copy, element by element -- the
// result needs no stylesheet at all. The PNG is that same SVG drawn onto a
// canvas.

const EXPORT_MIN_W = 200;   // smaller <svg>s are icons and legend swatches
const EXPORT_MIN_H = 90;

// What a standalone SVG needs carried over from the page's rules.
const SVG_STYLE_PROPS = [
  "fill", "fill-opacity", "fill-rule", "stroke", "stroke-width", "stroke-opacity", "stroke-dasharray",
  "stroke-linecap", "stroke-linejoin", "opacity", "display", "visibility", "font-family", "font-size",
  "font-weight", "font-style", "letter-spacing", "text-anchor", "dominant-baseline", "text-transform",
];

function inlineComputedStyles(src, clone) {
  const from = [src, ...src.querySelectorAll("*")];
  const to = [clone, ...clone.querySelectorAll("*")];
  for (let i = 0; i < from.length && i < to.length; i++) {
    const cs = getComputedStyle(from[i]);
    let css = "";
    for (const prop of SVG_STYLE_PROPS) {
      const v = cs.getPropertyValue(prop);
      if (v) css += `${prop}:${prop === "font-family" ? `${v}, sans-serif` : v};`;
    }
    to[i].setAttribute("style", css);
  }
}

function downloadBlob(filename, blob) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 0);
}

function exportSlug(name) {
  const s = (name || "figure").toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-+|-+$/g, "").slice(0, 70);
  return s || "figure";
}

// The name a downloaded file gets: the figure's own card title, prefixed by
// the network, so a folder of exports says what it is.
function exportName(node) {
  const card = node.closest(".card, .channel, aside");
  const title = card?.querySelector(".card-title, .chart-title, h3, h4")?.textContent
    || node.closest("[data-fig-name]")?.dataset.figName
    || document.querySelector(".page.active h1")?.textContent
    || "figure";
  const net = state.networkLabel ? `${state.networkLabel}-` : "";
  return exportSlug(`${net}${title.replace(/\s+/g, " ").trim()}`);
}

function serializeSvg(svg) {
  const clone = svg.cloneNode(true);
  const vb = svg.viewBox?.baseVal;
  const rect = svg.getBoundingClientRect();
  const [x, y, w, h] = vb && vb.width
    ? [vb.x, vb.y, vb.width, vb.height]
    : [0, 0, rect.width || 800, rect.height || 400];
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  clone.setAttribute("width", w);
  clone.setAttribute("height", h);
  clone.setAttribute("viewBox", `${x} ${y} ${w} ${h}`);
  inlineComputedStyles(svg, clone);
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", x); bg.setAttribute("y", y);
  bg.setAttribute("width", w); bg.setAttribute("height", h);
  bg.setAttribute("fill", "#ffffff");
  clone.insertBefore(bg, clone.firstChild);
  return new XMLSerializer().serializeToString(clone);
}

function exportSvgFile(svg, name) {
  const text = serializeSvg(svg);
  downloadBlob(`${name}.svg`, new Blob([text], { type: "image/svg+xml;charset=utf-8" }));
}

async function exportPngFile(svg, name, scale = 2) {
  const text = serializeSvg(svg);
  const url = URL.createObjectURL(new Blob([text], { type: "image/svg+xml;charset=utf-8" }));
  try {
    const img = new Image();
    await new Promise((res, rej) => { img.onload = res; img.onerror = () => rej(new Error("render failed")); img.src = url; });
    const vb = svg.viewBox?.baseVal, rect = svg.getBoundingClientRect();
    const w = (vb && vb.width) || rect.width || 800, h = (vb && vb.height) || rect.height || 400;
    const canvas = document.createElement("canvas");
    canvas.width = Math.round(w * scale);
    canvas.height = Math.round(h * scale);
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
    await new Promise(res => canvas.toBlob(b => { downloadBlob(`${name}.png`, b); res(); }, "image/png"));
  } finally {
    URL.revokeObjectURL(url);
  }
}

// --- CSV ---------------------------------------------------------------------
function csvCell(v) {
  const s = v === null || v === undefined ? "" : String(v);
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
function csvText(rows) {
  return rows.map(r => r.map(csvCell).join(",")).join("\r\n");
}
function exportCsv(rows, name) {
  // The BOM keeps Excel happy with the µ, ° and λ these tables are full of.
  downloadBlob(`${name}.csv`, new Blob(["﻿" + csvText(rows)], { type: "text/csv;charset=utf-8" }));
}

function tableRows(table) {
  return [...table.querySelectorAll("tr")].map(tr =>
    [...tr.children].map(td => {
      const input = td.querySelector("input, select");
      if (input) return input.type === "checkbox" ? (input.checked ? "closed" : "open") : input.value;
      return td.textContent.replace(/\s+/g, " ").trim();
    }));
}

// A chart's own numbers: one time column when every series shares it (the
// usual case), otherwise a time column per series.
function seriesRows(series) {
  const live = series.filter(s => s.t && s.t.length);
  if (!live.length) return [["(no data)"]];
  const t0 = live[0].t;
  const shared = live.every(s => s.t === t0 || (s.t.length === t0.length && s.t[0] === t0[0] && s.t[s.t.length - 1] === t0[t0.length - 1]));
  const label = s => s.name + (s.dash ? " (linearised)" : "");
  if (shared) {
    const rows = [["t", ...live.map(label)]];
    for (let i = 0; i < t0.length; i++) rows.push([t0[i], ...live.map(s => s.y[i])]);
    return rows;
  }
  const rows = [live.flatMap(s => [`t (${label(s)})`, label(s)])];
  const n = Math.max(...live.map(s => s.t.length));
  for (let i = 0; i < n; i++) rows.push(live.flatMap(s => (i < s.t.length ? [s.t[i], s.y[i]] : ["", ""])));
  return rows;
}

// --- Decoration ----------------------------------------------------------------
function exportToolbar(host, buttons) {
  if (getComputedStyle(host).position === "static") host.style.position = "relative";
  const bar = el(`<div class="fig-tools"></div>`);
  buttons.forEach(([label, title, fn]) => {
    const b = el(`<button type="button" title="${esc(title)}">${esc(label)}</button>`);
    b.addEventListener("click", async ev => {
      ev.stopPropagation(); ev.preventDefault();
      b.disabled = true;
      try { await fn(); } catch (e) { alert(`Export failed: ${e.message}`); } finally { b.disabled = false; }
    });
    bar.appendChild(b);
  });
  host.appendChild(bar);
  return bar;
}

function decorateFigures(root = document) {
  root.querySelectorAll("svg:not([data-exported])").forEach(svg => {
    if (svg.closest(".legend, .fig-tools, .brand, .net-zoom, .picker-pop, #tooltip, #tour-cursor, #tour-caption")) return;
    const box = svg.getBoundingClientRect();
    const vb = svg.viewBox?.baseVal;
    const w = box.width || vb?.width || 0, h = box.height || vb?.height || 0;
    if (w < EXPORT_MIN_W || h < EXPORT_MIN_H) return;
    svg.setAttribute("data-exported", "1");
    const chart = svg.closest(".chart");
    const host = chart || svg.parentElement;
    if (!host || host.querySelector(":scope > .fig-tools")) return;
    const name = () => exportName(svg);
    const buttons = [
      ["PNG", "Download this figure as a PNG image", () => exportPngFile(svg, name())],
      ["SVG", "Download this figure as a vector SVG", () => exportSvgFile(svg, name())],
    ];
    if (chart && chart.__chart) {
      buttons.push(["CSV", "Download the data behind this chart", () => exportCsv(seriesRows(chart.__chart.series), name())]);
    }
    exportToolbar(host, buttons);
  });

  root.querySelectorAll("table:not([data-exported])").forEach(table => {
    if (table.closest(".fig-tools, .picker-pop, #tooltip, #tour-caption")) return;
    if (table.rows.length < 2) return;
    table.setAttribute("data-exported", "1");
    const wrap = table.closest(".tablewrap");
    const host = (wrap && wrap.parentElement) || table.parentElement;
    if (!host || host.querySelector(":scope > .fig-tools")) return;
    exportToolbar(host, [["CSV", "Download this table as CSV", () => exportCsv(tableRows(table), exportName(table))]]);
  });
}

// Pages redraw constantly (a scrubber, a live EMT trace); decorating on a
// debounced MutationObserver keeps that out of every page's own code.
function watchForFigures() {
  let timer = null;
  const run = () => { timer = null; decorateFigures(); };
  new MutationObserver(() => { if (!timer) timer = setTimeout(run, 120); }).observe(document.body, { childList: true, subtree: true });
  decorateFigures();
}
