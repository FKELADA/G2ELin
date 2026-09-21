// My networks: save the network being edited under a name and get it back
// later, after the tab or the browser has been closed.
//
// Where it lives: this browser's localStorage. That is per user, needs no
// account or server, behaves the same on a laptop, on a colleague's machine and
// on the hosted app (where a server-side folder would be shared by everyone and
// erased by every redeploy). Its limits are the flip side: it belongs to one
// browser profile at one address (host:port), and clearing site data removes it
// -- hence "Export all", a JSON backup that can be imported anywhere.
//
// Two things are stored:
//   - the saved networks, by explicit "Save" (SavedNetworks.KEY);
//   - a draft of the work in progress, written automatically a moment after
//     every edit, so closing the tab before saving loses nothing. On the next
//     visit the page offers to restore it (SavedNetworks.DRAFT_KEY).
//
// A saved record: {id, name, savedAt, origin: {presetId, label}|null, network,
// positions}. `network` is the Network JSON the whole app already uses, and
// `positions` the diagram layout, so a network comes back exactly as it was
// left. The file "Export" writes is the same shape the page's own JSON export
// uses, and can be imported through "Import JSON".

const SavedNetworks = {
  KEY: "g2elin.saved.v1",
  DRAFT_KEY: "g2elin.draft.v1",
  BACKUP_FORMAT: "g2elin-saved-networks",
  mounted: false,
  _timer: null,

  // --- storage ---------------------------------------------------------------
  available() {
    try {
      localStorage.setItem("g2elin.probe", "1");
      localStorage.removeItem("g2elin.probe");
      return true;
    } catch { return false; }   // blocked, or a private window
  },
  items() {
    try {
      const v = JSON.parse(localStorage.getItem(this.KEY));
      const list = Array.isArray(v?.items) ? v.items : [];
      return list.filter(i => i && i.id && i.network && Array.isArray(i.network.buses))
        .sort((a, b) => String(b.savedAt).localeCompare(String(a.savedAt)));
    } catch { return []; }
  },
  // Returns null on success, else a sentence saying why not.
  persist(items) {
    try {
      localStorage.setItem(this.KEY, JSON.stringify({ version: 1, items }));
      return null;
    } catch (e) {
      return e && e.name === "QuotaExceededError"
        ? "The browser's storage for this site is full. Delete some saved networks, or export them first."
        : "The browser refused to store it (a private window, or storage blocked for this site). Use Export instead.";
    }
  },
  draft() {
    try {
      const d = JSON.parse(localStorage.getItem(this.DRAFT_KEY));
      return d && d.network && Array.isArray(d.network.buses) ? d : null;
    } catch { return null; }
  },
  clearDraft() { try { localStorage.removeItem(this.DRAFT_KEY); } catch { /* nothing to clear */ } },

  // --- what the current network is -------------------------------------------------
  clone: v => JSON.parse(JSON.stringify(v)),
  newId: () => `n${Date.now().toString(36)}${Math.random().toString(36).slice(2, 6)}`,
  describe(net) {
    const n = (k) => (net[k] || []).length;
    return `${n("buses")} buses · ${n("lines")} lines · ${n("der_units")} unit${n("der_units") === 1 ? "" : "s"}`;
  },
  when: iso => { const d = new Date(iso); return isNaN(d) ? "" : d.toLocaleString(undefined, { dateStyle: "medium", timeStyle: "short" }); },
  current() { return state.savedId ? this.items().find(i => i.id === state.savedId) || null : null; },
  // Edits that would be lost by replacing the network being edited.
  hasUnsavedWork() {
    const n = state.network;
    if (!n) return false;
    if (state.savedId) return JSON.stringify(n) !== state.savedBaseline;
    if (state.presetId) return isModified();
    return n.buses.length > 0;    // built from scratch, or imported
  },
  originOfCurrent() {
    const cur = this.current();
    if (cur) return cur.origin || null;
    if (state.presetId) {
      const p = state.presets.find(x => x.id === state.presetId);
      return { presetId: state.presetId, label: p ? p.name : state.presetId };
    }
    return null;
  },

  // --- actions ---------------------------------------------------------------------
  save(asNew = false) {
    const name = $("#sv-name").value.trim();
    if (!state.network) return;
    if (!name) { this.flash("Give it a name first.", "bad"); $("#sv-name").focus(); return; }
    let items = this.items();
    let target = !asNew && state.savedId ? items.find(i => i.id === state.savedId) || null : null;
    const clash = items.find(i => i.name.toLowerCase() === name.toLowerCase() && i !== target);
    if (clash) {
      if (!confirm(`A saved network called “${clash.name}” already exists (saved ${this.when(clash.savedAt)}).\n\nReplace it?`)) return;
      items = items.filter(i => i !== clash);
    }
    const rec = {
      id: target ? target.id : this.newId(),
      name, savedAt: new Date().toISOString(),
      origin: this.originOfCurrent(),
      network: this.clone(state.network),
      positions: this.clone(state.positions || {}),
    };
    items = items.filter(i => i.id !== rec.id);
    items.push(rec);
    const err = this.persist(items);
    if (err) { this.flash(err, "bad"); return; }
    state.savedId = rec.id;
    state.savedBaseline = JSON.stringify(state.network);
    state.networkLabel = name;
    this.clearDraft();          // what the draft held is now saved
    NetworkPage.renderSummary();
    this.render();
    this.flash(`Saved “${name}”.`, "ok");
  },

  load(id) {
    const item = this.items().find(i => i.id === id);
    if (!item) return;
    if (this.hasUnsavedWork() && !confirm(`Loading “${item.name}” replaces the network you are editing, which has unsaved changes.\n\nLoad it anyway? (Your unsaved work is kept as a draft you can restore.)`)) return;
    const positions = item.positions && Object.keys(item.positions).length ? this.clone(item.positions) : null;
    setNetwork(this.clone(item.network), { label: item.name, presetId: null, savedId: item.id, positions });
    if ($("#net-preset")) $("#net-preset").value = "";
    $("#sv-name").value = item.name;
    this.flash(`Loaded “${item.name}”.`, "ok");
  },

  rename(id) {
    const items = this.items(), item = items.find(i => i.id === id);
    if (!item) return;
    const name = (prompt("New name:", item.name) || "").trim();
    if (!name || name === item.name) return;
    if (items.some(i => i.id !== id && i.name.toLowerCase() === name.toLowerCase())) { alert(`There is already a saved network called “${name}”.`); return; }
    item.name = name;
    const err = this.persist(items);
    if (err) { alert(err); return; }
    if (state.savedId === id) { state.networkLabel = name; $("#sv-name").value = name; NetworkPage.renderSummary(); }
    this.render();
  },

  remove(id) {
    const items = this.items(), item = items.find(i => i.id === id);
    if (!item || !confirm(`Delete “${item.name}”?\n\nThis cannot be undone (export it first if in doubt).`)) return;
    const err = this.persist(items.filter(i => i.id !== id));
    if (err) { alert(err); return; }
    if (state.savedId === id) { state.savedId = null; state.savedBaseline = null; }   // the network stays open, now unsaved
    this.render();
  },

  // --- files ---------------------------------------------------------------------------
  download(filename, obj) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([JSON.stringify(obj, null, 2)], { type: "application/json" }));
    a.download = filename;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  },
  fileSlug: s => (s || "network").replace(/[^\w.-]+/g, "_"),
  exportOne(id) {
    const item = this.items().find(i => i.id === id);
    if (item) this.download(`${this.fileSlug(item.name)}.json`, { network: item.network, positions: item.positions || {} });
  },
  exportAll() {
    const items = this.items();
    if (!items.length) { this.flash("Nothing saved yet.", "bad"); return; }
    this.download(`g2elin-saved-networks-${new Date().toISOString().slice(0, 10)}.json`,
      { format: this.BACKUP_FORMAT, version: 1, exportedAt: new Date().toISOString(), items });
  },
  async importBackup(file) {
    if (!file) return;
    try {
      const data = JSON.parse(await file.text());
      const incoming = data.format === this.BACKUP_FORMAT && Array.isArray(data.items)
        ? data.items
        : (data.network || Array.isArray(data.buses))
          ? [{ id: this.newId(), name: (data.network || data).name || file.name.replace(/\.json$/i, ""), savedAt: new Date().toISOString(), origin: null, network: data.network || data, positions: data.positions || {} }]
          : null;
      if (!incoming) throw new Error("not a G2ELin backup or network file");
      const items = this.items(), have = new Set(items.map(i => i.id)), names = new Set(items.map(i => i.name.toLowerCase()));
      let added = 0, skipped = 0;
      for (const it of incoming) {
        if (!it || !it.network || !Array.isArray(it.network.buses)) { skipped++; continue; }
        if (have.has(it.id)) { skipped++; continue; }          // already here
        ["lines", "transformers", "loads", "der_units"].forEach(k => { it.network[k] ||= []; });
        let name = String(it.name || "Imported network");
        while (names.has(name.toLowerCase())) name += " (imported)";
        names.add(name.toLowerCase());
        items.push({ ...it, id: it.id || this.newId(), name, savedAt: it.savedAt || new Date().toISOString() });
        added++;
      }
      const err = this.persist(items);
      if (err) throw new Error(err);
      this.render();
      this.flash(`Imported ${added} network${added === 1 ? "" : "s"}${skipped ? ` (${skipped} skipped: already present or invalid)` : ""}.`, "ok");
    } catch (e) { this.flash(`Could not import: ${e.message}`, "bad"); }
    $("#sv-import-file").value = "";
  },

  // --- the draft ---------------------------------------------------------------------------
  autosaveSoon() {
    clearTimeout(this._timer);
    this._timer = setTimeout(() => this.autosave(), 700);
  },
  autosave() {
    if (!this.hasUnsavedWork()) return;
    try {
      localStorage.setItem(this.DRAFT_KEY, JSON.stringify({
        at: new Date().toISOString(), label: state.networkLabel, savedId: state.savedId,
        origin: this.originOfCurrent(), network: state.network, positions: state.positions || {},
      }));
    } catch { /* storage full or blocked: the explicit Save reports it */ }
    this.renderDraft();
  },
  restoreDraft() {
    const d = this.draft();
    if (!d) return;
    if (this.hasUnsavedWork() && !confirm("Restoring replaces the network you are editing, which has unsaved changes. Continue?")) return;
    const existing = d.savedId ? this.items().find(i => i.id === d.savedId) : null;
    const positions = d.positions && Object.keys(d.positions).length ? this.clone(d.positions) : null;
    setNetwork(this.clone(d.network), { label: d.label || "Recovered network", presetId: null, savedId: existing ? existing.id : null, positions });
    // Measured against what was last *saved*, so the edits still count as unsaved.
    if (existing) state.savedBaseline = JSON.stringify(existing.network);
    if ($("#net-preset")) $("#net-preset").value = "";
    $("#sv-name").value = existing ? existing.name : (d.label || "");
    // The draft stays until saved or discarded: restoring and closing again must not lose it.
    this.render();
  },
  discardDraft() { this.clearDraft(); this.renderDraft(); },

  // --- rendering ----------------------------------------------------------------------------
  flash(text, kind = "ok") {
    const box = $("#sv-status");
    if (!box) return;
    box.innerHTML = `<span class="${kind === "bad" ? "bad" : "ok"}">${esc(text)}</span>`;
    clearTimeout(this._flash);
    this._flash = setTimeout(() => this.renderStatus(), 4500);
  },
  renderStatus() {
    const box = $("#sv-status");
    if (!box) return;
    const cur = this.current();
    let html;
    if (!state.network) html = "";
    else if (cur && !this.hasUnsavedWork()) html = `<span class="ok">Saved as “${esc(cur.name)}”</span> · ${esc(this.when(cur.savedAt))}`;
    else if (cur) html = `<span class="warn">Unsaved changes</span> to “${esc(cur.name)}”`;
    else if (state.presetId && !isModified()) html = `<span class="muted">An unmodified preset — save it to keep your own copy.</span>`;
    else html = `<span class="warn">Not saved yet</span>`;
    box.innerHTML = html;
    $("#sv-save").textContent = cur ? "Save changes" : "Save";
    $("#sv-save-new").hidden = !cur;
    $("#sv-save").disabled = $("#sv-save-new").disabled = !state.network || this.usable === false;
  },
  renderDraft() {
    const box = $("#sv-draft");
    if (!box) return;
    const d = this.draft();
    // Nothing to offer when the draft is what is already open.
    if (!d || (state.network && JSON.stringify(d.network) === JSON.stringify(state.network))) { box.innerHTML = ""; return; }
    box.innerHTML = `<div class="notice warn-bg" style="display:flex;gap:0.8rem;align-items:center;flex-wrap:wrap;margin-bottom:0.9rem">
      <span style="flex:1;min-width:16rem"><b>Unsaved work was found</b> from ${esc(this.when(d.at))} — “${esc(d.label || "network")}” (${esc(this.describe(d.network))}).</span>
      <button class="small" data-sv="restore">Restore it</button><button class="ghost small" data-sv="discard">Discard</button></div>`;
  },
  renderList() {
    const box = $("#sv-list"), items = this.items();
    $("#sv-count").textContent = items.length ? `— ${items.length} saved` : "";
    if (!items.length) { box.innerHTML = `<p class="empty">Nothing saved yet. Edit a preset or build a network, give it a name above and press <b>Save</b>.</p>`; return; }
    box.innerHTML = `<div class="tablewrap" style="max-height:320px"><table><thead><tr><th>Name</th><th>Based on</th><th>Contents</th><th>Saved</th><th></th></tr></thead><tbody>${items.map(i => `
      <tr class="${i.id === state.savedId ? "hl" : ""}"><td class="name"><b>${esc(i.name)}</b>${i.id === state.savedId ? ' <span class="badge accent">open</span>' : ""}</td>
        <td class="name">${esc(i.origin && i.origin.label ? i.origin.label : "built from scratch")}</td>
        <td class="name">${esc(this.describe(i.network))}</td><td class="name">${esc(this.when(i.savedAt))}</td>
        <td style="white-space:nowrap"><button class="secondary small" data-sv="load" data-id="${esc(i.id)}">Load</button>
          <button class="ghost small" data-sv="rename" data-id="${esc(i.id)}">Rename</button>
          <button class="ghost small" data-sv="export" data-id="${esc(i.id)}" title="Download this network as a JSON file">Export</button>
          <button class="ghost small" data-sv="delete" data-id="${esc(i.id)}">Delete</button></td></tr>`).join("")}</tbody></table></div>`;
  },
  render() { this.renderStatus(); this.renderDraft(); this.renderList(); },

  // --- wiring -----------------------------------------------------------------------------------
  mount() {
    if (this.mounted) return;
    this.mounted = true;
    const root = $("#net-saved");
    const ok = this.usable = this.available();
    root.innerHTML = `
      <div class="card-title">My networks <span class="card-sub" id="sv-count"></span></div>
      <div id="sv-draft"></div>
      ${ok ? "" : `<div class="notice warn-bg" style="margin-bottom:0.9rem">Saving isn't available here — the browser is blocking local storage (a private window does this). Use <b>Export JSON</b> above to keep a network as a file.</div>`}
      <div class="controls" style="align-items:flex-end">
        <div class="field" style="min-width:min(24rem,100%);flex:1"><label for="sv-name">Save the network you are editing as</label>
          <input id="sv-name" type="text" maxlength="80" placeholder="e.g. WSCC with a grid-forming unit at bus 2" autocomplete="off" style="width:100%"${ok ? "" : " disabled"}></div>
        <button id="sv-save"${ok ? "" : " disabled"}>Save</button>
        <button class="secondary" id="sv-save-new" hidden title="Keep the saved version and store this one as a separate network">Save as new copy</button>
      </div>
      <p class="status-line" id="sv-status" style="margin:0.5rem 0 0.9rem"></p>
      <div id="sv-list"></div>
      <div class="controls" style="margin-top:0.8rem;align-items:center">
        <button class="ghost small" id="sv-export-all" title="Download every saved network as one JSON file">Export all (backup)</button>
        <button class="ghost small" id="sv-import" title="Add networks from a backup file, or a single network JSON">Import a backup</button>
        <input type="file" id="sv-import-file" accept="application/json,.json" hidden>
        <span class="muted" style="font-size:0.76rem">Stored in this browser only — export a backup to keep it safe or to open it on another computer.</span>
      </div>`;
    $("#sv-save").addEventListener("click", () => this.save(false));
    $("#sv-save-new").addEventListener("click", () => this.save(true));
    $("#sv-name").addEventListener("keydown", e => { if (e.key === "Enter") { e.preventDefault(); this.save(false); } });
    $("#sv-export-all").addEventListener("click", () => this.exportAll());
    $("#sv-import").addEventListener("click", () => $("#sv-import-file").click());
    $("#sv-import-file").addEventListener("change", e => this.importBackup(e.target.files[0]));
    root.addEventListener("click", e => {
      const b = e.target.closest("[data-sv]");
      if (!b) return;
      const id = b.dataset.id;
      ({ load: () => this.load(id), rename: () => this.rename(id), export: () => this.exportOne(id), delete: () => this.remove(id),
         restore: () => this.restoreDraft(), discard: () => this.discardDraft() })[b.dataset.sv]?.();
    });
    const cur = this.current();
    $("#sv-name").value = cur ? cur.name : (state.networkLabel || "");
    this.render();
  },
};

// Autosave watches every edit from the moment the app loads -- a breaker toggled
// on the Power Flow page counts as much as one on the Network page -- while the
// section itself is only drawn once the Network page is opened (mount()).
on("network:changed", () => { SavedNetworks.autosaveSoon(); if (SavedNetworks.mounted) SavedNetworks.renderStatus(); });
on("network:layout-moved", () => SavedNetworks.autosaveSoon());
// A newly loaded network offers its name as the starting point for a save.
on("network:loaded", () => {
  if (!SavedNetworks.mounted) return;
  const cur = SavedNetworks.current();
  $("#sv-name").value = cur ? cur.name : (state.networkLabel || "");
  SavedNetworks.render();
});
