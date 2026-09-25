// Model order: how much of each element's dynamics the models keep.
//
// Two levels of control, both backed by the same thing. A *level* is a named
// preset ("6th order", "Droop only"); the per-state-group selects underneath
// it are the actual setting, and a level is just a shortcut for filling them
// in. Picking a level fills the groups; changing a group leaves the level
// reading "Custom" -- there is no separate "advanced mode", because there is
// no separate mechanism.
//
// Nothing here knows the names of the levels or the state groups: the whole
// catalogue comes from GET /api/model-levels, so adding a level or a group
// server-side makes it appear here with no frontend change.

const MODE_LABEL = {
  dynamic: "Dynamic",
  algebraic: "Algebraic",
  frozen: "Frozen",
};
const MODE_HELP = {
  dynamic: "Integrated — the state has its own differential equation.",
  algebraic: "d/dt set to zero — the state follows the rest of the system instantly. For dynamics faster than the study.",
  frozen: "Held at its equilibrium value — the state does not move. For dynamics slower than the study.",
};
const MODE_CLASS = { dynamic: "ok", algebraic: "warn", frozen: "muted" };

const ModelOrder = {
  catalogue: null,      // {kind: ElementModelInfo}
  summary: null,
  _pending: null,

  async load() {
    if (this.catalogue) return this.catalogue;
    if (!this._pending) {
      this._pending = api("/api/model-levels").then(r => {
        this.catalogue = Object.fromEntries(r.elements.map(e => [e.kind, e]));
        return this.catalogue;
      });
    }
    return this._pending;
  },

  models() {
    const net = state.network;
    if (!net) return null;
    if (!net.models) net.models = {};
    return net.models;
  },

  // The mode of every group of one element type, level then overrides --
  // the same resolution order the server's ModelOptions.group_modes_for does.
  groupModes(kind, { level, overrides } = {}) {
    const e = this.catalogue && this.catalogue[kind];
    if (!e) return {};
    const m = this.models() || {};
    const lvlId = level !== undefined ? level : (m[`${kind}_level`] || e.default_level);
    const lvl = e.levels.find(l => l.id === lvlId) || e.levels.find(l => l.id === e.default_level);
    const out = { ...(lvl ? lvl.modes : {}) };
    Object.entries(overrides !== undefined ? overrides : (m[`${kind}_states`] || {})).forEach(([g, v]) => { out[g] = v; });
    return out;
  },

  // The state groups one unit really has. A machine's AVR and PSS groups
  // depend on which regulator models it carries, so the catalogue's default
  // groups are swapped for the chosen ones. Everything else — levels,
  // cascades, the mode controls — works off this, so none of it has to know
  // that regulators are selectable at all.
  // The groups a model option brings: one for a machine's regulators, a set
  // for a converter's control law, none at all for a model fitted as "none".
  optionGroups(opt) {
    if (!opt) return [];
    return (opt.groups && opt.groups.length) ? opt.groups : (opt.group ? [opt.group] : []);
  },

  groupsFor(kind, der = null) {
    const e = this.catalogue[kind];
    if (!e || !e.regulators || !e.regulators.length || !der) return e ? e.groups : [];
    // Every group any option of any slot could bring, and the ones the
    // chosen options actually do. What a slot owns is replaced wholesale, so
    // a law with two groups can stand where another had one, or none.
    const owned = new Set(), selected = [];
    e.regulators.forEach(slot => {
      slot.options.forEach(o => this.optionGroups(o).forEach(g => owned.add(g.id)));
      const chosen = slot.options.find(o => o.id === (der[slot.id] || slot.default));
      this.optionGroups(chosen).forEach(g => selected.push(g));
    });
    const byId = Object.fromEntries(selected.map(g => [g.id, g]));
    const out = [];
    let spliced = false;
    e.groups.forEach(g => {
      if (!owned.has(g.id)) { out.push(g); return; }
      // The chosen groups go where the default's were, keeping the panel's
      // reading order: plant, then control law, then angle.
      if (!spliced) { out.push(...selected); spliced = true; }
    });
    if (!spliced && selected.length) out.push(...selected);
    return out.filter((g, i) => out.findIndex(x => x.id === g.id) === i).map(g => byId[g.id] || g);
  },

  setRegulator(kind, slotId, modelId, der) {
    // The group ids are the same whichever model is chosen, so a level or an
    // override stays meaningful across a swap; only the states behind the
    // group change. Any override that the new model disallows is dropped.
    der[slotId] = modelId;
    const allowed = Object.fromEntries(this.groupsFor(kind, der).map(g => [g.id, g.allowed]));
    Object.keys(der.states || {}).forEach(gid => {
      if (allowed[gid] && !allowed[gid].includes(der.states[gid])) delete der.states[gid];
    });
    this.afterEdit();
  },

  // A regulator model's short name, for the summary table.
  regulatorLabel(kind, slotId, modelId) {
    const slot = (this.catalogue[kind]?.regulators || []).find(r => r.id === slotId);
    const opt = slot && slot.options.find(o => o.id === (modelId || slot.default));
    return opt ? opt.label.replace(/\s*\(.*\)\s*$/, "") : (modelId || "");
  },

  // Which named level a set of group modes is, or null for a custom one.
  matchingLevel(kind, modes, der = null) {
    const e = this.catalogue[kind];
    const groups = this.groupsFor(kind, der);
    const lvl = e.levels.find(l => groups.every(g => (l.modes[g.id] || g.default) === modes[g.id]));
    return lvl ? lvl.id : null;
  },

  // --- editing ---------------------------------------------------------------

  setLevel(kind, levelId, der = null) {
    // A level replaces everything, overrides included: picking "6th order"
    // must give the 6th-order model, not the 6th-order model plus whatever
    // was left over from the last custom edit.
    if (der) { der.level = levelId; der.states = {}; }
    else { this.models()[`${kind}_level`] = levelId; this.models()[`${kind}_states`] = {}; }
    this.afterEdit();
  },

  afterEdit() {
    networkChanged();
    this.refreshBlocks();
    this.refreshSummary();
  },

  // A control loop's integrator can only be made algebraic once what it
  // regulates is an unknown too (StateGroupInfo.requires — the server
  // refuses the rest). Rather than let the user pick a combination and be
  // told no, carry the dependency along: making a loop algebraic makes what
  // it needs algebraic, and putting one of those back dynamic puts the
  // loops that depend on it back too. Mutates `modes` in place.
  applyCascade(kind, modes, groupId, mode, der = null) {
    const groups = this.groupsFor(kind, der);
    const by = Object.fromEntries(groups.map(g => [g.id, g]));
    if (mode === "algebraic") {
      const pending = [...(by[groupId].requires || [])];
      while (pending.length) {
        const id = pending.pop();
        if (modes[id] === "algebraic" || !by[id] || by[id].locked) continue;
        modes[id] = "algebraic";
        pending.push(...(by[id].requires || []));
      }
      return;
    }
    // Going back the other way: anything that required this group can no
    // longer be algebraic either.
    let changed = true;
    while (changed) {
      changed = false;
      groups.forEach(g => {
        if (modes[g.id] !== "algebraic") return;
        if ((g.requires || []).some(r => modes[r] !== "algebraic")) {
          modes[g.id] = "dynamic";
          changed = true;
        }
      });
    }
  },

  setGroup(kind, groupId, mode, der = null) {
    const modes = der
      ? { ...this.groupModes(kind, { level: der.level, overrides: der.states || {} }), [groupId]: mode }
      : { ...this.groupModes(kind), [groupId]: mode };
    this.applyCascade(kind, modes, groupId, mode, der);
    // Store the result as a level plus the smallest set of overrides that
    // reproduces it, so the saved network says "6th order" where it can
    // rather than a list of nine groups that happens to mean that.
    const named = this.matchingLevel(kind, modes, der);
    const e = this.catalogue[kind];
    const groups = this.groupsFor(kind, der);
    if (der) {
      if (named) { der.level = named; der.states = {}; }
      else {
        const base = e.levels.find(l => l.id === (der.level || e.default_level)) || e.levels[0];
        der.states = Object.fromEntries(
          groups.filter(g => modes[g.id] !== (base.modes[g.id] || g.default)).map(g => [g.id, modes[g.id]])
        );
      }
    } else if (named) {
      this.models()[`${kind}_level`] = named;
      this.models()[`${kind}_states`] = {};
    } else {
      const base = e.levels.find(l => l.id === (this.models()[`${kind}_level`] || e.default_level)) || e.levels[0];
      this.models()[`${kind}_states`] = Object.fromEntries(
        e.groups.filter(g => modes[g.id] !== (base.modes[g.id] || g.default)).map(g => [g.id, modes[g.id]])
      );
    }
    this.afterEdit();
  },

  // --- rendering -------------------------------------------------------------

  // The level picker plus one row per state group. `der` non-null edits that
  // unit alone; otherwise it edits the network's default for `kind`.
  controlsHtml(kind, der = null) {
    const e = this.catalogue[kind];
    if (!e) return "";
    const modes = der
      ? this.groupModes(kind, { level: der.level, overrides: der.states || {} })
      : this.groupModes(kind);
    const current = this.matchingLevel(kind, modes, der);
    const scope = der ? `der:${der.id}` : "net";
    const groups = this.groupsFor(kind, der);
    const nDyn = groups.filter(g => modes[g.id] === "dynamic").reduce((n, g) => n + g.states.length, 0);
    const nAll = groups.reduce((n, g) => n + g.states.length, 0);
    // Regulator pickers, for one machine at a time: which exciter a unit has
    // is a property of that unit, not a network-wide default.
    const regs = (der && e.regulators) ? e.regulators : [];

    const lvl = e.levels.find(l => l.id === current);
    return `
      <div class="mo-block" data-mo-kind="${kind}" data-mo-scope="${scope}">
        <div class="mo-head">
          <select class="mo-level" data-mo-kind="${kind}" data-mo-scope="${scope}">
            ${e.levels.map(l => `<option value="${l.id}"${l.id === current ? " selected" : ""}>${esc(l.label)}</option>`).join("")}
            ${current === null ? `<option value="" selected>Custom</option>` : ""}
          </select>
          <span class="mo-count">${nDyn} of ${nAll} states</span>
        </div>
        ${lvl && lvl.note ? `<p class="mo-note">${esc(lvl.note)}</p>` : ""}
        ${current === null ? `<p class="mo-note warn">Custom combination — not one of the named models. Run the adequacy check below before trusting it.</p>` : ""}
        ${regs.length ? `<table class="mo-groups mo-regulators">
          ${regs.map(slot => {
            const chosen = der[slot.id] || slot.default;
            const opt = slot.options.find(o => o.id === chosen);
            const states = this.optionGroups(opt).flatMap(g => g.states);
            return `<tr><td><span class="mo-g">${esc(slot.label)}</span><span class="mo-s">${esc(states.length ? states.join(", ") : "no states")}</span></td>
              <td><select class="mo-regulator" data-mo-kind="${kind}" data-mo-scope="${scope}" data-mo-slot="${slot.id}" title="${esc(this.optionGroups(opt).map(g => g.note).filter(Boolean).join(" "))}">
                ${slot.options.map(o => `<option value="${o.id}"${o.id === chosen ? " selected" : ""}>${esc(o.label)}</option>`).join("")}
              </select></td></tr>`;
          }).join("")}
        </table>` : ""}
        <table class="mo-groups">
          ${groups.map(g => {
            const mode = modes[g.id];
            if (g.locked) {
              return `<tr class="locked"><td><span class="mo-g">${esc(g.label)}</span><span class="mo-s">${esc(g.states.join(", "))}</span></td>
                <td><span class="badge">always dynamic</span></td></tr>`;
            }
            const needs = (g.requires || []).map(r => groups.find(x => x.id === r)?.label || r);
            const title = [g.note, needs.length ? `Making this algebraic also makes ${needs.join(" and ")} algebraic.` : ""]
              .filter(Boolean).join(" — ");
            return `<tr><td><span class="mo-g">${esc(g.label)}</span><span class="mo-s">${esc(g.states.join(", "))}</span></td>
              <td><select class="mo-group" data-mo-kind="${kind}" data-mo-scope="${scope}" data-mo-group="${g.id}" title="${esc(title)}">
                ${g.allowed.map(m => `<option value="${m}"${m === mode ? " selected" : ""}>${MODE_LABEL[m]}</option>`).join("")}
              </select></td></tr>`;
          }).join("")}
        </table>
        <p class="mo-legend">${Object.entries(MODE_HELP).map(([m, h]) => `<span><b class="${MODE_CLASS[m]}">${MODE_LABEL[m]}</b> ${esc(h)}</span>`).join("")}</p>
      </div>`;
  },

  // Bind the selects inside `root` (idempotent — safe to call after each re-render).
  bind(root) {
    $$(".mo-level", root).forEach(sel => sel.addEventListener("change", () => {
      if (!sel.value) return;  // "Custom" is a readout, not a choice
      this.setLevel(sel.dataset.moKind, sel.value, this.scopeDer(sel.dataset.moScope));
    }));
    $$(".mo-group", root).forEach(sel => sel.addEventListener("change", () =>
      this.setGroup(sel.dataset.moKind, sel.dataset.moGroup, sel.value, this.scopeDer(sel.dataset.moScope))));
    $$(".mo-regulator", root).forEach(sel => sel.addEventListener("change", () => {
      const der = this.scopeDer(sel.dataset.moScope);
      if (der) this.setRegulator(sel.dataset.moKind, sel.dataset.moSlot, sel.value, der);
    }));
  },

  // Re-render every block in place after an edit: picking a level changes
  // the group selects under it, and changing a group changes the level
  // above it to "Custom", so neither can be left showing what it showed
  // before the edit. Blocks only — the <details> around them keep their
  // open state, and the panel keeps its scroll position.
  refreshBlocks() {
    $$(".mo-block").forEach(block => {
      const kind = block.dataset.moKind;
      const der = this.scopeDer(block.dataset.moScope);
      if (der === null && block.dataset.moScope !== "net") return;  // unit has gone
      const fresh = el(this.controlsHtml(kind, der));
      block.replaceWith(fresh);
      this.bind(fresh);  // the new block only -- binding its parent would
                         // double up the listeners on its siblings
    });
  },

  scopeDer(scope) {
    if (!scope || scope === "net") return null;
    const id = +scope.split(":")[1];
    return state.network.der_units.find(d => d.id === id) || null;
  },

  // --- the network-wide panel ------------------------------------------------

  panelHtml() {
    if (!this.catalogue) return `<p class="empty">Loading…</p>`;
    if (!state.network) return `<p class="empty">No network loaded.</p>`;
    const kinds = ["network", "sm", "gfm", "gfl"];
    const present = new Set((state.network.der_units || []).map(d => d.unit_type));
    const m = this.models();
    return `
      <div id="mo-summary" class="mo-summary"></div>
      ${kinds.filter(k => k === "network" || present.has(k)).map(k => `
        <details class="mo-section"${k === "network" ? " open" : ""}>
          <summary>${esc(this.catalogue[k].label)}${k === "network" ? "" : " — default for every unit of this type"}</summary>
          ${this.controlsHtml(k)}
          ${k === "network" ? `
            <div class="field" style="margin-top:0.6rem">
              <label for="mo-net-freq">Speed terms in the passive elements</label>
              <select id="mo-net-freq">
                <option value="frame"${(m.network_frequency || "frame") === "frame" ? " selected" : ""}>Follow the reference frame (EMT convention)</option>
                <option value="nominal"${m.network_frequency === "nominal" ? " selected" : ""}>Pinned to nominal (phasor convention)</option>
              </select>
              <span class="hint">Whether the ω·L and ω·C terms in lines, loads and bus capacitances use the frame's own speed or the nominal one. Phasor tools pin them; it matters for frequency studies.</span>
            </div>` : ""}
        </details>`).join("")}
      <div class="mo-adequacy">
        <div class="controls">
          <div class="field" style="max-width:190px"><label for="mo-band">Band to reproduce <span class="unit">(Hz)</span></label>
            <input id="mo-band" type="number" step="any" value="5" min="0.1"></div>
          <button class="secondary" id="mo-check">Check this model against the full one</button>
        </div>
        <p class="mo-note">Electromechanical studies want a few Hz; converter interactions want tens. The check compares this model against the same network at full order, over that band.</p>
        <div id="mo-adequacy-out"></div>
      </div>`;
  },

  // Which unit types the network has, which is what the panel's sections
  // are: one per type present, so adding the first converter or retyping
  // the last machine changes the panel itself and not just its numbers.
  presentKey() {
    return [...new Set((state.network?.der_units || []).map(d => d.unit_type))].sort().join(",");
  },

  async mount(box) {
    await this.load();
    this.box = box;
    this.rendered = this.presentKey();
    box.innerHTML = this.panelHtml();
    this.bind(box);
    const freq = $("#mo-net-freq", box);
    if (freq) freq.addEventListener("change", () => { this.models().network_frequency = freq.value; networkChanged(); });
    const btn = $("#mo-check", box);
    if (btn) btn.addEventListener("click", () => this.runAdequacy());
    this.refreshSummary();
  },

  async refreshSummary() {
    // A unit added, removed or retyped changes which sections belong here,
    // and the summary alone cannot express that -- rebuild the panel.
    if (this.box && this.rendered !== undefined && this.presentKey() !== this.rendered) {
      await this.mount(this.box);
      return;
    }
    const out = $("#mo-summary");
    if (!out || !state.network) return;
    const version = state.version;
    try {
      const r = await netPost("model");
      if (version !== state.version) return;
      this.summary = r;
      const cls = { EMT: "ok", RMS: "warn", Mixed: "muted" }[r.model_class] || "";
      const saved = r.n_states_full - r.n_states;
      out.innerHTML = `<div class="mo-badge"><span class="badge ${cls}">${esc(r.model_class)}</span>
        <span><b>${r.n_states}</b> states${saved > 0 ? ` <span class="muted">(${saved} fewer than the full model's ${r.n_states_full})</span>` : ""}</span></div>
        <p class="mo-note">${esc(MODEL_CLASS_NOTE[r.model_class] || "")}</p>
        ${r.units.length ? `<table class="mo-units"><thead><tr><th>Unit</th><th>Model</th><th>States</th></tr></thead><tbody>
          ${r.units.map(u => `<tr><td>${esc(u.label)}</td><td>${esc(this.levelLabel(u.unit_type, u.level))}${u.exciter ? `<span class="mo-s">${esc(this.regulatorLabel("sm", "exciter", u.exciter))} · ${esc(this.regulatorLabel("sm", "pss", u.pss))} · ${esc(this.regulatorLabel("sm", "governor", u.governor))}</span>`
          : u.controller ? `<span class="mo-s">${esc(this.regulatorLabel("gfm", "controller", u.controller))}</span>` : ""}</td><td>${u.n_states}</td></tr>`).join("")}
        </tbody></table>` : ""}`;
    } catch (e) {
      if (version !== state.version) return;
      out.innerHTML = errorHtml(e);
    }
  },

  levelLabel(kind, levelId) {
    if (levelId === null || levelId === undefined) return "Custom";
    const e = this.catalogue && this.catalogue[kind];
    const l = e && e.levels.find(x => x.id === levelId);
    return l ? l.label : levelId;
  },

  async runAdequacy() {
    const out = $("#mo-adequacy-out");
    const band = parseFloat($("#mo-band").value) || 5;
    setSpinner(out, "Linearising at full order and comparing…");
    const version = state.version;
    try {
      const r = await netPost("model/adequacy", { band_hz: band });
      if (version !== state.version) return;
      out.innerHTML = this.adequacyHtml(r);
    } catch (e) {
      if (version !== state.version) return;
      out.innerHTML = errorHtml(e);
    }
  },

  adequacyHtml(r) {
    const badge = { safe: "ok", check: "warn", unsafe: "crit", full_order: "" }[r.verdict] || "";
    const title = {
      safe: "Safe for this network",
      check: "Usable, with a caveat",
      unsafe: "Not safe for this network",
      full_order: "Full order",
    }[r.verdict] || r.verdict;
    const matched = r.modes.filter(m => m.matched);
    // A mode that vanished because it *was* the removed states is not a
    // defect, so it is not marked as one.
    const rows = r.modes.slice(0, 20).map(m => `<tr${m.matched || m.expected_loss ? "" : ' class="bad-row"'}>
      <td>${fmt(m.full_hz, 4)}</td><td>${fmt(m.full_damping_pct, 3)}</td>
      <td>${m.matched ? fmt(m.reduced_hz, 4) : "—"}</td><td>${m.matched ? fmt(m.reduced_damping_pct, 3) : "—"}</td>
      <td>${m.matched ? (m.d_hz * 1000).toFixed(1)
            : m.expected_loss ? `<span class="muted" title="${(m.removed_share * 100).toFixed(0)}% of this mode was the removed states">removed by design</span>`
            : `<b class="bad" title="only ${(m.removed_share * 100).toFixed(0)}% of this mode was the removed states">lost</b>`}</td></tr>`).join("");
    // The risk table is the screen, not the measurement: it says which
    // removed states had weight in the band, which is why a verdict is what
    // it is (and, when the verdict is good anyway, how close the call was).
    const risks = r.risks.slice(0, 12).map(x => `<tr><td>${esc(x.state)}</td><td class="muted">${esc(x.group)}</td>
      <td>${fmt(x.mode_hz, 3)}</td><td>${(x.participation * 100).toFixed(0)}%</td></tr>`).join("");
    return `
      <div class="notice ${r.verdict === "unsafe" ? "err-bg" : r.verdict === "check" ? "warn-bg" : ""}" style="margin-top:0.8rem">
        <span><b class="badge ${badge}">${esc(title)}</b></span></div>
      <ul class="mo-notes">${r.notes.map(n => `<li>${esc(n)}</li>`).join("")}</ul>
      ${r.modes.length ? `<h4 class="mo-h4">Modes below ${fmt(r.band_hz, 3)} Hz — full model vs. this one${r.modes.length > 20 ? ` (first 20 of ${r.modes.length})` : ""}</h4>
        <div class="tablewrap"><table><thead><tr><th>Full (Hz)</th><th>Damping (%)</th><th>Reduced (Hz)</th><th>Damping (%)</th><th>Shift (mHz) / fate</th></tr></thead>
        <tbody>${rows}</tbody></table></div>` : `<p class="empty">No modes in that band to compare.</p>`}
      ${risks ? `<h4 class="mo-h4">Removed states that carry weight in that band${r.risks.length > 12 ? ` (worst 12 of ${r.risks.length})` : ""}</h4>
        <div class="tablewrap"><table><thead><tr><th>State</th><th>Group</th><th>Mode (Hz)</th><th>Participation</th></tr></thead>
        <tbody>${risks}</tbody></table></div>`
        : (r.verdict === "full_order" ? "" : `<p class="status-line"><span class="ok">No removed state participates above 5% in the band</span> — the timescale separation this reduction assumes holds for this case.</p>`)}
      ${matched.length ? `<p class="status-line muted">${matched.length} of ${r.modes.length} mode(s) matched.</p>` : ""}`;
  },
};

const MODEL_CLASS_NOTE = {
  EMT: "Everything is integrated: the network's electromagnetic transients, the machines' stator flux and the converters' filters. A time-domain run is an EMT simulation, and a stiff one.",
  RMS: "The network is quasi-stationary and no unit keeps its fast electrical states. A time-domain run is an electromechanical (RMS) simulation — far cheaper, and blind to anything above a few tens of Hz.",
  Mixed: "Some elements keep their fast dynamics and others don't. A valid model, and the right one for isolating what a single element contributes — but its results are neither EMT nor RMS.",
};
