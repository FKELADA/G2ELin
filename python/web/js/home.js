// Home page: what the tool is, what it does, and a quick way in.

const ICONS = {
  network: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="5" cy="6" r="2.5"/><circle cx="19" cy="6" r="2.5"/><circle cx="12" cy="18" r="2.5"/><path d="M7.5 6h9M6.3 8.2l4.5 7.6M17.7 8.2l-4.5 7.6"/></svg>',
  pf: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M13 2 4 14h7l-1 8 9-12h-7z"/></svg>',
  modal: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M3 12h18M12 3v18"/><circle cx="7" cy="8" r="1.4" fill="currentColor"/><circle cx="7" cy="16" r="1.4" fill="currentColor"/><circle cx="5" cy="12" r="1.4" fill="currentColor"/><circle cx="9.5" cy="5.5" r="1.4" fill="currentColor"/><circle cx="9.5" cy="18.5" r="1.4" fill="currentColor"/></svg>',
  emt: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M2 12h3l2-6 3 12 3-9 2 5 2-2h5"/></svg>',
  docs: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M4 4h10l6 6v10H4z"/><path d="M14 4v6h6M8 13h8M8 17h6"/></svg>',
  open: '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><path d="M12 3v12m0 0-4-4m4 4 4-4M4 17v3h16v-3"/></svg>',
};

const HomePage = {
  inited: false,
  init() {
    if (this.inited) return;
    this.inited = true;
    $("#page-home").innerHTML = `
      <section class="hero">
        <div>
          <div class="eyebrow">Open-access · G2ELab · GPL-3</div>
          <h1>Small-signal & EMT analysis for <em>inverter-dominated</em> power systems</h1>
          <p>G2ELin builds the full small-signal model of a power network — synchronous machines, grid-forming and grid-following converters, lines, nodes and loads — linearises it around a solved power flow, and checks it against a nonlinear time-domain (EMT) simulation of the very same equations.</p>
          <div class="cta"><a class="primary" href="#/network">Open a network →</a><a class="outline" href="#/docs">Read the documentation</a></div>
          <div class="hero-stats" id="home-stats">
            <div><span class="n" id="home-n-presets">–</span><span class="l">Preset networks</span></div>
            <div><span class="n">4</span><span class="l">Unit models</span></div>
            <div><span class="n">7</span><span class="l">Modal views</span></div>
            <div><span class="n">6</span><span class="l">Power-flow solvers</span></div>
          </div>
        </div>
        <div class="hero-art" aria-hidden="true">${this.heroSvg()}</div>
      </section>

      <h2 class="section-title">What you can do</h2>
      <p class="section-sub">Each step builds on the previous one — the power flow is the operating point every dynamic study starts from.</p>
      <div class="feature-grid">
        ${this.feature("network", "#/network", "Build or load a network", "Start from one of the WSCC-9, CIGRE MV or SMIB presets, or draw your own on a canvas.",
          ["Drag-and-drop buses, units, loads, lines, transformers", "Edit every parameter in a side inspector", "Live validation, JSON import/export"])}
        ${this.feature("pf", "#/powerflow", "Power flow", "pandapower-based AC power flow with a choice of solvers and settings.",
          ["Newton-Raphson, Iwamoto, fast-decoupled, Gauss-Seidel, BFSW", "Results as a heatmap on the diagram", "Batch runs scaling loads and unit setpoints"])}
        ${this.feature("modal", "#/modal/eigenmap", "Modal analysis", "The whole network linearised into one closed-loop state matrix and eigen-decomposed.",
          ["Eigenvalue map & damping", "Participation, sensitivity, mode shapes", "Free-motion and step responses, multi-channel"])}
        ${this.feature("emt", "#/emt", "EMT simulation", "The nonlinear DAE model integrated in time after a state offset or an input step.",
          ["Coupled Newton solve at every step", "Live tracing of the solver", "Overlay of the linearised response"])}
        ${this.feature("docs", "#/docs", "Documentation", "The physics and architecture reference: equations, methodology, module by module, with worked examples.",
          ["Component equations", "Interconnection & elimination", "Notebook walkthroughs"])}
      </div>

      <h2 class="section-title">Workflow</h2>
      <div class="card"><div class="workflow">
        <div class="wf-step"><h4>Choose a network</h4><p>Pick a preset or build one; every page works on the same network.</p></div>
        <div class="wf-step"><h4>Solve the power flow</h4><p>Get the steady-state operating point — voltages, angles, flows.</p></div>
        <div class="wf-step"><h4>Linearise & analyse</h4><p>Eigenvalues, damping, participation — which control loop drives which mode.</p></div>
        <div class="wf-step"><h4>Verify in the time domain</h4><p>Run the nonlinear EMT model and compare it with the linear prediction.</p></div>
      </div></div>

      <h2 class="section-title">Component models</h2>
      <p class="section-sub">Symbolic nonlinear models, linearised analytically — the same equations drive modal analysis and EMT.</p>
      <div class="models">
        ${Object.keys(UNIT_LABEL).map(u => `<div class="model-chip"><span class="tag" style="background:${UNIT_COLOR[u]}">${UNIT_LABEL[u]}</span>${UNIT_NAME[u]}</div>`).join("")}
        <div class="model-chip"><span class="tag" style="background:#77766f">RL</span>Lines (dq dynamic)</div>
        <div class="model-chip"><span class="tag" style="background:#77766f">C</span>Nodes (shunt capacitance)</div>
        <div class="model-chip"><span class="tag" style="background:#77766f">Z</span>Constant-impedance loads</div>
      </div>

      <h2 class="section-title">Presets</h2>
      <p class="section-sub">Click one to load it on the Network page.</p>
      <div class="preset-grid" id="home-presets"><p class="empty">Loading…</p></div>

      <div class="home-foot"><span>G2ELin — developed by Fadi Kelada at G2ELab. Provided under the GNU GPL-3 licence, without warranty.</span>
        <span>Cite: F. Kelada, “G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool,” 2023.</span></div>`;
    on("presets:loaded", () => this.renderPresets());
    this.renderPresets();
  },

  feature(icon, href, title, text, bullets) {
    return `<a class="feature" href="${href}"><span class="ic">${ICONS[icon]}</span><h3>${title}</h3><p>${text}</p>
      <ul>${bullets.map(b => `<li>${b}</li>`).join("")}</ul><span class="go">Open →</span></a>`;
  },

  renderPresets() {
    const box = $("#home-presets");
    if (!box || !state.presets.length) return;
    $("#home-n-presets").textContent = state.presets.length;
    box.innerHTML = state.presets.map(p => `<button class="preset-card" data-preset="${esc(p.id)}">
      <div class="pc-fam">${esc(presetFamily(p.id))}</div><div class="pc-name">${esc(p.name)}</div>
      <div class="pc-meta"><span class="badge">${p.n_buses} buses</span>${p.der_unit_types.map(t => `<span class="badge" style="color:${UNIT_COLOR[t] || "inherit"}">${UNIT_LABEL[t] || t}</span>`).join("")}</div></button>`).join("");
    $$("#home-presets [data-preset]").forEach(b => b.addEventListener("click", async () => {
      await NetworkPage.loadPreset(b.dataset.preset);
      location.hash = "#/network";
    }));
  },

  // A small stylised one-line diagram with a damped response behind it.
  heroSvg() {
    const nodes = [[70, 190, "sm"], [185, 95, null], [300, 175, null], [410, 90, "gfm"], [250, 285, null], [390, 280, "gfl"], [130, 300, null]];
    const edges = [[0, 1], [1, 2], [2, 3], [2, 4], [4, 5], [4, 6], [6, 0], [1, 3]];
    let path = "";
    for (let x = 0; x <= 480; x += 4) {
      const t = x / 480;
      path += `${x ? "L" : "M"}${x},${(200 + 70 * Math.exp(-3.2 * t) * Math.sin(2 * Math.PI * 4.2 * t)).toFixed(1)}`;
    }
    return `<svg viewBox="0 0 480 380">
      <defs><radialGradient id="hg" cx="50%" cy="50%" r="50%"><stop offset="0" stop-color="#2a78d6" stop-opacity="0.35"/><stop offset="1" stop-color="#2a78d6" stop-opacity="0"/></radialGradient></defs>
      <circle cx="260" cy="190" r="190" fill="url(#hg)"/>
      <path d="${path}" fill="none" stroke="#8fb8ea" stroke-opacity="0.35" stroke-width="2"/>
      ${edges.map(([a, b]) => `<line x1="${nodes[a][0]}" y1="${nodes[a][1]}" x2="${nodes[b][0]}" y2="${nodes[b][1]}" stroke="#5d6b7d" stroke-width="2.5"/>`).join("")}
      ${nodes.map(([x, y, u]) => u
        ? `<circle cx="${x}" cy="${y}" r="21" fill="${UNIT_COLOR[u]}" stroke="#12161c" stroke-width="3"/><text x="${x}" y="${y + 4}" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="12" font-weight="700" fill="#fff">${UNIT_LABEL[u]}</text>`
        : `<circle cx="${x}" cy="${y}" r="9" fill="#c9ced6" stroke="#12161c" stroke-width="3"/>`).join("")}
      <polygon points="243,310 257,310 250,322" fill="#c9ced6"/><polygon points="123,325 137,325 130,337" fill="#c9ced6"/>
    </svg>`;
  },
};
