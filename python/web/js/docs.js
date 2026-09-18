// Documentation page: the built Sphinx site (served at /manual) embedded in
// place. Its table of contents is read from the built index page and listed
// in the app's own sidebar, under "Documentation"; Sphinx's own sidebars are
// hidden inside the frame so the pages read as part of the app.

const DOCS_ROOT = "/manual/";

const DocsPage = {
  inited: false,
  available: null,
  current: "index.html",
  toc: [],

  async init() {
    if (this.inited) return this.ready;
    this.inited = true;
    this.ready = (async () => {
      try {
        const res = await fetch(`${DOCS_ROOT}index.html`);
        if (!res.ok) throw new Error("not built");
        const doc = new DOMParser().parseFromString(await res.text(), "text/html");
        this.toc = this.parseToc(doc);
        this.available = true;
      } catch {
        this.available = false;
      }
      this.renderNav();
      this.renderFrame();
    })();
    return this.ready;
  },

  // Furo's sidebar tree: captions + toctree-l1 entries; the index page's
  // body toctree additionally lists each page's l2 sections.
  parseToc(doc) {
    const out = [{ caption: null, title: "Overview", href: "index.html", children: [] }];
    const tree = doc.querySelector(".sidebar-tree");
    if (!tree) return out;
    let caption = null;
    tree.childNodes.forEach(node => {
      if (node.nodeType !== 1) return;
      if (node.classList.contains("caption")) caption = node.textContent.trim();
      if (node.tagName === "UL") {
        node.querySelectorAll(":scope > li.toctree-l1 > a").forEach(a => {
          out.push({ caption, title: a.textContent.trim(), html: a.innerHTML, href: a.getAttribute("href"), children: [] });
          caption = null;
        });
      }
    });
    const body = doc.querySelector(".toctree-wrapper");
    if (body) {
      body.querySelectorAll("li.toctree-l1").forEach(li => {
        const a = li.querySelector(":scope > a");
        const entry = a && out.find(e => e.href === a.getAttribute("href"));
        if (!entry) return;
        li.querySelectorAll(":scope > ul > li.toctree-l2 > a").forEach(a2 => entry.children.push({ title: a2.textContent.trim(), href: a2.getAttribute("href") }));
      });
    }
    return out;
  },

  renderNav() {
    const box = $("#nav-docs-sub");
    if (!this.available) { box.innerHTML = `<span class="nav-caption">Not built yet</span>`; return; }
    box.innerHTML = this.toc.map(e => `${e.caption ? `<span class="nav-caption">${esc(e.caption)}</span>` : ""}
      <a href="#/docs/${encodeURIComponent(e.href)}" data-doc="${esc(e.href)}">${e.html || esc(e.title)}</a>
      <div class="nav-doc-children" data-parent="${esc(e.href)}" style="display:none">${e.children.map(c => `<a class="nav-sub-l2" href="#/docs/${encodeURIComponent(c.href)}" data-doc="${esc(c.href)}">${esc(c.title)}</a>`).join("")}</div>`).join("");
    this.highlight();
  },

  renderFrame() {
    const page = $("#page-docs");
    if (!this.available) {
      page.innerHTML = `<div class="docs-missing"><div class="page-head"><div class="crumb">Reference</div><h1>Documentation</h1></div>
        <div class="notice warn-bg"><span>The documentation hasn't been built yet. From the <code>python/</code> directory run <code>pip install -e ".[docs]"</code> then <code>python tools/build_docs.py</code>, and reload this page.</span></div></div>`;
      return;
    }
    page.innerHTML = `<iframe class="docs-frame" id="docs-frame" title="G2ELin documentation" data-src="${DOCS_ROOT}index.html"></iframe>`;
    const frame = $("#docs-frame");
    frame.addEventListener("load", () => this.onFrameLoad(frame));
    frame.src = DOCS_ROOT + this.current;
  },

  // Called by the router with the path after #/docs/.
  async show(path) {
    await this.init();
    if (!this.available) return;
    const target = path ? decodeURIComponent(path) : this.current;
    this.current = target;
    const frame = $("#docs-frame");
    const now = this.framePath(frame);
    if (now !== target) frame.src = DOCS_ROOT + target;
    this.highlight();
  },

  framePath(frame) {
    try {
      const loc = frame.contentWindow.location;
      if (!loc.pathname.startsWith(DOCS_ROOT)) return null;
      return loc.pathname.slice(DOCS_ROOT.length) + loc.hash;
    } catch { return null; }
  },

  onFrameLoad(frame) {
    try {
      const d = frame.contentDocument;
      if (d && !d.getElementById("g2elin-embed")) {
        const st = d.createElement("style");
        st.id = "g2elin-embed";
        st.textContent = `.sidebar-drawer, .mobile-header, .toc-overlay-icon, .sidebar-toggle, .theme-toggle-container { display: none !important; }
          .page { margin: 0 auto; } .main { justify-content: center; } .content { width: min(100%, 60em); padding: 0 2.5em; }
          @media (max-width: 67em) { .content { padding: 0 1.2em; } } .article-container { padding-top: 1.2rem; }`;
        d.head.appendChild(st);
        // Keep in-frame navigation reflected in the app URL + sidebar.
        d.addEventListener("click", e => {
          const a = e.target.closest("a[href]");
          if (a && a.origin === location.origin && !a.pathname.startsWith(DOCS_ROOT)) { e.preventDefault(); location.hash = a.hash || "#/home"; }
        });
        frame.contentWindow.addEventListener("hashchange", () => this.syncFromFrame(frame));
      }
    } catch { /* cross-origin (shouldn't happen: same server) */ }
    this.syncFromFrame(frame);
  },

  syncFromFrame(frame) {
    const p = this.framePath(frame);
    if (!p) return;
    this.current = p;
    const want = `#/docs/${encodeURIComponent(p)}`;
    if (location.hash !== want && location.hash.startsWith("#/docs")) history.replaceState(null, "", want);
    this.highlight();
  },

  highlight() {
    const page = this.current.split("#")[0];
    $$("#nav-docs-sub a[data-doc]").forEach(a => a.classList.toggle("active", a.dataset.doc === this.current || (a.dataset.doc === page && !this.current.includes("#")) || (a.dataset.doc === page && !$$(`#nav-docs-sub a[data-doc="${CSS.escape(this.current)}"]`).length && !a.classList.contains("nav-sub-l2"))));
    $$("#nav-docs-sub .nav-doc-children").forEach(div => { div.style.display = div.dataset.parent === page ? "" : "none"; });
  },
};
