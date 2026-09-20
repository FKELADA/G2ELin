// App shell: hash router, sidebar navigation, preset loading.
//   #/home  #/docs[/<page>]  #/network  #/powerflow  #/modal/<view>  #/emt

const PAGES = {
  home: { el: "page-home", init: () => HomePage.init() },
  docs: { el: "page-docs", init: () => DocsPage.init(), show: arg => DocsPage.show(arg) },
  network: { el: "page-network", init: () => NetworkPage.init(), show: () => NetworkPage.onShow() },
  powerflow: { el: "page-powerflow", init: () => PowerFlowPage.init(), show: () => PowerFlowPage.onShow() },
  modal: { el: "page-modal", init: () => ModalPage.init(), show: arg => ModalPage.show(arg || "eigenmap") },
  emt: { el: "page-emt", init: () => EmtPage.init(), show: () => EmtPage.onShow() },
};

function route() {
  hideTooltip();  // a hover tooltip must not outlive its page
  const raw = location.hash.replace(/^#\/?/, "") || "home";
  const slash = raw.indexOf("/");
  const name = slash < 0 ? raw : raw.slice(0, slash);
  const arg = slash < 0 ? "" : raw.slice(slash + 1);
  const page = PAGES[name] ? name : "home";
  $$(".page").forEach(p => p.classList.toggle("active", p.id === PAGES[page].el));
  PAGES[page].init();
  if (PAGES[page].show) PAGES[page].show(arg);

  $$(".nav-link[data-page]").forEach(a => a.classList.toggle("active", a.dataset.page === page && (page !== "modal" || !a.closest(".nav-group"))));
  // Expand a group when navigating *into* it; don't re-open one the user
  // just collapsed while staying on its page.
  if (page !== route.lastPage) $$(".nav-group").forEach(g => { if (g.dataset.group === page) g.classList.add("open"); });
  route.lastPage = page;
  $$("#nav-modal-sub a").forEach(a => a.classList.toggle("active", page === "modal" && a.dataset.view === (arg || "eigenmap")));
  $("#nav-modal").classList.toggle("active", page === "modal");
  $("#nav-docs").classList.toggle("active", page === "docs");
  $(".sidebar").classList.remove("open");
  if (page !== "docs") window.scrollTo(0, 0);
  document.title = `G2ELin — ${{ home: "Home", docs: "Documentation", network: "Network", powerflow: "Power Flow", modal: "Modal Analysis", emt: "EMT Simulation" }[page]}`;
}

function updateSidebarNetwork() {
  const box = $("#sidebar-network");
  if (!state.network) { box.innerHTML = `<span class="label">Network</span><span class="net-name">none loaded</span>`; return; }
  const n = state.network;
  box.innerHTML = `<span class="label">Current network</span><a href="#/network" class="net-name" style="text-decoration:none">${esc(state.networkLabel)}</a>
    <span class="net-meta">${n.buses.length} buses · ${n.der_units.length} units${isModified() ? ' · <span class="modified">modified</span>' : ""}</span>`;
}

async function boot() {
  // Group headers (Documentation, Modal Analysis): the chevron only toggles
  // the sub-pages; the label navigates to the page and expands it, or --
  // when that page is already open -- collapses/expands it again.
  $$(".nav-group > .nav-link").forEach(link => link.addEventListener("click", evt => {
    const g = link.parentElement;
    const onThisPage = location.hash.replace(/^#\/?/, "").split("/")[0] === g.dataset.group;
    if (evt.target.closest(".chev") || onThisPage) {
      evt.preventDefault();
      g.classList.toggle("open");
    } else {
      g.classList.add("open");
    }
  }));
  $("#menu-toggle").addEventListener("click", () => $(".sidebar").classList.toggle("open"));
  on("network:loaded", updateSidebarNetwork);
  on("network:changed", updateSidebarNetwork);
  updateSidebarNetwork();
  watchForFigures();   // PNG/SVG/CSV buttons on every figure and table (exports.js)
  window.addEventListener("hashchange", route);

  // Docs nav is filled lazily but should be browsable from any page.
  DocsPage.init();
  route();

  try {
    state.presets = await api("/api/presets");
    emit("presets:loaded");
    if (NetworkPage.inited) NetworkPage.fillPresetSelect();
    if (!state.network && state.presets.length) {
      const def = state.presets.find(p => p.id === "wscc9_3sm") || state.presets[0];
      await NetworkPage.loadPresetQuiet(def.id);
    }
  } catch (e) {
    console.error("failed to load presets", e);
  }
}

// Loads a preset without requiring the Network page to be initialised.
NetworkPage.loadPresetQuiet = async function (id) {
  const net = await api(`/api/presets/${encodeURIComponent(id)}/network`);
  const p = state.presets.find(x => x.id === id);
  setNetwork(net, { label: p ? p.name : id, presetId: id });
  if ($("#net-preset")) $("#net-preset").value = id;
};

document.addEventListener("DOMContentLoaded", boot);
