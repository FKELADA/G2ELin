// An exhaustive, scripted and captioned tour of the G2ELin web interface,
// recorded to WebM. See README.md in this directory for how to run it and
// why it is built the way it is.
//
// Every scene is independent: one that fails is logged and skipped, so the
// recording always completes and the log says what was missed.
//
// Live tracing stays on for the time-domain runs. Its pace is the solver's
// and nothing here can hurry it, so each live run is timestamped into
// segments.json and encode.py speeds exactly those stretches up.
const { chromium } = require("playwright-core");
const path = require("path");
const fs = require("fs");

const BASE = process.env.TOUR_BASE || "http://127.0.0.1:8765/";
const PRESET = "cigre_interconnected_1sm_1gfm_1gfl";
const W = 1600, H = 900;
const OUT_DIR = process.env.TOUR_OUT || path.join(__dirname, "recording");
// Playwright drives real Chrome; point at it if it is not installed here.
const CHROME = process.env.TOUR_CHROME || "C:/Program Files/Google/Chrome/Application/chrome.exe";
fs.rmSync(OUT_DIR, { recursive: true, force: true });

const OVERLAY = fs.readFileSync(path.join(__dirname, "overlay.js"), "utf8");
const LOGO = '<svg class="logo" viewBox="0 0 32 32"><rect width="32" height="32" rx="8" fill="#2a78d6"/><path d="M5 17h4l2.5-7 4 13 3.5-10 2 4H27" fill="none" stroke="#fff" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/></svg>';

(async () => {
  const browser = await chromium.launch({ executablePath: CHROME });
  const ctx = await browser.newContext({ viewport: { width: W, height: H }, recordVideo: { dir: OUT_DIR, size: { width: W, height: H } } });
  const t0 = Date.now();
  const p = await ctx.newPage();
  p.on("pageerror", e => console.log("pageerror", e.message));
  p.on("dialog", d => d.accept());
  const now = () => (Date.now() - t0) / 1000;
  const log = m => console.log(`${now().toFixed(1)}s ${m}`);

  // Stretches of the recording to speed up afterwards: the live traces,
  // which run at whatever pace the solver streams at.
  const fast = [];
  const faster = async (speed, fn) => {
    const a = now();
    try { return await fn(); }
    finally { fast.push({ start: a, end: now(), speed }); }
  };

  const wait = ms => p.waitForTimeout(ms);
  const moveTo = async (x, y, steps = 28) => { await p.mouse.move(x, y, { steps }); };
  const center = async sel => {
    const el = await p.waitForSelector(sel, { state: "visible", timeout: 60000 });
    await el.scrollIntoViewIfNeeded();
    const b = await el.boundingBox();
    return { x: b.x + b.width / 2, y: b.y + b.height / 2, b };
  };
  const hover = async (sel, steps) => { const c = await center(sel); await moveTo(c.x, c.y, steps); return c; };
  const click = async (sel, pause = 350) => { await hover(sel); await wait(180); await p.mouse.down(); await p.mouse.up(); await wait(pause); };
  const clickBox = async (box, pause = 400) => {
    await moveTo(box.x + box.width / 2, box.y + box.height / 2, 25);
    await wait(200); await p.mouse.down(); await p.mouse.up(); await wait(pause);
  };
  const caption = (k, t, s) => p.evaluate(([k, t, s]) => window.__tour.caption(k, t, s), [k, t, s]);
  const scrollTo = async (y, ms = 900) => {
    await p.evaluate(([y, ms]) => new Promise(res => {
      const y0 = window.scrollY, t0 = performance.now();
      const step = n => { const f = Math.min(1, (n - t0) / ms), e = f < 0.5 ? 2 * f * f : 1 - Math.pow(-2 * f + 2, 2) / 2;
        window.scrollTo(0, y0 + (y - y0) * e); if (f < 1) requestAnimationFrame(step); else res(); };
      requestAnimationFrame(step);
    }), [y, ms]);
  };
  const scrollToEl = async (sel, offset = 90) => {
    const y = await p.$eval(sel, (el, off) => el.getBoundingClientRect().top + window.scrollY - off, offset);
    await scrollTo(Math.max(0, y));
  };
  const chapter = async (title, sub) => {
    await p.evaluate(() => window.__tour.hideCaption());
    await p.evaluate(h => window.__tour.card(h), `<h2>${title}</h2><p>${sub}</p>`);
    await wait(2400);
    await p.evaluate(() => window.__tour.hideCard());
    await wait(600);
  };
  const scene = async (name, fn) => {
    const t = Date.now();
    try { await fn(); } catch (e) { console.log(`!! scene ${name} failed: ${e.message.split("\n")[0]}`); }
    // A caption belongs to its scene. Left standing it outlives what it
    // describes -- the RMS comparison's conclusion sat over the network page
    // while the next scene was already navigating.
    await p.evaluate(() => window.__tour.hideCaption()).catch(() => {});
    log(`${name} (${((Date.now() - t) / 1000).toFixed(1)}s)`);
  };
  const sweepChart = async (sel, ms = 2200) => {
    const b = (await center(sel)).b;
    await moveTo(b.x + b.width * 0.12, b.y + b.height * 0.5, 20);
    await moveTo(b.x + b.width * 0.85, b.y + b.height * 0.5, Math.round(ms / 40));
  };
  const breaker = async (host, tag) => {
    const el = await p.waitForSelector(`${host} [data-brk="${tag}"]`, { timeout: 30000 });
    await el.scrollIntoViewIfNeeded();
    return await el.boundingBox();
  };
  // The diagram's breaker hit areas overlap, so a click meant for one can
  // toggle its neighbour. The click stays (it is what the viewer sees); this
  // then states outright what should be out of service, so a stray toggle
  // cannot follow the tour into the next chapter.
  const setService = async ({ lines = [], transformers = [], loads = [], units = [] } = {}) => {
    await p.evaluate(open => {
      const n = state.network;
      n.lines.forEach((l, i) => { const o = open.lines.includes(i); l.from_closed = !o; l.to_closed = !o; });
      n.transformers.forEach((t, i) => { const o = open.transformers.includes(i); t.hv_closed = !o; t.lv_closed = !o; });
      n.loads.forEach((l, i) => { l.closed = !open.loads.includes(i); });
      n.der_units.forEach(d => { d.closed = !open.units.includes(d.id); });
      networkChanged();
    }, { lines, transformers, loads, units });
    await wait(700);
  };

  const loadPreset = async id => {
    await p.goto(BASE + "#/network");
    await wait(700);
    await p.selectOption("#net-preset", id);
    await p.waitForFunction(() => document.querySelector("#net-canvas svg"), null, { timeout: 60000 });
    await wait(1600);
  };
  // The signal pickers everywhere take element kind -> element -> signal.
  const pick = async (prefix, kind, el, signals) => {
    await p.selectOption(`${prefix}-el-kind`, kind); await wait(450);
    if (el != null) { await p.selectOption(`${prefix}-el`, String(el)); await wait(450); }
    for (const s of signals) {
      await p.evaluate(([pre, name]) => {
        const box = document.querySelector(`${pre}-sig`) || document;
        const hit = [...box.querySelectorAll("input[type=checkbox]")]
          .find(c => (c.value || c.dataset.sig || c.nextElementSibling?.textContent || "").includes(name));
        if (hit && !hit.checked) hit.click();
      }, [prefix, s]);
      await wait(200);
    }
  };

  await p.goto(BASE + "#/home");
  await p.addScriptTag({ content: OVERLAY });
  await p.evaluate(h => window.__tour.card(h), `${LOGO}<h1>G2ELin</h1><p>Small-signal &amp; EMT analysis for inverter-dominated power systems</p><div class="m">A complete tour of the web interface</div>`);
  await wait(3400);
  await p.evaluate(() => window.__tour.hideCard());
  await wait(600);

  // ============================================================ HOME
  await scene("home", async () => {
    await caption("Home", "Build, solve, linearise and simulate — in one place", "Everything is reachable from the navigation panel on the left");
    await moveTo(700, 360);
    await wait(2200);
    await scrollTo(620, 1400);
    await caption("Home", "Features, workflow and component models", "Synchronous machines · grid-forming and grid-following converters · infinite bus · lines · loads");
    await wait(2600);
    await scrollTo(1500, 1400);
    await caption("Home", "Preset networks, ready to load", "Single-machine cases · WSCC 9-bus · CIGRE MV · Kundur's two-area system · IEEE 14, 39 and 118 bus");
    await wait(3000);
  });

  // ======================================================= NEW NETWORKS
  await chapter("1 · The networks", "From one machine on an infinite bus to 118 buses");

  await scene("net-kundur", async () => {
    await loadPreset("kundur_two_area");
    await caption("Networks", "Kundur's two-area system", "The textbook inter-area case: four machines, two areas, a long tie — with the book's own exciter, stabiliser and machine data");
    await wait(3200);
    await caption("Networks", "Its published operating point, reproduced", "Every generator terminal lands within 0.1° of the book, angles included");
    await wait(2600);
  });

  await scene("net-ieee", async () => {
    await loadPreset("ieee118");
    await caption("Networks", "IEEE 118-bus", "Converted from the published case: a terminal bus and step-up per generator, 54 machines, 1837 states at full order");
    await wait(3400);
    await loadPreset("ieee39");
    await caption("Networks", "IEEE 39-bus — New England", "IEEE 14 and 39 are here too; every one of them linearises, simulates and sweeps like any other");
    await wait(2800);
  });

  // ============================================================ NETWORK
  await chapter("2 · The network", "Inspect it, edit it, or build one from scratch");

  await scene("network-plot", async () => {
    await loadPreset(PRESET);
    await caption("Network", "The CIGRE MV benchmark, interconnected", "Two feeders from the 110 kV grid, with a synchronous machine, a grid-forming and a grid-following converter");
    await wait(2800);
    await caption("Network", "Every unit is drawn as its own symbol", "An AC source for the machine; the converters carry a sine on the AC side if they form the voltage, an arrow if they follow it");
    await hover("#net-canvas");
    await wait(3000);
    await caption("Network", "Hover any element for its parameters", "Buses with their units and loads · lines · transformers");
    const b = (await center("#net-canvas svg")).b;
    await moveTo(b.x + b.width * 0.35, b.y + b.height * 0.45, 30);
    await wait(1700);
    await moveTo(b.x + b.width * 0.62, b.y + b.height * 0.6, 30);
    await wait(1900);
  });

  await scene("network-inspector", async () => {
    const gfmBus = await p.evaluate(() => state.network.der_units.find(u => u.unit_type === "gfm").bus);
    await click(`#net-canvas [data-bus-id="${gfmBus}"]`, 1200);
    await caption("Network", "Click an element to open its inspector", "Every editable field of that element — and of its unit and loads");
    await wait(2400);
    await scrollToEl("#insp-ctrl", 120);
    await caption("Network", "A unit's control and electrical parameters", "Filter, DC link, current and voltage loops — stored as overrides on top of the derived defaults");
    await wait(2800);
  });

  // =============================================== CONTROL LAWS (NEW)
  await scene("network-control-law", async () => {
    await caption("Network", "Which control law the converter runs is a choice", "Droop · droop behind a filter · dVOC · a virtual synchronous machine · matching control");
    await hover("select[data-regulator]");
    await wait(2200);
    for (const law of ["vsm", "matching", "dvoc", "droop"]) {
      await p.selectOption("select[data-regulator]", law).catch(() => {});
      await wait(1500);
    }
    await caption("Network", "Each law brings its own parameters", "A VSM has inertia and damping where droop has mp and nq — nothing is shared, and an override for the wrong law is refused");
    await p.selectOption("select[data-regulator]", "vsm").catch(() => {});
    await wait(3000);
    await caption("Network", "And they share one tuning", "Every law is written from the same droop equivalence, so the inertia can be set on a VSM that has no mp at all — and swapping a law keeps it");
    await hover('input[data-tuning="H"]').catch(() => {});
    await wait(2600);
    await p.fill('input[data-tuning="H"]', "6").catch(() => {});
    await p.dispatchEvent('input[data-tuning="H"]', "change").catch(() => {});
    await wait(2600);
    await p.selectOption("select[data-regulator]", "droop").catch(() => {});
    await wait(1800);
  });

  await scene("network-machine-regulators", async () => {
    const smBus = await p.evaluate(() => state.network.der_units.find(u => u.unit_type === "sm").bus);
    await click(`#net-canvas [data-bus-id="${smBus}"]`, 1200);
    await scrollToEl("#insp-ctrl", 120);
    await caption("Network", "A machine's regulators are chosen the same way", "Exciter, stabiliser and governor — each with its own parameters, under its own names");
    await wait(3000);
    const sel = 'select[data-regulator="exciter"]';
    await hover(sel).catch(() => {});
    await p.selectOption(sel, "kundur").catch(() => {});
    await wait(2400);
    await caption("Network", "Kundur's thyristor exciter, from the book", "A transducer, a gain and transient gain reduction — TR, KA, TA, TB, exactly as Fig. E12.9 labels them");
    await wait(3000);
    await caption("Network", "A stabiliser and a governor can simply be absent", "Not a gain turned down: the states, the parameters and the reduction group go with them");
    await p.selectOption('select[data-regulator="pss"]', "none").catch(() => {});
    await wait(1600);
    await p.selectOption('select[data-regulator="governor"]', "none").catch(() => {});
    await wait(2600);
    await p.selectOption('select[data-regulator="exciter"]', "g2elin").catch(() => {});
    await p.selectOption('select[data-regulator="pss"]', "g2elin").catch(() => {});
    await p.selectOption('select[data-regulator="governor"]', "g2elin").catch(() => {});
    await wait(1200);
  });

  await scene("network-tuner", async () => {
    // The tuner belongs to the converters' PI loops, so go back to one.
    const gfmBus = await p.evaluate(() => state.network.der_units.find(u => u.unit_type === "gfm").bus);
    await click(`#net-canvas [data-bus-id="${gfmBus}"]`, 1200);
    await wait(1200);
    const tuner = await p.$(".tuner");
    if (!tuner) return;
    await tuner.scrollIntoViewIfNeeded();
    await caption("Network", "Loop tuner: response time and damping ⟷ Kp / Ki", "Type a settling time and the gains follow, by the same pole placement the defaults use");
    await wait(3000);
  });

  await scene("network-breakers", async () => {
    await p.goto(BASE + "#/network"); await wait(1200);
    await caption("Network", "Breakers everywhere", "Both ends of every line and transformer, and every load and unit — filled = closed, hollow red = open");
    await scrollToEl("#net-canvas", 110);
    await wait(2600);
  });

  // ====================================== MODEL ORDER REDUCTION (NEW)
  await chapter("3 · Model order", "What makes the same case EMT or RMS");

  await scene("model-order", async () => {
    await scrollToEl("#net-model-order", 100);
    await caption("Model order", "Every element's dynamics are selectable", "Full detail is an EMT model; making the network algebraic and the machines lower-order turns the same case into an RMS study");
    await wait(3400);
    await hover("#mo-summary").catch(() => {});
    await caption("Model order", "The network: dynamic or quasi-stationary", "Dynamic integrates every L·di/dt and C·dv/dt — quasi-stationary is the admittance solve a phasor tool does");
    await wait(3000);
    await caption("Model order", "And a machine by its named order", "8th, 6th, 5th, 4th, 3rd, 2nd — or state group by state group: stator flux, each damper, the field winding");
    await scrollTo(await p.$eval("#net-model-order", e => e.getBoundingClientRect().top + window.scrollY + 260), 900);
    await wait(3600);
    await caption("Model order", "Algebraic and frozen are not the same thing", "One says the state settles instantly, the other that it never moves — the classical machine model needs the second, and choosing wrong is the usual way a reduced model goes quietly wrong");
    await wait(3600);
  });

  await scene("model-order-adequacy", async () => {
    await caption("Model order", "And the tool checks the reduction for you", "It builds the same network at full order and compares — on your network, not in general");
    await scrollToEl("#mo-check", 160);
    await wait(2000);
    await click("#mo-check", 600);
    await p.waitForFunction(() => (document.querySelector("#mo-adequacy-out")?.textContent || "").length > 40, null, { timeout: 300000 });
    await wait(1200);
    await caption("Model order", "Mode by mode, in the band you care about", "Which modes survived, how far they moved, and which were lost — with a verdict for the band");
    await scrollToEl("#mo-adequacy-out", 140);
    await wait(4200);
  });

  // ============================================================ POWER FLOW
  await chapter("4 · Power flow", "The operating point everything else is linearised about");

  const solve = async (timeout = 180000) => {
    await click("#pf-run", 200);
    await p.waitForSelector("#pf-show", { timeout });
    await wait(800);
  };

  await scene("pf-single", async () => {
    await click("a.nav-link[data-page='powerflow']", 900);
    await caption("Power Flow", "Pick the solver and its settings", "Newton-Raphson, Iwamoto, fast-decoupled, Gauss-Seidel, backward/forward sweep · iterations, tolerance, initialisation");
    await hover("#pf-algo").catch(() => {});
    await wait(2400);
    await caption("Power Flow", "Solve the initial operating point", "");
    await solve();
    await caption("Power Flow", "Results drawn on the network as a heatmap", "Buses: low voltage blue → high voltage red · branches shaded by the chosen flow");
    await scrollToEl("#pf-canvas", 230).catch(() => {});
    await wait(2600);
    const bus = await p.evaluate(() => {
      const c = document.querySelectorAll("#pf-canvas .nd-bus circle.core");
      if (!c.length) return null;
      const r = c[Math.min(6, c.length - 1)].getBoundingClientRect();
      return { x: r.x + r.width / 2, y: r.y + r.height / 2 };
    });
    if (bus) {
      await caption("Power Flow", "Hover an element for its own results", "Voltage, angle and injections on buses · flows, losses and current on branches");
      await moveTo(bus.x, bus.y, 30);
      await wait(2400);
    }
    await caption("Power Flow", "Full result tables", "Buses, lines, transformers, loads, generators, the slack — the “Show” selector sits right above the table");
    await scrollToEl("#pf-show", 150).catch(() => {});
    await wait(2800);
  });

  await scene("pf-open-line", async () => {
    await scrollToEl("#pf-canvas", 200).catch(() => {});
    await caption("Power Flow", "Breakers can be operated here too", "Open both ends of line #5 and solve again");
    for (const tag of ["line:4:from", "line:4:to"]) {
      const b = await breaker("#pf-canvas", tag).catch(() => null);
      if (b) await clickBox(b, 700);
    }
    await setService({ lines: [4] });
    await solve();
    await caption("Power Flow", "Line #5 is out of service", "It carries nothing, the feeder is fed the long way round, and the losses rise");
    await scrollToEl("#pf-canvas", 200).catch(() => {});
    await wait(3200);
  });

  await scene("pf-islanding", async () => {
    await caption("Power Flow", "Now open transformer #3, the grid infeed (bus 1 → 18)", "That disconnects the infinite bus — the slack itself");
    for (const tag of ["transformer:2:hv", "transformer:2:lv"]) {
      const b = await breaker("#pf-canvas", tag).catch(() => null);
      if (b) { await clickBox(b, 700); break; }
    }
    await setService({ lines: [4], transformers: [2] });
    await solve();
    await caption("Power Flow", "The network splits into two islands", "Each is solved against its own reference: the synchronous machine takes over the role the infinite bus had");
    await scrollToEl("#pf-canvas", 200).catch(() => {});
    await wait(3600);
    await caption("Power Flow", "Two references, one per island", "The machine now carries the MV network; the infinite bus sits alone on its side");
    await scrollToEl("#pf-show", 150).catch(() => {});
    await wait(3000);
    await setService();
    await caption("Power Flow", "Close them again", "");
    await solve();
    await wait(1600);
  });

  // ============================================================ MODAL
  await chapter("5 · Modal analysis", "Eigenvalues, participation, sensitivity, response");

  await scene("modal-eig", async () => {
    await p.goto(BASE + "#/modal/eigenmap"); await wait(1500);
    await p.waitForFunction(() => document.querySelector("#eig-map svg"), null, { timeout: 300000 }).catch(() => {});
    await caption("Modal analysis", "The whole network, linearised and eigen-decomposed", "Each view is its own sub-page in the navigation");
    await wait(2600);
    await caption("Modal analysis", "Modes grouped by what they are", "Synchronisation, control, unit electrical, network — a different shape each, so a glance says what kind of mode it is");
    await wait(3000);
    await sweepChart("#eig-map svg", 2000).catch(() => {});
    await wait(1800);
  });

  await scene("modal-views", async () => {
    await click("#nav-modal-sub a[data-view='participation']", 1400);
    await caption("Modal analysis", "Participation factors", "Which states drive which mode, for every mode at once");
    await wait(2600);
    await click("#nav-modal-sub a[data-view='single']", 1400);
    await caption("Modal analysis", "Single-mode participation", "The same, mode by mode, ranked");
    await wait(2600);
    await click("#nav-modal-sub a[data-view='shape']", 1400);
    await caption("Modal analysis", "Mode shapes", "Relative amplitude and phase of the states swinging in one mode");
    await wait(2600);
  });

  // =================================== PARAMETER SENSITIVITY (NEW)
  await scene("modal-sensitivity", async () => {
    await click("#nav-modal-sub a[data-view='sensitivity']", 1600);
    await p.waitForFunction(() => (document.querySelector("#sens-out")?.textContent || "").length > 40, null, { timeout: 300000 }).catch(() => {});
    await caption("Modal analysis", "Sensitivity: which entries of A a mode cares about", "∂λ/∂A_ij, ranked — and the heatmap over the whole matrix");
    await wait(3000);
    await caption("Modal analysis", "But A_ij is not a thing anyone can change", "It is an expression in the physical parameters, and that is the question actually being asked");
    await wait(3000);
    await scrollToEl("#param-run", 200).catch(() => {});
    await click("#param-run", 600);
    await p.waitForFunction(() => (document.querySelector("#param-out")?.textContent || "").length > 60, null, { timeout: 600000 });
    await caption("Modal analysis", "So it traces the mode back to the parameter", "Every unit parameter perturbed, ranked by how far it moves this mode — in millihertz and in damping points");
    await wait(4000);
    await scrollToEl("#param-out", 150).catch(() => {});
    await wait(3600);
    await caption("Modal analysis", "And it names the parameters inside each entry", "The “built from” column: which parameters that entry of A is actually made of — an electromechanical mode comes back to H, as the swing equation says it must");
    await scrollToEl("#sens-out", 140).catch(() => {});
    await wait(4000);
  });

  await scene("modal-free", async () => {
    await click("#nav-modal-sub a[data-view='free']", 1400);
    await caption("Modal analysis", "Free-motion response", "Release the system from a perturbed state and watch the modes decay");
    await wait(2400);
    await click("#ch-plot-all", 500).catch(() => {});
    await p.waitForSelector("#ch-list .chart svg", { timeout: 300000 }).catch(() => {});
    await wait(2600);
    await caption("Modal analysis", "Any plot can be windowed", "Pick X, Y or a box and drag — the band shows exactly what the release keeps; double-click goes back");
    const plot = await p.$(".chart svg");
    if (plot) {
      const b = await plot.boundingBox();
      await moveTo(b.x + b.width * 0.25, b.y + b.height * 0.5, 22);
      await p.mouse.down();
      await moveTo(b.x + b.width * 0.55, b.y + b.height * 0.5, 30);
      await p.mouse.up();
      await wait(2600);
      await p.mouse.dblclick(b.x + b.width * 0.5, b.y + b.height * 0.5);
      await wait(1600);
    }
  });

  await scene("modal-step", async () => {
    await click("#nav-modal-sub a[data-view='step']", 1400);
    await caption("Modal analysis", "Step response", "A 0.15 pu step in the grid-following converter's power reference");
    await wait(2600);
    await p.evaluate(() => {
      // A 0.15 pu step in the grid-following converter's power reference,
      // against every unit's electrical power and the PLL frequency.
      const ch = ModalPage.channels.step[0];
      if (!ch) return;
      const opt = (ModalPage.stepInputs || []).map(o => o.name || o);
      const inp = opt.find(n => /GFL/.test(n) && /(p_ref|Idc_ref)/.test(n));
      if (inp) ch.input = inp;
      ch.amplitude = 0.15;
      const outs = (ModalPage.stepOutputs || []).map(o => o.name || o);
      const pe = outs.filter(n => /^p_e/.test(n));
      const pll = outs.filter(n => /w_pll/.test(n));
      if (pe.length) ch.subplots = [pe, pll.length ? pll : pe.slice(0, 1)];
      ModalPage.renderChannels && ModalPage.renderChannels("step");
    }).catch(() => {});
    await wait(1200);
    await click("#ch-plot-all", 500).catch(() => {});
    await p.waitForSelector("#ch-list .chart svg", { timeout: 300000 }).catch(() => {});
    await wait(2800);
    await caption("Modal analysis", "Every unit's electrical power, and the PLL frequency", "Linearised, so this is the small-signal answer — the nonlinear one comes next");
    await wait(3000);
  });

  // ============================================== ROOT LOCUS
  await scene("modal-rootlocus", async () => {
    await caption("Modal analysis", "Root locus — sweep any parameter", "Network, bus, line, transformer, load or unit — including each control loop's response time");
    await click("#nav-modal-sub a[data-view='root']", 1400);
    await wait(1400);
    const row = "#rl-rows .rl-row:nth-child(1) ";
    await p.selectOption(row + "select[data-r='el']", "unit"); await wait(700);
    const gfmKey = await p.evaluate(() => String(state.network.der_units.find(d => d.unit_type === "gfm").id));
    await p.selectOption(row + "select[data-r='key']", gfmKey); await wait(900);
    await caption("Modal analysis", "The grid-forming converter's inner current loop", "Its response time, swept from 0.05 to 1.5 ms in 0.05 ms steps");
    await p.selectOption(row + "select[data-r='field']", "tune.cl.tr_ms"); await wait(900);
    await p.fill(row + "input[data-r='from']", "0.05");
    await p.fill(row + "input[data-r='to']", "1.5");
    await p.fill(row + "input[data-r='step']", "0.05");
    await hover("#rl-count");
    await wait(2000);
    await caption("Modal analysis", "Run the sweep", "Power flow and eigenvalues are re-solved at every value and streamed back");
    await click("#rl-run", 400);
    // Keep the plot in frame while it fills, rather than watching a button.
    await scrollToEl("#rl-main", 110);
    await p.waitForFunction(() => /done|stopped/.test(document.querySelector("#rl-progress")?.textContent || ""), null, { timeout: 900000 });
    await wait(1200);
    await caption("Modal analysis", "The loci, coloured by the parameter value", "Modes are tracked step to step, so each one draws as a single locus");
    await scrollToEl("#rl-main", 110);
    await wait(3000);
    await caption("Modal analysis", "Slowing that loop down drives a mode unstable", "The table ranks what moves most, and flags the value where it crosses over");
    await scrollToEl("#rl-table", 110);
    await wait(3000);

    // Pick the mode that actually moves: the one whose damping falls furthest
    // over the sweep. Mode 117 is that mode on this case; if the table has
    // moved on, fall back to whatever the table ranks first.
    const picked = await p.evaluate(() => {
      const want = document.querySelector('#rl-table tr[data-k="117"]');
      const row = want || document.querySelector("#rl-table tbody tr[data-k]");
      if (!row) return null;
      row.scrollIntoView({ block: "center" });
      row.click();
      return row.dataset.k;
    });
    await wait(1200);
    // Show the locus and the participation panel together, so the bars and
    // the moving ring on the plot are read as one thing.
    await p.evaluate(() => {
      const plot = document.querySelector("#rl-main"), part = document.querySelector("#rl-part");
      if (!plot || !part) return;
      const top = plot.getBoundingClientRect().top + window.scrollY;
      const bot = part.getBoundingClientRect().bottom + window.scrollY;
      window.scrollTo({ top: Math.max(0, (top + bot) / 2 - window.innerHeight / 2), behavior: "smooth" });
    });
    await wait(1400);
    await caption("Modal analysis", `Clicking a mode plays its participation along the sweep`, `Mode ${picked || "117"}: the states that make it up, moving as the loop is detuned — the ring on the locus marks the value shown`);
    await wait(11000);
    await caption("Modal analysis", "Both playbacks download as video", "The locus filling in, and the participation bars moving — straight from the page");
    await hover("#rl-video").catch(() => {});
    await wait(2600);
  });

  // ====================================================== TIME DOMAIN
  await chapter("6 · Time-domain simulation", "The full nonlinear model, integrated in time");

  // Run with the live trace on, and mark the stretch so the encoder speeds
  // it up: the pace is the solver's, and nothing in the browser can hurry it.
  const runLive = async (label, sub, ms = 900000) => {
    await caption("Time domain", label, sub);
    await click("#emt-run", 200);
    await faster(5, async () => {
      await p.waitForSelector("#emt-results .chart", { timeout: ms }).catch(() => {});
      await p.waitForFunction(() => {
        const b = document.querySelector("#emt-run");
        return b && !b.disabled;
      }, null, { timeout: ms }).catch(() => {});
      await wait(500);
    });
  };

  // The three scopes the runs are watched through: unit powers, bus
  // voltages, measured bus frequencies.
  const setScopes = (withPll) => p.evaluate(pll => {
    const opt = EmtPage.options.map(o => o.name);
    const pe = opt.filter(n => n.startsWith("o:p_e"));
    const w = pll ? opt.filter(n => n.includes("w_pll")) : [];
    const v = ["1", "3", "12"].map(b => `m:V_{bus${b}}`).filter(n => opt.includes(n));
    const f = ["1", "3", "12"].map(b => `m:f_{bus${b}}`).filter(n => opt.includes(n));
    EmtPage.scopes = [
      { ...EmtPage.newScope(pe.concat(w)), title: pll ? "Unit powers and PLL frequency" : "Unit powers" },
      { ...EmtPage.newScope(v), title: "Bus voltages (pu)" },
      { ...EmtPage.newScope(f), title: "Bus frequencies (Hz, measured)" },
    ];
    EmtPage.renderScopes();
  }, withPll);
  const liveOn = () => p.evaluate(() => { const c = document.querySelector("#emt-live"); if (c && !c.checked) c.click(); });

  await scene("emt-setup", async () => {
    // Whatever the power-flow scenes were doing with breakers, the
    // time-domain chapter starts from an intact network.
    await p.evaluate(() => {
      const n = state.network;
      n.lines.forEach(l => { l.from_closed = true; l.to_closed = true; });
      n.transformers.forEach(t => { t.hv_closed = true; t.lv_closed = true; });
      n.loads.forEach(l => { l.closed = true; });
      n.der_units.forEach(d => { d.closed = true; });
      networkChanged();
    });
    await wait(800);
    await click("a.nav-link[data-page='emt']", 700);
    await p.waitForFunction(() => !document.querySelector("#emt-perturb")?.disabled, null, { timeout: 300000 });
    await wait(1400);
    await caption("Time domain", "The same model, integrated rather than linearised", "A coupled Newton solve at every step, with the network's own electromagnetic transients in it");
    await wait(2800);
    await scrollToEl("#emt-model-class", 120).catch(() => {});
    await caption("Time domain", "The page says which kind of model this is", "EMT, RMS or mixed — read off the model order rather than asserted");
    await wait(3000);
    await scrollTo(0, 700);
    await caption("Time domain", "Disturb a state, step an input…", "Chosen the same way: element type, element, then the signal");
    await p.selectOption("#emt-el-kind", "unit"); await wait(600);
    const gflKey = await p.evaluate(() => { const e = elementCatalog().find(x => x.block.startsWith("GFL")); return e && e.key; });
    await p.selectOption("#emt-el", gflKey); await wait(800);
    await hover("#emt-perturb");
    await wait(1600);
    await p.selectOption("#emt-perturb", await p.evaluate(() =>
      [...document.querySelectorAll("#emt-perturb option")].map(o => o.value).find(v => v.includes("Idc_ref"))));
    await p.fill("#emt-amp", "0.15");
    await p.fill("#emt-tf", "0.4");
    await p.fill("#emt-dt", "0.002");
    await wait(900);
    await setScopes(true);
    await caption("Time domain", "Scopes: one panel per group of signals", "States, inputs, outputs and measurements together — measurements read like instruments");
    await wait(2600);
    await liveOn();
    await caption("Time domain", "Live tracing is on", "The solver streams its steps and the scopes fill as it goes, rather than waiting for the end");
    await hover("#emt-live").catch(() => {});
    await wait(2400);
  });

  await scene("emt-step-emt", async () => {
    await caption("Time domain", "The same 0.15 pu power step, on the nonlinear model", "Every unit's electrical power, some bus voltages, and the measured frequencies");
    await wait(2400);
    await runLive("Integrating — full EMT model", "Network dynamics, stator flux, filter currents: all of it is in this one");
    await scrollToEl("#emt-results .card:nth-child(2)", 120).catch(() => {});
    await caption("Time domain", "The converter takes up its step", "and the rest of the network answers — with the electromagnetic transient visible in the first cycles");
    await sweepChart("#emt-results .chart svg", 2400).catch(() => {});
    await wait(2600);
    await scrollToEl("#emt-results .card:nth-child(4)", 110).catch(() => {});
    await caption("Time domain", "Bus voltages and measured frequencies", "Read off the model exactly as a meter would, with a one-cycle filter on the frequency");
    await wait(3000);
  });

  // =============================== EMT vs RMS (NEW)
  await scene("emt-vs-rms", async () => {
    await caption("Time domain", "Now the same disturbance as an RMS study", "Network quasi-stationary, stator flux algebraic, transformer currents algebraic — the model a phasor tool carries");
    await wait(3400);
    await p.goto(BASE + "#/network"); await wait(1200);
    await scrollToEl("#net-model-order", 100);
    await caption("Model order", "One selector for the network…", "Quasi-stationary: every passive element becomes algebraic, and the interconnection's elimination becomes the admittance solve");
    const netLevel = '#net-model-order [data-mo-kind="network"] .mo-level';
    await p.selectOption(netLevel, "quasi_stationary").catch(() => {});
    await wait(2400);
    await caption("Model order", "…and one per unit type", "The machines to 6th order — stator flux algebraic, which is the standard stability model — and the converters to their outer loop");
    for (const [kind, level] of [["sm", "order6"], ["gfm", "droop"], ["gfl", "pll"]]) {
      await p.selectOption(`#net-model-order [data-mo-kind="${kind}"] .mo-level`, level).catch(() => {});
      await wait(900);
    }
    await wait(1800);
    await caption("Model order", "The state count falls with it", "Same network, same operating point, same disturbance — a fraction of the states");
    await hover("#mo-summary").catch(() => {});
    await wait(3000);
    await click("a.nav-link[data-page='emt']", 700);
    await p.waitForFunction(() => !document.querySelector("#emt-perturb")?.disabled, null, { timeout: 300000 });
    await wait(1200);
    await scrollToEl("#emt-model-class", 120).catch(() => {});
    await caption("Time domain", "The page now reads RMS", "Same network, same operating point, same 0.15 pu step");
    await wait(2600);
    await scrollTo(0, 700);
    await liveOn();
    await runLive("Integrating — reduced RMS model", "The same run, without the electromagnetic transients");
    await scrollToEl("#emt-results .card:nth-child(2)", 120).catch(() => {});
    await caption("Time domain", "The electromechanical answer is the same", "What has gone is the fast transient in the first cycles — which is exactly what the reduction claims to remove, and what the adequacy check quantifies");
    await wait(4800);
    // Back to the full model for the remaining scenes. The caption goes
    // first: what it says is about the plot that is on screen now.
    await p.evaluate(() => window.__tour.hideCaption());
    await wait(400);
    await p.goto(BASE + "#/network"); await wait(1000);
    await p.selectOption(netLevel, "full").catch(() => {});
    for (const kind of ["sm", "gfm", "gfl"]) {
      await p.selectOption(`#net-model-order [data-mo-kind="${kind}"] .mo-level`, "full").catch(() => {});
      await wait(400);
    }
    await wait(800);
  });

  await scene("emt-islanding", async () => {
    await click("a.nav-link[data-page='emt']", 700);
    await p.waitForFunction(() => !document.querySelector("#emt-perturb")?.disabled, null, { timeout: 300000 });
    await scrollTo(0, 800);
    await caption("Time domain", "A disturbance can also be a network event", "A breaker opening, a load step, or a phase jump — applied at T0");
    await click("[data-dist='event']", 900).catch(() => {});
    await hover("#emt-ev-kind").catch(() => {});
    await wait(2000);
    await caption("Time domain", "Trip the infinite bus itself", "The slack is a unit like any other now — the MV network is left islanded on its own machines");
    await p.selectOption("#emt-ev-breaker", await p.evaluate(() => {
      const ib = state.network.der_units.find(d => d.unit_type === "infinite_bus");
      return `unit:${ib.id}`;
    })).catch(() => {});
    await wait(2000);
    await setScopes(false);
    await liveOn();
    await runLive("Integrating the post-event network", "Same components, same setpoints — only the topology changed");
    await scrollToEl("#emt-results .card:nth-child(2)", 120).catch(() => {});
    await caption("Time domain", "The machine and the grid-former pick up the infeed", "The infinite bus's own trace stops at T0; the others take what it was carrying");
    await sweepChart("#emt-results .chart svg", 2400).catch(() => {});
    await wait(2400);
    await scrollToEl("#emt-results .card:nth-child(4)", 110).catch(() => {});
    await caption("Time domain", "And the island settles at its own frequency", "Nothing holds it at 50 Hz any more but the units left in it");
    await wait(3600);
  });

  await scene("emt-load", async () => {
    await scrollTo(0, 800);
    await caption("Time domain", "One more: a load disconnecting", "The network as loaded, and the load at bus 1 dropped at T0");
    await p.selectOption("#emt-ev-breaker", "load:0").catch(() => {});
    await wait(1800);
    await liveOn();
    await runLive("Integrating", "");
    await scrollToEl("#emt-results .card:nth-child(2)", 120).catch(() => {});
    await caption("Time domain", "The infeed reverses", "With its load gone, that bus sends power back up to the grid — the converters barely move");
    await sweepChart("#emt-results .chart svg", 2400).catch(() => {});
    await wait(2600);
    await scrollTo(0, 800);
    await caption("Time domain", "Also available", "A load step by a percentage, a phase jump, the linearised overlay for comparison, solver choice and fixed or variable step");
    await wait(3000);
  });

  // ============================================================ DOCS
  await chapter("7 · Documentation", "The manual, built in");

  await scene("docs", async () => {
    await p.goto(BASE + "#/docs"); await wait(2600);
    await caption("Documentation", "Every documentation page, in the navigation panel", "Getting started, the API, and a page per module");
    await wait(2800);
    await scrollTo(700, 1400);
    await caption("Documentation", "Equations, methodology and worked examples", "Component models, the reference frame, operating point, modal analysis, model order, EMT");
    await wait(3000);
  });

  await p.evaluate(() => window.__tour.hideCaption());
  await p.evaluate(h => window.__tour.card(h), `${LOGO}<h1>G2ELin</h1><p>Open-access power system linearization and EMT simulation</p><div class="m">github.com/FKELADA/G2ELin</div>`);
  await wait(4000);

  await ctx.close();
  await browser.close();
  fs.writeFileSync(path.join(__dirname, "segments.json"), JSON.stringify({ total: now(), fast }, null, 2));
  log(`done — ${fast.length} live stretches to speed up`);
})();
