
(() => {
  if (window.__tour) return;
  const css = document.createElement("style");
  css.textContent = `
    #tour-cursor { position: fixed; z-index: 99999; width: 22px; height: 22px; margin: -3px 0 0 -3px; pointer-events: none;
      transition: transform 0.08s; filter: drop-shadow(0 2px 3px rgba(0,0,0,0.35)); }
    .tour-ripple { position: fixed; z-index: 99998; width: 36px; height: 36px; margin: -18px 0 0 -18px; border-radius: 50%;
      border: 3px solid #2a78d6; pointer-events: none; animation: tourRipple 0.55s ease-out forwards; }
    @keyframes tourRipple { from { transform: scale(0.3); opacity: 1; } to { transform: scale(1.5); opacity: 0; } }
    #tour-caption { position: fixed; z-index: 99997; left: 50%; bottom: 28px; transform: translateX(-50%) translateY(20px); opacity: 0;
      background: rgba(18,22,28,0.94); color: #fff; border-radius: 14px; padding: 14px 26px; max-width: 1040px; text-align: center;
      font-family: "IBM Plex Sans", system-ui, sans-serif; box-shadow: 0 10px 40px rgba(0,0,0,0.35); transition: all 0.35s ease; pointer-events: none; }
    #tour-caption.on { opacity: 1; transform: translateX(-50%) translateY(0); }
    #tour-caption .t { font-size: 21px; font-weight: 700; letter-spacing: -0.01em; }
    #tour-caption .s { font-size: 15px; color: #c9ced6; margin-top: 3px; }
    #tour-caption .k { font-family: "IBM Plex Mono", monospace; font-size: 11px; letter-spacing: 0.12em; text-transform: uppercase; color: #8fb8ea; margin-bottom: 4px; }
    #tour-card { position: fixed; inset: 0; z-index: 99996; display: flex; flex-direction: column; align-items: center; justify-content: center;
      background: radial-gradient(900px 400px at 70% 10%, #24476f 0%, transparent 60%), linear-gradient(135deg, #12161c, #1b2330);
      color: #fff; font-family: "IBM Plex Sans", system-ui, sans-serif; opacity: 0; transition: opacity 0.6s; pointer-events: none; }
    #tour-card.on { opacity: 1; }
    #tour-card .logo { width: 84px; height: 84px; margin-bottom: 26px; }
    #tour-card h1 { font-size: 64px; margin: 0; letter-spacing: -0.03em; }
    #tour-card h2 { font-size: 44px; margin: 0; letter-spacing: -0.02em; }
    #tour-card p { font-size: 24px; color: #c9ced6; margin: 14px 0 0; }
    #tour-card .m { font-family: "IBM Plex Mono", monospace; font-size: 14px; letter-spacing: 0.14em; color: #8fb8ea; text-transform: uppercase; margin-top: 30px; }
  `;
  document.head.appendChild(css);
  const cur = document.createElement("div");
  cur.id = "tour-cursor";
  cur.innerHTML = '<svg viewBox="0 0 24 24" width="22" height="22"><path d="M3 2l7 19 2.6-7.4L20 11z" fill="#fff" stroke="#111" stroke-width="1.6" stroke-linejoin="round"/></svg>';
  cur.style.left = "800px"; cur.style.top = "450px";
  const cap = document.createElement("div"); cap.id = "tour-caption";
  const card = document.createElement("div"); card.id = "tour-card";
  document.body.append(cur, cap, card);
  window.addEventListener("mousemove", e => { cur.style.left = e.clientX + "px"; cur.style.top = e.clientY + "px"; }, true);
  window.addEventListener("pointermove", e => { cur.style.left = e.clientX + "px"; cur.style.top = e.clientY + "px"; }, true);
  window.addEventListener("mousedown", e => {
    const r = document.createElement("div"); r.className = "tour-ripple";
    r.style.left = e.clientX + "px"; r.style.top = e.clientY + "px";
    document.body.appendChild(r); setTimeout(() => r.remove(), 600);
    cur.style.transform = "scale(0.85)"; setTimeout(() => { cur.style.transform = ""; }, 120);
  }, true);
  window.__tour = {
    caption(kicker, title, sub) {
      cap.classList.remove("on");
      setTimeout(() => {
        cap.innerHTML = (kicker ? '<div class="k">' + kicker + '</div>' : "") + '<div class="t">' + title + '</div>' + (sub ? '<div class="s">' + sub + '</div>' : "");
        cap.classList.add("on");
      }, title ? 180 : 0);
    },
    hideCaption() { cap.classList.remove("on"); },
    card(html) { card.innerHTML = html; card.classList.add("on"); },
    hideCard() { card.classList.remove("on"); },
  };
})();