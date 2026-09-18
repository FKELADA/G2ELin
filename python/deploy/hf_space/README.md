---
title: G2ELin
emoji: ⚡
colorFrom: blue
colorTo: gray
sdk: docker
app_port: 7860
pinned: false
license: gpl-3.0
short_description: Power-system linearisation, modal analysis and EMT
---

# G2ELin — web app

Small-signal and EMT analysis for inverter-dominated power systems:
build or load a network, solve its power flow, linearise it for modal
analysis, and check the result against a nonlinear time-domain (EMT)
simulation of the same equations.

- Source code: <https://github.com/FKELADA/G2ELin>
- Developed by Fadi Kelada at G2ELab, under the GNU GPL-3 licence, without warranty.

This Space runs the `python/` implementation (FastAPI backend + web UI +
documentation). Networks are held in your browser only — nothing you edit
is stored on the server. The Space sleeps when idle, so the first request
after a pause takes a little longer while it restarts.

Cite: F. Kelada, "G2ELin: An Open-Access Power System Linearization and EMT
Simulation Tool," 2023.
