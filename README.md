# G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool

G2ELin is an open-access tool for power system linearization and time-domain simulation, developed to help researchers analyze small-signal stability and dynamics in systems with a high penetration of Inverter-Based Resources (IBRs).

Every element is modelled in full detail by default, which makes a run an Electromagnetic Transient (EMT) simulation. The model order of each element is also selectable — the network can be made quasi-stationary and the machines and converters lower-order, which turns the same case into an electromechanical (RMS) study for large networks. A built-in adequacy check compares the reduced model against the full one on *your* network and says whether the reduction is safe for it.

This repository contains two implementations:

- **[`matlab/`](matlab/README.md)** — the original MATLAB/Simulink implementation (model library + linearization/modal analysis toolbox).
- **[`python/`](python/README.md)** — a Python port of the same functionality.

## Documentation

**<https://fkelada-g2elin-docs.static.hf.space>** — the full manual (architecture,
equations, API reference, worked examples), in the Read the Docs layout.

## Try it

**[GETTING_STARTED.md](GETTING_STARTED.md) — three steps, no coding.** Install
Python, download this repository, then double-click `start-windows.bat`
(Windows) or run `./start-macos-linux.sh` (macOS/Linux): it sets everything
up in a private environment inside the folder and opens the web interface at
<http://127.0.0.1:8000>.

From a command line instead:

```bash
cd python
python -m venv .venv
source .venv/bin/activate                # Windows instead: .venv\Scripts\activate
pip install -e ".[api,docs]"
python tools/build_docs.py               # the in-app manual (optional)
python -m uvicorn g2elin_api.main:app --port 8000
```

[`python/docs/sphinx/installation.md`](python/docs/sphinx/installation.md)
has the full setup: the optional extras (`api`, `docs`, `dev`, `notebook`),
the test suite, the notebooks, troubleshooting and server deployment.

See each directory's README for requirements, installation, and usage instructions specific to that implementation.

## Copyright and License
This tool is developed by Fadi Kelada at G2ELab. It is provided under the GNU GPL-3 license.

## Usage and Warranty
G2ELin is provided as an open-access tool without any guarantees or warranties. Users are encouraged to explore, experiment, and contribute to its development. However, the developers do not assume any liability for issues arising from its use.

## Acknowledgements
I would like to extend my gratitude to Arshpreet Singh for providing the initial code for performing the linearization analysis, which served as the foundation for this tool. Additionally, I am grateful to the SYREL task force, a group of wonderful colleagues working on closely related research topics. Our discussions and meetings were instrumental in enhancing my understanding of many important aspects of my research and generating numerous interesting ideas. Special thanks to Sameh Betamony, Heitor Farias, and Marta Gomis for being part of this working group. Your insights were incredibly helpful and significantly contributed to advancing my research.

Additionally, I am grateful for the following citations, which motivated me to develop this tool and provided a foundational base:
- T. Qoria, T. Jouini, D. Gross, U. Markovic, G. Denis, and T. Prevost, "Data underlying the research of a 3 bus model for full inverter system - Migrate WP3." 4TU.Centre for Research Data, 2018. [Online]. Available: https://doi.org/10.4121/uuid:e5497fd2-f617-4573-b6d5-1202ebae411d.
- A. Tayyebi, D. Groß, and A. Anta. (2019). GridFormingConverters: Implementation of Grid-Forming Control Techniques in IEEE 9- Bus System. Git Repository. [Online]. Available: https://github.com/ATayebi/GridFormingConverters

## How to Cite
If you use G2ELin in your research, please cite it as follows:
F. Kelada, "G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool," 2023. https://github.com/FKELADA/G2ELin
