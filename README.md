# G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool

G2ELin is an open-access tool for power system linearization and Electromagnetic Transient (EMT) simulation, developed to help researchers analyze small-signal stability and dynamics in systems with a high penetration of Inverter-Based Resources (IBRs).

This repository contains two implementations:

- **[`matlab/`](matlab/README.md)** — the original MATLAB/Simulink implementation (model library + linearization/modal analysis toolbox).
- **[`python/`](python/README.md)** — a Python port of the same functionality.

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
