# G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool

## Introduction
G2ELin is an open-access MATLAB-based tool developed to address the need for comprehensive power system linearization and Electromagnetic Transient (EMT) simulations. This tool was developed to assist researchers in the field of power systems in understanding the development and application of small-signal models. It aims to facilitate the analysis of power system stability and dynamics, particularly in scenarios with a high penetration of Inverter-Based Resources (IBRs).

The tool is composed of two main components:

1. **Model Library (_G2ELib_)**: This library features comprehensive full-order Electromagnetic Transient (EMT) models of synchronous generators and Voltage Source Converters (VSC) in both Grid-Following (GFL) and Grid-Forming (GFM) modes. It includes various advanced GFM controllers, such as conventional droop-based controllers, Synchronverter, dispatchable Virtual Oscillator Controllers, and Matching Controller, enabling comparative analysis. These models are implemented in MATLAB Simulink (version 2020b).

2. **Linearization and Modal Analysis Toolbox (_G2ELin_)**: This toolbox linearizes the Simulink-simulated topology and provides a corresponding modal analysis summary.

## Requirements
To fully utilize G2ELin and all its features, you will need the following software licenses:
- MATLAB 2020b
- MATLAB Simulink 2020b
- Symbolic Math Toolbox

## Download and Installation
To download and run G2ELin, follow these steps:
1. Clone the repository from GitHub: git clone https://github.com/FKELADA/G2ELin
2. Open MATLAB and ensure you have the required licenses (MATLAB 2020b, Simulink 2020b, Symbolic Math Toolbox).
3. Run the "startup.m" file.
4. Further detailed instructions on how to use the tool are provided in the documentation file included in the repository. This documentation covers the full description of the tool, its components, and usage instructions, along with plans for future development.

## Copyright and License
This tool is developed by Fadi Kelada at G2ELab. It is provided under the GNU GPL-3 license. 

## Usage and Warranty
G2ELin is provided as an open-access tool without any guarantees or warranties. Users are encouraged to explore, experiment, and contribute to its development. However, the developers do not assume any liability for issues arising from its use.

## Acknowledgements
I would like to extend my gratitude to Arshpreet Singh for providing the initial code for performing the linearization analysis, which served as the foundation for this tool. Additionally, I am grateful to the SYREL task force, a group of wonderful colleagues working on closely related research topics. Our discussions and meetings were instrumental in enhancing my understanding of many important aspects of my research and generating numerous interesting ideas. Special thanks to Sameh Betamony, Heitor Farias, and Marta Gomis for being part of this working group. Your insights were incredibly helpful and significantly contributed to advancing my research.

Additionally, I am grateful for the following citations, which motivated me to develop this tool and provided a foundational base:
- T. Qoria, T. Jouini, D. Gross, U. Markovic, G. Denis, and T. Prevost, “Data underlying the research of a 3 bus model for full inverter system - Migrate WP3.” 4TU.Centre for Research Data, 2018. [Online]. Available: https://doi.org/10.4121/uuid:e5497fd2-f617-4573-b6d5-1202ebae411d.
- A. Tayyebi, D. Groß, and A. Anta. (2019). GridFormingConverters: Implementation of Grid-Forming Control Techniques in IEEE 9- Bus System. Git Repository. [Online]. Available: https://github.com/ATayebi/GridFormingConverters

## How to Cite
If you use G2ELin in your research, please cite it as follows: 
F. Kelada, “G2ELin: An Open-Access Power System Linearization and EMT Simulation Tool,” 2023. https://github.com/FKELADA/G2ELin
