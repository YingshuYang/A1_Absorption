# Multi-Layer Optical Laser Absorption Analysis

## Overview
`A1_Absorption` is a computational tool designed to analyze the absorption of laser light within multilayer systems using the Transfer Matrix Method (TMM) and the Poynting theorem. This software enables the calculation of reflection, transmission, and total absorption rates for multilayer structures. Additionally, it offers detailed insights into the absorption of each individual layer within the system.

## Features
- **Reflection, Transmission, and Absorption Calculation**: Quantifies the basic optical properties of multilayer systems.
- **Layer-specific Absorption Analysis**: Provides detailed absorption data for each layer, enhancing the understanding of laser-material interaction.

## Physical Background
The underlying physics and mathematical derivations implemented in this code are extensively discussed in Yingshu Yang's thesis:

- Yang, Yingshu. "Modeling terahertz wave propagation and production in out-of-equilibrium complex heterostructure devices." (2023).

## Citation
If you utilize `A1_Absorption` for your research, please cite the following publications to acknowledge the theoretical foundation and the software's contribution to your work:

1. Yang, Yingshu, Stefano Dal Forno, and Marco Battiato. "Removal of Spectral Distortion Due to Echo for Ultrashort THz Pulses Propagating Through Multilayer Structures with Thick Substrate." Journal of Infrared, Millimeter, and Terahertz Waves (2021): 1-11.
2. Yang, Yingshu, Stefano Dal Forno, and Marco Battiato. "Modeling spintronic terahertz emitters as a function of spin generation and diffusion geometry." Physical Review B 107.14 (2023): 144407.
3. Yang, Yingshu, Stefano Dal Forno, and Marco Battiato. "Transfer-matrix description of heterostructured spintronics terahertz emitters." Physical Review B 104.15 (2021): 155437.

## Instructions

`A1_Absorption` provides two main functionalities to cater to different research needs:

### Single Sample Analysis
- **Script**: `absorption_calculation_d_fixed_Main`
- **Purpose**: Calculates the absorption of the laser in a single sample with fixed layer thicknesses.
- **Outputs**: Prints the transmission, reflection, total absorption, and layer-specific absorption ratios.

### Multiple Samples Analysis
- **Script**: `absorption_calculation_d_change_Main`
- **Purpose**: Analyzes multiple samples by varying the thickness of one specific layer across samples.
- **Outputs**: Besides printing the transmission, reflection, total absorption, and layer-specific absorption ratios, it also plots the changes with respect to the varying thicknesses.

### Preparing for Calculations
- **Inputs**: Before running calculations, please refer to the comments within each main file for guidance on the required inputs.

## Future Developments
- We are committed to enhancing the usability of `A1_Absorption` and are currently working on developing a Graphical User Interface (GUI) to make the tool more accessible and user-friendly. Although the GUI is still under development, we aim to provide an intuitive interface that will simplify the process of configuring and running analyses. Stay tuned for updates on this exciting feature!
- We are now providing only single frequency calculation. We will update the files to multiple frequency calculation in the future for more accurate results.
