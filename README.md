# Blob Dynamics simulation - 3rd year project Physics Imperial College London

## Overview
This project explores turbulence and filament dynamics in the scrape-off layer (SOL) of a Tokamak using numerical simulations. The developed solver is based on drift-reduced plasma fluid equations and aims to model the evolution of turbulent filaments (or blobs) in a simplified 2D geometry. 

## Features
- 2D simulation of plasma blob dynamics
- Uses drift-reduced plasma fluid equations
- Implements E × B advection, magnetic curvature effects, and dissipative forces
- Numerical solver based on iterative methods and Runge-Kutta integration
- Visualization of results with ImageMagick-generated GIFs

## Requirements
To run the simulation, install the following dependencies:
- [CERN ROOT](https://root.cern/install/) (for executing the program)
- [ImageMagick](https://imagemagick.org/script/download.php) (for GIF creation)

## Installation
1. Install CERN ROOT:
   ```sh
   sudo apt-get install root-system
   ```
2. Install ImageMagick:
   ```sh
   sudo apt-get install imagemagick
   ```

## Usage
To execute the simulation, run:
```sh
root main.cpp
```
This will generate the necessary output data and visualizations.

## Simulation details
The solver simulates the evolution of a Gaussian filament in a Tokamak SOL, considering:
- Radial propagation driven by E × B advection
- Dissipative effects through perpendicular viscosity and diffusion
- Electrostatic potential interactions
- Filament deformation under magnetic curvature effects

## Results
Our simulations validate theoretical predictions of filament dynamics. The observed results align with expectations from plasma physics, reinforcing the use of sheath dissipation closures for SOL studies. Visualizations of the simulation results are provided as GIFs.

## Future work
Potential extensions include:
- 3D simulations incorporating parallel transport effects
- Inclusion of kinetic plasma effects
- Validation with experimental data from fusion devices

## Authors
- **Giancarlo Venturato**
- **Andrea Sanfilippo**
- **Supervisor:** Dr. Robert Kingham
- **Assessor:** Dr. Yasmin Andrew

## References
For a detailed explanation of the model and results, refer to our [project report]([https://github.com/ginca03/Blob-dynamics-simulation](https://github.com/ginca03/Blob-dynamics-simulation/blob/main/docs/project_report%20Y3%2C%20Giancarlo%20Venturato%2006033934%2C%20Andrea%20Sanfilippo%2006034359%2C%20Supervisor%20Dr%20Robert%20Kingham%2C%20Assessor%20Dr%20Yasmin%20Andrew%2C%20Exploring%20turbulence%20and%20filaments%20in%20Tokamak%20exhaust%20via%20simulations.pdf)).

---

This project contributes to understanding plasma behavior in Tokamak exhaust regions, aiding in the advancement of fusion energy research.

