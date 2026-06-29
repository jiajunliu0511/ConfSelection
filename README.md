# ConfSelection: Conformal Selection for Treatment Beneficiaries using Auxiliary External Data
This repository contains the R simulation code to support the methodology presented in "A Conformal Selection Framework for Individual Treatment Beneficiaries with Auxiliary External Data". The goal is to identify candidate patients who are likely to benefit from treatment while accounting for uncertainty in estimated individual treatment effects and False discovery rate (FDR) control.

## Repository Structure
The simulation pipeline is organized into data generation and statistical inference modules.
### 1. Data Generation
* `0_dataGen.R`: The main script to generate the synthetic datasets across various simulation settings explored in the study.
* `0_functions_DGP.R`: Contains the underlying functions for the Data Generating Process (DGP).
### 2. Simulation Modules
Both simulation scripts utilize helper functions defined in `0_functions_sim.R`.
* `1_sim_ConfP.R`: Performs the simulation for the **primary** analysis (Conformal Prediction).
* `1_sim_StochOrder.R`: Performs the simulation to check the **stochastic ordering condition**.
* `0_functions_sim.R`: Shared functions for executing simulation iterations, performance metrics, and data processing.
