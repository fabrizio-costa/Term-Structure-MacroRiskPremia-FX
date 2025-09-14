# Term Structure of Macroeconomic Risk Premia in Foreign Exchange Markets: A Bayesian Approach

This repository contains the code and data for my Master’s thesis:  
**“Term Structure of Macroeconomic Risk Premia in Foreign Exchange Markets: A Bayesian Approach.”**

The project is organized as follows:

## Repository Structure

### Root Directory
- **`Term-Structure-MacroRiskPremia-FX.Rproj`** – RStudio project file for reproducibility and environment management.  
- **`code/`** – all R scripts used for simulation, posterior sampling, and Gibbs samplers.  
- **`data/`** – datasets used in the analysis (currency and macroeconomic data).  

---

### Code Folder (`code/`)

- **`main.R`**  
  Entry point script. Running this file executes the full pipeline (data loading, estimation, and results generation).  

- **`matrix.notation.R`**  
  Contains functions to construct the matrices and vectors used in the model, matching the notation in the thesis (e.g. \( V_\rho, V_\eta, G, \bar{G}, V_r, B_r \)).  

- **`posterior.R`**  
  Implements conditional posterior distributions for all model parameters. Each function returns draws from the corresponding posterior, as derived in the thesis.  

- **`gibbs.sampler.R`**  
  Combines the posterior functions and matrix definitions into structured Gibbs samplers. Includes separate samplers for each block and a full sampler.  

- **`helpers.R`**  
  Utility functions used across scripts to keep the main code clean (e.g. diagnostics, stability checks).  

- **`checks/`**  
  Contains validation scripts and simulation code:  
  - **`1check.gibbs.sampler.R`**, **`2check.gibbs.sampler.R`**, **`3check.gibbs.sampler.R`**, **`4check.single.posterior.R`** – scripts to validate each Gibbs sampler or posterior individually (trace plots, density plots, etc.).  
  - **`simulation.data.R`** – script to simulate synthetic datasets used for testing the samplers.  
  - **`plot/`** – folder with PDF plots produced during the validation checks (trace plots, posterior densities, ACFs, etc.).  

---

### Data Folder (`data/`)

- **`currency_data/`** – datasets related to exchange rates and FX premia:  
  - `all_outright.rds`  
  - `df_long.RDS`, `df_short.RDS`  
  - `ds_data.xlsx`  

- **`macroeconomic_data/`** – datasets with macroeconomic series used in the thesis.  

---

## How to Use

1. Open the project in RStudio by clicking `Term-Structure-MacroRiskPremia-FX.Rproj`.  
2. Ensure required packages are installed (`MASS`, `MCMCpack`, `tidyverse`, etc.).  
3. Run `code/main.R` to reproduce the main results.  
4. Validation and diagnostic scripts can be run separately from `code/checks/`.  

---
