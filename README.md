# Salish Sea Mizer Ecosystem Model
# Created by: Viv Tulloch

**Supporting code for:**
*Tulloch et al. (in review). Multi-species ecosystem modelling to support conservation of Southern Resident Killer Whales and Pacific salmon. Ecological Applications.*

---

## Overview

This repository contains the code and workflows used to build, calibrate, and apply a **multi-species size-spectrum ecosystem model** of the Salish Sea using the [`mizer`](https://sizespectrum.org/mizer) framework. The model integrates **25 functional groups** (salmon, marine mammals, groundfish, forage fish, invertebrates) and is used to explore ecosystem responses to fishing, climate variability, and conservation strategies.

The workflow has three main components:

1. **Calibration and tuning (2011 baseline)**
   * Iterative adjustment of species parameters (feeding, reproduction, recruitment).
   * Ensures realistic biomass, diet, growth, and resilience.
  
2. **Historic runs (hindcast, 1986–2011)**
   * Spin-up and hindcast simulations with reconstructed fishing mortality.
   * Validation against survey and catch data.

3. **Scenario simulations (2011–2111)**
   * Forward projections under alternative management scenarios (e.g., reduced salmon fishing, pinniped changes, productivity shifts).
   * Trade-off analysis for Southern Resident Killer Whale (SRKW) recovery and salmon conservation.

---

## File Structure

* **Tuning Runs**
  * `ss_mizer_2011_250515.Rmd`: Interactive and manual parameter calibration.
  * Includes biomass calibration, yield curve analysis, and resilience adjustments.

* **Historic Runs**
  * `ss_mizer_2011_250813_runhistoric.Rmd`
  * Code to reconstruct fishing effort, hindcast 1986–2011, and validate against data.

* **Simulation Runs**
  * `Simulations.R`: Test script for scenario projections.
  * Includes functions for:
    * Model tuning (marine mammals, predators, ecosystem parameters).
    * Historical simulations and survey validation.
    * Future effort scenarios (status quo, reduced salmon, etc.).
    * Reference points (SSB at 50% unexploited).
    * Scenario comparisons and visualizations.
    * Export of outputs (biomass, abundance, yield, large fish indicator, etc.).

---

## Inputs

* **Parameter files:**

  * `ss_pars_2011_241212_2.csv` (species parameters)
  * `ss_gear_2011_241212.csv` (gear parameters)
  * `ss_mizer_interaction_241212.csv` (predator–prey interactions)

* **Effort data:**

  * `All_fisheries_effort_array_250521.csv` (time series of relative fishing effort, 1986–2011).

* **Survey data:**

  * Marine mammal and salmon biomass indices for validation.

---

## Dependencies

* R ≥ 4.0
* [`mizer`](https://github.com/sizespectrum/mizer) + extensions:

  * [`mizerExperimental`](https://github.com/sizespectrum/mizerExperimental)
  * [`mizerMR`](https://github.com/sizespectrum/mizerMR)
* `tidyverse`
* `plotly`
* `ggplot2`, `reshape2`, `abind`

Install with:

```r
remotes::install_github("sizespectrum/mizerExperimental")
remotes::install_github("sizespectrum/mizerMR")
install.packages(c("tidyverse", "ggplot2", "plotly", "reshape2", "abind"))
```

---

## Workflow

### 1. Historic Runs (1986–2011)

* Spin-up to equilibrium with no fishing.
* Apply reconstructed fishing effort by fleet/gear.
* Validate against observed biomass and catch data.

### 2. Calibration & Tuning (2011 baseline)

* Adjust feeding (`gamma`, `h`), reproduction (`erepro`, `R_max`), and recruitment dynamics.
* Use `tuneParams()` interface and yield curves to test resilience.
* Ensure realistic predator–prey dynamics and species coexistence.

### 3. Scenario Simulations (2011–2111)

* Construct future effort scenarios:

  * **Status quo**: Continue 2011 fishing levels.
  * **Reduced salmon fishing**: 50% reduction in salmon sectors, small increases in others, elimination of bycatch.
  * Additional climate/productivity scenarios (see manuscript).
* Project biomass, yield, and ecosystem indicators.
* Compare scenarios relative to reference points (50% unexploited SSB).

---

## Outputs

* **Model parameter sets** (`.rds`, `.csv`) at different stages (baseline, tuned, scenario-specific).
* **Time series outputs** (CSV):

  * Biomass, spawning stock biomass (SSB)
  * Yield and catch
  * Predation mortality (M2)
  * Large Fish Indicator (LFI)
* **Figures and animations**:

  * Spectra plots
  * Biomass trajectories
  * Scenario comparison plots
  * Validation vs survey data

---

## Notes

* File paths in scripts point to local directories (`/Dropbox/UBC/...`). Update paths before running externally.
* Effort matrices must be normalized relative to 2011 values before projecting scenarios.
* Species indices in parameter tuning (e.g., `19:21` = marine mammals) depend on the species ordering in your parameter files.

---

## Citation

If using this code or model, please cite:

Tulloch, V., Morzaria-Luna, H., Murray, C., & Martin, T. (in review).
*Multi-species ecosystem modelling to support conservation of Southern Resident Killer Whales and Pacific salmon.* Ecological Applications.

---

