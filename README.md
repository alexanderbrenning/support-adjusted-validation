# Support-adjusted validation metrics: code and data

This repository contains code and selected data for:

> Brenning, A. & Suesse, T. (2026). *Support-adjusted validation metrics in environmental prediction: a minimal-information framework*. Submitted manuscript.

The code reproduces the simulation study, post-hoc examples, interactive calculators, and manuscript/supplementary tables and figures. The repository is intended to be opened as an RStudio project from the `code_data` root. All paths are relative to this root.

## Folder structure

```text
R/core/          core support-adjustment functions
R/simulation/    simulation design, execution, summaries, and simulation figures
R/examples/      post-hoc example calculations and shared example functions
R/paper/         manuscript and supplementary table/figure generation
R/ancillary/     auxiliary checks and data-preparation utilities
R/apps/          launchers for the Shiny calculators
shiny/           Tier 1 and Tier 2/3 Shiny applications
html/            standalone browser versions of both interactive calculators
data/            redistributed input data where permitted or practically useful
results/         generated outputs; ignored/empty in the release except for .gitkeep files
```

## Main workflows

Complete reproduction of all released analysis outputs:

```r
source("R/run_all_outputs.R")
```

The main runner recalculates examples, simulation outputs, simulation figures, and paper/supplementary tables and figures.

Equivalent step-by-step workflow:

```r
source("R/examples/run_all_examples.R")
source("R/simulation/run_simulation.R")
source("R/simulation/summarise_simulation.R")
source("R/simulation/plot_simulation.R")
source("R/paper/run_paper_outputs.R")
```

The first four lines generate the example outputs and simulation figures. The final paper-output script generates the remaining manuscript/supplementary outputs, including adjustment-function figures and the simulation synopsis table.

Launch the interactive calculators directly:

```r
source("R/apps/launch_tier1_reference_error_app.R")
source("R/apps/launch_tier2_support_adjustment_app.R")
```

The Tier 1 app implements the reference-error adjustment. The Tier 2/3 app implements block-support adjustment either from direct `lambda_A,B` and `Delta_A,B` inputs or from a semivariogram support approximation. Its single example-input panel populates the validation fields and selects the relevant support-input mode.

Equivalent standalone browser applications are provided at:

```text
html/tier1_reference_error_app/index.html
html/tier2_support_adjustment_app/index.html
```

Open either HTML file directly in a modern web browser. Each is a self-contained, dependency-free page: all calculations, input synchronization, diagnostics, example settings, styling, and (for Tier 2/3) semivariogram plotting are implemented in embedded JavaScript and require neither R nor an internet connection.

## Package requirements

The scripts load the packages they use. Main non-base packages are `sf`, `sp`, `gstat`, `randomForest`, `ggplot2`, `ggrepel`, `terra` where needed, and `shiny` for the Shiny apps. Install missing packages in the usual way before running the corresponding workflow. The standalone HTML applications have no package requirements.

## Outputs

Examples write one CSV per example and a combined result table to `results/examples/`. Main-paper and supplementary figures/tables are written to `results/paper/` and `results/supplementary/` by the scripts in `R/paper/`. Supplementary output files use `suppl_` in their names.

## Data provenance and licences

See `data/DATA_PROVENANCE.md`. Code is released under GPL-3.0-or-later. Data files retain the rights/licences of their original providers and are redistributed only for reproducibility where this appears permissible or where the file is derived from open data. In particular, the original Vizcaino and Lavalle NO2--EU dataset is distributed through Mendeley Data under CC BY 4.0; the coordinate-augmented derivatives also incorporate historical AirBase/EEA metadata.
