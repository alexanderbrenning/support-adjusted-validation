# Reproduce all analysis outputs for the support-adjusted validation manuscript.
#
# Run from the RStudio project root:
#   source("R/run_all_outputs.R")
#
# Output files are regenerated in results/examples/, results/simulation/ and
# results/supplementary/. Interactive Shiny calculators are launched separately
# because shiny::runApp() blocks the R session.

VERBOSE <- 2

# Post-hoc example calculations and combined example tables.
source("R/examples/run_all_examples.R")

# Simulation study and simulation figures.
source("R/simulation/run_simulation.R")
source("R/simulation/summarise_simulation.R")
source("R/simulation/plot_simulation.R")

# Manuscript and supplementary tables/figures derived from the outputs above.
source("R/paper/run_paper_outputs.R")

# Interactive calculators, run separately when needed:
# source("R/apps/launch_tier1_reference_error_app.R")
# source("R/apps/launch_tier2_support_adjustment_app.R")

# Standalone HTML/JavaScript calculators, open separately when needed:
# utils::browseURL(normalizePath("html/tier1_reference_error_app/index.html"))
# utils::browseURL(normalizePath("html/tier2_support_adjustment_app/index.html"))
