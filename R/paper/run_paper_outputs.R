# Generate manuscript and supplementary tables/figures that are not produced
# directly by the example and simulation workflows.
#
# Run after:
#   source("R/examples/run_all_examples.R")
#   source("R/simulation/run_simulation.R")
#   source("R/simulation/summarise_simulation.R")
#   source("R/simulation/plot_simulation.R")

VERBOSE <- 2

source("R/paper/make_posthoc_summary_figure_and_table.R")
source("R/paper/make_suppl_posthoc_example_table.R")
source("R/paper/make_suppl_correction_function_plots.R")
source("R/paper/extract_suppl_simulation_synopsis.R")
