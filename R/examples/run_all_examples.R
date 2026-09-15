# Run and combine post-hoc examples
#
# Each script writes one example-level CSV to results/examples/.

OVERWRITE <- TRUE
VERBOSE <- 2

scripts <- c(
  "R/examples/example_biomass.R",
  "R/examples/example_smap_crns.R",
  "R/examples/example_bts.R",
  "R/examples/example_swe.R",
  "R/examples/example_no2_de.R",
  "R/examples/example_meuse.R",
  "R/examples/example_s1_sm.R",
  "R/examples/example_no2_eu.R"
)


for (s in scripts) {
  if (VERBOSE >= 1) message("\n--- Running ", s, " ---")
  source(s, local = new.env(parent = globalenv()))
}

source("R/examples/examples_functions.R")
combine_posthoc_examples(overwrite = TRUE, verbose = VERBOSE)


# Generate compact manuscript-facing table and quantitative summary figure.
source("R/paper/make_posthoc_summary_figure_and_table.R")
source("R/paper/make_suppl_posthoc_example_table.R")

