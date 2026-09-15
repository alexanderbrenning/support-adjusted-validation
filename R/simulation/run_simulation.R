# Run pedagogical support-adjustment simulation incrementally
#
# Existing scenario files are resumed and completed replicate-by-replicate unless OVERWRITE = TRUE.

source("R/simulation/simulation_functions.R")

OVERWRITE <- TRUE
VERBOSE <- 1

OUT_DIR <- file.path("results", "simulation")
SCENARIO_FILE <- file.path(OUT_DIR, "simulation_scenarios.csv")

# Simulation resolution and replication settings.
N_REPS <- 200
N_BLOCKS_SIDE <- 40
FINE_PER_BLOCK <- 5
N_SAMPLE <- 500
BASE_SEED <- 20260906

#' Construct simulation scenario table
#'
#' @return data.frame of scenarios.
make_scenarios <- function() {
  core <- expand.grid(
    q_epsilon = c(0, 0.10, 0.25, 0.40),
    range_block_ratio = c(2, 10, 50),
    microscale = "moderate",
    target_r2_block = c(0.30, 0.60, 0.80),
    prediction_relation = "block_targeted",
    scenario_type = "core",
    stringsAsFactors = FALSE
  )
  # Sparse microscale-sensitivity scenarios retain the standard block-targeted
  # prediction while varying the true microscale component.
  microscale_sensitivity <- expand.grid(
    q_epsilon = c(0, 0.25),
    range_block_ratio = 2,
    microscale = c("none", "large"),
    target_r2_block = 0.60,
    prediction_relation = "block_targeted",
    scenario_type = "targeted_microscale",
    stringsAsFactors = FALSE
  )
  # Separate failure-mode scenarios deliberately violate Cov(Y_B, delta) = 0.
  failure <- expand.grid(
    q_epsilon = c(0.10, 0.25),
    range_block_ratio = c(2, 10),
    microscale = "moderate",
    target_r2_block = 0.60,
    prediction_relation = "subgrid_informed",
    scenario_type = "failure_mode",
    stringsAsFactors = FALSE
  )
  out <- unique(rbind(core, microscale_sensitivity, failure))
  out$scenario_index <- seq_len(nrow(out))
  out$scenario_id <- sprintf(
    "sc%03d_q%03d_rb%02d_micro-%s_r2%02d_%s",
    out$scenario_index,
    round(100 * out$q_epsilon),
    out$range_block_ratio,
    out$microscale,
    round(100 * out$target_r2_block),
    out$prediction_relation
  )
  out$true_nugget <- vapply(out$microscale, microscale_fraction, numeric(1))
  out$partial_sill <- 1 - out$true_nugget
  out$true_sill <- 1
  out$sigma2_eps <- ifelse(out$q_epsilon == 0, 0, out$q_epsilon / (1 - out$q_epsilon))
  out$observed_nugget <- out$true_nugget + out$sigma2_eps
  out$observed_sill <- out$true_sill + out$sigma2_eps
  out$spherical_range_block_units <- out$range_block_ratio
  out$spherical_range_fine_units <- out$range_block_ratio * FINE_PER_BLOCK

  out[, c("scenario_index", "scenario_id", "scenario_type", "q_epsilon",
          "range_block_ratio", "microscale", "true_nugget", "partial_sill",
          "true_sill", "sigma2_eps", "observed_nugget", "observed_sill",
          "spherical_range_block_units", "spherical_range_fine_units",
          "target_r2_block", "prediction_relation")]
}

#' Run or resume all repetitions for one scenario
#'
#' Stores all repetitions for the scenario in one RDS file. The file is updated
#' after each completed repetition, so interrupted runs can resume from the last
#' completed replicate while keeping the number of files small.
#'
#' @param scenario One-row scenario data.frame.
#' @param n_reps Number of repetitions.
#' @return data.frame with all repetitions for the scenario.
run_scenario_incremental <- function(scenario, n_reps) {
  out_file <- file.path(OUT_DIR, paste0(scenario$scenario_id, "_reps.rds"))

  if (file.exists(out_file) && !OVERWRITE) {
    existing <- readRDS(out_file)
    done_reps <- sort(unique(existing$rep))
    if (length(done_reps) >= n_reps && all(seq_len(n_reps) %in% done_reps)) {
      if (VERBOSE >= 1) message("  all reps exist, skipped")
      return(existing)
    }
  } else {
    existing <- NULL
    done_reps <- integer(0)
  }

  rows <- if (is.null(existing) || OVERWRITE) list() else list(existing)

  for (rr in setdiff(seq_len(n_reps), done_reps)) {
    if (VERBOSE >= 2) message("  rep ", rr, "/", n_reps, ": running")

    res <- simulate_one_replicate(
      scenario = scenario,
      rep = rr,
      n_blocks_side = N_BLOCKS_SIDE,
      fine_per_block = FINE_PER_BLOCK,
      n_sample = N_SAMPLE,
      seed = BASE_SEED
    )

    rows[[length(rows) + 1]] <- res
    combined <- do.call(rbind, rows)
    saveRDS(combined, out_file)
  }

  if (VERBOSE >= 1) message("  wrote/resumed scenario file: ", out_file)
  readRDS(out_file)
}

scenarios <- make_scenarios()
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
write_csv_incremental(scenarios, SCENARIO_FILE, overwrite = TRUE, verbose = VERBOSE)

if (VERBOSE >= 1) {
  message("Scenarios: ", nrow(scenarios))
  message("Replicates per scenario: ", N_REPS)
  message("One RDS file per scenario; files are updated after each replicate.")
  message("Fine grid: ", N_BLOCKS_SIDE * FINE_PER_BLOCK, " x ", N_BLOCKS_SIDE * FINE_PER_BLOCK)
  message("Sample size: ", N_SAMPLE)
}

for (ii in seq_len(nrow(scenarios))) {
  sc <- scenarios[ii, , drop = FALSE]
  if (VERBOSE >= 1) message("Scenario ", ii, "/", nrow(scenarios), ": ", sc$scenario_id)
  run_scenario_incremental(sc, N_REPS)
}

message("Simulation run complete. Next run: source('R/simulation/summarise_simulation.R')")
