# Combine and summarize simulation scenario files
#
# Reads one RDS file per scenario, combines all repetitions, and writes compact
# you want to summarize.

OVERWRITE <- TRUE
VERBOSE <- 1

IN_DIR <- file.path("results", "simulation")
OUT_PREFIX <- file.path(IN_DIR, "simulation")

if (!dir.exists(IN_DIR)) stop("Input directory not found: ", IN_DIR)

files <- list.files(IN_DIR, pattern = "_reps\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No scenario replicate files found in ", IN_DIR)

if (VERBOSE >= 1) message("Reading ", length(files), " scenario files.")
res <- do.call(rbind, lapply(files, readRDS))

#' Mean with NA handling
mean_na <- function(x) mean(x, na.rm = TRUE)
#' Standard deviation with NA handling
sd_na <- function(x) stats::sd(x, na.rm = TRUE)
#' Quantile with NA handling
q_na <- function(x, p) as.numeric(stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE))

# Summarise each numeric output by scenario.
num_cols <- names(res)[vapply(res, is.numeric, logical(1))]
id_cols <- c(
  "scenario_index", "scenario_id", "scenario_type", "q_epsilon_target",
  "range_block_ratio", "microscale", "true_nugget", "partial_sill",
  "true_sill", "sigma2_eps", "observed_nugget", "observed_sill",
  "spherical_range", "spherical_range_block_units", "target_r2_block",
  "prediction_relation"
)
id_cols <- intersect(id_cols, names(res))

metric_cols <- setdiff(
  num_cols,
  c("scenario_index", "rep", "n_blocks_side", "fine_per_block", "n_sample")
)

split_res <- split(res, res$scenario_id)
summary <- do.call(rbind, lapply(split_res, function(d) {
  first <- d[1, id_cols, drop = FALSE]
  out <- first
  out$n_reps <- nrow(d)
  out$n_sample <- d$n_sample[1]
  out$n_blocks_side <- d$n_blocks_side[1]
  out$fine_per_block <- d$fine_per_block[1]
  for (nm in metric_cols) {
    out[[paste0(nm, "_mean")]] <- mean_na(d[[nm]])
    out[[paste0(nm, "_sd")]] <- sd_na(d[[nm]])
    out[[paste0(nm, "_q025")]] <- q_na(d[[nm]], 0.025)
    out[[paste0(nm, "_q975")]] <- q_na(d[[nm]], 0.975)
  }
  out
}))
rownames(summary) <- NULL

# Bias diagnostics against the realized block-support truth.
summary$R2_obs_bias_mean <- summary$R2_obs_mean - summary$R2_block_true_mean
summary$R2_ME_bias_mean <- summary$R2_ME_mean - summary$R2_block_true_mean
summary$R2_B_hat_bias_mean <- summary$R2_B_hat_mean - summary$R2_block_true_mean
summary$RMSE_obs_bias_mean <- summary$RMSE_obs_mean - summary$RMSE_block_true_mean
summary$RMSE_ME_bias_mean <- summary$RMSE_ME_mean - summary$RMSE_block_true_mean
summary$RMSE_B_hat_bias_mean <- summary$RMSE_B_hat_mean - summary$RMSE_block_true_mean
summary$lower_value_hold_rate <- vapply(split_res, function(d) mean(d$lower_value_holds, na.rm = TRUE), numeric(1))

# Compact table of semivariogram parameters for reporting.
svgm_cols <- intersect(
  c("scenario_id", "scenario_type", "q_epsilon_target", "microscale",
    "true_nugget", "partial_sill", "true_sill", "sigma2_eps",
    "observed_nugget", "observed_sill", "spherical_range_block_units",
    "range_block_ratio", "target_r2_block", "prediction_relation",
    "lambda_theoretical_mean", "sigma2_delta_theoretical_mean"),
  names(summary)
)
svgm_table <- summary[, svgm_cols, drop = FALSE]

saveRDS(res, paste0(OUT_PREFIX, "_replicates.rds"))
utils::write.csv(res, paste0(OUT_PREFIX, "_replicates.csv"), row.names = FALSE)
saveRDS(summary, paste0(OUT_PREFIX, "_summary.rds"))
utils::write.csv(summary, paste0(OUT_PREFIX, "_summary.csv"), row.names = FALSE)
utils::write.csv(svgm_table, paste0(OUT_PREFIX, "_semivariogram_parameters.csv"), row.names = FALSE)

if (VERBOSE >= 1) {
  message("Wrote: ", paste0(OUT_PREFIX, "_replicates.csv"))
  message("Wrote: ", paste0(OUT_PREFIX, "_summary.csv"))
  message("Wrote: ", paste0(OUT_PREFIX, "_semivariogram_parameters.csv"))
}

print(summary[, intersect(
  c("scenario_id", "n_reps", "true_nugget", "partial_sill", "sigma2_eps",
    "spherical_range_block_units", "R2_obs_mean", "R2_ME_mean",
    "R2_B_hat_mean", "R2_block_true_mean", "lower_value_hold_rate"),
  names(summary)
)])
