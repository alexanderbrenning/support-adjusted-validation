# Plot simulation results
#
# Produces paper-oriented figures and additional diagnostics.
#
# Plotting conventions:
# - Paper trend figures use only the core setup:
#     scenario_type == "core", prediction_relation == "block_targeted".
# - The main paper trend figures focus on target R_B^2 = 0.60 and
#   microscale = moderate, avoiding sparse diagnostic scenarios and zigzagging
#   lines across different target-skill levels.
# - The naive point-support metric is plotted as the dominant black baseline.
# - Sparse microscale-sensitivity and failure-mode scenarios are shown in diagnostics
#   without trend lines or with facetting that drops empty combinations.

VERBOSE <- 1

IN_DIR <- file.path("results", "simulation")
SUMMARY_FILE <- file.path(IN_DIR, "simulation_summary.csv")
REPLICATE_FILE <- file.path(IN_DIR, "simulation_replicates.csv")
PAPER_DIR <- file.path("results", "paper")
DIAG_DIR <- file.path("results", "supplementary")

library(ggplot2)

source("R/core/support_correction_core.R")

dir.create(PAPER_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)

summary <- utils::read.csv(SUMMARY_FILE, stringsAsFactors = FALSE)
replicates <- if (file.exists(REPLICATE_FILE)) {
  utils::read.csv(REPLICATE_FILE, stringsAsFactors = FALSE)
} else {
  if (VERBOSE >= 1) message("Replicate file not found; continuing with summary-only plots: ", REPLICATE_FILE)
  NULL
}

#' Save ggplot as PDF and PNG.
save_plot_pair <- function(p, file_stem, dir, width = 7, height = 4.5) {
  ggplot2::ggsave(file.path(dir, paste0(file_stem, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(dir, paste0(file_stem, ".png")), p, width = width, height = height, dpi = 300)
  if (VERBOSE >= 1) message("Wrote figure pair: ", file.path(dir, file_stem))
}

#' Convert selected summary columns into a long table.
make_long_metrics <- function(x, metric_prefixes) {
  out <- list()
  k <- 1
  for (mp in metric_prefixes) {
    col <- paste0(mp, "_mean")
    if (!col %in% names(x)) next
    tmp <- x
    tmp$metric_source_raw <- mp
    tmp$value <- as.numeric(x[[col]])
    out[[k]] <- tmp
    k <- k + 1
  }
  do.call(rbind, out)
}

metric_labels_r2 <- c(
  R2_obs = "Naive / current practice",
  R2_ME = "Tier 1: measurement error",
  R2_B_hat = "Tier 2: support adjustment",
  R2_block_true = "True block target"
)

metric_labels_rmse <- c(
  RMSE_obs = "Naive / current practice",
  RMSE_ME = "Tier 1: measurement error",
  RMSE_B_hat = "Tier 2: support adjustment",
  RMSE_block_true = "True block target"
)

# The baseline is deliberately dominant. Tier 2 and the true block target often
# nearly coincide; use strongly different line types, shapes, and colours.
metric_colours <- c(
  "Naive / current practice" = "black",
  "Tier 1: measurement error" = "grey55",
  "Tier 2: support adjustment" = "#0072B2",
  "True block target" = "#D55E00"
)
metric_linetypes <- c(
  "Naive / current practice" = "solid",
  "Tier 1: measurement error" = "dashed",
  "Tier 2: support adjustment" = "dotdash",
  "True block target" = "twodash"
)
metric_shapes <- c(
  "Naive / current practice" = 16,
  "Tier 1: measurement error" = 17,
  "Tier 2: support adjustment" = 15,
  "True block target" = 1
)
metric_linewidths <- c(
  "Naive / current practice" = 1.75,
  "Tier 1: measurement error" = 0.80,
  "Tier 2: support adjustment" = 0.95,
  "True block target" = 1.10
)

paper_core <- subset(summary, scenario_type == "core" & prediction_relation == "block_targeted")
if (nrow(paper_core) == 0) stop("No core settings available for paper plots.")

target_r2_focus <- 0.60
paper_focus <- subset(paper_core, abs(target_r2_block - target_r2_focus) < 1e-8)
if (nrow(paper_focus) == 0) paper_focus <- paper_core

# -----------------------------
# Paper figure 1: R2 adjustments
# -----------------------------

r2_long <- make_long_metrics(paper_focus, c("R2_obs", "R2_ME", "R2_B_hat", "R2_block_true"))
r2_long$metric_source <- factor(metric_labels_r2[r2_long$metric_source_raw], levels = unname(metric_labels_r2))
r2_long$range_block_ratio <- factor(r2_long$range_block_ratio, levels = sort(unique(r2_long$range_block_ratio)))
r2_other <- subset(r2_long, metric_source != "Naive / current practice")
r2_base <- subset(r2_long, metric_source == "Naive / current practice")

p_r2 <- ggplot2::ggplot(
  r2_other,
  ggplot2::aes(x = q_epsilon_target, y = value, colour = metric_source,
               linetype = metric_source, shape = metric_source,
               linewidth = metric_source, group = metric_source)
) +
  ggplot2::geom_line(na.rm = TRUE) +
  ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
  ggplot2::geom_line(data = r2_base, ggplot2::aes(group = metric_source),
                     colour = "black", linewidth = 1.9, linetype = "solid", na.rm = TRUE) +
  ggplot2::geom_point(data = r2_base, colour = "black", size = 3.1, shape = 16, na.rm = TRUE) +
  ggplot2::facet_wrap(~ range_block_ratio, nrow = 1, labeller = ggplot2::label_both) +
  ggplot2::scale_colour_manual(values = metric_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = metric_linetypes, drop = FALSE) +
  ggplot2::scale_shape_manual(values = metric_shapes, drop = FALSE) +
  ggplot2::scale_linewidth_manual(values = metric_linewidths, guide = "none", drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = expression("Measurement-error fraction " * q[epsilon]),
    y = expression(R^2), colour = "Metric", linetype = "Metric", shape = "Metric",
    title = expression("Simulation adjustments for " * R[B]^2 == 0.60 * " and moderate true microscale variation"),
    caption = "The heavy black line is the naive point-support validation metric and represents current practice."
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

save_plot_pair(p_r2, "fig_simulation_r2_corrections", PAPER_DIR, width = 8.5, height = 4.5)

# -----------------------------
# Paper figure 2: RMSE adjustments
# -----------------------------

rmse_long <- make_long_metrics(paper_focus, c("RMSE_obs", "RMSE_ME", "RMSE_B_hat", "RMSE_block_true"))
rmse_long$metric_source <- factor(metric_labels_rmse[rmse_long$metric_source_raw], levels = unname(metric_labels_rmse))
rmse_long$range_block_ratio <- factor(rmse_long$range_block_ratio, levels = sort(unique(rmse_long$range_block_ratio)))
rmse_other <- subset(rmse_long, metric_source != "Naive / current practice")
rmse_base <- subset(rmse_long, metric_source == "Naive / current practice")

p_rmse <- ggplot2::ggplot(
  rmse_other,
  ggplot2::aes(x = q_epsilon_target, y = value, colour = metric_source,
               linetype = metric_source, shape = metric_source,
               linewidth = metric_source, group = metric_source)
) +
  ggplot2::geom_line(na.rm = TRUE) +
  ggplot2::geom_point(size = 2.2, na.rm = TRUE) +
  ggplot2::geom_line(data = rmse_base, ggplot2::aes(group = metric_source),
                     colour = "black", linewidth = 1.9, linetype = "solid", na.rm = TRUE) +
  ggplot2::geom_point(data = rmse_base, colour = "black", size = 3.1, shape = 16, na.rm = TRUE) +
  ggplot2::facet_wrap(~ range_block_ratio, nrow = 1, labeller = ggplot2::label_both) +
  ggplot2::scale_colour_manual(values = metric_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = metric_linetypes, drop = FALSE) +
  ggplot2::scale_shape_manual(values = metric_shapes, drop = FALSE) +
  ggplot2::scale_linewidth_manual(values = metric_linewidths, guide = "none", drop = FALSE) +
  ggplot2::labs(
    x = expression("Measurement-error fraction " * q[epsilon]),
    y = "RMSE", colour = "Metric", linetype = "Metric", shape = "Metric",
    title = expression("Simulation adjustments for " * R[B]^2 == 0.60 * " and moderate true microscale variation"),
    caption = "The heavy black line is the naive point-support validation metric and represents current practice."
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

save_plot_pair(p_rmse, "fig_simulation_rmse_corrections", PAPER_DIR, width = 8.5, height = 4.5)

# ------------------------------------------------------
# Paper figure 3: spherical semivariograms for the core setup
# ------------------------------------------------------

make_svgm_curve <- function(range_block_ratio, true_nugget, partial_sill) {
  h <- seq(0, 55, length.out = 400)
  data.frame(
    h_block = h,
    gamma_true = spherical_semivariogram(h = h, true_nugget = true_nugget,
                                         partial_sill = partial_sill,
                                         range = range_block_ratio,
                                         sigma2_eps = 0, observed = FALSE)
  )
}

# For the paper, show the core-setup semivariogram settings only.
# The none/large microscale-sensitivity settings are sparse by design and are
# shown in the diagnostic figure below instead.
svgm_grid <- unique(paper_core[, c("range_block_ratio", "microscale", "true_nugget", "partial_sill")])
svgm_grid <- svgm_grid[order(svgm_grid$microscale, svgm_grid$range_block_ratio), ]

svgm_long <- do.call(rbind, lapply(seq_len(nrow(svgm_grid)), function(i) {
  s <- svgm_grid[i, , drop = FALSE]
  cc <- make_svgm_curve(s$range_block_ratio, s$true_nugget, s$partial_sill)
  cbind(s[rep(1, nrow(cc)), , drop = FALSE], cc, row.names = NULL)
}))
svgm_long$range_block_ratio <- factor(svgm_long$range_block_ratio, levels = sort(unique(svgm_long$range_block_ratio)))

p_svgm <- ggplot2::ggplot(
  svgm_long,
  ggplot2::aes(x = h_block, y = gamma_true, linetype = range_block_ratio,
               linewidth = range_block_ratio, group = range_block_ratio)
) +
  ggplot2::geom_line(na.rm = TRUE) +
  ggplot2::scale_linewidth_manual(values = c("2" = 1.2, "10" = 0.95, "50" = 0.75), guide = "none") +
  ggplot2::labs(
    x = "Lag distance / block side length",
    y = expression(gamma[true](h)),
    linetype = "Spherical range / block side",
    title = "Latent true spherical semivariograms used in the core simulation setup"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

save_plot_pair(p_svgm, "fig_simulation_semivariogram_parameters", PAPER_DIR, width = 6.8, height = 4.2)

# The previous diagnostic-only microscale-sensitivity semivariogram plot was removed
# from the publication bundle. It duplicated scenario metadata and could be
# misleading because the microscale-sensitivity scenarios are sparse and use only
# the short range/block-side setting.

# ------------------------------------------------------
# Diagnostic: all target block R2 levels, core setup only
# ------------------------------------------------------

r2_core_all <- make_long_metrics(paper_core, c("R2_obs", "R2_ME", "R2_B_hat", "R2_block_true"))
r2_core_all$metric_source <- factor(metric_labels_r2[r2_core_all$metric_source_raw], levels = unname(metric_labels_r2))
r2_core_all$target_r2_block <- factor(r2_core_all$target_r2_block)
r2_core_all$range_block_ratio <- factor(r2_core_all$range_block_ratio)

p_r2_all <- ggplot2::ggplot(
  r2_core_all,
  ggplot2::aes(x = q_epsilon_target, y = value, colour = metric_source,
               linetype = metric_source, shape = metric_source, group = metric_source)
) +
  ggplot2::geom_line(na.rm = TRUE) +
  ggplot2::geom_point(size = 1.8, na.rm = TRUE) +
  ggplot2::facet_grid(target_r2_block ~ range_block_ratio, labeller = ggplot2::label_both) +
  ggplot2::scale_colour_manual(values = metric_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = metric_linetypes, drop = FALSE) +
  ggplot2::scale_shape_manual(values = metric_shapes, drop = FALSE) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = expression("Measurement-error fraction " * q[epsilon]), y = expression(R^2),
    colour = "Metric", linetype = "Metric", shape = "Metric",
    title = "Diagnostic: all target block-skill levels in the core setup"
  ) +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(legend.position = "bottom")

save_plot_pair(p_r2_all, "diag_r2_all_core_target_r2_levels", DIAG_DIR, width = 8.5, height = 6.5)

# ------------------------------------------------------
# Diagnostic: microscale-sensitivity contrasts
# ------------------------------------------------------

microscale_sensitivity <- subset(summary, scenario_type == "targeted_microscale" & prediction_relation == "block_targeted")
if (nrow(microscale_sensitivity) > 0) {
  sensitivity_long <- make_long_metrics(microscale_sensitivity, c("R2_obs", "R2_ME", "R2_B_hat", "R2_block_true"))
  sensitivity_long$metric_source <- factor(metric_labels_r2[sensitivity_long$metric_source_raw], levels = unname(metric_labels_r2))

  p_sensitivity <- ggplot2::ggplot(
    sensitivity_long,
    ggplot2::aes(x = microscale, y = value, colour = metric_source, shape = metric_source)
  ) +
    ggplot2::geom_point(size = 2.4, position = ggplot2::position_dodge(width = 0.5), na.rm = TRUE) +
    ggplot2::facet_wrap(~ q_epsilon_target, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = metric_colours, drop = FALSE) +
    ggplot2::scale_shape_manual(values = metric_shapes, drop = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "True microscale variation", y = expression(R^2),
      colour = "Metric", shape = "Metric",
      title = "Microscale-sensitivity scenarios; sparse design shown without trend lines"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")

  # Retain the established output filename used by the Supplement TeX.
  save_plot_pair(p_sensitivity, "diag_targeted_microscale_r2", DIAG_DIR, width = 7.5, height = 4.5)
}

# ------------------------------------------------------
# Diagnostic: lower-value interpretation
# ------------------------------------------------------

# Core diagnostic: no empty panels because the core setup is factorially complete.
p_lower_core <- ggplot2::ggplot(
  paper_core,
  ggplot2::aes(x = factor(range_block_ratio), y = lower_value_hold_rate, shape = prediction_relation)
) +
  ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
  ggplot2::facet_grid(target_r2_block ~ q_epsilon_target, labeller = ggplot2::label_both) +
  ggplot2::labs(
    x = "Range/block-size ratio",
    y = "Proportion with Tier 1 r <= true block r",
    shape = "Prediction relation",
    title = "Diagnostic: lower-value condition in the core setup"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot_pair(p_lower_core, "diag_lower_value_hold_rate_core", DIAG_DIR, width = 8.5, height = 6.0)

# Additional diagnostics: show only existing microscale-sensitivity and
# failure-mode combinations with facet_wrap.
additional_scenarios <- subset(summary, scenario_type != "core")
if (nrow(additional_scenarios) > 0) {
  additional_scenarios$facet_label <- paste0("q=", additional_scenarios$q_epsilon_target,
                                              ", micro=", additional_scenarios$microscale)
  p_lower_sparse <- ggplot2::ggplot(
    additional_scenarios,
    ggplot2::aes(x = factor(range_block_ratio), y = lower_value_hold_rate,
                 shape = prediction_relation)
  ) +
    ggplot2::geom_point(size = 2.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~ facet_label, labeller = ggplot2::label_value) +
    ggplot2::labs(
      x = "Range/block-size ratio",
      y = "Proportion with Tier 1 r <= true block r",
      shape = "Prediction relation",
      title = "Diagnostic: lower-value condition in additional scenarios"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  save_plot_pair(p_lower_sparse, "diag_lower_value_hold_rate_targeted", DIAG_DIR, width = 8.5, height = 5.0)
}

# ------------------------------------------------------
# Diagnostic: theoretical and empirical lambda
# ------------------------------------------------------

lambda_long_core <- rbind(
  transform(paper_core, lambda_type = "theoretical", lambda_value = lambda_theoretical_mean),
  transform(paper_core, lambda_type = "empirical", lambda_value = lambda_empirical_mean)
)

p_lambda_core <- ggplot2::ggplot(
  lambda_long_core,
  ggplot2::aes(x = factor(range_block_ratio), y = lambda_value, shape = lambda_type)
) +
  ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
  ggplot2::facet_grid(target_r2_block ~ q_epsilon_target, labeller = ggplot2::label_both) +
  ggplot2::labs(
    x = "Range/block-size ratio",
    y = expression(lambda),
    shape = expression(lambda~"source"),
    title = "Diagnostic: theoretical versus empirical representativeness in the core setup"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")
save_plot_pair(p_lambda_core, "diag_lambda_by_scenario_core", DIAG_DIR, width = 8.5, height = 6.0)

lambda_sparse <- subset(summary, scenario_type != "core")
if (nrow(lambda_sparse) > 0) {
  lambda_sparse_long <- rbind(
    transform(lambda_sparse, lambda_type = "theoretical", lambda_value = lambda_theoretical_mean),
    transform(lambda_sparse, lambda_type = "empirical", lambda_value = lambda_empirical_mean)
  )
  lambda_sparse_long$facet_label <- paste0("q=", lambda_sparse_long$q_epsilon_target,
                                           ", micro=", lambda_sparse_long$microscale)
  p_lambda_sparse <- ggplot2::ggplot(
    lambda_sparse_long,
    ggplot2::aes(x = factor(range_block_ratio), y = lambda_value, shape = lambda_type)
  ) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~ facet_label) +
    ggplot2::labs(
      x = "Range/block-size ratio",
      y = expression(lambda),
      shape = expression(lambda~"source"),
      title = "Diagnostic: theoretical versus empirical representativeness in additional scenarios"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  save_plot_pair(p_lambda_sparse, "diag_lambda_by_scenario_targeted", DIAG_DIR, width = 8.5, height = 5.0)
}
