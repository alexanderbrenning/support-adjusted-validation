# Supplementary support-adjustment function plots
#
# Supplementary adjustment-function figures for the support-adjustment framework.
#
# The figures show relative changes implied by Tier 1 and Tier 2 adjustments as
# functions of q_epsilon (Tier 1) and lambda (Tier 2). Example results from
# results/examples/posthoc_examples_combined.csv are overlaid where the required
# columns can be identified.
#
# Outputs:
#   results/examples/suppl_correction_function_example_points.csv
#   results/examples/fig_correction_functions_four_panel.pdf/png
#   results/examples/fig_correction_functions_four_panel.tex

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

library(ggrepel)
has_ggrepel <- TRUE

`%||%` <- function(a, b) if (!is.null(a)) a else b

out_dir <- file.path("results", "examples")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

combined_file <- file.path(out_dir, "posthoc_examples_combined.csv")

# -----------------------------------------------------------------------------
# User-adjustable settings
# -----------------------------------------------------------------------------
# The adjustment-function plots are intended as rule-of-thumb figures. Restricting q to
# 0--0.4 keeps the displays in the practically relevant range used in the paper.
q_grid <- seq(0, 0.40, by = 0.0025)
lambda_grid <- seq(0.001, 1.00, by = 0.0025)

# RMSE adjustment curves depend on MSE_obs / sigma2_obs, equivalently on raw
# predictive R2 when that metric is available, not only on q_epsilon or lambda.
# These values index the background curves. Example points are assigned to the
# closest panel by their reported or reconstructed raw R2; the neutral facet
# heading is intentional because the overlay contains both R2 definitions.
rmse_reference_r2 <- c(0.25, 0.50, 0.75)

# Substantial-change thresholds: 5% relative R2 increase; 10% relative RMSE decrease.
r2_substantial_ratio <- 1.05
nrmse_substantial_ratio <- 0.90

q_xlim <- c(0, max(q_grid))
lambda_xlim <- c(0, 1)
lambda_ylim <- c(0, 1)

shade_col <- "grey88"
main_curve_col <- "#2166AC"
point_col <- "black"
leader_col <- "grey45"
grid_major <- "grey68"
grid_minor <- "grey82"

# Fixed plotting ranges for the simple Tier 1 panels. Use finite rectangles,
# not +/-Inf rectangles, to avoid device/version-specific clipping problems.
p1a_ylim <- c(1.00, 1.50)
p1b_ylim <- c(0.30, 1.00)


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------
first_existing <- function(nm, candidates) {
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1]] else NA_character_
}

as_num <- function(x) suppressWarnings(as.numeric(x))

safe_ratio <- function(num, den) {
  out <- as_num(num) / as_num(den)
  out[!is.finite(out)] <- NA_real_
  out
}

r2_ratio_tier1 <- function(q) 1 / (1 - q)
r2_ratio_tier2 <- function(lambda) 1 / lambda

# Tier 1 RMSE ratio after reference-error adjustment.
# raw_r2_pred is R2_pred,obs = 1 - MSE_obs / S_obs^2.
rmse_ratio_tier1 <- function(q, raw_r2_pred) {
  m <- 1 - raw_r2_pred
  out <- sqrt((m - q) / m)
  out[!is.finite(out) | (m - q) < 0] <- NA_real_
  out
}

# Tier 2 RMSE ratio using lambda = Var(block) / Var(observed reference).
# Since Delta + sigma2_eps = Sobs^2 - S_B^2 = Sobs^2(1-lambda), the relative
# RMSE adjustment depends on lambda and raw R2_pred, but not separately on q_eps.
rmse_ratio_tier2 <- function(lambda, raw_r2_pred) {
  m <- 1 - raw_r2_pred
  out <- sqrt((m - (1 - lambda)) / m)
  out[!is.finite(out) | (m - (1 - lambda)) < 0] <- NA_real_
  out
}

assign_r2_panel <- function(r2) {
  r2 <- as_num(r2)
  out <- ifelse(is.na(r2), NA_real_,
                ifelse(r2 < 0.375, 0.25,
                       ifelse(r2 <= 0.625, 0.50, 0.75)))
  out
}

save_plot <- function(p, basename, width, height) {
  ggsave(file.path(out_dir, paste0(basename, ".pdf")), p, width = width, height = height, units = "in")
  ggsave(file.path(out_dir, paste0(basename, ".png")), p, width = width, height = height, units = "in", dpi = 300)
}

label_layer <- function(data, mapping, size = 2.8) {
  if (!nrow(data)) return(NULL)
  if (has_ggrepel) {
    ggrepel::geom_label_repel(
      mapping = mapping,
      data = data,
      size = size,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.30,
      segment.colour = leader_col,
      segment.size = 0.22,
      label.size = 0.12,
      label.r = grid::unit(0.05, "lines"),
      label.padding = grid::unit(0.10, "lines"),
      fill = "white",
      alpha = 0.92,
      seed = 17,
      force = 1.7,
      force_pull = 0.35,
      show.legend = FALSE
    )
  } else {
    geom_label(mapping = mapping, data = data, size = size,
               nudge_y = 0.02, label.size = 0.12,
               check_overlap = TRUE, show.legend = FALSE)
  }
}

shorten_example_labels <- function(x) {
  x <- as.character(x)

  # Use typographic-neutral hyphens in this figure only.
  # Literal replacements use fixed = TRUE deliberately; only the compact
  # support/q label patterns below are regular expressions.
  x <- gsub("--", "-", x, fixed = TRUE)
  x <- gsub("–", "-", x, fixed = TRUE)
  x <- gsub("—", "-", x, fixed = TRUE)
  x <- gsub("\\textendash{}", "-", x, fixed = TRUE)
  x <- gsub("\\textendash", "-", x, fixed = TRUE)

  # Dataset-name replacements for this figure only.
  x <- gsub("Zandler et al. (2015)", "Biomass", x, fixed = TRUE)
  x <- gsub("Plattner et al. (2006)", "SWE", x, fixed = TRUE)
  x <- gsub("Brenning et al. (2005)", "BTS", x, fixed = TRUE)
  x <- gsub("Balenzano et al. (2021)", "S1-SM", x, fixed = TRUE)
  x <- gsub("Zhao et al. (2021)", "SMAP-CRNS", x, fixed = TRUE)
  x <- gsub("Vizcaino & Lavalle (2018)", "NO2-EU", x, fixed = TRUE)
  x <- gsub("NO2 Germany 2018", "NO2-DE", x, fixed = TRUE)
  x <- gsub("Meuse log zinc", "Meuse-Zn", x, fixed = TRUE)
  x <- gsub("Meuse log cadmium", "Meuse-Cd", x, fixed = TRUE)
  x <- gsub("Meuse--Zn", "Meuse-Zn", x, fixed = TRUE)
  x <- gsub("Meuse--Cd", "Meuse-Cd", x, fixed = TRUE)

  # Sensitivity-label compaction.
  x <- gsub("50/50 nugget split", "50/50 nugget", x, fixed = TRUE)
  x <- gsub("50% microscale, 50% meas. error", "50/50 nugget", x, fixed = TRUE)
  x <- gsub("50% microscale, 50% measurement error", "50/50 nugget", x, fixed = TRUE)
  x <- gsub("full nugget microscale", "full nugget", x, fixed = TRUE)
  x <- gsub("full nugget micro", "full nugget", x, fixed = TRUE)
  x <- gsub("100% nugget microscale", "full nugget", x, fixed = TRUE)
  x <- gsub("100% nugget micro", "full nugget", x, fixed = TRUE)

  # Remove comma in compact support labels, e.g. "SWE, 100 m" -> "SWE 100 m".
  # These three calls are intentionally regex-based.
  x <- gsub("SWE, *([0-9]+) *m", "SWE \\1 m", x)
  x <- gsub("SWE ([0-9]+) m, *q *= *", "SWE \\1 m q=", x)
  x <- gsub(", *q *= *", " q=", x)

  # Remove remaining double hyphens after replacements, just in case.
  x <- gsub("--", "-", x, fixed = TRUE)

  x
}

format_compact_q <- function(q) {
  out <- ifelse(is.finite(q), sprintf("%.3f", q), NA_character_)
  out <- sub("^0", "", out)
  out <- sub("\\.?0+$", "", out)
  blank <- !is.na(out) & out == ""
  out[blank] <- "0"
  out
}


collapse_label_set <- function(labs) {
  labs <- unique(labs[!is.na(labs) & nzchar(labs)])
  if (!length(labs)) return("")

  # Work only on labels already shortened for this figure. Literal string
  # replacements use fixed = TRUE; regex is used only for the compact support
  # patterns.
  labs <- gsub("--", "-", labs, fixed = TRUE)
  labs <- gsub("–", "-", labs, fixed = TRUE)
  labs <- gsub("—", "-", labs, fixed = TRUE)

  # Merge Meuse zinc/cadmium labels when their points coincide or nearly coincide.
  has_meuse_zn <- any(grepl("^Meuse-Zn($|[[:space:]])", labs))
  has_meuse_cd <- any(grepl("^Meuse-Cd($|[[:space:]])", labs))
  if (has_meuse_zn && has_meuse_cd) {
    labs <- labs[!grepl("^Meuse-(Zn|Cd)($|[[:space:]])", labs)]
    labs <- c("Meuse-Zn/Cd", labs)
  }

  # Merge SWE 10 m and 100 m labels when their points coincide or nearly coincide.
  # This intentionally drops the support-specific distinction only in this figure.
  has_swe_10 <- any(grepl("^SWE[[:space:]]+10[[:space:]]*m($|[[:space:]])", labs))
  has_swe_100 <- any(grepl("^SWE[[:space:]]+100[[:space:]]*m($|[[:space:]])", labs))
  if (has_swe_10 && has_swe_100) {
    labs <- labs[!grepl("^SWE[[:space:]]+(10|100)[[:space:]]*m($|[[:space:]])", labs)]
    labs <- c("SWE 10/100 m", labs)
  }

  paste(labs, collapse = " / ")
}

merge_near_labels <- function(data, x_col, y_col,
                              x_tol = 0.006,
                              y_tol = 0.018,
                              disambiguate_duplicate_labels = TRUE) {
  if (!nrow(data)) return(data)

  d <- data
  d <- d[is.finite(d[[x_col]]) & is.finite(d[[y_col]]), , drop = FALSE]
  if (!nrow(d)) return(d)

  d$.x_bin <- round(d[[x_col]] / x_tol)
  d$.y_bin <- round(d[[y_col]] / y_tol)
  d$.grp <- paste(d$.x_bin, d$.y_bin, sep = "_")

  groups <- split(d, d$.grp)
  merged <- lapply(groups, function(g) {
    out <- g[1, , drop = FALSE]
    out[[x_col]] <- mean(g[[x_col]], na.rm = TRUE)
    out[[y_col]] <- mean(g[[y_col]], na.rm = TRUE)
    if ("q_epsilon_plot" %in% names(g)) {
      q_mean <- mean(g$q_epsilon_plot, na.rm = TRUE)
      out$q_epsilon_plot <- if (is.nan(q_mean)) NA_real_ else q_mean
    }
    labs <- unique(g$example_label_plot[!is.na(g$example_label_plot) & nzchar(g$example_label_plot)])
    out$example_label_plot <- collapse_label_set(labs)
    out
  })
  d2 <- do.call(rbind, merged)
  rownames(d2) <- NULL
  d2$.x_bin <- d2$.y_bin <- d2$.grp <- NULL

  if (disambiguate_duplicate_labels && "q_epsilon_plot" %in% names(d2)) {
    dup <- duplicated(d2$example_label_plot) | duplicated(d2$example_label_plot, fromLast = TRUE)
    has_q <- is.finite(d2$q_epsilon_plot)
    d2$example_label_plot[dup & has_q] <- paste0(
      d2$example_label_plot[dup & has_q],
      " q=",
      format_compact_q(d2$q_epsilon_plot[dup & has_q])
    )
  }

  d2
}

substantial_q_r2 <- 1 - 1 / r2_substantial_ratio
# For Tier 1 RMSE: sqrt((m-q)/m)=0.90 -> q = m*(1-0.90^2).
substantial_q_rmse <- data.frame(
  raw_r2_reference = rmse_reference_r2,
  q_substantial = (1 - rmse_reference_r2) * (1 - nrmse_substantial_ratio^2)
)
# For Tier 2 R2: 1/lambda=1.05 -> lambda = 1/1.05.
substantial_lambda_r2 <- 1 / r2_substantial_ratio
# For Tier 2 RMSE: sqrt((m-(1-lambda))/m)=0.90 -> lambda=1-m*(1-0.90^2).
substantial_lambda_rmse <- data.frame(
  raw_r2_reference = rmse_reference_r2,
  lambda_substantial = 1 - (1 - rmse_reference_r2) * (1 - nrmse_substantial_ratio^2)
)

base_theme <- function() {
  theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(colour = grid_major, linewidth = 0.25),
      panel.grid.minor = element_line(colour = grid_minor, linewidth = 0.15),
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "grey95", colour = "grey70")
    )
}

# -----------------------------------------------------------------------------
# Theoretical Tier 1 data
# -----------------------------------------------------------------------------
tier1_r2 <- data.frame(
  q_epsilon = q_grid,
  ratio = r2_ratio_tier1(q_grid)
)
tier1_r2$substantial <- tier1_r2$ratio >= r2_substantial_ratio

tier1_rmse <- do.call(rbind, lapply(rmse_reference_r2, function(r2ref) {
  data.frame(
    q_epsilon = q_grid,
    raw_r2_reference = r2ref,
    ratio = rmse_ratio_tier1(q_grid, r2ref)
  )
}))
tier1_rmse$substantial <- tier1_rmse$ratio <= nrmse_substantial_ratio

# -----------------------------------------------------------------------------
# Theoretical Tier 2 data
# -----------------------------------------------------------------------------
grid2 <- expand.grid(q_epsilon = q_grid, lambda = lambda_grid)
# Feasible upper envelope when lambda is decomposed into measurement reliability
# and structural/support representativeness: lambda <= 1 - q_epsilon.
grid2$feasible_support_variance <- grid2$lambda <= (1 - grid2$q_epsilon)

tier2_r2 <- transform(grid2, ratio = r2_ratio_tier2(lambda))
tier2_r2$substantial <- tier2_r2$ratio >= r2_substantial_ratio

tier2_r2_panel <- do.call(rbind, lapply(rmse_reference_r2, function(r2ref) {
  out <- tier2_r2
  out$raw_r2_reference <- r2ref
  out
}))

tier2_rmse <- do.call(rbind, lapply(rmse_reference_r2, function(r2ref) {
  out <- grid2
  out$raw_r2_reference <- r2ref
  out$ratio <- rmse_ratio_tier2(out$lambda, r2ref)
  out$substantial <- out$ratio <= nrmse_substantial_ratio
  out
}))

# Filter finite contour data to avoid irrelevant stat_contour warnings.
tier2_r2_contour <- subset(tier2_r2_panel, is.finite(ratio))
tier2_rmse_contour <- subset(tier2_rmse, is.finite(ratio))

# -----------------------------------------------------------------------------
# Example points from the post-hoc results table
# -----------------------------------------------------------------------------
examples <- data.frame()
if (file.exists(combined_file)) {
  examples <- read.csv(combined_file, stringsAsFactors = FALSE, check.names = FALSE)
  nm <- names(examples)

  lab_col <- first_existing(nm, c("study", "source", "example_label", "example_id"))
  q_col <- first_existing(nm, c("q_epsilon", "q_eps", "q", "qeps"))
  lam_col <- first_existing(nm, c(
    "lambda_A_B", "lambda_p_B", "lambda", "lambda_center", "lambda_support",
    "lambda_tier2", "lambda_used", "representativeness_factor"
  ))
  r2_raw_col <- first_existing(nm, c("r2_raw", "r2_obs", "R2_raw", "R2_obs"))
  r2_t1_col <- first_existing(nm, c("r2_tier1", "r2_Tier1", "r2_lower", "R2_Tier1", "r2_ME", "R2_ME"))
  r2_t2_col <- first_existing(nm, c("r2_tier2", "r2_Tier2", "r2_upper", "R2_Tier2", "r2_B", "R2_B"))
  rmse_raw_col <- first_existing(nm, c("rmse_raw", "rmse_obs", "RMSE_raw", "RMSE_obs"))
  rmse_t1_col <- first_existing(nm, c("rmse_tier1", "rmse_Tier1", "rmse_measurement_error_only", "RMSE_Tier1", "rmse_ME", "RMSE_ME"))
  rmse_t2_col <- first_existing(nm, c("rmse_tier2", "rmse_Tier2", "rmse_block", "RMSE_Tier2", "rmse_B", "RMSE_B"))
  r2_def_col <- first_existing(nm, c("r2_definition", "R2_definition"))

  examples$example_label_plot <- if (!is.na(lab_col)) examples[[lab_col]] else examples[[1]]
  examples$q_epsilon_plot <- if (!is.na(q_col)) as_num(examples[[q_col]]) else NA_real_
  examples$lambda_plot <- if (!is.na(lam_col)) as_num(examples[[lam_col]]) else NA_real_
  examples$r2_raw_plot <- if (!is.na(r2_raw_col)) as_num(examples[[r2_raw_col]]) else NA_real_
  examples$r2_tier1_plot <- if (!is.na(r2_t1_col)) as_num(examples[[r2_t1_col]]) else NA_real_
  examples$r2_tier2_plot <- if (!is.na(r2_t2_col)) as_num(examples[[r2_t2_col]]) else NA_real_
  examples$rmse_raw_plot <- if (!is.na(rmse_raw_col)) as_num(examples[[rmse_raw_col]]) else NA_real_
  examples$rmse_tier1_plot <- if (!is.na(rmse_t1_col)) as_num(examples[[rmse_t1_col]]) else NA_real_
  examples$rmse_tier2_plot <- if (!is.na(rmse_t2_col)) as_num(examples[[rmse_t2_col]]) else NA_real_
  examples$r2_definition_plot <- if (!is.na(r2_def_col)) examples[[r2_def_col]] else NA_character_

  examples$r2_tier1_ratio <- safe_ratio(examples$r2_tier1_plot, examples$r2_raw_plot)
  examples$r2_tier2_ratio <- safe_ratio(examples$r2_tier2_plot, examples$r2_raw_plot)
  examples$rmse_tier1_ratio <- safe_ratio(examples$rmse_tier1_plot, examples$rmse_raw_plot)
  examples$rmse_tier2_ratio <- safe_ratio(examples$rmse_tier2_plot, examples$rmse_raw_plot)
  examples$raw_r2_reference <- assign_r2_panel(examples$r2_raw_plot)

  # Keep labels short enough for the manuscript sensitivity figure only.
  # Summary tables and example-result files keep their own labels unchanged.
  examples$example_label_plot <- shorten_example_labels(examples$example_label_plot)

  write.csv(examples, file.path(out_dir, "suppl_correction_function_example_points.csv"), row.names = FALSE)
} else {
  message("No combined post-hoc example file found at ", combined_file, "; figures will contain theoretical curves only.")
}

examples_t1_r2 <- subset(examples,
                         is.finite(q_epsilon_plot) & q_epsilon_plot >= q_xlim[1] & q_epsilon_plot <= q_xlim[2] &
                           is.finite(r2_tier1_ratio) & r2_tier1_ratio >= p1a_ylim[1] & r2_tier1_ratio <= p1a_ylim[2] &
                           r2_tier1_ratio > 1.01)
examples_t1_rmse <- subset(examples,
                           is.finite(q_epsilon_plot) & q_epsilon_plot >= q_xlim[1] & q_epsilon_plot <= q_xlim[2] &
                             is.finite(rmse_tier1_ratio) & rmse_tier1_ratio >= p1b_ylim[1] & rmse_tier1_ratio <= p1b_ylim[2] &
                             rmse_tier1_ratio < 0.99)
examples_t2_r2 <- subset(examples,
                         is.finite(q_epsilon_plot) & q_epsilon_plot >= q_xlim[1] & q_epsilon_plot <= q_xlim[2] &
                           is.finite(lambda_plot) & is.finite(r2_tier2_ratio) &
                           is.finite(raw_r2_reference))
examples_t2_rmse <- subset(examples,
                           is.finite(q_epsilon_plot) & q_epsilon_plot >= q_xlim[1] & q_epsilon_plot <= q_xlim[2] &
                             is.finite(lambda_plot) & is.finite(rmse_tier2_ratio) &
                             is.finite(raw_r2_reference))

# Label data for the manuscript panels. Points with coincident or nearly
# coincident locations are represented by one label wherever possible. This keeps
# the sensitivity plot readable when several SWE or Meuse scenarios imply almost
# the same relative adjustment.
examples_t1_r2_lab <- merge_near_labels(
  examples_t1_r2, "q_epsilon_plot", "r2_tier1_ratio",
  x_tol = 0.008, y_tol = 0.020,
  disambiguate_duplicate_labels = TRUE
)
examples_t1_rmse_lab <- merge_near_labels(
  examples_t1_rmse, "q_epsilon_plot", "rmse_tier1_ratio",
  x_tol = 0.008, y_tol = 0.025,
  disambiguate_duplicate_labels = TRUE
)

# -----------------------------------------------------------------------------
# Tier 1 figures
# -----------------------------------------------------------------------------
p1a <- ggplot(tier1_r2, aes(q_epsilon, ratio)) +
  annotate("rect", xmin = q_xlim[1], xmax = q_xlim[2],
           ymin = r2_substantial_ratio, ymax = p1a_ylim[2],
           fill = shade_col, alpha = 0.28) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_hline(yintercept = r2_substantial_ratio, linewidth = 0.25, linetype = "dashed") +
  geom_line(data = subset(tier1_r2, is.finite(ratio)), colour = main_curve_col, linewidth = 0.9) +
  geom_point(data = examples_t1_r2,
             aes(x = q_epsilon_plot, y = r2_tier1_ratio),
             colour = point_col, size = 1.8) +
  label_layer(examples_t1_r2_lab,
              aes(x = q_epsilon_plot, y = r2_tier1_ratio, label = example_label_plot),
              size = 2.9) +
  coord_cartesian(xlim = q_xlim, ylim = p1a_ylim, expand = FALSE) +
  labs(x = expression(q[epsilon]), y = expression(R^2~adjusted/raw),
       title = "Tier 1 reference-error adjustment for R²",
       subtitle = "Same relative change for predictive and squared-correlation R².") +
  base_theme()

p1b <- ggplot(tier1_rmse, aes(q_epsilon, ratio)) +
  annotate("rect", xmin = q_xlim[1], xmax = q_xlim[2],
           ymin = p1b_ylim[1], ymax = nrmse_substantial_ratio,
           fill = shade_col, alpha = 0.28) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_hline(yintercept = nrmse_substantial_ratio, linewidth = 0.25, linetype = "dashed") +
  geom_line(data = subset(tier1_rmse, is.finite(ratio)), colour = main_curve_col, linewidth = 0.9) +
  geom_point(data = examples_t1_rmse,
             aes(x = q_epsilon_plot, y = rmse_tier1_ratio),
             colour = point_col, size = 1.8) +
  label_layer(examples_t1_rmse_lab,
              aes(x = q_epsilon_plot, y = rmse_tier1_ratio, label = example_label_plot),
              size = 2.9) +
  facet_wrap(~ raw_r2_reference,
             labeller = as_labeller(function(x) paste0("Raw R² = ", x)), nrow = 1) +
  coord_cartesian(xlim = q_xlim, ylim = p1b_ylim, expand = FALSE) +
  labs(x = expression(q[epsilon]), y = "RMSE adjusted/raw",
       title = "Tier 1 reference-error adjustment for RMSE",
       subtitle = "Curves use raw predictive R²; points include both R² definitions.") +
  base_theme()


# -----------------------------------------------------------------------------
# Tier 2 figures
# -----------------------------------------------------------------------------
# With lambda = Var(block) / Var(observed reference), the relative Tier 2 R2
# adjustment is 1/lambda for both squared-correlation and predictive R2.
# The relative Tier 2 RMSE adjustment additionally depends on the raw predictive
# R2 because RMSE is an MSE-scale quantity.

p2a_ylim <- c(1.00, 3.00)
p2b_ylim <- c(0.30, 1.00)

# Use lambda as the x-axis. Keep the line data finite and inside the displayed
# range to avoid avoidable geom_line warnings.
tier2_r2_line <- data.frame(
  lambda = lambda_grid,
  ratio = r2_ratio_tier2(lambda_grid)
)
tier2_r2_line_plot <- subset(tier2_r2_line, is.finite(ratio))

# For Tier 2 RMSE, one curve is needed for each reference raw predictive R2.
tier2_rmse_line <- do.call(rbind, lapply(rmse_reference_r2, function(r2ref) {
  data.frame(
    lambda = lambda_grid,
    raw_r2_reference = r2ref,
    ratio = rmse_ratio_tier2(lambda_grid, r2ref)
  )
}))
tier2_rmse_line$substantial <- tier2_rmse_line$ratio <= nrmse_substantial_ratio
tier2_rmse_line_plot <- subset(tier2_rmse_line, is.finite(ratio))

# Thresholds for the line-graph version.
# R2: 1/lambda = 1.05 -> lambda = 1/1.05.
substantial_lambda_r2 <- 1 / r2_substantial_ratio
# RMSE: sqrt((lambda - r2_raw) / (1 - r2_raw)) = 0.90.
substantial_lambda_rmse <- data.frame(
  raw_r2_reference = rmse_reference_r2,
  lambda_substantial = rmse_reference_r2 +
    (1 - rmse_reference_r2) * nrmse_substantial_ratio^2
)

examples_t2_r2_line <- subset(examples_t2_r2,
                              is.finite(lambda_plot) & lambda_plot >= lambda_xlim[1] & lambda_plot <= lambda_xlim[2] &
                                is.finite(r2_tier2_ratio) & r2_tier2_ratio >= p2a_ylim[1] & r2_tier2_ratio <= p2a_ylim[2])
examples_t2_rmse_line <- subset(examples_t2_rmse,
                                is.finite(lambda_plot) & lambda_plot >= lambda_xlim[1] & lambda_plot <= lambda_xlim[2] &
                                  is.finite(rmse_tier2_ratio) & rmse_tier2_ratio >= p2b_ylim[1] & rmse_tier2_ratio <= p2b_ylim[2] &
                                  is.finite(raw_r2_reference))

examples_t2_r2_line_lab <- merge_near_labels(
  examples_t2_r2_line, "lambda_plot", "r2_tier2_ratio",
  x_tol = 0.018, y_tol = 0.050,
  disambiguate_duplicate_labels = FALSE
)
examples_t2_rmse_line_lab <- merge_near_labels(
  examples_t2_rmse_line, "lambda_plot", "rmse_tier2_ratio",
  x_tol = 0.018, y_tol = 0.035,
  disambiguate_duplicate_labels = FALSE
)

p2a <- ggplot(tier2_r2_line_plot, aes(lambda, ratio)) +
  annotate("rect", xmin = lambda_xlim[1], xmax = lambda_xlim[2],
           ymin = r2_substantial_ratio, ymax = p2a_ylim[2],
           fill = shade_col, alpha = 0.28) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_hline(yintercept = r2_substantial_ratio, linewidth = 0.25, linetype = "dashed") +
  geom_line(colour = main_curve_col, linewidth = 0.9) +
  geom_point(data = examples_t2_r2_line,
             aes(x = lambda_plot, y = r2_tier2_ratio),
             colour = point_col, size = 1.8) +
  label_layer(examples_t2_r2_line_lab,
              aes(x = lambda_plot, y = r2_tier2_ratio, label = example_label_plot),
              size = 2.9) +
  coord_cartesian(xlim = lambda_xlim, ylim = p2a_ylim, expand = FALSE) +
  labs(x = expression(lambda[A*','*B]), y = expression(R^2~adjusted/raw),
       title = "Tier 2 support adjustment for R²",
       subtitle = "Same relative change for predictive and squared-correlation R².") +
  base_theme()

p2b <- ggplot(tier2_rmse_line_plot, aes(lambda, ratio)) +
  annotate("rect", xmin = lambda_xlim[1], xmax = lambda_xlim[2],
           ymin = p2b_ylim[1], ymax = nrmse_substantial_ratio,
           fill = shade_col, alpha = 0.28) +
  geom_hline(yintercept = 1, linewidth = 0.3) +
  geom_hline(yintercept = nrmse_substantial_ratio, linewidth = 0.25, linetype = "dashed") +
  geom_line(colour = main_curve_col, linewidth = 0.9) +
  geom_point(data = examples_t2_rmse_line,
             aes(x = lambda_plot, y = rmse_tier2_ratio),
             colour = point_col, size = 1.8) +
  label_layer(examples_t2_rmse_line_lab,
              aes(x = lambda_plot, y = rmse_tier2_ratio, label = example_label_plot),
              size = 2.8) +
  facet_wrap(~ raw_r2_reference,
             labeller = as_labeller(function(x) paste0("Raw R² = ", x)), nrow = 1) +
  coord_cartesian(xlim = lambda_xlim, ylim = p2b_ylim, expand = FALSE) +
  labs(x = expression(lambda[A*','*B]), y = "RMSE adjusted/raw",
       title = "Tier 2 support adjustment for RMSE",
       subtitle = "Curves use raw predictive R²; points include both R² definitions.") +
  base_theme()


# -----------------------------------------------------------------------------
# Four-panel manuscript figure
# -----------------------------------------------------------------------------
# Combines:
#   a) Tier 1 R2 adjustment vs q_epsilon
#   b) Tier 1 RMSE adjustment vs q_epsilon
#   c) Tier 2 R2 adjustment vs lambda
#   d) Tier 2 RMSE adjustment vs lambda
# A small spacer column increases the distance between the left and right panels.

draw_four_panel <- function(p1, p2, p3, p4,
                            labels = c("a)", "b)", "c)", "d)"),
                            left_width = 1,
                            gap_width = 0.10,
                            right_width = 1,
                            top_height = 1,
                            bottom_height = 1) {
  grid::grid.newpage()

  lay <- grid::grid.layout(
    nrow = 2,
    ncol = 3,
    widths = grid::unit(c(left_width, gap_width, right_width), "null"),
    heights = grid::unit(c(top_height, bottom_height), "null")
  )
  grid::pushViewport(grid::viewport(layout = lay))

  plot_list <- list(p1, p2, p3, p4)
  positions <- data.frame(
    row = c(1, 1, 2, 2),
    col = c(1, 3, 1, 3)
  )

  for (i in seq_along(plot_list)) {
    vp <- grid::viewport(
      layout.pos.row = positions$row[i],
      layout.pos.col = positions$col[i]
    )
    print(plot_list[[i]], vp = vp)

    grid::pushViewport(vp)
    grid::grid.text(
      labels[i],
      x = grid::unit(0.015, "npc"),
      y = grid::unit(0.985, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontface = "bold", fontsize = 11)
    )
    grid::popViewport()
  }

  grid::popViewport()
  invisible(NULL)
}

save_four_panel <- function(filename_base,
                            width = 11.2,
                            height = 7.4,
                            dpi = 300) {
  pdf_file <- file.path(out_dir, paste0(filename_base, ".pdf"))
  png_file <- file.path(out_dir, paste0(filename_base, ".png"))

  grDevices::pdf(pdf_file, width = width, height = height)
  draw_four_panel(p1a, p1b, p2a, p2b)
  grDevices::dev.off()

  grDevices::png(
    png_file,
    width = width,
    height = height,
    units = "in",
    res = dpi
  )
  draw_four_panel(p1a, p1b, p2a, p2b)
  grDevices::dev.off()

  message("Wrote four-panel adjustment-function figure: ", pdf_file)
  message("Wrote four-panel adjustment-function figure: ", png_file)
}

make_four_panel_caption_file <- function(
    out_tex = file.path(out_dir, "fig_correction_functions_four_panel.tex"),
    figure_path = "../code_data/results/examples/fig_correction_functions_four_panel.pdf") {
  lines <- c(
    "% Automatically generated by R/paper/make_suppl_correction_function_plots.R",
    "\\begin{figure}[ht]",
    "\\centering",
    paste0("\\includegraphics[width=\\textwidth]{", figure_path, "}"),
    "\\caption{Relative support-adjustment functions and empirical example locations. Panels (a) and (b) show Tier~1 reference-error adjustment as a function of measurement-error fraction $q_\\varepsilon$; panel (a) applies to both $R^2_{\\mathrm{cor}}$ and $R^2_{\\mathrm{pred}}$, whereas panel (b) shows the RMSE ratio for three representative raw predictive $R^2$ levels. Panels (c) and (d) show Tier~2 support adjustment as a function of the representativeness factor $\\lambda_{A,B}$; panel (c) again applies to both $R^2$ definitions, and panel (d) shows the RMSE ratio for the same three raw predictive $R^2$ levels. The RMSE-panel facet headings use the neutral label raw $R^2$ because overlaid examples include both $R^2$ definitions; each point is assigned to the closest panel using its reported or reconstructed raw $R^2$, while the background curve retains its predictive-$R^2$ interpretation. Grey shading marks changes that meet or exceed the relevance thresholds used in the manuscript (at least 5\\% relative increase in $R^2$ or at least 10\\% relative decrease in RMSE). For legibility, Tier~1 example points are shown only when the relative adjustment exceeds 1\\%, i.e. $R^2_\\mathrm{adj}/R^2_\\mathrm{raw} > 1.01$ in panel (a) and $\\mathrm{RMSE}_\\mathrm{adj}/\\mathrm{RMSE}_\\mathrm{raw} < 0.99$ in panel (b).}",
    "\\label{fig:suppl-correction-functions-four-panel}",
    "\\end{figure}"
  )
  writeLines(lines, out_tex)
  invisible(out_tex)
}

save_four_panel("fig_correction_functions_four_panel")
make_four_panel_caption_file()

message("Wrote support-adjustment function figures to: ", out_dir)
message("Note: RMSE adjustment curves/surfaces require a reference raw predictive R2; see rmse_reference_r2 in this script.")
message("Note: Tier 2 panels use lambda on the x-axis. Once lambda = Var(block)/Var(observed reference) is fixed, q_epsilon is not a separate argument of the relative R2/RMSE ratios; its effect is contained in lambda and in the MSE-scale decomposition used to produce each example point.")
