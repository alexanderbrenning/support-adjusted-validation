# Plattner et al. (2006): snow water equivalent (SWE)
#
# Role in manuscript:
# - Tier 1--2 sensitivity example using the published ordinary-regression R2,
#   response standard deviation, and residual autocorrelation range.  The R2 is
#   treated in its original OLS sense, i.e. as a squared fitted-vs-observed
#   correlation proxy, not as an independent spatial-CV metric.
#
# Extracted quantities:
# - Ordinary regression R2 = 0.414.
# - SWE standard deviation = 164 mm.
# - Residual autocorrelation effective range approximately 250 m.
# - DEM/grid support = 10 m in the original analysis; an additional 100 m
#   target support is included as a sensitivity setting.
# - A nugget could not be safely estimated, but unresolved residual variability
#   below the 10 m DEM scale is described as considerable and nugget-like.
#
# Assumptions:
# - Exponential structured covariance with effective range 250 m, converted to a
#   practical 95%-range convention.
# - The visible/reported nugget-like short-scale component is represented by a
#   total observed nugget-equivalent SD of 40 mm.  This is a sensitivity value:
#   it is decomposed into measurement error plus true unresolved microscale
#   variability.  Thus measurement error is not added on top of the nugget; it
#   is treated as part of the observed nugget-like discontinuity.
# - Two measurement-error decompositions are used: 0 and 20 mm.  The 0 mm row
#   is a no-measurement-error baseline; the 20 mm row assigns one quarter of
#   the nugget-equivalent variance (20^2 / 40^2) to measurement error and the
#   remaining three quarters to true unresolved microscale variation.

source("R/examples/examples_functions.R")

OVERWRITE <- TRUE
VERBOSE <- 2

r2_obs <- 0.414
r_obs <- sqrt(r2_obs)
sigma_obs <- 164
sigma2_obs <- sigma_obs^2
mse_obs <- (1 - r2_obs) * sigma2_obs
rmse_obs <- sqrt(mse_obs)

effective_range <- 250
phi <- effective_range / 3  # exponential practical range convention

# Approximate total observed nugget-equivalent component.  This includes both
# possible measurement error and true unresolved microscale variability.
total_nugget_sd <- 40
total_nugget_observed <- total_nugget_sd^2

rows <- list()
for (cell_side in c(10, 100)) {
  # Compact two-row sensitivity set.  sigma_eps = 0 mm is the baseline in
  # which the full 40 mm nugget-equivalent component is treated as true
  # unresolved microscale variability.  sigma_eps = 20 mm assigns one quarter
  # of the nugget-equivalent variance to measurement error:
  #
  #   20^2 / 40^2 = 0.25.
  #
  # Both settings decompose the same total nugget-equivalent component;
  # measurement error is not added on top of the nugget-like term.
  for (sigma_eps in c(0, 20)) {
    sigma2_eps <- sigma_eps^2
    if (sigma2_eps > total_nugget_observed) {
      stop("Measurement-error variance exceeds total nugget-equivalent variance.")
    }

    true_nugget <- total_nugget_observed - sigma2_eps
    partial_sill <- sigma2_obs - total_nugget_observed
    sigma2_point_true <- true_nugget + partial_sill

    R <- equivalent_circle_radius(cell_side)
    gamma_true_R <- true_nugget + partial_sill * (1 - exp(-R / phi))

    terms <- equal_area_disk_linear_terms(
      sigma2_point_true = sigma2_point_true,
      sigma2_eps = sigma2_eps,
      gamma_true_R = gamma_true_R,
      true_nugget = true_nugget
    )

    res <- make_support_result(
      example_id = paste0("example_plattner_2006_swe_", cell_side, "m_eps", sigma_eps),
      source = "Plattner et al. (2006), Vernagtferner SWE regression",
      variable = "SWE",
      units = "mm",
      r_obs = r_obs,
      r2_obs = r2_obs,
      rmse_obs = rmse_obs,
      sigma2_obs = sigma2_obs,
      sigma2_eps = sigma2_eps,
      lambda = unname(terms["lambda"]),
      sigma2_delta = unname(terms["sigma2_delta"]),
      sigma2_block = unname(terms["sigma2_block"]),
      r2_definition = "squared_correlation",
      correction_type = "tier2_sensitivity",
      notes = paste(
        "Uses OLS R2=0.414 as squared-correlation proxy, SWE SD=164 mm,",
        "residual effective range=250 m, and a nugget-equivalent SD=40 mm.",
        "Measurement-error SD is either 0 or 20 mm and is treated as part of",
        "the nugget-equivalent component, not as an added variance term."
      )
    )

    res <- add_example_metadata(
      res,
      study = paste0("SWE, ", cell_side, " m"),
      prediction_support = paste0(cell_side, " m grid cell"),
      reference_support = "local snow survey site",
      information_tier = "Tier 2 sens.",
      available_information = "OLS R2, response SD, range, nugget-like term, two measurement-error decompositions",
      main_limitation = "OLS R2 used as squared-correlation proxy; nugget split into measurement error and true microscale variability is heuristic."
    )

    res$cell_side <- cell_side
    res$equivalent_radius <- R
    res$measurement_error_sd <- sigma_eps
    res$total_nugget_sd <- total_nugget_sd
    res$total_nugget_observed <- total_nugget_observed
    res$true_nugget <- true_nugget
    res$partial_sill <- partial_sill
    res$structural_covariance_model <- "Exponential"
    res$effective_range <- effective_range
    res$gamma_true_R <- gamma_true_R
    res$rho_R <- unname(terms["rho_R"])
    res$rho_struct_R <- unname(terms["rho_struct_R"])
    rows[[length(rows) + 1]] <- res
  }
}

out <- do.call(rbind, rows)
# Save one CSV for this multi-row sensitivity example.
dir.create("results/examples", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(out, "results/examples/example_plattner_2006_swe.csv", row.names = FALSE)
if (VERBOSE >= 1) message("Wrote: results/examples/example_plattner_2006_swe.csv")
print(out)
