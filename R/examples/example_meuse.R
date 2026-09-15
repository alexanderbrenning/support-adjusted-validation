# Meuse--Zn and Meuse--Cd: log-zinc and log-cadmium examples
#
# Role in manuscript/code supplement:
# - Reproducible implementation-check examples using the built-in Meuse data.
# - The soil observations are bulk samples with finite sampling support, but are
#   treated as point-support observations here to exercise the adjustment code.
#
# Data:
# - sp::meuse      station/sample observations
# - sp::meuse.grid interpolation grid and gridded covariates
#
# Response variables:
# - log(zinc)
# - log(cadmium)
#
# Supports:
# - Prediction target: 100 m Meuse interpolation-grid cell.
# - Reference support: treated as point-support soil measurement.
#
# Model/performance:
# - Simple linear regression models are fitted separately for log-zinc and
#   log-cadmium using the same available predictors.
# - Leave-one-out cross-validation (LOOCV) is used to compute predictive
#   RMSE and predictive R2:
#       R2_pred = 1 - MSE_LOOCV / Var(response).
# - For comparison and audit, the squared correlation between LOOCV predictions
#   and observations is also retained and support-adjusted through lambda:
#       R2_cor = cor(pred_LOOCV, obs)^2,
#       R2_cor,B = R2_cor / lambda,
#   with out-of-range values flagged rather than truncated.
#
# Correction status:
# - Measurement error is unknown and set to zero.
# - A fitted response semivariogram supplies the support-variance terms for a
#   100 m grid cell.
# - Because measurement error is not identified separately from any fitted
#   nugget-like component, these are Tier 3 support-sensitivity examples, not
#   fully identified Tier 2 adjustments.

source("R/examples/examples_functions.R")
library(sf)
library(sp)
library(gstat)

OVERWRITE <- TRUE
VERBOSE <- 2

library(sf)
library(sp)
library(gstat)

data(meuse, package = "sp")
data(meuse.grid, package = "sp")

meuse_sf <- sf::st_as_sf(meuse, coords = c("x", "y"), crs = 28992, remove = FALSE)
meuse_grid_sf <- sf::st_as_sf(meuse.grid, coords = c("x", "y"), crs = 28992, remove = FALSE)

candidate_predictors <- c("x", "y", "dist", "ffreq", "soil", "part.a", "part.b", "elev")
available_predictors <- intersect(candidate_predictors, intersect(names(meuse_sf), names(meuse_grid_sf)))
if (length(available_predictors) < 2) {
  available_predictors <- intersect(c("x", "y", "dist", "ffreq", "soil"), names(meuse_sf))
}

cell_side <- 100
R_equiv <- equivalent_circle_radius(cell_side)

safe_fit_variogram <- function(vg_emp, sigma2_obs) {
  starts <- list(
    gstat::vgm(psill = 0.50 * sigma2_obs, model = "Sph", range = 900,  nugget = 0.20 * sigma2_obs),
    gstat::vgm(psill = 0.70 * sigma2_obs, model = "Sph", range = 1100, nugget = 0.10 * sigma2_obs),
    gstat::vgm(psill = 0.40 * sigma2_obs, model = "Sph", range = 600,  nugget = 0.30 * sigma2_obs)
  )

  fits <- lapply(starts, function(v0) try(gstat::fit.variogram(vg_emp, v0), silent = TRUE))
  fits <- fits[!vapply(fits, inherits, logical(1), what = "try-error")]
  if (length(fits) == 0) return(NULL)

  valid <- vapply(fits, function(fit) {
    any(fit$model == "Sph") && all(is.finite(fit$psill)) && all(fit$psill >= 0)
  }, logical(1))
  fits <- fits[valid]
  if (length(fits) == 0) return(NULL)

  sse <- vapply(fits, function(fit) {
    x <- attr(fit, "SSErr")
    if (is.null(x) || !is.finite(x)) Inf else x
  }, numeric(1))
  fits[[which.min(sse)]]
}

run_meuse_metal <- function(metal,
                            study_label,
                            variable_label,
                            units_label = "log(mg kg^-1)") {
  if (!metal %in% names(meuse_sf)) stop("Column not found in meuse data: ", metal)

  dat <- meuse_sf
  response_col <- paste0("log_", metal)
  dat[[response_col]] <- log(dat[[metal]])

  for (v in intersect(c("ffreq", "soil"), names(dat))) {
    dat[[v]] <- as.factor(dat[[v]])
  }
  dat_df <- sf::st_drop_geometry(dat)

  form <- stats::as.formula(paste(response_col, "~", paste(available_predictors, collapse = " + ")))

  # Leave-one-out cross-validation for a simple linear model.
  pred <- rep(NA_real_, nrow(dat_df))
  for (i in seq_len(nrow(dat_df))) {
    fit_i <- stats::lm(form, data = dat_df[-i, , drop = FALSE])
    pred[i] <- stats::predict(fit_i, newdata = dat_df[i, , drop = FALSE])
  }

  obs <- dat_df[[response_col]]
  met <- validation_metrics(pred, obs)

  sigma2_obs <- stats::var(obs, na.rm = TRUE)
  mse_obs <- unname(met["mse"])
  r_obs <- unname(met["r"])
  r2_pred <- 1 - mse_obs / sigma2_obs
  r2_cor <- unname(met["r2"])

  # Fit a response semivariogram. This is a semivariogram of the response
  # variable itself, not of model residuals. It supplies the support-variance
  # ratio for reinterpretation of point-reference validation at 100 m support.
  vg_form <- stats::as.formula(paste(response_col, "~ 1"))
  vg_emp <- gstat::variogram(vg_form, dat)
  vg_fit <- safe_fit_variogram(vg_emp, sigma2_obs)

  gamma_R <- lambda <- sigma2_delta <- sigma2_block <- true_nugget <- partial_sill <- NA_real_
  rho_R <- rho_struct_R <- variogram_range <- nugget_to_sill <- NA_real_
  variogram_model <- NA_character_

  if (!is.null(vg_fit) && any(vg_fit$model == "Sph")) {
    true_nugget <- sum(vg_fit$psill[vg_fit$model == "Nug"], na.rm = TRUE)
    sph <- vg_fit[vg_fit$model == "Sph", ][1, ]
    partial_sill <- sph$psill
    variogram_range <- sph$range
    total_sill <- true_nugget + partial_sill
    nugget_to_sill <- if (is.finite(total_sill) && total_sill > 0) true_nugget / total_sill else NA_real_

    # In Tier 3, measurement error is unknown and set to zero. The fitted
    # nugget is therefore retained as true unresolved microscale variation.
    gamma_R <- true_nugget + partial_sill * spherical_shape(R_equiv, variogram_range)

    # Use the empirical response variance for normalization; the fitted
    # semivariogram supplies the nugget and short-lag structural increment.
    # A fitted-sill/empirical-variance difference is consequently retained as a
    # component that is effectively smooth at the target support.
    terms <- equal_area_disk_linear_terms(
      sigma2_point_true = sigma2_obs,
      sigma2_eps = 0,
      gamma_true_R = gamma_R,
      true_nugget = true_nugget
    )
    lambda <- unname(terms["lambda"])
    sigma2_delta <- unname(terms["sigma2_delta"])
    sigma2_block <- unname(terms["sigma2_block"])
    rho_R <- unname(terms["rho_R"])
    rho_struct_R <- unname(terms["rho_struct_R"])
    variogram_model <- paste(vg_fit$model, collapse = "+")
  }

  res <- make_support_result(
    example_id = paste0("example_meuse_log", metal, "_lm_100m"),
    source = "sp::meuse / sp::meuse.grid, linear model analysed with sf",
    variable = variable_label,
    units = units_label,
    r_obs = r_obs,
    r2_obs = r2_pred,
    rmse_obs = unname(met["rmse"]),
    sigma2_obs = sigma2_obs,
    sigma2_eps = 0,
    lambda = lambda,
    sigma2_delta = sigma2_delta,
    sigma2_block = sigma2_block,
    r2_definition = "predictive_skill",
    correction_type = "tier3_support_sensitivity_bulk_samples_as_points",
    notes = paste(
      "LOOCV linear model using predictors available in meuse.grid:",
      paste(available_predictors, collapse = ", "),
      ". Predictive R2 is computed as 1 - MSE_LOOCV / Var(response) and corrected on the MSE scale.",
      "For audit, squared-correlation R2 from LOOCV predictions is also stored and adjusted through lambda.",
      "Measurement error is unknown and set to zero. The observations are bulk soil samples with finite support but are treated as points here for illustration."
    )
  )

  res <- add_example_metadata(
    res,
    study = study_label,
    prediction_support = paste0(cell_side, " m grid cell"),
    reference_support = "treated as point soil sample",
    information_tier = ifelse(is.finite(lambda), "Tier 3", "Tier 0/3"),
    available_information = "LOOCV predictive R2, response variance and fitted response semivariogram",
    main_limitation = "Bulk soil samples treated as point support for illustration."
  )

  # Additional audit fields for comparing predictive R2 and squared-correlation R2.
  r2_cor_B_raw <- if (is.finite(r2_cor) && is.finite(lambda) && lambda > 0) r2_cor / lambda else NA_real_
  r2_cor_B_noninformative <- is.finite(r2_cor_B_raw) && (r2_cor_B_raw < 0 || r2_cor_B_raw > 1)
  r2_cor_B <- if (is.finite(r2_cor_B_raw) && !r2_cor_B_noninformative) r2_cor_B_raw else NA_real_

  res$cell_side <- cell_side
  res$equivalent_radius <- R_equiv
  res$gamma_true_R <- gamma_R
  res$r2_correlation <- r2_cor
  res$r2_correlation_tier1 <- r2_cor
  res$r2_correlation_tier2 <- r2_cor_B
  res$r2_correlation_tier2_raw <- r2_cor_B_raw
  res$r2_correlation_tier2_noninformative <- r2_cor_B_noninformative
  res$r2_predictive <- r2_pred
  res$variogram_model <- variogram_model
  res$structural_covariance_model <- if (is.finite(variogram_range)) "Spherical" else NA_character_
  res$variogram_range <- variogram_range
  res$true_nugget <- true_nugget
  res$partial_sill <- partial_sill
  res$nugget_to_sill <- nugget_to_sill
  res$rho_R <- rho_R
  res$rho_struct_R <- rho_struct_R
  res
}

rows <- list(
  run_meuse_metal("zinc", "Meuse--Zn", "Log zinc"),
  run_meuse_metal("cadmium", "Meuse--Cd", "Log cadmium")
)

out <- rbind_fill_base(rows)

# Save one CSV for both Meuse rows. The file name is kept compatible with the
# earlier log-zinc-only version so that rerunning with OVERWRITE=TRUE replaces
# the stale one-row output.
save_example_row(out, overwrite = OVERWRITE, verbose = VERBOSE)
print(out)
