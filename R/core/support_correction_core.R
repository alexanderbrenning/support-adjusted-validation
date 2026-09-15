# Shared support-adjustment core functions
#
# These functions are used by both compact published examples and the
# pedagogical simulation.

#' Equivalent circular radius for a square cell
#'
#' @param cell_side Numeric. Square-cell side length in arbitrary map units.
#' @return Equal-area disk radius, R = a / sqrt(pi).
equivalent_circle_radius <- function(cell_side) {
  stopifnot(is.numeric(cell_side), length(cell_side) == 1, cell_side > 0)
  cell_side / sqrt(pi)
}

#' Spherical semivariogram shape
#'
#' @param h Lag distance.
#' @param range Practical spherical range. The semivariogram reaches its sill at this distance.
#' @return Unitless semivariogram shape in [0,1].
spherical_shape <- function(h, range) {
  ifelse(
    h <= 0,
    0,
    ifelse(h < range, 1.5 * (h / range) - 0.5 * (h / range)^3, 1)
  )
}

#' Spherical semivariogram with separated true nugget and measurement error
#'
#' For h > 0, the observed semivariogram includes true microscale variation
#' and measurement error in the discontinuity at the origin. The latent true
#' semivariogram excludes measurement error.
#'
#' @param h Lag distance.
#' @param true_nugget True microscale variance, c0_true.
#' @param partial_sill Spatially structured partial sill, c1.
#' @param range Practical spherical range, a.
#' @param sigma2_eps Measurement-error variance. Included only when observed = TRUE.
#' @param observed Logical. Return observed semivariogram including measurement error?
#' @return Semivariance gamma(h).
spherical_semivariogram <- function(h,
                                    true_nugget,
                                    partial_sill,
                                    range,
                                    sigma2_eps = 0,
                                    observed = FALSE) {
  c0 <- true_nugget + if (observed) sigma2_eps else 0
  ifelse(h <= 0, 0, c0 + partial_sill * spherical_shape(h, range))
}

#' Spherical covariance of the latent true process
#'
#' Measurement error is not included in covariance for non-identical locations.
#'
#' @param h Lag distance.
#' @param true_nugget True microscale variance. Contributes only at h = 0.
#' @param partial_sill Spatially structured partial sill.
#' @param range Practical spherical range.
#' @return Latent true covariance C_true(h).
spherical_covariance_true <- function(h, true_nugget, partial_sill, range) {
  spatial <- partial_sill * (1 - spherical_shape(h, range))
  ifelse(h <= 0, true_nugget + partial_sill, spatial)
}

#' Numerical square-block covariance ingredients for spherical covariance
#'
#' Computes support-variance ingredients for a square block by
#' averaging covariance values over a fine regular grid within the block. The
#' support adjustment uses the block variance, not the covariance between a
#' fixed reference point and the block mean. True microscale/nugget variation
#' is treated as point-scale variation and does not contribute to a continuous
#' block average.
#'
#' @param fine_per_block Number of fine cells along each block side.
#' @param true_nugget True microscale variance c0_true.
#' @param partial_sill Spatially structured partial sill c1.
#' @param range Practical spherical range in fine-cell units.
#' @param sigma2_eps Measurement-error variance.
#' @return Named numeric vector with covariance ingredients and lambda.
square_block_terms_spherical <- function(fine_per_block,
                                         true_nugget,
                                         partial_sill,
                                         range,
                                         sigma2_eps = 0) {
  # Fine-cell centres inside one unit block in fine-cell coordinates.
  xy <- expand.grid(
    x = seq_len(fine_per_block) - 0.5,
    y = seq_len(fine_per_block) - 0.5
  )

  dx <- outer(xy$x, xy$x, "-")
  dy <- outer(xy$y, xy$y, "-")
  h <- sqrt(dx^2 + dy^2)

  # Continuous block averaging would remove pure nugget variation. We therefore
  # use only the spatially structured covariance component in block integrals.
  c_struct <- partial_sill * (1 - spherical_shape(h, range))

  # For the metric adjustment we need support variances, not the covariance
  # between a fixed reference point and the block mean.  True microscale
  # nugget variation contributes to point-support variance but not to the
  # continuous block average.
  cov_point_block <- NA_real_
  var_block <- mean(c_struct)
  sigma2_point_true <- true_nugget + partial_sill
  sigma2_obs <- sigma2_point_true + sigma2_eps
  sigma2_delta <- sigma2_point_true - var_block
  lambda <- var_block / sigma2_obs

  c(
    cov_point_block = cov_point_block,
    var_block = var_block,
    sigma2_point_true = sigma2_point_true,
    sigma2_obs_theoretical = sigma2_obs,
    sigma2_delta = sigma2_delta,
    lambda = lambda
  )
}

#' Safe correlation
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @return Pearson correlation, or NA if undefined.
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  if (stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  stats::cor(x[ok], y[ok])
}

#' Validation metrics from predictions and references
#'
#' @param pred Numeric vector of predictions.
#' @param obs Numeric vector of reference values.
#' @return Named numeric vector with bias, MSE, RMSE, r and R2_cor.
validation_metrics <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  err <- pred[ok] - obs[ok]
  r <- safe_cor(pred[ok], obs[ok])
  c(
    bias = mean(err),
    mse = mean(err^2),
    rmse = sqrt(mean(err^2)),
    r = r,
    r2 = r^2,
    n = sum(ok)
  )
}

#' Tier 1 measurement-error adjustment
#'
#' Adjusts observed point-support metrics for independent reference error.
#' These are latent point-support metrics, not full block-support adjustments.
#'
#' @param r_obs Observed correlation.
#' @param r2_obs Observed squared correlation or R2-like metric.
#' @param rmse_obs Observed RMSE.
#' @param sigma2_obs Observed reference variance.
#' @param sigma2_eps Reference-error variance.
#' @return Named numeric vector with q_eps and adjusted metrics.
tier1_measurement_error <- function(r_obs = NA_real_,
                                    r2_obs = NA_real_,
                                    rmse_obs = NA_real_,
                                    sigma2_obs,
                                    sigma2_eps,
                                    r2_definition = c("squared_correlation", "predictive_skill")) {
  r2_definition <- match.arg(r2_definition)
  if (!is.finite(sigma2_obs) || sigma2_obs <= 0) stop("sigma2_obs must be positive.")
  if (!is.finite(sigma2_eps) || sigma2_eps < 0) stop("sigma2_eps must be non-negative.")

  q_eps <- sigma2_eps / sigma2_obs
  if (q_eps >= 1) warning("q_eps >= 1; measurement-error adjustment undefined.")

  r_me <- if (is.finite(r_obs) && q_eps < 1) r_obs / sqrt(1 - q_eps) else NA_real_

  if (r2_definition == "squared_correlation") {
    r2_me_raw <- if (is.finite(r2_obs) && q_eps < 1) r2_obs / (1 - q_eps) else NA_real_
    if (is.finite(r2_me_raw) && (r2_me_raw < 0 || r2_me_raw > 1)) {
      warning("Tier 1 R2 adjustment is outside [0,1]; returning NA for R2_ME.")
      r2_me <- NA_real_
    } else {
      r2_me <- r2_me_raw
    }
  } else {
    if (is.finite(r2_obs) && is.finite(rmse_obs) && q_eps < 1 && rmse_obs^2 >= sigma2_eps) {
      mse_me <- rmse_obs^2 - sigma2_eps
      sigma2_me <- sigma2_obs - sigma2_eps
      r2_me <- 1 - mse_me / sigma2_me
    } else {
      r2_me <- NA_real_
    }
  }

  rmse_me <- if (is.finite(rmse_obs) && rmse_obs^2 >= sigma2_eps) {
    sqrt(rmse_obs^2 - sigma2_eps)
  } else {
    if (is.finite(rmse_obs)) warning("rmse_obs^2 < sigma2_eps; RMSE_ME undefined.")
    NA_real_
  }

  c(q_eps = q_eps, r_ME = r_me, R2_ME = r2_me, RMSE_ME = rmse_me)
}

#' Tier 2 support adjustment from a support-variance representativeness factor
#'
#' @param r_obs Observed reference-support correlation.
#' @param r2_obs Observed squared correlation or R2-like metric.
#' @param rmse_obs Observed reference-support RMSE.
#' @param lambda Support-variance representativeness factor, Var(Z_B) / Var(Z_obs_A).
#' @param sigma2_eps Reference-error variance.
#' @param sigma2_delta Support-variance difference, Delta_A,B = Var(Z_A) - Var(Z_B).
#' @return Named numeric vector with support-adjusted metrics and diagnostic flags.
tier2_support_correction <- function(r_obs = NA_real_,
                                     r2_obs = NA_real_,
                                     rmse_obs = NA_real_,
                                     sigma2_obs = NA_real_,
                                     lambda,
                                     sigma2_eps,
                                     sigma2_delta,
                                     sigma2_block = NA_real_,
                                     r2_definition = c("squared_correlation", "predictive_skill")) {
  r2_definition <- match.arg(r2_definition)
  if (!is.finite(lambda) || lambda <= 0) stop("lambda must be positive.")
  if (!is.finite(sigma2_eps) || sigma2_eps < 0) stop("sigma2_eps must be non-negative.")
  if (!is.finite(sigma2_delta) || sigma2_delta < 0) stop("sigma2_delta must be non-negative.")

  # The updated manuscript treats values beyond the feasible metric range as
  # diagnostic/non-informative, not as estimates truncated to the boundary.
  r_b_raw <- if (is.finite(r_obs)) r_obs / sqrt(lambda) else NA_real_
  r_noninformative <- is.finite(r_b_raw) && abs(r_b_raw) > 1
  r_b <- if (is.finite(r_b_raw) && !r_noninformative) r_b_raw else NA_real_

  mse_b_raw <- if (is.finite(rmse_obs)) rmse_obs^2 - sigma2_eps - sigma2_delta else NA_real_
  mse_noninformative <- is.finite(mse_b_raw) && mse_b_raw < 0
  rmse_b <- if (is.finite(mse_b_raw) && !mse_noninformative) sqrt(mse_b_raw) else NA_real_

  if (r2_definition == "squared_correlation") {
    r2_b_raw <- if (is.finite(r2_obs)) r2_obs / lambda else NA_real_
    r2_noninformative <- is.finite(r2_b_raw) && (r2_b_raw > 1 || r2_b_raw < 0)
    r2_b <- if (is.finite(r2_b_raw) && !r2_noninformative) r2_b_raw else NA_real_
  } else {
    if (is.finite(r2_obs) && is.finite(rmse_obs) && is.finite(sigma2_block) &&
        sigma2_block > 0) {
      r2_b_raw <- 1 - mse_b_raw / sigma2_block
      r2_noninformative <- is.finite(r2_b_raw) &&
        (mse_noninformative || r2_b_raw > 1 || r2_b_raw < 0)
      r2_b <- if (is.finite(r2_b_raw) && !r2_noninformative) r2_b_raw else NA_real_
    } else {
      r2_b_raw <- NA_real_
      r2_noninformative <- FALSE
      r2_b <- NA_real_
    }
  }

  if (isTRUE(mse_noninformative)) {
    warning("sigma2_eps + sigma2_delta exceeds rmse_obs^2; Tier 2 MSE/RMSE adjustment is non-informative for this row.")
  }
  if (isTRUE(r2_noninformative)) {
    warning("Tier 2 R2 adjustment is outside [0,1]; treating it as non-informative for this row.")
  }

  c(r_B = r_b,
    R2_B = r2_b,
    RMSE_B = rmse_b,
    r_B_raw = r_b_raw,
    R2_B_raw = r2_b_raw,
    RMSE_B_raw = if (is.finite(mse_b_raw) && mse_b_raw >= 0) sqrt(mse_b_raw) else NA_real_,
    mse_B_raw = mse_b_raw,
    mse_B = if (is.finite(mse_b_raw) && mse_b_raw >= 0) mse_b_raw else NA_real_,
    r_bound_truncated = as.numeric(r_noninformative),
    r2_bound_truncated = as.numeric(r2_noninformative),
    mse_bound_truncated = as.numeric(mse_noninformative),
    r_noninformative = as.numeric(r_noninformative),
    r2_noninformative = as.numeric(r2_noninformative),
    mse_noninformative = as.numeric(mse_noninformative),
    support_variance_exceeds_observed_mse = as.numeric(mse_noninformative))
}

#' Common result-row constructor
#'
#' @param example_id Short unique identifier.
#' @param source Human-readable source label.
#' @param variable Response-variable label.
#' @param units Response units.
#' @param r_obs Observed correlation.
#' @param r2_obs Observed R2-like metric.
#' @param rmse_obs Observed RMSE.
#' @param sigma2_obs Observed reference variance.
#' @param sigma2_eps Reference-error variance.
#' @param lambda Optional representativeness factor.
#' @param sigma2_delta Optional support-variance difference, Delta_A,B = Var(Z_A) - Var(Z_B).
#' @param correction_type Short adjustment-type label.
#' @param notes Free-text notes.
#' @return One-row data.frame.
make_support_result <- function(example_id,
                                source,
                                variable,
                                units = NA_character_,
                                r_obs = NA_real_,
                                r2_obs = NA_real_,
                                rmse_obs = NA_real_,
                                sigma2_obs,
                                sigma2_eps = 0,
                                lambda = NA_real_,
                                sigma2_delta = NA_real_,
                                sigma2_block = NA_real_,
                                r2_definition = c("squared_correlation", "predictive_skill"),
                                correction_type = "tier1_or_tier2",
                                notes = NA_character_) {
  r2_definition <- match.arg(r2_definition)
  if (r2_definition == "squared_correlation") {
    if (!is.finite(r2_obs) && is.finite(r_obs)) r2_obs <- r_obs^2
    if (!is.finite(r_obs) && is.finite(r2_obs)) r_obs <- sqrt(r2_obs)
  }

  t1 <- tier1_measurement_error(r_obs, r2_obs, rmse_obs, sigma2_obs, sigma2_eps,
                                r2_definition = r2_definition)

  if (is.finite(lambda) && is.finite(sigma2_delta)) {
    t2 <- tier2_support_correction(r_obs, r2_obs, rmse_obs,
                                   sigma2_obs = sigma2_obs,
                                   lambda = lambda,
                                   sigma2_eps = sigma2_eps,
                                   sigma2_delta = sigma2_delta,
                                   sigma2_block = sigma2_block,
                                   r2_definition = r2_definition)
  } else {
    t2 <- c(r_B = NA_real_, R2_B = NA_real_, RMSE_B = NA_real_,
            r_B_raw = NA_real_, R2_B_raw = NA_real_, RMSE_B_raw = NA_real_,
            mse_B_raw = NA_real_, mse_B = NA_real_,
            r_bound_truncated = 0, r2_bound_truncated = 0,
            mse_bound_truncated = 0,
            r_noninformative = 0, r2_noninformative = 0,
            mse_noninformative = 0,
            support_variance_exceeds_observed_mse = 0)
  }

  data.frame(
    example_id = example_id,
    source = source,
    variable = variable,
    units = units,
    correction_type = correction_type,
    r2_definition = r2_definition,
    r_obs = r_obs,
    r2_obs = r2_obs,
    rmse_obs = rmse_obs,
    sigma2_obs = sigma2_obs,
    sigma2_eps = sigma2_eps,
    q_eps = unname(t1["q_eps"]),
    r_ME = unname(t1["r_ME"]),
    R2_ME = unname(t1["R2_ME"]),
    RMSE_ME = unname(t1["RMSE_ME"]),
    lambda = lambda,
    sigma2_delta = sigma2_delta,
    sigma2_block = sigma2_block,
    r_B = unname(t2["r_B"]),
    R2_B = unname(t2["R2_B"]),
    RMSE_B = unname(t2["RMSE_B"]),
    r_B_raw = unname(t2["r_B_raw"]),
    R2_B_raw = unname(t2["R2_B_raw"]),
    RMSE_B_raw = unname(t2["RMSE_B_raw"]),
    mse_B_raw = unname(t2["mse_B_raw"]),
    mse_B = unname(t2["mse_B"]),
    r_bound_truncated = as.logical(unname(t2["r_bound_truncated"])),
    r2_bound_truncated = as.logical(unname(t2["r2_bound_truncated"])),
    mse_bound_truncated = as.logical(unname(t2["mse_bound_truncated"])),
    r_noninformative = as.logical(unname(t2["r_noninformative"])),
    r2_noninformative = as.logical(unname(t2["r2_noninformative"])),
    mse_noninformative = as.logical(unname(t2["mse_noninformative"])),
    support_variance_exceeds_observed_mse = as.logical(unname(t2["support_variance_exceeds_observed_mse"])),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

#' Write an object unless it already exists
#'
#' @param object R object to save.
#' @param file Path to RDS file.
#' @param overwrite Logical. Overwrite an existing file?
#' @param verbose Integer verbosity level.
#' @return TRUE if written, FALSE if skipped.
save_rds_incremental <- function(object, file, overwrite = TRUE, verbose = 1) {
  if (file.exists(file) && !overwrite) {
    if (verbose >= 2) message("Skipping existing file: ", file)
    return(FALSE)
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(object, file = file)
  if (verbose >= 1) message("Wrote: ", file)
  TRUE
}

#' Write CSV unless it already exists
#'
#' @param x data.frame.
#' @param file Path to CSV file.
#' @param overwrite Logical. Overwrite an existing file?
#' @param verbose Integer verbosity level.
#' @return TRUE if written, FALSE if skipped.
write_csv_incremental <- function(x, file, overwrite = TRUE, verbose = 1) {
  if (file.exists(file) && !overwrite) {
    if (verbose >= 2) message("Skipping existing file: ", file)
    return(FALSE)
  }
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(x, file = file, row.names = FALSE)
  if (verbose >= 1) message("Wrote: ", file)
  TRUE
}
