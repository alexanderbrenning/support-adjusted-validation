# Core functions for the pedagogical support-adjustment simulation
#
# The simulation follows the manuscript outline:
#   Z_true(x) = S(x) + M(x)
# S is a spatially correlated Gaussian field generated from a spherical
# covariance model. M is true point-scale/microscale variation. Coarse block
# means are computed from S only, representing the continuous-support block
# target; M is treated as point-scale variation that averages out at block
# support. Point observations are contaminated with measurement error.

source("R/core/support_correction_core.R")

#' Standardize numeric vector or matrix values
#'
#' @param x Numeric vector/matrix.
#' @return Object of same shape with mean 0 and variance 1.
standardize <- function(x) {
  z <- as.numeric(x)
  s <- stats::sd(z)
  if (!is.finite(s) || s == 0) stop("Cannot standardize a constant field.")
  (x - mean(z)) / s
}

#' Simulate an approximate Gaussian field with spherical covariance
#'
#' Uses a periodic circulant embedding on a regular grid. The target covariance
#' is spherical with zero nugget and unit partial sill; the simulated field is
#' subsequently standardized. Small negative spectral coefficients caused by
#' discretization/periodization are truncated to zero.
#'
#' @param n Number of grid cells in each dimension.
#' @param range Practical spherical range in grid-cell units.
#' @return n by n matrix with approximately mean 0 and variance 1.
simulate_spherical_field_fft <- function(n, range) {
  if (!is.finite(range) || range <= 0) stop("range must be positive.")

  # Periodic distances on the torus. This gives fast simulation while keeping
  # the covariance model explicit and interpretable.
  idx <- 0:(n - 1)
  dx <- pmin(idx, n - idx)
  h <- sqrt(outer(dx^2, dx^2, "+"))

  cov_kernel <- 1 - spherical_shape(h, range)
  cov_kernel[h >= range] <- 0
  cov_kernel[1, 1] <- 1

  eig <- Re(fft(cov_kernel))
  n_negative <- sum(eig < -1e-8, na.rm = TRUE)
  eig[eig < 0] <- 0

  z_complex <- matrix(stats::rnorm(n * n), n, n) + 1i * matrix(stats::rnorm(n * n), n, n)
  field <- Re(fft(sqrt(eig) * z_complex, inverse = TRUE)) / n
  out <- standardize(field)
  attr(out, "n_negative_eigenvalues") <- n_negative
  out
}

#' Aggregate a fine-grid matrix to block means
#'
#' @param z Fine-grid matrix.
#' @param fine_per_block Number of fine cells along one block side.
#' @return Matrix of block means.
aggregate_to_blocks <- function(z, fine_per_block) {
  n <- nrow(z)
  if (nrow(z) != ncol(z)) stop("z must be a square matrix.")
  if (n %% fine_per_block != 0) stop("fine grid size must be divisible by fine_per_block.")

  nb <- n / fine_per_block
  row_block <- rep(seq_len(nb), each = fine_per_block)
  col_block <- rep(seq_len(nb), each = fine_per_block)

  sums_row <- rowsum(z, group = row_block, reorder = FALSE)
  sums_block <- t(rowsum(t(sums_row), group = col_block, reorder = FALSE))
  sums_block / fine_per_block^2
}

#' Sample one fine-grid point inside selected blocks
#'
#' @param z_point_field Fine-grid matrix of latent true point-support values.
#' @param z_block_field Matrix of latent block means.
#' @param n_sample Number of sampled blocks.
#' @param fine_per_block Number of fine cells along one block side.
#' @return data.frame with block indices, point values, block means, and delta.
sample_points_in_blocks <- function(z_point_field, z_block_field, n_sample, fine_per_block) {
  nb <- nrow(z_block_field)
  all_blocks <- expand.grid(block_row = seq_len(nb), block_col = seq_len(nb))
  n_sample <- min(n_sample, nrow(all_blocks))
  pick <- all_blocks[sample(seq_len(nrow(all_blocks)), n_sample), , drop = FALSE]

  off_row <- sample(seq_len(fine_per_block), n_sample, replace = TRUE)
  off_col <- sample(seq_len(fine_per_block), n_sample, replace = TRUE)

  fine_row <- (pick$block_row - 1) * fine_per_block + off_row
  fine_col <- (pick$block_col - 1) * fine_per_block + off_col

  z_point <- z_point_field[cbind(fine_row, fine_col)]
  z_block <- z_block_field[cbind(pick$block_row, pick$block_col)]

  data.frame(
    block_row = pick$block_row,
    block_col = pick$block_col,
    fine_row = fine_row,
    fine_col = fine_col,
    z_point_true = z_point,
    z_block_true = z_block,
    delta = z_point - z_block
  )
}

#' Orthogonalize and standardize a vector against one or more columns
#'
#' @param x Numeric vector.
#' @param against Matrix/data.frame of variables to regress out.
#' @return Standardized residual vector.
orthogonal_standard <- function(x, against) {
  d <- data.frame(x = x, as.data.frame(against))
  fit <- stats::lm(x ~ ., data = d)
  standardize(stats::residuals(fit))
}

#' Generate predictions for validation sample
#'
#' @param z_block Latent block means at sampled blocks.
#' @param delta Point-to-block deviations at sampled locations.
#' @param target_r2_block Desired block-scale squared correlation for block-targeted predictions.
#' @param relation Either "block_targeted" or "subgrid_informed".
#' @return Numeric vector of block-support predictions.
generate_prediction_sample <- function(z_block,
                                       delta,
                                       target_r2_block,
                                       relation = c("block_targeted", "subgrid_informed")) {
  relation <- match.arg(relation)
  u <- as.numeric(standardize(z_block))

  if (relation == "block_targeted") {
    # Noise is orthogonal to the block mean, so the sample correlation with the
    # block target is approximately the requested value.
    eta <- orthogonal_standard(stats::rnorm(length(u)), data.frame(u = u))
    x <- sqrt(target_r2_block) * u + sqrt(1 - target_r2_block) * eta
    return(as.numeric(standardize(x)))
  }

  # Failure-mode scenario: the prediction contains partial information about
  # the sampled within-block deviation. This can break the conditional lower-
  # value interpretation of the Tier 1 corrected point-support correlation.
  v <- as.numeric(standardize(delta))
  eta <- orthogonal_standard(stats::rnorm(length(u)), data.frame(u = u, v = v))

  beta_delta <- 0.45
  beta_block <- sqrt(target_r2_block)
  beta_noise <- sqrt(max(0.05, 1 - beta_block^2 - beta_delta^2))

  x <- beta_block * u + beta_delta * v + beta_noise * eta
  as.numeric(standardize(x))
}

#' Microscale/nugget fraction from named level
#'
#' @param microscale Character level: none, moderate, or large.
#' @return True microscale variance fraction of the latent point-support sill.
microscale_fraction <- function(microscale) {
  switch(
    as.character(microscale),
    none = 0,
    moderate = 0.20,
    large = 0.50,
    stop("Unknown microscale level: ", microscale)
  )
}

#' Simulate one scenario replicate
#'
#' @param scenario One-row data.frame with simulation-factor values.
#' @param rep Integer replicate index.
#' @param n_blocks_side Number of coarse blocks along one domain side.
#' @param fine_per_block Fine cells along one block side.
#' @param n_sample Number of validation blocks sampled.
#' @param seed Integer base seed.
#' @return One-row data.frame of simulated true, observed and corrected metrics.
simulate_one_replicate <- function(scenario,
                                   rep,
                                   n_blocks_side,
                                   fine_per_block,
                                   n_sample,
                                   seed = 1) {
  scenario_index <- as.integer(scenario$scenario_index)
  set.seed(seed + 100000L * scenario_index + rep)

  n_fine <- n_blocks_side * fine_per_block
  block_size <- fine_per_block

  # Spherical semivariogram parameters in fine-grid units. The practical range
  # is explicit: gamma(h) reaches the structured sill at h = range.
  true_nugget <- microscale_fraction(scenario$microscale)
  partial_sill <- 1 - true_nugget
  true_sill <- true_nugget + partial_sill
  range <- scenario$range_block_ratio * block_size

  q_target <- scenario$q_epsilon
  sigma2_eps <- if (q_target == 0) 0 else q_target / (1 - q_target) * true_sill
  observed_nugget <- true_nugget + sigma2_eps
  observed_sill <- true_sill + sigma2_eps

  # Generate the structured field from the spherical covariance. True microscale
  # variation is generated separately and kept out of the block mean.
  s_unit <- simulate_spherical_field_fft(n = n_fine, range = range)
  n_negative_eigenvalues <- attr(s_unit, "n_negative_eigenvalues")

  s_field <- sqrt(partial_sill) * s_unit
  m_field <- if (true_nugget > 0) {
    matrix(stats::rnorm(n_fine * n_fine, sd = sqrt(true_nugget)), nrow = n_fine)
  } else {
    matrix(0, nrow = n_fine, ncol = n_fine)
  }

  z_point_field <- s_field + m_field

  # The target is an ideal continuous block mean. Under the pure-nugget model,
  # M(x) has zero covariance at non-zero separation and therefore contributes to
  # point support but averages out of the continuous block integral. The regular
  # fine grid numerically integrates only S(x); averaging the simulated m_field
  # over 25 cell centres would instead introduce an artificial discretization-
  # dependent nugget contribution to the target variance.
  z_block_field <- aggregate_to_blocks(s_field, fine_per_block = fine_per_block)

  smp <- sample_points_in_blocks(
    z_point_field = z_point_field,
    z_block_field = z_block_field,
    n_sample = n_sample,
    fine_per_block = fine_per_block
  )

  pred <- generate_prediction_sample(
    z_block = smp$z_block_true,
    delta = smp$delta,
    target_r2_block = scenario$target_r2_block,
    relation = scenario$prediction_relation
  )

  # Scale predictions to the realized block-target scale.
  pred <- mean(smp$z_block_true) + stats::sd(smp$z_block_true) * pred
  z_obs <- smp$z_point_true + stats::rnorm(nrow(smp), sd = sqrt(sigma2_eps))

  # Theoretical Tier-2 ingredients from the known spherical semivariogram.
  block_terms <- square_block_terms_spherical(
    fine_per_block = fine_per_block,
    true_nugget = true_nugget,
    partial_sill = partial_sill,
    range = range,
    sigma2_eps = sigma2_eps
  )

  # Empirical diagnostics from the simulated sample.
  sigma2_obs_emp <- stats::var(z_obs)

  # This sample covariance ratio is a design diagnostic, not the definition of
  # lambda. Under uniform random point sampling within blocks and independent
  # measurement error, its expectation equals Var(Z_B) / Var(Z_obs), the
  # variance-ratio representativeness factor used for the adjustment.
  lambda_emp <- stats::cov(z_obs, smp$z_block_true) / sigma2_obs_emp
  sigma2_delta_emp <- stats::var(smp$delta)

  m_obs <- validation_metrics(pred, z_obs)
  m_block <- validation_metrics(pred, smp$z_block_true)
  m_true_point <- validation_metrics(pred, smp$z_point_true)

  # Corrections use the known DGP ingredients to mimic an analyst with the
  # semivariogram and measurement-error variance available.
  t1 <- tier1_measurement_error(
    r_obs = unname(m_obs["r"]),
    r2_obs = unname(m_obs["r2"]),
    rmse_obs = unname(m_obs["rmse"]),
    sigma2_obs = observed_sill,
    sigma2_eps = sigma2_eps
  )

  t2 <- tier2_support_correction(
    r_obs = unname(m_obs["r"]),
    r2_obs = unname(m_obs["r2"]),
    rmse_obs = unname(m_obs["rmse"]),
    lambda = unname(block_terms["lambda"]),
    sigma2_eps = sigma2_eps,
    sigma2_delta = unname(block_terms["sigma2_delta"])
  )

  data.frame(
    scenario_index = scenario_index,
    scenario_id = scenario$scenario_id,
    scenario_type = scenario$scenario_type,
    rep = rep,
    q_epsilon_target = q_target,
    q_epsilon_theoretical = sigma2_eps / observed_sill,
    q_epsilon_empirical = sigma2_eps / sigma2_obs_emp,
    range_block_ratio = scenario$range_block_ratio,
    microscale = scenario$microscale,
    true_nugget = true_nugget,
    partial_sill = partial_sill,
    true_sill = true_sill,
    sigma2_eps = sigma2_eps,
    observed_nugget = observed_nugget,
    observed_sill = observed_sill,
    spherical_range = range,
    spherical_range_block_units = scenario$range_block_ratio,
    target_r2_block = scenario$target_r2_block,
    prediction_relation = scenario$prediction_relation,
    n_blocks_side = n_blocks_side,
    fine_per_block = fine_per_block,
    n_sample = nrow(smp),
    sigma2_obs_theoretical = observed_sill,
    sigma2_obs_empirical = sigma2_obs_emp,
    sigma2_delta_theoretical = unname(block_terms["sigma2_delta"]),
    sigma2_delta_empirical = sigma2_delta_emp,
    cov_point_block_theoretical = unname(block_terms["cov_point_block"]),
    var_block_theoretical = unname(block_terms["var_block"]),
    lambda_theoretical = unname(block_terms["lambda"]),
    lambda_empirical = lambda_emp,
    n_negative_eigenvalues = n_negative_eigenvalues,
    r_obs = unname(m_obs["r"]),
    R2_obs = unname(m_obs["r2"]),
    RMSE_obs = unname(m_obs["rmse"]),
    r_true_point = unname(m_true_point["r"]),
    R2_true_point = unname(m_true_point["r2"]),
    RMSE_true_point = unname(m_true_point["rmse"]),
    r_block_true = unname(m_block["r"]),
    R2_block_true = unname(m_block["r2"]),
    RMSE_block_true = unname(m_block["rmse"]),
    r_ME = unname(t1["r_ME"]),
    R2_ME = unname(t1["R2_ME"]),
    RMSE_ME = unname(t1["RMSE_ME"]),
    r_B_hat = unname(t2["r_B"]),
    R2_B_hat = unname(t2["R2_B"]),
    RMSE_B_hat = unname(t2["RMSE_B"]),
    lower_value_holds = unname(t1["r_ME"]) <= unname(m_block["r"]),
    stringsAsFactors = FALSE
  )
}
