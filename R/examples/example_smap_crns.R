# -----------------------------------------------------------------------------
# SMAP-CRNS: Zhao et al. (2021), soil moisture
# -----------------------------------------------------------------------------
# Purpose:
#   Partial Tier 1 reference-error calculation for ubRMSD only.
#
# Dataset label used in manuscript figures/tables:
#   SMAP-CRNS
#
# Source study:
#   Zhao et al. (2021), "The importance of subsurface processes in land
#   surface modeling over a temperate region: An analysis with SMAP, cosmic
#   ray neutron sensing and triple collocation analysis".
#
# Response variable:
#   Near-surface volumetric soil moisture [cm^3 cm^-3, numerically identical
#   to m^3 m^-3].
#
# Supports:
#   Prediction/product support:
#     SMAP enhanced soil-moisture product, approximately 9 km pixel support.
#   Reference support:
#     Cosmic-ray neutron sensing (CRNS) footprint.
#     This is an intermediate-area field-scale footprint, not a strict point
#     observation.
#
# Values reused from Zhao et al. (2021):
#   Conventional SMAP L3_SM_E_P--CRNS comparison:
#     r        = 0.699              # descriptive only; not adjusted here
#     RMSD     = 0.085 cm^3 cm^-3   # biased RMSD, retained as audit field
#     ubRMSD   = 0.056 cm^3 cm^-3   # unbiased RMSD corrected below
#
#   Triple-collocation result:
#     CRNS error standard deviations for the CLM-ParFlow triplet are entered
#     below as site-level values from the published TC summary table. The
#     aggregate CRNS reference-error variance used here is the arithmetic mean
#     of the site-specific variances, sigma_eps^2 = mean(sigma_eps_site^2),
#     giving sigma_eps approximately 0.026 cm^3 cm^-3.
#
# Quantities not reported in the tabulated information used here:
#   - variance of the CRNS reference series;
#   - q_epsilon = sigma_eps^2 / sigma_obs^2;
#   - point-scale, CRNS-footprint, or SMAP-pixel response semivariogram;
#   - lambda_A_B for the SMAP-pixel target support.
#
# Consequence:
#   We do not compute correlation, R2, or Tier 2 support adjustments. We only
#   subtract the TC-estimated CRNS reference-error variance from ubRMSD^2:
#
#     ubRMSD_adj = sqrt(ubRMSD_obs^2 - sigma_eps_CRNS^2)
#
#   This is a partial Tier 1 reference-error calculation for ubRMSD, not a full
#   support-corrected validation analysis.
# -----------------------------------------------------------------------------

source("R/examples/examples_functions.R")

OVERWRITE <- TRUE
VERBOSE <- 2

# Conventional SMAP--CRNS validation summary from Zhao et al. (2021).
r_obs <- 0.699
r2_obs <- r_obs^2
RMSD_obs <- 0.085
ubRMSD_obs <- 0.056

# Site-specific TC-estimated CRNS error standard deviations from Zhao et al.
# (2021), CLM-ParFlow triplet. Units are volumetric soil moisture.
sigma_eps_crns_sites <- c(
  0.021, 0.023, 0.011, 0.015, 0.028, 0.040, 0.030,
  0.039, 0.038, 0.021, 0.023, 0.010, 0.023
)

# Aggregate the site-specific CRNS error estimates on the variance scale.
sigma2_eps <- mean(sigma_eps_crns_sites^2)
sigma_eps <- sqrt(sigma2_eps)

# Partial Tier 1 ubRMSD adjustment.  We avoid max(0, .) here: if the published
# ubRMSD were smaller than the estimated reference-error SD, the result would be
# non-informative and should be flagged as NA rather than truncated to zero.
ubrmsd_me <- if (ubRMSD_obs^2 >= sigma2_eps) {
  sqrt(ubRMSD_obs^2 - sigma2_eps)
} else {
  NA_real_
}

res <- data.frame(
  example_id = "example_zhao_2021_tc",
  source = "Zhao et al. (2021), SMAP L3_SM_E_P vs CRNS; TC with CLM-ParFlow",
  variable = "Soil moisture",
  units = "cm^3 cm^-3",
  correction_type = "partial_tier1_reference_error_from_tc",
  r2_definition = NA_character_,

  # Published conventional validation quantities.
  r_obs = r_obs,
  r2_obs = NA_real_,
  rmse_obs = ubRMSD_obs,
  RMSD_obs_biased = RMSD_obs,
  ubRMSD_obs = ubRMSD_obs,

  # Reference-error component reused for partial Tier 1 adjustment.
  sigma2_obs = NA_real_,
  sigma2_eps = sigma2_eps,
  sigma_eps = sigma_eps,
  q_eps = NA_real_,

  # Tier 1 correlation/R2 adjustments cannot be computed because the CRNS
  # reference-series variance is not tabulated.
  r_ME = NA_real_,
  R2_ME = NA_real_,
  RMSE_ME = ubrmsd_me,
  ubRMSD_ME = ubrmsd_me,

  # Tier 2 support adjustment is not attempted because no point-scale or
  # footprint-scale response covariance model is used.
  lambda = NA_real_,
  sigma2_delta = NA_real_,
  sigma2_block = NA_real_,
  r_B = NA_real_,
  R2_B = NA_real_,
  RMSE_B = NA_real_,

  n_tc_sites = length(sigma_eps_crns_sites),
  tc_sigma_eps_min = min(sigma_eps_crns_sites),
  tc_sigma_eps_max = max(sigma_eps_crns_sites),
  correction_status =
    "partial Tier 1: TC-estimated CRNS reference-error variance subtracted from ubRMSD only",
  notes = paste(
    "SMAP-CRNS label. Uses published SMAP--CRNS r, RMSD, and ubRMSD plus",
    "TC-estimated CRNS error standard deviations from the CLM-ParFlow",
    "triplet. CRNS reference variance is not tabulated; q_epsilon,",
    "correlation, R2, and Tier 2 support adjustments are therefore not computed."
  ),
  stringsAsFactors = FALSE
)

res <- add_example_metadata(
  res,
  study = "SMAP-CRNS",
  prediction_support = "SMAP 9 km pixel",
  reference_support = "CRNS footprint",
  information_tier = "Partial Tier 1",
  available_information = "r, ubRMSD, TC-estimated CRNS error SD",
  main_limitation = "Reference variance not tabulated; only ubRMSD can be corrected."
)

res$source_study <- "Zhao et al. (2021)"
res$validation_metric_source <- "Zhao et al. (2021), conventional SMAP L3_SM_E_P--CRNS comparison"
res$reference_error_source <- "Zhao et al. (2021), triple-collocation CRNS error SDs, CLM-ParFlow triplet"
res$support_note <- "SMAP 9 km pixel; CRNS field-scale footprint; no support covariance model used"

save_example_row(res, overwrite = OVERWRITE, verbose = VERBOSE)
print(res)
