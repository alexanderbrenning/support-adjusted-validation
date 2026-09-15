# -----------------------------------------------------------------------------
# S1-SM: Balenzano et al. (2021), Sentinel-1 soil moisture
# -----------------------------------------------------------------------------
# Purpose:
#   Reported spatial-representativeness-error (SRE) adjustment example.
#
# Dataset label used in manuscript figures/tables:
#   S1-SM
#
# Source study:
#   Balenzano et al. (2021), "Sentinel-1 soil moisture at 1 km resolution:
#   a validation study".
#
# Response variable:
#   Surface volumetric soil moisture [m^3 m^-3], retrieved from Sentinel-1
#   synthetic-aperture radar observations.
#
# Supports:
#   Prediction/product support:
#     Sentinel-1 soil-moisture product at approximately 1 km resolution.
#   Reference support:
#     In situ point-scale soil-moisture measurements at 0.05 m depth.
#
# Values reused from Balenzano et al. (2021), Table 3, overall 1 km validation:
#   RMSE             = 0.088 m^3 m^-3
#   ubRMSE           = 0.085 m^3 m^-3
#   bias             = 0.021 m^3 m^-3
#   delta_SRE        = 0.053 m^3 m^-3
#   intrinsic RMSE   = 0.070 m^3 m^-3
#   intrinsic ubRMSE = 0.067 m^3 m^-3
#   OLS correlation  = 0.46          # descriptive audit field only
#   WLS correlation  = 0.54          # descriptive audit field only
#   number of pairs  = 15057
#   number outliers  = 82
#
# Calculation performed here:
#   The source study explicitly estimates the spatial representativeness error
#   caused by comparing a 1 km product against point-scale in situ data. We do
#   not fit a covariance model and do not re-adjust the reported correlations.
#   We only reproduce the reported intrinsic RMSE by subtracting the SRE
#   variance from the conventional RMSE:
#
#     RMSE_intrinsic = sqrt(RMSE_obs^2 - delta_SRE^2)
#
#   This is not our Tier 2 covariance formula; it is a direct reported-SRE
#   example showing that a support/representativeness-error term can materially
#   affect the RMSE scale.
# -----------------------------------------------------------------------------

OVERWRITE <- TRUE
VERBOSE <- 2

source("R/examples/examples_functions.R")

# Overall 1 km validation over the experimental sites, Table 3 in Balenzano et
# al. (2021). Units are volumetric soil moisture [m^3 m^-3].
rmse_obs <- 0.088
ubrmse_obs <- 0.085
bias_obs <- 0.021
sigma_sre <- 0.053
rmse_intrinsic_reported <- 0.070
ubrmse_intrinsic_reported <- 0.067
r_ols <- 0.46
r_wls <- 0.54
n_pairs <- 15057
n_outliers <- 82

# Reproduce the reported intrinsic RMSE by subtracting SRE variance. Do not use
# this as a covariance-model Tier 2 adjustment; Balenzano et al. model the
# heteroscedastic uncertainty structure directly.
rmse_intrinsic_recomputed <- if (rmse_obs^2 >= sigma_sre^2) {
  sqrt(rmse_obs^2 - sigma_sre^2)
} else {
  NA_real_
}

res <- data.frame(
  example_id = "example_balenzano_2021_sentinel1_soil_moisture",
  source = "Balenzano et al. (2021), Sentinel-1 soil moisture validation, Table 3",
  variable = "Soil moisture",
  units = "m^3 m^-3",
  correction_type = "reported_spatial_representativeness_error_correction",
  r2_definition = "not_used_reported_sre_example",

  # Reported correlations are retained as audit fields, not corrected.
  r_obs = NA_real_,
  r2_obs = NA_real_,
  r_ols = r_ols,
  r_wls = r_wls,

  # Conventional and reported-SRE RMSE quantities.
  rmse_obs = rmse_obs,
  ubrmse_obs = ubrmse_obs,
  bias_obs = bias_obs,
  sigma_sre = sigma_sre,
  sigma2_delta = sigma_sre^2,
  rmse_intrinsic_reported = rmse_intrinsic_reported,
  ubrmse_intrinsic_reported = ubrmse_intrinsic_reported,

  # No measurement-error-only Tier 1 adjustment is computed.
  sigma2_obs = NA_real_,
  sigma2_eps = 0,
  q_eps = NA_real_,
  r_ME = NA_real_,
  R2_ME = NA_real_,
  RMSE_ME = NA_real_,

  # Store the reproduced intrinsic RMSE in the Tier 2/RMSE_B column so that the
  # existing table and figure machinery can display the support-adjusted RMSE
  # value for reported-SRE examples.
  lambda = NA_real_,
  sigma2_block = NA_real_,
  r_B = NA_real_,
  R2_B = NA_real_,
  RMSE_B = rmse_intrinsic_recomputed,

  n_pairs = n_pairs,
  n_outliers = n_outliers,
  correction_status =
    "reported SRE: intrinsic RMSE reproduced by subtracting delta_SRE^2 from RMSE^2",
  notes = paste(
    "S1-SM label. Directly uses Balenzano et al. (2021), Table 3:",
    "RMSE, ubRMSE, bias, delta_SRE, intrinsic RMSE, intrinsic ubRMSE,",
    "and OLS/WLS correlations. No covariance model is fitted and R/R2 are",
    "not re-corrected because the source study treats heteroscedastic",
    "uncertainty in both retrieved and observed variables directly."
  ),
  stringsAsFactors = FALSE
)

res <- add_example_metadata(
  res,
  study = "S1-SM",
  prediction_support = "Sentinel-1 soil moisture product, approximately 1 km",
  reference_support = "in situ point-scale soil moisture at 0.05 m depth",
  information_tier = "Reported SRE",
  available_information = "reported RMSE, ubRMSE, bias, SRE, intrinsic RMSE, and OLS/WLS correlations",
  main_limitation = "Direct SRE example; not a covariance-model post-hoc adjustment."
)

res$source_study <- "Balenzano et al. (2021)"
res$table_source <- "Balenzano et al. (2021), Table 3, overall 1 km validation"
res$support_note <- "Prediction support approximately 1 km; reference support point-scale in situ measurements at 0.05 m depth"

save_example_row(res, overwrite = OVERWRITE, verbose = VERBOSE)
print(res)
