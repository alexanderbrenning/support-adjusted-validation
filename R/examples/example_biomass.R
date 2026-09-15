# -----------------------------------------------------------------------------
# Biomass: Zandler et al. (2015), remotely sensed dwarf-shrub biomass
# -----------------------------------------------------------------------------
# Purpose:
#   Tier 1 reference-error example for remotely sensed dwarf-shrub biomass.
#
# Dataset label used in manuscript figures/tables:
#   Biomass
#
# Source study:
#   Zandler et al. (2015), "Quantifying dwarf shrub biomass in an arid
#   environment: comparing empirical methods in a high dimensional setting".
#
# Response variable:
#   Total dwarf-shrub biomass (TB), kg ha^-1, Eastern Pamirs, Tajikistan.
#
# Supports:
#   Prediction/product support:
#     Image-pixel support. The source study evaluates predictors derived from
#     Landsat OLI, RapidEye and ASTER GDEM. The Landsat/ASTER predictors are
#     treated here as 30 m support; RapidEye predictors are at 5 m support.
#     We use the 30 m image-pixel support as the relevant map-prediction
#     support for this compact post-hoc example.
#
#   Reference support:
#     60 m x 60 m field site, summarized from two 4 m x 4 m subplots.
#     Zandler et al. describe 137 mapped homogeneous field plots, of which 122
#     were dwarf-shrub field sites. The 60 m side length was chosen from the
#     coarsest sensor resolution and geometric accuracy and increased to allow
#     for GPS deviations.
#
# Values taken from Zandler et al. (2015):
#   Best model / validation metric:
#     LASSO LSRE spatial cross-validation RMSE = 992 kg ha^-1.
#
#   Approximate response variance:
#     The paper states that cross-validation RMSE values from 992 to
#     1742 kg ha^-1 correspond to about 75% to 130% of the biomass standard
#     deviation. We therefore reconstruct the validation-reference SD for the
#     best model as
#       sigma_obs = 992 / 0.75 = 1322.667 kg ha^-1.
#     This is a rounded reconstruction, not a directly tabulated SD.
#
#   Reference-error components:
#     Allometric biomass model contribution:
#       sigma_allometric = 180 kg ha^-1
#       Derived by Zandler et al. from a plant-level allometric RMSE of about
#       70 g, an average of 34 plants per 4 m x 4 m plot, and two plots per
#       field site.
#
#     Within-site subplot-sampling contribution:
#       sigma_within_site = 518 kg ha^-1
#       Derived by Zandler et al. from an ANOVA of plot-level TB observations:
#       4 m x 4 m plot-scale TB varied within 60 m x 60 m field sites with
#       SD = 732 kg ha^-1, giving an approximate precision of 518 kg ha^-1
#       when estimating field-site TB from two plots.
#
# Assumptions made here:
#   - The allometric and within-site sampling components are treated as
#     independent reference-error components and are added on the variance
#     scale.
#   - No measurement-error/support decomposition beyond these reference-error
#     terms is attempted.
#   - No response semivariogram or point-to-block covariance model is available
#     from the source study, so Tier 2 support adjustment is not computed.
#   - The published spatial-CV RMSE is handled as predictive-skill R2 because
#     the metric is based on cross-validation prediction errors.
# -----------------------------------------------------------------------------

source("R/examples/examples_functions.R")

OVERWRITE <- TRUE
VERBOSE <- 2

# -----------------------------------------------------------------------------
# Published validation metric and reconstructed response variance
# -----------------------------------------------------------------------------

# Zandler et al. (2015), Table 3 / performance discussion:
# best-performing model = LASSO LSRE; spatial-CV RMSE = 992 kg ha^-1.
rmse_obs <- 992

# Zandler et al. state that this lower RMSE corresponds to about 75% of the
# biomass SD. Reconstruct the observed/reference variance from that rounded
# ratio; this value is therefore approximate.
rmse_to_sd_ratio <- 0.75
sigma_obs <- rmse_obs / rmse_to_sd_ratio
sigma2_obs <- sigma_obs^2

# Predictive-skill R2 reconstructed from MSE and observed-reference variance.
# This is not a reported ordinary-regression R2; it is the predictive R2 implied
# by RMSE = 0.75 * SD.
r2_obs <- 1 - rmse_obs^2 / sigma2_obs
r_obs <- sqrt(r2_obs)

# -----------------------------------------------------------------------------
# Reference-error variance from field-reference uncertainty components
# -----------------------------------------------------------------------------

# Zandler et al. expected allometric-model uncertainty at the 60 m field-site
# level to be about 180 kg ha^-1 for total biomass.
sigma_allometric <- 180

# Zandler et al. estimated the precision of field-site total-biomass estimates
# from two 4 m x 4 m plots as 518 kg ha^-1, based on within-site plot-scale
# variability of 732 kg ha^-1 within 60 m x 60 m field sites.
sigma_within_site <- 518

# Treat both components as independent reference-error terms and add on the
# variance scale.
sigma2_eps <- sigma_allometric^2 + sigma_within_site^2
sigma_eps <- sqrt(sigma2_eps)

# No semivariogram/covariance information is used in this example.
lambda <- NA_real_
sigma2_delta <- NA_real_
sigma2_block <- NA_real_

# -----------------------------------------------------------------------------
# Tier 1 calculation
# -----------------------------------------------------------------------------

res <- make_support_result(
  example_id = "example_zandler_2015_biomass",
  source = "Zandler et al. (2015), LASSO LSRE dwarf-shrub biomass model",
  variable = "Total biomass",
  units = "kg ha^-1",
  r_obs = r_obs,
  r2_obs = r2_obs,
  rmse_obs = rmse_obs,
  sigma2_obs = sigma2_obs,
  sigma2_eps = sigma2_eps,
  lambda = lambda,
  sigma2_delta = sigma2_delta,
  sigma2_block = sigma2_block,
  r2_definition = "predictive_skill",
  correction_type = "tier1_reference_error_only",
  notes = paste(
    "Tier 1 reference-error calculation only.",
    "Observed RMSE = 992 kg ha^-1 from the LASSO LSRE spatial-CV result.",
    "Observed biomass variance is reconstructed from the statement that this",
    "RMSE is about 75% of the biomass standard deviation.",
    "Reference-error variance combines 180^2 kg^2 ha^-2 allometric uncertainty",
    "and 518^2 kg^2 ha^-2 within-site subplot-sampling uncertainty.",
    "No response semivariogram is reported, so no Tier 2 support adjustment is attempted."
  )
)

res <- add_example_metadata(
  res,
  study = "Biomass",
  prediction_support = "30 m map pixel",
  reference_support = "60 m field site; two 4 m x 4 m subplots",
  information_tier = "Tier 1",
  available_information = paste(
    "Spatial-CV RMSE; approximate response SD reconstructed from RMSE/SD ratio;",
    "allometric and within-site reference-error components"
  ),
  main_limitation = paste(
    "Reference support is a field-site estimate from sparse subplots;",
    "response variance is reconstructed from a rounded RMSE/SD ratio;",
    "no semivariogram is available."
  )
)

# -----------------------------------------------------------------------------
# Audit fields specific to this example
# -----------------------------------------------------------------------------

# Manuscript-facing short label. The current table/figure scripts may still use
# the `study` field for citation mapping; this field is included for the shorter
# non-author label convention used in revised text.
res$dataset_label <- "Biomass"
res$display_label <- "Biomass"
res$source_study <- "Zandler et al. (2015)"

# Published/reconstructed validation quantities.
res$validation_model <- "LASSO LSRE"
res$validation_design <- "spatial cross-validation"
res$rmse_obs_source <- "Zandler et al. (2015), LASSO LSRE spatial-CV RMSE"
res$rmse_to_sd_ratio <- rmse_to_sd_ratio
res$sigma_obs_reconstructed <- sigma_obs
res$sigma2_obs_reconstructed <- sigma2_obs
res$response_variance_source <- "reconstructed from RMSE = 75% of biomass SD"

# Field-reference uncertainty components.
res$sigma_allometric <- sigma_allometric
res$sigma_within_site <- sigma_within_site
res$sigma_eps <- sigma_eps
res$reference_error_source <- paste(
  "Zandler et al. (2015) uncertainty discussion: allometric model error and",
  "within-site subplot-sampling precision"
)

# Support/covariance status.
res$n_field_sites_total <- 137
res$n_dwarf_shrub_sites <- 122
res$field_site_side_m <- 60
res$subplot_side_m <- 4
res$n_subplots_per_site <- 2
res$tier2_status <- "not computed: no response semivariogram or covariance model reported"

save_example_row(res, overwrite = OVERWRITE, verbose = VERBOSE)
print(res)
