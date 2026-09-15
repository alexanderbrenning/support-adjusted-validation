# Brenning et al. (2005): bottom temperature of the winter snow cover (BTS)
#
# Dataset label used in manuscript figures/tables:
#   BTS
#
# Role in manuscript:
# - Tier 2 example with published regression performance and a published
#   residual semivariogram for the Zermatt, Swiss Alps, BTS data set.
# - The measurement-error variance is not available for Zermatt. We therefore
#   transfer only the measurement-error component reported in Brenning et al.'s
#   Table 1 for related meso-scale BTS data from Jotunheimen, southern Norway,
#   originally from Oedegaard et al. (1999).
# - The published performance measure is ordinary-regression R2; this example is
#   therefore treated as a squared-correlation R2_cor example, not as predictive
#   R2_pred.
#
# Regional setting and response:
# - Response variable: bottom temperature of the winter snow cover (BTS), in deg C.
# - Main data set used for the validation metric and residual semivariogram:
#   Zermatt, Swiss Alps; high-mountain permafrost terrain.
# - BTS is used as an empirical proxy for subsurface thermal/permafrost
#   conditions, but the quantitative example uses BTS as the continuous response
#   variable, not a classified permafrost-presence response.
#
# Values taken from Brenning et al. (2005), Zermatt data set:
# - Regression model (2): BTS = 11.4 - 0.0061 Z + 0.12 PSWR.
# - Published ordinary-regression R2 = 0.451 and AIC = 485.99.
# - Residual semivariogram for model (2): range of spatial autocorrelation
#   = 180 m and nugget effect = 0.45 deg C^2.
# - Residual total sill used here = 3.3 deg C^2. This is the variance assigned
#   to the residual component in the support calculation.
# - Observed response variance used for the pooled validation metric = 4.19
#   deg C^2, obtained as the sample variance of the reference observations
#   from the original data set. This replaces the earlier reconstruction
#   sigma_obs^2 = 3.3 / (1 - 0.451).
# - Prediction target support = 30 m grid cell, consistent with the paper's
#   recommendation that meso-scale BTS prediction should not be finer than
#   about 20--30 m because shorter-distance variability is dominated by local
#   effects.
#
# Transferred measurement-error component:
# - Brenning et al. (2005), Table 1, reports for related meso-scale BTS data
#   from Jotunheimen, southern Norway, after Oedegaard et al. (1999):
#     local variability variance      = 0.41 deg C^2,
#     measurement-error variance      = 0.05 deg C^2,
#     range of autocorrelation        = 200 m.
# - Zermatt has local variability/nugget variance = 0.45 deg C^2 but no
#   separately reported measurement-error variance. We therefore use
#   sigma_eps^2 = 0.05 deg C^2 as a transferred component and treat the
#   remaining Zermatt nugget variance, 0.45 - 0.05 = 0.40 deg C^2, as true
#   unresolved microscale variability.
#
# Support-adjustment convention:
# - The updated support-adjustment framework uses support variances,
#   lambda = Var(Z_B) / Var(Z_A^obs), not a fixed point-to-block covariance.
# - The equal-area disk approximation is applied only to the spatially
#   structured residual covariance component. The true microscale/nugget
#   component contributes to reference-support variance and to Delta_{A,B},
#   but not to an ideal continuous 30 m block average.
# -----------------------------------------------------------------------------

source("R/examples/examples_functions.R")

OVERWRITE <- TRUE
VERBOSE <- 2

# Published ordinary-regression R2 for Brenning et al. (2005), model (2).
r2_obs <- 0.451
r_obs <- sqrt(r2_obs)

# Published/compiled Zermatt residual semivariogram quantities.
residual_sill_obs <- 3.3       # deg C^2; residual total sill used in this example
residual_nugget_obs <- 0.45    # deg C^2; Zermatt nugget/local variability
spherical_range <- 180         # m; range of spatial autocorrelation for model (2)

# Measurement-error component transferred from the Norway row in Table 1 of
# Brenning et al. (2005), originally from Oedegaard et al. (1999).
sigma2_eps <- 0.05             # deg C^2

# Decompose the observed Zermatt nugget into transferred measurement error and
# true unresolved microscale variation. Measurement error is not added on top of
# the Zermatt nugget; it is treated as one component of that nugget.
residual_nugget_true <- residual_nugget_obs - sigma2_eps
if (residual_nugget_true < 0) {
  stop("Transferred sigma2_eps exceeds the observed Zermatt residual nugget.")
}

# Spatially structured residual partial sill.
partial_sill <- residual_sill_obs - residual_nugget_obs
if (partial_sill <= 0) stop("Residual partial sill must be positive.")

# Target support: 30 m grid cell.
cell_side <- 30
R <- equivalent_circle_radius(cell_side)

# Observed RMSE is reconstructed from the residual total sill.
mse_obs <- residual_sill_obs
rmse_obs <- sqrt(mse_obs)

# Observed response variance used for q_epsilon and lambda. This is the sample
# variance of the reference observations extracted from the original Zermatt BTS
# data, while the residual semivariogram quantities remain those reported in
# Brenning et al. (2005).
sigma2_obs <- 4.19

# Smooth/fitted component variance implied by the observed response variance and
# residual sill.
# The trend component is treated as effectively smooth at 30 m support.
trend_var <- sigma2_obs - residual_sill_obs
if (trend_var < 0) stop("Implied trend variance is negative.")

# Structural semivariogram shape at the equal-area block radius.
shape_R <- spherical_shape(R, spherical_range)
gamma_resid_struct_R <- partial_sill * shape_R
rho_resid_struct_R <- 1 - gamma_resid_struct_R / partial_sill

# Full latent residual semivariance at R, including the true microscale nugget
# but excluding measurement error. This is retained as an audit quantity.
gamma_resid_true_R <- spherical_semivariogram(
  h = R,
  true_nugget = residual_nugget_true,
  partial_sill = partial_sill,
  range = spherical_range,
  sigma2_eps = 0,
  observed = FALSE
)
rho_resid_true_R <- 1 - gamma_resid_true_R / (residual_nugget_true + partial_sill)

# Equal-area disk linear approximation for the block support.
k_disk <- 128 / (45 * pi)
sigma2_delta_resid_struct <- k_disk * gamma_resid_struct_R

# Reference-to-block support variance difference for the residual component.
# It contains the true microscale nugget plus the block-averaging loss from the
# structural residual covariance.
sigma2_delta <- residual_nugget_true + sigma2_delta_resid_struct

# Residual block variance and total block variance. The true microscale nugget
# and measurement error do not contribute to an ideal continuous block average.
var_block_resid <- partial_sill - sigma2_delta_resid_struct
sigma2_block <- trend_var + var_block_resid
lambda <- sigma2_block / sigma2_obs

res <- make_support_result(
  example_id = "example_brenning_2005_bts",
  source = "Brenning et al. (2005), Zermatt BTS regression model; measurement-error component transferred from Norway BTS data in Table 1",
  variable = "BTS",
  units = "deg C",
  r_obs = r_obs,
  r2_obs = r2_obs,
  rmse_obs = rmse_obs,
  r2_definition = "squared_correlation",
  sigma2_obs = sigma2_obs,
  sigma2_eps = sigma2_eps,
  lambda = lambda,
  sigma2_delta = sigma2_delta,
  sigma2_block = sigma2_block,
  correction_type = "tier2_measurement_error_plus_semivariogram",
  notes = paste(
    "Zermatt, Swiss Alps BTS model (2): R2=0.451, observed response variance=4.19,",
    "residual sill=3.3, residual nugget=0.45 and range=180 m. Measurement-error variance=0.05",
    "is transferred from related meso-scale BTS data from Jotunheimen, southern",
    "Norway, reported in Brenning et al. (2005) Table 1 after Oedegaard et al. (1999).",
    "Target support is a 30 m grid cell."
  )
)

res <- add_example_metadata(
  res,
  study = "BTS",
  prediction_support = "30 m grid cell",
  reference_support = "BTS point observation",
  information_tier = "Tier 2",
  available_information = paste(
    "ordinary-regression R2, residual semivariogram, transferred",
    "measurement-error variance"
  ),
  main_limitation = paste(
    "Measurement-error variance is not Zermatt-specific; it is transferred",
    "from related Norway BTS data reported in the source study."
  )
)

# Audit fields documenting the provenance and decomposition of the calculation.
res$regional_setting <- "Zermatt, Swiss Alps"
res$measurement_error_regional_source <- "Jotunheimen, southern Norway; Oedegaard et al. (1999), as reported in Brenning et al. (2005) Table 1"
res$model_label <- "Brenning et al. (2005) model (2): BTS = 11.4 - 0.0061 Z + 0.12 PSWR"
res$cell_side <- cell_side
res$equivalent_radius <- R
res$residual_sill_obs <- residual_sill_obs
res$residual_nugget_obs <- residual_nugget_obs
res$residual_nugget_true <- residual_nugget_true
res$partial_sill <- partial_sill
res$structural_covariance_model <- "Spherical"
res$spherical_range <- spherical_range
res$sigma2_obs_source <- "sample variance of original Zermatt BTS reference observations"
res$trend_var <- trend_var
res$gamma_resid_struct_R <- gamma_resid_struct_R
res$rho_resid_struct_R <- rho_resid_struct_R
res$gamma_resid_true_R <- gamma_resid_true_R
res$rho_resid_true_R <- rho_resid_true_R
res$sigma2_delta_resid_struct <- sigma2_delta_resid_struct
res$var_block_resid <- var_block_resid

save_example_row(res, overwrite = OVERWRITE, verbose = VERBOSE)
print(res)
