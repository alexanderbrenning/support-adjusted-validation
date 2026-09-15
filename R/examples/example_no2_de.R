# NO2-DE: annual mean NO2 in Germany, 2018
#
# Manuscript source:
#   Brenning, A. & Suesse, T. (2026). Support-adjusted validation metrics
#   in environmental prediction: a minimal-information framework.
#
# Data sources:
#   Umweltbundesamt station metadata and annual NO2 observations for 2018.
#   See data/DATA_PROVENANCE.md for download/provenance notes.
#
# Role in the manuscript:
#   Tier 3 support-sensitivity example. A simple 10-fold cross-validated
#   station-level linear model supplies the raw reference-support metric. A
#   response semivariogram fitted to station annual means supplies the support
#   terms for a 2 km prediction support. Measurement error is set to zero and
#   the fitted nugget is retained as true microscale variation.

source("R/examples/examples_functions.R")
library(sf)
library(gstat)

OVERWRITE <- TRUE
VERBOSE <- 2

sample_file <- "data/no2_de/no2_2018.gpkg"
if (!file.exists(sample_file)) stop("NO2-DE file not found: ", sample_file)

prediction_cell_side_m <- 2000
prediction_radius_m <- equivalent_circle_radius(prediction_cell_side_m)

sample_sf <- sf::st_read(sample_file, quiet = TRUE)
needed <- c("annual_mean", "elevation", "urban", "population", "corine")
missing <- setdiff(needed, names(sample_sf))
if (length(missing)) stop("NO2-DE data are missing required columns: ", paste(missing, collapse = ", "))

sample_sf$.row_id <- seq_len(nrow(sample_sf))
sample_sf$z <- sample_sf$annual_mean
coords <- sf::st_coordinates(sample_sf)
sample_sf$x <- coords[, 1]
sample_sf$y <- coords[, 2]

predictor_vars <- c("elevation", "urban", "population", "corine")
dat <- sf::st_drop_geometry(sample_sf)
dat <- dat[stats::complete.cases(dat[, c("z", predictor_vars)]), , drop = FALSE]
form <- stats::as.formula(paste("z ~", paste(predictor_vars, collapse = " + ")))

# Simple random 10-fold CV. This is an illustrative validation calculation,
# not a reproduction of the TWCV case study.
set.seed(2018)
k <- 10
fold <- sample(rep(seq_len(k), length.out = nrow(dat)))
pred <- rep(NA_real_, nrow(dat))
for (j in seq_len(k)) {
  fit <- stats::lm(form, data = dat[fold != j, , drop = FALSE])
  pred[fold == j] <- stats::predict(fit, newdata = dat[fold == j, , drop = FALSE])
}

met <- validation_metrics(pred, dat$z)
sigma2_obs <- stats::var(dat$z, na.rm = TRUE)

sample_sf2 <- sample_sf[sample_sf$.row_id %in% dat$.row_id, , drop = FALSE]
sample_sf2 <- sample_sf2[match(dat$.row_id, sample_sf2$.row_id), , drop = FALSE]
vg_emp <- gstat::variogram(z ~ 1, sample_sf2)
vg0 <- gstat::vgm(psill = 0.5 * sigma2_obs, model = "Sph", range = 50000,
                  nugget = 0.1 * sigma2_obs)
vg_fit <- gstat::fit.variogram(vg_emp, vg0)

true_nugget <- sum(vg_fit$psill[vg_fit$model == "Nug"], na.rm = TRUE)
sph <- vg_fit[vg_fit$model == "Sph", ][1, ]
partial_sill <- sph$psill
sph_range <- sph$range
gamma_R <- true_nugget + partial_sill * spherical_shape(prediction_radius_m, sph_range)

# The empirical station variance remains the point-reference variance used to
# normalize q_epsilon and lambda. The fitted response semivariogram supplies the
# nugget and the short-lag structural increment. Any small difference between
# the empirical variance and the fitted total sill is therefore retained as a
# component that is effectively smooth at the 2 km target support.
terms <- equal_area_disk_linear_terms(
  sigma2_point_true = sigma2_obs,
  sigma2_eps = 0,
  gamma_true_R = gamma_R,
  true_nugget = true_nugget
)

res <- make_support_result(
  example_id = "example_no2_de",
  source = "UBA annual mean NO2 Germany 2018; simple 10-fold CV linear model",
  variable = "Annual mean NO2",
  units = "ug m^-3",
  r_obs = unname(met["r"]),
  r2_obs = unname(met["r2"]),
  rmse_obs = unname(met["rmse"]),
  sigma2_obs = sigma2_obs,
  sigma2_eps = 0,
  lambda = unname(terms["lambda"]),
  sigma2_delta = unname(terms["sigma2_delta"]),
  sigma2_block = unname(terms["sigma2_block"]),
  r2_definition = "squared_correlation",
  correction_type = "tier3_support_sensitivity_measurement_error_zero",
  notes = "NO2-DE uses UBA annual mean NO2 station data for 2018 and a simple 10-fold CV linear model. Measurement error is set to zero; the fitted response-semivariogram nugget is retained as true microscale variation."
)

res <- add_example_metadata(
  res,
  study = "NO2-DE",
  prediction_support = "2 km grid cell",
  reference_support = "monitoring station annual mean",
  information_tier = "Tier 3",
  available_information = "station data; simple 10-fold CV metrics; fitted response semivariogram; sigma_epsilon^2 set to 0",
  main_limitation = "Illustrative validation model; annual-mean measurement error set to zero; fitted nugget not decomposed."
)

res$data_source_observations <- "Umweltbundesamt annual mean NO2 observations, 2018"
res$data_source_stations <- "Umweltbundesamt station metadata, German regulatory monitoring network"
res$validation_design <- "simple random 10-fold cross-validation"
res$predictor_variables <- paste(predictor_vars, collapse = ";")
res$n_validation <- unname(met["n"])
res$cell_side <- prediction_cell_side_m
res$equivalent_radius <- prediction_radius_m
res$measurement_error_assumption <- "sigma_epsilon^2 = 0"
res$structural_covariance_model <- "Spherical"
res$variogram_model <- paste(vg_fit$model, signif(vg_fit$psill, 5), signif(vg_fit$range, 5), sep = ":", collapse = "; ")
res$variogram_nugget <- true_nugget
res$variogram_spherical_psill <- partial_sill
res$variogram_spherical_range_m <- sph_range
res$variogram_range <- sph_range
res$true_nugget <- true_nugget
res$partial_sill <- partial_sill
res$nugget_to_sill <- true_nugget / (true_nugget + partial_sill)
res$gamma_true_R <- gamma_R
res$gamma_struct_R <- unname(terms["gamma_struct_R"])
res$rho_R <- unname(terms["rho_R"])
res$rho_struct_R <- unname(terms["rho_struct_R"])
res$sigma2_delta_struct <- unname(terms["sigma2_delta_struct"])

save_example_row(res, overwrite = OVERWRITE, verbose = VERBOSE)
print(res)
