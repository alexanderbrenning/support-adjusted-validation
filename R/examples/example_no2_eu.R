# NO2-EU: European annual mean NO2, Vizcaino and Lavalle (2018)
#
# Manuscript source:
#   Brenning, A. & Suesse, T. (2026). Support-adjusted validation metrics
#   in environmental prediction: a minimal-information framework.
#
# Data sources:
#   Vizcaino and Lavalle (2018), Environmental Pollution, DOI:
#   10.1016/j.envpol.2018.03.075; accompanying Mendeley Data record DOI:
#   10.17632/kkss4z9yrs.1, licensed CC BY 4.0. Coordinates are recovered from
#   historical AirBase/European Environment Agency station metadata. See
#   data/DATA_PROVENANCE.md for attribution and derivative-data details.
#
# Role in the manuscript:
#   Tier 2 sensitivity example. Cross-validated station-level predictive
#   metrics are computed from the coordinate-augmented station table. A response
#   semivariogram is fitted from the recovered coordinates and used for a 100 m
#   prediction support. Two nugget partitions are reported.

source("R/examples/examples_functions.R")
library(sf)
library(gstat)
library(randomForest)

OVERWRITE <- TRUE
VERBOSE <- 2

input_file <- "data/no2_eu/NO2_2010_prox_with_coords.rds"
input_file_csv <- "data/no2_eu/NO2_2010_prox_with_coords.csv"
out_dir <- "results/examples"

response_col <- "statistic_"
lon_col <- "longitude_recovered"
lat_col <- "latitude_recovered"
cell_side <- 100
cv_k <- 10
cv_seed <- 2010
rf_ntree <- 500

predictor_vars <- c(
  "nb_luabloc", "nb_luAgloc2010", "nb_luFrloc", "nb_luGUloc",
  "nb_luIdloc2010", "nb_luUrloc2010", "nb_pop2010",
  "nb_proxAg_empl", "nb_proxAg_gva", "nb_proxId_empl", "nb_proxId_gva",
  "nb_proxUr_empl", "nb_proxUr_gva", "nb_rHloc", "nb_rLloc", "nb_rMloc",
  "nb_toploc", "re_coastal", "re_mdtall", "re_tempall", "re_windspa",
  "r_proxIf2010", "r_proxIn2010", "r_proxUr2010"
)

read_input <- function() {
  if (file.exists(input_file)) return(readRDS(input_file))
  if (file.exists(input_file_csv)) return(utils::read.csv(input_file_csv, stringsAsFactors = FALSE, check.names = FALSE))
  stop("NO2-EU coordinate-augmented input not found: ", input_file,
       " or ", input_file_csv,
       ". Run R/ancillary/recover_no2_eu_airbase_coordinates.R first if needed.")
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

d_raw <- read_input()
missing <- setdiff(c(response_col, lon_col, lat_col, predictor_vars), names(d_raw))
if (length(missing)) stop("NO2-EU input is missing required columns: ", paste(missing, collapse = ", "))

d_raw[[response_col]] <- as_num(d_raw[[response_col]])
d_raw[[lon_col]] <- as_num(d_raw[[lon_col]])
d_raw[[lat_col]] <- as_num(d_raw[[lat_col]])
for (v in predictor_vars) d_raw[[v]] <- as_num(d_raw[[v]])

# Validation data: complete response and predictors.
keep_model <- stats::complete.cases(d_raw[, c(response_col, predictor_vars), drop = FALSE])
d_model <- d_raw[keep_model, , drop = FALSE]
if (nrow(d_model) < 100) stop("Too few complete NO2-EU records for validation.")

form <- stats::as.formula(paste(response_col, "~", paste(predictor_vars, collapse = " + ")))
set.seed(cv_seed)
fold <- sample(rep(seq_len(cv_k), length.out = nrow(d_model)))
pred <- rep(NA_real_, nrow(d_model))
for (j in seq_len(cv_k)) {
  train <- d_model[fold != j, , drop = FALSE]
  test <- d_model[fold == j, , drop = FALSE]
  mtry <- max(1, floor(sqrt(length(predictor_vars))))
  fit <- randomForest::randomForest(
    x = train[, predictor_vars, drop = FALSE],
    y = train[[response_col]],
    ntree = rf_ntree,
    mtry = mtry,
    importance = FALSE
  )
  pred[fold == j] <- stats::predict(fit, newdata = test[, predictor_vars, drop = FALSE])
}

obs <- d_model[[response_col]]
met <- validation_metrics(pred, obs)
sigma2_obs <- stats::var(obs, na.rm = TRUE)
mse_obs <- unname(met["mse"])
rmse_obs <- unname(met["rmse"])
r_obs <- unname(met["r"])
r2_cor <- unname(met["r2"])
r2_pred <- 1 - mse_obs / sigma2_obs

# Response semivariogram from all coordinate-augmented stations with NO2 values.
keep_svgm <- is.finite(d_raw[[response_col]]) & is.finite(d_raw[[lon_col]]) & is.finite(d_raw[[lat_col]])
sf_all <- sf::st_as_sf(d_raw[keep_svgm, , drop = FALSE], coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
sf_all <- sf::st_transform(sf_all, 3035)
sf_all$z <- sf_all[[response_col]]

vg_emp <- gstat::variogram(z ~ 1, sf_all, cressie = TRUE, cutoff = 900000, width = 25000)
vg_emp <- vg_emp[is.finite(vg_emp$gamma) & is.finite(vg_emp$dist), , drop = FALSE]
vg0 <- gstat::vgm(psill = 0.4 * sigma2_obs, model = "Sph", range = 420000,
                  nugget = 0.6 * sigma2_obs)
vg_fit <- gstat::fit.variogram(vg_emp, vg0)

fitted_nugget <- sum(vg_fit$psill[vg_fit$model == "Nug"], na.rm = TRUE)
sph <- vg_fit[vg_fit$model == "Sph", ][1, ]
partial_sill <- sph$psill
variogram_range_m <- sph$range
R_equiv <- equivalent_circle_radius(cell_side)
shape_R <- spherical_shape(R_equiv, variogram_range_m)
gamma_struct_R <- partial_sill * shape_R
k_disk <- 128 / (45 * pi)

nugget_partitions <- data.frame(
  study = c("NO2-EU; 50/50 nugget split", "NO2-EU; full nugget microscale"),
  suffix = c("split", "full_nugget"),
  measurement_error_fraction = c(0.50, 0.00),
  microscale_fraction = c(0.50, 1.00),
  stringsAsFactors = FALSE
)

rows <- lapply(seq_len(nrow(nugget_partitions)), function(i) {
  part <- nugget_partitions[i, , drop = FALSE]
  sigma2_eps <- part$measurement_error_fraction * fitted_nugget
  sigma2_mu <- part$microscale_fraction * fitted_nugget
  sigma2_delta <- sigma2_mu + k_disk * gamma_struct_R

  # The empirical response variance is the denominator for q_epsilon and lambda.
  # The fitted semivariogram supplies the nugget partition and short-lag
  # structural increment used in Delta_A,B; it does not replace sigma2_obs when
  # its fitted total sill differs from the empirical response variance.
  sigma2_block <- sigma2_obs - sigma2_eps - sigma2_delta
  lambda <- sigma2_block / sigma2_obs

  res <- make_support_result(
    example_id = paste0("example_no2_eu_", part$suffix),
    source = "Vizcaino and Lavalle (2018) NO2 Europe station data; recovered AirBase/EEA coordinates",
    variable = "Annual mean NO2",
    units = "ug m^-3",
    r_obs = r_obs,
    r2_obs = r2_pred,
    rmse_obs = rmse_obs,
    sigma2_obs = sigma2_obs,
    sigma2_eps = sigma2_eps,
    lambda = lambda,
    sigma2_delta = sigma2_delta,
    sigma2_block = sigma2_block,
    r2_definition = "predictive_skill",
    correction_type = "tier2_response_variogram_nugget_partition",
    notes = paste(
      "RF 10-fold CV metrics computed from the released station table;",
      "100 m prediction support; fitted response semivariogram; nugget partition:",
      part$study
    )
  )

  res <- add_example_metadata(
    res,
    study = part$study,
    prediction_support = "100 m grid cell",
    reference_support = "AirBase monitoring station",
    information_tier = "Tier 2 sens.",
    available_information = "station data; recovered AirBase/EEA coordinates; response variance; CV RMSE; fitted response semivariogram; nugget partition scenario",
    main_limitation = "Does not reproduce pan-European raster map; fitted nugget decomposition is a sensitivity assumption."
  )

  res$data_source <- "Vizcaino and Lavalle (2018) Mendeley Data (CC BY 4.0); AirBase/EEA metadata"
  res$model_type <- "randomForest"
  res$validation_design <- paste0(cv_k, "-fold cross-validation")
  res$predictor_vars <- paste(predictor_vars, collapse = ";")
  res$n_total_raw <- nrow(d_raw)
  res$n_model_complete <- nrow(d_model)
  res$n_validation <- unname(met["n"])
  res$n_variogram <- nrow(sf_all)
  res$r2_predictive <- r2_pred
  res$r2_correlation <- r2_cor
  res$cell_side <- cell_side
  res$equivalent_radius <- R_equiv
  res$structural_covariance_model <- "Spherical"
  res$variogram_model <- "Sph + Nug"
  res$variogram_range <- variogram_range_m
  res$variogram_range_m <- variogram_range_m
  res$variogram_nugget <- fitted_nugget
  res$partial_sill <- partial_sill
  res$partial_sill_scaled <- partial_sill
  res$true_nugget <- sigma2_mu
  res$true_nugget_used_for_support <- sigma2_mu
  res$measurement_error_nugget <- sigma2_eps
  res$fitted_nugget_scaled <- fitted_nugget
  res$nugget_to_sill <- fitted_nugget / (fitted_nugget + partial_sill)
  res$nugget_measurement_error_fraction <- part$measurement_error_fraction
  res$nugget_support_fraction <- part$microscale_fraction
  res$gamma_struct_R <- gamma_struct_R
  res$sigma2_delta_struct <- k_disk * gamma_struct_R
  res$rho_struct_R <- 1 - gamma_struct_R / partial_sill
  res
})

res <- rbind_fill_base(rows)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(out_dir, "example_no2_eu.csv"), row.names = FALSE)
if (VERBOSE >= 1) message("Wrote: ", file.path(out_dir, "example_no2_eu.csv"))
print(res)
