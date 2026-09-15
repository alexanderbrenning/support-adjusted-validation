# Recover coordinates for Vizcaino & Lavalle (2018) NO2 station table
#
# Expected input supplied by user:
#   data/no2_eu/NO2_2010_prox.csv
#
# The station code column in this file is expected to be 'station_eu' or
# 'station__1', with AirBase/EoI-style codes such as 'AT0ILL1'.  The script
# attempts exact joins to EEA/AirBase metadata using these codes.
#
# Outputs:
#   data/no2_eu/NO2_2010_prox_with_coords.csv
#   data/no2_eu/NO2_2010_prox_with_coords.rds
#   results/examples/vizcaino_airbase_coordinate_match_report.csv
#   results/examples/vizcaino_airbase_coordinate_match_summary.csv
#   analysis/data/NO2_Europe_with_coordinates.rds        (if NO2_Europe.rds exists and can be joined)
#
# Notes:
# - The primary join key is AirQualityStationEoICode / station_european_code.
# - The script uses exact code matching first. No fuzzy matching is used by
#   default, because station names/cities are less auditable than EoI codes.
# - EEA services occasionally change URLs or field names. The code is written
#   defensively and keeps local cached copies of downloaded metadata.

# -----------------------------
# User settings
# -----------------------------

base_dir_candidates <- c("analysis", ".")
input_csv_rel <- file.path("data", "NO2_2010_prox.csv")
input_rds_rel <- file.path("data", "NO2_Europe.rds")

# Downloaded/cached metadata.
external_dir_rel <- file.path("data", "external_airbase_metadata")

# Current e-Reporting/AirBase bridge metadata. It is broad and may include
# post-2012 station metadata, but it contains AirQualityStationEoICode and
# coordinates for many historical AirBase stations.
paneuropean_metadata_url <- "https://ereporting.blob.core.windows.net/downloadservice/metadata.csv"

# Alternative documented URL used by some EEA examples. Keep as fallback.
paneuropean_metadata_url_alt <- "http://discomap.eea.europa.eu/map/fme/metadata/PanEuropean_metadata.csv"

# Historical AirBase station ArcGIS REST service, AirBase v8. This is useful
# for station-level metadata independent of pollutant time series.
airbase_station_query_url <- paste0(
  "https://air.discomap.eea.europa.eu/arcgis/rest/services/",
  "Airbase/AirBase_Stations_WM/MapServer/1/query"
)

# Set TRUE only if the blob metadata does not work and you want to try the
# ArcGIS REST station service. The blob metadata is usually enough.
try_arcgis_airbase_service <- TRUE

# -----------------------------
# Helpers
# -----------------------------

find_base_dir <- function(base_candidates, rel_path) {
  for (b in base_candidates) {
    p <- file.path(b, rel_path)
    if (file.exists(p)) return(normalizePath(b, mustWork = TRUE))
  }
  stop("Cannot find input file ", rel_path, " under: ", paste(base_candidates, collapse = ", "))
}

norm_code <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  gsub("[^A-Z0-9]", "", x)
}

pick_col <- function(x, candidates) {
  nm <- names(x)
  hit <- candidates[tolower(candidates) %in% tolower(nm)]
  if (length(hit) == 0) return(NA_character_)
  nm[match(tolower(hit[1]), tolower(nm))]
}

read_csv_flex <- function(path, ...) {
  # EEA metadata is comma-separated; user file may be comma-separated too.
  # read.csv handles quoted fields reasonably well.
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

safe_download <- function(url, destfile) {
  ok <- FALSE
  msg <- NULL
  tryCatch({
    utils::download.file(url, destfile, mode = "wb", quiet = TRUE)
    ok <- file.exists(destfile) && file.info(destfile)$size > 0
  }, error = function(e) {
    msg <<- conditionMessage(e)
  })
  if (!ok) warning("Download failed for ", url, if (!is.null(msg)) paste0(": ", msg))
  ok
}

collapse_station_metadata <- function(meta, code_cols, lon_col, lat_col, source_label) {
  if (is.na(lon_col) || is.na(lat_col)) return(NULL)
  out_list <- list()
  for (cc in code_cols) {
    if (!cc %in% names(meta)) next
    tmp <- meta[, unique(c(cc, lon_col, lat_col, names(meta))), drop = FALSE]
    names(tmp)[names(tmp) == cc] <- "station_code_raw"
    names(tmp)[names(tmp) == lon_col] <- "longitude"
    names(tmp)[names(tmp) == lat_col] <- "latitude"
    tmp$station_code_norm <- norm_code(tmp$station_code_raw)
    tmp$longitude <- suppressWarnings(as.numeric(tmp$longitude))
    tmp$latitude <- suppressWarnings(as.numeric(tmp$latitude))
    tmp$coord_source <- source_label
    tmp$key_source_column <- cc
    tmp <- tmp[nzchar(tmp$station_code_norm) & is.finite(tmp$longitude) & is.finite(tmp$latitude), , drop = FALSE]
    if (nrow(tmp) > 0) out_list[[length(out_list) + 1]] <- tmp[, c("station_code_norm", "station_code_raw", "longitude", "latitude", "coord_source", "key_source_column"), drop = FALSE]
  }
  if (length(out_list) == 0) return(NULL)
  out <- do.call(rbind, out_list)
  out <- unique(out)

  # Collapse duplicates by station code. Prefer codes with one unique coordinate;
  # otherwise average duplicate coordinates and mark as ambiguous.
  split_out <- split(out, out$station_code_norm)
  collapsed <- lapply(split_out, function(z) {
    lon_u <- unique(round(z$longitude, 7))
    lat_u <- unique(round(z$latitude, 7))
    data.frame(
      station_code_norm = z$station_code_norm[1],
      station_code_raw = z$station_code_raw[1],
      longitude = mean(z$longitude, na.rm = TRUE),
      latitude = mean(z$latitude, na.rm = TRUE),
      coord_source = paste(unique(z$coord_source), collapse = ";"),
      key_source_column = paste(unique(z$key_source_column), collapse = ";"),
      n_metadata_records = nrow(z),
      n_unique_coord_pairs = length(unique(paste(round(z$longitude, 7), round(z$latitude, 7)))),
      coordinate_ambiguous = length(unique(paste(round(z$longitude, 7), round(z$latitude, 7)))) > 1,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, collapsed)
}

fetch_arcgis_airbase_stations <- function(query_url, cache_path, page_size = 1000) {
  if (file.exists(cache_path)) {
    return(read_csv_flex(cache_path))
  }
  library(jsonlite)
  library(httr)

  all <- list()
  offset <- 0
  repeat {
    resp <- try(httr::GET(
      query_url,
      query = list(
        where = "1=1",
        outFields = "*",
        returnGeometry = "true",
        outSR = "4326",
        f = "json",
        resultOffset = offset,
        resultRecordCount = page_size
      ),
      httr::timeout(120)
    ), silent = TRUE)
    if (inherits(resp, "try-error") || httr::http_error(resp)) {
      warning("ArcGIS REST request failed at offset ", offset)
      break
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    obj <- jsonlite::fromJSON(txt, flatten = TRUE)
    if (!is.null(obj$error)) {
      warning("ArcGIS REST returned error: ", paste(unlist(obj$error), collapse = " | "))
      break
    }

    feats <- obj$features
    if (is.null(feats)) break

    # jsonlite::fromJSON(flatten = TRUE) may return either
    #   attributes.* / geometry.x columns, or
    #   a data frame with list/data-frame columns named attributes and geometry.
    # The previous version assumed feats$attributes was always a data.frame;
    # on some EEA responses it is NULL/list-like, which made nrow(dat) length 0.
    if (is.data.frame(feats) && any(startsWith(names(feats), "attributes."))) {
      dat <- feats[, startsWith(names(feats), "attributes."), drop = FALSE]
      names(dat) <- sub("^attributes\\.", "", names(dat))
      if (all(c("geometry.x", "geometry.y") %in% names(feats))) {
        dat$arcgis_longitude <- feats[["geometry.x"]]
        dat$arcgis_latitude <- feats[["geometry.y"]]
      }
    } else if (is.data.frame(feats) && "attributes" %in% names(feats)) {
      attrs <- feats$attributes
      if (is.data.frame(attrs)) {
        dat <- attrs
      } else if (is.list(attrs) && length(attrs) > 0) {
        dat <- tryCatch(do.call(rbind.data.frame, attrs), error = function(e) NULL)
      } else {
        dat <- NULL
      }
      geom <- if ("geometry" %in% names(feats)) feats$geometry else NULL
      if (!is.null(dat) && is.data.frame(geom)) {
        if (all(c("x", "y") %in% names(geom))) {
          dat$arcgis_longitude <- geom$x
          dat$arcgis_latitude <- geom$y
        }
      }
    } else {
      dat <- NULL
    }

    if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0) break
    all[[length(all) + 1]] <- dat
    if (nrow(dat) < page_size) break
    offset <- offset + page_size
  }

  if (length(all) == 0) return(NULL)
  res <- do.call(rbind, all)
  utils::write.csv(res, cache_path, row.names = FALSE)
  res
}

augment_by_station_code <- function(dat, meta_coords, station_col) {
  dat$station_code_norm <- norm_code(dat[[station_col]])
  m <- match(dat$station_code_norm, meta_coords$station_code_norm)
  dat$longitude_recovered <- meta_coords$longitude[m]
  dat$latitude_recovered <- meta_coords$latitude[m]
  dat$coord_source <- meta_coords$coord_source[m]
  dat$key_source_column <- meta_coords$key_source_column[m]
  dat$n_metadata_records <- meta_coords$n_metadata_records[m]
  dat$n_unique_coord_pairs <- meta_coords$n_unique_coord_pairs[m]
  dat$coordinate_ambiguous <- meta_coords$coordinate_ambiguous[m]
  dat$coordinate_match <- is.finite(dat$longitude_recovered) & is.finite(dat$latitude_recovered)
  dat
}

# -----------------------------
# Main workflow
# -----------------------------

base_dir <- find_base_dir(base_dir_candidates, input_csv_rel)
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "results", "examples")
ext_dir <- file.path(base_dir, external_dir_rel)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)

input_csv <- file.path(base_dir, input_csv_rel)
message("Reading station data: ", input_csv)
no2 <- read_csv_flex(input_csv)

station_col <- pick_col(no2, c("station_eu", "station__1", "station_european_code", "AirQualityStationEoICode"))
if (is.na(station_col)) stop("No station-code column found. Expected station_eu or station__1.")
message("Using station-code column: ", station_col)

country_col <- pick_col(no2, c("country_na", "Country", "country_name", "country"))

# ---- Source 1: Pan-European metadata CSV ----
metadata_cache <- file.path(ext_dir, "eea_paneuropean_metadata.csv")
if (!file.exists(metadata_cache)) {
  ok <- safe_download(paneuropean_metadata_url, metadata_cache)
  if (!ok) ok <- safe_download(paneuropean_metadata_url_alt, metadata_cache)
}

meta_coords_list <- list()
if (file.exists(metadata_cache)) {
  message("Reading EEA Pan-European metadata: ", metadata_cache)
  meta <- read_csv_flex(metadata_cache)
  code_cols <- names(meta)[tolower(names(meta)) %in% tolower(c(
    "AirQualityStationEoICode", "station_european_code", "AirQualityStation",
    "AirQualityStationNatCode", "SamplingPoint"
  ))]
  lon_col <- pick_col(meta, c("Longitude", "longitude", "lon", "Lon"))
  lat_col <- pick_col(meta, c("Latitude", "latitude", "lat", "Lat"))
  meta_coords_list[[length(meta_coords_list) + 1]] <- collapse_station_metadata(
    meta, code_cols = code_cols, lon_col = lon_col, lat_col = lat_col,
    source_label = "EEA PanEuropean metadata"
  )
} else {
  warning("No Pan-European metadata file available.")
}

# ---- Source 2: historical AirBase ArcGIS station service ----
if (isTRUE(try_arcgis_airbase_service)) {
  arcgis_cache <- file.path(ext_dir, "eea_airbase_v8_stations_arcgis.csv")
  arc <- fetch_arcgis_airbase_stations(airbase_station_query_url, arcgis_cache)
  if (!is.null(arc)) {
    message("Using AirBase ArcGIS station metadata.")
    code_cols <- names(arc)[tolower(names(arc)) %in% tolower(c(
      "station_european_code", "AirQualityStationEoICode", "EoICode",
      "AirQualityStation", "station_code", "station_local_code"
    ))]
    lon_col <- pick_col(arc, c("Longitude", "longitude", "lon", "Lon", "arcgis_longitude"))
    lat_col <- pick_col(arc, c("Latitude", "latitude", "lat", "Lat", "arcgis_latitude"))
    meta_coords_list[[length(meta_coords_list) + 1]] <- collapse_station_metadata(
      arc, code_cols = code_cols, lon_col = lon_col, lat_col = lat_col,
      source_label = "EEA AirBase v8 ArcGIS"
    )
  }
}

meta_coords_list <- Filter(Negate(is.null), meta_coords_list)
if (length(meta_coords_list) == 0) stop("No usable station-coordinate metadata could be obtained.")
meta_coords <- unique(do.call(rbind, meta_coords_list))

# Prefer non-ambiguous coordinate records and PanEuropean metadata if both exist.
meta_coords$source_rank <- ifelse(grepl("PanEuropean", meta_coords$coord_source), 1, 2)
meta_coords <- meta_coords[order(meta_coords$station_code_norm, meta_coords$coordinate_ambiguous, meta_coords$source_rank), ]
meta_coords <- meta_coords[!duplicated(meta_coords$station_code_norm), ]

# ---- Join ----
no2_aug <- augment_by_station_code(no2, meta_coords, station_col)

# Basic country-code consistency diagnostic, using first two characters of EoI code.
no2_aug$station_country_code <- substr(no2_aug$station_code_norm, 1, 2)
if (!is.na(country_col)) {
  no2_aug$country_name_original <- no2_aug[[country_col]]
}

# ---- Save augmented data ----
out_csv <- file.path(data_dir, "NO2_2010_prox_with_coords.csv")
out_rds <- file.path(data_dir, "NO2_2010_prox_with_coords.rds")
utils::write.csv(no2_aug, out_csv, row.names = FALSE)
saveRDS(no2_aug, out_rds)
message("Wrote: ", out_csv)
message("Wrote: ", out_rds)

# ---- Matching report ----
report_cols <- unique(c(
  station_col, "station_code_norm", country_col, "station_country_code",
  "longitude_recovered", "latitude_recovered", "coordinate_match",
  "coordinate_ambiguous", "coord_source", "key_source_column",
  "n_metadata_records", "n_unique_coord_pairs"
))
report_cols <- report_cols[!is.na(report_cols) & report_cols %in% names(no2_aug)]
report <- no2_aug[, report_cols, drop = FALSE]

report_path <- file.path(out_dir, "vizcaino_airbase_coordinate_match_report.csv")
utils::write.csv(report, report_path, row.names = FALSE)
message("Wrote: ", report_path)

# Summary by country, if available.
if (!is.na(country_col)) {
  country <- as.character(no2_aug[[country_col]])
} else {
  country <- no2_aug$station_country_code
}
summary <- aggregate(
  coordinate_match ~ country,
  data = data.frame(country = country, coordinate_match = no2_aug$coordinate_match),
  FUN = function(z) c(n = length(z), matched = sum(z), pct = 100 * mean(z))
)
summary <- data.frame(
  country = summary$country,
  n = summary$coordinate_match[, "n"],
  matched = summary$coordinate_match[, "matched"],
  pct_matched = round(summary$coordinate_match[, "pct"], 1),
  row.names = NULL
)
summary <- summary[order(-summary$pct_matched, summary$country), ]
summary_path <- file.path(out_dir, "vizcaino_airbase_coordinate_match_summary.csv")
utils::write.csv(summary, summary_path, row.names = FALSE)
message("Wrote: ", summary_path)

message("Overall coordinate matches: ", sum(no2_aug$coordinate_match), " / ", nrow(no2_aug),
        " (", round(100 * mean(no2_aug$coordinate_match), 1), "%)")

# ---- Optional: augment NO2_Europe.rds if it contains the same station-code key ----
rds_path <- file.path(base_dir, input_rds_rel)
if (file.exists(rds_path)) {
  message("Attempting to augment RDS: ", rds_path)
  d_rds <- readRDS(rds_path)
  station_col_rds <- pick_col(d_rds, c("station_eu", "station__1", "station_european_code", "AirQualityStationEoICode"))
  if (!is.na(station_col_rds)) {
    d_rds_aug <- augment_by_station_code(d_rds, meta_coords, station_col_rds)
    rds_out <- file.path(data_dir, "NO2_Europe_with_coordinates.rds")
    saveRDS(d_rds_aug, rds_out)
    message("Wrote: ", rds_out)
  } else {
    message("RDS has no station-code column; not augmented.")
  }
}
