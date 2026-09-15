# Shared helpers for compact post-hoc examples
#
# These helpers keep the published examples aligned with the manuscript
# terminology: raw reference-support metrics, Tier 1 reference-error adjustment,
# and Tier 2 support adjustment where available.

source("R/core/support_correction_core.R")

#' Safe numeric rounding for manuscript tables
fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), "", formatC(x, digits = digits, format = "f"))
}

#' Save one example row
save_example_row <- function(x,
                             out_dir = "results/examples",
                             overwrite = TRUE,
                             verbose = 1) {
  stopifnot(is.data.frame(x), nrow(x) >= 1)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(out_dir, paste0(x$example_id[1], ".csv"))
  write_csv_incremental(x, f, overwrite = overwrite, verbose = verbose)
  invisible(x)
}

#' Add manuscript-table fields to a core result row
#'
#' @param x Result row produced by make_support_result().
#' @param study Short study label.
#' @param prediction_support Prediction support label.
#' @param reference_support Reference support label.
#' @param information_tier Tier label.
#' @param available_information Short description of used information.
#' @param main_limitation Short limitation note.
#' @return data.frame with additional table fields.
add_example_metadata <- function(x,
                                 study,
                                 prediction_support,
                                 reference_support,
                                 information_tier,
                                 available_information,
                                 main_limitation = "") {
  x$study <- study
  x$prediction_support <- prediction_support
  x$reference_support <- reference_support
  x$information_tier <- information_tier
  x$available_information <- available_information
  x$main_limitation <- main_limitation
  x
}

#' Combine all example result CSVs and write a manuscript table
#' Row-bind data frames with different columns
#'
#' Some examples add audit fields that are specific to one study.  Base
#' rbind() fails when these columns differ.  This helper first creates the
#' union of all column names and fills absent columns with NA.
rbind_fill_base <- function(x) {
  stopifnot(is.list(x), length(x) > 0)
  all_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  x2 <- lapply(x, function(d) {
    missing <- setdiff(all_names, names(d))
    if (length(missing) > 0) {
      for (m in missing) d[[m]] <- NA
    }
    d[all_names]
  })
  do.call(rbind, x2)
}

#' Return a column or an NA vector of the required length
col_or_na <- function(d, nm) {
  if (nm %in% names(d)) d[[nm]] else rep(NA, nrow(d))
}

combine_posthoc_examples <- function(out_dir = "results/examples",
                                     overwrite = TRUE,
                                     verbose = 1) {
  files <- list.files(out_dir, pattern = "^example_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No example_*.csv files found in ", out_dir)

  if (verbose >= 2) {
    message("Combining ", length(files), " example result file(s):")
    message(paste("  -", basename(files)), sep = "\n")
  }

  rows <- rbind_fill_base(lapply(files, utils::read.csv, stringsAsFactors = FALSE))

  # Harmonized table columns for manuscript-facing summaries.
  tab <- data.frame(
    Study = col_or_na(rows, "study"),
    Variable = col_or_na(rows, "variable"),
    Tier = col_or_na(rows, "information_tier"),
    `Prediction support` = col_or_na(rows, "prediction_support"),
    `Reference support` = col_or_na(rows, "reference_support"),
    `R2 definition` = col_or_na(rows, "r2_definition"),
    `q_epsilon` = fmt_num(col_or_na(rows, "q_eps"), 3),
    `lambda` = fmt_num(col_or_na(rows, "lambda"), 3),
    `R2_raw` = fmt_num(col_or_na(rows, "r2_obs"), 3),
    `R2_Tier1` = fmt_num(col_or_na(rows, "R2_ME"), 3),
    `R2_Tier2` = fmt_num(col_or_na(rows, "R2_B"), 3),
    `RMSE_raw` = fmt_num(col_or_na(rows, "rmse_obs"), 3),
    `RMSE_Tier1` = fmt_num(col_or_na(rows, "RMSE_ME"), 3),
    `RMSE_Tier2` = fmt_num(col_or_na(rows, "RMSE_B"), 3),
    `R2_Tier2_raw` = fmt_num(col_or_na(rows, "R2_B_raw"), 3),
    `RMSE_Tier2_raw` = fmt_num(col_or_na(rows, "RMSE_B_raw"), 3),
    `R2_Tier2_noninformative` = as.character(as.logical(col_or_na(rows, "r2_noninformative"))),
    `RMSE_Tier2_noninformative` = as.character(as.logical(col_or_na(rows, "mse_noninformative"))),
    Limitation = col_or_na(rows, "main_limitation"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Keep rows stable and manuscript-readable: first by study, then by tier and
  # support setting.  Missing fields are pushed to the end by order().
  ord <- order(tab$Study, tab$Tier, tab$`Prediction support`, tab$q_epsilon, na.last = TRUE)
  rows <- rows[ord, , drop = FALSE]
  tab <- tab[ord, , drop = FALSE]

  combined_file <- file.path(out_dir, "posthoc_examples_combined.csv")
  table_file <- file.path(out_dir, "posthoc_examples_table.csv")
  tex_file <- file.path(out_dir, "posthoc_examples_table.tex")
  write_csv_incremental(rows, combined_file, overwrite = overwrite, verbose = verbose)
  write_csv_incremental(tab, table_file, overwrite = overwrite, verbose = verbose)

  writeLines(make_posthoc_examples_latex(tab), tex_file)
  if (verbose >= 1) message("Wrote: ", tex_file)

  invisible(list(rows = rows, table = tab))
}

#' Create LaTeX table for the post-hoc examples
make_posthoc_examples_latex <- function(tab) {
  esc <- function(x) {
    x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
    x <- gsub("_", "\\_", x, fixed = TRUE)
    x <- gsub("%", "\\%", x, fixed = TRUE)
    x
  }

  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Post-hoc examples illustrating how published or standard validation metrics can be reinterpreted under the proposed information tiers. Tier~1 uses measurement- or reference-error information only; Tier~2 additionally uses covariance or representativeness information. Empty cells indicate quantities that are not identifiable from the available summary information.}",
    "\\label{tab:posthoc-examples}",
    "\\small",
    "\\begin{tabular}{llllrrrrrrl}",
    "\\hline",
    "Study & Variable & Tier & Support & $q_\\varepsilon$ & $\\lambda$ & $R^2_{raw}$ & $R^2_{T1}$ & $R^2_{T2}$ & RMSE$_{raw}$ & RMSE$_{T1/T2}$ \\\\",
    "\\hline"
  )

  for (i in seq_len(nrow(tab))) {
    support <- paste0(tab$`Prediction support`[i], " / ", tab$`Reference support`[i])
    rmse_t <- paste0(
      tab$RMSE_Tier1[i],
      ifelse(tab$RMSE_Tier2[i] != "", paste0(" / ", tab$RMSE_Tier2[i]), "")
    )
    row <- paste(
      esc(tab$Study[i]),
      esc(tab$Variable[i]),
      esc(tab$Tier[i]),
      esc(support),
      tab$q_epsilon[i],
      tab$lambda[i],
      tab$R2_raw[i],
      tab$R2_Tier1[i],
      tab$R2_Tier2[i],
      tab$RMSE_raw[i],
      rmse_t,
      sep = " & "
    )
    lines <- c(lines, paste0(row, " \\\\"))
  }

  c(
    lines,
    "\\hline",
    "\\multicolumn{11}{p{0.95\\textwidth}}{\\footnotesize Note: The Meuse example treats bulk soil-sample observations as point-support observations for illustration, although they have finite sampling support.} \\\\",
    "\\end{tabular}",
    "\\end{table}"
  )
}

#' Closed-form equal-area disk approximation used in compact examples
#'
#' This is a transparent approximation for cases where only a short-lag
#' semivariogram value at the equivalent radius is available.
#' @return Named vector with lambda and sigma2_delta.
equal_area_disk_linear_terms <- function(sigma2_point_true,
                                         sigma2_eps,
                                         gamma_true_R,
                                         true_nugget = 0) {
  stopifnot(is.finite(sigma2_point_true), sigma2_point_true > 0)
  stopifnot(is.finite(sigma2_eps), sigma2_eps >= 0)
  stopifnot(is.finite(gamma_true_R), gamma_true_R >= 0)
  stopifnot(is.finite(true_nugget), true_nugget >= 0)

  # Closed-form equal-area disk approximation from the revised manuscript.
  # The linear short-lag approximation is applied only to the spatially
  # correlated structural covariance component.  A true microscale/nugget
  # component contributes to the point/reference-support variance and to
  # Delta_{p,B}, but averages out of an ideal continuous block mean.
  if (true_nugget > sigma2_point_true) {
    warning("true_nugget exceeds sigma2_point_true; truncating to sigma2_point_true.")
    true_nugget <- sigma2_point_true
  }

  sigma2_struct <- sigma2_point_true - true_nugget
  sigma2_obs <- sigma2_point_true + sigma2_eps
  q_eps <- sigma2_eps / sigma2_obs
  k <- 128 / (45 * pi)

  # gamma_true_R is allowed to be supplied as the full latent semivariance at R,
  # including true microscale variation.  For the linear approximation, remove
  # the microscale component and use only the structural short-lag semivariance.
  gamma_true_R <- min(gamma_true_R, sigma2_point_true)
  gamma_struct_R <- max(gamma_true_R - true_nugget, 0)
  gamma_struct_R <- min(gamma_struct_R, sigma2_struct)

  rho_struct_R <- if (sigma2_struct > 0) {
    max(min(1 - gamma_struct_R / sigma2_struct, 1), -1)
  } else {
    NA_real_
  }
  rho_R <- max(min(1 - gamma_true_R / sigma2_point_true, 1), -1)

  sigma2_delta_struct <- k * gamma_struct_R
  sigma2_block <- max(sigma2_struct - sigma2_delta_struct, 0)
  sigma2_delta <- true_nugget + sigma2_delta_struct
  lambda <- sigma2_block / sigma2_obs

  c(
    rho_R = rho_R,
    rho_struct_R = rho_struct_R,
    q_eps = q_eps,
    lambda = lambda,
    sigma2_delta = sigma2_delta,
    sigma2_delta_struct = sigma2_delta_struct,
    sigma2_block = sigma2_block,
    cov_point_block = NA_real_,
    true_nugget = true_nugget,
    partial_sill = sigma2_struct,
    gamma_struct_R = gamma_struct_R
  )
}
