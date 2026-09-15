# Generate supplementary post-hoc example overview table
#
# Numerical cells are populated from results/examples/posthoc_examples_combined.csv.
# The table is for the supplementary material and uses the same input terminology
# as the Tier 1 and Tier 2 Shiny apps.

VERBOSE <- 2

fmt <- function(x, digits = 3, bracket = FALSE) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("--")
  z <- formatC(x, digits = digits, format = "fg", flag = "#")
  if (bracket) paste0("[", z, "]") else z
}
fmt_fixed <- function(x, digits = 1, bracket = FALSE) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("--")
  z <- formatC(x, digits = digits, format = "f")
  if (bracket) paste0("[", z, "]") else z
}
fmt_int <- function(x, bracket = FALSE) fmt(x, digits = 0, bracket = bracket)
ni_or_fmt <- function(value, noninformative = FALSE, digits = 3) {
  if (isTRUE(noninformative)) return("n.i.")
  fmt(value, digits = digits)
}
tex <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "--"
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x
}
mc <- function(x) paste0("\\makecell[l]{", x, "}")
compact_prediction_support <- function(x) {
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return("--")

  # The supplementary overview table reports only the prediction/target support.
  # Reference-support descriptions are discussed in the text and in the example
  # scripts, but are intentionally suppressed here to keep the table compact.
  x <- strsplit(x, " / ", fixed = TRUE)[[1]][1]
  x <- gsub("approximately", "", x, ignore.case = TRUE)
  x <- gsub("about", "", x, ignore.case = TRUE)
  x <- gsub("grid cell", "", x, ignore.case = TRUE)
  x <- gsub("pixel", "", x, ignore.case = TRUE)
  x <- gsub("target support", "", x, ignore.case = TRUE)

  # Prefer the first explicit metric support, e.g. 30 m, 2 km, 100 m, 9 km.
  m <- regexpr("[0-9]+(?:\\.[0-9]+)?[[:space:]]*(?:km|m)", x, ignore.case = TRUE, perl = TRUE)
  if (m[1] > 0) {
    out <- regmatches(x, m)
    out <- gsub("[[:space:]]+", " ", out)
    return(trimws(out))
  }

  # Last-resort cleanup for unusual labels without an explicit metric unit.
  x <- gsub("^[,;:[:space:]]+|[,;:[:space:]]+$", "", x)
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

first <- function(d, names, default = NA) {
  for (nm in names) if (nm %in% names(d)) {
    v <- d[[nm]][which(!is.na(d[[nm]]) & d[[nm]] != "")]
    if (length(v)) return(v[1])
  }
  default
}
num <- function(d, names, default = NA_real_) suppressWarnings(as.numeric(first(d, names, default)))
logi <- function(d, names) {
  v <- first(d, names, FALSE)
  isTRUE(v) || identical(as.character(v), "TRUE") || identical(as.character(v), "1")
}

read_rows <- function(path = "results/examples/posthoc_examples_combined.csv") {
  if (!file.exists(path)) stop("Combined example result file not found: ", path)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

example_spec <- data.frame(
  col = c("Biomass", "BTS",
          "SWE 10 m q=0", "SWE 10 m q=.015",
          "SWE 100 m q=0", "SWE 100 m q=.015", "NO2--DE",
          "NO2--EU split", "NO2--EU full nugget", "Meuse--Zn", "Meuse--Cd",
          "S1--SM", "SMAP--CRNS"),
  pattern = c("Biomass", "BTS", "SWE", "SWE", "SWE", "SWE", "NO2-DE",
              "50/50", "full nugget", "Meuse--Zn", "Meuse--Cd", "S1--SM", "SMAP--CRNS"),
  support_short = c("30 m", "30 m", "10 m", "10 m", "100 m", "100 m", "2 km",
                    "100 m", "100 m", "100 m", "100 m", "1 km", "9 km"),
  stringsAsFactors = FALSE
)

pick_example <- function(dat, label) {
  st <- if ("study" %in% names(dat)) dat$study else dat$Study
  eid <- if ("example_id" %in% names(dat)) dat$example_id else rep("", nrow(dat))
  if (label == "SWE 10 m q=0") return(dat[grepl("SWE", st) & grepl("10m_eps0$", eid), , drop = FALSE])
  if (label == "SWE 10 m q=.015") return(dat[grepl("SWE", st) & grepl("10m_eps20$", eid), , drop = FALSE])
  if (label == "SWE 100 m q=0") return(dat[grepl("SWE", st) & grepl("100m_eps0$", eid), , drop = FALSE])
  if (label == "SWE 100 m q=.015") return(dat[grepl("SWE", st) & grepl("100m_eps20$", eid), , drop = FALSE])
  if (label == "NO2--EU split") return(dat[grepl("NO2-EU", st) & grepl("50/50", st), , drop = FALSE])
  if (label == "NO2--EU full nugget") return(dat[grepl("NO2-EU", st) & grepl("full nugget", st), , drop = FALSE])
  dat[st == label | gsub("-", "--", st, fixed = TRUE) == label, , drop = FALSE]
}

cell_for <- function(d, field) {
  if (!nrow(d)) return("--")
  d <- d[1, , drop = FALSE]
  br <- field %in% c("sigma2_eps", "sigma2_mu", "true_nugget", "partial_sill", "range", "lambda", "delta") && grepl("NO2-EU|SWE|BTS|Biomass|SMAP", first(d, "study", ""))
  switch(field,
    response = mc(gsub("Annual mean NO2", "annual mean\\\\NO$_2$", tex(first(d, "variable", "--")), fixed = TRUE)),
    support = mc(tex(compact_prediction_support(first(d, "prediction_support", "--")))),
    tier = mc(gsub("Tier ", "", tex(first(d, "information_tier", "--")), fixed = TRUE)),
    metric = {
      rdef <- first(d, "r2_definition", "--")
      mc(if (rdef == "squared_correlation") "$R^2_{\\mathrm{cor}}$" else if (rdef == "predictive_skill") "$R^2_{\\mathrm{pred}}$" else "--")
    },
    raw_r2 = mc(fmt(num(d, "r2_obs"), 3, bracket = grepl("Biomass|NO2-EU", first(d,"study","")))),
    raw_rmse = {
      value <- num(d, "rmse_obs")
      study <- first(d, "study", "")
      if (identical(study, "NO2-DE")) {
        mc(fmt_fixed(value, 1))
      } else {
        mc(fmt(value, ifelse(grepl("Meuse|S1|SMAP", study), 3, 2),
               bracket = grepl("Biomass|NO2-EU", study)))
      }
    },
    sigma_obs = mc(fmt(sqrt(num(d, "sigma2_obs")), 3, bracket = grepl("Biomass|NO2-EU|BTS", first(d,"study","")))),
    sigma2_obs = mc(fmt(num(d, "sigma2_obs"), 3, bracket = grepl("Biomass|NO2-EU|BTS", first(d,"study","")))),
    covariance_model = {
      model <- first(d, c("structural_covariance_model", "covariance_model", "variogram_model"), "--")
      model <- if (grepl("exp", model, ignore.case = TRUE)) {
        "Exponential"
      } else if (grepl("sph", model, ignore.case = TRUE)) {
        "Spherical"
      } else {
        model
      }
      mc(tex(model))
    },
    range = mc(fmt(num(d, c("effective_range", "variogram_range", "spherical_range", "variogram_spherical_range_m", "variogram_range_m")), 3, bracket = br)),
    sigma2_nugget_true = mc(fmt(num(d, c("residual_nugget_true", "true_nugget", "true_nugget_used_for_support", "variogram_nugget")), 3, bracket = br)),
    sigma2_eps = mc(fmt(num(d, "sigma2_eps"), 3, bracket = grepl("Biomass|NO2-EU|BTS|SWE|SMAP", first(d,"study","")))),
    sigma2_S = mc(fmt(num(d, c("partial_sill", "partial_sill_scaled", "variogram_spherical_psill")), 3, bracket = br)),
    q_eps = mc(fmt(num(d, "q_eps"), 3, bracket = grepl("NO2-EU|SWE", first(d,"study","")))),
    lambda = mc(fmt(num(d, "lambda"), 3, bracket = br)),
    delta = {
      value <- num(d, "sigma2_delta")
      if (identical(first(d, "study", ""), "NO2-DE")) {
        mc(fmt_fixed(value, 1))
      } else {
        mc(fmt(value, 3, bracket = br))
      }
    },
    tier1_r2 = mc(fmt(num(d, "R2_ME"), 3, bracket = grepl("NO2-EU|Biomass", first(d,"study","")))),
    tier1_rmse_ratio = {
      rr <- num(d, "RMSE_ME") / num(d, "rmse_obs")
      mc(fmt(rr, 2))
    },
    tier23_r2 = mc(ni_or_fmt(num(d, "R2_B"), logi(d, c("r2_noninformative", "support_variance_exceeds_observed_mse")), 3)),
    tier23_rmse_ratio = {
      rr <- num(d, "RMSE_B") / num(d, "rmse_obs")
      mc(ni_or_fmt(rr, logi(d, "mse_noninformative"), 2))
    },
    "--")
}

make_table <- function(dat) {
  labels <- example_spec$col
  n_table_cols <- length(labels) + 1
  section_row <- function(title, add_space = FALSE) {
    prefix <- if (add_space) "\\addlinespace[0.35em]" else ""
    paste0(prefix, "\\multicolumn{", n_table_cols,
           "}{@{}l}{\\textbf{", title, "}}\\\\")
  }
  rows <- list(
    c(section_row("Study and support setting")),
    c("Response", vapply(labels, function(l) cell_for(pick_example(dat,l), "response"), character(1))),
    c("Prediction support", vapply(labels, function(l) mc(tex(example_spec$support_short[example_spec$col == l])), character(1))),
    c("Tier", vapply(labels, function(l) cell_for(pick_example(dat,l), "tier"), character(1))),
    c(section_row("Validation inputs", add_space = TRUE)),
    c("Metric", vapply(labels, function(l) cell_for(pick_example(dat,l), "metric"), character(1))),
    c("Raw $R^2$", vapply(labels, function(l) cell_for(pick_example(dat,l), "raw_r2"), character(1))),
    c("Raw RMSE/ubRMSE", vapply(labels, function(l) cell_for(pick_example(dat,l), "raw_rmse"), character(1))),
    c("$\\sigma_{\\mathrm{obs}}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "sigma_obs"), character(1))),
    c("$\\sigma^2_{\\mathrm{obs}}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "sigma2_obs"), character(1))),
    c(section_row("Support inputs", add_space = TRUE)),
    c("Structural covariance model", vapply(labels, function(l) cell_for(pick_example(dat,l), "covariance_model"), character(1))),
    c("Range (m)", vapply(labels, function(l) cell_for(pick_example(dat,l), "range"), character(1))),
    c("True microscale variance $\\sigma^2_{\\mu}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "sigma2_nugget_true"), character(1))),
    c("Measurement-error variance $\\sigma^2_{\\varepsilon}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "sigma2_eps"), character(1))),
    c("Structural partial sill $\\sigma^2_S$", vapply(labels, function(l) cell_for(pick_example(dat,l), "sigma2_S"), character(1))),
    c("$q_{\\varepsilon}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "q_eps"), character(1))),
    c("$\\lambda_{A,B}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "lambda"), character(1))),
    c("$\\Delta_{A,B}$", vapply(labels, function(l) cell_for(pick_example(dat,l), "delta"), character(1))),
    c(section_row("Adjusted metrics", add_space = TRUE)),
    c("Tier 1 $R^2$", vapply(labels, function(l) cell_for(pick_example(dat,l), "tier1_r2"), character(1))),
    c("Tier 1 RMSE/raw", vapply(labels, function(l) cell_for(pick_example(dat,l), "tier1_rmse_ratio"), character(1))),
    c("Tier 2/3 $R^2$", vapply(labels, function(l) cell_for(pick_example(dat,l), "tier23_r2"), character(1))),
    c("Tier 2/3 RMSE/raw", vapply(labels, function(l) cell_for(pick_example(dat,l), "tier23_rmse_ratio"), character(1)))
  )
  lines <- c(
    "% Automatically generated by R/paper/make_suppl_posthoc_example_table.R",
    "% Requires: booktabs, adjustbox, array, makecell.",
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Numerical summary of the post-hoc examples. Both SWE measurement-error settings are reported for each support: $q_\\varepsilon=0$ and $q_\\varepsilon\\approx0.015$ ($\\sigma_\\varepsilon=20$~mm). Spherical ranges are practical ranges; the SWE exponential range is the reported effective range. Square brackets mark values reconstructed from reported summaries, transferred from related data, or specified as sensitivity assumptions. Dashes indicate quantities that were unavailable or not used. The entry n.i. denotes a non-informative adjusted endpoint, for example when the implied support adjustment would lie outside the feasible metric range.}",
    "\\label{tab:posthoc-examples-appendix}",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2.0pt}",
    "\\renewcommand{\\arraystretch}{1.07}",
    "\\begin{adjustbox}{width=\\linewidth,totalheight=0.90\\textheight,keepaspectratio}",
    paste0("\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{0.130\\linewidth}*{",
           length(labels), "}{>{\\raggedright\\arraybackslash}p{0.066\\linewidth}}@{}}"),
    "\\toprule",
    paste(c("Quantity", labels), collapse = " & "), "\\\\",
    "\\midrule"
  )
  for (r in rows) {
    if (length(r) == 1) lines <- c(lines, r)
    else lines <- c(lines, paste(r, collapse = " & "), "\\\\")
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{adjustbox}", "\\end{table}")
}

make_suppl_posthoc_example_table <- function(combined_csv = "results/examples/posthoc_examples_combined.csv",
                                             out_tex = "results/supplementary/suppl_posthoc_examples_overview_table.tex",
                                             out_csv = "results/supplementary/suppl_posthoc_examples_overview_table.csv") {
  dat <- read_rows(combined_csv)
  labels <- example_spec$col
  # also save CSV in wide display form for checking
  row_names <- c("Response", "Prediction support", "Tier", "Metric", "Raw R2", "Raw RMSE/ubRMSE",
                 "sigma_obs", "sigma2_obs", "Structural covariance model", "Range_m", "sigma2_mu", "sigma2_epsilon",
                 "sigma2_S", "q_epsilon", "lambda_A_B", "Delta_A_B",
                 "Tier 1 R2", "Tier 1 RMSE/raw", "Tier 2/3 R2", "Tier 2/3 RMSE/raw")
  fields <- c("response","support","tier","metric","raw_r2","raw_rmse","sigma_obs","sigma2_obs",
              "covariance_model","range","sigma2_nugget_true","sigma2_eps","sigma2_S","q_eps","lambda","delta",
              "tier1_r2","tier1_rmse_ratio","tier23_r2","tier23_rmse_ratio")
  wide <- data.frame(Quantity = row_names, stringsAsFactors = FALSE)
  for (l in labels) {
    wide[[l]] <- vapply(fields, function(f) {
      if (identical(f, "support")) return(example_spec$support_short[example_spec$col == l])
      gsub("\\\\makecell\\[l\\]\\{|\\}$", "", cell_for(pick_example(dat,l), f))
    }, character(1))
  }
  dir.create(dirname(out_tex), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(wide, out_csv, row.names = FALSE)
  writeLines(make_table(dat), out_tex)
  if (VERBOSE >= 1) message("Wrote: ", out_tex)
  invisible(wide)
}

if (file.exists("results/examples/posthoc_examples_combined.csv")) {
  make_suppl_posthoc_example_table()
}
