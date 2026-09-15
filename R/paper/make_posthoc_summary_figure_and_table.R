# Generate manuscript-facing post-hoc example table and summary figure
#
# This script replaces the earlier wide quantitative table by:
#   1. a compact overview table of settings, sources and information tiers;
#   2. a quantitative figure comparing raw reference-support metrics with
#      Tier 1 and Tier 2 adjustments.
#
# Input:
#   results/examples/posthoc_examples_combined.csv
#
# Output:
#   results/examples/posthoc_examples_overview_table.tex
#   results/examples/fig_posthoc_examples_summary.pdf
#   results/examples/fig_posthoc_examples_summary.png
#   results/examples/posthoc_examples_plot_data.csv
#
# The output paths are intentionally compatible with manuscript inclusion via:
#   \input{../analysis/results/examples/posthoc_examples_overview_table.tex}
#   \includegraphics[width=\textwidth]{../code_data/results/examples/fig_posthoc_examples_summary.pdf}

VERBOSE <- 2
OVERWRITE <- TRUE

# -----------------------------
# Small utilities
# -----------------------------

msg <- function(..., level = 1) {
  if (isTRUE(VERBOSE >= level)) message(...)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

latex_escape <- function(x) {
  # Escape text fields for LaTeX.  This character-by-character approach
  # avoids replacement-string backslash doubling in gsub().
  if (length(x) == 0) return(x)
  vapply(as.character(x), function(one) {
    if (is.na(one)) return(NA_character_)
    chars <- strsplit(one, "", fixed = TRUE)[[1]]
    mapped <- vapply(chars, function(ch) {
      switch(ch,
             "\\" = "\\textbackslash{}",
             "&"  = "\\&",
             "%"  = "\\%",
             "$"  = "\\$",
             "#"  = "\\#",
             "_"  = "\\_",
             "{"  = "\\{",
             "}"  = "\\}",
             "~"  = "\\textasciitilde{}",
             "^"  = "\\textasciicircum{}",
             ch)
    }, character(1))
    paste0(mapped, collapse = "")
  }, character(1), USE.NAMES = FALSE)
}
fmt_range <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) return("--")
  if (length(unique(round(x, digits))) == 1) return(formatC(x[1], format = "f", digits = digits))
  paste0(formatC(min(x), format = "f", digits = digits), "--", formatC(max(x), format = "f", digits = digits))
}

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & as.character(x) != ""]
  if (length(x) == 0) return(NA_character_)
  as.character(x[1])
}

source_label_latex <- function(study) {
  # Keep citation commands unescaped.  Add or adjust keys here if the BibTeX
  # database uses different names.
  if (grepl("^BTS$|Brenning", study)) return("\\citet{brenning2005bts}")
  if (grepl("^SWE|Plattner", study)) return("\\citet{plattner2006swe}")
  if (grepl("^Biomass$|Zandler", study)) return("\\citet{zandler2015dwarf}")
  if (grepl("SMAP-CRNS|Zhao", study)) return("\\citet{zhao2021subsurface}")
  if (grepl("S1-SM|Balenzano", study)) return("\\citet{balenzano2021sentinel1}")
  if (grepl("NO2-EU|Vizcaino|Lavalle", study)) return("\\citet{vizcaino2018no2}")
  if (grepl("Meuse", study)) return("Meuse data$^{a}$")
  if (grepl("NO2-DE|NO2|NO_2|NO\\$_2", study)) return("NO$_2$ Germany 2018")
  latex_escape(study)
}

make_example_group <- function(x) {
  # Collapses sensitivity rows into manuscript-facing examples.
  x <- as.character(x)
  x <- gsub(", 10 m$", "", x)
  x <- gsub(", 100 m$", "", x)
  x
}

support_size_m <- function(pred_support) {
  # Extract first number in metres from strings such as "100 m grid cell" or
  # "SMAP 9 km pixel". Returns NA for qualitative support descriptions.
  s <- tolower(as.character(pred_support))
  out <- rep(NA_real_, length(s))
  has_km <- grepl("km", s)
  num <- suppressWarnings(as.numeric(sub(".*?([0-9]+(?:\\.[0-9]+)?).*", "\\1", s)))
  has_num <- !is.na(num) & grepl("[0-9]", s)
  out[has_num] <- num[has_num]
  out[has_num & has_km] <- out[has_num & has_km] * 1000
  out
}

support_class <- function(size_m) {
  out <- rep("not reported", length(size_m))
  out[is.finite(size_m) & size_m <= 30] <- "small support"
  out[is.finite(size_m) & size_m > 30 & size_m <= 250] <- "medium support"
  out[is.finite(size_m) & size_m > 250] <- "large support"
  factor(out, levels = c("not reported", "small support", "medium support", "large support"))
}

# -----------------------------
# Overview table
# -----------------------------

make_posthoc_overview_table <- function(combined_csv = "results/examples/posthoc_examples_combined.csv",
                                        out_tex = "results/examples/posthoc_examples_overview_table.tex",
                                        overwrite = OVERWRITE,
                                        verbose = VERBOSE) {
  if (file.exists(out_tex) && !overwrite) {
    msg("Keeping existing table: ", out_tex, level = 1)
    return(invisible(out_tex))
  }
  dat <- utils::read.csv(combined_csv, stringsAsFactors = FALSE)
  dat$example_group <- make_example_group(dat$study)

  groups <- split(dat, dat$example_group)
  tab <- lapply(groups, function(d) {
    study <- first_nonmissing(d$example_group)
    qtxt <- fmt_range(d$q_eps, 3)
    ltxt <- fmt_range(d$lambda, 3)
    r2def <- first_nonmissing(d$r2_definition)
    avail <- first_nonmissing(d$available_information)
    if (is.na(avail)) avail <- "Published validation/error summaries"

    info <- paste0(latex_escape(avail),
                   "; $q_{\\varepsilon}$ = ", latex_escape(qtxt),
                   ", $\\lambda$ = ", latex_escape(ltxt),
                   "; ", latex_escape(r2def))

    data.frame(
      Source = source_label_latex(study),
      Variable = latex_escape(first_nonmissing(d$variable)),
      Support = latex_escape(paste0(first_nonmissing(d$prediction_support), " / ",
                                    first_nonmissing(d$reference_support))),
      Information = info,
      Tier = latex_escape(first_nonmissing(d$information_tier)),
      Limitation = latex_escape(first_nonmissing(d$main_limitation)),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, tab)

  body_rows <- paste0(tab$Source, " & ", tab$Variable, " & ", tab$Support, " & ",
                      tab$Information, " & ", tab$Tier, " & ", tab$Limitation, " \\\\")

  lines <- c(
    "% Automatically generated by R/30_make_posthoc_table_and_figure.R",
    "% Generated as a standard tabular environment; no tabularx dependency.",
    "\\begin{table}[t]",
    "\\centering",
    "\\footnotesize",
    "\\caption{Post-hoc examples used to illustrate the information tiers. The table records the source, support setting, definition of reported or reconstructed $R^2$, and information available for adjustment; quantitative changes in the metrics are displayed in Fig.~\\ref{fig:posthoc-examples-summary}. Empty or unavailable quantities in the underlying results indicate that the corresponding tier is not identifiable from the published information.}",
    "\\label{tab:posthoc-examples-overview}",
    "\\begin{tabular}{p{0.14\\textwidth} p{0.10\\textwidth} p{0.18\\textwidth} p{0.29\\textwidth} p{0.07\\textwidth} p{0.15\\textwidth}}",
    "\\hline",
    "Source & Variable & Prediction/reference support & Information used & Tier & Main limitation \\\\",
    "\\hline",
    body_rows,
    "\\hline",
    "\\multicolumn{6}{p{0.96\\textwidth}}{\\footnotesize $^{a}$For predictive/test-set $R^2$, adjustments are performed on the MSE scale and then transformed back to $R^2$ using the target-support variance. For the SWE sensitivity rows, the ordinary-regression $R^2$ is treated as a squared-correlation proxy; the nugget-like term is split heuristically into measurement error and true microscale variability. The S1-SM reported-SRE row uses the spatial representativeness error adjustment as reported in the original Sentinel-1 soil-moisture validation study. The SMAP--CRNS row uses the TC-estimated CRNS reference-error standard deviation for a partial Tier~1 ubRMSD adjustment only. For the NO2--EU rows, station coordinates recovered from AirBase metadata are used to fit a European semivariogram, with sensitivity rows that decompose the fitted nugget into measurement error and true microscale variation. The Meuse log-zinc and log-cadmium examples use LOOCV predictive $R^2$ and treat bulk soil-sample observations as point-support observations for illustration, although they have finite sampling support. If the inferred support-error variance exceeds the observed prediction-error variance, the Tier~2 RMSE adjustment is flagged as non-informative rather than truncated to zero.} \\\\",
    "\\end{tabular}",
    "\\end{table}"
  )
  ensure_dir(dirname(out_tex))
  writeLines(lines, out_tex, useBytes = TRUE)
  msg("Wrote overview table: ", out_tex, level = 1)
  invisible(out_tex)
}

# -----------------------------
# Plotting data and figure
# -----------------------------

make_posthoc_plot_data <- function(table_csv = "results/examples/posthoc_examples_table.csv",
                                   out_csv = "results/examples/posthoc_examples_plot_data.csv",
                                   overwrite = OVERWRITE,
                                   verbose = VERBOSE) {
  if (file.exists(out_csv) && !overwrite) {
    msg("Keeping existing plot data: ", out_csv, level = 1)
    return(utils::read.csv(out_csv, stringsAsFactors = FALSE))
  }

  dat <- utils::read.csv(table_csv, stringsAsFactors = FALSE)
  dat$row_id <- seq_len(nrow(dat))
  dat$study_group <- make_example_group(dat$Study)

  # Handle both check.names=TRUE and check.names=FALSE column names.
  get_col <- function(d, candidates) {
    nm <- names(d)
    hit <- candidates[candidates %in% nm]
    if (length(hit) == 0) stop("None of these columns found: ", paste(candidates, collapse = ", "))
    d[[hit[1]]]
  }
  get_col_or <- function(d, candidates, default) {
    nm <- names(d)
    hit <- candidates[candidates %in% nm]
    if (length(hit) == 0) return(rep(default, nrow(d)))
    d[[hit[1]]]
  }

  pred_support <- get_col(dat, c("Prediction.support", "Prediction support"))
  ref_support <- get_col(dat, c("Reference.support", "Reference support"))
  dat$support_m <- support_size_m(pred_support)
  dat$support_class <- support_class(dat$support_m)

  # Construct compact labels.  Keep SWE sensitivity rows separate because
  # they represent different supports and measurement-error scenarios.
  q <- suppressWarnings(as.numeric(get_col(dat, c("q_epsilon"))))
  dat$example_label <- dat$Study
  is_swe <- grepl("^SWE", dat$Study) | grepl("Plattner", dat$Study)
  dat$example_label[is_swe] <- paste0(dat$Study[is_swe], ", q=", formatC(q[is_swe], format = "f", digits = 3))

  # Keep a plain-character copy for robust ordering.  R factors can otherwise
  # carry stale levels through rbind()/split() operations, which can make the
  # plotted y-axis appear to follow file/read order rather than alphabetical
  # display-label order.
  dat$example_label_text <- as.character(dat$example_label)

  # Top-to-bottom alphabetical order in a horizontal ggplot requires the scale
  # limits to be the reverse of the alphabetic order, because the first discrete
  # y level is plotted at the bottom.
  example_levels_top <- sort(unique(dat$example_label_text), na.last = NA)
  example_levels_bottom <- rev(example_levels_top)
  dat$example_label <- factor(dat$example_label_text, levels = example_levels_bottom)

  r2_raw <- suppressWarnings(as.numeric(get_col(dat, c("R2_raw", "R2_naive"))))
  r2_t1 <- suppressWarnings(as.numeric(get_col(dat, c("R2_Tier1"))))
  r2_t2 <- suppressWarnings(as.numeric(get_col(dat, c("R2_Tier2"))))
  r2_t2_raw <- suppressWarnings(as.numeric(get_col_or(dat, c("R2_Tier2_raw"), NA_real_)))
  rmse_raw <- suppressWarnings(as.numeric(get_col(dat, c("RMSE_raw", "RMSE_naive"))))
  rmse_t1 <- suppressWarnings(as.numeric(get_col(dat, c("RMSE_Tier1"))))
  rmse_t2 <- suppressWarnings(as.numeric(get_col(dat, c("RMSE_Tier2"))))
  rmse_t2_raw <- suppressWarnings(as.numeric(get_col_or(dat, c("RMSE_Tier2_raw"), NA_real_)))
  r2_t2_noninformative <- suppressWarnings(as.logical(get_col_or(dat, c("R2_Tier2_noninformative", "R2_Tier2_truncated"), FALSE)))
  rmse_t2_noninformative <- suppressWarnings(as.logical(get_col_or(dat, c("RMSE_Tier2_noninformative", "RMSE_Tier2_truncated"), FALSE)))
  r2_t2_noninformative[is.na(r2_t2_noninformative)] <- FALSE
  rmse_t2_noninformative[is.na(rmse_t2_noninformative)] <- FALSE

  # Non-informative Tier 2 values are not interpreted as estimates, but the
  # raw diagnostic endpoint is useful for drawing an arrow to the panel edge.
  r2_t2_plot <- r2_t2
  r2_t2_plot[!is.finite(r2_t2_plot) & r2_t2_noninformative & is.finite(r2_t2_raw)] <-
    r2_t2_raw[!is.finite(r2_t2_plot) & r2_t2_noninformative & is.finite(r2_t2_raw)]
  r2_t2_plot[!is.finite(r2_t2_plot) & r2_t2_noninformative] <- 1

  rmse_t2_plot <- rmse_t2
  rmse_t2_plot[!is.finite(rmse_t2_plot) & rmse_t2_noninformative & is.finite(rmse_t2_raw)] <-
    rmse_t2_raw[!is.finite(rmse_t2_plot) & rmse_t2_noninformative & is.finite(rmse_t2_raw)]
  rmse_t2_plot[!is.finite(rmse_t2_plot) & rmse_t2_noninformative] <- 0

  make_long <- function(value, tier, panel, bound_noninformative = rep(FALSE, nrow(dat))) {
    data.frame(
      row_id = dat$row_id,
      Study = dat$Study,
      Variable = dat$Variable,
      Tier = tier,
      panel = panel,
      value = value,
      example_label = dat$example_label,
      example_label_text = dat$example_label_text,
      support_class = dat$support_class,
      support_m = dat$support_m,
      prediction_support = pred_support,
      reference_support = ref_support,
      bound_noninformative = as.logical(bound_noninformative),
      stringsAsFactors = FALSE
    )
  }

  r2_long <- rbind(
    make_long(r2_raw, "Raw reference support", "R2"),
    make_long(r2_t1, "Tier 1", "R2"),
    make_long(r2_t2_plot, "Tier 2", "R2", bound_noninformative = r2_t2_noninformative)
  )

  err_long <- rbind(
    make_long(rep(1, nrow(dat)), "Raw reference support", "Relative error"),
    make_long(rmse_t1 / rmse_raw, "Tier 1", "Relative error"),
    make_long(rmse_t2_plot / rmse_raw, "Tier 2", "Relative error", bound_noninformative = rmse_t2_noninformative)
  )

  plot_dat <- rbind(r2_long, err_long)
  plot_dat <- plot_dat[is.finite(plot_dat$value), ]
  plot_dat$Tier <- factor(plot_dat$Tier,
                          levels = c("Raw reference support", "Tier 1", "Tier 2"))
  # Tier 2 values that are outside the feasible metric range are diagnostic /
  # non-informative under the updated manuscript wording.  Keep their raw values
  # in the CSV, but omit the point itself from the figure scale and draw an arrow
  # to the relevant panel edge.
  plot_dat$bound_noninformative <- as.logical(plot_dat$bound_noninformative)
  plot_dat$bound_noninformative[is.na(plot_dat$bound_noninformative)] <- FALSE
  plot_dat$noninformative_endpoint <- plot_dat$Tier == "Tier 2" & plot_dat$bound_noninformative

  plot_dat$point_value <- plot_dat$value
  plot_dat$point_value[plot_dat$noninformative_endpoint] <- NA_real_

  # Axis limits are calculated from all finite displayed metric values;
  # non-informative Tier 2 endpoints are excluded because they are diagnostics,
  # not finite adjusted estimates.
  plot_dat$axis_range_value <- plot_dat$value
  plot_dat$axis_range_value[plot_dat$noninformative_endpoint] <- NA_real_

  # Facet labels are deliberately plain strings.  The Unicode superscript avoids
  # fragile expression labelling in generated TeX/graphics workflows.
  plot_dat$panel_label <- ifelse(plot_dat$panel == "R2", "R²", "RMSE / raw")
  plot_dat$panel_label <- factor(plot_dat$panel_label, levels = c("R²", "RMSE / raw"))

  # Symbol semantics: the raw metric is the current reference-support comparison
  # and therefore receives a fixed circular bullet.  Tier 1 is a latent
  # reference-support adjustment.  Tier 2 is block-/pixel-support adjusted and
  # uses support-size symbols.
  plot_dat$symbol_target <- as.character(plot_dat$support_class)
  plot_dat$symbol_target[plot_dat$Tier == "Raw reference support"] <- "raw reference-support metric"
  plot_dat$symbol_target[plot_dat$Tier == "Tier 1"] <- "latent reference-support metric"
  plot_dat$symbol_target <- factor(
    plot_dat$symbol_target,
    levels = c("raw reference-support metric", "latent reference-support metric",
               "small support", "medium support", "large support", "not reported")
  )

  ensure_dir(dirname(out_csv))
  utils::write.csv(plot_dat, out_csv, row.names = FALSE)
  msg("Wrote plot data: ", out_csv, level = 1)
  plot_dat
}

# A tiny infix helper; avoids importing rlang.
`%||%` <- function(a, b) if (!is.null(a)) a else b

make_posthoc_summary_figure <- function(table_csv = "results/examples/posthoc_examples_table.csv",
                                        out_pdf = "results/examples/fig_posthoc_examples_summary.pdf",
                                        out_png = "results/examples/fig_posthoc_examples_summary.png",
                                        overwrite = OVERWRITE,
                                        verbose = VERBOSE) {
  if (file.exists(out_pdf) && file.exists(out_png) && !overwrite) {
    msg("Keeping existing figure files.", level = 1)
    return(invisible(c(out_pdf, out_png)))
  }
  library(ggplot2)

  plot_dat <- make_posthoc_plot_data(table_csv = table_csv, overwrite = TRUE, verbose = verbose)

  # Enforce alphabetical order at the plotting stage as well, independent of
  # the order in which examples were read or combined.  This is deliberately
  # based on the character display label, not on an inherited factor order.
  y_levels_bottom <- rev(sort(unique(as.character(plot_dat$example_label_text)), na.last = NA))
  plot_dat$example_label <- factor(as.character(plot_dat$example_label_text), levels = y_levels_bottom)

  # Panel-specific plotting limits are based on all finite displayed metric
  # values. Non-informative Tier 2 diagnostics are excluded from scaling but
  # retained as arrows to the panel edge.
  finite_for_limits <- plot_dat[is.finite(plot_dat$axis_range_value), , drop = FALSE]
  panel_limits <- lapply(split(finite_for_limits, finite_for_limits$panel_label), function(d) {
    rng <- range(d$axis_range_value, finite = TRUE)
    if (!all(is.finite(rng))) rng <- c(0, 1)
    pad <- diff(rng) * 0.04
    if (!is.finite(pad) || pad == 0) pad <- max(abs(rng), 1) * 0.04
    data.frame(panel_label = d$panel_label[1], xmin_plot = rng[1] - pad, xmax_plot = rng[2] + pad)
  })
  panel_limits <- do.call(rbind, panel_limits)
  panel_limits$example_label <- factor(y_levels_bottom[1], levels = y_levels_bottom)

  get_panel_limit <- function(panel, which) {
    z <- panel_limits[as.character(panel_limits$panel_label) == as.character(panel), which]
    if (length(z) == 0 || !is.finite(z[1])) return(NA_real_)
    z[1]
  }

  # Scenario range from Tier 1 to Tier 2 in each panel.  The raw value is
  # plotted as the reference-support metric but is not part of the interval
  # segment.
  split_key <- paste(plot_dat$row_id, plot_dat$panel_label, sep = "__")
  pieces <- split(plot_dat, split_key)
  seg <- lapply(pieces, function(d) {
    tier1 <- d$value[d$Tier == "Tier 1"]
    tier2 <- d$value[d$Tier == "Tier 2"]
    if (length(tier1) == 0 || length(tier2) == 0 ||
        !is.finite(tier1[1]) || !is.finite(tier2[1])) return(NULL)

    tier2_row <- d[d$Tier == "Tier 2" & is.finite(d$value), , drop = FALSE][1, , drop = FALSE]
    # Draw arrows for non-informative Tier 2 endpoints rather than plotting
    # them as finite adjusted estimates.
    noninformative <- isTRUE(tier2_row$noninformative_endpoint[1])
    panel_chr <- as.character(d$panel_label[1])
    if (panel_chr == "R²") {
      xend_plot <- if (noninformative) get_panel_limit(d$panel_label[1], "xmax_plot") else tier2[1]
    } else {
      xend_plot <- if (noninformative) get_panel_limit(d$panel_label[1], "xmin_plot") else tier2[1]
    }

    data.frame(
      row_id = d$row_id[1],
      panel_label = d$panel_label[1],
      x = tier1[1],
      xend = tier2[1],
      xend_plot = xend_plot,
      noninformative_endpoint = noninformative,
      example_label = d$example_label[1],
      stringsAsFactors = FALSE
    )
  })
  seg <- do.call(rbind, seg)
  if (!is.null(seg) && nrow(seg) > 0) {
    seg$example_label <- factor(as.character(seg$example_label), levels = y_levels_bottom)
    seg$panel_label <- factor(as.character(seg$panel_label), levels = levels(plot_dat$panel_label))
  }

  baseline_df <- data.frame(panel_label = factor("RMSE / raw", levels = levels(plot_dat$panel_label)), xintercept = 1)

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = point_value, y = example_label)) +
    ggplot2::geom_blank(data = panel_limits,
                        ggplot2::aes(x = xmin_plot, y = example_label),
                        inherit.aes = FALSE) +
    ggplot2::geom_blank(data = panel_limits,
                        ggplot2::aes(x = xmax_plot, y = example_label),
                        inherit.aes = FALSE) +
    ggplot2::geom_vline(data = baseline_df,
                        ggplot2::aes(xintercept = xintercept),
                        inherit.aes = FALSE,
                        linewidth = 0.35, colour = "grey35") +
    ggplot2::geom_segment(data = seg[!seg$noninformative_endpoint, , drop = FALSE],
                          ggplot2::aes(x = x, xend = xend_plot, y = example_label, yend = example_label),
                          inherit.aes = FALSE,
                          linewidth = 0.8,
                          colour = "grey45") +
    ggplot2::geom_segment(data = seg[seg$noninformative_endpoint, , drop = FALSE],
                          ggplot2::aes(x = x, xend = xend_plot, y = example_label, yend = example_label),
                          inherit.aes = FALSE,
                          linewidth = 0.8,
                          colour = "grey45",
                          arrow = grid::arrow(length = grid::unit(0.09, "in"), type = "closed")) +
    # Draw the raw reference-support points first and larger. They mark the
    # reference-support baseline while allowing Tier 1/Tier 2 symbols to remain
    # visible on top when values overlap.
    ggplot2::geom_point(
      data = plot_dat[plot_dat$Tier == "Raw reference support" & is.finite(plot_dat$point_value), , drop = FALSE],
      ggplot2::aes(colour = Tier, fill = Tier, shape = symbol_target),
      size = 4.8,
      stroke = 0.95
    ) +
    ggplot2::geom_point(
      data = plot_dat[plot_dat$Tier != "Raw reference support" & is.finite(plot_dat$point_value), , drop = FALSE],
      ggplot2::aes(colour = Tier, fill = Tier,
                   shape = symbol_target, size = symbol_target),
      stroke = 0.85
    ) +
    ggplot2::facet_wrap(~ panel_label, scales = "free_x", nrow = 1) +
    ggplot2::scale_y_discrete(limits = y_levels_bottom, drop = FALSE) +
    ggplot2::scale_colour_manual(values = c("Raw reference support" = "black",
                                            "Tier 1" = "#D55E00",
                                            "Tier 2" = "#0072B2")) +
    ggplot2::scale_fill_manual(values = c("Raw reference support" = "black",
                                          "Tier 1" = "#D55E00",
                                          "Tier 2" = "#0072B2")) +
    ggplot2::scale_shape_manual(values = c("raw reference-support metric" = 16,
                                           "latent reference-support metric" = 21,
                                           "small support" = 22,
                                           "medium support" = 23,
                                           "large support" = 24,
                                           "not reported" = 4),
                                drop = FALSE) +
    ggplot2::scale_size_manual(values = c("raw reference-support metric" = 4.8,
                                          "latent reference-support metric" = 2.9,
                                          "small support" = 2.9,
                                          "medium support" = 3.6,
                                          "large support" = 4.3,
                                          "not reported" = 2.7),
                               drop = FALSE) +
    ggplot2::labs(x = NULL,
                  y = NULL,
                  colour = "Metric",
                  fill = "Metric",
                  shape = "Target/support",
                  size = "Target/support") +
    ggplot2::guides(
      fill = "none",
      colour = ggplot2::guide_legend(order = 1, override.aes = list(shape = 16, size = 4.8)),
      shape = ggplot2::guide_legend(order = 2),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(5.5, 8, 5.5, 5.5)
    )

  ensure_dir(dirname(out_pdf))
  ggplot2::ggsave(out_pdf, p, width = 8.2, height = 5.8, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(out_png, p, width = 8.2, height = 5.8, units = "in", dpi = 320)
  msg("Wrote figure: ", out_pdf, level = 1)
  msg("Wrote figure: ", out_png, level = 1)

  invisible(c(out_pdf, out_png))
}

make_posthoc_caption_file <- function(out_tex = "results/examples/posthoc_examples_summary_figure_caption.tex",
                                      overwrite = OVERWRITE) {
  if (file.exists(out_tex) && !overwrite) return(invisible(out_tex))
  lines <- c(
    "% Automatically generated by R/30_make_posthoc_table_and_figure.R",
    "\\begin{figure}[ht]",
    "\\centering",
    "\\includegraphics[width=\\textwidth]{../code_data/results/examples/fig_posthoc_examples_summary.pdf}",
    "\\caption{Post-hoc reinterpretation of published or example validation metrics under the proposed information tiers. Black: raw reference-support metric; orange: Tier~1 reference-error adjustment; blue: Tier~2 support adjustment using covariance or representativeness information. Grey segments connect Tier~1 and Tier~2 values; arrows indicate non-informative Tier~2 endpoints outside the feasible metric range. Left: reported or reconstructed $R^2$; predictive/test-set $R^2$ is adjusted on the MSE scale before back-transformation. Right: RMSE or ubRMSE divided by its raw value; the vertical line marks the raw reference-support baseline.}",
    "\\label{fig:posthoc-examples-summary}",
    "\\end{figure}"
  )
  ensure_dir(dirname(out_tex))
  writeLines(lines, out_tex)
  invisible(out_tex)
}

make_posthoc_outputs <- function(combined_csv = "results/examples/posthoc_examples_combined.csv",
                                 table_csv = "results/examples/posthoc_examples_table.csv",
                                 out_dir = "results/examples",
                                 overwrite = OVERWRITE,
                                 verbose = VERBOSE) {
  ensure_dir(out_dir)
  make_posthoc_overview_table(
    combined_csv = combined_csv,
    out_tex = file.path(out_dir, "posthoc_examples_overview_table.tex"),
    overwrite = overwrite,
    verbose = verbose
  )
  make_posthoc_summary_figure(
    table_csv = table_csv,
    out_pdf = file.path(out_dir, "fig_posthoc_examples_summary.pdf"),
    out_png = file.path(out_dir, "fig_posthoc_examples_summary.png"),
    overwrite = overwrite,
    verbose = verbose
  )
  make_posthoc_caption_file(
    out_tex = file.path(out_dir, "posthoc_examples_summary_figure.tex"),
    overwrite = overwrite
  )
  invisible(TRUE)
}

# Run when the script is sourced interactively from the project root.
if (file.exists("results/examples/posthoc_examples_combined.csv") &&
    file.exists("results/examples/posthoc_examples_table.csv")) {
  make_posthoc_outputs(overwrite = OVERWRITE, verbose = VERBOSE)
} else {
  msg("Post-hoc example summary files not found; run R/29_run_all_examples.R first.", level = 1)
}
