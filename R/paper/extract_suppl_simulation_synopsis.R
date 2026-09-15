# Extract compact manuscript-ready simulation summaries.
#
# Run from the manuscript or project root after the simulation has been
# summarized. Adjust RESULT_DIR if your analysis directory differs.

RESULT_DIR <- file.path("results", "simulation")
SUMMARY_FILE <- file.path(RESULT_DIR, "simulation_summary.csv")
OUT_CSV <- file.path(RESULT_DIR, "simulation_synopsis_for_ms.csv")
OUT_TEX <- file.path(RESULT_DIR, "simulation_synopsis_for_ms.tex")

stopifnot(file.exists(SUMMARY_FILE))

x <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE)

# Core settings used in the main paper figures and selected table.
core_mid <- subset(
  x,
  scenario_type == "core" &
    prediction_relation == "block_targeted" &
    abs(target_r2_block - 0.60) < 1e-12
)

selected <- subset(
  core_mid,
  (range_block_ratio == 2  & q_epsilon_target %in% c(0, 0.40)) |
    (range_block_ratio == 10 & q_epsilon_target == 0.25) |
    (range_block_ratio == 50 & q_epsilon_target %in% c(0, 0.40))
)

keep <- c(
  "q_epsilon_target", "range_block_ratio", "lambda_theoretical_mean",
  "R2_obs_mean", "R2_ME_mean", "R2_B_hat_mean", "R2_block_true_mean",
  "RMSE_obs_mean", "RMSE_ME_mean", "RMSE_B_hat_mean", "RMSE_block_true_mean"
)
selected <- selected[order(selected$range_block_ratio, selected$q_epsilon_target), keep]
write.csv(selected, OUT_CSV, row.names = FALSE)

fmt <- function(z) sprintf("%.3f", z)
rows <- apply(selected, 1, function(r) {
  paste(
    fmt(as.numeric(r["q_epsilon_target"])),
    as.integer(as.numeric(r["range_block_ratio"])),
    fmt(as.numeric(r["lambda_theoretical_mean"])),
    fmt(as.numeric(r["R2_obs_mean"])),
    fmt(as.numeric(r["R2_ME_mean"])),
    fmt(as.numeric(r["R2_B_hat_mean"])),
    fmt(as.numeric(r["RMSE_obs_mean"])),
    fmt(as.numeric(r["RMSE_ME_mean"])),
    fmt(as.numeric(r["RMSE_B_hat_mean"])),
    sep = " & "
  )
})

latex <- c(
  "\\begin{tabular}{rrrrrrrrr}",
  "\\toprule",
  "$q_\\varepsilon$ & range/block & $\\lambda$ & $R^2_{obs}$ & $R^2_{ME}$ & $\\hat R^2_B$ & RMSE$_{obs}$ & RMSE$_{ME}$ & $\\widehat{\\mathrm{RMSE}}_B$ \\\\",
  "\\midrule",
  paste0(rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(latex, OUT_TEX)

message("Wrote: ", OUT_CSV)
message("Wrote: ", OUT_TEX)
