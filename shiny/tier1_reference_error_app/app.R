# Tier 1 reference-error adjustment app
# Support-adjusted validation metrics in environmental prediction

library(shiny)

num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(NA_real_)
  suppressWarnings(as.numeric(x))
}

fmt <- function(x, digits = 3) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("--")
  if (!is.finite(x)) return("n.i.")
  formatC(x, digits = digits, format = "fg", flag = "#")
}

fmt_status <- function(ok, msg_ok = "finite", msg_bad = "non-informative") {
  if (isTRUE(ok)) msg_ok else msg_bad
}

metric_to_components <- function(metric_type, metric_value, sigma2_obs, rmse_obs) {
  out <- list(raw_r = NA_real_, raw_r2_cor = NA_real_, raw_r2_pred = NA_real_, mse_obs = NA_real_)
  if (is.na(metric_value)) return(out)
  if (metric_type == "correlation") {
    out$raw_r <- metric_value
    out$raw_r2_cor <- metric_value^2
  } else if (metric_type == "r2_cor") {
    out$raw_r2_cor <- metric_value
    out$raw_r <- if (metric_value >= 0) sqrt(metric_value) else NA_real_
  } else if (metric_type == "r2_pred") {
    out$raw_r2_pred <- metric_value
    if (!is.na(sigma2_obs)) out$mse_obs <- (1 - metric_value) * sigma2_obs
  }
  if (!is.na(rmse_obs)) out$mse_obs <- rmse_obs^2
  out
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      body { font-size: 13px; }\n      .container-fluid { max-width: 1320px; }\n      .panel { margin-bottom: 8px; }\n      .panel-heading { padding: 6px 10px; font-weight: 600; }\n      .panel-body { padding: 8px 10px; }\n      .form-group { margin-bottom: 7px; }\n      .radio { margin-top: 2px; margin-bottom: 2px; }\n      .sd-field .control-label { color: #1f5f99; font-weight: 600; }\n      .sd-field input { background-color: #eef6ff; border-color: #8dbbe5; }\n      .var-field .control-label { color: #7a4a00; font-weight: 600; }\n      .var-field input { background-color: #fff5df; border-color: #d4a64a; }\n      .metric-field .control-label { font-weight: 600; }\n      .equation { font-family: Georgia, serif; background: #f7f7f7; border-left: 3px solid #777; padding: 6px 8px; margin-bottom: 6px; }\n      .resultbox { background: #fbfbfb; border: 1px solid #ddd; border-radius: 4px; padding: 7px 9px; margin-bottom: 7px; }\n      .resultbox strong { display: inline-block; min-width: 170px; }\n      .diagnostic { font-size: 12px; color: #333; background: #fff8e8; border: 1px solid #e1c77a; padding: 7px; border-radius: 4px; }\n      .notice { font-size: 11px; color: #666; margin-top: 10px; border-top: 1px solid #ddd; padding-top: 6px; }\n    "))
  ),
  titlePanel("Tier 1 reference-error adjustment"),
  fluidRow(
    column(
      width = 6,
      wellPanel(
        tags$h5("Example inputs"),
        fluidRow(
          column(4, actionButton("ex_biomass", "Biomass", width = "100%")),
          column(4, actionButton("ex_bts", "BTS", width = "100%")),
          column(4, actionButton("ex_swe", "SWE", width = "100%"))
        ),
        fluidRow(
          column(4, actionButton("ex_no2eu", "NO2-EU split", width = "100%")),
          column(4, actionButton("ex_smap", "SMAP-CRNS", width = "100%")),
          column(4, actionButton("clear_inputs", "Reset", width = "100%"))
        ),

        radioButtons(
          "metric_type", "Reported validation metric",
          choices = c("correlation r" = "correlation", "R²_cor" = "r2_cor", "R²_pred" = "r2_pred"),
          selected = "r2_cor", inline = TRUE
        ),
        div(class = "metric-field", textInput("metric_value", HTML("Value of selected metric"), value = "0.50")),
        textInput("rmse_obs", HTML("Raw RMSE or ubRMSE (optional)"), value = ""),
        fluidRow(
          column(6, div(class = "sd-field", textInput("obs_sd", HTML("Observed reference SD &sigma;<sub>obs</sub>"), value = "1"))),
          column(6, div(class = "var-field", textInput("obs_var", HTML("Observed reference variance &sigma;<sup>2</sup><sub>obs</sub>"), value = "1")))
        ),
        fluidRow(
          column(6, div(class = "sd-field", textInput("eps_sd", HTML("Measurement-error SD &sigma;<sub>&epsilon;</sub>"), value = "0"))),
          column(6, div(class = "var-field", textInput("eps_var", HTML("Measurement-error variance &sigma;<sup>2</sup><sub>&epsilon;</sub>"), value = "0")))
        )
      ),
      div(class = "notice",
          HTML("&copy; Alexander Brenning and Thomas Suesse. License: GPL-3 or later.<br>Cite as: Brenning and Suesse, <em>Support-adjusted validation metrics in environmental prediction</em> (manuscript)."))
    ),
    column(
      width = 6,
      tags$h4("Equations"),
      div(class = "equation", HTML("q<sub>&epsilon;</sub> = &sigma;<sup>2</sup><sub>&epsilon;</sub> / &sigma;<sup>2</sup><sub>obs</sub>")),
      div(class = "equation", HTML("r<sub>ME</sub> = r<sub>obs</sub> / sqrt(1 - q<sub>&epsilon;</sub>)")),
      div(class = "equation", HTML("R²<sub>cor,ME</sub> = R²<sub>cor,obs</sub> / (1 - q<sub>&epsilon;</sub>)")),
      div(class = "equation", HTML("RMSE<sub>ME</sub> = sqrt(RMSE²<sub>obs</sub> - &sigma;<sup>2</sup><sub>&epsilon;</sub>)")),
      div(class = "equation", HTML("R²<sub>pred,ME</sub> = 1 - (MSE<sub>obs</sub> - &sigma;<sup>2</sup><sub>&epsilon;</sub>) / (&sigma;<sup>2</sup><sub>obs</sub> - &sigma;<sup>2</sup><sub>&epsilon;</sub>)")),
      tags$h4("Results"),
      uiOutput("results"),
      tags$h4("Diagnostics"),
      uiOutput("diagnostics")
    )
  )
)

server <- function(input, output, session) {
  set_inputs <- function(metric_type, metric_value, rmse, obs_var, eps_var) {
    updateRadioButtons(session, "metric_type", selected = metric_type)
    updateTextInput(session, "metric_value", value = as.character(metric_value))
    updateTextInput(session, "rmse_obs", value = as.character(rmse))
    updateTextInput(session, "obs_var", value = as.character(obs_var))
    updateTextInput(session, "obs_sd", value = as.character(signif(sqrt(as.numeric(obs_var)), 8)))
    updateTextInput(session, "eps_var", value = as.character(eps_var))
    updateTextInput(session, "eps_sd", value = as.character(signif(sqrt(as.numeric(eps_var)), 8)))
  }
  observeEvent(input$ex_biomass, set_inputs("r2_pred", 0.438, 992, 1749447, 300724))
  observeEvent(input$ex_bts, set_inputs("r2_cor", 0.451, 1.82, 4.19, 0.05))
  observeEvent(input$ex_swe, set_inputs("r2_cor", 0.414, 125.5, 26896, 400))
  observeEvent(input$ex_no2eu, set_inputs("r2_pred", 0.553, 9.5258, 203, 56.4))
  observeEvent(input$ex_smap, set_inputs("correlation", 0.699, 0.056, "", 0.000705))
  observeEvent(input$clear_inputs, set_inputs("r2_cor", 0.50, "", 1, 0))

  last_pair <- reactiveVal(NULL)
  observeEvent(input$obs_sd, { last_pair("obs_sd") }, ignoreInit = TRUE)
  observeEvent(input$obs_var, { last_pair("obs_var") }, ignoreInit = TRUE)
  observeEvent(input$eps_sd, { last_pair("eps_sd") }, ignoreInit = TRUE)
  observeEvent(input$eps_var, { last_pair("eps_var") }, ignoreInit = TRUE)

  sync_pair <- function(sd_id, var_id) {
    sd_react <- debounce(reactive(input[[sd_id]]), 900)
    var_react <- debounce(reactive(input[[var_id]]), 900)
    observeEvent(sd_react(), {
      if (!identical(last_pair(), sd_id)) return()
      x <- num_or_na(input[[sd_id]])
      if (!is.na(x) && x >= 0) updateTextInput(session, var_id, value = format(signif(x^2, 8), scientific = FALSE))
    }, ignoreInit = TRUE)
    observeEvent(var_react(), {
      if (!identical(last_pair(), var_id)) return()
      x <- num_or_na(input[[var_id]])
      if (!is.na(x) && x >= 0) updateTextInput(session, sd_id, value = format(signif(sqrt(x), 8), scientific = FALSE))
    }, ignoreInit = TRUE)
  }
  sync_pair("obs_sd", "obs_var")
  sync_pair("eps_sd", "eps_var")

  vals <- reactive({
    sigma2_obs <- num_or_na(input$obs_var)
    sigma2_eps <- num_or_na(input$eps_var)
    rmse_obs <- num_or_na(input$rmse_obs)
    metric_value <- num_or_na(input$metric_value)
    comp <- metric_to_components(input$metric_type, metric_value, sigma2_obs, rmse_obs)
    q_eps <- if (!is.na(sigma2_obs) && sigma2_obs > 0 && !is.na(sigma2_eps)) sigma2_eps / sigma2_obs else NA_real_
    ok_q <- !is.na(q_eps) && q_eps >= 0 && q_eps < 1

    r_me <- r2_cor_me <- r2_pred_me <- rmse_me <- NA_real_
    r_status <- r2cor_status <- r2pred_status <- rmse_status <- "not computed"

    if (ok_q && !is.na(comp$raw_r)) {
      r_me_raw <- comp$raw_r / sqrt(1 - q_eps)
      if (abs(r_me_raw) <= 1) { r_me <- r_me_raw; r_status <- "finite" } else r_status <- "n.i.: adjusted correlation exceeds [-1, 1]"
    }
    if (ok_q && !is.na(comp$raw_r2_cor)) {
      r2_raw <- comp$raw_r2_cor / (1 - q_eps)
      if (r2_raw >= 0 && r2_raw <= 1) { r2_cor_me <- r2_raw; r2cor_status <- "finite" } else r2cor_status <- "n.i.: adjusted R²_cor exceeds feasible range"
    }
    if (ok_q && !is.na(rmse_obs)) {
      mse_me <- rmse_obs^2 - sigma2_eps
      if (mse_me >= 0) { rmse_me <- sqrt(mse_me); rmse_status <- "finite" } else rmse_status <- "n.i.: measurement-error variance exceeds observed squared error"
    }
    if (ok_q && !is.na(comp$mse_obs) && !is.na(sigma2_obs)) {
      mse_me <- comp$mse_obs - sigma2_eps
      var_lat <- sigma2_obs - sigma2_eps
      r2_raw <- 1 - mse_me / var_lat
      if (mse_me >= 0 && var_lat > 0 && r2_raw <= 1) { r2_pred_me <- r2_raw; r2pred_status <- "finite" } else r2pred_status <- "n.i.: MSE or latent reference variance is not feasible"
    }

    list(sigma2_obs=sigma2_obs, sigma2_eps=sigma2_eps, q_eps=q_eps, ok_q=ok_q,
         raw=comp, rmse_obs=rmse_obs, r_me=r_me, r2_cor_me=r2_cor_me,
         r2_pred_me=r2_pred_me, rmse_me=rmse_me, r_status=r_status,
         r2cor_status=r2cor_status, r2pred_status=r2pred_status, rmse_status=rmse_status)
  })

  output$results <- renderUI({
    v <- vals()
    tagList(
      div(class="resultbox", HTML(paste0("<strong>q<sub>&epsilon;</sub></strong>", fmt(v$q_eps, 4)))),
      div(class="resultbox", HTML(paste0("<strong>Raw r</strong>", fmt(v$raw$raw_r)))),
      div(class="resultbox", HTML(paste0("<strong>Raw R²<sub>cor</sub></strong>", fmt(v$raw$raw_r2_cor)))),
      div(class="resultbox", HTML(paste0("<strong>Raw R²<sub>pred</sub></strong>", fmt(v$raw$raw_r2_pred)))),
      div(class="resultbox", HTML(paste0("<strong>Tier 1 r<sub>ME</sub></strong>", fmt(v$r_me), " &nbsp; <em>", v$r_status, "</em>"))),
      div(class="resultbox", HTML(paste0("<strong>Tier 1 R²<sub>cor,ME</sub></strong>", fmt(v$r2_cor_me), " &nbsp; <em>", v$r2cor_status, "</em>"))),
      div(class="resultbox", HTML(paste0("<strong>Tier 1 R²<sub>pred,ME</sub></strong>", fmt(v$r2_pred_me), " &nbsp; <em>", v$r2pred_status, "</em>"))),
      div(class="resultbox", HTML(paste0("<strong>Tier 1 RMSE<sub>ME</sub></strong>", fmt(v$rmse_me), " &nbsp; <em>", v$rmse_status, "</em>")))
    )
  })

  output$diagnostics <- renderUI({
    v <- vals()
    msg <- c()
    if (is.na(v$sigma2_obs) || v$sigma2_obs <= 0) msg <- c(msg, "Observed reference variance must be positive.")
    if (is.na(v$sigma2_eps) || v$sigma2_eps < 0) msg <- c(msg, "Measurement-error variance must be non-negative.")
    if (!is.na(v$q_eps) && v$q_eps >= 1) msg <- c(msg, "q_epsilon must be smaller than one. Otherwise the latent reference variance is not positive.")
    msg <- c(msg, "n.i. means non-informative: the supplied inputs imply an adjusted endpoint outside the feasible metric range. Such endpoints should not be truncated to R² = 1 or RMSE = 0.")
    div(class="diagnostic", HTML(paste(msg, collapse = "<br>")))
  })
}

shinyApp(ui, server)
