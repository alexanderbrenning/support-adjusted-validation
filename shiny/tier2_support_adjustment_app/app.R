# Tier 2 / Tier 3 reference-to-block support adjustment app
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

metric_to_components <- function(metric_type, metric_value, sigma2_obs, rmse_obs) {
  out <- list(raw_r = NA_real_, raw_r2_cor = NA_real_, raw_r2_pred = NA_real_, mse_obs = NA_real_)
  if (!is.na(metric_value)) {
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
  }
  if (!is.na(rmse_obs)) out$mse_obs <- rmse_obs^2
  out
}

rho_spherical <- function(h, a) {
  if (is.na(h) || is.na(a) || a <= 0) return(NA_real_)
  x <- h / a
  if (x >= 1) return(0)
  1 - (1.5 * x - 0.5 * x^3)
}

rho_exponential_practical <- function(h, a) {
  # a is interpreted as an effective/practical range at correlation about 0.05.
  if (is.na(h) || is.na(a) || a <= 0) return(NA_real_)
  exp(-3 * h / a)
}

support_from_svgm <- function(sigma2_obs, sigma2_eps, sigma2_mu, sigma2_S, L, range, model) {
  R <- L / sqrt(pi)
  if (is.na(sigma2_S) || sigma2_S < 0) {
    return(list(ok=FALSE, msg="negative structural partial sill", sigma2_S=sigma2_S,
                R=R, rho=NA_real_, delta=NA_real_, lambda=NA_real_, var_block=NA_real_))
  }
  rho <- if (model == "spherical") rho_spherical(R, range) else rho_exponential_practical(R, range)
  k <- 128 / (45 * pi)
  delta <- sigma2_mu + k * sigma2_S * (1 - rho)
  var_block <- sigma2_obs - sigma2_eps - delta
  lambda <- var_block / sigma2_obs
  remainder <- sigma2_obs - sigma2_eps - sigma2_mu - sigma2_S
  list(ok=TRUE, msg="finite support inputs", sigma2_S=sigma2_S, R=R, rho=rho,
       delta=delta, lambda=lambda, var_block=var_block, variance_remainder=remainder)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      body { font-size: 13px; }\n      .container-fluid { max-width: 1360px; }\n      .panel { margin-bottom: 8px; }\n      .panel-heading { padding: 6px 10px; font-weight: 600; }\n      .panel-body { padding: 8px 10px; }\n      .form-group { margin-bottom: 7px; }\n      .radio { margin-top: 2px; margin-bottom: 2px; }\n      .sd-field .control-label { color: #1f5f99; font-weight: 600; }\n      .sd-field input { background-color: #eef6ff; border-color: #8dbbe5; }\n      .var-field .control-label { color: #7a4a00; font-weight: 600; }\n      .var-field input { background-color: #fff5df; border-color: #d4a64a; }\n      .calc-field { background: #f0f0f0; border: 1px solid #ccc; padding: 5px 7px; border-radius: 4px; min-height: 34px; }\n      .calc-label { font-weight: 600; color: #555; }\n      .metric-field .control-label { font-weight: 600; }\n      .equation { font-family: Georgia, serif; background: #f7f7f7; border-left: 3px solid #777; padding: 6px 8px; margin-bottom: 6px; }\n      .resultbox { background: #fbfbfb; border: 1px solid #ddd; border-radius: 4px; padding: 7px 9px; margin-bottom: 7px; }\n      .resultbox strong { display: inline-block; min-width: 185px; }\n      .diagnostic { font-size: 12px; color: #333; background: #fff8e8; border: 1px solid #e1c77a; padding: 7px; border-radius: 4px; }\n      .notice { font-size: 11px; color: #666; margin-top: 10px; border-top: 1px solid #ddd; padding-top: 6px; }\n    "))
  ),
  titlePanel("Tier 2 / Tier 3 reference-to-block support adjustment"),
  fluidRow(
    column(
      width = 6,
      wellPanel(
        tags$h5("Example inputs"),
        fluidRow(
          column(3, actionButton("ex_bts", "BTS", width = "100%")),
          column(3, actionButton("ex_swe10", "SWE 10 m", width = "100%")),
          column(3, actionButton("ex_swe100", "SWE 100 m", width = "100%")),
          column(3, actionButton("ex_no2de", "NO2-DE", width = "100%"))
        ),
        fluidRow(
          column(3, actionButton("ex_no2eu_split", "NO2-EU split", width = "100%")),
          column(3, actionButton("ex_no2eu_full", "NO2-EU full", width = "100%")),
          column(3, actionButton("ex_meusezn", "Meuse-Zn", width = "100%")),
          column(3, actionButton("ex_meusecd", "Meuse-Cd", width = "100%"))
        ),

        radioButtons(
          "metric_type", "Reported validation metric",
          choices = c("correlation r" = "correlation", "R²_cor" = "r2_cor", "R²_pred" = "r2_pred"),
          selected = "r2_cor", inline = TRUE
        ),
        div(class = "metric-field", textInput("metric_value", HTML("Value of selected metric"), value = "0.451")),
        textInput("rmse_obs", HTML("Raw RMSE or ubRMSE (optional)"), value = "1.82"),
        fluidRow(
          column(6, div(class = "sd-field", textInput("obs_sd", HTML("Observed reference SD &sigma;<sub>obs</sub>"), value = "2.046948"))),
          column(6, div(class = "var-field", textInput("obs_var", HTML("Observed reference variance &sigma;<sup>2</sup><sub>obs</sub>"), value = "4.19")))
        ),
        fluidRow(
          column(6, div(class = "sd-field", textInput("eps_sd", HTML("Measurement-error SD &sigma;<sub>&epsilon;</sub>"), value = "0.2236068"))),
          column(6, div(class = "var-field", textInput("eps_var", HTML("Measurement-error variance &sigma;<sup>2</sup><sub>&epsilon;</sub>"), value = "0.05")))
        )
      ),
      wellPanel(
        tags$h5("Support inputs"),
        radioButtons(
          "support_mode", "Support information",
          choices = c("Enter λ_A,B and Δ_A,B directly" = "direct",
                      "Use semivariogram support approximation" = "svgm"),
          selected = "svgm"
        ),
        conditionalPanel(
          condition = "input.support_mode == 'direct'",
          textInput("lambda_direct", HTML("Representativeness factor &lambda;<sub>A,B</sub>"), value = "0.806"),
          div(class = "var-field", textInput("delta_direct", HTML("Support-variance difference &Delta;<sub>A,B</sub>"), value = "0.763")),
          helpText("Use this mode when λ_A,B and Δ_A,B come from numerical integration, a reported representativeness-error estimate, or another component-specific calculation. The example buttons use semivariogram mode whenever the corresponding parameter inputs are available.")
        ),
        conditionalPanel(
          condition = "input.support_mode == 'svgm'",
          fluidRow(
            column(6, div(class = "sd-field", textInput("mu_sd", HTML("True microscale SD &sigma;<sub>&mu;</sub>"), value = "0.6324555"))),
            column(6, div(class = "var-field", textInput("mu_var", HTML("True microscale variance &sigma;<sup>2</sup><sub>&mu;</sub>"), value = "0.40")))
          ),
          fluidRow(
            column(6, div(class = "sd-field", textInput("sigmaS_sd_in", HTML("Structural partial-sill SD &sigma;<sub>S</sub>"), value = "1.688194"))),
            column(6, div(class = "var-field", textInput("sigmaS_var_in", HTML("Structural partial sill &sigma;<sup>2</sup><sub>S</sub>"), value = "2.85")))
          ),
          fluidRow(
            column(6, textInput("block_L", HTML("Block side <em>L</em> (target resolution)"), value = "30")),
            column(6, textInput("range_a", HTML("Structural range <em>a</em>"), value = "180"))
          ),
          radioButtons("svgm_model", "Structural covariance model", choices = c("spherical" = "spherical", "exponential, practical range" = "exponential"), selected = "spherical", inline = TRUE),
          helpText("Enter the structural partial sill that belongs to the structural range. It may be the partial sill of a response semivariogram or of a residual/component semivariogram; it is not inferred from the observed reference variance.")
        )
      )
    ),
    column(
      width = 6,
      fluidRow(
        column(6,
          tags$h4("Equations"),
          div(class = "equation", HTML("&lambda;<sub>A,B</sub> = Var{Z<sub>B</sub>(B<sub>i</sub>)} / &sigma;<sup>2</sup><sub>obs</sub>")),
          div(class = "equation", HTML("R²<sub>cor,B</sub> &asymp; R²<sub>cor,obs</sub> / &lambda;<sub>A,B</sub>")),
          div(class = "equation", HTML("MSE<sub>B</sub> &asymp; MSE<sub>obs</sub> - &Delta;<sub>A,B</sub> - &sigma;<sup>2</sup><sub>&epsilon;</sub>")),
          div(class = "equation", HTML("R²<sub>pred,B</sub> &asymp; 1 - MSE<sub>B</sub> / Var{Z<sub>B</sub>(B<sub>i</sub>)}")),
          div(class = "equation", HTML("&Delta;<sub>p,B</sub> &asymp; &sigma;<sup>2</sup><sub>&mu;</sub> + 128/(45&pi;) &sigma;<sup>2</sup><sub>S</sub>{1 - &rho;<sub>S</sub>(R)}"))
        ),
        column(6,
          tags$h4("Semivariogram"),
          plotOutput("svgm_plot", height = "250px")
        )
      ),
      fluidRow(
        column(6,
          tags$h4("Results"),
          uiOutput("results")
        ),
        column(6,
          tags$h4("Diagnostics"),
          uiOutput("diagnostics"),
          div(class = "notice",
              HTML("&copy; Alexander Brenning and Thomas Suesse. License: GPL-3 or later.<br>Cite as: Brenning and Suesse, <em>Support-adjusted validation metrics in environmental prediction</em> (manuscript)."))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  set_common_example <- function(metric_type, metric_value, rmse, obs_var, eps_var) {
    updateRadioButtons(session, "metric_type", selected = metric_type)
    updateTextInput(session, "metric_value", value = as.character(metric_value))
    updateTextInput(session, "rmse_obs", value = as.character(rmse))
    updateTextInput(session, "obs_var", value = as.character(obs_var))
    updateTextInput(session, "obs_sd", value = ifelse(is.na(suppressWarnings(as.numeric(obs_var))), "", as.character(signif(sqrt(as.numeric(obs_var)), 8))))
    updateTextInput(session, "eps_var", value = as.character(eps_var))
    updateTextInput(session, "eps_sd", value = ifelse(is.na(suppressWarnings(as.numeric(eps_var))), "", as.character(signif(sqrt(as.numeric(eps_var)), 8))))
  }

  set_svgm_example <- function(metric_type, metric_value, rmse, obs_var, eps_var,
                               mu_var, sigmaS_var, block_L, range_a,
                               svgm_model = "spherical") {
    set_common_example(metric_type, metric_value, rmse, obs_var, eps_var)
    updateRadioButtons(session, "support_mode", selected = "svgm")
    updateTextInput(session, "mu_var", value = as.character(mu_var))
    updateTextInput(session, "mu_sd", value = ifelse(is.na(suppressWarnings(as.numeric(mu_var))), "", as.character(signif(sqrt(as.numeric(mu_var)), 8))))
    updateTextInput(session, "sigmaS_var_in", value = as.character(sigmaS_var))
    updateTextInput(session, "sigmaS_sd_in", value = ifelse(is.na(suppressWarnings(as.numeric(sigmaS_var))), "", as.character(signif(sqrt(as.numeric(sigmaS_var)), 8))))
    updateTextInput(session, "block_L", value = as.character(block_L))
    updateTextInput(session, "range_a", value = as.character(range_a))
    updateRadioButtons(session, "svgm_model", selected = svgm_model)
  }

  set_direct_example <- function(metric_type, metric_value, rmse, obs_var, eps_var, lambda, delta) {
    set_common_example(metric_type, metric_value, rmse, obs_var, eps_var)
    updateRadioButtons(session, "support_mode", selected = "direct")
    updateTextInput(session, "lambda_direct", value = as.character(lambda))
    updateTextInput(session, "delta_direct", value = as.character(delta))
  }

  # Examples with published, fixed or dynamically fitted semivariogram inputs.
  observeEvent(input$ex_bts, set_svgm_example("r2_cor", 0.451, 1.82, 4.19, 0.05, 0.40, 2.85, 30, 180, "spherical"))
  observeEvent(input$ex_swe10, set_svgm_example("r2_cor", 0.414, 125.5, 26896, 400, 1200, 25296, 10, 250, "exponential"))
  observeEvent(input$ex_swe100, set_svgm_example("r2_cor", 0.414, 125.5, 26896, 400, 1200, 25296, 100, 250, "exponential"))
  observeEvent(input$ex_no2eu_split, set_svgm_example("r2_pred", 0.553, 9.5258, 203, 56.4, 56.4, 70.2, 100, 342798, "spherical"))
  observeEvent(input$ex_no2eu_full, set_svgm_example("r2_pred", 0.553, 9.5258, 203, 0, 112.8, 70.2, 100, 342798, "spherical"))
  observeEvent(input$ex_no2de, set_svgm_example("r2_cor", 0.308, 11.2, 181, 0, 109, 68.3, 2000, 179324, "spherical"))
  observeEvent(input$ex_meusezn, set_svgm_example("r2_pred", 0.654, 0.42, 0.521, 0, 0.0507, 0.591, 100, 897, "spherical"))
  observeEvent(input$ex_meusecd, set_svgm_example("r2_pred", 0.711, 0.66, 1.50, 0, 0.548, 1.34, 100, 1149, "spherical"))

  last_pair <- reactiveVal(NULL)
  observeEvent(input$obs_sd, { last_pair("obs_sd") }, ignoreInit = TRUE)
  observeEvent(input$obs_var, { last_pair("obs_var") }, ignoreInit = TRUE)
  observeEvent(input$eps_sd, { last_pair("eps_sd") }, ignoreInit = TRUE)
  observeEvent(input$eps_var, { last_pair("eps_var") }, ignoreInit = TRUE)
  observeEvent(input$mu_sd, { last_pair("mu_sd") }, ignoreInit = TRUE)
  observeEvent(input$mu_var, { last_pair("mu_var") }, ignoreInit = TRUE)
  observeEvent(input$sigmaS_sd_in, { last_pair("sigmaS_sd_in") }, ignoreInit = TRUE)
  observeEvent(input$sigmaS_var_in, { last_pair("sigmaS_var_in") }, ignoreInit = TRUE)

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
  sync_pair("mu_sd", "mu_var")
  sync_pair("sigmaS_sd_in", "sigmaS_var_in")

  support_vals <- reactive({
    sigma2_obs <- num_or_na(input$obs_var)
    sigma2_eps <- num_or_na(input$eps_var)
    if (input$support_mode == "direct") {
      lambda <- num_or_na(input$lambda_direct)
      delta <- num_or_na(input$delta_direct)
      var_block <- if (!is.na(lambda) && !is.na(sigma2_obs)) lambda * sigma2_obs else NA_real_
      return(list(mode="direct", lambda=lambda, delta=delta, var_block=var_block,
                  sigma2_S=NA_real_, sigmaS=NA_real_, R=NA_real_, rho=NA_real_, msg="direct support quantities"))
    }
    sigma2_mu <- num_or_na(input$mu_var)
    sigma2_S <- num_or_na(input$sigmaS_var_in)
    L <- num_or_na(input$block_L)
    range <- num_or_na(input$range_a)
    sv <- support_from_svgm(sigma2_obs, sigma2_eps, sigma2_mu, sigma2_S, L, range, input$svgm_model)
    sv$mode <- "svgm"
    sv$sigmaS <- if (!is.na(sv$sigma2_S) && sv$sigma2_S >= 0) sqrt(sv$sigma2_S) else NA_real_
    sv
  })

  vals <- reactive({
    sigma2_obs <- num_or_na(input$obs_var)
    sigma2_eps <- num_or_na(input$eps_var)
    rmse_obs <- num_or_na(input$rmse_obs)
    metric_value <- num_or_na(input$metric_value)
    comp <- metric_to_components(input$metric_type, metric_value, sigma2_obs, rmse_obs)
    sup <- support_vals()
    q_eps <- if (!is.na(sigma2_obs) && sigma2_obs > 0 && !is.na(sigma2_eps)) sigma2_eps / sigma2_obs else NA_real_

    r2_cor_B <- r2_pred_B <- rmse_B <- NA_real_
    r2cor_status <- r2pred_status <- rmse_status <- "not computed"
    if (!is.na(comp$raw_r2_cor) && !is.na(sup$lambda) && sup$lambda > 0) {
      z <- comp$raw_r2_cor / sup$lambda
      if (z >= 0 && z <= 1) { r2_cor_B <- z; r2cor_status <- "finite" } else r2cor_status <- "n.i.: adjusted R²_cor exceeds feasible range"
    }
    if (!is.na(comp$mse_obs) && !is.na(sup$delta) && !is.na(sigma2_eps) && !is.na(sup$var_block)) {
      mse_B <- comp$mse_obs - sup$delta - sigma2_eps
      z <- 1 - mse_B / sup$var_block
      if (mse_B >= 0 && sup$var_block > 0 && z <= 1) { r2_pred_B <- z; r2pred_status <- "finite" } else r2pred_status <- "n.i.: MSE-scale endpoint outside feasible range"
    }
    if (!is.na(rmse_obs) && !is.na(sup$delta) && !is.na(sigma2_eps)) {
      mse_B <- rmse_obs^2 - sup$delta - sigma2_eps
      if (mse_B >= 0) { rmse_B <- sqrt(mse_B); rmse_status <- "finite" } else rmse_status <- "n.i.: adjusted MSE is negative"
    }
    list(sigma2_obs=sigma2_obs, sigma2_eps=sigma2_eps, q_eps=q_eps, raw=comp,
         support=sup, rmse_obs=rmse_obs, r2_cor_B=r2_cor_B, r2_pred_B=r2_pred_B,
         rmse_B=rmse_B, r2cor_status=r2cor_status, r2pred_status=r2pred_status,
         rmse_status=rmse_status)
  })

  output$results <- renderUI({
    v <- vals(); s <- v$support
    tagList(
      div(class="resultbox", HTML(paste0("<strong>q<sub>&epsilon;</sub></strong>", fmt(v$q_eps, 4)))),
      div(class="resultbox", HTML(paste0("<strong>&lambda;<sub>A,B</sub></strong>", fmt(s$lambda, 4)))),
      div(class="resultbox", HTML(paste0("<strong>&Delta;<sub>A,B</sub></strong>", fmt(s$delta, 4)))),
      div(class="resultbox", HTML(paste0("<strong>Var{Z<sub>B</sub>(B<sub>i</sub>)}</strong>", fmt(s$var_block, 4)))),
      div(class="resultbox", HTML(paste0("<strong>Tier 2/3 R²<sub>cor,B</sub></strong>", fmt(v$r2_cor_B), " &nbsp; <em>", v$r2cor_status, "</em>"))),
      div(class="resultbox", HTML(paste0("<strong>Tier 2/3 R²<sub>pred,B</sub></strong>", fmt(v$r2_pred_B), " &nbsp; <em>", v$r2pred_status, "</em>"))),
      div(class="resultbox", HTML(paste0("<strong>Tier 2/3 RMSE<sub>B</sub></strong>", fmt(v$rmse_B), " &nbsp; <em>", v$rmse_status, "</em>")))
    )
  })

  output$diagnostics <- renderUI({
    v <- vals(); s <- v$support
    msg <- c()
    if (is.na(v$sigma2_obs) || v$sigma2_obs <= 0) msg <- c(msg, "Observed reference variance must be positive.")
    if (is.na(v$sigma2_eps) || v$sigma2_eps < 0) msg <- c(msg, "Measurement-error variance must be non-negative.")
    if (!is.na(v$q_eps) && v$q_eps >= 1) msg <- c(msg, "q_epsilon must be smaller than one.")
    if (identical(s$mode, "svgm")) {
      msg <- c(msg, paste0("Entered structural partial sill σ²_S = ", fmt(s$sigma2_S, 4), "."))
      if (!is.na(s$sigma2_S) && s$sigma2_S < 0) msg <- c(msg, "Structural partial sill must be non-negative.")
      if (!is.na(s$variance_remainder) && abs(s$variance_remainder) > 1e-6) {
        msg <- c(msg, paste0("Variance remainder σ²_obs - σ²_ε - σ²_μ - σ²_S = ", fmt(s$variance_remainder, 4), ". A non-zero value is acceptable when the semivariogram is residual- or component-specific."))
      }
      msg <- c(msg, "Use direct λ and Δ entry when these support quantities have already been computed by numerical integration, a reported SRE, or the example-specific R workflow.")
    }
    msg <- c(msg, "n.i. means non-informative: the supplied inputs imply an adjusted endpoint outside the feasible metric range. Do not truncate such endpoints to R² = 1 or RMSE = 0.")
    div(class="diagnostic", HTML(paste(msg, collapse = "<br>")))
  })

  output$svgm_plot <- renderPlot({
    if (input$support_mode != "svgm") {
      plot.new(); text(0.5, 0.5, "Semivariogram plot shown in semivariogram mode", cex = 0.9); return()
    }
    sigma2_obs <- num_or_na(input$obs_var)
    sigma2_eps <- num_or_na(input$eps_var)
    sigma2_mu <- num_or_na(input$mu_var)
    sigma2_S <- num_or_na(input$sigmaS_var_in)
    L <- num_or_na(input$block_L)
    range <- num_or_na(input$range_a)
    sv <- support_from_svgm(sigma2_obs, sigma2_eps, sigma2_mu, sigma2_S, L, range, input$svgm_model)
    if (any(is.na(c(sigma2_obs, sigma2_eps, sigma2_mu, sigma2_S, L, range))) || range <= 0 || L <= 0 || is.na(sv$sigma2_S) || sv$sigma2_S < 0) {
      plot.new(); text(0.5, 0.5, "Enter valid semivariogram inputs", cex = 0.9); return()
    }
    h <- seq(0, max(range, L) * 1.1, length.out = 200)
    rho <- if (input$svgm_model == "spherical") vapply(h, rho_spherical, numeric(1), a = range) else vapply(h, rho_exponential_practical, numeric(1), a = range)
    gamma_latent <- sigma2_mu + sv$sigma2_S * (1 - rho)
    gamma_obs <- sigma2_eps + gamma_latent
    plot(h, gamma_obs, type = "l", lwd = 2, xlab = "lag distance h", ylab = "semivariance", ylim = range(c(0, gamma_obs), finite = TRUE))
    lines(h, gamma_latent, lwd = 2, lty = 2)
    abline(h = sigma2_eps, lty = 3)
    abline(h = sigma2_eps + sigma2_mu, lty = 3)
    abline(v = sv$R, lty = 3)
    legend("bottomright", bty = "n", cex = 0.8,
           legend = c("observed semivariogram", "latent true semivariogram", "σ²_ε", "σ²_ε + σ²_μ", "R = L / √π"),
           lty = c(1, 2, 3, 3, 3), lwd = c(2, 2, 1, 1, 1))
  })
}

shinyApp(ui, server)
