## Numerical check for support-variance approximation
## Unit square, structural sill = 1, no nugget.
## Block variance = E{ C(H) }, where H is the distance between two
## independently and uniformly sampled points in the support.

## ------------------------------------------------------------
## Covariance models
## ------------------------------------------------------------

corr_fun <- function(model, practical_range) {
  ## practical_range:
  ## - spherical: compact range
  ## - others: distance where correlation is 0.05
  
  if (model == "spherical") {
    a <- practical_range
    return(function(h) {
      z <- h / a
      ifelse(h >= a, 0, 1 - 1.5 * z + 0.5 * z^3)
    })
  }
  
  if (model == "exponential") {
    a <- practical_range / log(20)   # exp(-h/a) = 0.05 at h = practical_range
    return(function(h) exp(-h / a))
  }
  
  if (model == "matern32") {
    t05 <- uniroot(function(t) (1 + t) * exp(-t) - 0.05,
                   interval = c(1e-8, 50))$root
    a <- practical_range / t05
    return(function(h) {
      t <- h / a
      (1 + t) * exp(-t)
    })
  }
  
  if (model == "matern52") {
    t05 <- uniroot(function(t) (1 + t + t^2 / 3) * exp(-t) - 0.05,
                   interval = c(1e-8, 50))$root
    a <- practical_range / t05
    return(function(h) {
      t <- h / a
      (1 + t + t^2 / 3) * exp(-t)
    })
  }
  
  if (model == "gaussian") {
    a <- practical_range / sqrt(log(20))  # exp(-(h/a)^2) = 0.05
    return(function(h) exp(-(h / a)^2))
  }
  
  stop("Unknown model: ", model)
}

## ------------------------------------------------------------
## Distance-density integration
## ------------------------------------------------------------

square_block_var <- function(corr) {
  ## Difference density for two independent uniform points in unit square:
  ## f(dx, dy) = (1 - |dx|)(1 - |dy|), dx,dy in [-1,1]^2.
  ## Use symmetry: 4 * integral_0^1 integral_0^1 ...
  
  inner <- function(x) {
    integrate(
      f = function(y) corr(sqrt(x^2 + y^2)) * (1 - x) * (1 - y),
      lower = 0,
      upper = 1,
      rel.tol = 1e-9,
      abs.tol = 1e-11
    )$value
  }
  
  4 * integrate(
    f = Vectorize(inner),
    lower = 0,
    upper = 1,
    rel.tol = 1e-8,
    abs.tol = 1e-10
  )$value
}

disk_distance_density <- function(h, R) {
  ## Density of the distance between two independent uniform points
  ## in a disk of radius R.
  out <- numeric(length(h))
  ok <- h >= 0 & h <= 2 * R
  hh <- h[ok]
  
  out[ok] <-
    (4 * hh / (pi * R^2)) * acos(hh / (2 * R)) -
    (2 * hh^2 / (pi * R^3)) * sqrt(pmax(0, 1 - hh^2 / (4 * R^2)))
  
  out
}

disk_block_var <- function(corr, R) {
  integrate(
    f = function(h) corr(h) * disk_distance_density(h, R),
    lower = 0,
    upper = 2 * R,
    rel.tol = 1e-9,
    abs.tol = 1e-11,
    subdivisions = 1000
  )$value
}

linear_disk_block_var <- function(corr, R) {
  ## Manuscript approximation:
  ## Var_B approx 1 - E(H)/R * {1 - rho(R)}
  ## with E(H)/R = 128/(45*pi) for a disk.
  rho_R <- corr(R)
  1 - 128 / (45 * pi) * (1 - rho_R)
}

## ------------------------------------------------------------
## Run comparison
## ------------------------------------------------------------

models <- c("spherical", "exponential", "matern32", "matern52", "gaussian")
model_labels <- c(
  spherical   = "Spherical",
  exponential = "Exponential",
  matern32    = "Matern nu=3/2",
  matern52    = "Matern nu=5/2",
  gaussian    = "Gaussian"
)

range_ratios <- c(1, 2, 10, 50)

L <- 1
R <- L / sqrt(pi)  # equal-area disk radius for unit square

res <- do.call(
  rbind,
  lapply(range_ratios, function(rr) {
    do.call(
      rbind,
      lapply(models, function(m) {
        cf <- corr_fun(m, practical_range = rr)
        
        exact_square <- square_block_var(cf)
        exact_disk   <- disk_block_var(cf, R)
        linear_disk  <- linear_disk_block_var(cf, R)
        
        data.frame(
          range_block_side = rr,
          model = model_labels[[m]],
          exact_square = exact_square,
          exact_disk = exact_disk,
          linear_disk = linear_disk,
          err_disk_vs_square_pct =
            100 * (exact_disk / exact_square - 1),
          abs_err_disk_vs_square_pct =
            100 * abs(exact_disk / exact_square - 1),
          err_linear_vs_square_pct =
            100 * (linear_disk / exact_square - 1),
          abs_err_linear_vs_square_pct =
            100 * abs(linear_disk / exact_square - 1)
        )
      })
    )
  })
)

## Full table: more transparent than the current compact manuscript table
print(
  within(res, {
    exact_square <- round(exact_square, 4)
    exact_disk <- round(exact_disk, 4)
    linear_disk <- round(linear_disk, 4)
    err_disk_vs_square_pct <- round(err_disk_vs_square_pct, 2)
    err_linear_vs_square_pct <- round(err_linear_vs_square_pct, 2)
  })[, c(
    "range_block_side", "model",
    "exact_square", "exact_disk", "linear_disk",
    "err_disk_vs_square_pct", "err_linear_vs_square_pct"
  )],
  row.names = FALSE
)

## Compact manuscript-style worst-case summary
summary_tab <- aggregate(
  cbind(abs_err_disk_vs_square_pct, abs_err_linear_vs_square_pct) ~ range_block_side,
  data = res,
  FUN = max
)

summary_tab$abs_err_disk_vs_square_pct <-
  round(summary_tab$abs_err_disk_vs_square_pct, 1)
summary_tab$abs_err_linear_vs_square_pct <-
  round(summary_tab$abs_err_linear_vs_square_pct, 1)

names(summary_tab) <- c(
  "Range / block side",
  "Max abs. error: exact disk vs exact square (%)",
  "Max abs. error: linear disk vs exact square (%)"
)

print(summary_tab, row.names = FALSE)
