# Load Packages -----------------------------------------------------------
library(blendR)
library(survHE)
library(survHEhmc)
library(rstan)
options(mc.cores = parallel::detectCores())
# Example -----------------------------------------------------------------
## data
data("TA174_FCR", package = "blendR")

## externally estimated data
data_sim <- blendR::ext_surv_sim(t_info = 144,
                                 S_info = 0.5,
                                 T_max = 360)
## data fit
obs_Surv <- fit.models(formula = Surv(death_t, death) ~ 1,
                       data = dat_FCR, distr = "exponential",
                       method = "hmc")

ext_Surv <- fit.models(formula = Surv(time, event) ~ 1,
                       data = data_sim, distr = "weibull",
                       method = "hmc")

blend_interv <- list(min = 36, max = 180)
beta_params <- list(alpha = 3, beta = 3)

ble_Surv <- blendR::blendsurv(
  obs_Surv, ext_Surv,
  blend_interv,
  beta_params
)

plot(ble_Surv)

blend_df <- data.frame(
  time = ble_Surv$times,
  surv = ble_Surv$S$S,
  lower = ble_Surv$S$low,
  upper = ble_Surv$S$upp
)

km_fit <- survfit(Surv(death_t, death) ~ 1, data = dat_FCR)
km_df <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper
)

ggplot() +
  # Kaplan-Meier curve (stepwise)
  geom_step(data = km_df, aes(x = time, y = surv), color = "#7CAE00", linetype = "dashed", size = 1) +
  geom_ribbon(data = km_df, aes(x = time, ymin = lower, ymax = upper), alpha = 0.1) +
  # Blended curve (smoothed)
  geom_line(data = blend_df, aes(x = time, y = surv), color = "#F8766D", size = 1.2) +
  geom_ribbon(data = blend_df, aes(x = time, ymin = lower, ymax = upper), alpha = 0.05) +
  # Blending interval as a shaded rectangle
  geom_rect(aes(xmin = blend_interv$min, xmax = blend_interv$max, ymin = 0, ymax = 1),
            fill = "lightblue", alpha = 0.2) +
  labs(title = "Calibration Plot: Blended vs Kaplan–Meier Survival Curve",
       x = "Time", y = "Survival Probability") +
  theme_minimal()

# Log-Likelihood ----------------------------------------------------------
## Observed CDF and PDF Survival functions --------------------------------
obs_Obj <- obs_Surv$models$Exponential # extract observed fit
obs_posterior <- rstan::extract(obs_Obj) # extract posterior data
obs_lambda <- obs_posterior$rate # extract posterior parameters
ext_Obj <- ext_Surv$models$Exponential # extract external fit
ext_posterior <- rstan::extract(ext_Obj) # extract posterior data
ext_lambda <- ext_posterior$rate # extract posterior parameters

blend_loglik <- function(beta_params, time, event,
                         blend_interv, obs_lambda,
                         ext_lambda) {
  alpha <- beta_params[[1]]
  beta  <- beta_params[[2]]
  int_min <- blend_interv[[1]]
  int_max <- blend_interv[[2]]

  M <- length(obs_lambda)  # posterior sample size
  N <- length(time)        # number of observations

  # Normalize time
  u <- (time - int_min) / (int_max - int_min)
  u <- pmin(pmax(u, 1e-8), 1 - 1e-8)  # numerical stability

  # Blending weights
  w <- pbeta(u, alpha, beta)
  w_prime <- dbeta(u, alpha, beta) / (int_max - int_min)

  # Compute matrix of survival and density values: M x N
  S_obs_mat <- exp(-outer(obs_lambda, time))             # exp(-λt)
  f_obs_mat <- sweep(S_obs_mat, 1, obs_lambda, "*")      # λ * exp(-λt)

  S_ext_mat <- exp(-outer(ext_lambda, time))
  f_ext_mat <- sweep(S_ext_mat, 1, ext_lambda, "*")

  # Blend each row of posterior samples (M rows, N columns)
  S_blend_mat <- sweep(S_obs_mat, 2, w, "*") +
    sweep(S_ext_mat, 2, (1 - w), "*")

  f_blend_mat <- sweep(S_obs_mat - S_ext_mat, 2, -w_prime, "*") +
    sweep(f_obs_mat, 2, w, "*") +
    sweep(f_ext_mat, 2, (1 - w), "*")

  # Avoid log(0)
  f_blend_mat <- pmax(f_blend_mat, 1e-15)
  S_blend_mat <- pmax(S_blend_mat, 1e-15)

  # Log-likelihood for each posterior sample
  loglik_mat <- sweep(log(f_blend_mat), 2, event, "*") +
    sweep(log(S_blend_mat), 2, 1 - event, "*")

  # Now average likelihood over posterior samples
  loglik_pointwise <- log(colMeans(exp(loglik_mat)))  # log of mean likelihood

  return(-sum(loglik_pointwise))  # negative log-likelihood for minimization
}

blend_loglik(
  beta_params = list(alpha = 2, beta = 2),
  time = dat_FCR$death_t,
  event = dat_FCR$death,
  blend_interv = list(min = 12, max = 300),
  obs_lambda = obs_posterior$rate,
  ext_lambda = ext_posterior$rate
)

# 2 -----------------------------------------------------------------------
get_posterior_params <- function(model_obj, distr = NULL) {
  models_available <- names(model_obj$models)
  # Auto-detect if distr is not supplied
  if (is.null(distr)) {
    if (length(models_available) != 1) {
      stop("Multiple models found. Please specify the distribution.")
    }
    matched_name <- models_available[[1]]
  } else {
    distr <- tolower(distr)
    match_idx <- grep(distr, tolower(models_available))

    if (length(match_idx) == 0) {
      stop(paste("No model found matching distribution:", distr))
    } else if (length(match_idx) > 1) {
      stop(paste("Multiple models match distribution:", distr,
                 "\nMatches:", paste(models_available[match_idx],
                                     collapse = ", ")))
    }
    matched_name <- models_available[match_idx]
  }

  stanfit <- model_obj$models[[matched_name]]
  post <- rstan::extract(stanfit)

  if (grepl("exponential", tolower(matched_name))) {
    return(list(rate = post$rate))
  } else if (grepl("weibull", tolower(matched_name))) {
    return(list(shape = post$alpha, scale = post$scale))
  } else if (grepl("gompertz", tolower(matched_name))) {
    return(list(shape = post$shape, rate = post$rate))
  } else {
    stop("Unsupported or unrecognized model structure.")
  }
}

get_posterior_params(obs_Surv)
get_posterior_params(ext_Surv)


get_surv_dens_fns <- function(posterior, distr, time) {
  if (distr == "exponential") {
    lambda <- posterior$rate
    S <- exp(-outer(lambda, time, "*"))
    f <- sweep(S, 1, lambda, "*")

  } else if (distr == "weibull") {
    shape <- posterior$shape
    scale <- posterior$scale
    S <- outer(1:length(shape), time, function(i, t) {
      exp(-(t / scale[i])^shape[i])
    })
    f <- outer(1:length(shape), time, function(i, t) {
      (shape[i] / scale[i]) * (t / scale[i])^(shape[i] - 1) *
        exp(-(t / scale[i])^shape[i])
    })

  } else if (distr == "gompertz") {
    shape <- posterior$shape
    rate <- posterior$rate
    S <- outer(1:length(shape), time, function(i, t) {
      exp(-(exp(shape[i] * t) - 1) / shape[i] * rate[i])
    })
    f <- outer(1:length(shape), time, function(i, t) {
      rate[i] * exp(shape[i] * t) *
        exp(-(exp(shape[i] * t) - 1) / shape[i] * rate[i])
    })

  } else {
    stop("Unsupported distribution.")
  }

  list(S = S, f = f)
}

blend_loglik_posterior <- function(par, time, event, blend_interv,
                                   obs_posterior, ext_posterior, distr) {
  alpha <- par[1]
  beta <- par[2]

  int_min <- blend_interv$min
  int_max <- blend_interv$max
  u <- pmin(pmax((time - int_min) / (int_max - int_min), 0), 1)
  w <- pbeta(u, alpha, beta)
  w_prime <- dbeta(u, alpha, beta) / (int_max - int_min)

  obs_vals <- get_surv_dens_fns(obs_posterior, distr, time)
  ext_vals <- get_surv_dens_fns(ext_posterior, distr, time)

  # Take average over all MCMC samples
  S_obs <- colMeans(obs_vals$S)
  f_obs <- colMeans(obs_vals$f)
  S_ext <- colMeans(ext_vals$S)
  f_ext <- colMeans(ext_vals$f)

  # Blended model
  S_blend <- w * S_obs + (1 - w) * S_ext
  f_blend <- -w_prime * (S_obs - S_ext) + w * f_obs + (1 - w) * f_ext

  # Stability
  S_blend <- pmax(S_blend, 1e-15)
  f_blend <- pmax(f_blend, 1e-15)

  loglik <- event * log(f_blend) + (1 - event) * log(S_blend)
  return(-sum(loglik))  # Negative for minimization
}
rstan::extract(obs_Surv$models)
# Example: Using Weibull
distr <- "Exponential"
str(obs_Surv$models$Exponential, max.level = 2)



