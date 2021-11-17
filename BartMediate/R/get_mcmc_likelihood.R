get_mcmc_likelihood <- function(fit_y, Y, fit_m, M, design_y, design_m) {
  N <- nrow(design_y)

  ## M weights ----

  m_hats <- bart_machine_get_posterior(fit_m, design_m)[["y_hat_posterior_samples"]]
  
  if(fit_m$pred_type == "regression") {
    rez <- t(M - m_hats)
    sigma_sq_m <- get_sigsqs(fit_m)
    log_weight_m <- -0.5 * rez^2 / sigma_sq_m - 0.5 * log(2 * pi * sigma_sq_m)  
  } else if(fit_m$pred_type == "classification") {
    m_levels <- levels(M)
    log_weight_m <- t((M == m_levels[2]) * log(m_hats) + (M == m_levels[1]) * log(1 - m_hats))
  }

  ## Y weights ----

  y_hats <- bart_machine_get_posterior(fit_y, design_y)[["y_hat_posterior_samples"]]
  
  if(fit_y$pred_type == "regression") {
    rez <- t(Y - y_hats)
    sigma_sq_y <- get_sigsqs(fit_y)
    log_weight_y <- -0.5 * rez^2 / sigma_sq_y - 0.5 * log(2 * pi * sigma_sq_y)  
  } else if(fit_m$pred_type == "classification") {
    y_levels <- levels(Y)
    log_weight_y <- t((Y == y_levels[2]) * log(y_hats) + (M == y_levels[1]) * log(1 - y_hats))
  }
  
  log_weight <- log_weight_m + log_weight_y
  
  return(log_weight)

}