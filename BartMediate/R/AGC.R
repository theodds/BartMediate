#' Inference using accelerated g-computation
#'
#' Combines the results of analysis on multiple simulated datasets. Associated
#' to each simulated dataset is an estimate of the parameters of interest.
#' Datasets are indexed by (i) the iteration of the bootstrap/MCMC scheme used
#' to generate the simulated data and (ii) an id indicating which of K
#' simulations the dataset corresponds to for a particular iteration.
#'
#' @param results A data.frame with columns Iteration giving the iteration of
#' the sample parameter value, CopyID giving the copy number,and Param giving the
#' estimate of the desired causal parameter, and ParamName giving the name of
#' the parameter.
#'
#' @param level The confidence level to be used for constructing confidence
#' intervals and performing hypothesis tests.
#'
#' @return A data.frame object with the following columns
#' \itemize{
#'   \item ParamName - The name of the parameter
#'   \item estimate  - The aggregated estimate of the parameter
#'   \item standard_error - The estimated standard error of the estimate
#'   \item between_var - The "between" component of the variance
#'   \item within_var - The "within" component of the variance
#'   \item var_est  - The estimated variance (square of standard_error)
#'   \item dof - The degrees-of-freedom of the Sattherthwaite approximation
#'   \item Z - Test statistic for testing if the parameter is zero
#'   \item p_value - The (two-sided) p-value for the test that the parameter is zero.
#' }
AGCCombine <-function(results, level = 0.95) {
  ParamName <- NULL; Iteration <- NULL; Param <- NULL

  estimate <- mean(results$Param)
  M        <- length(unique(results$Iteration))
  K        <- length(unique(results$CopyID))
  unique_params <- unique(results$ParamName)

  out <- data.frame(ParamName = character(0), estimate = numeric(0),
                    standard_error = numeric(0),
                    between_var = numeric(0),
                    within_var = numeric(0),
                    var_est = numeric(0),
                    dof = numeric(0),
                    Z = numeric(0),
                    p_value = numeric(0))

  for(name in unique_params) {
    result_sub <- results %>% filter(ParamName == name)
    within_var <- result_sub %>%
      group_by(Iteration) %>%
      summarise(within = var(Param) / K) %>%
      pluck("within") %>% mean()
    between_var <- result_sub %>%
      group_by(Iteration) %>%
      summarise(within_est = mean(Param)) %>%
      pluck("within_est") %>% var()
    between_var <- between_var
    estimate <- mean(result_sub$Param)
    var_est <- (1 + 1/M) * between_var - within_var

    dof <- var_est ^ 2 / (
      (1 + 1/M)^2 * between_var^2 / (M - 1) +
        within_var^2 / M / (K - 1)
    )

    Z <- estimate / sqrt(var_est)
    p <- 2 * pt(abs(Z), df = dof, lower.tail = FALSE)
    q <- 1 - (1 - level) / 2
    Lower <- estimate - qt(p = q, df = dof) * sqrt(var_est)
    Upper <- estimate + qt(p = q, df = dof) * sqrt(var_est)

    out <- rbind(out, data.frame(
      ParName = name,
      estimate = estimate, standard_error = sqrt(var_est),
      between_var = between_var, within_var = within_var, var_est = var_est,
      dof = dof, Z = Z, p_value = p,
      Lower = Lower, Upper = Upper
    ))
  }

  return(out)
}
