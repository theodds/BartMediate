mediate_bart_sensitivity <- function(fit_y, fit_m, design_y, design_m, 
                                     trt_name, med_name, 
                                     med_pred_name = NULL,
                                     max_iter = NULL, 
                                     iters = NULL,
                                     num_copy = 2,
                                     effect = avg_mediation_effects) {
  
  N <- nrow(design_y)
  pred_type_y <- fit_y$pred_type
  pred_type_m <- fit_m$pred_type
  # stopifnot(pred_type_y %in% c("regression", "classification"))
  # stopifnot(pred_type_m %in% c("regression", "classification"))
  stopifnot(pred_type_y == "regression")
  stopifnot(pred_type_m == "regression")
  
  stopifnot(!is.null(max_iter) | !is.null(iters))
  stopifnot(is.null(max_iter) | is.null(iters))
  if(!is.null(max_iter)) iters <- 1:max_iter
  n_iter <- length(iters)

  design_m_0 <- design_m; design_m_0[[trt_name]] <- 0
  design_m_1 <- design_m; design_m_1[[trt_name]] <- 1
  design_y_00 <- design_y; design_y_00[[trt_name]] <- 0
  design_y_01 <- design_y; design_y_01[[trt_name]] <- 0
  design_y_10 <- design_y; design_y_10[[trt_name]] <- 1
  design_y_11 <- design_y; design_y_11[[trt_name]] <- 1

  sigma_m <- NULL
  sigma_y <- NULL
  if(pred_type_y == "regression") sigma_y <- sqrt(get_sigsqs(fit_y))
  if(pred_type_m == "regression") sigma_m <- sqrt(get_sigsqs(fit_m))
  M_0 <- M_1 <- Y_00 <- Y_01 <- Y_10 <- Y_11 <- NULL


  out_list <- list()

  for(j in 1:n_iter) {
    i <- iters[j]

    omega   <- as.numeric(MCMCpack::rdirichlet(1, rep(1,N)))
    mu_m_0  <- predict_iteration(fit_m, design_m_0, i)
    mu_m_1  <- predict_iteration(fit_m, design_m_1, i)
    
    theta_a <- sum(omega * (mu_m_1 - mu_m_0))

    for(k in 1:num_copy) {

      epsilon <- rnorm(N)
      if(pred_type_m == "regression") {
        M_0 <- mu_m_0 + sigma_m[i] * epsilon
        M_1 <- mu_m_1 + sigma_m[i] * epsilon
      } else if(pred_type_m == "classification") {
        ## NOTE: For classification, mu_m the probability of being at level 2
        m_levels <- levels(design_y[[med_name]])
        M_0 <- ifelse(qnorm(mu_m_0) + epsilon > 0, m_levels[2], m_levels[1])
        M_1 <- ifelse(qnorm(mu_m_1) + epsilon > 0, m_levels[2], m_levels[1])
      } else stop()

      design_y_00[[med_name]] <- M_0
      design_y_01[[med_name]] <- M_1
      design_y_10[[med_name]] <- M_0
      design_y_11[[med_name]] <- M_1

      epsilon_y <- rnorm(N)
      mu_y_00   <- predict_iteration(fit_y, design_y_00, i)
      mu_y_01   <- predict_iteration(fit_y, design_y_01, i)
      mu_y_10   <- predict_iteration(fit_y, design_y_10, i)
      mu_y_11   <- predict_iteration(fit_y, design_y_11, i)
      if(pred_type_y == "regression") {
        Y_00      <- mu_y_00 + sigma_y[i] * epsilon_y
        Y_01      <- mu_y_01 + sigma_y[i] * epsilon_y
        Y_10      <- mu_y_10 + sigma_y[i] * epsilon_y
        Y_11      <- mu_y_11 + sigma_y[i] * epsilon_y
      } else if(pred_type_y == "classification") {
        Y_00      <- ifelse(qnorm(mu_y_00) + epsilon_y > 0, 1, 0)
        Y_01      <- ifelse(qnorm(mu_y_01) + epsilon_y > 0, 1, 0)
        Y_10      <- ifelse(qnorm(mu_y_10) + epsilon_y > 0, 1, 0)
        Y_11      <- ifelse(qnorm(mu_y_11) + epsilon_y > 0, 1, 0)
      } else stop()

      out_list[[length(out_list) + 1]] <-
        data.frame(
          Iteration = i,
          CopyID    = k,
          Param     = c(effect(Y_00, Y_01, Y_10, Y_11, omega), theta_a),
          ParamName = c("delta_0", "delta_1", "zeta_0", "zeta_1", "tau", "theta_a")
        )

    }

    if(j %% 10 == 0) cat(paste("Finishing iteration", j, "\n"))

  }

  return(do.call(rbind, out_list))
  
}