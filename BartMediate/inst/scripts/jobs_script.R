## Load packages ----

library(mediation)
library(gridExtra)
library(tidyverse)
options(java.parameters = "-Xmx5g")
library(BartMediate)

set.seed(digest::digest2int("jobs_analysis"))

## Load data ----

jobs <- mediation::jobs

jobs_y <- jobs %>%
  select(treat, sex, age, econ_hard, 
         marital, nonwhite, educ, income, depress1, job_seek) %>%
  mutate(educ = as.numeric(educ), income = as.numeric(income))

jobs_m <- jobs_y %>% select(-job_seek)

## Fit mediator models ----

jobs_m_0 <- jobs_m %>% filter(treat == 0) %>% select(-treat)
jobs_m_1 <- jobs_m %>% filter(treat == 1) %>% select(-treat)
jobs_m_a <- jobs_m %>% select(-treat)

fit_m_0 <- bartMachine(X = jobs_m_0, y = jobs$job_seek[jobs$treat == 0],
                       num_burn_in = 5000,
                       num_iterations_after_burn_in = 10000,
                       do_ard = TRUE, do_prior = TRUE, 
                       seed = sample.int(1E6, 1))

fit_m_1 <- bartMachine(X = jobs_m_1, y = jobs$job_seek[jobs$treat == 1],
                       num_burn_in = 5000,
                       num_iterations_after_burn_in = 10000,
                       do_ard = TRUE, do_prior = TRUE, 
                       seed = sample.int(1E6, 1))

## Get clever covariates ----

jobs_y_ric <- jobs_y %>% 
  mutate(M_hat_0 = predict(fit_m_0, jobs_m_a), 
         M_hat_1 = predict(fit_m_0, jobs_m_a))

## Fit outcome models ----

jobs_y_0 <- jobs_y_ric %>% filter(treat == 0) %>% select(-treat)
jobs_y_1 <- jobs_y_ric %>% filter(treat == 1) %>% select(-treat)

fit_y_0    <- bartMachine(X = jobs_y_0, y = jobs$depress2[jobs$treat == 0],
                          num_burn_in = 5000,
                          num_iterations_after_burn_in = 10000,
                          do_ard = TRUE, do_prior = TRUE,
                          seed = sample.int(1E6, 1))

fit_y_1    <- bartMachine(X = jobs_y_1, y = jobs$depress2[jobs$treat == 1],
                          num_burn_in = 5000,
                          num_iterations_after_burn_in = 10000,
                          do_ard = TRUE, do_prior = TRUE,
                          seed = sample.int(1E6, 1))

## Doing the mediation ----

mediated_bart_ric <- 
  mediate_bart_stratify_sensitivity(fit_y_0, fit_y_1, fit_m_0, fit_m_1,
                                 design_y = jobs_y_ric %>% select(-treat),
                                 design_m = jobs_m_a,
                                 trt_name = "treat", med_name = "job_seek",
                                 iters = seq(from = 1, to = 10000, by = 10),
                                 num_copy = 2)

wide_bart_ric <- mediated_bart_ric %>% filter(CopyID == 1) %>%
  pivot_wider(names_from = "ParamName", values_from = "Param") %>%
  select(-Iteration, -CopyID)

print(AGCCombine(mediated_bart_ric))

## Setting up for sensitivity ----

shift_effect <- function(lambda, df) {
  df[,"delta_0"] <- df[,"delta_0"] - lambda * df[,"theta_a"]
  df[,"delta_1"] <- df[,"delta_1"] - lambda * df[,"theta_a"]
  df[,"zeta_0"]  <- df[,"zeta_0"] + lambda * df[,"theta_a"]
  df[,"zeta_1"]  <- df[,"zeta_1"] + lambda * df[,"theta_a"]
  df[,"lambda"]  <- lambda
  return(df)
}

lambda_grid <- seq(from = -0.5, to = 0.5, by = 0.01)
lambda_bart <- do.call(rbind, lapply(lambda_grid, shift_effect, df = wide_bart_ric))

## Making SA Figure ----

sa_bart <- lambda_bart %>% mutate(iter = rep(1:1000, length(lambda_grid)))
sa_bart <- sa_bart %>% pivot_longer(cols = c("delta_0", "delta_1", "zeta_0", "zeta_1", "tau", "theta_a"), 
                            names_to = "ParName", values_to = "Param")

sa_bart_sum <- sa_bart %>% group_by(lambda, ParName) %>% summarise(
  mean = mean(Param),
  lcl  = quantile(Param, c(0.025)), 
  ucl = quantile(Param, 0.975)
)

p_1 <- 
  sa_bart_sum %>% filter(ParName == "delta_0") %>% 
  ggplot(aes(x = lambda, y = mean)) + geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = scales::muted("gray", 70, 70)) + 
  geom_line(size = 2, linetype = "dashed") + theme_bw() + xlab("$\\lambda$") + ylab("$\\delta(0)$") + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + ylim(-0.25, 0.25)

p_2 <- 
  sa_bart_sum %>% filter(ParName == "zeta_1") %>% 
  ggplot(aes(x = lambda, y = mean)) + geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = scales::muted("gray", 70, 70)) + 
  geom_line(size = 2, linetype = "dashed") + theme_bw() + xlab("$\\lambda$") + ylab("$\\zeta(1)$") + 
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0)  + ylim(-0.25, 0.25)

grid.arrange(p_1, p_2)

## Treatment-as-covariate Approach ----

fit_m <- bartMachine(X = jobs_m, y = jobs$job_seek,
                     num_burn_in = 5000,
                     num_iterations_after_burn_in = 10000,
                     do_ard = TRUE, do_prior = TRUE, 
                     seed = sample.int(1E6, 1))

## Get clever covariates ----

jobs_y_ric_2 <- jobs_y %>%
  mutate(M_hat_0 = predict(fit_m, jobs_m),
         M_hat_1 = predict(fit_m, jobs_m))

## Fit outcome models ----

fit_y    <- bartMachine(X = jobs_y_ric_2, y = jobs$depress2,
                        num_burn_in = 5000,
                        num_iterations_after_burn_in = 10000,
                        do_ard = TRUE, do_prior = TRUE,
                        seed = sample.int(1E6, 1))

## Doing the mediation ----

mediated_bart_ric_cc <-
  mediate_bart_sensitivity(fit_y, fit_m,
                           design_y = jobs_y_ric_2,
                           design_m = jobs_m,
                           trt_name = "treat", med_name = "job_seek",
                           iters = seq(from = 1, to = 10000, by = 10),
                           num_copy = 2)

## Tidy Bayes Figure ----

library(tidybayes)

mediated_bart_2 <- rbind(mediated_bart_ric %>% mutate(Method = "Stratified"), 
                         mediated_bart_ric_cc %>% mutate(Method = "CC"))
mediated_bart_2$Variable <- factor(
  paste0("$\\", mediated_bart_2$ParamName, "$"), 
  levels = rev(c("$\\delta_0$", "$\\delta_1$", "$\\zeta_0$", "$\\zeta_1$", 
                 "$\\tau$", "$\\theta_a$")))

mediated_bart_2 %>%
ggplot(aes(y = Variable, x = Param)) + stat_halfeyeh() + theme_bw() + 
xlab("Value") + facet_wrap(~Method, nrow = 2)

## One More Figure ----

mediated_bart_cc_wide <- pivot_wider(mediated_bart_ric_cc, 
                                     names_from = "ParamName", 
                                     values_from = "Param") %>% 
filter(CopyID == 1) %>%
select(-Iteration, -CopyID) %>% 
mutate(Method = "CC", Delta = delta_1 - delta_0)

wide_bart_ric_2 <- wide_bart_ric %>% mutate(Method = "Stratified", Delta = delta_1 - delta_0)

wide_bart_comb <- rbind(mediated_bart_cc_wide, wide_bart_ric_2)

ggplot(wide_bart_comb, aes(x = Delta, fill = Method)) + geom_histogram(bins = 50, color = 'white') + 
facet_wrap(~Method, scales = "free_y") + theme_bw() + theme(legend.position = "none") + 
xlab("$\\delta(1) - \\delta(0)$") + ylab("Frequency")
