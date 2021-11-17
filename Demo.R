## This script reproduces some components of Figure 4 from the paper. 
## Note: dartMachine is a hacked version of the bartMachine package of 
##       Kapelner and Bleich. I cannot promise that there are no bugs.
##       There are several limitations of the package
##         1. I have not checked that the code is correct when categorical
##            predictors are used.
##         2. The conjugate update from the paper is always used, ignoring
##            the issues mentioned when this update needs modification.
##         3. The prior on alpha is the Beta(0.5, 1) prior on alpha/(alpha+P)
##            This cannot be changed without changing the source code. 
##         4. Technically, the prior on alpha/(alpha+P) is a discretized
##            version of the prior claimed in the paper. This is just to make
##            the update easy, as one does not really need to nail alpha to 
##            get good performance.

options(java.parameters = "-Xmx2g")
library(dartMachine)

Fried <- function(X) {
  10 * sin(pi * X[, 1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
}

FriedSamp <- function(n, dim_x, sigma=1) {
  X <- matrix(runif(n * dim_x), n, dim_x)
  mu <- Fried(X)
  Y <- rnorm(n, mu, sigma)

  return(list(X=X, Y=Y))
}

RMSE <- function(x,y) sqrt(mean((x - y)^2))

set.seed(123)
train <- FriedSamp(250, 100, sqrt(10))
test  <- FriedSamp(1000, 100, sqrt(10))

## Fit bart
bart      <- bartMachine(X = as.data.frame(train$X), y = train$Y)
bart_pred <- predict(bart, as.data.frame(test$X))
rm(bart); gc()

## Fit dart, with alpha / (alpha + P) ~ Beta(0.5, 1)
dart         <- bartMachine(X = as.data.frame(train$X), y = train$Y,
                    do_ard = TRUE,
                    num_trees = 200,
                    num_burn_in = 5000,
                    num_iterations_after_burn_in = 5000,
                    do_prior = TRUE)
dart_pred    <- predict(dart, as.data.frame(test$X))
split_counts <- get_var_counts_over_chain(dart, type = "splits")
s_samples    <- get_cov_prior_select(dart)

rm(dart); gc()

## Fit dart, with alpha = 1
dart_no_prior   <- bartMachine(X = as.data.frame(train$X), y = train$Y,
                             do_ard = TRUE,
                             num_trees = 200,
                             num_burn_in = 5000, 
			     num_iterations_after_burn_in = 5000,
                             alpha_0 = 1)
dart_np_pred    <- predict(dart_no_prior, as.data.frame(test$X))
split_counts_np <- get_var_counts_over_chain(dart_no_prior, type = "splits")
s_samples_np    <- get_cov_prior_select(dart_no_prior)
rm(dart_no_prior); gc()

## RMSE
RMSE(bart_pred, Fried(test$X))
RMSE(dart_pred, Fried(test$X))
RMSE(dart_np_pred, Fried(test$X))

## Reproduction of figure in paper, comparing prior to no-prior
par(mfrow = c(2,2))
plot(colMeans(split_counts > 0), cex=.2)
abline(h = 0.5)
plot(colMeans(s_samples))

plot(colMeans(split_counts_np > 0), cex=.2)
abline(h = 0.5)
plot(colMeans(s_samples_np))
