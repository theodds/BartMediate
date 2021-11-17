log_sum_exp <- function(x) {
  M <- max(x)
  M + log(sum(exp(x - M)))
}
