avg_mediation_effects <- function(Y_00, Y_01, Y_10, Y_11, omega) {
  out <- c(sum(omega * (Y_01 - Y_00)), sum(omega * (Y_11 - Y_10)), 
           sum(omega * (Y_10 - Y_00)), sum(omega * (Y_11 - Y_01)),
           sum(omega * (Y_11 - Y_00)))
  
  return(out)
}