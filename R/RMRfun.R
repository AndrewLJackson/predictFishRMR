#' Calculate RMR for a fish of given mass, body temperature and thermal strategy
#'
#' @param m body mass in g, typically provided as a vector but could be scalar.
#' @param Tm body temperature in C
#' @param meso a logical value where TRUE indicates mesotherm and FALSE indicates ectotherm
#' @param pars a vector of length 4 containing the coefficients of the corresponding fitted linear regression
#'
#' @returns a vector the same length as input m containing the corresponding estimates of resting metabolic rate RMR
#' @export
#'
RMRfun <-
function(m, Tm, meso, pars) {
  
  # determines the coefficient associate with mesothermy or ectothermy
  ifelse(meso, b_meso <- pars[4],  b_meso <- 0)

  # calculate RMR based on the regression parameters
  RMR <- exp(pars[1] + pars[2]*log(m) + pars[3]*Tm + b_meso)
  
  # return output
  return(RMR)
  
  
}
