#' Calculate RMR for a fish of given mass, body temperature and thermal strategy
#'
#' @param m body mass in g
#' @param Tm body temperature in C
#' @param meso a logical value where TRUE indicates mesotherm and FALSE indicates ectotherm
#' @param pars a vector of length 4 containing the coefficients of the corresponding fitted linear regression
#'
#' @returns a vector the same length as input m containing the corresponding estimates of resting metabolic rate RMR
#' @export
#'
RMRfun <-
function(m, Tm, meso, pars) {
  
  # determines the coefficients for 
  ifelse(meso, b_meso <- pars[4],  b_meso <- 0)

  
  RMR <- 10 ^ (pars[1]*log10(m) + pars[2]*Tm + pars[3] + b_meso)
  
  return(RMR)
  
  
}
