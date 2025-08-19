#' Esimate Routine Metabolic Rate (RMR) from ambient water temperature (Ta)
#'
#' @param m mass of fish in kg
#' @param Ta temperature of the ambient water
#' @param meso a logical value where TRUE indicates mesotherm and FALSE
#'   indicates ectotherm
#' @param pars pars a vector of length 4 containing the coefficients of the
#'   corresponding fitted linear regression
#' @param kpars a vector contain "a" the prefactor of the exponential scaling
#'   equation and "b" the scaling exponent. Supplying NULL, the default, uses
#'   default parameter. See `?kfun` for details.
#'
#' @returns a vector the same length as `m` with corresponding RMR values
#' @export
#'
#' @examples RMRfromTa(m = c(100, 1000), Ta = 15, meso = T, pars = c(2.95, 0.83, 0.08, 1.23), kpars = c(a = 0.00201, b = -0.617))
RMRfromTa <- function(m, Ta, meso, pars, kpars = NULL) {
  
  # convert Ta to Tm
  Tm <- Ta + T0fun(m, Ta, meso = meso, pars = pars) / kfun(m, kpars)
  
  # calculate RMR
  RMR_values <- RMRfun(m, Tm, meso, pars)
  
  
  return(RMR_values)
}
