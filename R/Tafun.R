#' Calculates the ambient water temperature associated with heat balance for a given body temperature
#'
#' @param m body mass in g
#' @param Tm body temperature in C
#' @param meso a switch parameter that is either TRUE for mesotherm (or regional
#'   endotherm) or FALSE for ectotherm
#' @param pars the coefficients of the fitted regression model
#' @param kpars a vector contain "a" the prefactor of the exponential scaling
#'   equation and "b" the scaling exponent.
#'
#' @returns the ambient water temperature in C
#' @export
#'
Tafun <-
function(m, Tm, meso, pars, kpars = NULL){
  # calculates the ambient water temperature associated with heat balance
  # over a range of body temperature values Tm for a 
  # fish of mass m, meso = T or F as a logical toggle for whether to 
  # indicate a regional endotherm (T) or an ectotherm (F), and parameters
  # of the corresponding linear model.
  return(Tm - T0fun(m, Tm, meso = meso, pars = pars) / kfun(m, kpars = kpars))
  
}
