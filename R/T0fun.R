#' Estimate metabolic heat production T0 using a fitted linear model
#'
#' @param m body mass of the fish in g
#' @param Tm body temperature of the fish
#' @param meso a switch parameter that is either TRUE for mesotherm (or regional
#'   endotherm) or FALSE for ectotherm
#' @param pars the coefficients of the fitted regression model
#' @param omega the coefficient of mass to convert metabolic rate into temperature. Defaults to 10^0.5945
#'
#' @returns a vector the same length as input m containing the corresponding estimated metabolic heat production parameter T0 
#' @export
#'
T0fun <-
function(m, Tm, meso, pars){
  # calculates the rate of heat generation T0 
  # over a range of body temperature values Tm for a 
  # fish of mass m, meso = T or F as a logical toggle for whether to 
  # indicate a regional endotherm (T) or an ectotherm (F), and parameters
  # of the corresponding linear model.
  
  # use the default package supplied value for omega if omitted
  # if(is.null(omega)){omega <- omega}
  
  # Equation 6 in the paper
  return(RMRfun(m, Tm,meso = meso, pars = pars) / (m*omega))
  
}
