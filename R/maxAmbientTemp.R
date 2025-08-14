#' Calculate the maximum value of water temperature that yields real numbers
#'
#' @param m body mass of the fish in g
#' @param meso a switch parameter that is either TRUE for mesotherm (or regional
#'   endotherm) or FALSE for ectotherm
#' @param pars the coefficients of the fitted regression model
#' @param preK is a prefactor that is applied to the estimated cooling
#'   coefficient K that allows uncertainty in this parameter to be explored
#' @param kpars a vector contain "a" the prefactor of the exponential scaling
#'   equation and "b" the scaling exponent. If k is provided, kpars is ignored.
#' @param omega the coefficient of mass to convert metabolic rate into temperature. Defaults to 10^0.5945
#'
#' @returns a vector of length m with corresponding max water temperatures in C
#' @export
#' 
maxAmbientTemp <-
function(m, meso, pars, preK = 1, kpars = NULL){
  
  # calculate the cooling rate k from allometric relationship with mass
  # using provided scaling parameters if supplied. 
  # if(is.null(kpars))  kk <- kfun(m) * preK
  # if(!is.null(kpars)) kk <- kfun(m, kpars) * preK
  ifelse(is.null(kpars),
         kk <- kfun(m) * preK,
         kk <- kfun(m, kpars) * preK)
  
  # kk <- kfun(m) * preK
  
  
  # extract the parameters to pass onwards
  gg <- pars[1] # gamma, the intercept
  aa <- pars[2] # alpha, the coefficient of log10(mass) in RMR
  bb <- pars[3] # beta, the coefficient of body temperature
  ifelse(meso, pp <- pars[4],  pp <- 0) # psi, the effect of mesothermy
  
  # omega is an internal data object loaded and available to files in R/
  # use the default package supplied value for omega if omitted
  # if(is.null(omega)){omega <- omega}
  om <- omega
  
  # this approximation involves complex numbers owing to the branching of the 
  # Lambert function, but only the real component is meaningful for the two 
  # branches that concern us.
  # x_max <- (0.434294481903252 * 
  #             (log((0.434294481903252 * kk * m^(1 - aa)) / (bb)) - 
  #                         (2.30258509299405 * (gg + pp -omega)) - 1 )) / bb
  x_max <- (-gg - pp + log(kk*om*m^(1-aa) /  bb) - 1) / bb
  
  return(x_max)
  
}
