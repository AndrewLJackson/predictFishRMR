#' Calculate the maximum value of water temperature that yields real numbers
#'
#' @param m body mass of the fish in g
#' @param meso a switch parameter that is either TRUE for mesotherm (or regional
#'   endotherm) or FALSE for ectotherm
#' @param pars the coefficients of the fitted regression model
#' @param preK is a prefactor that is applied to the estimated cooling
#'   coefficient K that allows uncertainty in this parameter to be explored
#'
#' @returns a vector of length m with corresponding max water temperatures in C
#' @export
#' 
maxAmbientTemp <-
function(m, meso, pars, preK = 1){
  
  # calculate the cooling rate k from allometric relationship with mass
  kk <- kfun(m) * preK
  
  # extract the parameters to pass onwards
  aa <- pars[1] # alpha, the coefficient of log10(mass) in RMR
  bb <- pars[2] # beta, the coefficient of body temperature
  gg <- pars[3] # gamma, the intercept
  ifelse(meso, pp <- pars[4],  pp <- 0) # psi, the effect of mesothermy
  om <- pars[5] # omega, the exponent of the multiplier of mass to convert RMR to T0
  
  # this approximation involves complex numbers owing to the branching of the 
  # Lambert function, but only the real component is meaningful for the two 
  # branches that concern us.
  x_max <- (0.434294481903252 * 
              (log((0.434294481903252 * kk * m^(1 - aa)) / (bb)) - 
                          (2.30258509299405 * (gg + pp -om)) - 1 )) / bb
  
  return(x_max)
  
}
