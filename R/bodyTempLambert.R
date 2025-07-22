#' Solve the model for body temperature as a function of ambient water temp
#'
#' @param x ambient water temperature in C
#' @param m body mass in g
#' @param k cooling coefficient k. Defaults to NULL which causes it to be
#'   estimated from m via allometric scaling kfun()
#' @param meso a switch parameter that is either TRUE for mesotherm (or regional
#'   endotherm) or FALSE for ectotherm
#' @param pars the coefficients of the fitted regression model
#'
#' @returns body temperature in C
#' @export
#'
bodyTempLambert <-
function(x, m, k = NULL,meso, pars){
  
  # calculate the cooling rate k from allometric relationship with mass if not
  # provided
  kk <- ifelse(is.null(k), kk <- kfun(m), k)
  
  # extract the parameters to pass onwards
  aa <- pars[1] # alpha, the coefficient of log10(mass) in RMR
  bb <- pars[2] # beta, the coefficient of body temperature
  gg <- pars[3] # gamma, the intercept
  ifelse(meso, pp <- pars[4],  pp <- 0) # psi, the effect of mesothermy
  om <- pars[5] # omega, the exponent of the multiplier of mass to convert RMR to T0
  
  # evaluate the two arms of the Lambert function that gives Real solutions
  y0 <- x - 
    (lamW::lambertW0( - ((10^(bb*x+gg+pp-om)) * (m^(aa-1)) * bb * log(10)) / kk )) / (bb* log(10))
  y1 <- x - 
    (lamW::lambertWm1(- ((10^(bb*x+gg+pp-om)) * (m^(aa-1)) * bb * log(10)) / kk )) / (bb* log(10))
  
  
  # lmbd <- gg + pp -om
  
  # x_crit <- (0.43429 * (log((0.43429 * kk * m^(1 - aa)) / (bb)) - 
  #                         (2.3026 * (gg + pp -om)) - 1 )) / bb
  # 
  # return the data
  return(data.frame(x = c(x, NA, x), y = c(y0, NA, y1), 
                    branch = c(rep(0, length(y0)), NA,rep(-1, length(y1)))))
}
