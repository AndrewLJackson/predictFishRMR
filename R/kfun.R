#' Estimate k from Nakamura's exponential scaling
#'
#' @param m body mass of the fish in g
#' @param kpars a vector contain "a" the prefactor of the exponential scaling
#'   equation and "b" the scaling exponent
#'
#' @returns an estimate of the cooling coefficient k
#' @export
#'
#' @examples
#' kfun(1000)
#' kfun(m = c(1, 10, 100, 1000), a = 0.0018, b = -0.63)
kfun <-
function(m, kpars = c(a = 0.0018, b = -0.63)){
  # Uses the reported allometric scaling of rate of cooling 
  # with body mass. Defaults to Nakamura values. 
  return( kpars[1] * m^(kpars[2]) )
}

