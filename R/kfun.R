#' Estimate k from Nakamura's exponential scaling
#'
#' @param m body mass of the fish in g
#' @param a the prefactor of the exponential scaling equation
#' @param b the scaling exponent 
#'
#' @returns an estimate of the cooling coefficient k
#' @export
#'
#' @examples
#' kfun(1000)
#' kfun(m = c(1, 10, 100, 1000), a = 0.00161, b = -0.62)
kfun <-
function(m, a = 0.00161, b = -0.62){
  # Uses the reported allometric scaling of rate of cooling 
  # with body mass. Defaults to Nakamura values. 
  return( a * m^(b) )
}
