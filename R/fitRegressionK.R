#' Fit a log-log regression model to estimate K from body mass
#' 
#' Fits a random intercept model of cooling coefficient by body mass with 
#' species as random factor
#'
#' @param dd a data.frame or tibble object comprising three columns: "species", 
#' "M"(log mass in ) and cooling coefficient "K". 
#' Defaults to the included dataset Kmass.
#'
#' @returns a model object of a log-log regression
#' @export
#' 
fitRegressionK <- 
function(dd = NULL){
  
  if(is.null(dd)) {utils::data("Kmass"); dd <- Kmass}
  
  out <- nlme::lme(K ~ M, 
                   random = ~ 1 | species, 
                   data = dd, method="ML")
  
  return(out)
}