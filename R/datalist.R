#' Cooling coefficients by body mass and species.
#'
#' Taken from Nakamuru XYZH REF!
#'
#' @docType data
#'
#' @usage data(Kmass)
#'
#' @format An object of class `"data.frame"` containing three variables. 
#' The first column is a character vector of species names labelled "species". 
#' The second is labelled "M" is a numeric vector of body mass values in
#' log10(tons). The third is labelled "K" represents cooling coefficients. 
#'
#' @keywords datasets
#' @author Andrew Jackson
"Kmass"


#' Fitted RMR model.
#'
#'
#'
#' @docType data
#'
#' @usage data(fittedModelRMR)
#'
#' @format An object of class "`brms`" containing the full model fit for the RMR
#'   regression by mass, temperature and species. This is a large file but is
#'   used in a vignette to detail the analysis of model convergence and
#'   predictive power.
#'
#' @keywords datasets
#' @author Andrew Jackson
"fittedModelRMR"

