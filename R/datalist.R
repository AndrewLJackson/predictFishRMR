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
#' @format A list containing two matrix arrays. $fixef contains a summary of the
#'   fixed effects and $ranef the random effects by species. These have been
#'   obtained by using the brms pacakge to extract these from a fitted model.
#'   The actual brms model is several tens of Mbs and so we elected to only
#'   provide the estimated coefficients. In subsequent analyses, only the means
#'   of the fixed effect coefficients are used.
#'
#' @keywords datasets
#' @author Andrew Jackson
"fittedModelRMR"

