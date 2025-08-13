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
#' @format A list object obcontaining the fixed effects estimates derived from
#'   the brms model fit of RMR by mass, temperature and species. The entire
#'   model object is not provided as it is tens of Mb in size and is not
#'   strictly required.
#'
#' @keywords datasets
#' @author Andrew Jackson
"fittedModelRMR"

