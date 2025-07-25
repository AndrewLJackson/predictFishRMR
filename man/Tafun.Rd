% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tafun.R
\name{Tafun}
\alias{Tafun}
\title{Calculates the ambient water temperature associated with heat balance for a given body temperature}
\usage{
Tafun(m, Tm, meso, pars)
}
\arguments{
\item{m}{body mass in g}

\item{Tm}{body temperature in C}

\item{meso}{a switch parameter that is either TRUE for mesotherm (or regional
endotherm) or FALSE for ectotherm}

\item{pars}{the coefficients of the fitted regression model}
}
\value{
the ambient water temperature in C
}
\description{
Calculates the ambient water temperature associated with heat balance for a given body temperature
}
