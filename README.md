
# predictFishRMR

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/1008301182.svg)](https://doi.org/10.5281/zenodo.16409669)
<!-- badges: end -->

The goal of predictFishRMR is to estimate the Routine Metabolic Rate of free swimming fish based on allometric relationships with mass and their thermal strategy: either ectothermy or mesothermy. 


## Installation


### Preferred installation route

You can install the development version of predictFishRMR from [GitHub](https://github.com/AndrewLJackson/predictFishRMR) after installing the package `pak`:

``` r
# install.packages("pak")
pak::pak("AndrewLJackson/predictFishRMR", build_vignettes = TRUE)
```

or you can install a particular version

``` r
# install.packages("pak")
pak::pak("AndrewLJackson/predictFishRMR@v0.1.2-beta", build_vignettes = TRUE)
```

Alternatively you can install this package easily by cloning the repository to your local computer or downloading the zip file and extracting it to your local working directory. Then open the `predictFishRMR.Rproj` to load Rstudio and point it to this package. Then you can install it along with the vingettes using:

``` r
# install.packages("devtools")
devtools::install(build_vignettes = TRUE)
```

## Notes for Nacho and Nick

+ The file `vignettes/fit-RMR-regression.Rmd` contains the code to run the phylogenetic regression estimating RMR from mass and thermal strategy. This file needs to be updated with the latest dataset and re-run. The last chunk is not evaluated by default, but we need to then save the fitted model object and put it in the package folder `data/`. This will then allow the other functions to look up this model fit and pull the coefficients for the rest of the functions. `vignette("fit-RMR-regression", package = "predictFishRMR")`
+ The file `vignettes/regression-K-from_mass.Rmd` similarly runs the regression estimating K from mass and produces a figure with the confidence intervals for this regression. The model is fitted by the function `fitRegression()` using the raw data save as a dataset loaded by `data(Kmass)`. There is no need to update this model fitting file or data unless something new has arisen. `vignette("regression-K-from-mass", package = "predictFishRMR")`
+ The file `vignettes/predict-Ta-crit-sensitivity.Rmd` generates some figures that we discussed adding to the paper including a sensitivity analysis of how the estimated Ta that permits thermal equilibrium varies with mass and thermal strategy for double and half K estimates.`vignette("predict-Ta-crit-sensitivity", package = "predictFishRMR")`
+ The file `vignettes/worked-example.Rmd` illustrates how you can calculate RMR for a fish of given size and thermal strategy at a range of water temperatures. `vignette("worked-example", package = "predictFishRMR")`
+ When it comes to resubmitting the paper, we can create a new release and give it a stable link to install that particular version of the package along with a doi. We currently have a doi for the version `v0.1.1-beta` [GitHub](https://github.com/AndrewLJackson/predictFishRMR)


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(predictFishRMR)
## basic example code TBC
```

