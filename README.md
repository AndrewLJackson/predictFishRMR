
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

or you can install a particular version such as the latest stable release

``` r
# install.packages("pak")
pak::pak("AndrewLJackson/predictFishRMR@v1.0.0", build_vignettes = TRUE)
```

### Alternative installation route

You can install this package easily by cloning the repository to your local computer or downloading the zip file and extracting it to your local working directory. Then open the `predictFishRMR.Rproj` to load Rstudio and point it to this package. Then you can install it along with the vingettes using:

``` r
# install.packages("devtools")
devtools::install(build_vignettes = TRUE)
```

## Notes for co-authors

+ The file `vignettes/predict-Ta-crit-sensitivity.Rmd` generates some figures that we discussed adding to the paper including a sensitivity analysis of how the estimated Ta that permits thermal equilibrium varies with mass and thermal strategy for double and half K estimates.`vignette("predict-Ta-crit-sensitivity", package = "predictFishRMR")`
+ The file `vignettes/worked-example.Rmd` illustrates how you can calculate RMR for a fish of given size and thermal strategy at a range of water temperatures. `vignette("worked-example", package = "predictFishRMR")`



## Example

This is a basic example which shows you how to solve a common problem: to be implemented. 

``` r
library(predictFishRMR)
## basic example code TBC
```

