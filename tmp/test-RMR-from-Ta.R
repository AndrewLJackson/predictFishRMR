
library(predictFishRMR)


# load the fitted regression model and extract the parameters
data("fittedModelRMR")
pars <- fittedModelRMR$fixef

Ta = 15

meso = FALSE

m = c(100, 1000, 2000, 5000)


# convert Ta to Tm
Tm <- Ta + T0fun(m, Ta, meso = meso, pars = pars) / kfun(m)

RMR_values <- RMRfun(m, Tm, meso, pars)

# Initialize a vector to store the RMR values
# RMR_values <- numeric(length(size_range))



# Calculate RMR for each size in the range
for (i in 1:length(size_range)) {
  m <- size_range[i]
  RMR_values[i] <- RMRfun(m, Tm, meso, pars)
}