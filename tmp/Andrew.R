
# order of the parameters is:
# 1 = alpha the coefficient of log10(mass)
# 2 = beta the coefficient of Body temperature
# 3 = Intercept (gamma)
# 4 = the effect of mesothermy (psi)
# 5 = omega, the exponent of the multiplier of mass for converting RMR to T0
pars <- c(0.8161589, 0.05, 0.98, 0.68, 5.945)


# ------------------------------------------------------------------------------
# apply our already fitted linear model of how RMR scales
# with temperature, mass and physiology.
# ------------------------------------------------------------------------------
RMRfun <- function(m, Tm, meso, pars) {
  # calculates RMR over a range of body temperature values Tm for a 
  # fish of mass m, meso = T or F as a logical toggle for whether to 
  # indicate a regional endotherm (T) or an ectotherm (F), and parameters
  # of the corresponding linear model.
  
  # determines the coefficints for 
  ifelse(meso, b_meso <- pars[4],  b_meso <- 0)

  
  RMR <- 10 ^ (pars[1]*log10(m) + pars[2]*Tm + pars[3] + b_meso)
  
  return(RMR)
  
  
}


# ------------------------------------------------------------------------------
# estimate metabolic heat production T0 using our linear model
# ------------------------------------------------------------------------------
T0fun <- function(m, Tm, meso, pars){
  # calculates the rate of heat generation T0 
  # over a range of body temperature values Tm for a 
  # fish of mass m, meso = T or F as a logical toggle for whether to 
  # indicate a regional endotherm (T) or an ectotherm (F), and parameters
  # of the corresponding linear model.
  
  return(RMRfun(m, Tm,meso = meso, pars = pars) / (m*10^pars[5]))
  
}

# ------------------------------------------------------------------------------
## Estimate k from Nakamura scaling
# ------------------------------------------------------------------------------
kfun <- function(m, a = 0.00161, b = -0.62){
  # Uses the reported allometric scaling of rate of cooling 
  # with body mass. Defaults to Nakamura values. 
  return( a * m^(b) )
}
  
# ------------------------------------------------------------------------------
## calculate DeltaT = (Tm - Ta)
# ------------------------------------------------------------------------------
Tafun <- function(m, Tm, meso, pars){
  # calculates the ambient water temperature associated with heat balance
  # over a range of body temperature values Tm for a 
  # fish of mass m, meso = T or F as a logical toggle for whether to 
  # indicate a regional endotherm (T) or an ectotherm (F), and parameters
  # of the corresponding linear model.
  return(Tm - T0fun(m, Tm, meso = meso, pars = pars) / kfun(m))
  
}


# ------------------------------------------------------------------------------
# Solve the model for body temperature as a function of ambient water temp
# ------------------------------------------------------------------------------
bodyTempLambert <- function(x, m, meso, pars){
  
  # calculate the cooling rate k from allometric relationship with mass
  kk <- kfun(m)
  
  # extract the parameters to pass onwards
  aa <- pars[1] # alpha, the coefficient of log10(mass) in RMR
  bb <- pars[2] # beta, the coefficient of body temperature
  gg <- pars[3] # gamma, the intercept
  ifelse(meso, pp <- pars[4],  pp <- 0) # psi, the effect of mesothermy
  om <- pars[5] # omega, the exponent of the multiplier of mass to convert RMR to T0
  
  # evaluate the two arms of the Lambert function that gives Real solutions
  y0 <- x - 
    (lambertW0( - ((10^(bb*x+gg+pp-om)) * (m^(aa-1)) * bb * log(10)) / kk )) / (bb* log(10))
  y1 <- x - 
    (lambertWm1(- ((10^(bb*x+gg+pp-om)) * (m^(aa-1)) * bb * log(10)) / kk )) / (bb* log(10))
  
  
  # lmbd <- gg + pp -om
  
  # x_crit <- (0.43429 * (log((0.43429 * kk * m^(1 - aa)) / (bb)) - 
  #                         (2.3026 * (gg + pp -om)) - 1 )) / bb
  # 
  # return the data
  return(data.frame(x = c(x, NA, x), y = c(y0, NA, y1), 
                    branch = c(rep(0, length(y0)), NA,rep(-1, length(y1)))))
}

# ------------------------------------------------------------------------------
# Calculate the maximum value of water temperature that yields real numbers
# ------------------------------------------------------------------------------
maxAmbientTemp <- function(m, meso, pars){
  
  # calculate the cooling rate k from allometric relationship with mass
  kk <- kfun(m)
  
  # extract the parameters to pass onwards
  aa <- pars[1] # alpha, the coefficient of log10(mass) in RMR
  bb <- pars[2] # beta, the coefficient of body temperature
  gg <- pars[3] # gamma, the intercept
  ifelse(meso, pp <- pars[4],  pp <- 0) # psi, the effect of mesothermy
  om <- pars[5] # omega, the exponent of the multiplier of mass to convert RMR to T0
  
  # this approximation involves complex numbers owing to the branching of the 
  # Lambert function, but only the real component is meaningful for the two 
  # branches that concern us.
  x_max <- (0.43429 * (log((0.43429 * kk * m^(1 - aa)) / (bb)) - 
                          (2.3026 * (gg + pp -om)) - 1 )) / bb
  
  return(x_max)
  
}

# ------------------------------------------------------------------------------
# Generate a sequence of Ta values up to the critical value
# NB NOT IMPLEMENTED PROPERLY YET
# ------------------------------------------------------------------------------
seqPlausibleTa <- function(lwr, upp, m, meso, pars){
  # Returns a sequence of Ta values suitable for estimating and plotting a 
  # smooth curve.
  # Values that fall outside the maximum possible Ta will be truncated.
  # lower (lwr) and upper (upp) desired values are requested but are 
  # truncated to the maximum possible. The sequence is scaled logarithmically
  # to ensure more fine scale gradient close to the maximum.
  
  # calculate the maximum ambient water temperature
  x_max <- maxAmbientTemp(m, meso, pars)
  
  # generate the sequence
  a <- lwrv
  b <- min(x_max, upp)
  
  # create a logarithmic scale
  x_log <- abs(a) * (exp(seq(0,1, length.out = 10)) ) / (exp(1) + abs(b))
  
}

# ------------------------------------------------------------------------------
# Get RMR using fixed Ta
# ------------------------------------------------------------------------------
calcRMRForSizeRange <- function(size_range, Ta = Ta, meso, pars) {
  # This function calculates the RMR for a range of fish sizes given a fixed ambient temperature (Ta)
  # size_range: a vector of fish sizes
  # Ta: fixed ambient temperature (default is 20 degrees)
  # meso: boolean indicating whether the fish is mesothermic (TRUE) or ectothermic (FALSE)
  # pars: parameters of the linear model
  
  # Initialize a vector to store the RMR values
  RMR_values <- numeric(length(size_range))
  
  # Calculate Tm from Ta
  Tm <- Ta + T0fun(m, Ta, meso = meso, pars = pars) / kfun(m)
  
  # Calculate RMR for each size in the range
  for (i in 1:length(size_range)) {
    m <- size_range[i]
    RMR_values[i] <- RMRfun(m, Tm, meso, pars)
  }
  
  return(RMR_values)
}


# ------------------------------------------------------------------------------
# Get Tm above Ta across size
# ------------------------------------------------------------------------------
calcElevationTmAboveTa <- function(size_range, Ta, meso, pars) {
  # This function calculates the elevation of Tm above Ta for a range of fish sizes given Ta
  # size_range: a vector of fish sizes
  # Ta: ambient temperature
  # meso: boolean indicating whether the fish is mesothermic (TRUE) or ectothermic (FALSE)
  # pars: parameters of the linear model
  
  # Initialize a vector to store the elevation of Tm above Ta
  elevation_values <- numeric(length(size_range))
  
  # Calculate the elevation of Tm above Ta for each size in the range
  for (i in 1:length(size_range)) {
    m <- size_range[i]
    Tm <- Ta + T0fun(m, Ta, meso = meso, pars = pars) / kfun(m)
    elevation_values[i] <- Tm - Ta
  }
  
  return(elevation_values)
}



