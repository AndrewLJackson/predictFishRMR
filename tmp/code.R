


# High fuel demands and overheating risk challenge mesothermic fishes in warming oceans

# Nicholas L. Payne*, Edward P. Snelling, Ignacio Peralta-Maraver, Dave E. Cade, Taylor K. Chapple, Alex G. McInturf, Yuuki Y. Watanabe, David W. Sims, 
# Nuno Queiroz, Ivo da Costa, Lara L. Sousa, Jeremy A. Goldbogen, Haley R. Dolton & Andrew L. Jackson

# *email: paynen@tcd.ie


# --- libraries
library(phytools) # phylogenetic analysis
library(nlme) # pgls analysis
library(geiger)
library(plyr) # organise dataset
library(yarrr) # transparent color
library(viridis) # color palettes
library(brms) # fit baysian models
library(parallel) # run models in parallel
library(rstan) # run bayes
library(bayestestR)#Bayesian p-values
library(geiger) # covariance matrix
library(raster) # plot maps
library(sp)
library(maps) # plot maps





# SECTION 1: HOUSEKEEPING

# --- open data
rmr_data <- read.csv("MandT.csv")

# this code is used to fit models using Tc data
#rmr_data2 <- read.xlsx("T_core.xlsx")
#rmr_data[1:10,]$T <- rmr_data2$Tc
#rmr_data[1:10,]$logRMR <- rmr_data2$logRMR

rmr_data$spp <- gsub("_", " ", rmr_data$spp) #remove underscore

# some corrections are needed
rmr_data$spp[rmr_data$spp == "Scomberesocidae sp."] <- "Scomberesocidae" # correct this species
rmr_data$spp[rmr_data$spp == "Ambassis sp"] <- "Ambassis sp." # correct this species
rmr_data$spp[rmr_data$spp == "Hypoatherina sp"] <- "Hypoatherina sp." # correct this species
rmr_data$spp[rmr_data$spp == "Monacanthidae sp"] <- "Monacanthidae" # correct this species
rmr_data$spp[rmr_data$spp == "Monacanthidae sp."] <- "Monacanthidae" # correct this species
rmr_data$spp[rmr_data$spp == "Carcharodon carcharius"] <- "Carcharodon carcharias" # correct this species
rmr_data$spp[rmr_data$spp == "Stegostoma tigrinum"] <- "Stegostoma fasciatum"
rmr_data <- subset(rmr_data, spp != "Caulophrynidae sp.") # this family is not included in the phylo tree
rmr_data <- rmr_data[order(rmr_data$spp),]

sp_names <- read.csv("sp_list.csv", sep = ";")
sp_names$spp <- gsub("_", " ", sp_names$spp) #remove underscore
sp_names$original.name <- gsub("_", " ", sp_names$original.name) #remove underscore
sp_names$sp_all <- gsub("_", " ", sp_names$sp_all) #remove underscore
sp_names$original.name[sp_names$original.name == "Carcharodon carcharius"]<- "Carcharodon carcharias" # correct this species
sp_names$original.name[sp_names$original.name == "Stegostoma tigrinum"]<- "Stegostoma fasciatum" # correct this species

#sp_names[sp_names$original.name == "Carcharhinus acronotus",]
tree <- read.tree("shark_fish.trees")
tree$tip.label <- gsub("_", " ", tree$tip.label) #remove underscore
tree$tip.label <- gsub("-", ",", tree$tip.label) #remove underscore
tree$tip.label[tree$tip.label == "Pseudocaranx dentex ,formaly Longirostrum delicatissimus"] <- "Pseudocaranx dentex (formaly Longirostrum delicatissimus) "
tree$tip.label[tree$tip.label == "Sardinops sagax ,formaly S. caerulea"] <- "Sardinops sagax (formaly S. caerulea)"

# Info for Stegostoma_tigrinum comes from Stegostoma_fasciatum
sp_tree <- unique(tree$tip.label)
sp <- unique(rmr_data$spp)

# check point: diff must be 0
setdiff(rmr_data$spp, tree$tip.label)
setdiff(rmr_data$spp, sp_names$original.name)

# --- housekeeping
rmr_data$method <- factor(rmr_data$method)
rmr_data$ram <- factor(rmr_data$ram)
rmr_data$therm <- factor(rmr_data$therm)
rmr_data <- rmr_data[complete.cases(rmr_data),]

#rownames(rmr_data) <- rmr_data$spp
#spp <- rownames(rmr_data) 

# --- check for potential cofounding effects between 

library(MatchIt)

# calculate propensity scores considering method as the treatment and logM, T, therm as covariates:

# fit a propensity score model using nearest neighbor matching
ps_model <- matchit(method ~ logM, data = rmr_data, method = "nearest")
summary(ps_model)

# extract the matched dataset and fit the model again
matched_data <- match.data(ps_model)
mod_matched <- lm(logRMR ~ logM + T + therm + method, data = matched_data)
summary(mod_matched) # method does not have sig. effect.


# SECOND SECTION: BAYESIAN LINEAR MIXED MODEL 

# construct covariance matrix of species (Hadfield & Nakagawa, 2010).
A <- ape::vcv.phylo(compute.brlen(tree,method="Grafen"))

# adding species as a factor into our dataset
rmr_data$spp_name<-factor(rmr_data$spp)

# Many species include serveral masurements. Thus, we use a bayesian models.
# General STAN specifications 
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

		
		
		# MODEL
		# We use a bayesian approach to study the relationship
		# across species (phylogeny) while taking into considration
		# intraspaecific variation and multiple observations for
		# several species.



# set model priors
priors =c(
  prior(normal(0,10), "b"),
  prior(normal(0,50), "Intercept"),
  prior(student_t(3,0,20), "sd"),
  prior(student_t(3,0,20), "sigma")
)



# run Bayesian model including phylogenetic tree
mod <- brm((logRMR)  ~ logM + T + therm # fixed coefficients
								+ (1 | gr(spp, cov = A)) + (1| spp_name), # random factors
  								family = gaussian(), 
  								data2 = list(A = A),
  								prior = priors,
  								sample_prior = TRUE, 
  								chains = 2, 
  								cores = 8, # cores to be used 8 in mac pro
  								iter = 3000, 
  								warmup = 500,
  								control = list(adapt_delta = 0.99, max_treedepth = 15), # Improve convergence
  								save_pars = save_pars(all = TRUE),
  								data = rmr_data)
  								
# check model convergence 
plot(mod)
pp_check(mod)

# extract posteriors
posterior_samples <- posterior_samples(mod)

# get coefficients
summary(mod) 

# get bayesian p-value
p_map(mod)

# ---
# assign color for the two states

map_data <- ddply(rmr_data, c("spp", "therm"), summarise, logM = mean(logM))[,c(1,2)]
row.names(map_data) <- map_data$spp

#Vector with information of the caracter
a <-  tree$tip.label[order(tree$tip.label)]
b <-  map_data[map_data$spp %in% tree$tip.label, ]$therm

info_caracter <- setNames(c(b), a)
ancestor <- make.simmap(tree, info_caracter,model="ER",nsim=50)
ecomorph<-as.factor(getStates(ancestor,"tips"))

cols<-setNames(c("blue","red"), levels(ecomorph))

# coef mod
summary(mod)
coef <- c(0.98, 0.81, 0.05, 0.68)


# --- Figure 1. 
# A persistently higher energetic cost to fuel mesothermy in fishes. Left panel shows the phylogeny used in the regression of mass, 
# body temperature and thermal strategy (red for mesotherms and blue for ectotherms) on fish routine metabolic rate, shown in right panel. 
# Data with black outline indicate RMR estimated by our new heat production method, with all other data derived from respirometry. 


	quartz("Fig 1", 9,4.2)
	par(mfrow = c(1,2))

	plot(ancestor[[1]],colors=sapply(cols,make.transparent,alpha=0.1),type = "fan",add=F, lwd=4, ftype="off", fsize=0.6,offset=0.5)
	for(i in 1:length(ancestor))
	plot(ancestor[[i]],colors=sapply(cols,make.transparent,alpha=0.1),type = "fan",add=T,lwd=4,ftype="off",fsize=0.6,offset=0.5)
	par(mar = c(5.1, 5.1, 2.1,2.1))
	rmr_data$therm <- factor(rmr_data$therm)
	pt.cols <- setNames(c("blue", "red"), levels(rmr_data$therm))

	plot((logRMR-T*coef[3]) ~ logM, data = rmr_data, pch = c(21,24)[as.numeric(rmr_data$method)], cex = c(0.9,1.3)[as.numeric(rmr_data$method)], bg = pt.cols[rmr_data$therm], col = "white", lwd = 0.5, ylab = "", xlab = "", xaxt = "n", yaxt = "n")
	
	mtext(side = 2, expression("Routine metabolic rate (g O"[2]*" h"^-1*" )"), line = 3.5)

	axis(1, at = log10(c(0.00001, 0.001, 0.1,10, 1000)), labels = c("10 mg", "1 g", "100 g", "10 Kg", "1Ton"))
	axis(2,las = 1, at = log10(c(0.0001,0.001, 0.01,0.1, 1, 10, 100, 1000, 10000)), labels = c(expression("10"^-7*""), expression("10"^-6*""), expression("10"^-5*""), expression("10"^-4*""), 	expression("10"^-3*""), 
	expression("10"^-2*""), expression("10"^-1*""), expression("1"), expression("10")))    
	
	legend("bottomright", c("Respirometry", "Heat production", "", "Mesotherms", "Ectotherms"), pch = c(1,2, NA, 21,21), col = c("black", "black", NA, "red", "blue"), pt.bg = c(NA,NA,NA,"red",  "blue"), bty = "n", pt.cex = 1.2, cex = 0.8)

	# adding endotherms fit
	x = seq(-1,3.2, length = 100)
	y = coef[1]+coef[4] + x* coef[2] 
	lines(x,y, lwd = 4, col = "white"); lines(x,y, lwd = 2, col = "red")

	# adding ectotherm fit
	x = seq(-5,3, length = 100)
	y = coef[1] + x* coef[2] 
	lines(x,y, lwd = 4, col = "white"); lines(x,y, lwd = 2, col = "blue")
	mtext(side = 1, "Body weight", line = 3)
	
	points((logRMR-T*coef[3]) ~ logM, data = rmr_data, pch = c(NA,24)[as.numeric(rmr_data$method)], cex = c(1,1.3)[as.numeric(rmr_data$method)], bg = pt.cols[rmr_data$therm], col = "white")

box()


# --- Figure 2
# Scaling of energy demand and heat balance limits in large-bodied fishes. (A) Routine metabolic rate RMR is 7-8 fold higher for mesotherms than ectotherms,
# and fluctuates widely for the mesotherms as body temperature varies. Plotted are RMR estimates for fish swimming in water ranging from 20-25°C, 
# with body temperature estimated based on equation 10. (B) The empirically derived gigantothermy phenomenon in fish shows how bigger
#  fish – particularly mesotherms – have elevated body temperature due to a scaling mismatch between heat production and loss. 
# Shown is estimated Ta when Tm (body temperature) is 20°C. (C & D) Theoretical heat balance thresholds (Ta) arise from an inability 
# of larger fish to balance heat production and loss in warm water. As fish mass increases, there is a threshold water temperature in which 
# heat produced by metabolism will be balanced by heat loss to the water (y axis); above those Ta the fish will continue to heat unless
#  they reduce heat production (e.g. swim slower), increase heat loss by adjusting convective cooling, or both. Larger fish have lower threshold temperatures,
# especially mesotherms given higher rates of heat production, theoretically restricting them to cooler climates unless they engage physiological
# or behavioural adjustments. Dashed line represents xy unity.

source("Andrew.R")

n = 100
m = size_range = seq(1,2000, length = n)

## A meso - expected Ta and RMR
meso <- T  # indicating whether it is mesothermic (TRUE) or ectothermic (FALSE)
RMR_meso_20 <- calcRMRForSizeRange(size_range, Ta = 20, meso, pars)
RMR_meso_25 <- calcRMRForSizeRange(size_range, Ta = 25, meso, pars)

## A endo - expected Ta and RMR
meso <- F  
RMR_ecto_20 <- calcRMRForSizeRange(size_range, Ta = 20, meso, pars)
RMR_ecto_25 <- calcRMRForSizeRange(size_range, Ta = 25, meso, pars)

# --- plot

# Pannel A
quartz("Fig2", 9,7)

layout(matrix(c(1,3,0,0,2,4,0,5,0,0), 2), widths = c(1,0.4,1,0.4, 0.2))
	par(oma = c(2,5,2,1), mar = c(3.5,0,3.5,0))

	m <-  size_range
	plot(m, RMR_meso_20/1000, type = "n", ylim = c(0,450), xlim = c(0,2000), las = 1, ylab = "", xlab = "", col = "red", xaxs = "i",yaxs = "i", xaxt = "n")
	axis(1, at = c(1,500,1000,2000), labels = c("1 Kg", "500 Kg", "1Ton", "2 Ton"))
	
	mtext(side = 2, line = 3, expression("Routine metabolic rate (g O"[2]*" h"^-1*" )"))
	mtext(side = 1, line = -25, expression("Body mass"), outer = T, adj = 0.39)
	
	polygon(c(m, rev(m)), c(RMR_meso_20/1000, rev(RMR_meso_25/1000)), border = F, col = transparent("black", trans.val = 0.4))
	lines(m, RMR_meso_20/1000,  col = "black", lty = 2)
	lines(m, RMR_meso_25/1000,  col = "black")

	polygon(c(m, rev(m)), c(RMR_ecto_20/1000, rev(RMR_ecto_25/1000)), border = F, col = transparent("grey", trans.val = 0.5))
	lines(m, RMR_ecto_20/1000,  col = "black", lty = 2)
	lines(m, RMR_ecto_25/1000,  col = "black")

		text(2800, 310, expression("20"*degree*"C"))
		text(2800, 660, expression("25"*degree*"C"))

		legend("topleft", c("Mesotherms", "Ectotherms"), pch = 21, col = c("black", "black"), pt.bg = c("grey20",  "grey"), bty = "n", pt.cex = 2, cex = 1)


# Panel B
Ta = 20

		meso <- TRUE  # indicating the fish is ectothermic (FALSE)
		elevation_values <- calcElevationTmAboveTa(size_range, Ta, meso, pars)

		plot(size_range ,elevation_values, type = "l", col = "black", yaxs = "i", las = 1, xlim = c(0,2000), ylim = c(0,9.5), xaxs = "i", ylab = "", xaxt = "n", xlab = "", lwd = 1.5)
axis(1, at = c(1,500,1000,2000), labels = c("1 Kg", "500 Kg", "1Ton", "2 Ton"))
	
		meso <- F  # indicating the fish is ectothermic (FALSE)
		elevation_values <- calcElevationTmAboveTa(size_range, Ta, meso, pars)

		lines(size_range ,elevation_values, type = "l", col = "grey", lwd = 1.5)
		mtext(side = 2, line = 3, expression("Elevation of T"["m"]*" above T"["a"]*" ("*degree*"C)"))


# Panel C 
 # Initialize a vector to store the Ta values
 Tm <- seq(1, 35, by = 0.5)
 Ta <- numeric(length(Tm))
  m = seq(1,2000, length = 2500)


  # Calculate Ta for each Tm value
  meso = F
   col <-(plasma(length(m)))
 
  plot(Tm, Tm, xlim = c(1,35), ylim = c(1,36), type = "n", xaxs = "i",  yaxs = "i", xlab = "", ylab = "", las = 1)
 for(i in 1:length(m)){lines(Tm, Tafun(m[i], Tm, meso, pars), col =col[i])}
abline(0,1, lty = 2)

mtext(side = 2, expression("Heat balanced T"["a"]*" ("*degree*" C)") , line = 3)
box()


  # Calculate Ta for each Tm value
  meso = T
  plot(Tm, Tm, xlim = c(1,35), ylim = c(1,36), type = "n", xaxs = "i",  yaxs = "i", xlab = "", ylab = "", las = 1)
for(i in 1:length(m)){lines(Tm, Tafun(m[i], Tm, meso, pars), col =col[i])}
abline(0,1, lty = 2)

mtext(side = 1, expression("Body temperature T"["m"]*" ("*degree*" C)"), line = 0, outer = T, adj = 0.4)

polygon(c(c(20.5,36), c(36,20.5)), c(c(9.5,24), c(-100, -100)), col = "white", border = F)
box()

par(mar = c(5,1,5,6))

library(plotrix)

plot(NA, NA, ylim = c(1,2000), xlim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "", xaxs = "i")
axis(4, at = c(1, 500,1000,2000), labels = c("1 Kg", "500 Kg", "1Ton", "2 Ton"), las = 1)

gradient.rect(0,-85,1, 2085, col=plasma(length(m)),border=NA,gradient="y")
box()
mtext(side = 4, "Body mass", line = 4)

segments(21,10,35,20)


# Figure 3. 
# Energy landscape and theoretical heat balance thresholds for current and future large-bodied fishes. 
# Routine metabolic rate increases with water temperature and body size, and is higher for mesothermic fishes. 
# In all panels RMR is estimated for current and future Sea Surface Temperature scenarios for a model 1 ton ectotherm 
# and 300 kg mesotherm via equations 10 and 5, with white regions indicating SSTs that exceed the theoretical thresholds 
# beyond which those fish cannot balance heat production and loss unless engaging behavioural 
# (e.g. diving to cooler depths as seen in salmon sharks; 33) or physiological mechanisms 
# (such as swimming more slowly or altering blood circulation patterns, as seen in bluefin tuna; 39). 
# These heat-balance dynamics could help explain why large fishes migrate seasonally (winter vs summer panels); 
# why contemporary mesotherms are distributed toward higher latitudes than ectotherms (top vs bottom panels); 
# and could help predict range shifts for large fishes under future ocean warming (“Today” vs “2080-2100” SST scenarios). 
#Scenario SSP4-6.0 was used to generate future warming scenarios, with winter and summer reflecting minimum and maximum annual SST estimates.


	quartz("Fig3", 24,8)
	layout(matrix(c(1,2,0,0,3,4,0,0,5,6,0,0,7,8,9,9),2), widths = c(1, 0.08,1, 0.15, 1, 0.08, 1, 0.4))
	par(oma = c(4,6,5,0), mar = c(2,0,1,0), cex.axis = 2)

#panel A
size <- m <- 1000
meso = F

	world_ta1 <- raster("thetao_baseline_2010_2019_depthsurf_min.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta1) <- ifelse(values(world_ta1)>x_max, NA, values(world_ta1) )
	Ta <- world_ta1
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n", yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim =c(1500, 89469.57), col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	mtext(line = 1, "Today", cex = 1.5)
	axis(2, las = 1)
	axis(1, labels = F)
	box()
	
#panel B
size <- m <- 300
meso = T

	world_ta1 <- raster("thetao_baseline_2010_2019_depthsurf_min.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta1) <- ifelse(values(world_ta1)>x_max, NA, values(world_ta1) )
	Ta <- world_ta1
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n", yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim =c(1500, 89469.57), col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	axis(2, las = 1)
	axis(1, labels = T)
	box()
	
#panel C
size <- m <- 1000
meso = F
	world_ta2 <- raster("thetao_ssp460_2090_2100_depthsurf_min.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta2) <- ifelse(values(world_ta2)>x_max, NA, values(world_ta2) )
	Ta <- world_ta2
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n",yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim = c(1500, 89469.57),col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	mtext(line = 1, "2080-2100", cex = 1.5)
	axis(1, labels = F)
    axis(2, labels = F)
	
#panel D
size <- m <- 300
meso = T
	world_ta2 <- raster("thetao_ssp460_2090_2100_depthsurf_min.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta2) <- ifelse(values(world_ta2)>x_max, NA, values(world_ta2) )
	Ta <- world_ta2
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n",yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim = c(1500, 89469.57),col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	axis(1)
    axis(2, labels = F)
		
mtext("Winter", line = 2, outer = T, adj = 0.21, cex = 2)
	
#panel E
size <- m <- 1000
meso = F

	world_ta1 <- raster("thetao_baseline_2010_2019_depthsurf_max.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta1) <- ifelse(values(world_ta1)>x_max, NA, values(world_ta1) )
	Ta <- world_ta1
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n", yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim =c(1500, 89469.57), col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	mtext(line = 1, "Today", cex = 1.5)
	axis(1, labels = F)
    axis(2, labels = F)
	
#panel F
size <- m <- 300
meso = T

	world_ta1 <- raster("thetao_baseline_2010_2019_depthsurf_max.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta1) <- ifelse(values(world_ta1)>x_max, NA, values(world_ta1) )
	Ta <- world_ta1
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, yaxt = "n", xaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim =c(1500, 89469.57), col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	axis(1)
    axis(2, labels = F)
	
#panel G
size <- m <- 1000
meso = F
	world_ta2 <- raster("thetao_ssp585_2090_2100_depthsurf_max.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta2) <- ifelse(values(world_ta2)>x_max, NA, values(world_ta2) )
	Ta <- world_ta2
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n",yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim = c(1500, 89469.57),col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	mtext(line = 1, "2080-2100", cex = 1.5)
	axis(1, labels = F)
    axis(2, labels = F)

#panel H
size <- m <- 300
meso = T
	world_ta2 <- raster("thetao_ssp585_2090_2100_depthsurf_max.nc")
	x_max <- maxAmbientTemp(m, meso, pars)
	values(world_ta2) <- ifelse(values(world_ta2)>x_max, NA, values(world_ta2) )
	Ta <- world_ta2
    Tm <- Ta + T0fun(size, Ta, meso, pars) / kfun(size)
    RMR_values <- RMRfun(m, Tm, meso, pars)
	plot(NA,NA,  xlim = c(-160,160), ylim = c(-75,75), las = 1, xaxt = "n",yaxt = "n", ylab = "", xlab = "")
	image((RMR_values), zlim = c(1500, 89469.57),col = magma(50), add = T)
	map("world", fill=TRUE, col="grey", add=TRUE, lwd=0.5)
	axis(1)
    axis(2, labels = F)
mtext("Summer", line = 2, outer = T, adj = 0.71, cex = 2)	
	
par(mar = c(8,1,8,6))
library(plotrix)
plot(NA, NA, ylim = c(1500, 89469.57), xlim = c(0,1), xaxt = "n", yaxt = "n", ylab = "", xlab = "", xaxs = "i",yaxs = "i")
axis(4, las = 1, at = c(10,20,30,40,50,60,70,80,90)*1000, labels = c(10,20,30,40,50,60,70, 80, 90))

gradient.rect(0,1500, 1, 89469.57, col= magma(50),border=NA, gradient="y")
box()
mtext(side = 4, expression("Routine metabolic rate (g O"[2]*" h"^-1*" )"), line = 3.5, cex = 1.5)

mtext(side = 1, outer = T, expression("Longitude ("*degree*")"), line = 3, adj = 0.46, cex = 1.5)
mtext(side = 2, outer = T, expression("Latitude ("*degree*")"), line = 3, adj = 0.475, cex = 1.5)

