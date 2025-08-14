# High fuel demands and overheating risk challenge mesothermic fishes in warming oceans

# Nicholas L. Payne*, Edward P. Snelling, Ignacio Peralta-Maraver, Dave E. Cade, Taylor K. Chapple, 
# Alex G. McInturf, Yuuki Y. Watanabe, David W. Sims, Nuno Queiroz, Ivo da Costa, Lara L. Sousa, 
# Jeremy A. Goldbogen, Haley R. Dolton & Andrew L. Jackson

# *email: paynen@tcd.ie

# --- libraries
library(phytools)        # pylogenetic analysis
library(nlme)            # pgls analysis
library(geiger)
library(plyr)            # organise dataset
library(yarrr)           # transparent color
library(viridis)         # color palettes
library(brms)            # fit baysian models
library(parallel)        # run models in parallel
library(rstan)           # run bayes
library(bayestestR)      # Bayesian p-values
library(geiger)          # covariance matrix
library(raster)          # plot maps
library(sp)
library(maps)            # plot maps

# SECTION 1: HOUSEKEEPING

# --- open data
rmr_data <- read.csv("tmp/nacho_RMR/MandT_indiv_whalefixed_raw.csv")
tree <- read.tree("tmp/nacho_RMR/fish_shark_combined.tre")

# Compare species between dataset and phylogeny
unique(rmr_data$spp)
tree$tip.label

tree$tip.label <- gsub("_", " ", tree$tip.label)

tree$tip.label[tree$tip.label == "Pseudocaranx dentex -formaly Longirostrum delicatissimus"] <- "Pseudocaranx dentex (formaly Longirostrum delicatissimus) "
tree$tip.label[tree$tip.label == "Ranzania laevis- Ostracion boops stage"] <- "Ranzania laevis, Ostracion boops stage"
tree$tip.label[tree$tip.label == "Coregonus lavaretus -larvae"] <- "Coregonus lavaretus (larvae)"
tree$tip.label[tree$tip.label == "Sardinops sagax -formaly S. caerulea"] <- "Sardinops sagax (formaly S. caerulea)"
tree$tip.label[tree$tip.label == "Ophichthus gomesii- leptocephalus larvae"] <- "Ophichthus gomesii, leptocephalus larvae"
tree$tip.label[tree$tip.label == "Paraconger caudilimbatus- leptocephalus larvae"] <- "Paraconger caudilimbatus, leptocephalus larvae"
tree$tip.label[tree$tip.label == "Ariosoma balearicum- leptocephalus larvae"] <- "Ariosoma balearicum, leptocephalus larvae"
tree$tip.label[tree$tip.label == "Gymnothorax saxicola- leptocephalus larvae"] <- "Gymnothorax saxicola, leptocephalus larvae"

rmr_data$spp[rmr_data$spp == "Negaprion brevirostris\xa0"] <- "Negaprion brevirostris"
rmr_data$spp[rmr_data$spp == "Mustelus antarcticus\xa0"] <- "Mustelus antarcticus"
rmr_data$spp[rmr_data$spp == "A. fimbria"] <- "Anoplopoma fimbria"
rmr_data$spp[rmr_data$spp == "Oncorhynchus\xa0tshawytscha"] <- "Oncorhynchus tshawytscha"
rmr_data$spp[rmr_data$spp == "Trachurus capensis\xa0"] <- "Trachurus capensis"
rmr_data$spp[rmr_data$spp == "Hypoatherina sp"] <- "Hypoatherina sp."
rmr_data$spp[rmr_data$spp == "Monacanthidae sp"] <- "Monacanthidae"
rmr_data$spp[rmr_data$spp == "Monacanthidae sp."] <- "Monacanthidae"
rmr_data$spp[rmr_data$spp == "Ambassis sp"] <- "Ambassis sp."
rmr_data$spp[rmr_data$spp == "Caulophrynidae sp."] <- "Caulophrynidae"
rmr_data$spp[rmr_data$spp == "Scomberesocidae sp."] <- "Scomberesocidae"
rmr_data$spp[rmr_data$spp == "Carcharodon carcharius"] <- "Carcharodon carcharias"
rmr_data$spp[rmr_data$spp == "Thunnys maccoyyi"] <- "Thunnus maccoyii"

# Especies en el dataset que faltan en el árbol
missing_in_tree <- setdiff(unique(rmr_data$spp), tree$tip.label)
missing_in_tree

# Especies en el árbol que no están en el dataset
extra_in_tree <- setdiff(tree$tip.label, unique(rmr_data$spp))
extra_in_tree

# Podar el árbol para eliminar especies extra
tree <- drop.tip(tree, extra_in_tree)

# Especies comunes entre ambos (dataset y árbol)
common_species <- intersect(unique(rmr_data$spp), tree$tip.label)
common_species

# Check if all dataset species are present in phylogeny
all_match <- all(sort(unique(rmr_data$spp)) %in% sort(tree$tip.label))
all_match  # TRUE if all dataset species are included in the tree

# --- housekeeping
rmr_data$method <- factor(rmr_data$method)
rmr_data$therm <- factor(rmr_data$therm)
rmr_data <- rmr_data[complete.cases(rmr_data),]

#rownames(rmr_data) <- rmr_data$spp
#spp <- rownames(rmr_data)

## Andrew code added 11 Aug 2025
# calculate natural log Ln(M) and Ln(RMR)
rmr_data$logM <- log(rmr_data$mass.kg)
rmr_data$logRMR <- log(rmr_data$RMR.mg.O2.h)

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
A <- ape::vcv.phylo(compute.brlen(tree, method = "Grafen"))

# adding species as a factor into our dataset
rmr_data$spp_name <- factor(rmr_data$spp)

# Many species include serveral masurements. Thus, we use a bayesian models.
# General STAN specifications 
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel::detectCores()

# MODEL
# We use a bayesian approach to study the relationship
# across species (phylogeny) while taking into considration
# intraspaecific variation and multiple observations for
# several species.

# set model priors
priors = c(
  prior(normal(0,10), "b"),
  prior(normal(0,50), "Intercept"),
  prior(student_t(3,0,20), "sd"),
  prior(student_t(3,0,20), "sigma")
)

# run Bayesian model including phylogenetic tree
set.seed(23); mod <- brm((logRMR)  ~ logM + T + therm # fixed coefficients
  + (1 | gr(spp, cov = A)) + (1| spp_name), # random factors
  family = gaussian(), 
  data2 = list(A = A),
  prior = priors,
  sample_prior = TRUE, 
  chains = 2, 
  cores = 16, # cores to be used 8 in mac pro
  iter = 3000, 
  warmup = 500,
  control = list(adapt_delta = 0.99, max_treedepth = 15), # Improve convergence
  save_pars = save_pars(all = TRUE),
  data = rmr_data)

# get coefficients
summary(mod) 

# check model convergence 

# par(oma = c(4,4,4,4))
# plot(mod)
# pp_check(mod)

# get bayesian p-value
# p_map(mod)

natural_mod <- mod
save(natural_mod, file = "tmp/nacho_RMR/natural_model.rda", compress = "xz")


