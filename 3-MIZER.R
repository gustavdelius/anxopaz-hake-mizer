
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~   Hake's MIZER model  ~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(mizerMR)
library(TMB)

source( './scripts/aux_functions.R')


# Biological parameters ----------------------------

load( './input/Bio_Pars.RData')     # './1-Bio_Pars.R' with Biological Parameters


# Fishing Mortality --------------------------

load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data


# SSB ----------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

aver_y

ss_biomass <- assessment$SSB[assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

species_params(bio_pars)$biomass_observed <- sum(ss_biomass)/length(aver_y)
species_params(bio_pars)$biomass_cutoff <- lwf(4,a,b)   # SS smallest size is 4 cm

species_params(bio_pars)$biomass_observed/1e6
species_params(bio_pars)$biomass_cutoff

bio_pars <- setBevertonHolt( bio_pars,        # Rdd = Rdi * (Rmax/(Rdi+Rmax))
                reproduction_level = 0.001)   # rep_level = Rdd/Rmax (density dependance degree)



# MIZER model --------------------------------

hake_model <- bio_pars |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
  calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()

gear_names <- Catch$fleet

### Double sigmoid selectivity initial parameters

gear_params( hake_model) <- data.frame(
  gear = gear_names, species = "Hake", catchability = 1,
  sel_func = "double_sigmoid_length",
  l50 = c(       28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
  l25 = c(       23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
  l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, 35.6, 58.0, 54.4, 27.3),
  l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, 45.9, 67.9, 60.8, 28.4))

### Catch by gear

gear_params( hake_model)$yield_observed <- corLFDs$catch   # == Catch$catch; != LFDs$catch

gear_params( hake_model)$yield_observed/1e6
sum(gear_params( hake_model)$yield_observed/1e6)


### Initial effort

initial_effort( hake_model) <- .25

hake_model <- matchYield( hake_model)
hake_model <- steady( hake_model)



# Fit ----------------------

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
source( './scripts/MIZER.R')

# nofixed <- c('h','n','ks','p','k','f0','alpha')
nofixed <- c('h','ks','p','k','f0','alpha','n')

# hake_model_fitted <- MIZER( model = hake_model, catch = corLFD, nofixed = nofixed)
hake_model_fitted <- MIZER( model = hake_model, catch = corLFD)


## Check ---------------------------------------

plot_lfd( hake_model_fitted, corLFD)
plot_lfd_gear( hake_model_fitted, corLFD)

getYield( hake_model_fitted)
sum( hake_model_fitted@gear_params$yield_observed)
getYield( hake_model_fitted)/sum( hake_model_fitted@gear_params$yield_observed)

getBiomass( hake_model_fitted)
hake_model_fitted@species_params$biomass_observed
getBiomass( hake_model_fitted)/hake_model_fitted@species_params$biomass_observed

# model_vector <- as.numeric(species_params(hake_model_fitted)[nofixed])
# pre_vector <- as.numeric(species_params(hake_model)[nofixed])
# 
# biopars <- rbind( model_vector, pre_vector)
# colnames(biopars) <- nofixed
# biopars
# 
# hake_model_fitted@gear_params
# 
# plotSpectra( hake_model_fitted, power = 2) + theme_bw() 
# 
# 
# ### HERE!!!! ---------------
# 
# ener <- getEnergy( hake_model_fitted, return_df = TRUE)
# getEnergy( hake_model_fitted, log = FALSE)
# 
# sp <- species_params( hake_model_fitted)
# wt <- w( hake_model_fitted)
# 
# 
# ## No update??
# 
# getMetabolicRate(hake_model_fitted)
# getMetabolicRate(hake_model_fitted)/(hake_model@species_params$ks*(wt^hake_model@species_params$p))
# getMetabolicRate(hake_model_fitted)/(sp$ks*(wt^sp$p))
# getMetabolicRate(hake_model_fitted)/getMetabolicRate(hake_model)
# getMetabolicRate(bio_pars)/getMetabolicRate(hake_model)
# 
# 
# 
# ## Then: old ---------------
# 
# getEnergy( hake_model)
# getEnergy( hake_model, log = TRUE)
# 
# sp <- species_params( hake_model)
# wt <- w( hake_model)
# 
# emetab <- sp$ks * wt^sp$p
# eactiv <- sp$k * wt
# etotal <- sp$alpha*sp$f0*(sp$h*wt^sp$n)
# ereproandgrowth <- etotal - emetab - eactiv
# repp <- (wt/sp$w_max)^(1-sp$n)
# phi <-  (wt/sp$w_max)^(1-sp$n) * (1/(1+(wt/sp$w_mat)^(-sp$U)))
# egrowth <- ereproandgrowth * (1-phi)
# erepro <- ereproandgrowth * (phi)
# 
# getReproductionProportion(hake_model)/repp
# getEGrowth(hake_model)/egrowth
# getMetabolicRate(hake_model)/emetab
# getEReproAndGrowth(hake_model)/ereproandgrowth
# getERepro(hake_model)/erepro


# Background ---------------

hake_mizer <- scaleDownBackground( hake_model_fitted, 2e-7)

plotSpectra( hake_mizer, power = 2) + theme_bw() 



# # Old MIZER fit --------------------------------
# 
# ngear <- length(gear_names)
# 
# catch_lengths <- data.frame( species = "Hake",gear = rep(gear_names,each=bins_no), 
#   length = rep(LFDc$length,ngear), dl = 1, weight = rep(w(hake_model),ngear), 
#   dw = rep(dw(hake_model),ngear), catch = c(LFD$number))
# 
# # tuneParams( hake_mizer, catch = catch_lengths)
# # 
# # hake_model <- readParams("./output/hake_model.rds")
# # hake_model <- scaleDownBackground( hake_model, 1/8000000)



# Save ----------------------

save.image( './output/hake_model.RData')

