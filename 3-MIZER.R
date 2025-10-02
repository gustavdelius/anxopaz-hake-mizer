
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

ng <- length( gear_names)

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

initial_effort( hake_model) <- 1

hake_model <- matchYield( hake_model)
hake_model <- steady( hake_model)



# Fit ----------------------

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
source( './scripts/MIZER.R')

# nofixed <- c( 'a', 'b','beta', 'sigma', 'inter_HR', 'inter_HH', 'gamma', 'q', 'h', 'n', 
#               'ks', 'p', 'k', 'alpha', 'U', 'w_mat', 'w_min', 'w_max', 'M', 'd')

hake_model_fitted <- MIZER( model = hake_model, catch = corLFD)

hake_model_fitted_nf <- MIZER( model = hake_model, catch = corLFD, 
  fixed_sel = F, nofixed = c('h'))

hake_model_fitted_nf2 <- MIZER( model = hake_model_fitted, catch = corLFD, 
  fixed_sel = T, nofixed = c('h'))

hake_model_fitted_nf3 <- MIZER( model = hake_model_fitted, catch = corLFD, 
  fixed_sel = T, nofixed = c( 'gamma','q','h','n','ks','p','k','alpha','U','M','d'))

hake_model_fitted_nf4 <- MIZER( model = hake_model_fitted, catch = corLFD, 
  fixed_sel = T, nofixed = c( 'gamma','q','h','n','ks','p','k','alpha','U','M','d','kappa','lambda'))

hake_model_fitted@species_params$h
hake_model_fitted_nf@species_params$h
hake_model_fitted_nf2@species_params$h
hake_model_fitted_nf3@species_params$h
hake_model_fitted_nf4@species_params$h


## Check ---------------------------------------

getYield( hake_model_fitted)/sum( hake_model_fitted@gear_params$yield_observed)
getYield( hake_model_fitted_nf)/sum( hake_model_fitted@gear_params$yield_observed)
getYield( hake_model_fitted_nf2)/sum( hake_model_fitted@gear_params$yield_observed)
getYield( hake_model_fitted_nf3)/sum( hake_model_fitted@gear_params$yield_observed)
getYield( hake_model_fitted_nf4)/sum( hake_model_fitted@gear_params$yield_observed)

nofixed = c( 'gamma','q','h','n','ks','p','k','alpha','U','M','d')

model_vector1 <- as.numeric(species_params(hake_model_fitted)[nofixed])
model_vector2 <- as.numeric(species_params(hake_model_fitted_nf2)[nofixed])
model_vector3 <- as.numeric(species_params(hake_model_fitted_nf3)[nofixed])

biopars <- rbind( model_vector1, model_vector2, model_vector3)
colnames(biopars) <- nofixed
biopars

res_pars <- c('kappa','lambda')
model_vector1 <- as.numeric(resource_params(hake_model_fitted)[res_pars])
model_vector2 <- as.numeric(resource_params(hake_model_fitted_nf4)[res_pars])

respars <- rbind( model_vector1, model_vector2)
colnames(respars) <- res_pars
respars


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

modelo <- hake_mizer
save(modelo, file = "./fit.RData")
