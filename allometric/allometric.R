
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

install_github("gustavdelius/mizerEcopath")
library(mizerEcopath)

source( './scripts/aux_functions.R')

# Load data ----

# Biological parameters
load( './input/Bio_Pars.RData')     # './1-Bio_Pars.R' with Biological Parameters

# Fishing Mortality
load( './input/Catch.RData')   # './2-Catch.R' with Catch and LFD data

# SSB / Bio
load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

quantity <- 'biomass' # 'SSB'

ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm

# Allometric params ----
sp <- bio_pars@species_params |>
    select(species, w_mat, age_mat, w_max, a, b)
sp$biomass_observed <- obs_q
sp$biomass_cutoff <- b_min

hake_model <- newAllometricParams(sp)
plotSpectra(hake_model, resource = FALSE)

# Gears ----

gear_names <- Catch$fleet

### Double sigmoid selectivity initial parameters

ng <- length( gear_names)

gear_params(hake_model) <- data.frame(
    gear = gear_names, species = "Hake", catchability = 1,
    sel_func = "double_sigmoid_length",
    l50 = c(       28.6, 30.7, 29.8, 14.8, 27.5, 30.3, 51.2, 54.9, 16.1),
    l25 = c(       23.8, 28.5, 27.4, 13.0, 26.6, 28.1, 47.5, 51.3, 13.5),
    l50_right = c( 38.3, 33.6, 42.0, 20.6, 33.1, 35.6, 58.0, 54.4, 27.3),
    l25_right = c( 43.3, 45.0, 47.8, 27.0, 38.9, 45.9, 67.9, 60.8, 28.4))


### Catch by gear

gear_params(hake_model)$yield_observed <- corLFDs$catch   # == Catch$catch; != LFDs$catch

gear_params( hake_model)$yield_observed/1e6
sum(gear_params( hake_model)$yield_observed/1e6)
