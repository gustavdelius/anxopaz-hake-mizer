
rm(list = ls())

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
# Transform landings data to required shape
catch <- corLFD |>
    group_by(length) |>
    summarise(catch = sum(number))
catch$dl = 1
catch$species = "Hake"
catch$gear = "Total"

# SSB / Bio
load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

# Species params ----

sp <- bio_pars@species_params |>
    select(species, w_mat, age_mat, w_max, a, b)

# Biomass
quantity <- 'biomass' # 'SSB'
ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)
obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm

sp$biomass_observed <- obs_q
sp$biomass_cutoff <- b_min

# max size
l_max = 1.1 * max(catch$length + catch$dl, na.rm = TRUE)
sp$w_max = lwf(l_max, sp$a, sp$b)


# Allometric params ----
hake_model <- newAllometricParams(sp) |>
    matchGrowth() |>
    matchBiomasses()
plotSpectra(hake_model, power = 2, resource = FALSE)

# Fishing ----

gp <- data.frame(
    gear = "Total", species = "Hake", catchability = 1,
    sel_func = "sigmoid_length",
    l50 = 30,
    l25 = 28,
    yield_observed = sum(corLFDs$catch)
)

gear_params(hake_model) <- gp
initial_effort(hake_model) <- 1

# Set catchability to get the observed yield
yield <- getYield(hake_model)
gp$catchability <- gp$yield_observed / yield
gear_params(hake_model) <- gp

hake_model <- matchCatch(hake_model, catch = catch)

plotYieldVsSize(hake_model, x_var = "Length", catch = catch)

hake_model@species_params$erepro

plotGrowthCurves(hake_model)
