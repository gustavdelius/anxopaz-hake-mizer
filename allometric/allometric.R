# This script sets up a Hake model 
rm(list = ls())

library(dplyr)
library(ggplot2)
library(mizerExperimental)

# We will use functions from the mizerEcopath package. 
# That package is still undergoing rapid change
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

# Species params ----

sp <- bio_pars@species_params |>
    select(species, w_mat, age_mat, w_max, a, b, n,
           pred_kernel_type, beta, sigma)

# Biomass
quantity <- 'biomass' # 'SSB'
ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)
obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm
sp$biomass_observed <- obs_q
sp$biomass_cutoff <- b_min

# Max size (from landings)
l_max <- 1.001 * max(corLFD$length + 1, na.rm = TRUE)
sp$w_max = lwf(l_max, sp$a, sp$b)

# Mortality exponent (from SS)
sp$d <- -0.1217

# Allometric params ----
# We now set up a Hake model with power law encounter rate and
# power law mortality rate.
hake_model <- newAllometricParams(sp)
# This is not yet calibrated

# Fishing ----
# We work with only one "Total" gear
gp <- data.frame(
    gear = "Total", species = "Hake", catchability = 1,
    sel_func = "sigmoid_length",
    l50 = 30,
    l25 = 28,
    yield_observed = sum(corLFDs$catch)
)
gear_params(hake_model) <- gp
initial_effort(hake_model) <- 1

# Set initial catchability to get the observed yield
# This will be calibrated properly in `matchCatch()` below.
yield <- getYield(hake_model)
gp$catchability <- gp$yield_observed / yield
gear_params(hake_model) <- gp

# Transform landings data to required shape, adding landings from
# all gears to give total landings
catch <- corLFD |>
    group_by(length) |>
    summarise(catch = sum(number))
catch$dl = 1
catch$species = "Hake"
catch$gear = "Total"

# Calibrate to total yield and landings size distribution
hake_model <- matchCatch(hake_model, catch = catch)

plotYieldVsSize(hake_model, x_var = "Length", catch = catch)

# Add density dependencies ----
# We now have a model whose steady state matches the landings data
# But we still need to calibrate its sensitivity to changes away from
# the steady state, for example its sensitivity to changes in fishing.
# We don't have a good way to choose the following three parameters yet,
# so we'll just make up values for now.

# Set reproduction level
hake_model <- setBevertonHolt(hake_model, reproduction_level = 0.6)

# Set feeding level
hake_model <- setFeedingLevel(hake_model, 0.6)

# Set resource level
hake_model <- alignResource(hake_model) |>
    setResourceInteraction(resource_dynamics = "resource_semichemostat")
resource_level(hake_model) <- 1/2

# Turn on cannibalism ----
# I'll assume below that 17% of the total diet comes from cannibalism
# You will get a warning if you try to increase this and we should discuss this.
diet_matrix <- matrix(c(0.17, 0.83), ncol = 2,
                      dimnames = list(predator = "Hake",
                                      prey = c("Hake", "other")))

cannibal_hake_model <- matchDiet(hake_model, diet_matrix)
plotDietX(cannibal_hake_model)

# Check that the models are in steady state ----
sim <- project(hake_model, t_max = 8)
plotBiomass(sim)
sim <- project(cannibal_hake_model, t_max = 8)
plotBiomass(sim)

# Run simulations with increased fishing effort ----
sim12 <- project(hake_model, effort = 1.2, t_max = 12, t_save = 0.2)
sim_cannibal12 <- project(cannibal_hake_model, effort = 1.2, t_max = 12, t_save = 0.2)

# Compare biomasses
bio12 <- melt(getBiomass(sim12))
bio12$Cannibalism <- "Off"
bio_cannibal12 <- melt(getBiomass(sim_cannibal12))
bio_cannibal12$Cannibalism <- "On"
bio <- rbind(bio12, bio_cannibal12)
ggplot(bio) +
    geom_line(aes(x = time, y = value, 
                  colour = Cannibalism)) 
