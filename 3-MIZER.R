
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


### Resource dynamics

# hake_model@resource_params$kappa <- 5000000000 
# hake_model@species_params$gamma <- 0.000001 

### Steady state

hake_model <- matchYield( hake_model)
hake_model <- steady( hake_model)



# Fit ----------------------

pd_name <- paste0( getwd(), '/output/plots/base/')
dir.create( path = pd_name, showWarnings = TRUE, recursive = TRUE)

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
source( './scripts/MIZER.R')

# nofixed <- c( 'a', 'b','beta', 'sigma', 'inter_HR', 'inter_HH', 'gamma', 'q', 'h', 'n', 
#               'ks', 'p', 'k', 'alpha', 'U', 'w_mat', 'w_min', 'w_max', 'M', 'd')

only_sel_mod <- MIZER( model = hake_model, catch = corLFD,
  plot = T, plot_dir = paste0( pd_name, 'only_sel/'))

sel_and_gamma_mod <- MIZER( model = hake_model, catch = corLFD,
  fixed_sel = F, nofixed = c('gamma'),
  plot = T, plot_dir = paste0( pd_name, 'sel_and_gamma/'))

only_gamma_mod <- MIZER( model = only_sel_mod, catch = corLFD,
  fixed_sel = T, nofixed = c('gamma'),
  plot = T, plot_dir = paste0( pd_name, 'only_gamma/'))

only_others_mod <- MIZER( model = only_sel_mod, catch = corLFD,
  fixed_sel = T, nofixed = c( 'gamma','q','n','ks','p','k','alpha'),
  plot = T, plot_dir = paste0( pd_name, 'only_others/'))

only_others_resources_mod <- MIZER( model = only_sel_mod, catch = corLFD,
  fixed_sel = T, nofixed = c( 'gamma','q','n','ks','p','k','alpha','kappa','lambda'),
  plot = T, plot_dir = paste0( pd_name, 'only_others_resources/'))

only_resources_mod <- MIZER( model = only_sel_mod, catch = corLFD,
  fixed_sel = T, nofixed = c( 'kappa','lambda'),
  plot = T, plot_dir = paste0( pd_name, 'only_resources/'))

all_fitted_mod <- MIZER( model = hake_model, catch = corLFD,
  fixed_sel = F, nofixed = c( 'gamma','q','n','ks','p','k','alpha','kappa','lambda'),
  plot = T, plot_dir = paste0( pd_name, 'all_fitted/'))

hake_model_fitted <- only_sel_mod


## Comparison --------------------

mod_list <- list( only_sel=only_sel_mod, sel_and_gamma=sel_and_gamma_mod, 
  only_gamma=only_gamma_mod, only_others=only_others_mod, only_others_resources=only_others_resources_mod,
  only_resources=only_resources_mod, all_fitted=all_fitted_mod)

nofixed <- c( 'gamma','q','n','ks','p','k','alpha')
res_pars <- c('kappa','lambda')
all_pars <- c( 'inc_SSB', 'inc_Y', nofixed, res_pars)

pars_table <- matrix( NA, nrow = length(mod_list), ncol = length(all_pars), 
  dimnames = list( names(mod_list), all_pars))

for(i in names(mod_list)){
  pars_table[i,] <- c( 
    signif( getBiomass( mod_list[[i]])/hake_model_fitted@species_params$biomass_observed,3),
    signif( (getYield( mod_list[[i]])/sum( hake_model_fitted@gear_params$yield_observed)),3),
    signif( as.numeric(species_params(mod_list[[i]])[nofixed]),3),
    signif( as.numeric(resource_params(mod_list[[i]])[res_pars]),3))
}

pars_table

spectradf <- NULL

for(i in names(mod_list)){
  idf <- plotSpectra( mod_list[[i]], power = 2)
  idf <- idf$data[which(idf$data$Species=='Hake'),c('w','value')]; idf$model <- i
  spectradf <- rbind( spectradf, idf)}

ggplot( spectradf, aes( x = w, y = value, color = model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparison.jpg'), width = 12, height = 7)

ggplot( spectradf, aes( x = w, y = value, color = model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + 
  scale_y_log10( limits = c( 100000000, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparisonx10.jpg'), width = 12, height = 7)


  
## Background ---------------

plotSpectra( hake_model_fitted, power = 2) + theme_bw() 
ggsave( paste0( pd_name, 'base_rr.jpg'), width = 12, height = 7)

hake_mizer <- scaleDownBackground( hake_model_fitted, 2e-7)

plotSpectra( hake_mizer, power = 2) + theme_bw() 
ggsave( paste0( pd_name, 'base_pluskappa.jpg'), width = 12, height = 7)

pars_table <- rbind( pars_table, rescaled_res = c( 
  signif( getBiomass( hake_mizer)/hake_model_fitted@species_params$biomass_observed,3),
  signif( getYield( hake_mizer)/sum( hake_model_fitted@gear_params$yield_observed),3),
  signif( as.numeric(species_params(hake_mizer)[nofixed]),3),
  signif( as.numeric(resource_params(hake_mizer)[res_pars]),3)))

pars_table


# # Old MIZER fit
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



# Cannibalism --------------------------------------

load( './data/Diet.RData')

cannibal <- hake_mizer
# cannibal2 <- hake_model_fitted


## Turn on cannibalism

pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% aver_y)])

interaction_matrix( cannibal)[] <- pcann
# interaction_matrix( cannibal2)[] <- pcann

ext_mort( cannibal) <- ext_mort( cannibal) - getPredMort( cannibal)
# ext_mort( cannibal2) <- ext_mort( cannibal2) - getPredMort( cannibal2)

# p_test <- steadySingleSpecies( cannibal)
# all.equal( initialN(cannibal), initialN(p_test), tolerance = 1e-5)

cdir <- paste0( getwd(), '/output/cannibal/')
dir.create( path = cdir, showWarnings = TRUE, recursive = TRUE)


## Fit

cannibal_hake <- MIZER( model = cannibal, catch = corLFD, plot = T, plot_dir = cdir)
# cannibal_hake2 <- MIZER( model = cannibal2, catch = corLFD, plot = T, plot_dir = paste0(cdir,'/2'))


## Add to Comparison

pars_table <- rbind( pars_table, cannibal = c( 
  signif( getBiomass( cannibal_hake)/hake_model_fitted@species_params$biomass_observed,3),
  signif( getYield( cannibal_hake)/sum( hake_model_fitted@gear_params$yield_observed),3),
  signif( as.numeric(species_params(cannibal_hake)[nofixed]),3),
  signif( as.numeric(resource_params(cannibal_hake)[res_pars]),3)))

pars_table

plotDiet( hake_mizer)
plotDiet( cannibal_hake)
# plotDiet( cannibal_hake2)


s1 <- plotSpectra2(hake_mizer, name1 = "Base Model", cannibal_hake, name2 = "With cannibalism",
  power = 2, resource = FALSE)

s2 <- plotSpectra2(cannibal, name1 = "No refit", hake_mizer, name2 = "Base Model",
  power = 2, resource = FALSE)

stot <- rbind( s1$data, s2$data[which(s2$data$Model=='No refit'),])

ggplot( stot, aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + 
  scale_y_log10( limits = c( 100000000, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( cdir, 'spectra_comparison.jpg'), width = 12, height = 7)


p1 <- plotDeath( hake_mizer, return_data = T)
p2 <- plotDeath( cannibal_hake, return_data = T)

p1$Model = 'No cannibalism'
p2$Model = 'With cannibalism'

ppt <- rbind( p1, p2)
ppt <- ppt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))

pdeath <- ggplot(ppt, aes(x = w, y = value, fill = Cause)) +
  geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
  facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
  scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                     labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))

pdeath
ggsave( paste0( cdir, 'death_comparison.jpg'), width = 12, height = 7)
# ggsave( paste0(cdir,'death_poster.jpg'), width = 5, height = 3)


l1 <- plot_lfd( hake_mizer, corLFD, return_df = T)
l2 <- plot_lfd( cannibal_hake, corLFD, return_df = T)

l1$Model <- 'No cannibalism'
l2$Model <- 'With cannibalism'

ldf1 <- rbind( l1, subset(l2, Type == 'Estimated'))

ldfp1 <- ggplot(ldf1, aes(x = Length, y = Density)) +
  geom_bar(data = subset(ldf1, Type != 'Estimated'), aes( fill = Type), stat = "identity", position = "dodge", alpha = 0.6) +
  geom_line(data = subset(ldf1, Type == 'Estimated'), aes( color = Model), linewidth = 1) +
  scale_fill_manual(values = c("Observed" = "yellowgreen")) +
  scale_color_manual(values = c( "No cannibalism" = "#E41A1C", "With cannibalism" = "#377EB8")) +
  theme_bw() + labs( x = "Size [cm]", y = "Normalised number density [1/cm]", fill = NULL, color = NULL)

ldfp1
ggsave( paste0( cdir,'ldf.jpg'), width = 9, height = 7)
# ggsave( paste0( cdir,'ldf_poster.jpg'), width = 5, height = 3)



# Time periods --------------------------------

# Loop for time periods

tpdir <- paste0( getwd(), '/output/time_periods/')
dir.create( path = tpdir, showWarnings = TRUE, recursive = TRUE)

tp_mods <- list()

years

for( i in 1:length(years)){
  
  ny <- names(years)[i]
  vy <- years[[i]]
  ychar <- paste0( vy[1],'-',vy[length(vy)])
  
  ssbio_tp <- sum(assessment$SSB[assessment$Year %in% vy]*1e6)/length(vy)
  
  itp_mod <- hake_mizer
  
  species_params(itp_mod)$biomass_observed <- ssbio_tp
  
  itp_mod <- itp_mod |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
  
  gear_params( itp_mod)$yield_observed <- corLFD_sum[[i]]$catch
  
  tp_mods[[ychar]] <- MIZER( model =  itp_mod, catch = corLFD_list[[i]],
     plot = T, plot_dir = paste0(ychar,'/'))
  
  plotSpectra( tp_mods[[ychar]], power = 2) + theme_bw()
  ggsave(paste0(tpdir,ychar,'/spectra.jpg'), width = 9, height = 7) 
  
}


spdf <- plotSpectra( tp_mods[[1]], power = 2)
spdf <- subset( spdf$data, Legend == 'Resource')[,c(1,2)]
spdf$Model <- 'Resource'

for(i in names(tp_mods)){
  ispdf <- plotSpectra( tp_mods[[i]], power = 2)
  ispdf <- subset(ispdf$data, Legend != 'Resource')[,c(1,2)]  
  ispdf$Model <- i
  spdf <- rbind( spdf, ispdf)
}

ggplot( spdf, aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + 
  scale_y_log10( limits = c( 100000000, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Period')

ggsave( paste0( tpdir, 'spectra_comparison_res.jpg'), width = 12, height = 7)

spdf <- subset(spdf, Model!= 'Resource')

ggplot( spdf, aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + 
  scale_y_log10( limits = c( 100000000, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Period')

ggsave( paste0( tpdir, 'spectra_comparison.jpg'), width = 12, height = 7)


tp_table <- matrix( NA, nrow = length(tp_mods), ncol = length(all_pars), 
                      dimnames = list( names(tp_mods), all_pars))

for(i in names(tp_mods)){
  tp_table[i,] <- c( 
    signif( getBiomass( tp_mods[[i]])/tp_mods[[i]]@species_params$biomass_observed,3),
    signif( (getYield( tp_mods[[i]])/sum( tp_mods[[i]]@gear_params$yield_observed)),3),
    signif( as.numeric(species_params(tp_mods[[i]])[nofixed]),3),
    signif( as.numeric(resource_params(tp_mods[[i]])[res_pars]),3))
}

tp_table


tpcdir <- paste0( tpdir, 'cannibal/')
dir.create( path = tpcdir, showWarnings = TRUE, recursive = TRUE)

tp_cann <- list()
tpc_table <- NULL

for(i in 1:length(tp_mods)){
  
  icann <- tp_mods[[i]]
  iname <- names(tp_mods)[[i]]
  
  vy <- years[[i]]
  
  pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% vy)])
  interaction_matrix( icann)[] <- pcann
  ext_mort( icann) <- ext_mort( icann) - getPredMort( icann)
  
  tp_cann[[iname]] <- icann <- MIZER( model = icann, catch = corLFD_list[[i]], 
    plot = T, plot_dir = paste0(tpcdir,iname,'/'))
  
  tpc_table <- rbind( tpc_table, c( 
    signif( getBiomass( icann)/icann@species_params$biomass_observed,3),
    signif( getYield( icann)/sum( icann@gear_params$yield_observed),3),
    signif( as.numeric(species_params(icann)[nofixed]),3),
    signif( as.numeric(resource_params(icann)[res_pars]),3)))
  
  s1 <- plotSpectra2( tp_mods[[i]], name1 = "No cannibalism", tp_cann[[i]], name2 = "With cannibalism",
    power = 2, resource = FALSE)
  
  ggplot( s1$data, aes( x = w, y = value, color = Model)) + 
    geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + 
    theme_bw() +
    labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')
  
  ggsave( paste0( tpcdir,iname, '-spectra_comparison.jpg'), width = 12, height = 7)
  
  ip1 <- plotDeath( tp_mods[[i]], return_data = T)
  ip2 <- plotDeath( tp_cann[[i]], return_data = T)
  
  ip1$Model = 'No cannibalism'
  ip2$Model = 'With cannibalism'
  
  ippt <- rbind( ip1, ip2)
  ippt <- ippt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))
  
  ipdeath <- ggplot(ippt, aes(x = w, y = value, fill = Cause)) +
    geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
    facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
    scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                       labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))
  
  ipdeath

  ggsave( paste0( tpcdir,iname, '-death_comparison.jpg'), width = 12, height = 7)
  
}

rownames(tpc_table) <- paste0('cann_',names(tp_mods))

alltp_tab <- rbind( tp_table, tpc_table)
alltp_tab

save.image( './output/alldata.RData')




