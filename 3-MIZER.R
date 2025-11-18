
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


# SSB / Bio ----------------

load( './input/Hake_SS_Data.RData')   # './scripts/WGBIE24.R' WGBIE assessment results

quantity <- 'biomass' # 'SSB'

ss_biomass <- assessment[,quantity][assessment$Year %in% aver_y]*1e6  # SS biomass (tonnes to grams)

obs_q <- sum(ss_biomass)/length(aver_y); obs_q/1e6 
b_min <- lwf(4,a,b); b_min   # SS smallest size is 4 cm

species_params(bio_pars)$biomass_observed <- obs_q
species_params(bio_pars)$biomass_cutoff <- b_min

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

pd_name <- paste0( getwd(), '/output/OMs/')
dir.create( path = pd_name, showWarnings = TRUE, recursive = TRUE)

# TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
source( './scripts/MIZER.R')

# nofixed <- c( 'a', 'b','beta', 'sigma', 'inter_HR', 'inter_HH', 'gamma', 'q', 'h', 'n',
#               'ks', 'p', 'k', 'alpha', 'U', 'w_mat', 'w_min', 'w_max', 'M', 'd')


## OMs ---------------------------------

compiler <- !file.exists( "./TMB/fit.o")

only_sel_mod <- MIZER( model = hake_model, catch = corLFD, compiler = compiler,
                       plot = T, plot_dir = paste0( pd_name, 'only_sel/'))

sel_and_gamma_mod <- MIZER( model = hake_model, catch = corLFD, compiler = compiler,
                            fixed_sel = F, nofixed = c('gamma'),
                            plot = T, plot_dir = paste0( pd_name, 'sel_and_gamma/'))

only_gamma_mod <- MIZER( model = only_sel_mod, catch = corLFD, compiler = compiler,
                         fixed_sel = T, nofixed = c('gamma'),
                         plot = T, plot_dir = paste0( pd_name, 'only_gamma/'))

only_others_mod <- MIZER( model = only_sel_mod, catch = corLFD, compiler = compiler,
                          fixed_sel = T, nofixed = c( 'gamma','q','n','ks','p','k','alpha'),
                          plot = T, plot_dir = paste0( pd_name, 'only_others/'))

only_others_resources_mod <- MIZER( model = only_sel_mod, catch = corLFD, compiler = compiler,
                                    fixed_sel = T, nofixed = c( 'gamma','q','n','ks','p','k','alpha','kappa','lambda'),
                                    plot = T, plot_dir = paste0( pd_name, 'only_others_resources/'))

only_resources_mod <- MIZER( model = only_sel_mod, catch = corLFD, compiler = compiler,
                             fixed_sel = T, nofixed = c( 'kappa','lambda'),
                             plot = T, plot_dir = paste0( pd_name, 'only_resources/'))

all_fitted_mod <- MIZER( model = hake_model, catch = corLFD, compiler = compiler,
                         fixed_sel = F, nofixed = c( 'gamma','q','n','ks','p','k','alpha','kappa','lambda'),
                         plot = T, plot_dir = paste0( pd_name, 'all_fitted/'))


### Comparison --------------------

mod_list <- list( only_sel = only_sel_mod, sel_and_gamma = sel_and_gamma_mod,
                  only_gamma = only_gamma_mod, only_others = only_others_mod, 
                  only_others_resources = only_others_resources_mod,
                  only_resources = only_resources_mod, all_fitted = all_fitted_mod)

nofixed <- c( 'gamma', 'q', 'n', 'ks', 'p', 'k', 'alpha')
res_pars <- c( 'kappa', 'lambda')
all_pars <- c( 'inc_SSB', 'inc_Y', nofixed, res_pars)

pars_table <- matrix( NA, nrow = length(mod_list), ncol = length(all_pars), 
                      dimnames = list( names(mod_list), all_pars))

spectradf <- NULL

for(i in names(mod_list)){
  pars_table[i,] <- parsf( mod_list[[i]])
  idf <- spf(mod_list[[i]], name = i)
  spectradf <- rbind( spectradf, idf)}

pars_table

ggplot( spectradf, aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparison.jpg'), width = 12, height = 7)

ggplot( spectradf %>% filter(w>10), aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparisonx10.jpg'), width = 12, height = 7)

ggplot( spectradf, aes( x = w, y = value2, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparison2.jpg'), width = 12, height = 7)

ggplot( spectradf %>% filter(w>10), aes( x = w, y = value2, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( pd_name, 'spectra_comparison2x10.jpg'), width = 12, height = 7)



# OM selection ------------

hake_model_fitted <- only_sel_mod  


## Background ---------------

plotSpectra( hake_model_fitted, power = 2) + theme_bw() 
ggsave( paste0( pd_name, 'base_rr.jpg'), width = 12, height = 7)

hake_mizer <- scaleDownBackground( hake_model_fitted, 2e-8)

plotSpectra( hake_mizer) + theme_bw() 
ggsave( paste0( pd_name, 'base_pluskappa.jpg'), width = 12, height = 7)

plotSpectra( hake_mizer, power = 2) + theme_bw() 
ggsave( paste0( pd_name, 'base_pluskappa2.jpg'), width = 12, height = 7)

rescaled <- parsf( hake_mizer)
pars_table <- rbind( pars_table, rescaled = rescaled); pars_table


## Background rescaled fit ---------------------

hake_mizer2 <- MIZER( model = hake_mizer, catch = corLFD, compiler = compiler,
                      plot = T, plot_dir = paste0( pd_name, 'after_res/'))

hake_mizer3 <- MIZER( model = hake_mizer2, catch = corLFD, compiler = compiler,
                      nofixed = c('gamma', 'p'),
                      plot = T, plot_dir = paste0( pd_name, 'after_res_gamma/'))

hake_mizer4 <- hake_model_fitted
hake_mizer4@resource_params$kappa <- hake_mizer@resource_params$kappa

hake_mizer4 <- MIZER( model = hake_mizer4, catch = corLFD, compiler = compiler,
                      nofixed = c('gamma'),
                      plot = T, plot_dir = paste0( pd_name, 'after_res_gamma2/'))

after_rescaled <- parsf( hake_mizer2)
pars_table <- rbind( pars_table, after_rescaled = after_rescaled); pars_table

after_rescaled_gamma <- parsf( hake_mizer3)
pars_table <- rbind( pars_table, after_rescaled_gamma = after_rescaled_gamma); pars_table

after_kappa_gamma <- parsf( hake_mizer4)
pars_table <- rbind( pars_table, after_kappa_gamma = after_kappa_gamma); pars_table

hake_mizer5 <- hake_model
hake_mizer5@resource_params$kappa <- hake_mizer@resource_params$kappa

hake_mizer5 <- MIZER( model = hake_mizer5, catch = corLFD, nofixed = c('gamma'), compiler = compiler,
                      plot = T, plot_dir = paste0( pd_name, 'after_res_gamma3/'))

after_kappa_gamma2 <- parsf( hake_mizer5)
pars_table <- rbind( pars_table, after_kappa_gamma2 = after_kappa_gamma2); pars_table


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


## Save base fit----------------------

save.image( './output/hake_model.RData')

modelo <- hake_mizer
save( modelo, file = "./fit.RData")



# Cannibalism --------------------------------------

load( './data/Diet.RData')

cannibal <- cannibal2 <- hake_mizer

pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% aver_y)]); pcann

interaction_matrix( cannibal)[] <- interaction_matrix( cannibal2)[] <- pcann

ext_mort( cannibal) <- ext_mort( cannibal) - getPredMort( cannibal)

p_test <- steadySingleSpecies( cannibal)
all.equal( initialN(cannibal), initialN(p_test), tolerance = 1e-5)

cdir <- paste0( getwd(), '/output/cannibal/')
dir.create( path = cdir, showWarnings = TRUE, recursive = TRUE)


## Fit

cannibal_hake <- MIZER( model = cannibal, catch = corLFD, compiler = compiler, plot = T, plot_dir = cdir)

compiler2 <- !file.exists( dynlib("./TMB/fit_cann"))

cannibal_hake2 <- MIZER( model = cannibal, catch = corLFD, compiler = compiler2,
                         plot = T, plot_dir = paste0( cdir, '2/'), cannibalism = 'diff')


## Add to Comparison

base_model <- pars_table['rescaled',]

cpars_table <- rbind( base_model = base_model, cannibal = parsf( cannibal_hake))
cpars_table        # equal

plotDiet( hake_mizer) + scale_x_log10( limits = c( 10, NA))
plotDiet( cannibal_hake) + scale_x_log10( limits = c( 10, NA))
plotDiet( cannibal_hake2) + scale_x_log10( limits = c( 10, NA))

bmsp <- spf( hake_mizer, name = 'Base Model')
csp1 <- spf( cannibal, name = 'No Refit')
csp2 <- spf( cannibal_hake, name = 'With Cannibalism')
csp3 <- spf( cannibal_hake2, name = 'With Cannibalism (rare)')

stot <- rbind( bmsp, csp1, csp2, csp3) %>% filter( w > 10)

ggplot( stot, aes( x = w, y = value, color = Model)) + theme_bw() +
  geom_line( linewidth = .8) + scale_x_log10() + 
  labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')

ggsave( paste0( cdir, 'spectra_comparison.jpg'), width = 12, height = 7)

ggplot( stot, aes( x = w, y = value2, color = Model)) + theme_bw() +
  geom_line( linewidth = .8) + scale_x_log10() + 
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( cdir, 'spectra_comparison2.jpg'), width = 12, height = 7)


p1 <- plotDeath( hake_mizer, return_data = T)
p2 <- plotDeath( cannibal_hake, return_data = T)
p3 <- plotDeath( cannibal_hake2, return_data = T)

p1$Model = 'Base Model'
p2$Model = 'With Cannibalism'
p3$Model = 'With Cannibalism (rare)'

ppt <- rbind( p1, p2, p3)
ppt <- ppt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))
ppt <- subset( ppt, w > 1)

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

l1$Model <- 'Base Model'
l2$Model <- 'With Cannibalism'

ldf1 <- rbind( l1, subset(l2, Type == 'Estimated'))

ldfp1 <- ggplot(ldf1, aes(x = Length, y = Density)) +
  geom_bar(data = subset(ldf1, Type != 'Estimated'), aes( fill = Type), stat = "identity", position = "dodge", alpha = 0.6) +
  geom_line(data = subset(ldf1, Type == 'Estimated'), aes( color = Model), linewidth = 1) +
  scale_fill_manual(values = c("Observed" = "#00BFC4")) +
  scale_color_manual(values = c( "Base Model" = "#F8766D", "With Cannibalism" = "#7CAE00")) +
  theme_bw() + labs( x = "Size [cm]", y = "Normalised number density [1/cm]", fill = NULL, color = NULL)

ldfp1

ggsave( paste0( cdir,'ldf.jpg'), width = 9, height = 7)



# Natural Mortality ----------------------------------

tpnmdir <- paste0( getwd(), '/output/natural_mortality/')
dir.create( path = tpnmdir, showWarnings = TRUE, recursive = TRUE)

natmmod <- natmmod2 <- natmmod3 <- hake_mizer

load( './data/Natural_Mortality.RData')   # './scripts/Natural_Mortality.R' results

mu0_lm <- NatM_pars['lm','mu0']
d_lm <- NatM_pars['lm','d']
lm_mort <- mu0_lm*w(hake_mizer)^(d_lm)

mu0_nls <- NatM_pars['nls','mu0']
d_nls <- NatM_pars['nls','d']
nls_mort <- mu0_nls*w(hake_mizer)^(d_nls)

mu0_DB <- NatM_pars['ad_hoc','mu0']
d_DB <- NatM_pars['ad_hoc','d']
DB_mort <- mu0_DB*w(hake_mizer)^(d_DB)

ext_mort( natmmod) <- array( lm_mort, dim=c(1,bins_no))
ext_mort( natmmod2) <- array( DB_mort, dim=c(1,bins_no))
ext_mort( natmmod3) <- array( nls_mort, dim=c(1,bins_no))

species_params( natmmod)[c('d','M')] <- c( d_lm, mu0_lm)
species_params( natmmod2)[c('d','M')] <- c( d_DB, mu0_DB)
species_params( natmmod3)[c('d','M')] <- c( d_nls, mu0_nls)

natmmod <- MIZER( model = natmmod, catch = corLFD, compiler = compiler)
natmmod2 <- MIZER( model = natmmod2, catch = corLFD, compiler = compiler)
natmmod3 <- MIZER( model = natmmod3, catch = corLFD, compiler = compiler)

msp1 <- spf( natmmod, name = 'Linear Model')
msp2 <- spf( natmmod2, name = 'ad-hoc Model')
msp3 <- spf( natmmod3, name = 'NLS Model')

stotnm <- rbind( msp1, msp2, msp3) %>% filter( w > 10)

ggplot( stotnm, aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')

ggsave( paste0( tpnmdir,'spectra_comparison.jpg'), width = 9, height = 7)

ggplot( stotnm, aes( x = w, y = value2, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10( limits = c( 10, NA)) + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')

ggsave( paste0( tpnmdir,'spectra_comparison2.jpg'), width = 9, height = 7)


# Time periods --------------------------------

# Loop for time periods

tpdir <- paste0( getwd(), '/output/time_periods/')
dir.create( path = tpdir, showWarnings = TRUE, recursive = TRUE)

tp_mods <- list()

years <- list( years[[2]], years[[3]], years[[4]])
for(i in 1:length(years)) names(years)[[i]] <- paste0( years[[i]][1],' - ', years[[i]][5]) 

for( i in 1:length(years)){
  
  ny <- names(years)[i]
  vy <- years[[i]]
  ychar <- paste0( vy[1],'-',vy[length(vy)])
  
  ssbio_tp <- sum(assessment[,quantity][assessment$Year %in% vy]*1e6)/length(vy)
  
  itp_mod <- hake_mizer
  
  species_params(itp_mod)$biomass_observed <- ssbio_tp
  
  itp_mod <- itp_mod |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
  
  gear_params( itp_mod)$yield_observed <- corLFD_sum[[i]]$catch
  
  tp_mods[[ychar]] <- MIZER( model =  itp_mod, catch = corLFD_list[[i]], compiler = compiler,
                             plot = T, plot_dir = paste0(tpdir,ychar,'/'))
  
  plotSpectra( tp_mods[[ychar]]) + theme_bw()
  ggsave(paste0(tpdir,ychar,'/spectra.jpg'), width = 9, height = 7) 
  
  plotSpectra( tp_mods[[ychar]], power = 2) + theme_bw()
  ggsave(paste0(tpdir,ychar,'/spectra2.jpg'), width = 9, height = 7) 
  
}


tpdf <- NULL
for(i in names(tp_mods)){ tpdf <- rbind( tpdf, spf( tp_mods[[i]], name = i))}

ggplot( tpdf %>% filter(w>10), aes( x = w, y = value, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Period')

ggsave( paste0( tpdir, 'spectra_comparison.jpg'), width = 12, height = 7)
ggsave( paste0( tpdir, 'spectra_comparison_ppt.jpg'), width = 6, height = 4)

ggplot( tpdf %>% filter(w>10), aes( x = w, y = value2, color = Model)) + 
  geom_line( linewidth = .8) + scale_x_log10() + scale_y_log10() + theme_bw() +
  labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Period')

ggsave( paste0( tpdir, 'spectra_comparison2.jpg'), width = 12, height = 7)
ggsave( paste0( tpdir, 'spectra_comparison2_ppt.jpg'), width = 6, height = 4)


tp_table <- matrix( NA, nrow = length(tp_mods), ncol = length(all_pars), 
                    dimnames = list( names(tp_mods), all_pars))

for(i in names(tp_mods)){ tp_table[i,] <- parsf( tp_mods[[i]])}
tp_table



## NM and cannibalism for TP -------------

tp_cann1 <- list()
tp_cann2 <- list()
tpc_table1 <- NULL
tpc_table2 <- NULL

tp_nm1 <- list()
tp_nm2 <- list()
tp_nm3 <- list()
tpnm_table1 <- NULL
tpnm_table2 <- NULL
tpnm_table3 <- NULL

for(i in 1:length(tp_mods)){
  
  inm1 <- inm2 <- inm3 <- icann <- icann2 <- tp_mods[[i]]
  iname <- names(tp_mods)[[i]]
  
  ifiles <- list.files( list.files(tpdir, full.names = TRUE)[i], full.names = TRUE)
  
  # NLS power law fit
  
  mu0_lm <- NatM_pars['lm','mu0']; d_lm <- NatM_pars['lm','d']
  lm_mort <- mu0_lm*w(tp_mods[[i]])^(d_lm)
  
  mu0_DB <- NatM_pars['ad_hoc','mu0']; d_DB <- NatM_pars['ad_hoc','d']
  DB_mort <- mu0_DB*w(tp_mods[[i]])^(d_DB)
  
  mu0_nls <- NatM_pars['nls','mu0']; d_nls <- NatM_pars['nls','d']
  nls_mort <- mu0_nls*w(tp_mods[[i]])^(d_nls)
  
  ext_mort( inm1) <- array( lm_mort, dim=c(1,bins_no))
  ext_mort( inm2) <- array( DB_mort, dim=c(1,bins_no))
  ext_mort( inm3) <- array( nls_mort, dim=c(1,bins_no))
  
  species_params( inm1)[c('d','M')] <- c( d_lm, mu0_lm)
  species_params( inm2)[c('d','M')] <- c( d_DB, mu0_DB)
  species_params( inm3)[c('d','M')] <- c( d_nls, mu0_nls)
  
  vy <- years[[i]]
  
  pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% vy)])
  interaction_matrix( icann1)[] <- pcann
  interaction_matrix( icann2)[] <- pcann
  
  ext_mort( icann1) <- ext_mort( icann1) - getPredMort( icann1)
  
  tp_cann1[[iname]] <- icann1 <- MIZER( model = icann1, catch = corLFD_list[[i]], compiler = compiler, 
                                        plot = T, plot_dir = paste0(cdir,iname,'/minus/'))
  
  tp_cann2[[iname]] <- icann2 <- MIZER( model = icann2, catch = corLFD_list[[i]], compiler = compiler, 
                                        plot = T, plot_dir = paste0(cdir,iname,'/plus/'))
  
  tp_nm1[[iname]] <- inm1 <- MIZER( model = inm1, catch = corLFD_list[[i]], compiler = compiler, 
                                    plot = T, plot_dir = paste0(tpnmdir,iname,'/lm/'))
  
  tp_nm2[[iname]] <- inm2 <- MIZER( model = inm2, catch = corLFD_list[[i]], compiler = compiler, 
                                    plot = T, plot_dir = paste0(tpnmdir,iname,'/adhoc/'))
  
  tp_nm3[[iname]] <- inm3 <- MIZER( model = inm3, catch = corLFD_list[[i]], compiler = compiler, 
                                    plot = T, plot_dir = paste0(tpnmdir,iname,'/nls/'))
  
  file.copy( ifiles, paste0( tpnmdir,iname,'/'), recursive = TRUE)
  file.copy( ifiles, paste0( cdir,iname,'/'), recursive = TRUE)
  
  tpc_table1 <- rbind( tpc_table1, parsf(icann1))
  tpc_table2 <- rbind( tpc_table2, parsf(icann2))
  
  tpnm_table1 <- rbind( tpnm_table1, parsf(inm1))
  tpnm_table2 <- rbind( tpnm_table2, parsf(inm2))
  tpnm_table3 <- rbind( tpnm_table3, parsf(inm3))
  
  tpb <- spf( tp_mods[[i]], name = 'Base Model') 
  ctp1 <- spf( tp_cann1[[i]], name = "With Cannibalism")
  ctp2 <- spf( tp_cann2[[i]], name = "Cannibalism + NM")
  
  mtp1 <- spf( tp_nm1[[i]], name = "Linear Model")
  mtp2 <- spf( tp_nm2[[i]], name = "ad-hoc Model")
  mtp3 <- spf( tp_nm3[[i]], name = "NLS Model")
  
  cfdf <- rbind( tpb, ctp1, ctp2) %>% filter(w>10)
  mfdf <- rbind( mtp1, mtp2, mtp3) %>% filter(w>10)
  
  ggplot( cfdf, aes( x = w, y = value, color = Model)) + 
    geom_line( linewidth = .8) + scale_x_log10() + theme_bw() +
    labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')
  
  ggsave( paste0( cdir, iname, '/spectra_comparison.jpg'), width = 6, height = 4)
  
  ggplot( cfdf, aes( x = w, y = value2, color = Model)) + 
    geom_line( linewidth = .8) + scale_x_log10() + theme_bw() +
    labs ( x = 'Weigth [g]', y = 'Biomass density [g]', color = 'Model')
  
  ggsave( paste0( cdir,iname, '/spectra_comparison2.jpg'), width = 6, height = 4)
  
  ggplot( mfdf, aes( x = w, y = value, color = Model)) + 
    geom_line( linewidth = .8) + scale_x_log10() + theme_bw() +
    labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')
  
  ggsave( paste0( tpnmdir,iname, '/spectra_comparison.jpg'), width = 6, height = 4)
  
  ggplot( mfdf, aes( x = w, y = value2, color = Model)) + 
    geom_line( linewidth = .8) + scale_x_log10() + theme_bw() +
    labs ( x = 'Weigth [g]', y = 'Biomass density', color = 'Model')
  
  ggsave( paste0( tpnmdir,iname, '/spectra_comparison2.jpg'), width = 6, height = 4)
  
  
  cd0 <- plotDeath( tp_mods[[i]], return_data = T)
  cd1 <- plotDeath( tp_cann1[[i]], return_data = T)
  cd2 <- plotDeath( tp_cann2[[i]], return_data = T)
  
  md1 <- plotDeath( tp_nm1[[i]], return_data = T)
  md2 <- plotDeath( tp_nm2[[i]], return_data = T)
  md3 <- plotDeath( tp_nm3[[i]], return_data = T)
  
  cd0$Model = 'Base Model'
  cd1$Model = 'With Cannibalism'
  cd2$Model = 'Cannibalism + NM'
  
  md1$Model = 'Linear Model'
  md2$Model = 'ad-hoc Model'
  md3$Model = 'NLS Model'
  
  cippt <- rbind( cd0, cd1, cd2) %>% filter(w>1)
  mippt <- rbind( md1, md2, md3) %>% filter(w>1)
  
  cippt <- cippt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))
  mippt <- mippt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))
  
  cippt$Model <- factor( cippt$Model, levels = c( "Base Model","With Cannibalism","Cannibalism + NM"))
  
  cdeath <- ggplot( cippt, aes(x = w, y = value, fill = Cause)) +
    geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
    facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
    scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                       labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))
  
  cdeath
  
  ggsave( paste0( cdir,iname, '/death_comparison.jpg'), width = 6, height = 4)
  
  mdeath <- ggplot( mippt, aes(x = w, y = value, fill = Cause)) +
    geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
    facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
    scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                       labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))
  
  mdeath
  
  ggsave( paste0( tpnmdir,iname, '/death_comparison.jpg'), width = 6, height = 4)
  
  wv <- as.numeric(icann@w)
  
  tpb$Mort <- as.numeric(getMort(tp_mods[[i]])); tpb$Growth <- as.numeric(getEGrowth(tp_mods[[i]]))
  ctp1$Mort <- as.numeric(getMort(icann1)); ctp1$Growth <- as.numeric(getEGrowth(icann1))
  ctp2$Mort <- as.numeric(getMort(icann2)); ctp2$Growth <- as.numeric(getEGrowth(icann2))
  
  mtp1$Mort <- as.numeric(getMort(inm1)); mtp1$Growth <- as.numeric(getEGrowth(inm1))
  mtp2$Mort <- as.numeric(getMort(inm2)); mtp2$Growth <- as.numeric(getEGrowth(inm2))
  mtp3$Mort <- as.numeric(getMort(inm3)); mtp3$Growth <- as.numeric(getEGrowth(inm3))
  
  
  cffdf <- rbind( tpb, ctp1, ctp2) %>% filter(w>10)
  mffdf <- rbind( mtp1, mtp2, mtp3) %>% filter(w>10)
  
  cffdf$Model <- factor(cffdf$Model, levels = c("Base Model","With Cannibalism","Cannibalism + NM"))
  
  ggplot( cffdf, aes(x = w, y = Mort, color = Model)) +
    geom_line( linewidth = 0.8) + scale_x_log10() +   theme_bw() + 
    labs( x = "Weight [g]", y = "Mortality")
  
  ggsave( paste0( cdir,iname, '/mortality_comparison.jpg'), width = 6, height = 4)
  
  ggplot( cffdf, aes(x = w, y = Growth, color = Model)) +
    geom_line( linewidth = 0.8) + scale_x_log10() +   theme_bw() +
    labs( x = "Weight [g]")
  
  ggsave( paste0( cdir,iname, '/growth_comparison.jpg'), width = 6, height = 4)
  
  ggplot( mffdf, aes(x = w, y = Mort, color = Model)) +
    geom_line( linewidth = 0.8) + scale_x_log10() +   theme_bw() +
    labs( x = "Weight [g]", y = "Mortality")
  
  ggsave( paste0( tpnmdir,iname, '/mortality_comparison.jpg'), width = 6, height = 4)
  
  ggplot( mffdf, aes(x = w, y = Growth, color = Model)) +
    geom_line( linewidth = 0.8) + scale_x_log10() +   theme_bw() +
    labs( x = "Weight [g]")
  
  ggsave( paste0( tpnmdir,iname, '/growth_comparison.jpg'), width = 6, height = 4)
  
  
}



# Save all ----------------------

save.image( './output/alldata.RData')




