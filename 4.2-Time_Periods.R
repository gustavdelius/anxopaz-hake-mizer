
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model by time periods  ~~~~~~ #
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

# load( './output/hake_model.RData')   # Results from "3-MIZER.R"
load( './output/cannibal_hake.RData')   # Results from "4.1-Cannibalism.R"


# Loop for time periods

tp_mods <- list()

years

for( i in c(5)){
  
  ny <- names(years)[i]
  vy <- years[[i]]
  
  ssbio_tp <- sum(assessment$SSB[assessment$Year %in% vy]*1e6)/length(vy)
  
  tp_mods[[ny]] <- hake_mizer
  
  species_params(tp_mods[[ny]])$biomass_observed <- ssbio_tp
  
  tp_mods[[ny]] <- tp_mods[[ny]] |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> steady()
  
  gear_params( tp_mods[[ny]])$yield_observed <- corLFD_sum[[ny]]$catch   # == Catch$catch; != LFDs$catch
  
  tp_mods[[ny]] <- MIZER( model =  tp_mods[[ny]], catch = corLFD_list[[ny]])
  tp_mods[[ny]] <- MIZER( model =  tp_mods[[ny]], catch = corLFD_list[[ny]])
  tp_mods[[ny]] <- MIZER( model =  tp_mods[[ny]], catch = corLFD_list[[ny]])
  
  plot_lfd_gear( tp_mods[[ny]], corLFD_list[[ny]])
  
}


hake_mizer2 <- tp_mods$aver_y2
vy2 <- years$aver_y2
lfd2 <- corLFD_list$aver_y2

plot_lfd( hake_mizer2, lfd2)
plot_lfd_gear( hake_mizer2, lfd2)

getYield( hake_mizer2); sum( hake_mizer2@gear_params$yield_observed)
getBiomass( hake_mizer2); hake_mizer2@species_params$biomass_observed


# Cannibalism ----------------------

cannibal2 <- hake_mizer2

pcann2 <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% vy2)])
interaction_matrix( cannibal2)[] <- pcann2

ext_mort( cannibal2) <- ext_mort( cannibal2) - getPredMort( cannibal2)

cannibal_hake2 <- MIZER( model = cannibal2, catch = lfd2)
cannibal_hake2 <- MIZER( model = cannibal_hake2, catch = lfd2)


# Comparison -----------------

plotDiet( hake_mizer)
plotDiet( cannibal_hake)
plotDiet( hake_mizer2)
plotDiet( cannibal_hake2)

plotDeath( hake_mizer) + theme_bw()
plotDeath( cannibal_hake) + theme_bw()
plotDeath( hake_mizer2) + theme_bw()
plotDeath( cannibal_hake2) + theme_bw()

p3 <- plotDeath( hake_mizer2, return_data = T)
p4 <- plotDeath( cannibal_hake2, return_data = T)

plot_lfd( hake_mizer, corLFD)
plot_lfd( cannibal_hake, corLFD)
plot_lfd( hake_mizer2, lfd2)
plot_lfd( cannibal_hake2, lfd2)

l3 <- plot_lfd( hake_mizer2, lfd2, return_df = T)
l4 <- plot_lfd( cannibal_hake2, lfd2, return_df = T)

plot_lfd_gear( hake_mizer, corLFD)
plot_lfd_gear( cannibal_hake, corLFD)
plot_lfd_gear( hake_mizer2, lfd2)
plot_lfd_gear( cannibal_hake2, lfd2)



dir_name2 <- paste0( getwd(), '/plots/cannibalism/', vy2[1], '-', vy2[length(vy2)], '/')
dir.create( path = dir_name2, showWarnings = TRUE, recursive = TRUE)

pspectra2 <- plotSpectra2(hake_mizer2, name1 = "No cannibalism", cannibal_hake2, name2 = "With cannibalism",
                         power = 2, resource = FALSE, wlim = c(10, NA)) + theme_bw()

pspectra2
ggsave( paste0(dir_name2,'spectra.jpg'), width = 9, height = 7)
ggsave( paste0(dir_name2,'spectra_poster.jpg'), width = 5, height = 3)

plotSpectra2( hake_mizer2, name1 = "No cannibalism", cannibal_hake2, name2 = "With cannibalism",
             power = 2, resource = FALSE, wlim = c(10, NA)) + 
  scale_y_log10("Biomass density [g]", limits = c( 1000, 10000000000)) + theme_bw()

ggsave( paste0(dir_name2,'spectra_poster2.jpg'), width = 5, height = 3)


p3$Model = 'No cannibalism'
p4$Model = 'With Cannibalism'

ppt2 <- rbind( p3, p4)
ppt2 <- ppt2 %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))

pdeath2 <- ggplot(ppt2, aes(x = w, y = value, fill = Cause)) +
  geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
  facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
  scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                     labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))

pdeath2
ggsave( paste0(dir_name2,'death.jpg'), width = 9, height = 7)
ggsave( paste0(dir_name2,'death_poster.jpg'), width = 5, height = 3)


l3$Model <- 'No cannibalism'
l4$Model <- 'With Cannibalism'

ldf2 <- rbind( l3, subset(l4, Type == 'Estimated'))

ldfp2 <- ggplot(ldf2, aes(x = Length, y = Density)) +
  geom_bar(data = subset(ldf2, Type != 'Estimated'), aes( fill = Type), stat = "identity", position = "dodge", alpha = 0.6) +
  geom_line(data = subset(ldf2, Type == 'Estimated'), aes( color = Model), linewidth = 1) +
  scale_fill_manual(values = c("Observed" = "yellowgreen")) +
  scale_color_manual(values = c( "No cannibalism" = "#E41A1C", "With Cannibalism" = "#377EB8")) +
  theme_bw() + labs( x = "Size [cm]", y = "Normalised number density [1/cm]", fill = NULL, color = NULL)

ldfp2
ggsave( paste0(dir_name2,'ldf.jpg'), width = 9, height = 7)
ggsave( paste0(dir_name2,'ldf_poster.jpg'), width = 5, height = 3)


save.image( './output/tp_hake.RData')
