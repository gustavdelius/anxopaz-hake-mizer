
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~   Hake's MIZER model with cannibalism  ~~~~~ #
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

load( './output/hake_model.RData')
load( './data/Diet.RData')

plotSpectra( hake_mizer, power = 2)

cannibal <- hake_mizer2 <- hake_mizer

ext_mort( hake_mizer2)[1:60] <- 0

# Turn on cannibalism ----------------------------

pcann <- mean( cannibal_byyear$Percentage[ which(cannibal_byyear$Year %in% aver_y)])

interaction_matrix( cannibal)[] <- pcann
interaction_matrix( hake_mizer2)[] <- pcann

getPredMort( cannibal) - getPredMort( hake_mizer2)

# Reduce external mortality because Hake predation is now modeled explicitly
# cannibal2 <- cannibal
ext_mort( cannibal) <- ext_mort( cannibal) - getPredMort( cannibal)

# Check that steady state has not changed
p_test <- steadySingleSpecies( cannibal)
all.equal( initialN(cannibal), initialN(p_test), tolerance = 1e-5)

plotSpectra( cannibal, power = 2, total = TRUE)
# plotSpectra( cannibal2, power = 2, total = TRUE)

cannibal <- setBevertonHolt( cannibal, reproduction_level = 0.01)
# cannibal2 <- setBevertonHolt( cannibal2, reproduction_level = 0.01)

plotYieldVsF( cannibal, species="Hake",F_max=15)
# plotYieldVsF( cannibal2, species="Hake",F_max=15)


## Fit ------------

cannibal_hake <- MIZER( model = cannibal, catch = corLFD)
# cannibal_hake2 <- MIZER( model = cannibal2, catch = corLFD)


# Comparison -----------------

plotDiet( hake_mizer)
plotDiet( cannibal_hake)
# plotDiet( cannibal_hake2)

p1 <- plotDeath( hake_mizer) + theme_bw(); p1
p2 <- plotDeath( cannibal_hake) + theme_bw(); p2
# p3 <- plotDeath( cannibal_hake2) + theme_bw(); p3

plot_lfd( hake_mizer, corLFD)
plot_lfd( cannibal_hake, corLFD)
# plot_lfd( cannibal_hake2, corLFD)

plot_lfd_gear( hake_mizer, corLFD)
plot_lfd_gear( cannibal_hake, corLFD)
# plot_lfd_gear( cannibal_hake2, corLFD)

ggpubr::ggarrange(p1, p2, nrow=1, common.legend = TRUE, legend="bottom")
# ggpubr::ggarrange(p1, p2, p3, nrow=1, common.legend = TRUE, legend="bottom")

# plotSpectra2(hake_mizer, name1 = "Base Model",
#              cannibal_hake2, name2 = "With cannibalism",
#              power = 2, resource = FALSE) + theme_bw()

plotSpectra2(hake_mizer, name1 = "Base Model",
             cannibal_hake, name2 = "With cannibalism",
             power = 2, resource = FALSE) + theme_bw()

plotSpectra2(cannibal, name1 = "No refit",
             hake_mizer, name2 = "Base Model",
             power = 2, resource = FALSE) + theme_bw()

plotSpectra2(cannibal, name1 = "No refit",
             cannibal_hake, name2 = "Cannibal refit",
             power = 2, resource = FALSE) + theme_bw()

getYield( hake_mizer)
getYield( cannibal_hake)
# getYield( cannibal_hake2)


pdf("./plots/cannibalism.pdf", width = 10, height = 6, onefile = TRUE)

ggpubr::ggarrange(p1, p2, nrow=1, common.legend = TRUE, legend="bottom")

plotSpectra2(hake_mizer, name1 = "Base Model",
             cannibal_hake, name2 = "With cannibalism",
             power = 2, resource = FALSE) + theme_bw()

plotSpectra2(cannibal, name1 = "No refit",
             cannibal_hake, name2 = "Cannibal refit",
             power = 2, resource = FALSE) + theme_bw()
dev.off()


dir_name <- paste0( getwd(), '/plots/cannibalism/', aver_y[1], '-', aver_y[length(aver_y)], '/')
dir.create( path = dir_name, showWarnings = TRUE, recursive = TRUE)

pspectra <- plotSpectra2(hake_mizer, name1 = "Base Model", cannibal_hake, name2 = "With cannibalism",
  power = 2, resource = FALSE, wlim = c(10, NA)) + theme_bw()

pspectra
ggsave( paste0(dir_name,'spectra.jpg'), width = 9, height = 7)

p1 <- plotDeath( hake_mizer, return_data = T)
p2 <- plotDeath( cannibal_hake, return_data = T)

p1$Model = 'Base Model'
p2$Model = 'With Cannibalism'

ppt <- rbind( p1, p2)
ppt <- ppt %>% mutate(Cause = factor( Cause, levels = c("External", "Fishing", "Hake"), labels = c("Natural", "Fishing", "Cannibalism")))

pdeath <- ggplot(ppt, aes(x = w, y = value, fill = Cause)) +
  geom_area(position = "stack") + scale_x_log10() +   theme_bw() +
  facet_wrap(~ Model) + labs( x = "Weight [g]", y = "Proportion of Death", fill = "Mortality Cause") +
  scale_fill_manual( values = c("Natural" = "#66c2a5", "Fishing" = "#fc8d62", "Cannibalism" = "#8da0cb"),
                     labels = c("Natural" = "Natural", "Fishing" = "Fishing", "Cannibalism" = "Cannibalism"))

pdeath
ggsave( paste0(dir_name,'death.jpg'), width = 9, height = 7)


l1 <- plot_lfd( hake_mizer, corLFD, return_df = T)
l2 <- plot_lfd( cannibal_hake, corLFD, return_df = T)

l1$Model <- 'Base Model'
l2$Model <- 'With Cannibalism'

ldf1 <- rbind( l1, subset(l2, Type == 'Estimated'))

ldfp1 <- ggplot(ldf1, aes(x = Length, y = Density)) +
  geom_bar(data = subset(ldf1, Type != 'Estimated'), aes( fill = Type), stat = "identity", position = "dodge", alpha = 0.6) +
  geom_line(data = subset(ldf1, Type == 'Estimated'), aes( color = Model), linewidth = 1) +
  scale_fill_manual(values = c("Observed" = "yellowgreen")) +
  scale_color_manual(values = c( "Base Model" = "#E41A1C", "With Cannibalism" = "#377EB8")) +
  theme_bw() + labs( x = "Size [cm]", y = "Normalised number density [1/cm]", fill = NULL, color = NULL)

ldfp1
ggsave( paste0(dir_name,'ldf.jpg'), width = 9, height = 7)

save.image( './output/cannibal_hake.RData')
