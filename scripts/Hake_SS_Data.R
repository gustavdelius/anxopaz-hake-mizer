######## WGBIE24 southern hake data from ICES ##########

rm( list=ls())

assessment <- icesSAG::getSAG( 'hke.27.8c9a', 2024)

ss_stk <- r4ss::SSgetoutput( dirvec = "./data/WGBIE24", getcovar = F, verbose = FALSE)[[1]]$timeseries
biomass <- ss_stk$Bio_all[ which(ss_stk$Yr%in%assessment$Year & ss_stk$Seas == 1)]

assessment <- assessment[,-c(14:ncol(assessment))]

for ( i in 1:nrow(assessment)) 
  assessment$catches[i] <- assessment$landings[i] + ifelse(is.na(assessment$discards[i]),0,assessment$discards[i])

assessment$biomass <- biomass

save( assessment, file = './input/Hake_SS_Data.RData')

