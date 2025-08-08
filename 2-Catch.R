####    CATCH DATA         ####
####  for southern Hake    ####

rm(list=ls())

library(dplyr)
library(ggplot2)
library(plotly)
library(patchwork)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(mizerMR)
library(tidyr)

source( './scripts/aux_functions.R')

plotdir <- paste0( './plots/data/catch/')
dir.create( path = plotdir, showWarnings = TRUE, recursive = TRUE)


# Biological parameters ---------------------
# 'Bio_Pars.R' script

load( './input/bio_pars.RData')


# Catch data ----------------------

# catch1 <- read.csv("./data/catch/catch 1948-1981.csv", header=T)
# catch2 <- read.csv("./data/catch/catch 1982-1993.csv", header=T)
# catch3 <- read.csv("./data/catch/catch 1994-2023.csv", header=T)
# Index <- read.csv("./data/catch/indices 1982-2023.csv", header=T)

catch <- read.csv("./data/catch/catch 1994-2023.csv", header=T)
LFD3 <- read.csv("./data/catch/LFDs 1994-2023.csv", header=T)

LFD3 <- LFD3[,-(4:5)]
names(LFD3)[4] <- "length"
LFD3$length <- as.numeric(gsub("len","",LFD3$length))
LFD3$fleet <- gsub("_fem","",LFD3$fleet)
LFD3$fleet <- gsub("_mal","",LFD3$fleet)
LFD3$fleet <- gsub("_ind","",LFD3$fleet)

LFD3$weight <- lwf(LFD3$length, a, b)


# Years for average -------------

aver_y <- 2014:2023

years <- list( aver_y = aver_y,
               tp1 = 1996:2000,
               tp2 = 2008:2012,
               tp3 = 2019:2023,
               aver_y2 = aver_y-20)


# Ordered gears --------------------

gear_names <- c("Art","ptArt","bakka","cdTrw","ptTrw","pair","palangre","vol","disc","PtSurv","SpSurv","cdSurv")
surv <- c( "PtSurv", "SpSurv", "cdSurv")


# Loop ---------------------

catch_list <- catch_plot <- LFD_list <- LFDc_list <- LFD_plot <- LFD_sum <- catch_comp <- catch_table <- 
  corLFD_list <- corLFDc_list <- corLFD_plot <- corLFD_sum <- LFD_list2 <- LFD_plot2 <- LFD_sum2 <- LFDc_list2 <- 
  list()


for( i in 1:length(years)){
  
  iy <- years[[i]]
  cy <- names(years)[i]
  
  cat(cy,'=', iy[1], '-',iy[length(iy)],'\n' )
  
  ## Catch -------------------
  
  icatch <- catch %>% filter(year %in% iy)
  
  icatch <- icatch %>% mutate( fleet = recode( fleet, "baka" = "bakka", "volanta" = "vol", "pairTrw" = "pair")) %>% 
    arrange( factor( fleet, levels = gear_names)) 
  
  icatch_summary <- icatch %>% group_by(year, fleet) %>%
    summarise(total_amount = sum(amount, na.rm = TRUE), .groups = "drop")
  
  catch_list[[cy]] <- icatch_summary %>% group_by(fleet) %>%
    summarise(catch = mean(total_amount, na.rm = TRUE)*1e3, .groups = "drop")
  
  catch_plot[[cy]] <- catch_list[[cy]] %>% ggplot( aes(x = fleet, y = catch/1e+6, fill = fleet)) +
    geom_bar(stat = "identity", position = "stack") + guides(fill="none") +
    labs( title = 'Catch by fleet', y='average yield (t)') + theme_bw(); catch_plot[[cy]]
  
  
  
  ## LFD -------------
  
  iLFD <- LFD3 %>% filter(year %in% iy, !fleet %in% surv)
  
  iLFD_summary <- iLFD %>% group_by(year, fleet, length, weight) %>%
    summarise(total_number = sum(number, na.rm = TRUE), .groups = "drop")
  
  LFD_list[[cy]] <- iLFD_summary %>% group_by(fleet, length, weight) %>%
    summarise( number = mean(total_number, na.rm = TRUE), .groups = "drop")
  
  LFDc_list[[cy]] <- LFD_list[[cy]] %>% filter(fleet %in% gear_names) %>% pivot_wider( names_from = fleet, values_from = number)
  
  LFD_plot[[cy]] <- LFD_list[[cy]] %>% ggplot( aes( x = length, color = fleet)) + 
    stat_density( aes( weight = number), adjust = 0.2, geom = "line", position = "identity") +
    theme_bw() + scale_x_continuous( n.breaks = 20) + 
    labs(title="Length distribution",x="length (cm)"); LFD_plot[[cy]]
  
  LFD_sum[[cy]] <- LFD_list[[cy]] %>% group_by(fleet) %>%
    summarise(catch = sum(weight * number, na.rm = TRUE), .groups = "drop") %>% 
    filter( !fleet %in% surv)
  
  
  ## Plots ----------
  
  chn <- paste0( iy[1], '-',iy[length(iy)])
  
  pdf(paste0( plotdir, "Catch-LDs_",chn,".pdf"), width = 10, height = 6, onefile = TRUE)
  print(catch_plot[[cy]])
  print(LFD_plot[[cy]])
  dev.off()
  
  
  ## Check C vs LFD -----------------
  
  catch_comp[[cy]] <- data.frame( fleet = catch_list[[cy]]$fleet, catch_t = catch_list[[cy]]$catch/1e6, 
    catch_lfd = LFD_sum[[cy]]$catch/1e6, diff = catch_list[[cy]]$catch/LFD_sum[[cy]]$catch)
  
  catch_table[[cy]] <- rbind( catch_comp[[cy]], data.frame( fleet = 'TOTAL (t)', catch_t = sum(catch_list[[cy]]$catch)/1e6, 
    catch_lfd = sum(LFD_sum[[cy]]$catch)/1e6, diff = sum(catch_list[[cy]]$catch)/sum(LFD_sum[[cy]]$catch)))
  
  print(catch_table[[cy]])
  
  
  ## Corrected LFD (N_LFD = N_Catch) ---------------------------
  
  iLFD <- LFD3 %>% filter(year %in% iy, !fleet %in% surv) %>%
    left_join(catch_comp[[cy]] %>% select(fleet, diff), by = "fleet") %>%
    mutate(number = number * diff) %>%
    select(-diff)
  
  iLFD_summary <- iLFD %>% group_by(year, fleet, length, weight) %>%
    summarise(total_number = sum(number, na.rm = TRUE), .groups = "drop")
  
  corLFD_list[[cy]] <- iLFD_summary %>% group_by(fleet, length, weight) %>%
    summarise( number = mean(total_number, na.rm = TRUE), .groups = "drop")
  
  corLFDc_list[[cy]] <- corLFD_list[[cy]] %>% 
    filter(fleet %in% gear_names) %>% pivot_wider( names_from = fleet, values_from = number)
  
  corLFD_plot[[cy]] <- corLFD_list[[cy]] %>% ggplot( aes( x = length, color = fleet)) + 
    stat_density( aes( weight = number), adjust = 0.2, geom = "line", position = "identity") +
    theme_bw() + scale_x_continuous( n.breaks = 20) + 
    labs(title="Length distribution",x="length (cm)"); LFD_plot[[cy]]
  
  corLFD_sum[[cy]] <- corLFD_list[[cy]] %>% group_by(fleet) %>%
    summarise(catch = sum(weight * number, na.rm = TRUE), .groups = "drop") %>% 
    filter( !fleet %in% surv)
  
  
  ## Only surveys -------------
  
  iLFD2 <- LFD3 %>% filter(year %in% iy, fleet %in% surv)
  
  iLFD_summary2 <- iLFD2 %>% group_by(year, fleet, length, weight) %>%
    summarise(total_number = sum(number, na.rm = TRUE), .groups = "drop")
  
  LFD_list2[[cy]] <- iLFD_summary2 %>% group_by(fleet, length, weight) %>%
    summarise( number = mean(total_number, na.rm = TRUE), .groups = "drop")
  
  LFDc_list2[[cy]] <- LFD_list2[[cy]] %>% filter(fleet %in% gear_names) %>% pivot_wider( names_from = fleet, values_from = number)
  
  LFD_plot2[[cy]] <- LFD_list2[[cy]] %>% ggplot( aes( x = length, color = fleet)) + 
    stat_density( aes( weight = number), adjust = 0.2, geom = "line", position = "identity") +
    theme_bw() + scale_x_continuous( n.breaks = 20) + 
    labs(title="Length distribution",x="length (cm)"); LFD_plot2[[cy]]
  
  LFD_sum2[[cy]] <- LFD_list2[[cy]] %>% group_by(fleet) %>%
    summarise(catch = sum(weight * number, na.rm = TRUE), .groups = "drop") %>% 
    filter( fleet %in% surv)
  
  pdf(paste0( plotdir, "Surv-Catch-LDs_",chn,".pdf"), width = 10, height = 6, onefile = TRUE)
  print(LFD_plot2[[cy]])
  dev.off()
  
  cat('\n \n')
  
}



# Main Catch & LFD -----------

Catch <- catch_list[['aver_y']]
LFD <- LFD_list[['aver_y']]
LFDc <- LFDc_list[['aver_y']]
LFDs <- LFD_sum[['aver_y']]
corLFD <- corLFD_list[['aver_y']]
corLFDc <- corLFDc_list[['aver_y']]
corLFDs <- corLFD_sum[['aver_y']]


pdf(paste0( plotdir, "Catch.pdf"), width = 10, height = 6, onefile = TRUE)
print(catch_plot[['aver_y']])
dev.off()

pdf(paste0( plotdir, "LFD.pdf"), width = 10, height = 6, onefile = TRUE)
print(LFD_plot[['aver_y']])
dev.off()

pdf(paste0( plotdir, "Surveys-LFD.pdf"), width = 10, height = 6, onefile = TRUE)
print(LFD_plot2[['aver_y']])
dev.off()


data.frame( fleet = Catch$fleet, catch_t = Catch$catch/1e6, 
            catch_lfd = corLFDs$catch/1e6, diff = Catch$catch/corLFDs$catch)
  

# Save --------------------

save( catch_list, Catch, LFD_list, LFDc_list, LFDs, LFD_sum, LFD, LFDc, catch_comp, catch_table, 
      aver_y, years, corLFD_list, corLFDc_list, corLFDs, corLFD_sum, corLFD, corLFDc, 
      file = './input/Catch.RData')

