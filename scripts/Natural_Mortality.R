######## Mortality data from SS to MIZER ##########

rm(list=ls())

# Load libraries --------------------------------------------------------------

library(r4ss)
library(dplyr)
library(ggplot2)
library(tidyr)

source( './scripts/aux_functions.R')

dir.create( path = paste0( getwd(), '/plots/data/natural_mortality'), showWarnings = TRUE, recursive = TRUE)


# SS data --------------------------------------------------------------

ss_stk <- r4ss::SSgetoutput( dirvec = "./data/WGBIE24", getcovar = F, verbose = FALSE)[[1]]
sspars <- ss_stk$parameters
grpars <- ss_stk$Growth_Parameters


## Growth ---------------------

Linf_f <- grpars[1,'Linf']
Linf_m <- Linf_f*exp(sspars['L_at_Amax_Mal_GP_1','Value'])
Linf <- (Linf_f+Linf_m)/2

Kvb <- grpars$K[1]
latamin <- grpars$L_a_A1[1]
al0_f <- grpars$A_a_L0[1]
al0_m <- al0_f*exp(0.6)   # from SS control file

al0 <- (al0_f+al0_m)/2

a <- sspars['Wtlen_1_Fem_GP_1','Value']*1000   # kg to g
b <- sspars['Wtlen_2_Fem_GP_1','Value'] 


# Plot check

r4ss::SSplotBiology( ss_stk, subplots = 1)
lines( 0:15, ss_stk$endgrowth$Len_Beg[1:16], col='black')
lines( seq( 0, 15, length.out=500), alf( seq( 0, 15, length.out=500), Linf_f, Kvb, al0_f), col = 'green')
lines( laf( seq( 0, Linf_f, length.out=500), Linf_f, Kvb, al0_f), seq( 0, Linf_f, length.out=500), col='darkgreen', lty = 'dashed')
lines( laf( seq( 0, Linf_m, length.out=500), Linf_m, Kvb, al0_m), seq( 0, Linf_m, length.out=500), col='black', lty = 'dashed')
lines( laf( seq( 0, Linf, length.out=500), Linf, Kvb, al0), seq( 0, Linf, length.out=500), col='black')


## Natural Mortality ----------------------

M <- ss_stk$Natural_Mortality_Bmark

M_female <- subset( M, Seas==1 & Settlement==1 & Sex==1)[-(1:4)]   # Females (Season 1 Settlement 1)
M_male <- subset( M, Seas==1 & Settlement==1 & Sex==2)[-(1:4)]   # Males (Season 1 Settlement 1)
NatM <- ( M_female + M_male)/2    # Mean of both sexes

ages <- as.numeric( names( NatM))
NatM <- as.numeric( NatM)

lengths <- alf( ages, Linf, Kvb, al0)
lengths[1] <- 4

weights <- lwf( lengths,a,b)

names(lengths) <- names(weights) <- names(NatM) <- ages


mdf <- rbind( 
  data.frame( Ages = ages, Length = lengths, Weight = weights, Sex = 'Female', Mortality = as.numeric(M_female)),
  data.frame( Ages = ages, Length = lengths, Weight = weights, Sex = 'Male', Mortality = as.numeric(M_male)))

ggplot( mdf, aes( x = Ages, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()

ggsave( "./plots/data/natural_mortality/natural_mortality_SS.jpg", width = 5, height = 4)


ggplot( mdf, aes( x = Length, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()

ggplot( mdf, aes( x = Weight, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()


M_ext <- data.frame( Ages = ages, Length = lengths, Weight = weights, Sex = 'Average', Mortality = as.numeric(NatM))

fmdf <- rbind( mdf, M_ext)

ggplot( fmdf, aes( x = Ages, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()

ggsave( "./plots/data/natural_mortality/natural_mortality_SS_aver.jpg", width = 5, height = 4)

ggplot( fmdf, aes( x = Length, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()

ggplot( fmdf, aes( x = Weight, y = Mortality, color = Sex)) + theme_bw() + 
  geom_line() + geom_point()



# #- Eval. of straight lines between values at x(age)=0,1,5,15
# len0_1 <- seq( floor(lengths['0']), floor(lengths['1']), by =0.1)
# len1_5 <- ceiling(lengths['1']):floor(lengths['5'])
# len5_15 <- ceiling(lengths['5']):ceiling(lengths['15'])
# 
# NatM0_1 <- linef( len0_1, lengths['0'], lengths['1'], NatM['0'], NatM['1'])
# NatM1_5 <- linef( len1_5, lengths['1'], lengths['5'], NatM['1'], NatM['5'])
# NatM5_15 <- linef( len5_15, lengths['5'], lengths['15'], NatM['5'], NatM['15'])
# 
# NatMt <- c( NatM0_1, NatM1_5, NatM5_15)
# lvec <- c( len0_1, len1_5, len5_15)
# 
# M_ext <- data.frame( length = lvec, weight = a*(lvec)^b, NatM = NatMt)

M_ext %>% ggplot( aes( x = Length, y = Mortality)) +
  geom_line() + ylab("Mortality") + xlab("Length (cm)") + theme_bw()

ggsave( "./plots/data/natural_mortality/natural_mortality_SS_l.jpg", width = 4, height = 2)

M_ext %>% ggplot( aes( x = Weight, y = Mortality)) +
  geom_line() + ylab("Mortality") + xlab("Weigth (g)") + theme_bw()

ggsave( "./plots/data/natural_mortality/natural_mortality_SS_w.jpg", width = 4, height = 2)


### Power law fitting: mu(w)=mu0*w^d ----------

fit <- nls( Mortality ~ I( mu0 * Weight^d), data = M_ext, start = list( mu0=1.5, d=-0.05))

fitlm <- lm( log10(Mortality) ~ log10(Weight), data = M_ext)

lmdf <- data.frame( NatM = NatM, weight = weights)

mu0_nls <- coef(fit)['mu0'] 
d_nls <- coef(fit)['d']

mu0_lm <- 10^coef(fitlm)['(Intercept)']
d_lm <- coef(fitlm)['log10(Weight)']

NatM <- M_ext %>% mutate( M_fit.db = .7629334*(weights)^(-0.1094125)) %>% 
  mutate( M_fit = mu0_nls*(weights)^d_nls) %>%
  mutate( M_fit.lm = mu0_lm*(weights)^d_lm)

NatM_pars <- matrix( c( .7629334, -0.1094125, mu0_nls, d_nls, mu0_lm, d_lm),
                     byrow=T, 3, 2, dimnames = list( c('ad_hoc','nls','lm'), c('mu0','d')))

NatM_pars

M_ext_plot <- NatM %>% pivot_longer( cols = Mortality:M_fit.lm, names_to = "Model", values_to = "Mortality")

M_ext_plot <- M_ext_plot %>% mutate( Model = case_when(
  Model == "Mortality" ~ "SS model (input)", Model == "M_fit" ~ "Non-parametric model",
  Model == "M_fit.lm" ~ "Linear model", Model == "M_fit.db" ~ "ad hoc model", TRUE ~ Model))

mcols = c( "SS model (input)" = "black", "Non-parametric model" = "#1b9e77", "Linear model" = "#d95f02", "ad hoc model" = "#7570b3")
mlines <- c( "SS model (input)" = 1.2, "Non-parametric model" = 0.8, "Linear model" = 0.8, "ad hoc model" = 0.8)

lnmplot <- M_ext_plot %>% ggplot( aes( x = Length, y = Mortality, col = Model, linewidth = Model)) + 
  geom_line() + theme_bw() + labs( x = "Length (cm)", color= "Model") +
  scale_color_manual( values = mcols) + scale_linewidth_manual( values = mlines)
lnmplot

ggsave( "./plots/data/natural_mortality/natural_mortality_fit_l.jpg", width = 6, height = 4)

wnmplot <- M_ext_plot %>% ggplot( aes( x = Weight, y = Mortality, col = Model, linewidth = Model)) + 
  geom_line() + theme_bw() + labs( x = "Weight (g)", color= "Model") +
  scale_color_manual( values = mcols) + scale_linewidth_manual( values = mlines)
wnmplot

ggsave( "./plots/data/natural_mortality/natural_mortality_fit_w.jpg", width = 6, height = 3)

lognmplot <- M_ext_plot %>% ggplot( aes( x = Weight, y = Mortality, col = Model, linewidth = Model)) + 
  scale_x_log10() + scale_y_log10() + geom_line() + theme_bw() + 
  labs( y = "log-Mortality", x = "log-Weight (g)", color= "Model") +
  scale_color_manual( values = mcols) + scale_linewidth_manual( values = mlines)
lognmplot

ggsave( "./plots/data/natural_mortality/natural_mortality_fit_logw.jpg", width = 6, height = 3)

pdf("./plots/data/natural_mortality/natural_mortality_fit.pdf", width = 10, height = 6, onefile = TRUE)
print(lnmplot)
print(wnmplot)
print(lognmplot)
dev.off()


save( NatM, NatM_pars, file = './data/Natural_Mortality.RData')


