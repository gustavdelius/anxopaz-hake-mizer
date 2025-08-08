######## Mortality data from SS to MIZER ##########

rm(list=ls())

# Load libraries --------------------------------------------------------------

library(r4ss)
library(dplyr)
library(ggplot2)
library(tidyr)

source( './scripts/aux_functions.R')


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

M_female <- subset( M, Seas==1 & Settlement==1 & Sex==1)   # Females (Season 1 Settlement 1)
M_male <- subset( M, Seas==1 & Settlement==1 & Sex==2)   # Males (Season 1 Settlement 1)
NatM <- (M_female[-(1:4)] + M_male[-(1:4)])/2    # Mean of both sexes

ages <- as.numeric( names( NatM))
NatM <- as.numeric( NatM)

lengths <- alf( ages, Linf, Kvb, al0)
lengths[1] <- 4
names(lengths) <- names(NatM) <- ages


#- Eval. of straight lines between values at x(age)=0,1,5,15
len0_1 <- seq( floor(lengths['0']), floor(lengths['1']), by =0.1)
len1_5 <- ceiling(lengths['1']):floor(lengths['5'])
len5_15 <- ceiling(lengths['5']):ceiling(lengths['15'])

NatM0_1 <- linef( len0_1, lengths['0'], lengths['1'], NatM['0'], NatM['1'])
NatM1_5 <- linef( len1_5, lengths['1'], lengths['5'], NatM['1'], NatM['5'])
NatM5_15 <- linef( len5_15, lengths['5'], lengths['15'], NatM['5'], NatM['15'])

NatMt <- c( NatM0_1, NatM1_5, NatM5_15)
lvec <- c( len0_1, len1_5, len5_15)

M_ext <- data.frame( length = lvec, weight = a*(lvec)^b, NatM = NatMt)

M_ext %>% ggplot( aes( x = length, y = NatM)) +
  geom_line() + ylab("Mortality") + xlab("Length (cm)") + theme_bw()

M_ext %>% ggplot( aes( x = weight, y = NatM)) +
  geom_line() + ylab("Mortality") + xlab("Weigth (g)") + theme_bw()


### Power law fitting: mu(w)=mu0*w^d ----------

# Nonlinear least squares (NLS) fit

fit <- nls( NatM ~ I( mu0 * weight^d), data=M_ext, start=list(mu0=2,d=0))

fitlm <- lm( log10(NatM) ~ log10(weight), data = M_ext)

mu0 <- coef(fit)['mu0'] 
d <- coef(fit)['d']

mu0_lm <- 10^coef(fitlm)['(Intercept)']
d_lm <- coef(fitlm)['log10(weight)']

NatM <- M_ext %>% mutate( M_fit.db = .7629334*(weight)^(-0.1094125)) %>% 
  mutate( M_fit = mu0*(weight)^d) %>%
  mutate( M_fit.lm = mu0_lm*(weight)^d_lm)

NatM_pars <- matrix( c( .7629334, -0.1094125,
                        coef(fit)['mu0'], coef(fit)['d'],
                        10^coef(fitlm)['(Intercept)'], coef(fitlm)['log10(weight)']), 
                     byrow=T, 3, 2, dimnames = list( c('DB','nls','lm'), c('mu0','d')))

M_ext_plot <- NatM %>% pivot_longer( cols = NatM:M_fit.lm, names_to = "Model", values_to = "Mortality")

lnmplot <- M_ext_plot %>% ggplot( aes( x = length, y = Mortality, col = Model)) + 
  geom_line( ) + geom_point( data = subset( M_ext_plot, Model == 'NatM')) + 
  theme_bw() + labs( x = "Length (cm)", color= "Model")
lnmplot

wnmplot <- M_ext_plot %>% ggplot( aes( x = weight, y = Mortality, col = Model)) + 
  geom_line( ) + geom_point( data = subset( M_ext_plot, Model == 'NatM')) + 
  theme_bw() + labs( x = "Weight (g)", color= "Model")
wnmplot

lognmplot <- M_ext_plot %>% ggplot( aes( x = weight, y = Mortality, col = Model)) + 
  scale_x_log10() + scale_y_log10() +
  geom_line( ) + geom_point( data = subset( M_ext_plot, Model == 'NatM')) + 
  theme_bw() + labs( y = "log( mortality)", x = "log( weight (g))", color= "Model")
lognmplot


pdf("./plots/data/natural_mortality_fit.pdf", width = 10, height = 6, onefile = TRUE)
print(lnmplot)
print(wnmplot)
print(lognmplot)
dev.off()

save( NatM, NatM_pars, file = './data/Natural_Mortality.RData')


