#### BIOLOGICAL PARAMETERS ####
####  for southern Hake    ####

rm(list=ls())

library(dplyr)
library(ggplot2)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(mizerMR)

source( './scripts/aux_functions.R')


# SS WGBIE24 hake data --------------------------------

replist <- r4ss::SSgetoutput( dirvec = "./data/WGBIE24", getcovar = F, verbose = FALSE)[[1]]

sspars <- replist$parameters

grpars <- replist$Growth_Parameters     ## no males, then:
grpars <- readLines("./data/WGBIE24/Report.sso")
pos <- grep("^Growth_Parameters report:45", grpars)
grpars <- grpars[(pos + 1):(pos + 4 + 1)]
grpars <- read.table(text = paste(grpars, collapse = "\n"), header = TRUE)


# Model parameters ----------------------------

## Growth -----------------

a <- grpars$WtLen1[1] * 1e3; a   # 0.00377 (kg to g)
b <- grpars$WtLen2[1]; b        # 3.168 

Kvb <- grpars$K[1]; Kvb

Linf_f <- grpars$Linf[1]; Linf_f
Linf_m <- grpars$Linf[3]; Linf_m     # Linf_f*exp(sspars['L_at_Amax_Mal_GP_1','Value'])
Linf <- (Linf_f+Linf_m)/2; Linf

al0_f <- grpars$A_a_L0[1]; al0_f
al0_m <- grpars$A_a_L0[3]; al0_m     # from SS report file (report:45)
al0 <- (al0_f+al0_m)/2; al0

r4ss::SSplotBiology(replist, subplots = 1)
lines( 0:15, alf(0:15, Linf_f, Kvb, al0_f), col='black')
lines( 0:15, alf(0:15, Linf_m, Kvb, al0_m), col='black')
lines( 0:15, alf(0:15, Linf, Kvb, al0), col='black')

bins_no <- as.numeric(colnames(replist$natlen)[ncol(replist$natlen)]) - 1  # ceiling(Linf)


## Maturity -------------------

load( './data/Maturity_Size.RData')   # './scripts/Maturity.R' results

L50_f <- MatSize['L50','Females']; L50_f    # sspars['Mat50%_Fem_GP_1','Value'] 
L50_m <- MatSize['L50','Males']; L50_m    # aprox L50_m=L50_f*exp(sspars['L_at_Amax_Mal_GP_1','Value'])
L50 <- (L50_m+L50_f)/2; L50

L25_f <- MatSize['L25','Females']; L25_f
L25_m <- MatSize['L25','Males']; L25_m
L25 <- (L25_m+L25_f)/2; L25

a50 <- laf( L50, Linf, Kvb, al0)
a25 <- laf( L25, Linf, Kvb, al0)

w50 <- lwf( L50, a, b)
w25 <- lwf( L25, a, b)

kmat_f <- MatSize['k','Females']
kmat_m <- MatSize['k','Males']
kmat <- (kmat_m+kmat_f)/2; kmat

U <- kmat*L50/b

w_max <- lwf( bins_no, a, b)


## Predation -----------------------------

# beta: preferred predator-press mass ratio (PPMR); beta = 2^b = 8.99 (100 defalult)

load( './data/Predation.RData')

beta <- beta
sigma <- sigma


## Mizer pars -------------

h <- 4.75 * Kvb * Linf^0.75; h   # max. consumption rate; h(w) = h_1*w^n; 
# h_1 = 4.75*Kvb*linf^n; n = 0.75 (def)

bio_pars <- newSingleSpeciesParams(
  species_name="Hake", no_w=bins_no, w_max=w_max, w_mat=w50, lambda=2, h=h, beta=beta, sigma=sigma)


# lambda: exponent of the spectrum's power law; N(w) = kappa * w^(-lambda); expected to be around 2

species_params(bio_pars)$w_mat25 <- w25
species_params(bio_pars)$U <- U
species_params(bio_pars)$a <- a
species_params(bio_pars)$b <- b
species_params(bio_pars)$age_mat <- a50


## Natural Mortality -----------------------------

load( './data/Natural_Mortality.RData')   # './scripts/Natural_Mortality.R' results
# NLS power law fit

mu0 <- NatM_pars['nls','mu0']
d <- NatM_pars['nls','d']
mort <- mu0*w(bio_pars)^(d)    

ext_mort(bio_pars) <- array( mort, dim=c(1,bins_no))

species_params(bio_pars)$d <- d
species_params(bio_pars)$M <- mu0


# Save ----------------

rm_functions <- function( envir = globalenv()) {
  objs <- ls( envir = envir)
  funs <- objs[sapply(objs, function(x) is.function(get(x, envir = envir)))]
  rm(list = funs, envir = envir)
}

rm_functions()

save.image( './input/Bio_Pars.RData')

