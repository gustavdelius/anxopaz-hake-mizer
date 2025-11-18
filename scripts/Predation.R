######## Hake's Predation fit for MIZER ##########

rm(list=ls())

library(dplyr)
library(ggplot2)
library(plotly)
library(reshape)
library(sm)
library(mizer)
library(mizerExperimental)
library(mizerMR)


# SS WGBIE24 hake data --------------------------------

replist <- r4ss::SSgetoutput( dirvec = "./data/WGBIE24", getcovar = F, verbose = FALSE)[[1]]
sspars <- replist$parameters
grpars <- replist$Growth_Parameters


# Model parameters ----------------------------

## Growth -----------------

a <- sspars['Wtlen_1_Fem_GP_1','Value']*1000; a   # 0.00377 (kg to g)
b <- sspars['Wtlen_2_Fem_GP_1','Value']; b        # 3.168 


# Predation -------------------------------

Lpred <- 0:90

alpha <- 1.775
beta <- 0.438

adif <- 1.725
bdif <- 0.082


## L(prey) vs L(pred) -------------------

Lprey <- alpha + beta*(Lpred+0.5)     # lm from './data/Cannibalism/South Hake GADGET Cannibalism for CS meeting.ppt'
Lpreyu <- (alpha+adif) + (beta+bdif)*(Lpred+0.5)
Lpreyl <- (alpha-adif) + (beta-bdif)*(Lpred+0.5)

plot( Lpred, Lprey, type='l', ylim=c(0,40), xlim=c(5,85), xaxt = "n", yaxt = "n",
      xlab = 'Predator length (cm)', ylab = 'Prey length (cm)')
axis( 1, at = seq(5,85,by=10))
axis( 2, at = seq(0,40,by=5))
lines( Lpred, Lpreyu, lty=2)
lines( Lpred, Lpreyl, lty=2)
abline( h = seq(0,40,by=5), col='lightgrey', lty=3, lwd=0.2)


## W(prey) vs W(pred) -------------------
## Considering W=a+L^b

Wprey <- c( 0:400)
Wpred <- a * (((Wprey/a)^(1/b) - alpha - 0.5*beta)/ beta)^b
Wpredu <- a * (((Wprey/a)^(1/b) - (alpha+adif) - 0.5*(beta+bdif))/ (beta+bdif))^b
Wpredl <- a * (((Wprey/a)^(1/b) - (alpha-adif) - 0.5*(beta-bdif))/ (beta-bdif))^b

plot( Wpred, Wprey, type='l', xaxt = "n", yaxt = "n",
      xlab = 'Predator weigth (g)', ylab = 'Prey weigth (g)')
axis(1, at = seq(0,4500,by=500))
axis(2, at = seq(0,400,by=50))
lines( Wpredu, Wprey, lty=2)
lines( Wpredl, Wprey, lty=2)
abline( h = seq(0,400,by=50), col='lightgrey', lty=3, lwd=0.2)


# Fit ---------------------------

preypred <- c( mean = as.numeric( 1/coef( lm(Wprey~Wpred-1))), 
               lower = as.numeric( 1/coef( lm( Wprey~Wpredl-1))), 
               upper = as.numeric( 1/coef( lm(Wprey~Wpredu-1))))

preypred

sigma <- (preypred['upper']/preypred['mean'] +  preypred['mean']/preypred['lower'])/2
beta <- preypred['mean']

betaupp <- beta/sigma
betalow <- beta*sigma


## Check ---------------------------

predsize <- 538
predsize2 <- 1030

intercepts <- c( (1/beta)*predsize, (1/betaupp)*predsize, (1/betalow)*predsize)
intercepts2 <- c( (1/beta)*predsize2, (1/betaupp)*predsize2, (1/betalow)*predsize2)


lines( (1/beta)*1:4500, col='red')
lines( (1/betaupp)*1:4500, col='red', lty=2)
lines( (1/betalow)*1:4500, col='red', lty=2)

abline( h = intercepts, col='blue', lty=c(1,2,2), lwd=0.2)
abline( v = predsize, col='blue', lwd=0.2)

abline( h = intercepts2, col='green', lty=c(1,2,2), lwd=0.2)
abline( v = predsize2, col='green', lwd=0.2)


# Save -------------------------------------------------

save( beta, sigma, file = './data/Predation.RData')



# After reading pars in MIZER

load( './input/Bio_Pars.RData')

species_params(bio_pars)$beta <- beta
species_params(bio_pars)$sigma <- sigma

pred_kernel <- melt( getPredKernel(bio_pars)[, as.character(predsize), , drop = FALSE])
pred_kernel2 <- melt( getPredKernel(bio_pars)[, as.character(predsize2), , drop = FALSE])

intercepts <- c( predsize, (1/beta)*predsize, (1/betaupp)*predsize, (1/betalow)*predsize)
intercepts2 <- c( predsize2, (1/beta)*predsize2, (1/betaupp)*predsize2, (1/betalow)*predsize2)

ggplot(pred_kernel) + geom_line(aes(x = w_prey, y = value)) + xlim(c(0, predsize)) + 
  geom_vline( xintercept = intercepts, linetype=c(1,1,2,2), color=c('blue4','red', 'red','red')) + 
  annotate( "rect", xmin = (1/betalow)*predsize, xmax = (1/betaupp)*predsize, ymin = 0, ymax = 1, fill = "red", alpha = .05) + 
  theme_bw() + 
  labs(title=paste0('Prey size for a ', predsize,' g Predator'), x = 'Prey weight (g)', y = 'Kernel density')

ggplot(pred_kernel2) + geom_line(aes(x = w_prey, y = value)) + xlim(c(0, predsize2)) + 
  geom_vline( xintercept = intercepts2, linetype=c(1,1,2,2), color=c('blue4', 'red', 'red','red')) + theme_bw() + 
  annotate( "rect", xmin = (1/betalow)*predsize2, xmax = (1/betaupp)*predsize2, ymin = 0, ymax = 1, fill = "red", alpha = .05) + 
  labs(title=paste0('Prey size for a ', predsize2,' g Predator'), x = 'Prey weight (g)', y = 'Kernel density')



