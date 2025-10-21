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


# Data -------------

load("./fit.RData")
source( './scripts/aux_functions.R')

newcatch <- approx(x = catch$weight[1:129], y = catch$number[1:129], xout = wmizer, rule = 2)$y
plot( catch$weight[1:129], y = catch$number[1:129], type='l')
plot( catch$weight[1:129], y = catch$number[1:129], type='l', xlim=c(0,500))
plot( catch$weight[1:129], y = catch$number[1:129], type='l', xlim=c(0,800))
lines( w, newcatch, col='red')
lines( wmizer, newcatch, col='red')

modelo <- cannibal_hake

plotSpectra( modelo, power = 2) + theme_bw() 

sp <- species_params(modelo)
gp <- gear_params(modelo)
N <- modelo@initial_n
eff <- modelo@initial_effort

plotSpectra( modelo, power = 2)


## Parameters --------------------

a <- sp$a
b <- sp$b
h <- sp$h 
ks <- sp$ks 
alpha <- sp$alpha
beta <- sp$beta
gamma <- sp$gamma
sigma <- sp$sigma
n <- sp$n 
p <- sp$p
k <- sp$k
m <- sp$m
U <- sp$U
w_mat <- sp$w_mat
w_max <- sp$w_max
d <- sp$d
M <- sp$M
q <- sp$q


## Weights ----------------

l <- 1:129
dl <- 1
lb <- seq( l[1]-dl/2, l[length(l)]+dl/2, by=dl) 

w <- unique(catch$weight)
logw_min_pp <- log(modelo@w_full[1])
logw_max    <- min(log(w))
wp <- c( exp(seq(logw_min_pp, logw_max, length.out = 288-length(w)+1)),w[-1])

dw <- diff(w); dw <- c( dw, 0)
dwp <- diff(wp); dwp <- c( dwp, 0)

wmizer <- modelo@w
wpmizer <- modelo@w_full


## Resources ----------------------

Nres <- function( wp, kappa, lambda){ kappa*wp^(-lambda)}

reso <- plotSpectra( modelo); reso
resod <- reso$data$value[which(reso$data$Species=='Resource')]
resow <- reso$data$w[which(reso$data$Species=='Resource')]

Nr <- Nres( wp, modelo@resource_params$kappa, modelo@resource_params$lambda)

plot( wp, Nr, log = 'x', type='l')
plot( wp, log(Nr), type='l', log='x')
lines( resow, log( resod / (resow * modelo@w_full[2] / modelo@w_full[1])), col='red')



# Predation ------------------------

interM <- matrix( c( 0, 1, 0, pcann), nrow=2, dimnames = list( pred = c('resource', 'Hake'), prey = c('resource', 'Hake')))

pred_kernel <- function( w, wp, beta, sigma){
  
  predk <- array( NA, c( length(w), length(wp)), dimnames = list( pred = signif(w,4), prey = signif(wp,4)))
  
  for(i in 1:length(w)) for(j in 1:length(wp)) 
    predk[i,j] <- exp( -( (log(w[i]/wp[j]) - log(beta))^2) / (2 * sigma^2)) 
  
  return(predk)
}

predk <- pred_kernel( w, wp, beta, sigma)
pkt <- getPredKernel( modelo)

plot( wp, predk[length(w),], type='l', xlim=c(0,w[length(w)]*1.1),xlab='Weight (g)', ylab='Predation Kernel') 
abline(v=c(w[length(w)],w[length(w)]/beta), lty=2)
lines( wpmizer, pkt[1,length(w),], col='red')


# Energy --------------

Eavailable <- function( w, wp, dw, dwp, beta, sigma, interM, N, Nr){
  
  n_w <- length(w)
  enc <- numeric(n_w)
  
  for(i in 1:n_w){
    
    enc_re <- sum( predk[i,] * interM["Hake", "resource"] * Nr * wp * dwp)
    enc_os <- sum( predk[i,as.character(signif(w,4))] * interM["Hake", "Hake"] * N * w * dw)
    enc[i] <- enc_re + enc_os
    
  }
  
  enc <- as.numeric(enc)
  
  return(enc)
  
}

Eavail <- Eavailable( w, wp, dw, dwp, beta, sigma, interM, N, Nr)
plot( w, Eavail, log = 'x', type='l', xlim=c(1,w_max))


Eencounter <- function( w, Eavail, gamma, q){ Eavail * gamma * w^q}

Eenc <- Eencounter( w, Eavail, gamma, q)
plot( w, Eenc, log = 'x', type='l', xlim=c(1,w_max))
lines( wmizer, getEncounter(modelo), col='red')
range(Eenc/getEncounter(modelo))           # 1 1


feeding_level <- function( w, Eenc, h, n){ Eenc/(Eenc + h*w^n)}

feed_level <- feeding_level( w, Eenc, h, n)
plot( w, feed_level, log = 'x', type='l', ylim=c(0.56,0.62), xlim=c(1,w_max))
lines( w, getFeedingLevel(modelo), col='red')
range(feed_level/getFeedingLevel(modelo))         # 1 1


max_consumption <- function( w, h, n){ h*w^n}

max_cons <- max_consumption( w, h, n)
plot( w, max_cons, type='l')
lines( w, getMaxIntakeRate(modelo), col='red')
range(max_cons/getMaxIntakeRate(modelo))    # 1 1


Emetabolism <- function( w, ks, p){  ks*w^p}

Emet <- Emetabolism( w, ks, p)
plot( w, Emet, type='l')
lines( w, getMetabolicRate(modelo), col='red')
range(Emet/getMetabolicRate(modelo))    # 1 1


Emovement <- function( w, k){  k*w}

Emov <- Emovement( w, k)
plot( w, Emov, type='l')


Ereproandgrowth <- function( w, h, n, ks, p, k, alpha, feed_level){
  
  emet <- ks*w^p
  emov <- k*w
  maxc <- h*w^n
  
  erg <- alpha * feed_level * maxc - emet - emov
  for(i in 1:length(erg)) erg[i] <- max(0, erg[i])
  
  return(erg)
  
}

Erepgro <- Ereproandgrowth( w, h, n, ks, p, k, alpha, feed_level)
plot( w, Erepgro, type='l')
lines( w, getEReproAndGrowth(modelo), col='red')
range(Erepgro/getEReproAndGrowth(modelo))          # 0.9999999 1.0000001


repro_prop<- function( w, wmax, n){ (w/wmax)^(1-n)}

repp <- repro_prop( w, w_max, n)
plot( w, repp, type='l')
lines( w, getReproductionProportion(modelo), col='red')
range(repp-getReproductionProportion(modelo))    # -2.754523e-06  8.084945e-02


Maturity <- function( w, wmat, U){ (1+(w/wmat)^(-U))^(-1)}

matt <- Maturity( w, w_mat, U)
plot( w, matt, log='x', type='l', xlim=c(1,w_max))
abline( h=0.5, v=w_mat, lty='dashed')
lines( w, getMaturityProportion(modelo), col='red')
range(matt-getMaturityProportion(modelo))    # -0.01525598  0.01520397



psi <- function( w, wmax, wmat, U, n){
  repp <- (w/wmax)^(1-n)          # repro_prop( w, wmax, n)
  mat <- (1+(w/wmat)^(-U))^(-1)   # Maturity( w, wmat, U)
  fi <- repp*mat
  return(fi)
}

fi <- psi( w, w_max, w_mat, U, n)
plot( w, fi, log='x', type='l', xlim=c(1,w_max))


Egrowth <- function( Erg, phi){ Erg*(1-phi)}

Egro <- Egrowth( Erepgro, fi)
plot( w, Egro, log='x', type='l', xlim=c(1,w_max))
lines( w, getEGrowth(modelo), col='red')
range(Egro-getEGrowth(modelo))          # -5.743581  2.155383


Erepro <- function( Erg, phi){ Erg*(phi)}

Erep <- Erepro( Erepgro, fi)
plot( w, Erep, log='x', type='l', xlim=c(1,w_max))
lines( w, getERepro(modelo), col='red')
range(Erep-getERepro(modelo))          # -2.155354  5.743630




# Mortality --------------------

## Fishing ----------------------

FishM <- function( l, gp, eff = rep(0.25, nrow(gp))){
  
  l50 <- gp$l50; l25 <- gp$l25; l50r <- gp$l50_right; l25r <- gp$l25_right
  q <- gp$catchability; eff <- rep(0.25, length(q))
  
  sr <- l50 - l25
  s1 <- l50 * log(3) / sr
  s2 <- s1 / l50;
  
  srr <- l50r - l25r
  s1r <- l50r*log(3)/srr
  s2r <- s1r /l50r
  
  sel <- FM <- array( NA, c(length(l),length(l50)), dimnames = list( l, gp$gear))
  
  for(i in 1:ncol(sel)){
    sel[,i] <- (1/(1+exp(s1[i]-s2[i]*l)))*(1/(1+exp(s1r[i]-s2r[i]*l)))
    FM[,i] <- q[i]*eff[i]*sel[,i]
  }
  
  return(FM) 
  
}

FM <- FishM( l, gp, eff)

par( mfrow = c(3,3))
for(i in 1:ncol(FM)) plot( l, FM[,i], type='l', xlab ='Length (cm)', ylab='Selectivity', main=gp$gear[i])
par( mfrow = c(1,1))

par( mfrow = c(3,3))
for(i in 1:ncol(FM)) plot( lwf(l,a,b), FM[,i], xlim=c(1,w_max), log = 'x', type='l', 
                           xlab ='Weight (g)', ylab='Selectivity', main=gp$gear[i])
par( mfrow = c(1,1))

par(mfrow=c(3,3))
for(i in 1:ncol(FMw)) {
  plot(w, FMw[,i], log="x", type="l", xlim=c(1, w_max),
       xlab="Weight (g)", ylab="Selectivity", main=gpw$gear[i])
  lines(lwf(l,a,b), FM[,i], col="red") 
}

par(mfrow=c(1,1))

par(mfrow=c(3,3))
for(i in 1:ncol(FMw)) {
  plot( w, FMw[,i], log="x", type="l", xlim=c(1, w_max),
        xlab="Weight (g)", ylab="Selectivity", main=gpw$gear[i])
  lines( lwf(l,a,b), FM[,i], col="red") 
  lines( w, FMeq[,i], col="green", lty='dashed')
}

par(mfrow=c(1,1))




## Background --------------------

NatM <- function( w, M, d){ M*w^d}

NM <- NatM( w, M, d)

plot( w, NM, type = 'l', log = 'x', xlim=c(1,w_max), xlab = 'Weight (g)', ylab = 'Natural Mortality')
lines( w, getExtMort(modelo), col='red')
lines( w, getExtMort(modelo) + getPredMort(modelo), col='green')
lines( w, getExtMort(hake_mizer), col='blue')


## Predation -------------------------------

PredM <- function(w, dw, interM, gamma, q, N, feed_level, predk) {
  
  pmort <- numeric(length(w))
  
  for(i in 1:length(w)) { for(j in 1:length(w)) {
    pmort[i] <- pmort[i] + 
      interM["Hake","Hake"] * predk[as.character(signif(w[j],4)), as.character(signif(w[i],4))] * 
      (1 - feed_level[j]) * gamma * w[j]^q * N[j] * dw[j]
  }}
  
  return(pmort)
}


predm <- PredM( w, dw, interM, gamma, q, N, feed_level, predk)
pmizer <- as.numeric(getPredMort(modelo))

plot( w, predm, type = 'l', xlab = 'Weight (g)', ylab = 'Predation')
lines( w, pmizer, col='red')

range(predm-pmizer)      # -7.642208e-09  2.741028e-08

