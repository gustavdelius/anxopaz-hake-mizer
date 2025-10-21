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

w <- modelo@w
lw <- a*l^b
lwb <- a*lb^b

logw_min <- log(sp$w_min)
logw_max <- log(sp$w_max)
w <- exp(seq(logw_min, logw_max, length.out = 129))

logw_min_pp <- log(modelo@w_full[1])
logw_max    <- max(log(w))
wp <- exp(seq(logw_min_pp, logw_max, length.out = 288))

dw <- diff(w); dw <- c( dw, 0)
dwp <- diff(wp); dwp <- c( dwp, 0)

Nres <- function( wp, kappa, lambda){ kappa*wp^(-lambda)}
Nr <- Nres( wp, modelo@resource_params$kappa, modelo@resource_params$lambda)

interM <- matrix( c( 0, 1, 0, pcann), nrow=2, dimnames = list( pred = c('resource', 'Hake'), prey = c('resource', 'Hake')))

pred_kernel <- function( w, wp, beta, sigma){
  
  predk <- array( NA, c( length(w), length(wp)), dimnames = list( pred = signif(w,4), prey = signif(wp,4)))
  
  for(i in w) for(j in wp) 
    predk[as.character(signif(i,4)),as.character(signif(j,4))] <- 
      exp( -( (log(i/j) - log(beta))^2) / (2 * sigma^2)) 
  
  return(predk)
}

predk <- pred_kernel( w, wp, beta, sigma); dim(predk)

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

Eencounter <- function( w, Eavail, gamma, q){ Eavail * gamma * w^q}
Eenc <- Eencounter( w, Eavail, gamma, q)

feeding_level <- function( w, Eenc, h, n){ Eenc/(Eenc + h*w^n)}
feed_level <- feeding_level( w, Eenc, h, n)

max_consumption <- function( w, h, n){ h*w^n}
max_cons <- max_consumption( w, h, n)

Emetabolism <- function( w, ks, p){  ks*w^p}
Emet <- Emetabolism( w, ks, p)

Emovement <- function( w, k){  k*w}
Emov <- Emovement( w, k)

Ereproandgrowth <- function( w, h, n, ks, p, k, alpha, feed_level){
  
  emet <- ks*w^p
  emov <- k*w
  maxc <- h*w^n
  
  erg <- alpha * feed_level * maxc - emet - emov
  for(i in 1:length(erg)) erg[i] <- max(0, erg[i])
  
  return(erg)
  
}

Erepgro <- Ereproandgrowth( w, h, n, ks, p, k, alpha, feed_level)

repro_prop<- function( w, wmax, n){ (w/wmax)^(1-n)}
repp <- repro_prop( w, w_max, n)

Maturity <- function( w, wmat, U){ (1+(w/wmat)^(-U))^(-1)}
matt <- Maturity( w, w_mat, U)

psi <- function( w, wmax, wmat, U, n){
  repp <- (w/wmax)^(1-n)          # repro_prop( w, wmax, n)
  mat <- (1+(w/wmat)^(-U))^(-1)   # Maturity( w, wmat, U)
  fi <- repp*mat
  return(fi)
}
fi <- psi( w, w_max, w_mat, U, n)

Egrowth <- function( Erg, phi){ Erg*(1-phi)}
Egro <- Egrowth( Erepgro, fi)

Erepro <- function( Erg, phi){ Erg*(phi)}
Erep <- Erepro( Erepgro, fi)

FishMw <- function(w, gpw, eff = rep(0.25, nrow(gpw))) {
  
  w50  <- gpw$w50;  w25  <- gpw$w25
  w50r <- gpw$w50_right; w25r <- gpw$w25_right
  q <- gpw$catchability; eff <- rep(0.25, length(q))
  
  s2  <- log(3) / (w50  - w25)
  s1  <- s2 * w50
  s2r <- log(3) / (w50r - w25r)
  s1r <- s2r * w50r
  
  sel <- FM <- array(NA, c(length(w), length(w50)), dimnames = list(w, gpw$gear))
  
  for(i in 1:ncol(sel)) {
    sel[, i] <- (1 / (1 + exp(s1[i]  - s2[i]  * w))) *
      (1 / (1 + exp(s1r[i] - s2r[i] * w)))
    FM[, i]  <- q[i] * eff[i] * sel[, i]
  }
  
  return(FM)
}

FM <- FishMw(w, gpw, eff)

NatM <- function( w, M, d){ M*w^d}
NM <- NatM( w, M, d)

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
