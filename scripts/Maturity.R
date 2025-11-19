
rm(list=ls())

library(tidyr)
library(ggplot2)
library(INLA)

fdata <- read.csv("./data/Maturity/mat_data.csv", header = T, check.names = FALSE, sep = ";" ,dec = ".", stringsAsFactors = F)
fdata$year_mat <- as.factor(fdata$year_mat)
fdata$sex <- as.factor(fdata$sex)
fdata$lab <- as.factor(fdata$lab)
fdata$month <- as.factor(fdata$month)

last <- unique(fdata$year)

MatSize <- matrix( NA, 3, 2, dimnames = list( c('L50','L25', 'k'), c('Males','Females')))

for( i in c(1,2)){
  
  data <- subset( fdata, fdata$sex == i)
  head(data)
  
  ind <- is.na(data$mat)
  ind <- which(ind==TRUE)
  data <- data[-ind,]
  
  cutoff_lengths <- c(seq(min(data$lt),19,by=1),seq(from=20, to=40, by=1), seq(from=42, to=70, by=2),seq(71,max(data$lt),by=4))
  data$bin <- cut(data$lt, cutoff_lengths, labels = cutoff_lengths[-1])
  
  aux <- subset(data,data$lt<21)[,c(3,5)]
  ind <- which(aux$mat==1)
  
  data$mat[data$lt < 21 ] <- 0
  
  NLbins<-c(seq(from=20, to=40, by=1),seq(from=42, to=70, by=2)) # Desired bins (SS model) 67
  l_b <- length(NLbins)
  
  len <- data$lt
  l_len <- length(len); aux <- rep(0,l_len)
  
  years <- (min(as.numeric(as.character(data$year_mat))):max(as.numeric(as.character(data$year_mat))))
  
  data_ieo <- subset(data,data$lab=="ieo")
  data_ipma <- subset(data,data$lab=="ipma")
  data <- rbind(data_ieo,data_ipma)
  
  ind_ieo <- which(data$lab=="ieo"); ind_ipma <- which(data$lab=="ipma")
  len <- length(data$lab)
  
  len_ieo <- length(ind_ieo)
  len_ipma <- length(ind_ipma)
  
  YCombined <- matrix(NA, nrow = len, ncol = 2)
  YCombined[1:len_ieo, 1] <- (data$mat[ind_ieo])
  YCombined[(len_ieo+1):(len_ipma+len_ieo), 2] <- (data$mat[ind_ipma])
  
  data$Gyear_mat <- as.character(data$year_mat)
  ind <- which(as.numeric(as.character(data$year_mat))<2001)
  data$Gyear_mat[ind] <- "1980_2000"
  
  ind <- which(as.numeric(as.character(data$year_mat))>2016)
  data$Gyear_mat[ind] <- "2017-2019"
  
  data$Gyear_mat <- as.factor(data$Gyear_mat)
  
  f3 <-  YCombined ~ 1 + lt + f( Gyear_mat, model = "iid")
  
  I3 <- inla( f3, control.compute = list(config=TRUE, dic = TRUE, cpo=TRUE),
              family = c("binomial","binomial"), data = data, 
              control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads = 1)
  
  intercept <- I3$summary.fixed[1, 1]; coef_lt <- I3$summary.fixed[2, 1]
 
  MatSize['L50',i] <- - intercept/coef_lt
  MatSize['L25',i] <- (log(0.25 / (1 - 0.25)) - intercept) / coef_lt
  MatSize['k',i] <- coef_lt
  
}


MatSize  

source( './scripts/aux_functions.R')

replist <- r4ss::SSgetoutput( dirvec = "./data/WGBIE24", getcovar = F, verbose = FALSE)[[1]]

sspars <- replist$parameters

grpars <- replist$Growth_Parameters     ## no males, then:
grpars <- readLines("./data/WGBIE24/Report.sso")
pos <- grep("^Growth_Parameters report:45", grpars)
grpars <- grpars[(pos + 1):(pos + 4 + 1)]
grpars <- read.table(text = paste(grpars, collapse = "\n"), header = TRUE)

a <- grpars$WtLen1[1] * 1e3; a   # 0.00377 (kg to g)
b <- grpars$WtLen2[1]; b        # 3.168 

Kvb <- grpars$K[1]; Kvb

Linf_f <- grpars$Linf[1]; Linf_f
Linf_m <- grpars$Linf[3]; Linf_m     # Linf_f*exp(sspars['L_at_Amax_Mal_GP_1','Value'])
Linf <- (Linf_f+Linf_m)/2; Linf

al0_f <- grpars$A_a_L0[1]; al0_f
al0_m <- grpars$A_a_L0[3]; al0_m     # from SS report file (report:45)
al0 <- (al0_f+al0_m)/2; al0

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

lengths <- 0:129
weights <- lwf( lengths, a, b)

# males <- 1/(1+exp(log(3)*((MatSize['L50','Males']-lengths)/(MatSize['L50','Males']-MatSize['L25','Males']))))
# females <- 1/(1+exp(log(3)*((MatSize['L50','Females']-lengths)/(MatSize['L50','Females']-MatSize['L25','Females']))))

males <- 1/(1+exp(-kmat_m*(lengths-L50_m)))
females <- 1/(1+exp(-kmat_f*(lengths-L50_f)))

meanmf <-  1/(1+exp(-kmat*(lengths-L50)))
mizer <- 1/(1+(weights/w50)^(-U))

msplot <- data.frame( Length = lengths, Weight = weights, 
  Male = males, Female = females, Mean = meanmf, Mizer = mizer)

msplot <- msplot %>% pivot_longer( cols = Male:Mizer, names_to = "Model", values_to = "Maturity")

lmplot <- msplot %>% ggplot( aes( x = Length, y = Maturity, col = Model)) + 
  geom_line( ) + theme_bw() + labs( x = "Length (cm)", color= "Sex")
lmplot

wmplot <- msplot %>% ggplot( aes( x = Weight, y = Maturity, col = Model)) + 
  geom_line( ) + theme_bw() + labs( x = "Weight (g)", color= "Sex")
wmplot

logwmplot <- msplot[-c(1:2),] %>% ggplot( aes( x = Weight, y = Maturity, col = Model)) + scale_x_log10() +
  geom_line( ) + theme_bw() + labs( x = "log( weight (g))", color= "Sex") 
logwmplot

lmplot + geom_vline( xintercept = c(L50_f, L50_m, L50), linetype = c( 'dashed', 'dashed', 'solid')) +
  # geom_vline( xintercept = c(L25_f, L25_m, L25), linetype = c( 'dashed', 'dashed', 'solid'), color = 'red') +
  geom_hline( yintercept = c(0.25,0.5), linetype = 'dashed')

logwmplot + geom_vline( xintercept = c(lwf( L50_m, a, b), lwf( L50_f, a, b), w50), linetype = c( 'dashed', 'dashed', 'solid')) +
  # geom_vline( xintercept = c(lwf( L25_m, a, b), lwf( L25_f, a, b), w25), linetype = c( 'dashed', 'dashed', 'solid'), color = 'red') +
  geom_hline( yintercept = c(0.25, 0.5), linetype = 'dashed')


pdf("./data/plots/maturity.pdf", width = 10, height = 6, onefile = TRUE)
print(lmplot)
print(wmplot)
print(logwmplot)
dev.off()
  
save( MatSize, file = './data/Maturity_Size.RData')
  
