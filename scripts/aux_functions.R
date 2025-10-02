### Von Bertalanffy

alf <- function( a, linf, K, al0)  linf * ( 1 - exp( -K * (a - al0)))
laf <- function( l, linf, K, al0)  al0 - log(1 - l / linf) / K


### Weight - Length

lwf <- function( l, a, b)  a*l^b
wlf <- function( w, a, b) (w/a)^(1/b)


### Logistic function: f(x) = 1 / (1 + exp(-b (x - c)))

logf <- function(x, s, a50)  1 / (1 + exp( -s * (x - a50)))
logf_shift <- function(x, s, a50, shift)  1 / (1 + exp( -s * (x - a50 - shift)))


### Line equation

linef <- function( x, L0, L1, M0, M1)  M0 + (M1-M0)*(x-L0)/(L1-L0)


### Power Low function

powerlow <- lwf


### Double Normal Selectivity

double_normal_length <- function( w, params, species_params) {
  
  # 'params' vector elements:
  # 1: location of ascending peak
  # 2: logistic scale width of plateau
  # 3: ascending slope log scale value (slope = exp(p[3]))
  # 4: descending slope log scale value (slope = exp(p[4]))
  # 5: sel value at initial bin;
  # 6: sel value at final bin;
  
  # Preemptive checks
  assert_that(is.numeric(w) && is.numeric(params))
  assert_that(params[1] > 0)
  assert_that(params[5]>=0, params[5]<=1, params[6]>=0, params[6]<=1)
  
  # Weight-length conversion
  a <- species_params[["a"]]
  b <- species_params[["b"]]
  # if (is.null(a) || is.null(b)) {
  #   stop("The selectivity function needs the weight-length parameters ", 
  #        "`a` and `b` to be provided in the species_params data frame.")
  # }
  l <- (w/a)^(1/b)
  
  # Extract parameters
  peak1 <- params[1]
  upselex <- exp(params[3])
  downselex <- exp(params[4])
  
  # Compute derived parameters and points
  peak2 <- peak1 + (0.99 * max(l) - peak1) / (1 + exp(-params[2]))
  point1 <- ifelse(params[5] > 0, params[5], NA)
  point2 <- ifelse(params[6] > 0, params[6], NA)
  
  # Precompute initial scaling factors if points are defined
  t1min <- if (!is.na(point1)) exp(-((min(l) - peak1)^2) / upselex) else NA
  t2min <- if (!is.na(point2)) exp(-((max(l) - peak2)^2) / downselex) else NA
  
  # Vectorized computation for asc and dsc selectivity across x
  t1 <- l - peak1
  t2 <- l - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  
  # Ascending and descending selectivity calculations
  asc <- exp(-t1^2 / upselex)
  dsc <- exp(-t2^2 / downselex)
  
  # Scale asc and dsc selectivity if points are defined
  asc_scl <- if (!is.na(point1)) point1 + (1 - point1) * (asc - t1min) / (1 - t1min) else asc
  dsc_scl <- if (!is.na(point2)) 1 + (point2 - 1) * (dsc - 1) / (t2min - 1) else dsc
  
  # Compute final selectivity using vectorized operations
  sel <- asc_scl * (1 - join1) + join1 * (1 - join2 + dsc_scl * join2)
  
  # Plot and return
  plot(l, sel, col = "red", type = "l", ylab = "Selectivity", xlab = "Length")
  return(sel)
}


double_logistic <- function(w, pars) {
  sel_left  <- 1 / (1 + exp(-log(19) * (w - pars[1]) / (pars[1] - pars[2])))
  sel_right <- 1 / (1 + exp( log(19) * (w - pars[3]) / (pars[4] - pars[3])))
  sel <- sel_left * sel_right
  return(sel)
}


double_normal_sel <- function( w, params) {
  
  assert_that(is.numeric(w) && is.numeric(params))
  assert_that(params[1] > 0)
  assert_that(params[5]>=0, params[5]<=1, params[6]>=0, params[6]<=1)
  
  peak1 <- params[1]
  upselex <- exp(params[3])
  downselex <- exp(params[4])

  peak2 <- peak1 + (0.99 * max(w) - peak1) / (1 + exp(-params[2]))
  point1 <- ifelse(params[5] > 0, params[5], NA)
  point2 <- ifelse(params[6] > 0, params[6], NA)

  t1min <- if (!is.na(point1)) exp(-((min(w) - peak1)^2) / upselex) else NA
  t2min <- if (!is.na(point2)) exp(-((max(w) - peak2)^2) / downselex) else NA
  
  # Vectorized computation for asc and dsc selectivity across x
  t1 <- w - peak1
  t2 <- w - peak2
  join1 <- 1 / (1 + exp(-(20 / (1 + abs(t1))) * t1))
  join2 <- 1 / (1 + exp(-(20 / (1 + abs(t2))) * t2))
  
  # Ascending and descending selectivity calculations
  asc <- exp(-t1^2 / upselex)
  dsc <- exp(-t2^2 / downselex)
  
  # Scale asc and dsc selectivity if points are defined
  asc_scl <- if (!is.na(point1)) point1 + (1 - point1) * (asc - t1min) / (1 - t1min) else asc
  dsc_scl <- if (!is.na(point2)) 1 + (point2 - 1) * (dsc - 1) / (t2min - 1) else dsc
  
  # Compute final selectivity using vectorized operations
  sel <- asc_scl * (1 - join1) + join1 * (1 - join2 + dsc_scl * join2)
  
  # Plot and return
  plot(w, sel, col = "red", type = "l", ylab = "Selectivity", xlab = "Length")
  return(sel)
}


### Energies

getEnergy <- function( model, return_df = FALSE, log = FALSE){
  
  sp <- species_params(model)
  wt <- w(model)
  
  merepro <- as.numeric(getERepro(model))
  mrepp <- as.numeric(getReproductionProportion(model))
  mgrowth <- as.numeric(getEGrowth(model))
  mmetab <- as.numeric(getMetabolicRate(model))
  merg <- as.numeric(getEReproAndGrowth(model))
    
  metab <- sp$ks * wt^sp$p
  activ <- sp$k * wt
  totale <- sp$alpha*sp$f0*(sp$h*wt^sp$n)
  ereproandgrowth <- totale-metab-activ
  repp <- ((wt/sp$w_max)^(1-sp$n))
  phi <- (1/(1+(wt/sp$w_mat)^(-sp$U))) * (wt/sp$w_max)^(1-sp$n)
  
  egrowth <- ereproandgrowth * (1-phi)
  erepro <- ereproandgrowth * (phi)
  
  
  dafr <- rbind( data.frame( ERepro = erepro, ReproProp = repp, EGrowth = egrowth, EReproandGrowth = ereproandgrowth,  
                             Weight = wt, Type = 'By Hand'),
                 data.frame( ERepro = merepro, ReproProp = mrepp, EGrowth = mgrowth, EReproandGrowth = merg, 
                             Weight = wt, Type = 'MIZER'))
  
  dafr_long <- dafr %>% pivot_longer( cols = c(ERepro, ReproProp, EGrowth, EReproandGrowth),
      names_to = "Variable", values_to = "Value")
  
  pl <- ggplot( dafr_long, aes(x = Weight, y = Value, color = Type)) +
    geom_line() + facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
    theme_bw() + theme( legend.title = element_blank(), strip.text = element_text(face = "bold")) +
    labs(y = NULL) + guides(color = guide_legend(title = NULL, override.aes = list(alpha = 1)))
  
  if(log) pl <- pl + scale_x_log10()
  
  print(pl)
  
  if( return_df) return(dafr)
  
}


### Plots LFDs

plot_lfd <- function( model, catch, return_df = FALSE) {
  
  heights <- aggregate(number ~ length, data = catch, FUN = sum)$number / sum(catch$number)
  
  lengths <- (model@w / model@species_params$a)^(1/model@species_params$b)
  
  model_catch <- model@initial_n * getFMort(model)
  model_catch <- model_catch  / sum(model_catch * model@dw)
  model_catch <- model_catch * model@species_params$b * model@w / lengths
  
  df <- rbind( data.frame( Length = unique(catch$length), Density = heights, Type = 'Observed'),
               data.frame( Length = lengths, Density = as.numeric(model_catch), Type = 'Estimated'))
  
  pl <- ggplot( df, aes(x = Length, y = Density, color = Type, fill = Type)) + 
    geom_bar(data = subset(df, Type != 'Estimated'), stat = "identity", position = "dodge", alpha = 0.6) + 
    geom_line(data = subset(df, Type == 'Estimated'), linewidth = 1) + theme_bw() +
    labs( x = "Size [cm]", y = "Normalised number density [1/cm]", color = NULL, fill = NULL)
  
  print(pl)
  
  if(return_df) return(df)
  
}

plot_lfd_gear <- function( model, catch, return_df = FALSE){
  
  params <- validParams(model)
  
  s <- params@species_params['species']
  a <- as.numeric(params@species_params["a"])
  b <- as.numeric(params@species_params["b"])
  w <- params@w
  l <- wlf(w,a,b)
  gears <- unique(catch$fleet)
  glength <- length(gears)
  
  df <- NULL
  
  for( i in 1:glength){
    
    igear <- gears[i]
    
    f_mort <- getFMortGear(params)[i,1,]
    
    catch_w <- f_mort * params@initial_n[1,]
    catch_w <- catch_w/sum(catch_w * params@dw)
    catch_l <- catch_w * b * w/l
    
    df <- rbind(df, data.frame( Weight=w, Length=l, Catch_w=catch_w, Catch_l=catch_l, 
                                Gear=igear, Type="Estimated"))
    
    cind <- which(catch$fleet==igear)
    
    len <- catch$length[cind]
    catch_l <- catch$number[cind]
    catch_l <- catch_l/sum(catch_l)
    wei <- catch$weight[cind]
    catch_w <- catch_l/b * len/wei
    
    df <- rbind(df, data.frame(Weight=wei, Length=len, Catch_w=catch_w, Catch_l=catch_l, 
                               Gear=igear, Type = "Observed"))
    
  }
  
  pl <-   ggplot( df, aes(x = Length, y = Catch_l, color = Type, fill = Type)) + 
    geom_bar(data = subset(df, Type != 'Estimated'), stat = "identity", position = "dodge", alpha = 0.6) + 
    geom_line(data = subset(df, Type == 'Estimated'), linewidth = 1) + theme_bw() +
    facet_wrap( ~Gear) + 
    labs( x = "Size [cm]", y = "Normalised number density [1/cm]", color = NULL, fill = NULL)
  
  print(pl)
  
  if(return_df) return(df)
  
}



plot_wfd <- function( model, catch, return_df = FALSE) {
  
  Nw <- as.numeric(model@initial_n)
  w <- model@w
  gp <- gear_params(model)
  
  Cw_total <- rep(0, length(w))
  
  for (i in 1:nrow(gp)) {
    Fw <- double_logistic(w, as.numeric(gp[i, c("w50", "w25", "w50_right", "w25_right")]))
    Cw <- Nw * Fw
    Cw_total <- Cw_total + Cw
  }
  
  Cw_total <- Cw_total / sum(Cw_total * model@dw)
  df <- data.frame( Weight = w, Catch_w = Cw_total, Type = "Estimated")
    
  obs_w <- aggregate(number ~ weight, data = catch, FUN = sum)
  obs_w$number <- obs_w$number / sum(obs_w$number)
  
  df <- rbind(df, data.frame( Weight = obs_w$weight, Catch_w = obs_w$number, Type = "Observed"))
  
  pl <- ggplot(df, aes(x = Weight, y = Catch_w, color = Type, fill = Type)) +
    geom_bar(data = subset(df, Type == "Observed"),
             stat = "identity", position = "dodge", alpha = 0.6) +
    geom_line(data = subset(df, Type == "Estimated"), linewidth = 1) +
    theme_bw() + scale_x_log10(limits = c(1,max(w))) + 
    labs(x = "Weight [g]", y = "Normalised number density [1/g]", color = NULL, fill = NULL)
  
  print(pl)
  
  if (return_df) return(df)
  
}

plot_wfd_gear <- function( model, catch, return_df = FALSE){

  w <- model@w
  gp <- gear_params(model)
  Nw <- as.numeric(model@initial_n)
  
  df <- NULL
  
  for (i in 1:nrow(gp)) {
    
    Fw <- double_logistic(w, as.numeric(gp[i, c("w50", "w25", "w50_right", "w25_right")]))
    Cw <- Nw * Fw
    Cw <- Cw / sum(Cw * model@dw)
    
    df <- rbind(df, data.frame( Weight = w, Catch_w = Cw, Gear = gp$gear[i], Type = "Estimated"))
    
    cind <- which(catch$fleet == gp$gear[i])
    obs_w <- aggregate(number ~ weight, data = catch[cind, ], FUN = sum)
    obs_w$number <- obs_w$number / sum(obs_w$number)
    
    df <- rbind(df, data.frame( Weight = obs_w$weight, Catch_w = obs_w$number, Gear = gp$gear[i], Type = "Observed"))
  }
  
  pl <- ggplot( df, aes(x = Weight, y = Catch_w, color = Type, fill = Type)) + facet_wrap(~Gear) +
    geom_bar(data = subset(df, Type != 'Estimated'), stat = "identity", position = "dodge", alpha = 0.6) + 
    geom_line(data = subset(df, Type == 'Estimated'), linewidth = 1) + theme_bw() + scale_x_log10(limits = c(1,max(w))) + 
    labs( x = "Weight [g]", y = "Normalised number density [1/cm]", color = NULL, fill = NULL)
  
  print(pl)
  
  if(return_df) return(df)
  
}

