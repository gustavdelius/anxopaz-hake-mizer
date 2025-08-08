
MIZER <- function( model, catch, dl = 1, nofixed = NULL, fixed_q = FALSE, yield_lambda = 1e7, biomass_lambda = 1e7) {
  
  # Prefit
  
  params <- validParams(model)
  sp <- params@species_params
  gp <- gear_params(model)
  
  gears <- unique(catch$fleet)
  n_g <- length(gears)
  
  lengths <- sort(unique(catch$length))
  weights <- sp$a * lengths^sp$b
  l_bins <- seq( lengths[1]-dl/2, lengths[length(lengths)]+dl/2, by=dl) 
  w_bins <- sp$a * l_bins^sp$b
  w_bin_widths <- diff(w_bins)
  
  counts <- catch %>% pivot_wider( names_from = fleet, values_from = number, values_fill = 0)
  counts <- as.matrix(counts[,gears])
  
  data_list <- list( counts = counts, bin_widths = w_bin_widths, weight = weights, blength = lengths, 
    yield = gp$yield_observed, biomass = sp$biomass_observed, min_len = min(l_bins), max_len = max(l_bins), 
    yield_lambda = yield_lambda, biomass_lambda = biomass_lambda, n_g = n_g)
  
  pars <- list(
    logit_l50 = qlogis((gp$l50 - min(l_bins))/(max(l_bins) - min(l_bins))),
    log_ratio_left = log((gp$l50 - gp$l25)/gp$l50),
    log_l50_right_offset = log(pmax(1e-3, gp$l50_right - gp$l50)),
    log_ratio_right = log((gp$l25_right - gp$l50_right)/gp$l50_right),
    log_catchability = log(gp$catchability))
  
  lower_bounds <- upper_bounds <- NULL
  
  if( fixed_q){ for(i in names(pars)){ 
    ipars <- as.numeric(pars[[i]])
    names(ipars) <- paste0(i,1:n_g)
    lower_bounds <- c(lower_bounds,ipars)} 
    upper_bounds <- lower_bounds
    
  } else { for(i in names(pars)){
    lower_bounds[paste0(i,1:n_g)] <- rep(-Inf,n_g)
    upper_bounds[paste0(i,1:n_g)] <- rep(Inf,n_g)}
  }
  
  parvec <- c('U','M','d','h','n','ks','p','k','f0','alpha','w_mat')
  
  if( sp['k'] <= 0) sp['k'] <- 1e-7
  
  for( i in parvec){
    
    idata <- as.numeric(sp[i])
    
    if( i == 'd'){
      iname <- i; pars[[iname]] <- idata
    } else {
      iname <- paste0('log_',i); pars[[iname]] <- ifelse(i %in% c('f0','alpha','n'), qlogis(idata), log(idata))}
    
    if( i %in% nofixed){ lower_bounds[iname] <- -Inf; upper_bounds[iname] <- Inf} else {
      lower_bounds[iname] <- upper_bounds[iname] <- pars[[iname]]}
    
  }
  
  # Fit
  
  dyn.load( dynlib("./TMB/fit"))
  
  obj <- MakeADFun( data = data_list, parameters = pars, DLL = "fit")
  
  optim_result <- nlminb( obj$par, obj$fn, obj$gr, lower = lower_bounds, upper = upper_bounds,
                          control = list( eval.max = 10000, iter.max = 10000))
  
  # Update model
  
  newpars <- optim_result$par
  
  lmin <- min(lengths)
  lmax <- max(lengths)
  
  sp <- model@species_params
  gp <- model@gear_params
  
  gplist <- list()
  gpnames <- c('logit_l50','log_ratio_left','log_l50_right_offset','log_ratio_right','log_catchability')
  for (i in gpnames) gplist[[i]] <- as.numeric(newpars[grep(i, names(newpars))])
  
  l50 <- lmin + (lmax - lmin) * plogis(gplist$logit_l50)
  l25 <- l50 * (1 - exp(gplist$log_ratio_left))
  l50_right <- l50 + exp(gplist$log_l50_right_offset)
  l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
  catchability <- exp(gplist$log_catchability)
  
  gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)
  
  gp[,'l50'] <- gp_res$l50
  gp[,'l25'] <- gp_res$l25
  gp[,'l50_right'] <- gp_res$l50_right
  gp[,'l25_right'] <- gp_res$l25_right
  gp[,'catchability'] <- gp_res$catchability
  
  gear_params(model) <- gp
  
  for(i in c('M','U','h','k','ks','p','alpha','w_mat')) sp[i] <- exp(newpars[paste0('log_',i)])
  for(i in c('f0','alpha','n')) sp[i] <- plogis(newpars[paste0('log_',i)])
  sp['d'] <- newpars['d']
  
  model@species_params <- sp
  
  model <- steadySingleSpecies(model)
  
  totalb <- sum(model@initial_n * model@w * model@dw)  # == getBiomass(model)
  diffb <- sp$biomass_observed / totalb
  model@initial_n <- model@initial_n * diffb
  
  totaly <- getYield(model)
  diffy <- sum(gp$yield_observed) / totaly  # != gp$yield_observed / totaly
  model@initial_effort <- model@initial_effort * diffy
  
  return(model)
  
}

