
MIZER <- function( model, catch, w_from_catch = TRUE, compiler = FALSE,
                   nofixed = NULL, fixed_sel = FALSE, low_bounds = NULL, upp_bounds = NULL,
                   yield_lambda = 1e7, biomass_lambda = 1e7, cannibalism = 'add',
                   plot = F, plot_dir = getwd()) {
  
  # Prefit
  
  params <- validParams(model)
  sp <- params@species_params
  gp <- gear_params(model)
  
  gears <- unique(catch$fleet)
  n_g <- length(gears)
  
  if( w_from_catch){
    
    lengths <- sort(unique(catch$length))
    w <- lwf( lengths, sp$a, sp$b)
    
    counts <- catch %>% pivot_wider(names_from = fleet, values_from = number, values_fill = 0) %>% 
      select(all_of(gears)) %>% as.matrix()
    
  } else {
    
    w <- exp( seq( log(sp$w_min), log(sp$w_max), length.out = length(model@w)))
    lengths <- wlf( w, sp$a, sp$b)
    
    counts2 <- catch %>% group_by(fleet) %>%
      summarise(approx_counts = list( approx(x = weight, y = number, xout = w, rule = 2, ties = mean)$y)) %>%
      ungroup()
    
    counts <- do.call(cbind, counts2$approx_counts)
    colnames(counts) <- counts2$fleet
    rownames(counts) <- signif(w, 4)
    
  }
  
  lmin <- min(lengths)
  lmax <- max(lengths)
  
  wp <- c( exp(seq(log(model@w_full[1]), min(log(w)), length.out = 288-length(w)+1)), w)
  
  data_list <- list( counts = counts, w = w, wp = wp, yield = gp$yield_observed, biomass = sp$biomass_observed,
                     yield_lambda = yield_lambda, biomass_lambda = biomass_lambda, n_g = n_g)
  
  log_eff <- log(as.numeric(initial_effort(model)))
  
  pars <- list(
    logit_l50 = qlogis((gp$l50 - lmin)/(lmax - lmin)),
    log_ratio_left = log((gp$l50 - gp$l25)/gp$l50),
    log_l50_right_offset = log(pmax(1e-3, gp$l50_right - gp$l50)),
    log_ratio_right = log((gp$l25_right - gp$l50_right)/gp$l50_right),
    log_catchability = log(gp$catchability),
    log_effort = log_eff)
  
  lower_bounds <- upper_bounds <- NULL
  
  if( fixed_sel){ for(i in names(pars)){ 
    ipars <- as.numeric(pars[[i]])
    names(ipars) <- paste0(i,1:n_g)
    lower_bounds <- c(lower_bounds,ipars)} 
    upper_bounds <- lower_bounds
    
  } else { for(i in names(pars)){
    lower_bounds[paste0(i,1:n_g)] <- rep(-Inf,n_g)
    upper_bounds[paste0(i,1:n_g)] <- rep(Inf,n_g)}
  }
  
  for(i in 1:n_g) lower_bounds[paste0('log_effort',i)] <- upper_bounds[paste0('log_effort',i)] <- log_eff[i]
  
  for( i in c('kappa','lambda')){ 
    pars[[i]] <- model@resource_params[[i]]
    if( i %in% nofixed){ lower_bounds[i] <- -Inf; upper_bounds[i] <- Inf} else {
      lower_bounds[i] <- upper_bounds[i] <- pars[[i]]}
  }
  
  sp$inter_HR <- sp$interaction_resource
  sp$inter_HH <- as.numeric(model@interaction)
  
  parvec <- c( 'a', 'b','beta', 'sigma', 'inter_HR', 'inter_HH', 'gamma', 'q', 'h', 'n', 
               'ks', 'p', 'k', 'alpha', 'U', 'w_mat', 'w_min', 'w_max', 'M', 'd')
  
  parvecname <- parvec
  
  if( sp['k'] <= 0) sp['k'] <- 1e-7
  
  for( i in parvec){
    
    idata <- as.numeric(sp[i])
    
    if( i %in% c('d','inter_HR', 'inter_HH')){
      iname <- i; pars[[iname]] <- idata
    } else {
      iname <- paste0('log_',i)
      pars[[iname]] <- ifelse(i %in% c('alpha','n'), qlogis(idata), log(idata))}
    
    if( i %in% nofixed){ 
      lower_bounds[iname] <- -Inf; upper_bounds[iname] <- Inf
      if(i %in% c('inter_HH','inter_HR')){lower_bounds[iname] <- 0; upper_bounds[iname] <- 1}
    } else {
      lower_bounds[iname] <- upper_bounds[iname] <- pars[[iname]]}
    
  }
  
  if(!is.null(nofixed)){
    
    if(!is.null(upp_bounds)){
      bpars <- names(upp_bounds)
      for(i in bpars){ 
        if( i %in% c('d','inter_HR', 'inter_HH')){ upper_bounds[i] <- upp_bounds[[i]]
        } else { upper_bounds[paste0('log_',i)] <- ifelse(i %in% c('alpha','n'), qlogis(upp_bounds[[i]]), log(upp_bounds[[i]]))}
      }
    }
    
    if(!is.null(low_bounds)){
      bpars <- names(low_bounds)
      for(i in bpars){ 
        if( i %in% c('d','inter_HR', 'inter_HH')){ lower_bounds[i] <- low_bounds[[i]]
        } else { lower_bounds[paste0('log_',i)] <- ifelse(i %in% c('alpha','n'), qlogis(low_bounds[[i]]), log(low_bounds[[i]]))}
      }
    }
    
  }
  
  # Fit
  
  if( cannibalism == 'add'){
    
    if( compiler == TRUE){
      suppressWarnings(file.remove("./TMB/fit.o", "./TMB/fit.so"))
      TMB::compile("./TMB/fit.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
    }
    
    dyn.load( dynlib("./TMB/fit"))
    obj <- MakeADFun( data = data_list, parameters = pars, DLL = "fit")
    optim_result <- nlminb( obj$par, obj$fn, obj$gr, lower = lower_bounds, upper = upper_bounds,
                            control = list( eval.max = 1000, iter.max = 1000))
  } else {
    
    if( compiler == TRUE){
      suppressWarnings(file.remove("./TMB/fit_cann.o", "./TMB/fit_cann.so"))
      TMB::compile("./TMB/fit_cann.cpp", flags = "-Og -g", clean = TRUE, verbose = TRUE)
    }
    
    dyn.load( dynlib("./TMB/fit_cann"))
    obj <- MakeADFun( data = data_list, parameters = pars, DLL = "fit_cann")
    optim_result <- nlminb( obj$par, obj$fn, obj$gr, lower = lower_bounds, upper = upper_bounds,
                            control = list( eval.max = 1000, iter.max = 1000))
  }
  
  
  # Update model
  
  newpars <- optim_result$par
  
  wmin <- min(w)
  wmax <- max(w)
  
  sp <- model@species_params
  gp <- model@gear_params
  
  gplist <- list()
  gpnames <- c( 'logit_l50', 'log_ratio_left', 'log_l50_right_offset', 'log_ratio_right',
                'log_catchability', 'log_effort')
  
  for (i in gpnames) gplist[[i]] <- as.numeric(newpars[grep(i, names(newpars))])
  
  l50 <- lmin + (lmax - lmin) * plogis(gplist$logit_l50)
  l25 <- l50 * (1 - exp(gplist$log_ratio_left))
  l50_right <- l50 + exp(gplist$log_l50_right_offset)
  l25_right <- l50_right * (1 + exp(gplist$log_ratio_right))
  catchability <- exp(gplist$log_catchability)
  effort <- exp(gplist$log_effort)
  
  gp_res <- data.frame( l50 = l50, l25 = l25, l50_right = l50_right, l25_right = l25_right, catchability = catchability)
  
  gp[,'l50'] <- gp_res$l50
  gp[,'l25'] <- gp_res$l25
  gp[,'l50_right'] <- gp_res$l50_right
  gp[,'l25_right'] <- gp_res$l25_right
  gp[,'catchability'] <- gp_res$catchability
  
  gear_params(model) <- gp
  
  initial_effort(model) <- effort
  
  for(i in c( 'a','b','beta','sigma','gamma','q','M','U','h','k','ks','p','alpha','w_mat','w_min', 'w_max')) 
    sp[i] <- exp(newpars[paste0('log_',i)])
  
  for(i in c( 'alpha','n')) sp[i] <- plogis(newpars[paste0('log_',i)])
  
  sp['d'] <- newpars['d']
  
  model@interaction[1,1] <- newpars['inter_HH']
  
  for(i in c( 'kappa','lambda')) model@resource_params[[i]] <- newpars[i]
  
  model@species_params <- sp
  
  model <- steadySingleSpecies(model)
  
  if( w_from_catch){ model@initial_n <- model@initial_n * (sp$biomass_observed / sum(model@initial_n * model@w * model@dw))
  } else { model@initial_n[1,] <- obj$report(optim_result$par)$final_B / (w * c(diff(w),1))}
  
  plot_lfd( model, catch)
  if(plot == T){
    dir.create( path = plot_dir, showWarnings = TRUE, recursive = TRUE)
    ggsave( paste0( plot_dir, 'LFD.jpg'), width = 9, height = 7)
  } 
  
  plot_lfd_gear( model, catch)
  if(plot == T){
    ggsave( paste0( plot_dir, 'LFD_gear.jpg'), width = 9, height = 7)
  }
  
  return(model)
  
}

