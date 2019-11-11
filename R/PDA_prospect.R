#' Run PDA on prospect
#'
#' @export

PDA_prospect <- function(settings,Spectrum_leaf_data,df_PFT,wl.min,wl.max,use_meta.analysis=TRUE,nrChains=4,nrIter=10000){
    
  iter <- nrIter
  settings_MCMC <- list(iterations = iter, nrChains = nrChains)
  
  PDA <- list()
  best_set<- best_run <- list()
  
  # Uniform prior
  Prospect_param_names <- c("Ca","Cb","Car","Cw","Cm","Nlayers","ssigma")
  pft_lowers <- c(chlorophyll_a = 0,chlorophyll_b = 0, carotenoids = 0,Cw = 0, SLA = 1, Nlayers = 1, ssigma = 0)
  pft_uppers <-  c(chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1, SLA = 100, Nlayers = 5, ssigma = 1)
  # prior <- BayesianTools::createUniformPrior(lower = pft_lowers, upper = pft_uppers)
  
  dis2find <- c('chlorophyll_a','chlorophyll_b','carotenoids','Cw','SLA','Nlayers',"ssigma")
  
  N=10000
  
  # Define likelihood
  create_likelihood <- function(observed, waves) {
    function(params) {
      
      # PROSPECT parameters
      Ca <- params[1]
      Cb <- params[2]
      Cab <- Ca + Cb
      Car <- params[3]
      Cw <- params[4]
      Cm <- params[5]
      Nlayers <- params[6]
      
      ssigma <- params[7]
      
      optical_param <- c(Nlayers,Cab,Car,Cw,Cm)
      names(optical_param) <- c('Nlayers','chlab','carotenoids','Cw','Cm')
      
      # Call RTM
      result <- tryCatch(
        PEcAnRTM::prospect(optical_param, version = "5"),
        error = function(e) NULL)
      if (is.null(result)) return(-1e20)
      reflectance <- result[,"reflectance"]
      if (any(!is.finite(reflectance))) return(-1e20)
      if (any(reflectance < 0)) return(-1e20)
      if (any(reflectance > 1)) return(-1e20)
      
      reflectance_waves <- interp1(x = PEcAnRTM::wavelengths(reflectance),
                                   y = as.vector(matrix(reflectance)),
                                   waves)
      # Calculate likelihood
      ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
      return(ll)
    }
  }
  
  
  #############################################################################################
  # Loop over PFTs
   
  for (ipft in seq(df_PFT$names)){
    current_pft <- as.character(df_PFT$names[ipft])
    
    # Define priors
  
    # From posterior
    postfile <- file.path(settings$outdir,"pft",current_pft,"post.distns.Rdata")
    if (file.exists(postfile) & use_meta.analysis){
      load(postfile)
      distns <- post.distns
      Distributions <- rownames(post.distns)
    } else {
      priorfile <- file.path(settings$outdir,"pft",current_pft,"prior.distns.Rdata")
      load(priorfile)
      distns <- prior.distns
      Distributions <- rownames(prior.distns)
    }
  
    sampler <- matrix(NA,N,length(dis2find))
    lower <- upper <- best <- rep(NA,length(dis2find))
    
    for (idis in seq(dis2find)){
      current_dis = dis2find[idis]
      pos <- which(Distributions == current_dis)
      
      if (isempty(pos)){
        unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
        sampling <- get.sample(unif_distribution,n=N)
        
      } else {
        sampling <- get.sample(distns[pos,],n=N)
      }
      
      if (current_dis == 'SLA'){
        sampling = 1/sampling/10 # unit conversions
      }
      sampler[,idis] <- sampling
      lower[idis] <- min(sampler[,idis])
      upper[idis] <- max(sampler[,idis])
      best[idis]  <- median(sampler[,idis])
    }
    
    # 1/as.numeric(config_temp$SLA)/10*2
    prior <- createPriorDensity(sampler, method = "multivariate", eps = 1e-10,
                       lower = lower, upper = upper, best = best)
    
    # Read Leaf spectra 
    temp <-Spectrum_leaf_data %>% filter(pft==current_pft & wavelength>wl.min & wavelength<wl.max) %>% select(c('wavelength','Reflectance'))
    
    observation <- as.vector(temp$Reflectance)
    waves <- temp$wavelength
    
    likelihood <- create_likelihood(observation, waves)
    
    # Run inversion
    setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
    samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
    samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)
    
    for (ichain in seq(nrChains)){
      samples[[ichain]][["setup"]][["names"]] <- Prospect_param_names
    }

    
    PDA[[current_pft]] <- samples
    
  }
  
  return(PDA)
}
