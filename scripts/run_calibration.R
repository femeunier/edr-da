## library(redr)
devtools::load_all(".")

# library(devtools)
# install_github("ashiklom/pecan", subdir = "modules/rtm",force = TRUE,dependencies = FALSE)
# install_github("PecanProject/pecan", subdir = "models/ed",force = TRUE,dependencies = FALSE)
# devtools::install_github("ashiklom/rrtm",dependencies = FALSE)
# /kyukon/home/gent/425/vsc42558/R/x86_64-pc-linux-gnu-library/3.5

library(PEcAnRTM)
library(PEcAn.priors)
library(rrtm)
library(dplyr)
library(ggplot2)
library(tidyr)

load_local <- function(file) {
  menv <- new.env(parent = baseenv())
  load(file, envir = menv)
  as.list(menv)
}

Ntests <- 10
Niter <- 100
nrChains <-  2
ncore <- 2
# srcdir <- "/data/gent/vo/000/gvo00074/felicien/IC/BCI"
srcdir <- "/home/femeunier/Documents/data/BCI"
outdir <- "./outputs/"
plot.figure <- TRUE

#####################################################################################

# params2retrieve <- round(runif(n = Nparameters))
# True set of parameters
names_parameters <- c("ssigma","soil_moisture","direct_sky_frac","czen","b1leaf",
                    "b2leaf","N","Cab","Car","Cw","Cm","orient_factor","clumping_factor")
params2retrieve <- c(1,0,0,0,
                     0,0,0,0,0,0,0,1,1)
Nparameters <- length(params2retrieve) 
params2retrieve[c(3,4)] <- 0
params2retrieve[c(1)] <- 1
Nparams2retrieve <- sum(params2retrieve)

params2retrieve_pos <- which(params2retrieve==1)

#####################################################################################
# All uniformed initially
pft_lowers <- rep(c(
  b1leaf = 0, b2leaf = -5, N = 1, Cab = 0, Car = 0, Cw = 0, Cm = 0,
  orient_factor = -0.5, clumping_factor = 0.4), npft)
pft_uppers <- rep(c(
  b1leaf = 3, b2leaf = 5,N = 5, Cab = 100, Car = 40, Cw = 0.1, Cm = 0.1,
  orient_factor = 0.5, clumping_factor = 1), npft)
lowers <- c(
  "ssigma" = 0, "soil_moisture" = 0, "direct_sky_frac" = 0, "czen" = 0,
  pft_lowers
)
uppers <- c(
  "ssigma" = 1, "soil_moisture" = 1, "direct_sky_frac" = 1, "czen" = 1,
  pft_uppers
)
#####################################################################################
cl <- parallel::makeCluster(ncore)

#####################################################################################
# Forest inventory
site.file <- file.path(srcdir,"ForestGEO.2012.lat9.1543lon-79.8461.css")
site_dat <- PEcAn.ED2::read_css(site.file) %>% filter(patch == 1)
patch_simple <- merge_cohorts(patch = site_dat,
                              nplant.col = "n",
                              dbh_diff = 5)

site_dat <- patch_simple[seq(dim(patch_simple)[1],1),]

dbh <- site_dat[["dbh"]]
pft <- site_dat[["pft"]] - max( site_dat[["pft"]]) + 1
nplant <- site_dat[["n"]]
npft <- length(unique(pft))
ncohort <- length(dbh)

# Define priors, either uniformed or very informed according to params2retrieve
prior <- BayesianTools::createUniformPrior(lower = lowers, upper = uppers)

########################################################################################
run_ED_RTM <- function(params_model,waves){
  # Define constants
  ## direct_sky_frac <- 0.9
  ## czen <- 1
  # Pull out site-specific params
  soil_moisture <- params_model[1]
  direct_sky_frac <- params_model[2]
  czen <- params_model[3]
  
  # Remaining params are pft-specific
  # Construct a parameter matrix -- nparam x npft
  # Each column contains the parameters in the order:
  # b1leaf, b2leaf, sla,
  # N, Cab, Car, Cw, Cm,
  # clumping_factor, orient_factor
  pft_params <- matrix(params_model[-(1:3)], ncol = npft)
  
  # Calculate allometries
  b1leaf <- pft_params[1, pft]
  b2leaf <- pft_params[2, pft]
  sla <- 1/(10*pft_params[7,pft])
  bleaf <- size2bl(dbh, b1leaf, b2leaf)
  lai <- nplant * bleaf * sla
  
  wai <-  rep(0, ncohort)
  cai <- rep(1, ncohort)
  
  # Extract remaining parameters
  N <- pft_params[3, ]
  Cab <- pft_params[4, ]
  Car <- pft_params[5, ]
  Cw <- pft_params[6, ]
  Cm <- pft_params[7, ]
  orient_factor <- pft_params[8, ]
  clumping_factor <- pft_params[9, ]
  
  # Call RTM
  result <- tryCatch(
    edr_r_v2(pft, lai, wai, cai,
             N, Cab, Car, Cw, Cm,
             orient_factor, clumping_factor,
             soil_moisture,
             direct_sky_frac,
             czen,
             wavelengths = waves),
    error = function(e) NULL)
  if (is.null(result)) return(NULL)
  albedo <- result[["albedo"]]
  return(albedo)
  
}

# Define likelihood
create_likelihood <- function(observed, waves, pft, dbh, nplant) {
  ncohort <- length(dbh)
  npft <- max(pft)
  stopifnot(
    length(pft) == ncohort,
    length(dbh) == ncohort,
    length(nplant) == ncohort
  )
  function(params){
    
    params_all <- ref_parameters
    params_all[params2retrieve_pos] <- params  
    ssigma <- params_all[1]
    params_model <- params_all[-c(1)]
    albedo <- run_ED_RTM(params_model,waves)
    
    if (is.null(albedo)) return(-1e20)
    if (any(!is.finite(albedo))) return(-1e20)
    if (any(albedo < 0)) return(-1e20)
    if (any(albedo > 1)) return(-1e20)
    
    # Calculate likelihood
    ll <- sum(dnorm(albedo, observed, ssigma, log = TRUE))
    ll
  }
}

for (itest in seq(Ntests)){
  
  # Create observation
  waves = 400:2500
  observation <- run_ED_RTM(params = ref_parameters[-c(1)],waves)
  # plot(waves,observation,type = 'l')
  
  # True set of parameters
  ref_parameters <- runif(n = length(prior$lower),prior$lower,prior$upper)
  names(ref_parameters) <- names_parameters
  
  N <- 10000
  sampler <- matrix(NA,N,Nparameters)
  lower <- upper <- best <- rep(NA,Nparameters)
  
  for (iparam in seq(1,Nparameters)){
    if (params2retrieve[iparam]){
      # Uninformed prior
      sampling <- runif(n=N,lowers[iparam],uppers[iparam])
    } else { 
      # Informed prior
      param_mean <- ref_parameters[iparam]
      param_sd <- abs(param_mean)/50
      norm_distribution <- list(distn="norm",parama=param_mean,paramb=param_sd)
      sampling <- get.sample(norm_distribution,n=N)
    }
    
    sampler[,iparam] <- sampling
    lower[iparam] <- min(sampler[,iparam])
    upper[iparam] <- max(sampler[,iparam])
    best[iparam]  <- median(sampler[,iparam])
  }
  
  prior_params <- BayesianTools::createPriorDensity(sampler[,params2retrieve_pos], method = "multivariate", eps = 1e-10,
                                             lower = lower[params2retrieve_pos], upper = upper[params2retrieve_pos], best = best[params2retrieve_pos])
  
  likelihood <- create_likelihood(observation, waves, pft, dbh, nplant)
  
  # Run inversion
  setup <- BayesianTools::createBayesianSetup(likelihood, prior_params, parallel = ncore)
  
  settings_MCMC <- BayesianTools::applySettingsDefault(settings = NULL, sampler = "DEzs", check = FALSE)
  settings_MCMC$iterations=Niter
  settings_MCMC$nrChains=nrChains
  settings_MCMC$startValue = setup$prior$sampler(ncore)
  
  samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
  samples <- BayesianTools::runMCMC(samples,settings = settings_MCMC)
  
  # Results
  Posteriors <- BayesianTools::getSample(samples,coda=FALSE,numSamples=2500)
  
  best_params <- BayesianTools::MAP(samples)$parametersMAP
  names(best_params) <- names(ref_parameters[params2retrieve_pos])
  
  waves_best = seq(400,2500,10)
  best_params_all <- ref_parameters
  best_params_all[params2retrieve_pos] <- best_params
  best_params_model <- best_params_all[-1]
  albedo <- run_ED_RTM(best_params_model,waves_best)
  
  data.opt <- rbind(as.data.frame(cbind(waves = waves_best,albedo)) %>% mutate(type = "best run"),
                    as.data.frame(cbind(waves,as.vector(observation))) %>% rename(albedo = V2) %>% mutate(type = "obs"))
  
  albedo_all <- run_ED_RTM(best_params_model,waves)
  data.opt_all <- rbind(as.data.frame(cbind(waves = waves,albedo = albedo_all)) %>% mutate(type = "best run"),
                        as.data.frame(cbind(waves,as.vector(observation))) %>% rename(albedo = V2) %>% mutate(type = "obs"))
  
  params2plot <- names(ref_parameters[params2retrieve_pos])[-1]
  
  Priors <- sampler[,params2retrieve_pos]
  colnames(Priors) <- colnames(Posteriors) <- names(ref_parameters[params2retrieve_pos])
  
  df_ref <- data.frame(name = (params2plot), value = as.matrix(ref_parameters[names(ref_parameters) %in% params2plot])) %>% mutate(type = "reference")
  df_bestrun <- data.frame(name = (params2plot), value = as.matrix(best_params[names(best_params) %in% params2plot])) %>% mutate(type = "best_set")
  
  all.distributions <-
    rbind(as.data.frame(Priors) %>% dplyr::select(params2plot) %>% pivot_longer(cols = params2plot) %>% mutate(type = 'prior'),
          as.data.frame(Posteriors) %>% dplyr::select(params2plot) %>% pivot_longer(cols = params2plot) %>% mutate(type = 'posterior'))
  
  all.OP <- rbind(all.distributions,
                  df_ref,df_bestrun)
  saveRDS(file = file.path(outdir,paste0("results",itest,".RDS")),object = all.OP)
}