rm(list = ls())

devtools::load_all(".")

library(redr)
library(PEcAnRTM)
library(rrtm)
library(dplyr)
library(ggplot2)
library(pracma)
library(PEcAn.priors)

load_local <- function(file) {
  menv <- new.env(parent = baseenv())
  load(file, envir = menv)
  as.list(menv)
}

#####################################################################################################
outdir <- "~/Documents/R/edr-da/data/"
srcdir <- "/home/femeunier/Documents/data/gigante"
Niter <- 25000
select = "Marvin" # Marvin Foster Kalacska Sanchez

use_leaf_PDA <- TRUE
use_meta.analysis <- TRUE
use_prior <- TRUE

#####################################################################################################

df_PFT <- data.frame(names = c("Tree_optical","Liana_optical"),
                     nums = c(1,2),
                     colors = c("darkgreen","darkblue"))

# Forest inventory
site.file <- c(file.path(srcdir,"Gigante_control.lat9.000lon-79.000.css"),
               file.path(srcdir,"Gigante_removal.lat9.000lon-79.000.css"))

names <- c("Control","Removal")
patch.select <- c(1,3)
census <- data.frame()

for (isite in seq_along(site.file)){
  site_dat <- PEcAn.ED2::read_css(site.file[isite]) %>% filter(patch == patch.select[isite]) 
  
  dbh <- site_dat[["dbh"]]
  pft <-  as.factor(site_dat[["pft"]])
  nplant <- site_dat[["n"]]
  patch <- site_dat[["patch"]]
  
  is_liana <- (pft == 17)
  npft <- length(levels(pft))
  levels(pft) <- 1:npft
  
  href <- 61.7;  b1Ht <- 0.035; b2Ht <- 0.69;
  h <- pmin(35,(href*(1 -exp(-b1Ht*(dbh**b2Ht)))))
  href <- 61.7;  b1Ht <- 0.11; b2Ht <- 0.87;
  h[is_liana] <- pmin(35,(href*(1 -exp(-b1Ht*(dbh[is_liana]**b2Ht)))))
  h[dbh>3 & is_liana] <- max(h) + 0.5
  temp <- sort(h, decreasing = TRUE, index.return = TRUE)
  
  # reorganize by height
  h <- temp$x
  dbh <- dbh[temp$ix]
  pft <- pft[temp$ix]
  nplant <- nplant[temp$ix]
  patch <- patch[temp$ix]
  is_liana <- is_liana[temp$ix]
  
  patch <- data.frame(dbh = dbh, h = h, pft = pft, nplant = nplant, patch = patch, type = names[isite],is_liana = is_liana)
  patch_simple <- merge_cohorts(patch)
  patch.order <- patch_simple[seq(dim(patch_simple)[1],1),]
  census <- rbind(census,
                  patch.order)
}

# Merge cohorts

npft <- length(unique(census$pft))
# patch.num <- nrow(census %>% group_by(type,patch) %>% add_count() %>% summarise(n = mean(n)))
# ncohort <- sum(census %>% group_by(type,patch) %>% add_count() %>% summarise(n = mean(n)) %>% pull(n))

run_ED_RTM <- function(params_model,patch,waves){
  # Define constants
  ## direct_sky_frac <- 0.9
  ## czen <- 1
  # Pull out site-specific params
  nplant <- patch %>% pull(nplant)
  ncohort <- length(nplant)
  dbh <- patch %>% pull(dbh)
  pft <- patch %>% pull(pft)
  
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
  sla <- pft_params[3, pft]
  bleaf <- size2bl(dbh, b1leaf, b2leaf)
  lai <- nplant * bleaf * sla
  
  wai <-rep(0, ncohort)
  cai <- rep(1, ncohort)
  
  # Extract remaining parameters
  N <- pft_params[4, ]
  Cab <- pft_params[5, ]
  Car <- pft_params[6, ]
  Cw <- pft_params[7, ]
  Cm <- pft_params[8, ]
  orient_factor <- pft_params[9, ]
  clumping_factor <- pft_params[10, ]
  
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
create_likelihood <- function(observed, patches) {
  
  function(params){
    
    ssigma <- params[1]
    params_model <- params[-c(1)]
    
    patch <- patches %>% filter(type == "Control")
    waves_C <- observed %>% filter(type == "Control") %>% pull(waves)
    albedo_C <- run_ED_RTM(params_model,patch,waves = waves_C)
    
    patch <- patches %>% filter(type == "Removal")
    waves_R <- observed %>% filter(type == "Removal") %>% pull(waves)
    albedo_R <- run_ED_RTM(params_model,patch,waves = waves_R)
    
    if (is.null(albedo_C) | is.null(albedo_R)) return(-1e20)
    if (any(!is.finite(albedo_C)) | any(!is.finite(albedo_R))) return(-1e20)
    if (any(albedo_C < 0) | any(albedo_R < 0)) albedo_C(-1e20)
    if (any(albedo_C > 1) | any(albedo_R > 1)) return(-1e20)
    
    
    observed_C <- observed %>% filter(type == "Control") %>% pull(reflectance)
    observed_R <- observed %>% filter(type == "Removal") %>% pull(reflectance)
    
    # Calculate likelihood
    ll <- sum(dnorm(albedo_C, observed_C, ssigma, log = TRUE)) +
      sum(dnorm(albedo_R, observed_R, ssigma, log = TRUE))
    ll
  }
}

#######################################################################################
# Observations
observation <- readRDS("./data/All_canopy_spectra.rds") %>% filter(wavelength > 400, wavelength < 2500)

observed <-
  observation %>% filter(ref == select) %>% ungroup() %>% rename(type = scenario) %>% mutate(type = case_when(type == "low" ~ "Removal",
                                                                                                          type == "high" ~ "Control")) %>% dplyr::select(type, reflectance,wavelength) %>%
  rename(waves = wavelength)

#######################################################################################
# Priors
prospect_PDA_file <- file.path(outdir,"PDA_all_leaf.RDS")
if (file.exists(prospect_PDA_file) & use_leaf_PDA) {PDA_results <- readRDS(prospect_PDA_file)} else {PDA_results <- NULL}

# Import Posterior/Optimized distributions
dis2find <- c('b1Bl_large','b2Bl_large','SLA','Nlayers',
              'Cab','Car','Cw','Cm','orient_factor','clumping_factor')
dis2find_prim <- c('b1Bl','b2Bl','SLA','Nlayers',
                   'Cab','Car','Cw','Cm','orient_factor','clumping_factor')

pft_lowers <- c(b1Bl_large = 0.001, b2Bl_large = 1 , SLA = 1, orient_factor = -0.5,
                clumping_factor = 0.4, b1Ca = 0.1, b2Ca = 0.5, Nlayers = 1,
                Cab = 0, Car = 0,Cw = 0,Cm = 0.)
pft_uppers <- c(b1Bl_large = 0.08, b2Bl_large = 2.5, SLA = 100, orient_factor = 0.5,
                clumping_factor = 0.9, b1Ca = 2, b2Ca = 3, Nlayers = 5,
                Ca = 150, Car = 50,Cw = 0.1,Cm = 0.1)

N <- 10000
Ndist <- length(dis2find)
sampler_all <- matrix(NA,N,(Ndist)*npft)
lower_all <- upper_all <- best_all <- rep(NA,(Ndist)*npft)

for (ipft in seq(npft)){
  current_pft <- as.character(df_PFT$names[ipft])
  
  postfile <- file.path(outdir,"pft",current_pft,"post.distns.Rdata")
  priorfile <- file.path(outdir,"pft",current_pft,"prior.distns.Rdata")
  
  if (file.exists(prospect_PDA_file) & use_leaf_PDA) {
    names_dis_PDA <- names(BayesianTools::MAP(PDA_results[[current_pft]])$parametersMAP)
  } else {
    names_dis_PDA <- NULL
  }
  
  if (file.exists(postfile) & use_meta.analysis){
    load(postfile)
    distns <- post.distns
    Distributions <- rownames(post.distns)
  } else if (file.exists(priorfile) & use_prior){
    load(priorfile)
    distns <- prior.distns
    Distributions <- rownames(prior.distns)
  } else { 
    distns <- NULL
    Distributions <- NULL
  }
  
  sampler <- matrix(NA,N,length(dis2find))
  lower <- upper <- best <- rep(NA,length(dis2find))
  
  for (idis in seq(dis2find)){
    current_dis <- dis2find[idis]
    if (current_dis %in% names_dis_PDA){
      pos <- which(current_dis == names_dis_PDA)
      sampling <- as.vector(BayesianTools::getSample(PDA_results[[current_pft]],whichParameters=pos,numSamples = N))[1:N]
    } else {
      pos <- which(Distributions == current_dis)
      
      if (isempty(pos)){
        unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
        sampling <- get.sample(unif_distribution,n=N)
        
      } else {
        sampling <- get.sample(distns[pos,],n=N)
      }
    }
    
    sampler[,idis] <- sampling
    lower[idis] <- min(sampler[,idis])
    upper[idis] <- max(sampler[,idis])
    best[idis]  <- median(sampler[,idis])
  }
  
  pos <- ((ipft-1)*(Ndist)+1):(((ipft)*(Ndist)))
  
  sampler_all[,pos]<-sampler
  lower_all[pos]<-lower
  upper_all[pos]<-upper
  best_all[pos]<-best
  
  names(lower_all)[pos] <- names(upper_all)[pos] <- names(best_all)[pos] <- dis2find_prim
}

colnames(sampler_all) <- rep(dis2find_prim,npft)

ssigma_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
soil_moisture_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
direct_sky_frac_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
czen_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)

lower_all <- c(0.001,0.001,0.001,0.001,lower_all)
upper_all <- c(1,1,1,1,upper_all)
best_all <- c(0.5,0.5,0.5,0.5,best_all)
sampler_all <- cbind(ssigma_sampling,soil_moisture_sampling,direct_sky_frac_sampling,czen_sampling,
                     sampler_all)

prior <- BayesianTools::createPriorDensity(sampler_all, method = "multivariate", eps = 1e-10,
                            lower = lower_all, upper = upper_all, best = best_all)

###################################################################################
patches <- census
likelihood <- create_likelihood(observed = observed,patches = census)

## debug(edr_r)
## likelihood(prior$sampler())

# We test first
params <- (upper_all+lower_all)/2
time0 <- Sys.time()
test <- likelihood(params = lower_all)
Sys.time() - time0

print(test)

# Settings
settings_MCMC <- BayesianTools::applySettingsDefault(settings = NULL, sampler = "DEzs", check = FALSE)
settings_MCMC$iterations=Niter

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
samples <- BayesianTools::runMCMC(samples,settings = settings_MCMC)

saveRDS(samples, file.path(getwd(),"outputs","Marvin_edr.rds"))

# BayesianTools::gelmanDiagnostics(samples)
# summary(samples)

# Results
waves <- 400:2500
best_params <- BayesianTools::MAP(samples)$parametersMAP
best_params_model <- best_params[-c(1)]
names(best_params_model) <- names(params[-c(1)])

patch <- census %>% filter(type == "Control")
albedo_C <- run_ED_RTM(params_model = best_params_model,patch,waves = waves)

patch <- census %>% filter(type == "Removal")
albedo_R <- run_ED_RTM(best_params_model,patch,waves = waves)

best_run <- as.data.frame(rbind(data.frame(waves = waves,reflectance = albedo_C,type = "Control"),
                                data.frame(waves = waves,reflectance = albedo_R,type = "Removal")))



ggplot() +
  geom_point(data = observed,aes(x = waves, y = reflectance, color = type)) +
  geom_line(data = best_run,aes(x = waves, y = reflectance, color = type)) +
  theme_bw()
