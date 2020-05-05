## library(redr)
devtools::load_all(".")

# library(devtools)
# install_github("ashiklom/pecan", subdir = "modules/rtm",force = TRUE,dependencies = FALSE)
# install_github("PecanProject/pecan", subdir = "models/ed",force = TRUE,dependencies = FALSE)
# devtools::install_github("ashiklom/rrtm",dependencies = FALSE)

library(PEcAnRTM)
library(rrtm)
library(dplyr)
library(ggplot2)

load_local <- function(file) {
  menv <- new.env(parent = baseenv())
  load(file, envir = menv)
  as.list(menv)
}


# Forest inventory
site.file <- "/home/femeunier/Documents/data/BCI/ForestGEO.2012.lat9.1543lon-79.8461.css"
site_dat <- PEcAn.ED2::read_css(site.file) %>% filter(patch == 1) %>% top_n(n = 20)

dbh <- site_dat[["dbh"]]
pft <- site_dat[["pft"]] - max( site_dat[["pft"]]) + 1
nplant <- site_dat[["n"]]
npft <- length(unique(pft))
ncohort <- length(dbh)

# True set of parameters
ref_parameters <- c( "soil_moisture"= 0.5,
                     "direct_sky_frac" = 0.5,
                     "czen" = 0.5,
                     "b1leaf"=0.1,
                     "b2leaf"=1.6,
                     "sla"=12,
                     "b1wood"=0,
                     "b2wood" =0,
                     "N"=2,
                     "Cab"=50,
                     "Car"=20,
                     "Cw"=0.01,
                     "Cm"=0.01,
                     "orient_factor"= 0,
                     "clumping_factor"=0.7)

run_ED_RTM <- function(params_model){
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
    # b1wood, b2wood,
    # N, Cab, Car, Cw, Cm,
    # clumping_factor, orient_factor
    pft_params <- matrix(params_model[-(1:3)], ncol = npft)
    
    # Calculate allometries
    b1leaf <- pft_params[1, pft]
    b2leaf <- pft_params[2, pft]
    sla <- pft_params[3, pft]
    bleaf <- size2bl(dbh, b1leaf, b2leaf)
    lai <- nplant * bleaf * sla
    
    b1wood <- pft_params[4, pft]
    b2wood <- pft_params[5, pft]
    wai <- wai_allometry(dbh, nplant, b1wood, b2wood)
    cai <- rep(1, ncohort)
    
    # Extract remaining parameters
    N <- pft_params[6, ]
    Cab <- pft_params[7, ]
    Car <- pft_params[8, ]
    Cw <- pft_params[9, ]
    Cm <- pft_params[10, ]
    orient_factor <- pft_params[11, ]
    clumping_factor <- pft_params[12, ]
    
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
   
   ssigma <- params[1]
   params_model <- params[-c(1)]
   albedo <- run_ED_RTM(params_model)
     
   if (is.null(albedo)) return(-1e20)
   if (any(!is.finite(albedo))) return(-1e20)
   if (any(albedo < 0)) return(-1e20)
   if (any(albedo > 1)) return(-1e20)
   
   # Calculate likelihood
   ll <- sum(dnorm(albedo, observed, ssigma, log = TRUE))
   ll
 }
}

# Create observation
waves = 400:2500
observation <- run_ED_RTM(params = ref_parameters)
plot(waves,observation,type = 'l')


# Define priors
## prior_mvtraits <- load_local("priors/mvtraits_priors.RData")

pft_lowers <- rep(c(
  b1leaf = 0, b2leaf = -5, sla = 0,
  b1wood = 0, b2wood = -5,
  N = 1, Cab = 0, Car = 0, Cw = 0, Cm = 0,
  orient_factor = -0.5, clumping_factor = 0.4
), npft)
pft_uppers <- rep(c(
  b1leaf = 3, b2leaf = 5, sla = 100,
  b1wood = 3, b2wood = 5,
  N = 5, Cab = 100, Car = 40, Cw = 0.1, Cm = 0.1,
  orient_factor = 0.5, clumping_factor = 1
), npft)
lowers <- c(
  "ssigma" = 0, "soil_moisture" = 0, "direct_sky_frac" = 0, "czen" = 0,
  pft_lowers
)
uppers <- c(
  "ssigma" = 1, "soil_moisture" = 1, "direct_sky_frac" = 1, "czen" = 1,
  pft_uppers
)
prior <- BayesianTools::createUniformPrior(lower = lowers, upper = uppers)
likelihood <- create_likelihood(observation, aviris_use_wl, pft, dbh, nplant)

## debug(edr_r)
## likelihood(prior$sampler())

# We test first
params <- (uppers+lowers)/2
time0 <- Sys.time()
test <- likelihood(params)
Sys.time() - time0

print(test)


# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup)
samples <- BayesianTools::runMCMC(samples)
BayesianTools::gelmanDiagnostics(samples)
summary(samples)

# Results
best_params <- BayesianTools::MAP(samples)$parametersMAP
best_params_model <- best_params[-c(1)]
names(best_params_model) <- names(ref_parameters)
albedo <- run_ED_RTM(best_params_model)

plot(log(best_params_model))
lines(log(ref_parameters),type='p',col='red')
plot(ref_parameters,best_params_model,log='xy')

data.opt <- rbind(as.data.frame(cbind(waves,albedo)) %>% mutate(type = "best run"),
                  as.data.frame(cbind(waves,as.vector(observation))) %>% rename(albedo = V2) %>% mutate(type = "obs"))

ggplot() +
  geom_line(data = data.opt %>% filter(type == "best run"),aes(waves,albedo)) +
  geom_point(data = data.opt %>% filter(type == "obs"),aes(waves,albedo)) +
  theme_bw()


ggplot(data = data.opt %>% group_by(waves) %>% summarise(diff = albedo[type == "best run"]-albedo[type == "obs"],
                                                         obs = albedo[type == "obs"])) +
  geom_point(aes(waves,diff)) +
  theme_bw()

ggplot(data = data.opt %>% group_by(waves) %>% summarise(diff = albedo[type == "best run"]-albedo[type == "obs"],
                                                         obs = albedo[type == "obs"])) +
  geom_point(aes(waves,100*diff/obs)) +
  theme_bw()

## Local Variables:
## ess-r-package--project-cache: (redr . "/Users/shik544/Box Sync/Projects/edr_pda/edr-da/")
## End:
