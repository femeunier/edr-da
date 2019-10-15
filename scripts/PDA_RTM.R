# Argument functions
use_meta.analysis <- FALSE
use_leaf_PDA <- FALSE

alpha_frac <- 0.8

PFTselect = 17

crown_mod <- 1

# Data 
observation <- Spectrum_canopy_data

npft <- length(df_PFT$names)

# Test params
# params <- c(0.1,0.5,0.9,1,
#             rep(c(0.02, 1.847, 22, 0.01,2, 0.7, 1.12,2,40,20,15,0.01,0,0.7),npft))

h5file <- "/home/carya/output/PEcAn_99000000002/out/SA-median/RTM/history-S-2004-01-01-120000-g01.h5"

#########################################################################################################
# Function
#########################################################################################################
# hdf5 file reading

hfile <- hdf5r::H5File$new(h5file)

dbh <- readDataSet(hfile[["DBH"]])
nplant <- readDataSet(hfile[["NPLANT"]])
hite <- readDataSet(hfile[["HITE"]])
pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
PACO_N <- readDataSet(hfile[["PACO_N"]])

hfile$close_all()

PFTselect <- match(PFTselect,df_PFT$PFTnum)

#########################################################################################################
# likelihood function
create_likelihood <- function(observed,dbh,nplant,hite,pft,PACO_N,PFTselect,crown_mod) {
  
  Ncohort <- length(dbh)
  npft <- max(as.numeric(pft))
  
  Npatch <- length(PACO_N)
  PAN <- rep(1:Npatch,PACO_N)
  
  stopifnot(
    length(pft) == Ncohort,
    length(dbh) == Ncohort,
    length(nplant) == Ncohort,
    length(hite) == Ncohort
  )
  
  function(params){
    
    # print("here")
    # Model Parameters
    ssigma <- params[1]
    soil_moisture <- params[2]
    direct_sky_frac <- params[3]
    czen <- params[4]
    
    pft_params <- matrix(params[-(1:4)], ncol = npft)
    
    b1leaf <- pft_params[1, pft]
    b2leaf <- pft_params[2, pft]
    sla <- pft_params[3, pft]
    b1wood <- pft_params[4, pft]
    b2wood <- pft_params[5, pft]
    b1Ca <- pft_params[6, pft]
    b2Ca <- pft_params[7, pft]
    N <- pft_params[8, ]
    Ca <- pft_params[9, ]
    Cb <- pft_params[10, ]
    Car <- pft_params[11, ]
    Cw <- pft_params[12, ]
    orient_factor <- pft_params[13, ]
    clumping_factor <- pft_params[14, ]
    
    Cm <- 1/pft_params[3,]/10
    Cab <- Ca + Cb
    # Calculate allometries
    bleaf <- redr::size2bl(dbh, b1leaf, b2leaf)
    lai <- nplant * bleaf * sla
    wai <- wai_allometry(dbh, nplant, b1wood, b2wood)
    
    if (crown_mod == 0){
      cai <- rep(1,Ncohort)
    } else {
      cai <- cai_allometry(dbh,nplant,b1leaf,b2leaf,sla,b1Ca,b2Ca)
    }
    
    # Compute Crown occupancy index
    COI <-rep(NA,Npatch)
    
    for (ipa in seq(1,Npatch)){
      
      pos <- (PAN == ipa)
      
      hiteconow   <- hite[pos]
      dbhconow    <- dbh[pos]
      pftconow    <- pft[pos]
      laiconow    <- lai[pos]
      nplantconow <- nplant[pos]
      CAconow <- cai[pos]
    
      # Compute COI
      LAItot <- sum(laiconow)
      LAI_alpha <- (1-alpha_frac)*LAItot
      h_alpha <- interp1(x = cumsum(laiconow), y = hiteconow, xi = ifelse(LAI_alpha<laiconow[1],laiconow[1],LAI_alpha))
      
      COI[ipa] <- max(0,min(1,sum(laiconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/LAI_alpha))
    }
    
    # classifify them
    patch_class <- list()
    patch_class[["low"]] <- which(COI < 0.25)
    patch_class[["high"]] <- which(COI > 0.5)
    # crash if one is empty
    
    if (any(sapply(patch_class,isempty))) return(-1e20)
    
    patch2simulate <- as.vector(unlist(patch_class))
    Npatch_simulated <- length(patch2simulate)
    
    scenarios <- as.character(unique(observed[["scenario"]]))
    
    ll <- rep(NA,length(scenarios))
    
    # ED-RTM
    plot(NA,NA,xlim=c(300,2500),ylim=c(0,0.6))
    C<-c('black','blue')
    
    for (iscenar in seq(scenarios)){
      temp <- observed %>% filter(scenario == scenarios[iscenar])
      
      # Observations
      waves <- temp %>% dplyr::select(wavelength) %>% pull()
      observed_Reflectance <- temp %>% dplyr::select(reflectance) %>% pull()
        
      # Simulations
      patches_scenar <- patch_class[[scenarios[iscenar]]]
      Npatch2simulate <- length(patches_scenar)
      
      # Simulation output
      output_RTM <- matrix(NA,length(c(par.wl,nir.wl)),Npatch2simulate)
      
      for (ipa in seq(1,Npatch2simulate)){
        
        pos <- (PAN == patches_scenar[ipa])
        
        pftconow    <- pft[pos]
        laiconow    <- lai[pos]*20
        waiconow    <- lai[pos]
        caiconow    <- cai[pos]
        
        output_RTM[,ipa] <- (edr_r(
          pft=pftconow,
          lai= laiconow,
          wai= waiconow, 
          cai = caiconow,
                N, Cab, Car, Cw, Cm,
                orient_factor, clumping_factor,
                soil_moisture,
                direct_sky_frac,
                czen,
                wavelengths = c(par.wl,nir.wl)))[["albedo"]]
      } # patch loop
      
      if (any(!is.finite(output_RTM))) return(-1e20)
      if (any(output_RTM < 0)) return(-1e20)
      if (any(output_RTM > 1)) return(-1e20)
      
      simulated_reflectance <- apply(output_RTM,1,median)
      simulated_reflectance_Waves <- approxExtrap(x = c(par.wl,nir.wl),
                                                  y = simulated_reflectance,
                                                  xout = waves)
      # lines(waves,simulated_reflectance_Waves$y,lty=1,col=C[iscenar]) # solid line = simulated
      # lines(waves,observed_Reflectance,lty=2,col=C[iscenar])
      
      ll[iscenar] <- sum(dnorm(simulated_reflectance_Waves$y, observed_Reflectance, ssigma, log = TRUE))
    } # scenario loop
    return(sum(ll))
  } # end function
} # end create_likelihood

###################################################################################
# Define priors
## prior_mvtraits <- load_local("priors/mvtraits_priors.RData")
prospect_PDA_file <- file.path(settings$outdir,"pfts","all","Leaf_spectra_PDA.RDS")
if (file.exists(prospect_PDA_file) & use_leaf_PDA) {PDA_results <- load_rds(prospect_PDA_file)} else {PDA_results <- NULL}


dis2find <- c('b1Bl_large','b2Bl_large','SLA','b1wood','b2wood','b1Ca','b2Ca',
              'Nlayers','chlorophyll_a','chlorophyll_b','carotenoids','Cw',
              'orient_factor','clumping_factor')

pft_lowers <- c(b1Bl_large = 0.001, b2Bl_large = 1 , SLA = 1, b1wood = 0.001, b2wood = -5, b1Ca = 0.1, b2Ca = 1,
                Nlayers = 1,chlorophyll_a = 1,chlorophyll_b = 1, carotenoids = 1,Cw = 0.001, orient_factor = -0.5,
                clumping_factor = 0.001)
pft_uppers <- c(b1Bl_large = 0.08, b2Bl_large = 2.5, SLA = 100, b1wood = 3, b2wood = -2,  b1Ca = 1.2, b2Ca = 2,
                Nlayers = 5, chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1, orient_factor = 0.5,
                clumping_factor = 1)


# Import Posterior/Optimized distributions
N <- 10000
Ndist <- length(dis2find)
sampler_all <- matrix(NA,N,Ndist*npft)
lower_all <- upper_all <- best_all <- rep(NA,Ndist*npft)

for (ipft in seq(npft)){
  current_pft <- as.character(df_PFT$names[ipft])
  
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
    current_dis <- dis2find[idis]
    if (current_dis %in% PDA_results[[current_pft]]$setup$names){
      pos <- which(current_dis == PDA_results[[current_pft]]$setup$names)
      sampling <- as.vector(getSample(PDA_results[[current_pft]],whichParameters=pos,numSamples = N))[1:N]
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
  
  pos <- ((ipft-1)*Ndist+1):(((ipft)*Ndist))
  sampler_all[,pos]<-sampler
  lower_all[pos]<-lower
  upper_all[pos]<-upper
  best_all[pos]<-best
  
  names(lower_all)[pos] <- 
    names(upper_all)[pos] <- names(best_all)[pos] <- dis2find
}
colnames(sampler_all) <- rep(dis2find,npft)

ssigma_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
soil_moisture_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
direct_sky_frac_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
czen_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)

lower_all <- c(0.001,0.001,0.001,0.001,lower_all)
upper_all <- c(1,1,1,1,upper_all)
best_all <- c(0.5,0.5,0.5,0.5,best_all)
sampler_all <- cbind(ssigma_sampling,soil_moisture_sampling,direct_sky_frac_sampling,czen_sampling,
                     sampler_all)

prior <- createPriorDensity(sampler_all, method = "multivariate", eps = 1e-10,
                            lower = lower_all, upper = upper_all, best = best_all)

###################################################################################
# Likelihood function

likelihood <- create_likelihood(observation,dbh,nplant,hite,pft,PACO_N,PFTselect,crown_mod) 
# Test likelihood function
# params <- sampler_all[20,]
# likelihood(params)

###################################################################################
# Run inversion

iter <- 100
settings_MCMC <- list(iterations = iter, consoleUpdates = 1)

setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
samples <- BayesianTools::runMCMC(samples,settings = settings_MCMC)
summary(samples)

samples$setup$names<-c("ssigma","soil.moisture","direct_sky_frac","czen",rep(dis2find,npft))
MAP_samples <- MAP(samples)$parametersMAP
params <- MAP_samples

par(mar=c(5,4,4,2))
par(mar=c(0,0,0,0))
samples_small <- getSample(samples,parametersOnly=FALSE,whichParameters = 1:3,coda = TRUE,includesProbabilities=TRUE)
correlationPlot(samples_small)

plot(samples_small)
marginalPlot(samples_small, prior = samples$setup$prior)

## Local Variables:
## ess-r-package--project-cache: (redr . "/Users/shik544/Box Sync/Projects/edr_pda/edr-da/")
## End:
