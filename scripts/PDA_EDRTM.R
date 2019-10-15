# Argument function
use_meta.analysis <- FALSE
use_leaf_PDA <- FALSE

# Data 
observation <- Spectrum_canopy_data

alpha_frac <- 0.8
PFTselect <- 17
crown_mod <- 0

h5file <- "/home/carya/output/PEcAn_99000000002/out/SA-median/RTM/history-S-2004-01-01-120000-g01.h5"

npft <- length(df_PFT$names)

rundir <- "/home/carya/output/PEcAn_99000000002/run/SA-median"
outdir <- file.path(settings$modeloutdir,"PDA_EDRTM")
outdir <- "/home/carya/output/PEcAn_99000000002/out/PDA_EDRTM"

ed2in <- PEcAn.ED2::read_ed2in(file.path(rundir,"ED2IN"))
ed2in$FRQFAST <- 86400
ed2in$IFOUTPUT = 0
ed2in$IDOUTPUT = 0
ed2in$IMOUTPUT = 0
ed2in$IQOUTPUT = 0
ed2in$IYOUTPUT = 0
ed2in$ITOUTPUT = 0
ed2in$IOOUTPUT = 0
ed2in$ISOUTPUT = 0

ed2in$DTLSM = 900
ed2in$RADFRQ = 900

edr_ed2in <- PEcAnRTM::setup_edr(ed2in, outdir, date,TRUE)

ED2IN_file <- file.path(outdir,"ED2IN")
ed2in <- PEcAn.ED2::read_ed2in(ED2IN_file)
h5file <- paste0(paste(ed2in$SFILIN,"S",ed2in$IYEARH,sprintf('%02d',ed2in$IMONTHH),
                       sprintf('%02d',ed2in$IDATEH),paste0(ed2in$ITIMEH,"00"),sep='-'),"-g01.h5")

hfile <- hdf5r::H5File$new(h5file)

dbh <- readDataSet(hfile[["DBH"]])
nplant <- readDataSet(hfile[["NPLANT"]])
Npatch <- readDataSet(hfile[["NPATCHES_GLOBAL"]])
hite <- readDataSet(hfile[["HITE"]])
pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
PACO_N <- readDataSet(hfile[["PACO_N"]])

hfile$close_all()

PFTselect <- match(PFTselect,df_PFT$PFTnum)

Ncohort <- length(dbh)
npft <- max(as.numeric(pft))

Npatch <- length(PACO_N)
PAN <- rep(1:Npatch,PACO_N)

inventory <- list(dbh = dbh,pft = pft,hite = hite,nplant = nplant, PAN = PAN,PFTselect = PFTselect,
                  PACO_N = PACO_N, Ncohort = Ncohort, Npatch = Npatch)

if(!dir.exists(outdir)) dir.create(outdir)

##############################################################################
# Define likelihood
create_likelihood <- function(observed,inventory,crown_mod,rundir,outdir,plot = FALSE){
  
  function(params) {
      
    ssigma <- params[1]
    
    outputs <- run_ED_RTM(rundir,outdir,params[-1],crown_mod,inventory)
       
    COI <- outputs[["COI"]]
    output_RTM <- outputs[["output_RTM"]]
   
    # classifify them
    patch_class <- list()
    patch_class[["low"]] <- which(COI < 0.25)
    patch_class[["high"]] <- which(COI > 0.5)
    
    if (any(sapply(patch_class,isempty))) return(-1e20)
    
    # Get outputs
    scenarios <- as.character(unique(observed[["scenario"]]))
    
    ll <- rep(NA,length(scenarios))
    
    if (plot) {
      plot(NA,NA,xlim=c(300,2500),ylim=c(0,0.6))
      C<-c('black','blue')
    }

    for (iscenar in seq(scenarios)){
      temp <- observed %>% filter(scenario == scenarios[[iscenar]])
      
      waves <- temp %>% dplyr::select(wavelength) %>% pull()
      observed_Reflectance <- temp %>% dplyr::select(reflectance) %>% pull()
      
      if (length(patch_class[[scenarios[[iscenar]]]])>1){
        simulated_reflectance <- apply(output_RTM[patch_class[[scenarios[[iscenar]]]],],2,median)
      } else {
        simulated_reflectance <- output_RTM[patch_class[[scenarios[[iscenar]]]],]
      }
      simulated_reflectance_Waves <- approxExtrap(x = c(par.wl,nir.wl),
                                                  y = simulated_reflectance,
                                                  xout = waves)
      if (plot){
        lines(waves,simulated_reflectance_Waves$y,lty=1,col=C[iscenar]) # solid line = simulated
        lines(waves,observed_Reflectance,lty=2,col=C[iscenar]) # solid line = simulated
      }
      # lines(observed_Reflectance,simulated_reflectance_Waves$y,type='p',col=C[iscenar])
      
      # Calculate likelihood
      ll[iscenar] <- sum(dnorm(simulated_reflectance_Waves$y, observed_Reflectance, ssigma, log = TRUE))
    }  
    return(sum(ll))
  }
} # end likewihood

##############################################################################
# Define priors

if (crown_mod == 0){
  dis2find <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                'chlorophyll_a','chlorophyll_b','carotenoids','Cw')
  dis2find_prim <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                     'Cab','carotenoids','Cw')
} else {
  dis2find <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','b1Ca','b2Ca','Nlayers',
                'chlorophyll_a','chlorophyll_b','carotenoids','Cw')
  dis2find_prim <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','b1Ca','b2Ca','Nlayers',
                     'Cab','carotenoids','Cw')
}


pft_lowers <- c(b1Bl_large = 0.001, b2Bl_large = 1 , SLA = 1, orient_factor = -0.5,
                clumping_factor = 0.4, b1Ca = 0.1, b2Ca = 1, Nlayers = 1,
                chlorophyll_a = 0,chlorophyll_b = 0, carotenoids = 0,Cw = 0)
pft_uppers <- c(b1Bl_large = 0.08, b2Bl_large = 2.5, SLA = 100, orient_factor = 0.5,
                 clumping_factor = 0.9, b1Ca = 1.2, b2Ca = 2, Nlayers = 5,
                 chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1)

prospect_PDA_file <- file.path(settings$outdir,"pfts","all","Leaf_spectra_PDA.RDS")
if (file.exists(prospect_PDA_file) & use_leaf_PDA) {PDA_results <- load_rds(prospect_PDA_file)} else {PDA_results <- NULL}


# Import Posterior/Optimized distributions
N <- 10000
Ndist <- length(dis2find)
sampler_all <- matrix(NA,N,(Ndist-1)*npft)
lower_all <- upper_all <- best_all <- rep(NA,(Ndist-1)*npft)

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
  
  # Merge Ca and Cb
  colnames(sampler) <- dis2find
  names(lower) <- names(upper) <- names(best) <- dis2find
  
  sampler[,"chlorophyll_a"] <- sampler[,"chlorophyll_a"] + sampler[,"chlorophyll_b"]
  sampler <- sampler[,!(colnames(sampler) %in% c("chlorophyll_b"))]
  colnames(sampler)[(colnames(sampler) %in% c("chlorophyll_a"))] <- "Cab"
  
  lower["chlorophyll_a"] <- lower["chlorophyll_a"] + lower["chlorophyll_b"]
  lower <- lower[!(names(lower) %in% c("chlorophyll_b"))]
  names(lower)[(names(lower) %in% c("chlorophyll_a"))] <- "Cab"
  
  upper["chlorophyll_a"] <- upper["chlorophyll_a"] + upper["chlorophyll_b"]
  upper <- upper[!(names(upper) %in% c("chlorophyll_b"))]
  names(upper)[(names(upper) %in% c("chlorophyll_a"))] <- "Cab"
  
  best["chlorophyll_a"] <- best["chlorophyll_a"] + best["chlorophyll_b"]
  best <- best[!(names(best) %in% c("chlorophyll_b"))]
  names(best)[(names(best) %in% c("chlorophyll_a"))] <- "Cab"
  
  pos <- ((ipft-1)*(Ndist-1)+1):(((ipft)*(Ndist-1)))
  
  sampler_all[,pos]<-sampler
  lower_all[pos]<-lower
  upper_all[pos]<-upper
  best_all[pos]<-best
  
  names(lower_all)[pos] <- names(upper_all)[pos] <- names(best_all)[pos] <- dis2find_prim
}

colnames(sampler_all) <- rep(dis2find_prim,npft)

ssigma_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)
soil_moisture_sampling <- get.sample(list(distn = "unif",parama = 0,paramb = 1),n=N)

lower_all <- c(0.001,0.001,lower_all)
upper_all <- c(1,1,upper_all)
best_all <- c(0.5,0.5,best_all)
sampler_all <- cbind(ssigma_sampling,soil_moisture_sampling,sampler_all)

prior <- createPriorDensity(sampler_all, method = "multivariate", eps = 1e-10,
                            lower = lower_all, upper = upper_all, best = best_all)

##############################################################################
# likelihood function

likelihood <- create_likelihood(observation,inventory,crown_mod,rundir,outdir,plot = FALSE)

# Test likelihood function
params<-sampler_all[1,]
time0 <- Sys.time()
likelihood(params)
Sys.time() - time0

# Run inversion
iter <- 10000
settings_MCMC <- list(iterations = iter, consoleUpdates = 1)

setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
samples <- BayesianTools::runMCMC(samples,settings = settings_MCMC)
BayesianTools::gelmanDiagnostics(samples)

samples$setup$names<-c("ssigma","soil.moisture",rep(dis2find_prim,npft))
summary(samples)

MAP_samples <- MAP(samples)$parametersMAP
params <- MAP_samples
likelihood <- create_likelihood(observation,inventory,crown_mod,rundir,outdir,plot = TRUE)
likelihood(params)

# saveRDS(samples,file.path(outdir,"samples.RDS"))
# samples <- load_rds(file.path(outdir,"samples.RDS"))