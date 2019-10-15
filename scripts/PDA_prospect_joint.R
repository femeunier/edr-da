

best_set<- best_run <- list()

# Uniform prior
Prospect_param_names <- c("Ca","Cb","Car","Cw","Cm","Nlayers","ssigma")
pft_lowers <- c(chlorophyll_a = 0,chlorophyll_b = 0, carotenoids = 0,Cw = 0, SLA = 1, Nlayers = 1, ssigma = 0)
pft_uppers <-  c(chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1, SLA = 100, Nlayers = 5, ssigma = 1)

npft <- length(df_PFT$names)

# Read Leaf spectra 


#############################################################################################
# Loop over PFTs
dis2find <- c('chlorophyll_a','chlorophyll_b','carotenoids','Cw','SLA','Nlayers',"ssigma")

N=10000
Ndist <- length(dis2find)
sampler_all <- matrix(NA,N,Ndist*npft)
lower_all <- upper_all <- best_all <- rep(NA,Ndist*npft)

for (ipft in seq(npft)){
  current_pft <- as.character(df_PFT$names[ipft])
  
  # From posterior
  load(file.path(settings$outdir,"pft",current_pft,"post.distns.Rdata"))
  Distributions <- rownames(post.distns)
  
  observation[as.character(observation$pft)==current_pft,"pft_num"]=ipft
  
  sampler <- matrix(NA,N,Ndist)
  lower <- upper <- best <- rep(NA,Ndist)
  
  for (idis in seq(dis2find)){
    current_dis = dis2find[idis]
    pos <- which(Distributions == current_dis)
    
    if (isempty(pos)){
      unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
      sampling <- get.sample(unif_distribution,n=N)
      
    } else {
      sampling <- get.sample(post.distns[pos,],n=N)
    }
    
    if (current_dis == 'SLA'){
      sampling = 1/sampling/10*2 # unit conversions
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
}

prior <- createPriorDensity(sampler_all, method = "multivariate", eps = 1e-10,
                            lower = lower_all, upper = upper_all, best = best_all)
  

# Define likelihood
create_likelihood <- function(observed_all,npft) {
  function(params){
  
      simulated_all <- observed_all_values <- c()
      observed <- reflectance_waves <- waves <- list()
      
      for (ipft in seq(1,npft)){
        
        observed[[ipft]] <- observed_all %>% filter(pft_num==ipft) %>% select(reflectance) %>% pull()
        waves[[ipft]]    <- observed_all %>% filter(pft_num==ipft) %>% select(wavelength) %>% pull()
        
        # PROSPECT parameters
        Ca <- params[(ipft-1)*7+1]
        Cb <- params[(ipft-1)*7+2]
        Cab <- Ca + Cb
        Car <- params[(ipft-1)*7+3]
        Cw <- params[(ipft-1)*7+4]
        Cm <- params[(ipft-1)*7+5]
        Nlayers <- params[(ipft-1)*7+6]
        
        # print(Cm)
        
        ssigma <- params[(ipft-1)*7+7]
        
        optical_param <- c(Nlayers,Cab,Car,Cw,Cm)
        names(default) <- c('Nlayers','chlab','carotenoids','Cw','Cm')
        
        # Call RTM
        result <- tryCatch(
          PEcAnRTM::prospect(optical_param, version = "5"),
          error = function(e) NULL)
        if (is.null(result)) return(-1e20)
        reflectance <- result[,"reflectance"]
        if (any(!is.finite(reflectance))) return(-1e20)
        if (any(reflectance < 0)) return(-1e20)
        if (any(reflectance > 1)) return(-1e20)
        
        reflectance_waves[[ipft]] <- interp1(x = PEcAnRTM::wavelengths(reflectance),
                                     y = as.vector(matrix(reflectance)),
                                     waves[[ipft]])
        simulated_all <- c(simulated_all,reflectance_waves[[ipft]])
        observed_all_values <- c(observed_all_values,observed[[ipft]])
      }
      
      # Add differences
      for (ipft in seq(1,npft)){
        if (ipft<npft){
          for (ipft_prim in seq(ipft+1,npft)){
            waves_interp <- waves[[ipft]]
            
            reflectance_waves_interp <- approxExtrap(x = waves[[ipft_prim]],
                                                y = reflectance_waves[[ipft_prim]],
                                                xout = waves_interp)
            
            observed_interp   <- approxExtrap(x = waves[[ipft_prim]],
                                                y = observed[[ipft]],
                                                xout = waves_interp)
            
            simulated_all <- c(simulated_all,reflectance_waves_interp$y-reflectance_waves[[ipft]])
            
            observed_all_values <- c(observed_all_values,observed_interp$y-observed[[ipft]])
          }
        }
      }

      # Calculate likelihood
      ll <- sum(dnorm(simulated_all, observed_all_values, ssigma, log = TRUE))
      return(ll)
  }
}
  
likelihood <- create_likelihood(observation, npft)
# likelihood(rep(c(50,20,10,0.01,0.02,2,0.1),2))

# Run inversion
setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
samples <- BayesianTools::runMCMC(setup)
samples <- BayesianTools::runMCMC(samples)

samples$setup$names <- paste(rep(Prospect_param_names,npft),rep(df_PFT$names,7),sep='_')
BayesianTools::gelmanDiagnostics(samples)
summary(samples)


# correlationPlot(samples)
# plot(samples)
# marginalPlot(samples)

posteriorMat <- getSample(samples, parametersOnly = TRUE)
posteriorMat_df <- melt(posteriorMat) %>%select(Var2,value) %>% rename(param = Var2)
marginalPlot <- 
  ggplot(posteriorMat_df) +
  geom_density_ridges_gradient(aes(x = value, y = as.factor(param), fill = ..x..),
                               alpha = 0.5,rel_min_height=0.01,scale = 3) + 
  theme_minimal() +
  scale_fill_viridis()

MAP_samples <- MAP(samples)$parametersMAP

plot(x=NA,y=NA,xlim=c(0,2500),ylim=c(0,0.6),ylab= "Leaf reflectance",xlab="Wavelength")

for (ipft in seq(npft)){
  
  current_pft <- as.character(df_PFT$names[ipft])
  
  best_set[[current_pft]]<- c(MAP_samples[paste('Nlayers',current_pft,sep='_')],
                              MAP_samples[paste('Ca',current_pft,sep='_')]+MAP_samples[paste('Cb',current_pft,sep='_')],
                              MAP_samples[paste('Car',current_pft,sep='_')],
                              MAP_samples[paste('Cw',current_pft,sep='_')],
                              MAP_samples[paste('Cm',current_pft,sep='_')])
  best_run[[current_pft]] <- PEcAnRTM::prospect(best_set[[current_pft]], version = "5")
  
  temp <-Spectrum_leaf_data %>% filter(pft==current_pft) %>% select(c('wavelength','reflectance'))
  
  observation <- as.vector(temp$reflectance)
  waves <- temp$wavelength
  
  C <- Colors[ipft]
  lines(PEcAnRTM::wavelengths(best_run[[current_pft]]),best_run[[current_pft]][,1],col=C,type='l')
  lines(waves,observation,type='p',col=C,pch=19,lwd=0.1)
}


