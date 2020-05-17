rm(list = ls())

library(pracma)
library(dplyr)
library(PEcAnRTM)
library(ggplot2)
library(tidyr)
library(BayesianTools)
library(purrr)
library(ggridges)
library(ggpointdensity)

wv.min = 400
wv.max = 2500

directory <- "/home/carya/data/RTM/data/Guzman_all"
metadata.file <- file.path(directory,"meta_data_all2.csv")
metadata <- read.csv(metadata.file) %>% mutate(GF = case_when(substr(Code,1,1) == "L" ~ "Liana",
                                                              substr(Code,1,1) == "T" ~ "Tree"),
                                               Species = (substring(sub("\\_.*", "", Code),2)),
                                               Ind = as.numeric(sub(".*\\_", "", Code)),
                                               N  = 1 + To.ASD.spectra - From.ASD.spectra,
                                               data.exist = FALSE,
                                               site = case_when(substr(Site.date,1,1) == "P" ~ "PNM",
                                                                substr(Site.date,1,1) == "F" ~ "FTS"))

# metadata %>% group_by(GF) %>% summarise(N = sum(N))

files <- c('PNM050307.csv','PNM060307.csv','FS070307.csv')
data.files <- list()

for (ifile in seq(files)){
  data.files[[tools::file_path_sans_ext(files[ifile])]] <- read.csv(file.path(directory,files[ifile])) 
}

ncol(data.files[[1]])-1+ncol(data.files[[2]])-1+ncol(data.files[[3]])-1
All.cols <- unique(c(colnames(data.files[[1]]),colnames(data.files[[2]]),colnames(data.files[[3]])))

data.all <- list()
data.all[["Liana"]] <- list()
data.all[["Tree"]] <- list()

df.data.all <- data.frame()

for (i in seq(1,nrow(metadata))){
  
  print(i/nrow(metadata))
  
  temp <- metadata[i,]
  wv <- data.files[[as.character(temp[["Site.date"]])]][["Wavelength"]]
  columns <- paste0(as.character(temp$ASD.code..prefix.),".",sprintf("%03d",seq(temp$From.ASD.spectra,temp$To.ASD.spectra)))
  if (any(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))){
    
    columns <- columns[which(columns %in% colnames(data.files[[as.character(temp[["Site.date"]])]]))]

    Reflectance <- data.files[[as.character(temp[["Site.date"]])]][,columns]
    metadata[["data.exist"]][i] <- TRUE
    
    temp.mat <- cbind(wv,Reflectance)
    colnames(temp.mat) <- c("wv",1,2,3)
    df.single <- temp.mat %>% pivot_longer(c("1","2","3")) %>% filter(wv >= wv.min,wv <= wv.max)
   
    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] <- list() 
    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] [["spectrum"]] <- df.single
    data.all[[metadata[["GF"]][i]]] [[metadata[["Species"]][i]]] [[metadata[["Ind"]][i]]] [["site"]] <- metadata[["site"]][i]
    
    df.single_small <- df.single %>% filter(wv %in% seq(wv.min,wv.max,5))
    df.data.all <- rbind(df.data.all,
                         df.single_small %>% mutate(GF = metadata[["GF"]][i],
                                                    Species = metadata[["Species"]][i],
                                                    Ind = metadata[["Ind"]][i],
                                                    site = metadata[["site"]][i]))
    
  }
}

summary <- metadata %>% filter(data.exist) %>% ungroup() %>% group_by(GF,site) %>% summarise(Nspectra = sum(N),
                                                                                             Nspecies = length(unique(Species)),
                                                                                             max_Nleaf = length(unique(Ind)))

sum(summary$Nspectra)

df.data.all_sum <- df.data.all %>% group_by(site,GF,wv) %>% summarise(R = mean(value),
                                                                      Rmin = quantile(value,0.025),
                                                                      Rmax = quantile(value,0.975))

ggplot(data = df.data.all_sum) +
  geom_ribbon(aes(x = wv,ymin=Rmin,ymax = Rmax,fill = GF),color = NA,alpha = 0.2) +
  geom_line(aes(x = wv,y=R,color = GF)) +
  facet_wrap(~ site) +
  theme_bw()


run_prospect <- function(params){
  # PROSPECT parameters
  Nlayers <- params[1]
  Cab <- params[2]
  Car <- params[3]
  Cw <- params[4]
  Cm <- params[5]
  
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
  return(reflectance_waves)
}

# test optimization
create_likelihood <- function(observed, waves) {
  function(params) {
    
    ssigma <- params[6]
    reflectance_waves <- run_prospect(params)
    # Calculate likelihood
    ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    return(ll)
  }
}

# Uniform prior
Prospect_param_names <- c("Nlayers","Cab","Car","Cw","Cm","ssigma")
pft_lowers <- c(Nlayers = 1, Cab = 0 , Car = 0,Cw = 0, Cm = 0, ssigma = 0)
pft_uppers <-  c(Nlayers = 5, Cab = 100, Car = 50,Cw = 0.1, Cm = 0.1, ssigma = 1)

compt <- 1
for (pft in seq(names(data.all))){
  for (species in seq(names(data.all[[pft]]))){
    for (ind in seq(length(data.all[[pft]][[species]]))) {
        
      print(compt)
      data.temp <- data.all[[pft]][[species]][[ind]]
      
      if (!is.null(data.temp)){
        temp <- data.temp[["spectrum"]]
        observed <- as.vector(t(temp[,"value"]))
        waves <- as.vector(t(temp[,"wv"]))
        
        prior <- createUniformPrior(pft_lowers, pft_uppers)
        likelihood <- create_likelihood(observed=observed, waves)
        settings_MCMC <- list(iterations = 10000, nrChains = 1)
        
        # Run inversion
        setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
        samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
        samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)
        
        best_param <- BayesianTools::MAP(samples)$parametersMAP
        names(best_param) <- Prospect_param_names
        
        best_run <- run_prospect(best_param)
        
        data.all[[pft]][[species]][[ind]] [["param"]] <- best_param
        data.all[[pft]][[species]][[ind]] [["best_run"]] <- best_run
        
        LM <- lm(formula = y ~ x, data = data.frame(x = observed, y = best_run))
        data.all[[pft]][[species]][[ind]] [["r2"]] <- summary(LM)$r.squared
        data.all[[pft]][[species]][[ind]] [["rmse"]] <- sqrt(mean(LM$residuals^2))
        
        # plot(waves,observed,type='l')
        # lines(waves,best_run,col='red')
        
        print(c(summary(LM)$r.squared, 100*sqrt(mean(LM$residuals^2))/mean(observed)))
        
        compt <- compt +1
      }
    }
  }
}

# saveRDS(file = "Guzman_allfits.RDS",object = data.all)
data.all <- readRDS(file = "./scripts/Guzman_allfits.RDS")

##############################################################
GFs <- names(data.all)
stats.all <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  r2 <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'r2'))
  rmse <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'rmse'))
  site <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'site'))
  if (is.null(r2)){
    return(NULL)
    } else{
    data.frame(r2 = r2, rmse = rmse, site = site,species = i, GF = GFs[iGF])
  }
  }))
}))

ggplot(stats.all %>% pivot_longer(c(r2,rmse))) +
  geom_boxplot(aes(x = GF, y = value,fill = GF),alpha = 0.5) +
  facet_wrap(site~name,scales = "free") + 
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw()

params.all <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  params <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'param'))
  site <- unlist(sapply(data.all[[GFs[iGF]]][[i]], '[', 'site'))
  if (is.null(params)){
    return(NULL)
  } else{
    data.frame(value = params, params = sub(".*\\.", "", names(params)), site = site, species = i, GF = GFs[iGF])
  }
}))
}))

params.all2 <- params.all %>% filter(!(params == "ssigma"))
params.all2 %>% group_by(GF,params) %>% summarise(m = mean(value)) %>% arrange(params)

ggplot(data = params.all2,aes(x = value, y = GF, fill = GF)) +
  geom_density_ridges(alpha= 0.5) +
  facet_grid(site ~ params,scales = "free") +
  scale_color_manual(values = c("darkblue","darkgreen")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  theme_bw() 


wv_select <- seq(wv.min,wv.max,50)

data.all_formated <- do.call(rbind,map(1:length(GFs),function(iGF) {do.call(rbind,map(1:length(data.all[[GFs[iGF]]]),function(i){
  print(c(iGF,i))
  
  clist <- data.all[[GFs[iGF]]][[i]]
  clist[sapply(clist, is.null)] <- NULL
  waves <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,1]}))
  num <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,2]}))
  Reflectance <- unlist(lapply(sapply(clist, '[', 'spectrum'),function(x){x[,3]}))
  best_run <- unlist(sapply(clist, '[', 'best_run'))
  site <- unlist(sapply(clist, '[', 'site'))
  
  select <- which(waves %in% wv_select)
    
  if (is.null(best_run)){
    return(NULL)
  } else{
    data.frame(wv = waves[select], num = num[select], obs = Reflectance[select], site = site,sim = best_run[select], species = i ,GF = GFs[iGF])
  }
}))
}))

ggplot(data.all_formated) +
  geom_point(aes(x=obs,y=sim,color = GF, shape = site),alpha = 0.2,size = 0.5) +
  geom_abline(slope = 1, color = "black",size = 1.3) +
  theme_bw()
# 
ggplot(data.all_formated, aes(x=obs, y=sim)) +
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density") +
  facet_wrap(site~GF,scales = "free") +
  geom_abline(slope = 1,color = "black") +
  theme_bw()
#   #+ geom_point(shape = '.')

ggplot(data.all_formated, aes(x=obs, y=sim)) +
  geom_pointdensity() + scale_color_viridis_c() + 
  geom_abline(slope = 1, color = "black",size = 1.3) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(site~GF,scales = "free") +
  theme_bw()

data.all_formated %>% group_by(GF,site) %>% summarise(r2=summary(lm(formula = sim ~ obs))$r.squared)

