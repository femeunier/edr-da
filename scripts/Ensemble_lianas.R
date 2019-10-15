#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#

# Library
library(redr)
# library(PEcAn.all)
library(rlist)
library(hdf5r)
library(ggplot2)
library(viridis)
library(ggridges)
library(pracma)
library(dplyr)
library(tidyr)
library(reshape2)
library(purrr)
library(zoo)
library(stringr)
library(BayesianTools)
library(Hmisc)

# Inputs
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou_PFT3.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou_PFT3_CM0.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_BCI_PFT3.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_BCI_PFT3_CM0.xml"
pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Gigante.xml"

edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt_2") 
rerun <- TRUE       # Rerun ED2 simulation workflow?
run_all <- FALSE    # run all simulations or just the first one (median)?
rerun_RTM <- TRUE   # Rerun RTM simulations?

use_meta.analysis <- FALSE # use meta.analysis posteriors for PDA?
use_leaf_PDA <- FALSE

par.wl = 400:2499 
nir.wl = 2500

wl.min <- 400    # Min. of observations to be considered 
wl.max <- 2500   # Max. of observations to be considered

crown_mod = 0    # 0 or 1

alpha = 0.05        # Confidence intervals

ispatch <- TRUE     # All patches or just the last one?

Npatch_max <- 0
Ncohort_max <- 200

Colors <- c("#137300","#1E64C8")

# PDA parameters
nrChains <- 4
nrIter <- 10000
Nrun_prospect <- 500

############################################################################################################
# Output
output_RTM <- list ()

settings <- PEcAn.settings::read.settings(pecanxmlfile)

# PFT table
df_PFT <- data.frame(PFTnum = as.numeric(unlist(lapply(lapply(settings$pfts, "[[","constants"),"[[","num"))),
                     names = as.character(unlist(lapply(settings$pfts, "[[","name"))))
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

npft <- length(df_PFT$names)
mapdf <- data.frame(pfts=df_PFT$names,col=Colors)


# Lambdas()
Lambdas <- c(par.wl,nir.wl)



############################################################################################################
############################################################################################################
# Step 1
############################################################################################################
############################################################################################################

settings$model$ed2in_tags[["MAXPATCH"]] <- Npatch_max
settings$model$ed2in_tags[["MAXCOHORT"]] <- Ncohort_max

############################################################################################################
# Run PEcAn initial workflow simulation
if (rerun){
  # Clean output folders
  system(paste("rm -rf",paste0(settings$outdir,"/*")))
  redr::run_workflow_pecan(settings,run_all)
}

# Output directories
if (!dir.exists(file.path(settings$outdir,"pfts"))){
  dir.create(file.path(settings$outdir,"pfts"))
  dir.create(file.path(settings$outdir,"pfts","all"))
}

############################################################################################################
# Read updated settings
outdir <- settings$outdir
pecanxml <- file.path(outdir,"pecan.CONFIGS.xml")
settings <- PEcAn.settings::read.settings(pecanxml)

# Inputs
start.date <- settings$run$start.date
date <- ISOdatetime(
  lubridate::year(start.date),
  lubridate::month(start.date),
  lubridate::mday(start.date),
  12, 00, 00, tz = "UTC")

runfile <- file.path(settings$rundir,"runs.txt")
runs <- readLines(runfile)

############################################################################################################
# Meta-analysis optical traits
# Prospect version 5
# Version 5 N, Cab, Car, Cw, Cm
# Default optical parameters
default <- c(1.4,73.2,15.6,0.016,0.0105)
names(default) <- c('Nlayers','chlab','carotenoids','Cw','Cm')
optical_param_default <- trait_values_default <- list()

for (i in seq(df_PFT$PFTnum)){
  PFT_name <- as.character(df_PFT$names[i])
  optical_param_default[[PFT_name]] <- default
  trait_values_default[[PFT_name]] <- list()
}

#######################################################################
# Prospect runs

spectra_list_all <- trait_values_all <- list()

for (irun in seq(runs)){
  
  # directories
  rundir <- file.path(settings$rundir,runs[irun])
  outdir <- file.path(settings$outdir,"out",runs[irun])

  # read_config
  optical_param <- optical_param_default
  trait_values  <- trait_values_default
  
  config_file <- file.path(rundir,"config.xml")
  xml <- XML::xmlParse(config_file)
  config <- XML::xmlToList(xml)
  as.numeric(unlist(lapply(config, "[[","num")))
  
  config_temp <- config[["pft"]]
  while (!is.null(config_temp)){
    pft_num <- config_temp$num
    pft_name <- as.character(df_PFT$names[df_PFT$PFTnum==pft_num])
    
    if (!is.null(config_temp$chlorophyll_a) & !is.null(config_temp$chlorophyll_a)){
      optical_param[[pft_name]]['chlab'] <- as.numeric(config_temp$chlorophyll_a) + as.numeric(config_temp$chlorophyll_b) 
    }
    
    if (!is.null(config_temp$carotenoids)){
      optical_param[[pft_name]]['carotenoids']  <- as.numeric(config_temp$carotenoids)
    }
    
    if (!is.null(config_temp$Cw)){
      optical_param[[pft_name]]['Cw']  <- as.numeric(config_temp$Cw)
    }
    
    if (!is.null(config_temp$Nlayers)){
      optical_param[[pft_name]]['Nlayers']  <- as.numeric(config_temp$Nlayers)
    }
    
    if (!is.null(config_temp$SLA)){
      optical_param[[pft_name]]['Cm']  <- 1/as.numeric(config_temp$SLA)/10*2
    }
    
    if ("SLA" %in% names(config_temp)){
      config_temp$SLA <- as.character(as.numeric(config_temp$SLA)*0.48) # To go back to PEcAn values
    }
    if ("Vcmax" %in% names(config_temp)){
      config_temp$Vcmax <- as.character(as.numeric(config_temp$Vcmax)*2.4) # To go back to PEcAn values
    }
    
    trait_values[[pft_name]] <- config_temp
    
    config[["pft"]] <- NULL
    config_temp <- config[["pft"]]
  }
  
  trait_values_all[[runs[irun]]] <- trait_values
  
  spectra_list <- list()
  
  for (i in seq(df_PFT$PFTnum)){
    PFT_name <- as.character(df_PFT$names[i])
    spectra_list[[PFT_name]] <-  PEcAnRTM::prospect(optical_param[[PFT_name]], version = "5")   
    spectra_list[[PFT_name]][,] <- na.approx(spectra_list[[PFT_name]][,],rule =2)
  }

  spectra_list_all[[runs[irun]]] <- spectra_list

}

saveRDS(spectra_list_all,file = file.path(settings$outdir,"pfts","all","Leaf_spectra.RDS"))

##################################################################################
# Boxplots traits
sample_file <- file.path(settings$outdir,"samples.Rdata")
load(sample_file)
PEcAn_traits <- ensemble.samples
PEcAn_traits$env <- NULL

names_col <- sapply(ensemble.samples,names)
common_traits <- names(which(table(unlist(names_col[df_PFT$names]))==length(df_PFT$names)))

data_traits <- data.frame()
p_values <- rep(NA,length(common_traits))
names(p_values) <- common_traits

for (trait in common_traits){
  for (pft in df_PFT$names){
    value <- PEcAn_traits[[pft]][,trait]
    data_traits <- rbind(data_traits,data.frame(traits = trait,pfts = pft,values = value))
  }
  
  data_aov <- data_traits %>% filter(traits == trait)
  AoV <- summary(aov(data = data_aov, values ~ pfts))
  p_values[trait] <- AoV[[1]][1,5]
}

data_traits <- data_traits %>% mutate(p_values = p_values[match(data_traits$traits,names(p_values) )]) %>%
  mutate(signif = case_when(
    p_values <0.001 ~ '***',
    p_values <0.01 ~ '**',
    p_values <0.05 ~ '*',
    TRUE ~ ''))

data_traits <- data_traits %>% mutate(Trait = traits) %>%  unite("traits_signif",c(traits,signif),sep='')

plot_traits <- ggplot(data_traits, aes(x=factor(pfts), y=values,fill=factor(pfts))) + 
  geom_boxplot() + 
  scale_fill_manual(values = df_PFT$Col) +
  theme_bw()+
  facet_wrap(~traits_signif,scales = "free")

ggsave(filename = file.path(settings$outdir,"pfts","all","traits_ensemble.png"),dpi = 300,
       plot = plot_traits,
       width = 10, height = 8)

traits2remove <- c('b1Ca','b2Ca','clumping_factor','orient_factor',
                   'orient_factor_shifted','b1Bl','b2Bl',
                   'b1Bl_large','b2Bl_large','b1Bl_small','b2Bl_small')

plot_traits_select <- ggplot(data_traits %>% filter(!(Trait %in% traits2remove)), aes(x=factor(pfts), y=values,fill=factor(pfts))) + 
  geom_boxplot() + 
  scale_fill_manual(values = df_PFT$Col) +
  theme_bw() +
  facet_wrap(~traits_signif,scales = "free")

ggsave(filename = file.path(settings$outdir,"pfts","all","traits_ensemble_select.png"),dpi = 300,
       plot = plot_traits_select,
       width = 10, height = 8)

##################################################################################
# Leaf spectra

ens.runs <- runs[unlist(map(.x = runs,.f = function(run){
  runtype <- readLines(file.path(settings$outdir,"run",run, "README.txt"),n = 1)
  return((grepl("ensemble", runtype)) )
}),use.names = FALSE)]

sa.runs <- runs[!(runs %in% ens.runs)]

spectra_in <- spectra_list_all[ens.runs]
data_spectra <- data.frame()

for (run in ens.runs){
  for (pft in df_PFT$names){
    spectrum_temp <- spectra_in[[run]][[pft]]
    data_spectra <- rbind(data_spectra,data.frame(runs = run,
                                                  pfts = pft,
                                                  reflectance = spectrum_temp[,1],
                                                  transmittance = spectrum_temp[,2],
                                                  Lambda = Lambdas))
  }
}


data_spectra_select <- data_spectra %>% group_by(pfts,Lambda) %>% summarise(rmin = min(reflectance,na.rm = TRUE),
                                                                         rmax = max(reflectance,na.rm = TRUE),
                                                                         alphamin = quantile(reflectance,alpha/2,na.rm = TRUE),
                                                                         alphamax = quantile(reflectance,1-alpha/2,na.rm = TRUE),
                                                                         median = median(reflectance,na.rm = TRUE))
# Load leaf spectra data
Spectrum_liana_data <- load_rds("~/data/RTM/Spectrum_liana_data.R")

Spectrum_leaf_data <-
  rbind(cbind(as.data.frame(Spectrum_liana_data[[1]]),pft=as.character(df_PFT$names[2])),
        cbind(as.data.frame(Spectrum_liana_data[[2]]),pft=as.character(df_PFT$names[1])))

# Leaf spectra ensemble
leaf_spectra <- ggplot(data_spectra_select,
                       aes(x = Lambda,y = median,colour = pfts)) +
  geom_line(size = 1, alpha = 0.5,linetype = 2) +
  geom_ribbon(aes(ymin = alphamin, ymax = alphamax,fill=pfts),alpha = 0.5, size=0.5, linetype=0) + 
  theme_bw() +
  scale_color_manual(values = df_PFT$Col) +
  scale_fill_manual(values = df_PFT$Col) +
  xlab('Wavelength (nm)') + 
  ylab('Leaf reflectance (-)') +
  theme(panel.grid=element_blank())+
  geom_line(data = Spectrum_leaf_data, aes(x = wavelength, y = reflectance,colour = as.factor(pft)),size=1)

ggsave(filename = file.path(settings$outdir,"pfts","all","leaf_spectra_MA.png"),dpi = 300,
       plot = leaf_spectra,
       width = 10, height = 5)

############################################################################################################
# PDA

PDA_results <- PDA_prospect(settings,Spectrum_leaf_data,df_PFT,wl.min,wl.max,use_meta.analysis,
                            nrChains=nrChains,nrIter=nrIter)
saveRDS(PDA_results,file = file.path(settings$outdir,"pfts","all","Leaf_spectra_PDA.RDS"))

# Plot results of PDA

best_set <- best_run <- list()

rsquare <- rmse <- rep(NA,npft)
prospect_performance <- ensemble_posterior <- data.frame()
posterior_dis <- data.frame()

for (ipft in seq(npft)){
  
  current_pft <- as.character(df_PFT$names[ipft])
  
  temp <-Spectrum_leaf_data %>% filter(pft==current_pft & wavelength>wl.min & wavelength<wl.max) %>% select(c('wavelength','reflectance'))
  
  observation <- as.vector(temp$reflectance)
  waves <- temp$wavelength
   
  samples <- PDA_results[[current_pft]]

  MAP_samples <- MAP(samples)$parametersMAP
  best_set[[current_pft]]<- c(MAP_samples['Nlayers'],MAP_samples['Ca']+MAP_samples['Cb'],MAP_samples['Car'],MAP_samples['Cw'],MAP_samples['Cm'])
  best_run[[current_pft]] <- PEcAnRTM::prospect(best_set[[current_pft]], version = "5")
  
  best_run_interp <- interp1(x=PEcAnRTM::wavelengths(best_run[[current_pft]]),
                             y = as.vector(best_run[[current_pft]][,1]),
                             xi = waves)
  
  C <- Colors[ipft]
  
  current_model <- data.frame(obs = observation, mod = best_run_interp , pft = current_pft)
  
  lm_pft <- lm(data = current_model,formula = mod ~ obs)
  rsquare[ipft] <- summary(lm_pft)$adj.r.squared
  rmse[ipft] <- sqrt(c(crossprod(lm_pft$residuals))/(length(lm_pft$residuals)-1))
  coef_pft <- coef(lm_pft)

  prospect_performance <- rbind(prospect_performance,current_model)
  
  posteriorMat <- getSample(samples, numSamples = Nrun_prospect,
                                           parametersOnly = TRUE)
  
  posteriorMat_rel <- posteriorMat
  posteriorMat_rel[,"Ca"] <-  posteriorMat_rel[,"Ca"] +   posteriorMat_rel[,"Cb"]
  colnames(posteriorMat_rel)[colnames(posteriorMat_rel)=="Ca"] = "Cab"
  posteriorMat_rel <- posteriorMat_rel[,-which(colnames(posteriorMat_rel)=="Cb")]
  
  posterior_dis <- rbind(posterior_dis, 
                         (melt(posteriorMat_rel) %>% rename(Param = Var2) %>% select(c(Param,value)) %>% mutate(pft = current_pft)))
  
  temp_ensemble_prospect <- matrix(NA,nrow(posteriorMat),nrow(best_run[[current_pft]]))
  colnames(temp_ensemble_prospect) <- c(par.wl,nir.wl)
  for (irun in seq( nrow(posteriorMat))){
    current_parameter_set <- c(posteriorMat[irun,'Nlayers'],posteriorMat[irun,'Ca']+posteriorMat[irun,'Cb'],posteriorMat[irun,'Car'],posteriorMat[irun,'Cw'],posteriorMat[irun,'Cm'])
    current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
    temp_ensemble_prospect[irun,] <- current_model_output[,1]
  }
 
  ensemble_posterior <- rbind(ensemble_posterior,
   as.data.frame(melt(temp_ensemble_prospect) %>% select(Var2,value) %>%rename(wavelength = Var2, reflectance = value) %>% 
    group_by(wavelength) %>% summarise(rmin = min(reflectance,na.rm=TRUE),
                                    rmax = max(reflectance,na.rm=TRUE),
                                    alphamin = quantile(reflectance,alpha/2,na.rm=TRUE),
                                    alphamax = quantile(reflectance,1-alpha/2,na.rm=TRUE),
                                    median = median(reflectance,na.rm=TRUE))) %>% mutate(pft = current_pft))
  
}

performance_prospect_plot <-
  ggplot(prospect_performance,aes(x=obs,y=mod,colour = as.factor(pft),fill = as.factor(pft))) +
  geom_smooth(method = "lm",level = 1-alpha) + 
  scale_color_manual(values = df_PFT$Col) +
  scale_fill_manual(values = df_PFT$Col) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,colour = 'black',linetype=2) + 
  theme_bw()

ggsave(filename = file.path(settings$outdir,"pfts","all","Performance_PDA.png"),dpi = 300,
       width = 10, height = 5, plot = performance_prospect_plot)

spectra_post <- ggplot(ensemble_posterior,
                  aes(x = wavelength,y = median, colour = as.factor(pft))) +
  geom_line(size = 1, alpha = 0.5,linetype = 2) +
  geom_ribbon(aes(ymin = alphamin, ymax = alphamax,fill=as.factor(pft)),alpha = 0.5, size=0.5, linetype=0) + 
  theme_bw() +
  scale_color_manual(values = df_PFT$Col) +
  scale_fill_manual(values = df_PFT$Col) +
  xlab('Wavelength (nm)') + 
  ylab('Leaf reflectance (-)') +
  theme(panel.grid=element_blank()) + 
  geom_point(data = Spectrum_leaf_data,aes(x = wavelength,y=reflectance,colour=as.factor(pft)))

ggsave(filename = file.path(settings$outdir,"pfts","all","leaf_spectra_PDA.png"),dpi = 300,
       width = 10, height = 5, plot = spectra_post)


# Plot of marginal distirbutions
# Remove 5 of extreme
posterior_dis_95 <- posterior_dis %>% group_by(Param,pft) %>% filter(value <= quantile(value,0.975) & value >= quantile(value,0.025))

marginalPlots <-
  ggplot(posterior_dis %>% group_by(Param) %>% mutate(rel_value = (value - min(value))/(max(value)-min(value))), 
         aes(x = rel_value, y = Param, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_fill_manual(values = df_PFT$Col) +
  labs(y="") +
  theme_ridges(font_size = 13) + 
  theme_bw()

ggsave(filename = file.path(settings$outdir,"pfts","all","leaf_marginalPlots.png"),dpi = 300,
       width = 10, height = 5, plot = marginalPlots)

############################################################################################################
############################################################################################################
# Step 2
############################################################################################################
############################################################################################################

############################################################################################################
# RTM runs

nruns <- length(runs)
pb <- txtProgressBar(min = 0, max = nruns, style = 3, title = "RTM simulations")
pbi <- 0

for (irun in seq(runs)){
  
  # Porgress bar
  setTxtProgressBar(pb, pbi)
  
  #######################################################################
  # Directories
  rundir <- file.path(settings$rundir,runs[irun])
  outdir <- file.path(settings$outdir,"out",runs[irun])
    
  output_dir <- file.path(outdir,"RTM")
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  if (rerun_RTM){  
    ED2IN_file <- file.path(rundir,"ED2IN")
    ed2in <- PEcAn.ED2::read_ed2in(ED2IN_file)
    
    if (!run_all){
      if (irun == 1){
        ed2in_sfilout_ref <- ed2in$SFILOUT
      } else {
        ed2in$SFILOUT <- ed2in_sfilout_ref 
      }
    }
    history.path <- dirname(ed2in$SFILOUT)
    ed2in$MAXPATCH <- 0
    ed2in$MAXCOHORT <- 0
    
    # crown mod config !
    # ed2in$CROWN_MOD <- crown_mod
    PFT_ref <- ed2in$INCLUDE_THESE_PFT
    
    
    #######################################################################
    # RTM run
    ed2in$IBRANCH_THERMO=1
    
    edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, date,patches=ispatch)
    
    soil_reflect <- hapke_soil(rep(rand(1),2101))
    soil_file <- file.path("/tmp","soil_reflect_par.dat")
    writeLines(text=paste(as.character(soil_reflect),collapse=" "),con = soil_file)
    
    # wood_reflect_path_default <- system.file("extdata", 
    #                                 "wood_reflect_par.dat", package = "PEcAnRTM")
    # wood_reflect_path <- scan(wood_reflect_path_default)
    
    output_RTM[[runs[irun]]] <- 
      PEcAnRTM::EDR(img_path = NULL,
                    ed2in_path = edr_ed2in,
                    spectra_list = spectra_list_all[[runs[irun]]],
                    trait.values = trait_values_all[[runs[irun]]],
                    soil_reflect_path = soil_file,
                    edr_exe_path =  edr_exe_path,
                    par.wl = par.wl, 
                    nir.wl = nir.wl,
                    patches = ispatch)
    
    # PEcAn standard
    model2netcdf.EDR(outdir,
                     sitelat = settings$run$site$lat,
                     sitelon = settings$run$site$lon,
                     start_date = start.date,
                     par.wl = par.wl, 
                     nir.wl = nir.wl,
                     patches = ispatch) 
    
    # Plot patch structure
    ED2IN_file <- file.path(settings$outdir,"out",runs[irun],"RTM","ED2IN")
    ed2in <- PEcAn.ED2::read_ed2in(ED2IN_file)
    h5file <- paste0(paste(ed2in$SFILIN,"S",ed2in$IYEARH,sprintf('%02d',ed2in$IMONTHH),
                           sprintf('%02d',ed2in$IDATEH),paste0(ed2in$ITIMEH,"00"),sep='-'),"-g01.h5")
    plot_patch_vertical_distribution(h5file,ispatch)
    fraction_LAI(h5file,ispatch,PFTselect = 17,alpha_frac = 0.8)
    
    # fractionLAI<-as.numeric(readLines(file.path(output_dir,"LAI_fraction.dat")))
    # COI<-as.numeric(readLines(file.path(output_dir,"COI.dat")))
    # COI2<-as.numeric(readLines(file.path(output_dir,"COI2.dat")))
    
  } else{
    output_RTM[[runs[irun]]] <- read.EDR.output(output_dir,par.wl,
                                                nir.wl, ispatch)
  }
  
  pbi <- pbi + 1
}

settings$sensitivity.analysis$ensemble.id <- NULL
settings$ensemble$ensemble.id <- NULL
settings$sensitivity.analysis$perpft <- TRUE

saveRDS(output_RTM,file = file.path(settings$outdir,"pfts","all","Canopy_spectra.RDS"))

##################################################################################
# Plotting results
# VD plot all pfts 

variables <- c('All_leaf','PAR_leaf','NIR_leaf',                
               'Spectrum','Vis','Red','PAR','NIR') # Leaf spectrum and Canopy spectrum!
parameter2remove <- c(rep(list(traits2remove),3),
                      rep(list(c("")),5))

for (ivar in seq(variables)){
  VDP_allPFTs(variables[ivar],mapdf = mapdf,parameter2remove[[ivar]])
}

#####################################################################################
# Ensemble runs

# Data canopy spectra
Spectrum_canopy_data <- load_rds("~/data/RTM/Spectrum_canopy_data.R")
Spectrum_canopy_data <-
  rbind(cbind(as.data.frame(Spectrum_canopy_data[[2]]),scenario='low'),
        cbind(as.data.frame(Spectrum_canopy_data[[1]]),scenario='high'))


#########################################################"
#2 change!!!
Liana_fraction <- as.numeric(readLines(file.path(settings$modeloutdir,runs[1],"RTM","LAI_fraction.dat")))
scenarios <- c(0,0.1,0.5)
scenar_name <- c('low','high')

data_spectra_canopy <- data.frame()
for (iscenar in seq(scenar_name)){
  for (run in ens.runs){
    if (ispatch){
      pos <- which(Liana_fraction > scenarios[iscenar] & Liana_fraction < scenarios[iscenar+1])
      canopy_spectrum_temp <- output_RTM[[run]][pos,]
      rownames(canopy_spectrum_temp) <- pos
      colnames(canopy_spectrum_temp) <- c(par.wl,nir.wl)
    } else {
      pos <- 1
      canopy_spectrum_temp <- output_RTM[[run]]
      rownames(canopy_spectrum_temp) <- pos
      colnames(canopy_spectrum_temp) <- c(par.wl,nir.wl)
    }
    
    canopy_spectrum_temp <- melt(canopy_spectrum_temp) %>% rename(Patch = Var1, Lambda = Var2, Reflectance = value)
    
    canopy_spectrum_temp <- canopy_spectrum_temp %>% ungroup() %>%
      group_by(Lambda) %>% summarise(rmin = min(Reflectance,na.rm=TRUE),
                                     rmax = max(Reflectance,na.rm=TRUE),
                                     alphamin = quantile(Reflectance,alpha/2,na.rm=TRUE),
                                     alphamax = quantile(Reflectance,1-alpha/2,na.rm=TRUE),
                                     median = median(Reflectance,na.rm=TRUE))
    
    data_spectra_canopy <- rbind(data_spectra_canopy,
                                 data.frame(runs = run,
                                            reflectance = c(canopy_spectrum_temp$median,canopy_spectrum_temp$rmin,canopy_spectrum_temp$rmax),
                                            scenario = scenar_name[iscenar],
                                            Lambda = rep(canopy_spectrum_temp$Lambda,3)))
  }
}

data_spectra_canopy_select <- data_spectra_canopy %>% 
  group_by(Lambda,scenario) %>% summarise(rmin = min(reflectance),
                                       rmax = max(reflectance),
                                       alphamin = quantile(reflectance,alpha/2),
                                       alphamax = quantile(reflectance,1-alpha/2),
                                       median = median(reflectance))

# Ecosystem spectra ensemble
spectra <- ggplot(data_spectra_canopy_select,
       aes(x = Lambda,y = median, colour = as.factor(scenario))) +
  geom_line(size = 1, alpha = 0.5,linetype = 2) +
  geom_ribbon(aes(ymin = rmin, ymax = rmax,fill=as.factor(scenario)),alpha = 0.5, size=0.5, linetype=0) + 
  theme_bw() +
  scale_color_manual(values = Colors) +
  scale_fill_manual(values = Colors) +
  xlab('Wavelength (nm)') + 
  ylab('Canopy reflectance (-)') +
  theme(panel.grid=element_blank()) + 
  geom_point(data = Spectrum_canopy_data,aes(x = wavelength,y=reflectance,colour=as.factor(scenario)))
  
ggsave(filename = file.path(settings$outdir,"pfts","all","spectra.png"),dpi = 300,
       width = 10, height = 5, plot = spectra)

################################################################################
# Canopy/Leaf spectra sensitivity analysis

load(file.path(settings$outdir,"samples.Rdata"))
sa.samples[["env"]] <- NULL
pft.names <- names(sa.samples)

pecandir <- settings$outdir
outdir <- settings$modeloutdir

data_sa <- data.frame()
for (pft.name in pft.names){
  traits <- colnames(sa.samples[[pft.name]])
  quantiles <- rownames(sa.samples[[pft.name]])

  start.year <- settings$sensitivity.analysis$start.year
  end.year <- settings$sensitivity.analysis$end.year  
  
  sa.run.ids <- runs.samples$sa
  
  data_pft <- read.sa.output.spectrum(traits, quantiles, pecandir, outdir, variable = "Spectrum", sa.run.ids, par.wl,nir.wl)
  
  data_pft_leaf <- read.sa.output.leafspectrum(traits, quantiles, pecandir, outdir, variable, sa.run.ids,par.wl,nir.wl,spectra_list_all,pft.name)
  
  for (ipatch in seq(unique(data_pft$patch))){
    patch_plot <- ggplot(data_pft %>% filter(patch == ipatch), aes(x = lambda, y = Reflectance,group = quantiles,linetype = quantiles)) +
      geom_line(size = 0.2, color = mapdf$col[mapdf$pfts == pft.name]) +
      theme_bw() +
      facet_wrap(~traits) +
      theme(panel.grid=element_blank())
    
    path_patches <- file.path(settings$outdir,"pft",pft.name,"patches")
    if (!dir.exists(path_patches)) dir.create(path_patches)
    ggsave(filename = file.path(path_patches,paste0("spectra",ipatch,".png")),dpi = 300,
           width = 10, height = 5, plot = patch_plot)
    
    
  }
  
  patch_plot_leaf <- ggplot(data_pft_leaf %>% filter(var == 'reflectance'), aes(x = lambda, y = Reflectance,group = quantiles,linetype = quantiles)) +
    geom_line(size = 0.2, color = mapdf$col[mapdf$pfts == pft.name]) +
    theme_bw() +
    facet_wrap(~traits) +
    theme(panel.grid=element_blank())
  
  ggsave(filename = file.path(path_patches,paste0("leaf_spectra",".png")),dpi = 300,
         width = 10, height = 5, plot = patch_plot_leaf)
}  


#####################################################################################
# GSA

tmppath <- settings$outdir
PFT <- as.character(df_PFT$names)
vars<- c('NIR','Spectrum','Vis','Red')
mapdf <- data.frame(pfts=c(PFT),col=Colors)

plot_GSA(tmppath,vars,mapdf)

