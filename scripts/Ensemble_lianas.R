#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#

# Library
library(redr)
library(PEcAn.all)
library(rlist)
library(hdf5r)
library(ggplot2)
library(pracma)
library(dplyr)
library(tidyr)
library(reshape2)
library(purrr)

# Inputs
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou_PFT3.xml"
pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_BCI_PFT3.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_BCI.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou.xml"

edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt") 
rerun <- TRUE   # Rerun ED2 simulation workflow?
run_all <- FALSE  # run all simulations or just the first one (median)?
rerun_RTM <- TRUE   # Rerun RTM simulations?

par.wl = 400:2499 
nir.wl = 2500

ispatch <- TRUE

# Output
output_RTM <- list ()

############################################################################################################
settings <- PEcAn.settings::read.settings(pecanxmlfile)

# PFT table
df_PFT <- data.frame(PFTnum = as.numeric(unlist(lapply(lapply(settings$pfts, "[[","constants"),"[[","num"))),names = as.character(unlist(lapply(settings$pfts, "[[","name"))))
df_PFT <- df_PFT %>% arrange(PFTnum)

# Lambdas()
Lambdas <- c(par.wl,nir.wl)
############################################################################################################
# Run ED2 simulation
if (rerun){
  # Clean output folders
  system(paste("rm -rf",paste0(settings$outdir,"/*")))
  redr::run_workflow_pecan(pecanxmlfile,run_all)
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

############################################################################################################
# RTM runs

spectra_list_all <- list()

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
  
  ED2IN_file <- file.path(rundir,"ED2IN")
  ed2in <- PEcAn.ED2::read_ed2in(ED2IN_file)
  
  history.path <- dirname(ed2in$SFILOUT)
  ed2in$MAXPATCH <- 0
  ed2in$MAXCOHORT <- 0
  
  # crown mod config !
  ed2in$CROWN_MOD <- 0
  PFT_ref <- ed2in$INCLUDE_THESE_PFT
  
  #######################################################################
  # Meta-analysis
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
    
    trait_values[[pft_name]] <- config_temp
    
    config[["pft"]] <- NULL
    config_temp <- config[["pft"]]
  }
  
  spectra_list <- list()
  
  for (i in seq(df_PFT$PFTnum)){
    PFT_name <- as.character(df_PFT$names[i])
    spectra_list[[PFT_name]] <-  PEcAnRTM::prospect(optical_param[[PFT_name]], version = "5")    
  }
  
  spectra_list_all[[runs[irun]]] <- spectra_list
    
  #######################################################################
  # RTM run
  if (rerun_RTM){
    
    edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, date,patches=ispatch)
    
    output_RTM[[runs[irun]]] <- 
      PEcAnRTM::EDR(img_path = NULL,
                    ed2in_path = edr_ed2in,
                    spectra_list = spectra_list,
                    trait.values = trait_values,
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
    
  } else{
    output_RTM[[runs[irun]]] <- read.EDR.output(output_dir,par.wl,
                                                nir.wl, ispatch)
  }
  
  pbi <- pbi + 1
}

settings$sensitivity.analysis$ensemble.id <- NULL
settings$ensemble$ensemble.id <- NULL
settings$sensitivity.analysis$perpft <- TRUE

##################################################################################
# Plotting results
# VD plot all pfts

Colors <- c("#137300","#1E64C8")
mapdf <- data.frame(pfts=c("Tree_optical","Liana_optical"),col=Colors)
variables <- c('NIR','PAR','Vis','Red','Spectrum')

for (ivar in seq(variables)){
  VDP_allPFTs(variables[ivar],mapdf = mapdf)
}

#################################################################################
# per.pft SA
settings$sensitivity.analysis$variable <- "Spectrum"
runModule.run.sensitivity.analysis(settings)

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

data_traits <- data_traits %>% unite("traits_signif",c(traits,signif),sep='')

plot_traits <- ggplot(data_traits, aes(x=factor(pfts), y=values,fill=factor(pfts))) + 
  geom_boxplot() + 
  scale_fill_manual(values = Colors) +
  facet_wrap(~traits_signif,scales = "free")

ggsave(filename = file.path(settings$outdir,"pfts","all","traits_ensemble.png"),dpi = 300,
       plot = plot_traits,
       width = 10, height = 8)

##################################################################################
# Leaf spectra

ens.runs <- runs[unlist(map(.x = runs,.f = function(run){
  runtype <- readLines(file.path(settings$outdir,"run",run, "README.txt"),n = 1)
  return((grepl("ensemble", runtype)) )
}),use.names = FALSE)]

sa.runs <- runs[!(runs %in% ens.runs)]

spectra_in <- spectra_list_all[ens.runs]
data_spectra <- data_spectra_eco <- data.frame()

for (run in ens.runs){
  for (pft in df_PFT$names){
    spectrum_temp <- spectra_in[[run]][[pft]]
    data_spectra <- rbind(data_spectra,data.frame(runs = run,
                                                  pfts = pft,
                                                  reflectance = spectrum_temp[,1],
                                                  transmittance = spectrum_temp[,2],
                                                  Lambda = Lambdas))
  }
  
  if (ispatch){
    eco_spectrum_temp <- apply(output_RTM[[run]],2,mean)
  } else {
    eco_spectrum_temp <- output_RTM[[run]]
  }

  data_spectra_eco <- rbind(data_spectra_eco,data.frame(runs = run,
                                                        reflectance = eco_spectrum_temp,
                                                        Lambda = Lambdas))
}

alpha = 0.05

data_spectra_select <- data_spectra %>% group_by(pfts,Lambda) %>% mutate(rmin = min(reflectance),
                                                                         rmax = max(reflectance),
                                                                         alphamin = quantile(reflectance,alpha/2),
                                                                         alphamax = quantile(reflectance,1-alpha/2),
                                                                         median = median(reflectance))

data_spectra_eco_select <- data_spectra_eco %>% group_by(Lambda) %>% mutate(rmin = min(reflectance),
                                                                            rmax = max(reflectance),
                                                                            alphamin = quantile(reflectance,alpha/2),
                                                                            alphamax = quantile(reflectance,1-alpha/2),
                                                                            median = median(reflectance))

# Leaf spectra ensemble
leaf_spectra <- ggplot(data_spectra_select,
       aes(x = Lambda,y = median, ymin = alphamin, ymax = alphamax,fill = pfts,colour = pfts,linetype = pfts)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.5, size=0.5, linetype=0) + 
  theme_bw() +
  scale_color_manual(values = Colors) +
  scale_fill_manual(values = Colors) +
  scale_linetype_manual(values=c("solid", "longdash")) +
  xlab('Wavelength (nm)') + 
  ylab('Leaf reflectance (-)') +
  theme(panel.grid=element_blank())

ggsave(filename = file.path(settings$outdir,"pfts","all","leaf_spectra.png"),dpi = 300,
       plot = leaf_spectra,
       width = 10, height = 5)

# Ecosystem spectra ensemble
spectra <- ggplot(data_spectra_eco_select,
       aes(x = Lambda,y = median, ymin = alphamin, ymax = alphamax)) +
  geom_line(size = 1, alpha = 0.5) +
  geom_ribbon(alpha = 0.5, size=0.5, linetype=0) + 
  theme_bw() +
  scale_color_manual(values = Colors) +
  xlab('Wavelength (nm)') + 
  ylab('Canopy reflectance (-)') +
  theme(panel.grid=element_blank())

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


# GSA
