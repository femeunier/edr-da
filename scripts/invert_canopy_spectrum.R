rm(list = ls())

# Library
library(redr)
library(PEcAn.all)
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

Colors <- c("#137300","#1E64C8")

par.wl = 400:2499 
nir.wl = 2500

edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt")

# Data
All_canopy_spectra <- readRDS(file= "~/data/RTM/All_canopy_spectra.rds")

refs <- unique(All_canopy_spectra$ref)

pecanxmlfile <- "~/output/PEcAn_99000000002/pecan.CONFIGS.xml"
settings <- PEcAn.settings::read.settings(pecanxmlfile)

modeloutdir <- file.path(settings$modeloutdir,"PDA_EDRTM")
settings$modeloutdir <- modeloutdir
settings$sensitivity.analysis$perpft <- TRUE
settings$sensitivity.analysis$ensemble.id <- NULL

Nensemble <- 500

alpha <- 0.05

# SA analysis
quantiles <- c(0.159,0.841)
quantiles_all <- c(0.159,0.5,0.841) # Add median
vars<- c("PAR_liana",'NIR','Spectrum','Vis','Red',"PAR")

# PFT table
df_PFT <- data.frame(PFTnum = as.numeric(unlist(lapply(lapply(settings$pfts, "[[","constants"),"[[","num"))),
                     names = as.character(unlist(lapply(settings$pfts, "[[","name"))))
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)
PFT <- as.character(df_PFT$names)
mapdf <- data.frame(pfts=c(as.character(df_PFT$names),"soil"),col=c(df_PFT$Col,"#000000"))
Values <- c("Liana_optical" = as.character(mapdf$col[1]),
            "Tree_optical" =as.character(mapdf$col[2]),
            "soil" = as.character(mapdf$col[3]))

use_meta.analysis <- FALSE
use_leaf_PDA <- FALSE
use_prior <- FALSE

alpha_frac <- 0.8
PFTselect <- 17
crown_mod <- 0

h5file <- "/home/carya/output/PEcAn_99000000002/out/SA-median/RTM/history-S-2004-01-01-120000-g01.h5"

npft <- length(df_PFT$names)

if (crown_mod == 0){
 dis2find_prim <- c('b1Bl_small','b2Bl_small','b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                     'Cab','carotenoids','Cw','Cm')
} else {
  dis2find_prim <- c('b1Bl_small','b2Bl_small','b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
                     'Cab','carotenoids','Cw','Cm','b1Ca','b2Ca')
}
dis.names <- c("ssigma","soil.moisture",rep(dis2find_prim,npft))

rundir <- "/home/carya/output/PEcAn_99000000002/run/SA-median"


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

start.date <- settings$run$start.date
date <- ISOdatetime(
  lubridate::year(start.date),
  lubridate::month(start.date),
  lubridate::mday(start.date),
  12, 00, 00, tz = "UTC")

edr_ed2in <- PEcAnRTM::setup_edr(ed2in, modeloutdir, date,TRUE)

ed2in <- PEcAn.ED2::read_ed2in(file.path(modeloutdir,"ED2IN"))
PEcAn.ED2::write_ed2in(ed2in, file.path(modeloutdir,"ED2IN"))

ED2IN_file <- file.path(modeloutdir,"ED2IN")
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

if(!dir.exists(modeloutdir)) dir.create(modeloutdir)

##########################################################################
# Loop over papers

best_sets <- c("~/data/RTM/current_samples_Kalacska.rds",
               "~/data/RTM/current_samples_marvin.rds",
               "~/data/RTM/current_samples_sanchez.rds")

# data.frames to save between loops
model_ensemble_all <- parameter_dis <- model_sensitivities_all <- p.values_all <-
  performance_all <- data.frame()

for (iref in seq(1,length(best_sets))){
  
  current_ref <- refs[iref]
  
  Spectrum_canopy_data <- All_canopy_spectra %>% select(c(ref,scenario,wavelength,Reflectance_median)) %>%
    rename(reflectance = Reflectance_median) %>% filter (ref == current_ref)
  
  observed <- observation <- Spectrum_canopy_data
  
  best_set <- readRDS(best_sets[iref])
  
  ensemble <- getSample(best_set,numSamples = Nensemble)
  ensemble <- rbind(MAP(best_set)$parametersMAP,ensemble)
  
  ll_all <- rep(NA,Nensemble)
  
  ensemble_results <- data.frame()
  print('... Ensemble runs ...')
  for (iensemble in seq(1,nrow(ensemble))){
    print(iensemble/nrow(ensemble))
    
    params <- ensemble[iensemble,]
    
    ssigma <- params[1]
    
    outputs <- run_ED_RTM(rundir,modeloutdir,params,crown_mod,inventory,par.wl,nir.wl)
    
    COI <- outputs[["COI"]]
    output_RTM <- outputs[["output_RTM"]]
    
    # classifify them
    patch_class <- list()
    patch_class[["low"]] <- which(COI < 0.25)
    patch_class[["high"]] <- which(COI > 0.5)
    
    if (any(sapply(patch_class,isempty))){
      ll_all[iensemble] <- -1e20
      next
    }
    
    # Get outputs
    scenarios <- as.character(unique(observed[["scenario"]]))
    ll <- rep(NA,length(scenarios))
    
    for (iscenar in seq(scenarios)){
      temp <- observed %>% filter(scenario == scenarios[[iscenar]]) %>%ungroup()
      
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
      
      ensemble_results <- rbind(ensemble_results,
                                data.frame(run = iensemble,
                                           wavelength = waves,
                                           scenar = scenarios[iscenar],
                                           sim = simulated_reflectance_Waves$y,
                                           obs = observed_Reflectance))
      
      ll[iscenar] <- sum(dnorm(simulated_reflectance_Waves$y, observed_Reflectance, ssigma, log = TRUE))
    }
    
    ll_all[iensemble] <- sum(ll)
  }
  
  ensemble_results_select <- ensemble_results %>% group_by(scenar,wavelength) %>% summarise(rmin = min(sim,na.rm = TRUE),
                                                                                            rmax = max(sim,na.rm = TRUE),
                                                                                            alphamin = quantile(sim,alpha/2,na.rm = TRUE),
                                                                                            alphamax = quantile(sim,1-alpha/2,na.rm = TRUE),
                                                                                            median = median(sim,na.rm = TRUE))
  best_run <- which.max(ll_all)
  ensemble_results_select[["rbest"]] <- ensemble_results %>% filter(run == best_run) %>% pull(sim)
  
  model_ensemble_all <- rbind(model_ensemble_all,
                              ensemble_results_select %>% ungroup() %>% mutate(ref = current_ref))
  
  # ggplot(ensemble_results_select,
  #        aes(x = wavelength,y = median,colour = as.factor(scenar))) +
  #   geom_line(size = 1, alpha = 0.5,linetype = 2) +
  #   geom_ribbon(aes(ymin = alphamin, ymax = alphamax,fill=as.factor(scenar)),alpha = 0.5, size=0.5, linetype=0) + 
  #   theme_bw() +
  #   scale_color_manual(values = df_PFT$Col) +
  #   scale_fill_manual(values = df_PFT$Col) +
  #   xlab('Wavelength (nm)') + 
  #   ylab('Canopy reflectance (-)') +
  #   theme(panel.grid=element_blank())+
  #   geom_point(data = observation, aes(x = wavelength, y = reflectance,colour = as.factor(scenario)),size=1)

  data_comp <- ensemble_results %>% group_by(scenar,wavelength) %>% summarise(sim = median(sim),
                                                                              obs = median(obs))
  
  data_comp[["simbest"]] <- ensemble_results %>% filter(run == best_run) %>% pull(sim)
  
  summary(lm(data = data_comp,
             formula = obs ~ sim))$adj.r.squared
  
  performance_all <- rbind(performance_all,
                           data.frame(ref = current_ref,
                                      wave = data_comp$wavelength,
                                      scenar = data_comp$scenar,
                                      sim = data_comp$sim,
                                      obs = data_comp$obs,
                                      sim_best = data_comp$simbest))

  
  #########################################################
  # Univariate sensitivity analysis
  print('... SA runs ...')
  
  sampling <- getSample(best_set,parametersOnly = TRUE,numSamples = 10000)
  median_bs <- apply(sampling,2,quantile,0.5)
  
  # Median simulations
  OP_dir <- file.path("/home/carya/output/PEcAn_99000000002/out/PDA_EDRTM","SA-median")
  if (!dir.exists(OP_dir)){
    dir.create(OP_dir)
    dir.create(file.path(OP_dir,"patches"))
  }
  OP <- run_ED_RTM(rundir,modeloutdir,params=median_bs,crown_mod,inventory,par.wl,nir.wl,
                   temp = FALSE,outputdir = OP_dir)
  model2netcdf.EDR(outdir = OP_dir,
                   sitelat = settings$run$site$lat,
                   sitelon = settings$run$site$lon,
                   start_date = start.date,
                   par.wl = par.wl,
                   nir.wl = nir.wl,
                   patches = TRUE,
                   PFTselect = PFTselect)
  
  Nparams <- length(median_bs)
  Names_params <- dis.names
  
  settings$modeloutdir <- modeloutdir
  
  # SA simulations
  PFT_names <- as.character(df_PFT$names)
  PFT_all <- c("","",rep(PFT_names[1],(length(median_bs)-2)/2),rep(PFT_names[2],(length(median_bs)-2)/2))
  
  for (iquantile in seq(quantiles)){
    quantiles_bs <- apply(sampling,2,quantile,quantiles[iquantile])
    
    for (iparam in seq(1,Nparams)){
      temp_params <- median_bs
      temp_params[iparam] <- quantiles_bs[iparam]
      
      simname <- paste0("SA-",ifelse(!identical(PFT_all[iparam],""),paste0(PFT_all[iparam],"-"),""),Names_params[iparam],"-",quantiles[iquantile])
      OP_dir <- file.path(modeloutdir,simname)
      if (!dir.exists(OP_dir)){
        dir.create(OP_dir)
        dir.create(file.path(OP_dir,"patches"))
      }
      
      OP <- run_ED_RTM(rundir,modeloutdir,params=temp_params,crown_mod,inventory,par.wl,nir.wl,
                       temp = FALSE,outputdir = OP_dir)

      model2netcdf.EDR(OP_dir,sitelat = settings$run$site$lat,sitelon = settings$run$site$lon,
                       start_date = start.date,par.wl = par.wl,nir.wl = nir.wl,patches = TRUE)
    }
  }
  
  # Variables to save
  sa.run.ids <- list()
  sa.samples <- list(env = data.frame())
  pft.names <- PFT_names
  trait.names <- trait.samples <- list()
  
  delta_param = 0
  for (pft in (PFT_names)){
    mat <- runs <- array(NA,c(length(quantiles_all),(Nparams-2)/2))
    rownames(mat) <-  rownames(runs) <- quantiles_all*100
    colnames(mat) <- colnames(runs) <- Names_params[3:(((Nparams-2)/2)+2)]
    
    for (iparam in seq(3,(((Nparams-2)/2)+2))){
      for (iquantile in seq(quantiles_all)){
        simname <- paste0("SA-",pft,"-",Names_params[iparam],"-",quantiles_all[iquantile])
        quantiles_bs <- apply(sampling,2,quantile,quantiles_all[iquantile])
        
        temp_params <- median_bs
        temp_params[delta_param+iparam] <- quantiles_bs[delta_param+iparam]
        
        mat[iquantile,iparam-2] <- temp_params[delta_param+iparam]
        if (quantiles_all[iquantile] != 0.5){
          runs[iquantile,iparam-2] <- simname
        } else {
          runs[iquantile,iparam-2] <- "SA-median"       
        }
      }
      trait.samples[[pft]][[Names_params[iparam]]] <- sampling[,delta_param+iparam]
    }
    
    sa.run.ids[[pft]] <- runs
    sa.samples[[pft]] <- as.data.frame(mat)
    trait.names[[pft]] <- Names_params[seq(3,((Nparams-2)/2)+2)]
    
    delta_param <- delta_param + ((Nparams-2)/2)
  }
  
  temp <- array(as.character(c("SA-soil.moisture-0.159","SA-median","SA-soil.moisture-0.841")),
        c(3,1))
  rownames(temp) <- quantiles_all*100
  colnames(temp) <- "soil_moisture"
  sa.run.ids[["soil"]] <- temp
  
  pft.names <- c(pft.names,"soil")
  
  trait.names[["soil"]] <- "soil_moisture"
  trait.samples[["soil"]][["soil_moisture"]] <- sampling[,2]
  
  temp <- as.data.frame(quantile(sampling[,2],quantiles_all))
  rownames(temp) <- quantiles_all*100
  colnames(temp) <- "soil_moisture"
  sa.samples[["soil"]] <- temp
  
  settings$pfts[[3]] <- list(name = "soil",outdir = file.path(settings$outdir,"pft","soil"))
  
  save(sa.run.ids,pft.names,trait.names,trait.samples,sa.samples, file = file.path(settings$outdir,"sensitivity.samples.NOENSEMBLEID.Rdata"))
  save(trait.names,pft.names,trait.samples, file = file.path(settings$outdir,"samples.Rdata"))
  
  # Acutal senstivitiy analysis
  for (ivar in seq(1,length(vars))){
    OP_all <- VDP_allPFTs(vars[ivar],mapdf = mapdf)
    model_sensitivities_all <-
      rbind(model_sensitivities_all,
          data.frame(ref = current_ref,OP_variable = vars[ivar],
                     pft = OP_all$PFT_all,CV=OP_all$coef.vars,
                     sens = OP_all$sensitivities,variance = OP_all$variances,
                     par.var = OP_all$partial.variances,
                     variable = names(OP_all$elasticities)))
  }
  
  if (!dir.exists(file.path(settings$outdir,"pfts",current_ref))) dir.create(file.path(settings$outdir,"pfts",current_ref))
  
  # Save pdf and OP files
  system2("cp",args=list("-r",file.path(settings$outdir,"pfts","all","*pdf"),file.path(settings$outdir,"pfts",current_ref)))
  system2("cp",args=list(file.path(settings$outdir,"*Rdata"),file.path(settings$outdir,"pfts",current_ref)))
  
  # load("/home/carya/output/PEcAn_99000000002/sensitivity.output.NOENSEMBLEID.NIR.2004.2004.Rdata")
 
  colnames(ensemble) <- paste(dis.names,PFT_all,sep='^')
  # Merge dataframes
  
  param_dis_current <-
    melt(ensemble) %>% select(c(Var2,value)) %>% rename(name = Var2) %>% 
      mutate(param = sub("\\^.*", "", name),
             pft = sub(".*\\^","",name),
             ref = current_ref) %>% mutate(pft = case_when(
               pft == "" ~ "soil",
               TRUE ~ pft)) %>% select(-c(name))
  
  parameter_dis <- rbind(parameter_dis,param_dis_current)
  
  # p-values
  ensemble_small <- getSample(best_set,numSamples = 30)
  colnames(ensemble_small) <- paste(dis.names,PFT_all,sep='^')
  
  param_dis_current_small <-
    melt(ensemble_small) %>% select(c(Var2,value)) %>% rename(name = Var2) %>% 
    mutate(param = sub("\\^.*", "", name),
           pft = sub(".*\\^","",name),
           ref = current_ref) %>% mutate(pft = case_when(
             pft == "" ~ "soil",
             TRUE ~ pft)) %>% select(-c(name)) %>% filter(pft != 'soil')
  
  params2compare <- unique(param_dis_current_small$param) 
  Nparams2compare <- length(params2compare)
  p_values <- rep(NA,Nparams2compare)
  names(p_values) <- params2compare
  
  for (trait in params2compare){
    
    data_aov <- param_dis_current_small %>% filter(param == trait)
    
    Test <- kruskal.test(formula = value ~ as.factor(pft), data = data_aov)
    p_values[trait] <- Test$p.value
  }
  
  
  p.values_all <- rbind(p.values_all,
                        data.frame(param = names(p_values),
                                   p.values = p_values,
                                   ref = current_ref))
}

save.image(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")

######################################################################
# Figures

# Model prediction vs inputs
OP_plot <- ggplot(data = model_ensemble_all,
       aes(x = wavelength,
           y = rbest,
           ymin = alphamin,
           ymax = alphamax,
           fill = scenar,
           color = scenar)) +
  geom_ribbon(alpha = 0.5,size = 0.5,linetype=0) +
  geom_line() +
  geom_point(data = All_canopy_spectra,
            aes(x = wavelength,
                y = Reflectance_median,
                ymin = Reflectance_min,
                ymax = Reflectance_max,
                fill = scenario,
                color = scenario),linetype = 2) + 
  facet_wrap(ref ~ .,scales = "free") +
  scale_color_manual(values = c("#137300","#1E64C8")) +
  scale_fill_manual(values = c("#137300","#1E64C8")) +
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","PDA_EDRTM_outputs.png"),dpi = 300,
       plot = OP_plot,
       width = 20, height = 10,units = "cm")

  
# Model performance all together

performance_plot <- ggplot(performance_all,aes(x=obs,y=sim_best,colour = as.factor(scenar),fill = as.factor(scenar),
                           shape = as.factor(ref))) +
  scale_color_manual(values = c("#137300","#1E64C8")) +
  scale_fill_manual(values = c("#137300","#1E64C8")) +
  scale_shape_manual(values=seq(0,length(unique(performance_all$ref)))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,colour = 'black',linetype=2) + 
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","PDA_EDRTM_performance.png"),dpi = 300,
       plot = performance_plot,
       width = 20, height = 15,units = "cm")

# All posterior sensitivities

model_sensitivities_all <- model_sensitivities_all %>% mutate(variable.pft = paste(variable,pft,sep = "."))

ready.to.plot <- model_sensitivities_all %>% filter(100*par.var >= 0.)
ready.to.plot[5,"par.var"] <- 0.5
ready.to.plot <- ready.to.plot %>%  group_by(OP_variable,ref) %>%
  top_n(5, abs(par.var)) %>% ungroup() %>% arrange(OP_variable,ref,(par.var)) %>%
  mutate(order = row_number())

sensitivity_plot<-ggplot(data = ready.to.plot %>% filter(OP_variable %in% c("PAR_liana","PAR")),
       aes(x = order,y = 100*par.var,color = as.factor(pft), fill = as.factor(pft))) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7) +
  facet_wrap(ref~OP_variable,scales = "free",nrow = length(unique(ready.to.plot$ref))) +
  theme_bw() + 
  scale_fill_manual(values = Values) +
  scale_color_manual(values = Values) +
  scale_x_continuous(
    breaks = ready.to.plot$order,
    labels = ready.to.plot$variable,
    expand = c(0.01,0.01)
  ) +
  coord_flip() +
  ylab("Partial variance (%)") + 
  xlab("") 

ggsave(filename = file.path("~/data/RTM/","PDA_EDRTM_sensitivity.png"),dpi = 300,
       plot = sensitivity_plot,
       width = 15, height = 12,units = "cm")


# Parameter posterior distribution
temp_max <- parameter_dis %>% filter(pft!= "soil")%>% group_by(param) %>%
  summarise(value_max = max(value)) %>% ungroup()
p.values_all <- p.values_all %>% left_join(temp_max,by = "param") %>% mutate(signif = paste(case_when(
  p.values < 0.001 ~ '***',
  p.values < 0.01 ~ '**',
  p.values < 0.05 ~ '*',
  TRUE ~ 'N.S.')))

marginal_plot <- ggplot(data = parameter_dis %>% filter(pft!= "soil"),
       aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  geom_text(data = p.values_all,
            aes(x = 1.1*value_max, y = ref, label = signif),color = 'black',
            nudge_y = 0.5,
            size = 2) + 
  labs(y="") +
  theme_ridges(font_size = 13) +
  facet_wrap(param ~ .,scales = "free_x") +
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","PDA_EDRTM_marginalD.png"),dpi = 300,
       plot = marginal_plot,
       width = 30, height = 20,units = "cm")
