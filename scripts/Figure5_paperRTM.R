rm(list = ls())

library(dplyr)
library(hdf5r)
library(redr)
library(pracma)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(BayesianTools)
library(zoo)

load(file = "~/data/RTM/Inverse_canopy_spectrum_N250.Rdata")

Select = "Marvin"  # "Kalacska" "Marvin"   "Sanchez" 

Nensemble = 250
alpha = 0.05
crown_mod = 0
par.wl = c(400:2499)
nir.wl = c(2500)
PFTselect = 1
alpha_frac = 0.8

#############################################################################################################
# Prepare BCI runs

h5file <- "~/data/BCI/bci-S-2004-01-01-000000-g01.h5"

Colors <- c("#137300","#1E64C8")
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)  

Lindex <- calc_liana_index(h5file,PFTselect,alpha_frac)

par(mfrow=c(1,4)) 

boxplot(Lindex$COI)
boxplot(Lindex$COI2)
plot(Lindex$COI,Lindex$COI2)
plot(Lindex$fraction_LAI,Lindex$COI)

rundir <- "/home/carya/output/PEcAn_99000000002/run/SA-median"
modeloutdir <- "/home/carya/output/PEcAn_99000000002/out/PDA_EDRTM"
edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt")

# Modify ED2IN
ed2in$SFILIN <- file.path(dirname(h5file), sub("\\-.*", "",basename(h5file)))
ed2in$RUNTYPE = "HISTORY"
ed2in$POI_LAT = 9.2
ed2in$POI_LON = -79.8

PEcAn.ED2::write_ed2in(ed2in, file.path(modeloutdir,"ED2IN"))

ED2IN_file <- file.path(modeloutdir,"ED2IN")

hfile <- hdf5r::H5File$new(h5file)

dbh <- readDataSet(hfile[["DBH"]])
nplant <- readDataSet(hfile[["NPLANT"]])
Npatch <- readDataSet(hfile[["NPATCHES_GLOBAL"]])
hite <- readDataSet(hfile[["HITE"]])
pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
PACO_N <- readDataSet(hfile[["PACO_N"]])
hfile$close_all()

Ncohort <- length(dbh)
npft <- max(as.numeric(pft))

Npatch <- length(PACO_N)
PAN <- rep(1:Npatch,PACO_N)

inventory <- list(dbh = dbh,pft = pft,hite = hite,nplant = nplant, PAN = PAN,PFTselect = PFTselect,
                  PACO_N = PACO_N, Ncohort = Ncohort, Npatch = Npatch)

OP_dir <- modeloutdir
h5file_paste <- file.path(modeloutdir,"history-S-2004-01-01-000000-g01.h5")
system2("cp",args = list(h5file,h5file_paste))

#############################################################################################################
current_ref <- Select

Spectrum_canopy_data <- All_canopy_spectra %>% dplyr::select(c(ref,scenario,wavelength,Reflectance_median)) %>%
  rename(reflectance = Reflectance_median) %>% filter (ref == current_ref)

observed <- observation <- Spectrum_canopy_data

set_select <- which(grepl(Select,best_sets))
best_set <- readRDS(best_sets[set_select])

ensemble <- getSample(best_set,numSamples = Nensemble,start = 3500)
ensemble <- rbind(MAP(best_set)$parametersMAP,
                  ensemble)

ensemble_results <- data.frame()
ll_all <- rep(NA,Nensemble)

print('... Ensemble runs ...')
for (iensemble in seq(1,nrow(ensemble))){
  print(iensemble/nrow(ensemble))

  params <- ensemble[iensemble,]

  ssigma <- params[1]

  outputs <- run_ED_RTM(rundir,
                        outdir = modeloutdir,params,crown_mod,inventory,par.wl,nir.wl)

  # COI <- outputs[["COI"]]
  COI <- Lindex$fraction_LAI
  output_RTM <- outputs[["output_RTM"]]

  # classifify them
  patch_class <- list()
  patch_class[["low"]] <- which(COI < 0.05)
  patch_class[["high"]] <- which(COI > 0.4)

  if (any(sapply(patch_class,isempty)) | all(is.na(output_RTM))){
    next
  }

  # Get outputs
  scenarios <- as.character(unique(observed[["scenario"]]))

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
                                         wavelength = c(par.wl,nir.wl),
                                         scenar = scenarios[iscenar],
                                         sim = simulated_reflectance))

    ll[iscenar] <- sum(dnorm(simulated_reflectance_Waves$y, observed_Reflectance, ssigma, log = TRUE))
  }
  ll_all[iensemble] <- sum(ll)
}

# load("~/data/RTM/BCI.RData")



# sA <-
#   ggplot(p_values) +
#   geom_abline(slope=0, intercept = 0.05) +
#   geom_point(aes(x = wavelength,
#                  y = p_value))

####################################################################################

rm(list = ls())

library(PEcAnRTM)
library(dplyr)
library(hdf5r)
library(redr)
library(pracma)
library(ggplot2)
library(Hmisc)
library(cowplot)
library(BayesianTools)
library(zoo)
library(reshape2)
library(RColorBrewer)

load("~/data/RTM/BCI.RData")

best_run <- which.max(ll_all)
runs <- seq(1,nrow(ensemble))

print(length(runs[ll_all > (0.9)*ll_all[best_run] & !is.na(ll_all)]))
ensemble_results_filtered <- ensemble_results %>% filter(run %in% runs[ll_all > (0.7)*ll_all[best_run] & !is.na(ll_all)])

ensemble_results_select <- cbind(ensemble_results_filtered %>% group_by(scenar,wavelength) %>% summarise(rmin = min(sim,na.rm = TRUE),
                                                                                                         rmax = max(sim,na.rm = TRUE),
                                                                                                         alphamin = quantile(sim,alpha,na.rm = TRUE),
                                                                                                         alphamax = quantile(sim,1-alpha,na.rm = TRUE),
                                                                                                         median = median(sim,na.rm = TRUE)),
                                 rbest = ensemble_results_filtered %>% group_by(scenar,wavelength) %>% filter(run == best_run) %>% pull(sim))


# p_values <- 
#   ensemble_results_filtered %>% group_by(wavelength) %>% 
#   summarise(p_value = kruskal.test(formula = sim ~ as.factor(scenar))$p.value)


p_values <-
  ensemble_results_filtered %>% group_by(wavelength) %>%
  summarise(p_value = kruskal.test(formula = sim ~ as.factor(scenar))$p.value)


modelled.spectra <- ensemble_results_select %>% ungroup() %>% mutate(scenar = case_when(
  scenar == "low" ~'Low',
  scenar == "high" ~'High'))

ensemble_results_diff <- cbind(ensemble_results_filtered %>% filter(scenar == "low"),
                               sim_h = ensemble_results_filtered %>% filter(scenar == "high") %>% mutate(sim_h = sim) %>% pull(sim_h)) %>%
  mutate(diff = sim - sim_h)
ensemble_results_diff_select <- ensemble_results_diff %>% group_by(wavelength) %>% summarise(diff_min = min(diff,na.rm = TRUE),
                                                                                             diff_max = max(diff,na.rm = TRUE),
                                                                                             diff_alphamin = wilcox.test(diff,conf.int=TRUE)$conf.int[1],
                                                                                             diff_alphamax = wilcox.test(diff,conf.int=TRUE)$conf.int[2],
                                                                                             diff_median = median(diff,na.rm = TRUE))

p_values_diff <- 
  ensemble_results_diff %>% group_by(wavelength) %>% 
  summarise(p.val = wilcox.test(diff)$p.value)




sensor.list <- c("identity","aviris.ng", "aviris.classic",
                 "hyperion", "chris.proba", "landsat5", "landsat7",
                 "landsat8", "modis", "viirs", "avhrr")

sensor.proper <- c("Full spectrum","AVIRIS NG", "AVIRIS Classic",
                   "Hyperion", "CHRIS-Proba", "Landsat 5", "Landsat 7",
                   "Landsat 8", "MODIS", "VIIRS", "AVHRR")

good.runs <- unique(ensemble_results_filtered %>% pull(run))
scenars <- unique(ensemble_results_filtered %>% pull(scenar))

sensor.rsr <- PEcAnRTM:::generate.rsr.all()
sensor.rsr[["identity"]] <- cbind(400:2500,diag(x=1,nrow=length(400:2500)))
colnames(sensor.rsr[["identity"]]) <- c('index',1:2101)
rownames(sensor.rsr[["identity"]]) <- 400:2500

data.sensors <- data.frame()
sensors.bands <- data.frame()
tol <- 1e-5
for (sensor in sensor.list){
  rsr <- sensor.rsr[[sensor]][,-1]
  data.sensors <- rbind(data.sensors,
                        data.frame(melt(rsr) %>% rename(wl = Var1,
                                                        band = Var2,
                                                        cv = value) %>% mutate(sensor = sensor)))
  
  bands <- ncol(rsr)
  bands.names <- colnames(rsr)
  for (iband in seq(1,bands)){
    pos <- abs(rsr[,iband])>tol
    if( sensor == "identity"){
      all_wl <- as.numeric((sensor.rsr[[sensor]][pos,1]))
    } else {
      all_wl <- as.numeric(names(sensor.rsr[[sensor]][pos,1]))
    }

    sensors.bands <- rbind(sensors.bands,
                           data.frame(sensor = sensor,
                                      band = as.character(bands.names[iband]),
                                      wl =  c(min(all_wl),max(min(all_wl)+1,max(all_wl)))))
                                      
  }
}

remote.sensing <- data.frame()

Nmax <- 30
for (i in good.runs[1:Nmax]){
  cat(i)
  for (scenario in scenars){
    
    spec <- ensemble_results_filtered %>% filter(run == i & scenar == scenario)
    spectrum <- spectra(spectra = spec %>% pull(sim), 
                        wavelengths = spec %>% pull(wavelength)) 
    
    for (sensor in sensor.list){
      SR <- spectral.response(spec=spectrum, sensor=sensor)
      rsr <- sensor.rsr[[sensor]][,-1]
      if (sensor == "identity"){
        remote.sensing <- rbind(remote.sensing,
                                data.frame(run = i,
                                           scenar = scenario,
                                           sensor = sensor,
                                           band = as.numeric(colnames(rsr)),
                                           data = SR))
      } else{
        remote.sensing <- rbind(remote.sensing,
                                data.frame(run = i,
                                           scenar = scenario,
                                           sensor = sensor,
                                           band = colnames(rsr),
                                           data = SR))
      }

    }
  }
}



# p.value_RS <- remote.sensing %>% group_by(sensor,band) %>% summarise(p_value = kruskal.test(formula = data ~ as.factor(scenar))$p.value) %>%
#   mutate(signif = paste(case_when(
#     p_value < 0.001 ~ '***',
#     p_value < 0.01 ~ '**',
#     p_value < 0.05 ~ '*',
#     TRUE ~ 'N.S.'))) %>% ungroup %>% mutate(sensor = as.character(sensor))


p.value_RS <- cbind(remote.sensing %>% group_by(sensor,band) %>% filter(scenar == "low"),
                    data_h=remote.sensing %>% group_by(sensor,band) %>% filter(scenar == "high") %>% mutate(data_h = data) %>% pull(data)) %>%
  mutate(diff = data_h - data) %>% ungroup() %>% group_by(sensor,band) %>% summarise(p_value = wilcox.test(diff)$p.value) %>%
  mutate(signif = paste(case_when(
    p_value < 0.001 ~ '***',
    p_value < 0.01 ~ '**',
    p_value < 0.05 ~ '*',
    TRUE ~ 'N.S.'))) %>% ungroup %>% mutate(sensor = as.character(sensor)) %>% ungroup()


remote.sensing_diff <-
  cbind(remote.sensing %>% filter(sensor == "identity") %>% filter(scenar == "low") %>% rename(low = data),
      high = remote.sensing %>% filter(sensor == "identity") %>% filter(scenar == "high") %>% rename(high = data) %>% pull(high)) %>%
  mutate(diff = low - high,
         wavelength = as.numeric(band) +399)


remote.sensing_diff_CI <- remote.sensing_diff %>% 
  group_by(wavelength) %>% summarise(diff_median = median(diff,na.rm=TRUE),
                                     diff_alphamin = wilcox.test(diff,conf.int=TRUE)$conf.int[1],
                                     diff_alphamax = wilcox.test(diff,conf.int=TRUE)$conf.int[2])

# p.value_RS %>% filter(sensor == "landsat5")
# p.value_RS %>% filter(sensor == "landsat7")
# p.value_RS %>% filter(sensor == "identity")
# 
# band_s <- sort(as.numeric(p.value_RS[p.value_RS$sensor == "identity",] %>% pull(band)),index.return=TRUE)
# p.value_RS[p.value_RS$sensor == "identity",] <- 
# 
# p.value_RS <- p.value_RS %>% group_by(sensor) %>% arrange(band) %>%ungroup()

subplot.C_bis <- subplot.C <- data.frame()

compt=length(sensor.list)
for (sens in sensor.list){
  rsr <- sensor.rsr[[sens]][,-1]

  bands <- ncol(rsr)
  bands.names <- colnames(rsr)
  
  for (iband in seq(1,bands)){
    p.value_RS_sel <- p.value_RS %>% filter(band == bands.names[iband] & sensor == sens)
    sensors.bands_sel <- sensors.bands %>% filter(band == bands.names[iband] & sensor == sens)
    if (p.value_RS_sel %>% pull( p_value) < 0.05){
      subplot.C <- rbind(subplot.C,
                              data.frame(rbind(sensors.bands_sel,sensors.bands_sel[2,]),
                                         y_min = c(compt,compt,NA),y_max=c(compt,compt,NA)+0.25))
    } 
    subplot.C_bis <- rbind(subplot.C_bis,
                       data.frame(rbind(sensors.bands_sel,sensors.bands_sel[2,]),
                                  y_min = c(compt,compt,NA),y_max=c(compt,compt,NA)+0.25))
  }
  compt = compt-1
}

sC <- ggplot() + 
  geom_ribbon(data = subplot.C_bis,mapping = aes(x = wl,ymin = y_min, ymax = y_max,group =as.factor(sensor)),show.legend = FALSE,
              alpha = 0.3,linetype=0) +
  geom_ribbon(data = subplot.C,mapping = aes(x = wl,ymin = y_min, ymax = y_max,group =as.factor(sensor)),show.legend = FALSE,
              alpha = 1,linetype=0) +
  labs(y = "",
       x = "Wavelength [nm]") +
  scale_y_continuous(breaks = seq(length(sensor.proper),compt+1),labels = sensor.proper) +
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

final.subplotB <-  remote.sensing_diff_CI %>% mutate(scenar = 'test')

sB <-
  ggplot(data =final.subplotB,
         mapping = aes(x = wavelength,
                       ymin = diff_alphamin,
                       ymax = diff_alphamax,
                       fill = as.factor(scenar),
                       colour = as.factor(scenar),
                       y = diff_median)) +
  geom_ribbon(alpha = 0.5,size = 0.5,linetype=0,show.legend = FALSE) +
  geom_line(linetype = 1) +
  scale_y_continuous(limits = c(-0.005,0.02)) +
  scale_x_continuous() +
  theme_bw() + 
  labs(x = "Wavelength [nm]",
       y = "Reflectance difference [-]",
       colour = "") +
  scale_colour_manual(labels = "Low - High",values = "black") +
  scale_fill_manual(values = "darkgrey") +
  geom_abline(slope=0,intercept=0,linetype=3,colour='red',size=0.4) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.position = c(0.9, 0.75))

sA <-
  ggplot() +
  # geom_ribbon(data = modelled.spectra,
  #             aes(x = wavelength,
  #                 ymin = rbest - (median-alphamin),
  #                 ymax = rbest +  (alphamax-median),
  #                 fill = scenar,
  #                 color = scenar),alpha = 0.5,size = 0.5,linetype=0) +
  geom_line(data = modelled.spectra,
            aes(x = wavelength,
                y = rbest,
                color = scenar),linetype = 1) +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  scale_y_continuous(limits = c(0,0.45)) +
  scale_x_continuous() +
  theme_bw() + 
  labs(x = "Wavelength [nm]",
       y = "Canopy reflectance [-]",
       fill = "Liana infestation",
       colour = "Liana infestation") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        legend.position = c(0.9, 0.8))


plot_grid(sA,sB,sC,
          align = c("hv"),
          nrow = 3)

ggsave(plot = last_plot(),
       dpi = 300,
       width = 30,
       height = 25,
       units = "cm",
       file = "~/data/RTM/Figure5.png")

