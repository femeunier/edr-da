rm(list = ls())

######################################################################################################
# Leaf reflectance
######################################################################################################
library(redr)

load(file = "~/data/RTM/Inverse_leaf_spectrum2.Rdata")

posterior <- ensemble_posterior_all %>% dplyr::select(wavelength,alphamin,alphamax,median,ref,pft) %>% 
  rename(posterior_alphamin=alphamin,
         posterior_median=median,
         posterior_alphamax=alphamax)

ready.for.table <- All_leaf_spectra %>% mutate(wavelength = as.integer(round(wavelength))) %>%left_join(posterior,by=c("ref","pft","wavelength"))

####################################################################################################
references <- unique(All_leaf_spectra$ref)

all_data <- all_models <- data.frame()
for (reference in references){
  
  
  ##################################################################
  # Data 
  temp_T <- All_leaf_spectra %>% filter(ref == reference &
                                   pft == "Tree_optical")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(Reflectance_median)
  
  
  temp_L <- All_leaf_spectra %>% filter(ref == reference &
                                   pft == "Liana_optical")
  wavelengths_L <- temp_L %>% pull(wavelength)
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L <- temp_L %>% pull(Reflectance_median)
  
  R_L_out <- approxExtrap(x = wavelengths_L,
                     y = R_L,
                     xout = wavelengths_L_out)$y
  
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
               y = R_T,
               xout = wavelengths_L_out)$y
  
  all_data <- rbind(all_data,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out))
  
  ##################################################################
  # Models 
  temp_T <- posterior %>% filter(ref == reference &
                                          pft == "Tree_optical")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(posterior_median)
  
  
  temp_L <- posterior %>% filter(ref == reference &
                                          pft == "Liana_optical")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(posterior_median)
  
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L_out <- approxExtrap(x = wavelengths_L,
                          y = R_L,
                          xout = wavelengths_L_out)$y
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
                 y = R_T,
                 xout = wavelengths_L_out)$y
  
  all_models <- rbind(all_models,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out))
}

# Data
Delta_PAR <- all_data %>% filter(wv <=800)
metrics_delta_PAR <- redr::calc_metrics(obs=Delta_PAR %>% pull(R_L),
                                  sim=Delta_PAR %>% pull(R_T))

Delta_NIR <- all_data %>% filter(wv >800)
metrics_delta_NIR <- redr::calc_metrics(Delta_NIR %>% pull(R_L),
                                  Delta_NIR %>% pull(R_T))

# Models
Delta_PAR_m <- all_models %>% filter(wv <=800)
metrics_delta_PAR_m <- redr::calc_metrics(Delta_PAR_m %>% pull(R_L),
                                    Delta_PAR_m %>% pull(R_T))

Delta_NIR_m <- all_models %>% filter(wv >800)
metrics_delta_NIR_m <- redr::calc_metrics(Delta_NIR_m %>% pull(R_L),
                                    Delta_NIR_m %>% pull(R_T))

#####################################################################################################
# Calculation
PAR_tree <- ready.for.table %>% filter(pft == "Tree_optical" & wavelength<=800)
metrics_PAR_tree <- redr::calc_metrics(obs=PAR_tree%>%pull(Reflectance_median),
                                 sim=PAR_tree%>%pull(posterior_median))

PAR_liana <- ready.for.table %>% filter(pft == "Liana_optical" & wavelength<=800)
metrics_PAR_liana <- redr::calc_metrics(obs=PAR_liana%>%pull(Reflectance_median),
                                 sim=PAR_liana%>%pull(posterior_median))

NIR_tree <- ready.for.table %>% filter(pft == "Tree_optical" & wavelength>800)
metrics_NIR_tree <- redr::calc_metrics(obs=NIR_tree%>%pull(Reflectance_median),
                                 sim=NIR_tree%>%pull(posterior_median))

NIR_liana <- ready.for.table %>% filter(pft == "Liana_optical" & wavelength>800)
metrics_NIR_liana <- redr::calc_metrics(obs=NIR_liana%>%pull(Reflectance_median),
                                 sim=NIR_liana%>%pull(posterior_median))


############################################################################################################
# Canopy reflectance
############################################################################################################

load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")

model.results <- model_ensemble_all %>% dplyr::select(scenar,wavelength,rbest,alphamin,alphamax,ref) %>% rename(scenario = scenar,median = rbest)
data <- All_canopy_spectra %>% dplyr::select(ref,wavelength,scenario,Reflectance_median,Reflectance_alphamin,Reflectance_alphamax)

data_canopy <- data %>% left_join(model.results,by=c("ref","scenario","wavelength")) 

# differences 
references <- unique(data_canopy$ref)

all_data <- all_models <- data.frame()
for (reference in references){
  
  #######################################################
  # Data
  temp_T <- data_canopy %>% filter(ref == reference &
                                     scenario == "low")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(Reflectance_median)
  
  
  temp_L <- data_canopy %>% filter(ref == reference &
                                     scenario == "high")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(Reflectance_median)
  
  wavelengths_T_out <- seq(round(min(wavelengths_T))+1,round(max(wavelengths_T))-1,1)
  R_T_out <- approxExtrap(x = wavelengths_T,
                          y = R_T,
                          xout = wavelengths_T_out)$y
  
  
  R_L_extrap <- 
    approxExtrap(x = wavelengths_L,
                 y = R_L,
                 xout = wavelengths_T_out)$y
  
  all_data <- rbind(all_data,
                    data.frame(wv = wavelengths_T_out,
                               R_T = R_T_out,
                               R_L = R_L_extrap))
  
  #######################################################
  # Models
  temp_T <- model.results %>% filter(ref == reference &
                                     scenario == "low")
  wavelengths_T <- temp_T %>% pull(wavelength)
  R_T <- temp_T %>% pull(median)
  
  
  temp_L <- model.results %>% filter(ref == reference &
                                     scenario == "high")
  wavelengths_L <- temp_L %>% pull(wavelength)
  R_L <- temp_L %>% pull(median)
  
  wavelengths_L_out <- seq(round(min(wavelengths_L))+1,round(max(wavelengths_L))-1,1)
  R_L_out <- approxExtrap(x = wavelengths_L,
                          y = R_L,
                          xout = wavelengths_L_out)$y
  
  
  R_T_extrap <- 
    approxExtrap(x = wavelengths_T,
                 y = R_T,
                 xout = wavelengths_L_out)$y
  
  all_models <- rbind(all_models,
                    data.frame(wv = wavelengths_L_out,
                               R_T = R_T_extrap,
                               R_L = R_L_out))
}

# Data
Delta_PAR_canopy <- all_data %>% filter(wv <=700)
metrics_delta_PAR_canopy <- redr::calc_metrics(Delta_PAR_canopy %>% pull(R_L),
                                         Delta_PAR_canopy %>% pull(R_T))

Delta_NIR_canopy <- all_data %>% filter(wv >800) 
metrics_delta_NIR_canopy <- redr::calc_metrics(Delta_NIR_canopy %>% pull(R_L),
                                         Delta_NIR_canopy %>% pull(R_T))

# Models
Delta_PAR_canopy_m <- all_models %>% filter(wv <=700)
metrics_delta_PAR_canopy_m <- redr::calc_metrics(Delta_PAR_canopy_m %>% pull(R_L),
                                           Delta_PAR_canopy_m %>% pull(R_T))

Delta_NIR_canopy_m <- all_models %>% filter(wv >800) 
metrics_delta_NIR_canopy_m <- redr::calc_metrics(Delta_NIR_canopy_m %>% pull(R_L),
                                           Delta_NIR_canopy_m %>% pull(R_T))

#############################################################################################################
# Calculations

PAR_tree_canopy <- data_canopy %>% filter(scenario == "low" & wavelength<=800)
metrics_PAR_tree_canopy <- redr::calc_metrics(obs=PAR_tree_canopy%>%pull(Reflectance_median),
                                 sim=PAR_tree_canopy%>%pull(median))

PAR_liana_canopy <- data_canopy %>% filter(scenario == "high" & wavelength<=800)
metrics_PAR_liana_canopy <- redr::calc_metrics(obs=PAR_liana_canopy%>%pull(Reflectance_median),
                                 sim=PAR_liana_canopy%>%pull(median))

NIR_tree_canopy <- data_canopy %>% filter(scenario == "low" & wavelength>800)
metrics_NIR_tree_canopy <- redr::calc_metrics(obs=NIR_tree_canopy%>%pull(Reflectance_median),
                                 sim=NIR_tree_canopy%>%pull(median))

NIR_liana_canopy <- data_canopy %>% filter(scenario == "high" & wavelength>800)
metrics_NIR_liana_canopy <- redr::calc_metrics(obs=NIR_liana_canopy%>%pull(Reflectance_median),
                                 sim=NIR_liana_canopy%>%pull(median))

table3 <- 
  rbind(
    signif(
      rbind(cbind(rbind(metrics_PAR_tree,metrics_PAR_liana),rbind(metrics_NIR_tree,metrics_NIR_liana)),
            rbind(c(metrics_delta_PAR,metrics_delta_NIR),
                  c(metrics_delta_PAR_m,metrics_delta_NIR_m))),3),
    signif(
      rbind(cbind(rbind(metrics_PAR_tree_canopy,metrics_PAR_liana_canopy),rbind(metrics_NIR_tree_canopy,metrics_NIR_liana_canopy)),
            rbind(c(metrics_delta_PAR_canopy,metrics_delta_NIR_canopy),
                  c(metrics_delta_PAR_canopy_m,metrics_delta_NIR_canopy_m))),3))

write.csv(x= table3,file = "~/data/RTM/table3.csv")
