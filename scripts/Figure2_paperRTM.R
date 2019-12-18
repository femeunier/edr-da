rm(list = ls())

library(dplyr)
library(ggplot2)
library(pracma)
library(Hmisc)
library(cowplot)

load(file = "~/data/RTM/Inverse_leaf_spectrum2.Rdata")
All_leaf_spectra <- readRDS(file= "~/data/RTM/All_leaf_spectra.rds") %>% mutate(wavelength = as.integer(round(wavelength)))

# Subplot a 
Select = "Castro_PNM"           #  "Guzman"      "Kalacska"    "Sanchez_PNM" "Sanchez_FTS" ""  "Castro_FTS" "Castro_PNM"

subplota <-  ensemble_posterior_all %>% ungroup %>% filter(ref == Select) %>% mutate(pft = case_when(
  pft == "Liana_optical" ~ "Liana",
  pft == "Tree_optical" ~ "Tropical tree"))

data_subplota <- All_leaf_spectra %>% ungroup %>% filter(ref == Select) %>% mutate(pft = case_when(
  pft == "Liana_optical" ~ "Liana",
  pft == "Tree_optical" ~ "Tropical tree"))

subplotA <-
  ggplot(data =subplota,
         aes(x = wavelength,
             y = median,
             ymin = alphamin,
             ymax = alphamax,
             fill = pft,
             color = pft)) +
  geom_ribbon(alpha = 0.5,linetype = 0) +
  # geom_line() +
  geom_line(data = data_subplota,
            aes(x = wavelength,
                y = Reflectance_median,
                ymin = Reflectance_min,
                ymax = Reflectance_max,
                fill = pft,
                color = pft),linetype = 2) + 
  labs(x = "Wavelength [nm]",
       y = "Leaf reflectance [-]",
       fill = "Growth form",
       colour = "Growth form") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  scale_y_continuous(limits = c(0,0.65),expand = c(0.001, 0.001),breaks = seq(0,0.6,0.2)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        legend.title=element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

# Subplot b 

posterior <- ensemble_posterior_all %>% dplyr::select(wavelength,alphamin,alphamax,median,ref,pft) %>% 
  rename(posterior_alphamin=alphamin,
         posterior_median=median,
         posterior_alphamax=alphamax)

ready.for.plot <- All_leaf_spectra %>%left_join(posterior,by=c("ref","pft","wavelength"))

ready.for.plot_filt <- ready.for.plot %>% mutate(group_wd = round(wavelength*2/100)-8) %>%
  group_by(ref,group_wd) %>% sample_n(1) %>% mutate(pft = case_when(
    pft == "Liana_optical" ~ "Liana",
    pft == "Tree_optical" ~ "Tropical tree")) %>% ungroup %>% mutate(ref = case_when(
      ref == "Sanchez_PNM" ~ "Sanchez (PNM)",
      ref == "Sanchez_FTS" ~ "Sanchez (FTS)",
      ref == "Castro_PNM" ~ "Castro (PNM)",
      ref == "Castro_FTS" ~ "Castro (FTS)",
      ref ==  "Guzman" ~  "Guzman",
      ref ==  "Kalacska" ~  "Kalacska"))


model.ind.fit <- ready.for.plot %>% group_by(ref,pft) %>% summarise(slope = coef(lm(formula = posterior_median ~ Reflectance_median))[2],
                                                   r2 = summary(lm(formula = posterior_median ~ Reflectance_median))$adj.r.squared,
                                                   SSE = c(crossprod(lm(formula = posterior_median ~ Reflectance_median)$residuals))/length(lm(formula = posterior_median ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))
             
model.ind.fit_VIS <- ready.for.plot %>% filter(wavelength < 800) %>%group_by(ref,pft) %>% summarise(slope = coef(lm(formula = posterior_median ~ Reflectance_median))[2],
                                                                    r2 = summary(lm(formula = posterior_median ~ Reflectance_median))$adj.r.squared,
                                                                    SSE = c(crossprod(lm(formula = posterior_median ~ Reflectance_median)$residuals))/length(lm(formula = posterior_median ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))

model.ind.fit_NIR <- ready.for.plot %>% filter(wavelength > 800) %>%group_by(ref,pft) %>% summarise(slope = coef(lm(formula = posterior_median ~ Reflectance_median))[2],
                                                                                                    r2 = summary(lm(formula = posterior_median ~ Reflectance_median))$adj.r.squared,
                                                                                                    SSE = c(crossprod(lm(formula = posterior_median ~ Reflectance_median)$residuals))/length(lm(formula = posterior_median ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))      

subplotB <-
  ggplot(ready.for.plot_filt,aes(x=Reflectance_median,y=posterior_median,
                                 colour = as.factor(pft),
                                 shape = as.factor(ref))) +
  geom_errorbar(aes(ymin = posterior_alphamin,ymax = posterior_alphamax)) +
  geom_errorbarh(aes(xmin = Reflectance_alphamin,xmax = Reflectance_alphamax)) +
  geom_point(size = 2) +
  # geom_smooth(method = "lm",level = 1-alpha) + 
  labs(x = "Observed leaf reflectance [-]",
       y = "Simulated leaf reflectance [-]",
       shape = "Reference",
       color = "Growth form") +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  scale_y_continuous(limits = c(0,0.65),expand = c(0.001, 0.001),breaks = seq(0,0.6,0.2)) +
  scale_x_continuous(limits = c(0,0.65),expand = c(0.001, 0.001),breaks = seq(0,0.6,0.2)) +
  scale_shape_manual(values=c(0,1,2,15,16,17)) +
  geom_abline(slope = 1, intercept = 0,colour = 'black',linetype=2) + 
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

####################################################################################################
# Subplot C

load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")
# load(file = "~/data/RTM/Inverse_canopy_spectrum_CM1.Rdata")
# load(file = "~/data/RTM/Inverse_canopy_spectrum_kalackaonly.Rdata")

All_canopy_spectra <- readRDS(file= "~/data/RTM/All_canopy_spectra.rds")
All_canopy_spectra <- All_canopy_spectra %>% group_by(ref,scenario,wavelength) %>% summarise(Reflectance_min = min(reflectance),
                                                                                             Reflectance_max = max(reflectance),
                                                                                             Reflectance_median = median(reflectance),
                                                                                             Reflectance_alphamin = min(reflectance),
                                                                                             Reflectance_alphamax = max(reflectance)) %>% ungroup()

Select = "Marvin"  # "Kalacska" "Marvin"   "Sanchez" 

subplotc <- model_ensemble_all %>% filter(ref == Select) %>% ungroup() %>% mutate(scenar = case_when(
  scenar == "low" ~'Low',
  scenar == "high" ~'High'))

data_subplotc <- All_canopy_spectra %>% filter(ref == Select) %>% ungroup() %>% mutate(scenario = case_when(
  scenario == "low" ~'Low',
  scenario == "high" ~'High')) 


scenarios <- unique(data_subplotc$scenario)
for (scenario in scenarios){
  wavelengths <- subplotc %>% filter(scenar == scenario) %>% pull(wavelength)  
  
  for (iwl in seq(1,length(wavelengths)-1)){
    if ((wavelengths[iwl + 1]-wavelengths[iwl]) > 100){
      subplotc <- rbind(subplotc,
                        data.frame(scenar = scenario,
                                   wavelength = c(wavelengths[iwl] + 1,wavelengths[iwl + 1] -1),
                                   rmin = NA,
                                   rmax = NA,
                                   alphamin = NA,
                                   alphamax = NA,
                                   median = NA,
                                   rbest = NA,
                                   ref = Select))
    }
  }
  
  wavelengths <- data_subplotc %>% filter(scenario == scenario) %>% pull(wavelength)  
  
  for (iwl in seq(1,length(wavelengths)-1)){
    if ((wavelengths[iwl + 1]-wavelengths[iwl]) > 100){
      data_subplotc <- rbind(data_subplotc,
                             data.frame(scenario = scenario,
                                        wavelength = c(wavelengths[iwl] + 1,wavelengths[iwl + 1] -1),
                                        Reflectance_min = NA,
                                        Reflectance_max = NA,
                                        Reflectance_alphamin = NA,
                                        Reflectance_alphamax = NA,
                                        Reflectance_median = NA,
                                        ref = Select))
    }
  }
  
}

subplotC <-
  ggplot() +
  geom_ribbon(data = subplotc,
              aes(x = wavelength,
                  ymin = rbest - (median-alphamin),
                  ymax = rbest +  (alphamax-median),
                  # ymin = alphamin,
                  # ymax = alphamax,
                  fill = scenar,
                  color = scenar),alpha = 0.5,size = 0.5,linetype=0) +
  geom_line(data = subplotc,
            aes(x = wavelength,
                y = rbest,
                color = scenar),linetype = 1) +
  geom_line(data = data_subplotc,
            aes(x = wavelength,
                y = Reflectance_median,
                color = scenario),linetype = 2) + 
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  scale_y_continuous(limits = c(0,0.65),expand = c(0.001, 0.001)) +
  scale_x_continuous(expand = c(0.001, 0.001)) +
  theme_bw() + 
  labs(x = "Wavelength [nm]",
       y = "Canopy reflectance [-]",
       fill = "Liana infestation",
       colour = "Liana infestation") +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


###############################################################################################
# Subplot d

model.results <- model_ensemble_all %>% dplyr::select(scenar,wavelength,median,alphamin,alphamax,ref,rbest) %>% rename(scenario = scenar)
model.results <- model.results %>% mutate(alphamin = rbest - (median - alphamin),
                                          alphamax = rbest + (alphamax - median),
                                          median = rbest) 

raw_data <- All_canopy_spectra %>% dplyr::select(ref,wavelength,scenario,Reflectance_median,Reflectance_alphamin,Reflectance_alphamax)

data_Marvin_sd <- read.csv(file = "~/data/RTM/marvin_sd.csv",header = TRUE)
data_Marvin <- raw_data %>% filter(ref == "Marvin")

scenarios <- unique(data_Marvin$scenario)
N <- c(150,40)
compt = 1
for (scenar in scenarios){
  wavelengths <- data_Marvin %>% filter(scenario == scenar) %>% pull(wavelength)
  R <-  data_Marvin %>% filter(scenario == scenar) %>% pull(Reflectance_median)
  data_temp <- data_Marvin_sd %>% filter(scenar == scenar)
  sd_temp <- abs(R - approxExtrap(data_temp %>% pull(wl),data_temp %>% pull(Rsd),wavelengths)$y)
  data_Marvin[data_Marvin$scenario == scenar,c("Reflectance_alphamin","Reflectance_alphamax")] <- 
    #cbind(R - sd_temp, R + sd_temp)
    cbind(R - 1.96*(sd_temp)/sqrt(N[compt]),R + 1.96*(sd_temp)/sqrt(N[compt]))
  compt = compt + 1
}

data <- rbind(raw_data %>% filter(ref != "Marvin"),data_Marvin)

subplot.d <- data %>% left_join(model.results,by=c("ref","scenario","wavelength")) %>% mutate(group_wd = round(wavelength/100)-3) %>%
  group_by(ref,group_wd) %>% sample_n(1) %>% mutate(scenario = case_when(
    scenario == "low" ~ "Low",
    scenario == "high" ~ "High"))


SSB <- data %>% left_join(model.results,by=c("ref","scenario","wavelength")) 

model.ind.fit <- SSB %>% group_by(ref,scenario) %>% summarise(slope = coef(lm(formula =  rbest ~ Reflectance_median))[2],
                                                              m = mean(Reflectance_median),
                                                                    r2 = summary(lm(formula =  rbest ~ Reflectance_median))$adj.r.squared,
                                                                    SSE = c(crossprod(lm(formula =  rbest ~ Reflectance_median)$residuals))/length(lm(formula = rbest ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))

model.ind.fit_VIS <- SSB %>% filter(wavelength < 800) %>%group_by(ref,scenario) %>% summarise(slope = coef(lm(formula =  rbest ~ Reflectance_median))[2],
                                                                                                    r2 = summary(lm(formula =  rbest ~ Reflectance_median))$adj.r.squared,
                                                                                                    SSE = c(crossprod(lm(formula =  rbest ~ Reflectance_median)$residuals))/length(lm(formula = rbest ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))

model.ind.fit_NIR <- SSB %>% filter(wavelength > 800) %>%group_by(ref,scenario) %>% summarise(slope = coef(lm(formula = rbest ~ Reflectance_median))[2],
                                                                                                    r2 = summary(lm(formula = rbest ~ Reflectance_median))$adj.r.squared,
                                                                                                    SSE = c(crossprod(lm(formula = rbest ~ Reflectance_median)$residuals))/length(lm(formula = rbest ~ Reflectance_median)$residuals)) %>%
  mutate(RMSE = sqrt(SSE))      


subplotD <-
  ggplot(subplot.d,aes(x=Reflectance_median,
                       y=median,
                       colour = as.factor(scenario),
                       shape = as.factor(ref))) +
  geom_errorbar(aes(ymin = alphamin,ymax = alphamax)) +
  geom_errorbarh(aes(xmin = Reflectance_alphamin,xmax = Reflectance_alphamax)) +
  geom_point(size = 2) +
  labs(x = "Observed canopy reflectance [-]",
       y = "Simulated canopy reflectance [-]",
       shape = "Reference",
       color = "Liana infestation") +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  scale_y_continuous(limits = c(0,0.65),expand = c(0.001, 0.001),breaks = seq(0,0.6,0.2)) +
  scale_x_continuous(limits = c(0,0.65),expand = c(0.001, 0.001),breaks = seq(0,0.6,0.2)) +
  scale_shape_manual(values=c(15,16,17)) +
  geom_abline(slope = 1, intercept = 0,colour = 'black',linetype=2) + 
  theme_bw() + 
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

plot_grid(subplotA,subplotB,subplotC,subplotD,
          align = c("hv"),
          nrow = 2)

ggsave(plot = last_plot(),
       dpi = 300,
       width = 30,
       height = 15,
       units = "cm",
       file = "~/data/RTM/Figure2.png")

