rm(list = ls())

library(cowplot)

###########################################################################################################
# Subplot A
load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")

# dis2find_prim <- c('b1Bl_large','b2Bl_large','SLA','orient_factor','clumping_factor','Nlayers',
#                    'Cab','Car','Cw','Cm')
# npft <- 2
# Names <- c("ssigma","soil.moisture",rep(dis2find_prim,npft))

raw_data <- Pall %>% filter(Param!= "ssigma") %>% ungroup() %>% dplyr::select(-c(sample,rel_value))

# Integrate Kalacska posterior
Kalacska_post <- readRDS(file = "~/data/RTM/current_samples_Kalacska_b.rds")
Nensemble <- raw_data %>% filter(ref ==  "Kalacska" & pft == "Liana_optical" & Param == "Cab") %>% nrow()

samples_Kalacska <- getSample(Kalacska_post,numSamples = Nensemble)

samples_Kalacska_L <- melt(samples_Kalacska[,3:12]) %>% rename(Param = Var2) %>% dplyr::select(Param,value) %>% mutate(pft = "Liana_optical",
                                                                                                                     ref = "Kalacska") %>%
filter(Param == "Cab" | Param == "Car" | Param == "Cw" | Param == "Cm" | Param == "Nlayers")

samples_Kalacska_T <- melt(samples_Kalacska[,13:22]) %>% rename(Param = Var2) %>% dplyr::select(Param,value) %>% mutate(pft = "Tree_optical",
                                                                                                                ref = "Kalacska") %>%
filter(Param == "Cab" | Param == "Car" | Param == "Cw" | Param == "Cm" | Param == "Nlayers")

raw_data <- rbind(raw_data %>% filter(ref !=  "Kalacska"),
                samples_Kalacska_L,
                samples_Kalacska_T) 

raw_data <- raw_data %>% mutate(pft = case_when(
pft == "Liana_optical" ~ "Liana",
pft == "Tree_optical" ~ "Tropical tree")) %>% mutate(Param = case_when(
  Param == "Cab" ~ "Cab [µg/cm²]",
  Param == "Car" ~ "Car [µg/cm²]",
  Param == "Cw"  ~ "Cw [cm]",
  Param == "Cm" ~ "Cm [µg/cm²]",
  Param == "Nlayers" ~ "Nlayers [-]")) %>% mutate(ref = case_when(
    ref == "Sanchez_PNM" ~ "Sanchez (PNM)",
    ref == "Sanchez_FTS" ~ "Sanchez (FTS)",
    ref == "Castro_PNM" ~ "Castro (PNM)",
    ref == "Castro_FTS" ~ "Castro (FTS)",
    ref ==  "Guzman" ~  "Guzman",
    ref ==  "Kalacska" ~  "Kalacska"))

# Recalculate p-value (homoegeneize)
references <- unique(raw_data$ref)
PFTs <- unique(raw_data$pft)
params <- unique(raw_data$Param)
data.filtered <- p.values_all <- data.frame()

N = 30
for (reference in references){
  print(reference)
  for (param in params){
    
    p.value <- c()
    
    # we calculate 100 times because it depends on the sampling
    for (i in seq(1,100)){
      
      param_temp <- data.frame()
      for (PFT in PFTs){
  
        temp <- raw_data %>% filter(ref == reference & pft == PFT & Param == param)
        temp.filtered <- temp %>% filter(value <= (mean(value) + 1.5*IQR(value)) & value >= (mean(value) - 1.5*IQR(value)))
        
        if (i == 1){
          data.filtered <- rbind(data.filtered,
                                 temp.filtered)
        }
        
        param_temp <- rbind(param_temp,
                            temp.filtered %>% sample_n(N) %>% dplyr::select(Param,value,pft,ref))
        
      }
      Test <- kruskal.test(formula = value ~ as.factor(pft), data = param_temp)
      p.value <- c(p.value,Test$p.value)
    }
    
    p.values_all <- rbind(p.values_all,
                          data.frame(Param = param,
                                     p.values = mean(p.value),
                                     ref = reference,
                                     pft = "Liana",
                                     val_max = max(param_temp %>% pull(value))))
  }
}

temp <- p.values_all %>% group_by(Param) %>% summarise(value_max_all = max(val_max))

p.values_all <- p.values_all %>% mutate(signif = paste(case_when(
p.values < 0.001 ~ '***',
p.values < 0.01 ~ '**',
p.values < 0.05 ~ '*',
TRUE ~ 'N.S.'))) %>% left_join(temp)


subplotA <-
  ggplot(data = data.filtered ,
       aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="",
       x = "",
       fill = "Growth form") +
  theme_ridges(font_size = 13) + 
  geom_text(data = p.values_all,
            aes(x = 1.05*value_max_all, y = ref, label = signif),color = 'black',
            nudge_y = 0.25,
            size = 2) +
  facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
  theme_bw() +
  theme(axis.text = element_text(size=12))

###########################################################################################################
# Subplot B  

Nensemble = 1000
best_sets <- c("~/data/RTM/current_samples_Kalacska_b.rds",
               "~/data/RTM/current_samples_Marvin_b.rds",
               "~/data/RTM/current_samples_Sanchez_b.rds")

references <- c("Kalacska","Marvin","Sanchez")

ensemble_all <- data.frame()

for (iref in seq(1,length(best_sets))){
  
  best_set <- readRDS(best_sets[iref])
  ensemble <- getSample(best_set,numSamples = Nensemble)
  samples_L <- melt(ensemble[,3:12]) %>% rename(Param = Var2) %>% dplyr::select(Param,value) %>% mutate(pft = "Liana_optical",
                                                                                                        ref = references[iref]) %>%
    filter(!(Param == "Cab" | Param == "Car" | Param == "Cw" | Param == "Cm" | Param == "Nlayers"))
  
  samples_T <- melt(ensemble[,13:22]) %>% rename(Param = Var2) %>% dplyr::select(Param,value) %>% mutate(pft = "Tree_optical",
                                                                                                        ref = references[iref]) %>%
    filter(!(Param == "Cab" | Param == "Car" | Param == "Cw" | Param == "Cm" | Param == "Nlayers"))
  
  ensemble_all <- rbind(ensemble_all,
                        samples_L,samples_T)
                                                                                                                         
}


# Calculate p.values 
references <- unique(ensemble_all$ref)
PFTs <- unique(ensemble_all$pft)
params <- unique(ensemble_all$Param)
data.filtered <- p.values_all <- data.frame()

N = 30
for (reference in references){
  print(reference)
  for (param in params){
    
    p.value <- c()
    
    # we calculate 100 times because it depends on the sampling
    for (i in seq(1,100)){
      
      param_temp <- data.frame()
      for (PFT in PFTs){
        
        temp <- ensemble_all %>% filter(ref == reference & pft == PFT & Param == param)
        temp.filtered <- temp %>% filter(value <= (mean(value) + 1.5*IQR(value)) & value >= (mean(value) - 1.5*IQR(value)))
        
        if (i == 1){
          data.filtered <- rbind(data.filtered,
                                 temp.filtered)
        }
        
        param_temp <- rbind(param_temp,
                            temp.filtered %>% sample_n(N) %>% dplyr::select(Param,value,pft,ref))
        
      }
      Test <- kruskal.test(formula = value ~ as.factor(pft), data = param_temp)
      p.value <- c(p.value,Test$p.value)
    }
    
    p.values_all <- rbind(p.values_all,
                          data.frame(Param = param,
                                     p.values = mean(p.value),
                                     ref = reference,
                                     pft = "Liana",
                                     val_max = max(param_temp %>% pull(value))))
  }
}

temp <- p.values_all %>% group_by(Param) %>% summarise(value_max_all = max(val_max))

p.values_all <- p.values_all %>% mutate(signif = paste(case_when(
  p.values < 0.001 ~ '***',
  p.values < 0.01 ~ '**',
  p.values < 0.05 ~ '*',
  TRUE ~ 'N.S.'))) %>% left_join(temp) %>% mutate(Param = case_when(
    Param == "b1Bl_large" ~ "b1Bl",
    Param == "b2Bl_large" ~ "b2Bl",
    Param == "SLA" ~ "SLA [m²/kg]",
    Param == "orient_factor" ~ "ω [-]",
    Param == "clumping_factor" ~ "Ω [-]"))

subplot.b <- ensemble_all %>% mutate(pft = case_when(
  pft == "Liana_optical" ~ "Liana",
  pft == "Tree_optical" ~ "Tropical tree")) %>% mutate(Param = case_when(
    Param == "b1Bl_large" ~ "b1Bl",
    Param == "b2Bl_large" ~ "b2Bl",
    Param == "SLA" ~ "SLA [m²/kg]",
    Param == "orient_factor" ~ "ω [-]",
    Param == "clumping_factor" ~ "Ω [-]"))

subplotB <- 
  ggplot(data = subplot.b ,
       aes(x = value, y = ref, fill = as.factor(pft))) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="",
       x = "",
       fill = "Growth form") +
  theme_ridges(font_size = 13) + 
  geom_text(data = p.values_all,
            aes(x = 1.08*value_max_all, y = ref, label = signif),color = 'black',
            nudge_y = 0.25,
            size = 2) +
  facet_wrap(Param ~ .,scales = "free_x",nrow = 1) +
  theme_bw() + 
  theme(axis.text = element_text(size=12))

plot_grid(subplotA,subplotB,
          align = c("hv"),
          nrow = 2)

ggsave(plot = last_plot(),
       dpi = 300,
       width = 30,
       height = 15,
       units = "cm",
       file = "~/data/RTM//Figure3.png")

