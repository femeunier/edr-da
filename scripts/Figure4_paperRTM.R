rm(list = ls())

library(dplyr)

ymax <- 105

load(file = "~/data/RTM/Inverse_leaf_spectrum2.Rdata")
UA.step1 <- sensitivities_all %>% dplyr::select(ref,pft,param,par.var,OP_variable) %>% 
  mutate(ref = case_when(
    ref == "Sanchez_PNM" ~ "Sanchez (PNM)",
    ref == "Sanchez_FTS" ~ "Sanchez (FTS)",
    ref == "Castro_PNM" ~ "Castro (PNM)",
    ref == "Castro_FTS" ~ "Castro (FTS)",
    ref ==  "Guzman" ~  "Guzman",
    ref ==  "Kalacska" ~  "Kalacska")) %>% filter(OP_variable %in% c("PAR","NIR"))

load(file = "~/data/RTM/Inverse_canopy_spectrum_N250.Rdata")
# load(file = "~/data/RTM/Inverse_canopy_spectrum_CM1.Rdata")

UA.step2 <- model_sensitivities_all %>% dplyr::select(ref,OP_variable,pft,par.var,variable) %>%
  rename(param = variable) %>% mutate(OP_variable = as.character(OP_variable)) %>% mutate(OP_variable = case_when(
    OP_variable == "nir_liana" ~ "NIR_Liana",
    OP_variable == "nir_tree" ~ "NIR_Tree",
    OP_variable == "PAR_liana" ~ "PAR_Liana",
    OP_variable == "PAR_tree" ~ "PAR_Tree")) %>% filter(OP_variable %in% c("NIR_Liana","NIR_Tree","PAR_Liana","PAR_Tree")) %>%
  mutate(ref = case_when(
    ref == "Sanchez_PNM" ~ "Sanchez (PNM)",
    ref == "Sanchez_FTS" ~ "Sanchez (FTS)",
    ref == "Castro_PNM" ~ "Castro (PNM)",
    ref == "Castro_FTS" ~ "Castro (FTS)",
    ref ==  "Guzman" ~  "Guzman",
    ref ==  "Kalacska" ~  "Kalacska"))  

#######################################################################################################
# Subplot A

subplot.a <-
  UA.step1 %>% group_by(OP_variable,param,pft) %>% mutate(par.var = 100*par.var) %>%
  summarise(par.var_m = mean(par.var),
            par.var_sd = sd(par.var),
            par.var_se = sd(par.var)/sqrt(length(par.var))) %>% filter(param != "ssigma") %>% mutate(
              pft = case_when(
                pft == "Liana_optical" ~ "Liana",
                pft == "Tree_optical"  ~ "Tree")) %>% mutate(OP_var = case_when(
                  OP_variable == "NIR" ~ "NIR [801-2500nm]",
                  OP_variable == "PAR" ~ "PAR [400-800nm]")) %>% mutate(col = 1)


subplotA <-
  ggplot(data = subplot.a,
       aes(x = param,
           y = par.var_m,
           ymin = 0.95*par.var_m,
           ymax = par.var_m + par.var_sd,
           color = as.factor(pft),
           fill = as.factor(pft))) +
  geom_errorbar(position=position_dodge(width=0.8),width = 0.5) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7) +
  facet_grid(OP_var ~ col,scales = "fixed") +
  coord_flip() +
  theme_bw() +     
  scale_y_continuous(limits = c(0,ymax),
                     breaks = seq(0,100,20),
                     expand = c(0.001,0.001)) +
  scale_fill_manual(values =  as.character(df_PFT$Col)) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  labs(y = "Partial variance [%]",
       x = "",
       fill = "Growth form",
       colour = "Growth form") +
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        strip.background.x = element_rect(colour=NA, fill= NA),
        strip.text.x = element_text(colour = NA))

#######################################################################################################
# Subplot B

subplot.b <- UA.step2 %>% group_by(OP_variable,pft,param) %>% mutate(par.var = 100*par.var) %>%
    summarise(par.var_m = mean(par.var),
              par.var_sd = sd(par.var),
              par.var_se = sd(par.var)/sqrt(length(par.var))) %>% ungroup() %>% filter(param != "ssigma") %>% mutate(
                pft = case_when(
                  pft == "Liana_optical" ~ "Liana",
                  pft == "Tree_optical"  ~ "Tree",
                  pft == "soil" ~ "Soil")) %>% mutate(OP_var = sub("\\_.*", "", OP_variable),
                                                               scenario = sub(".*\\_", "", OP_variable))
  
names(Values) <- c("Liana","Tree","Soil")  
Values["Soil"] <- "dimgray"

subplot.b2 <- subplot.b %>% filter(!(param %in% c("clumping_factor","SLA","b1Bleaf","b2Bleaf"))) %>%
  mutate(scenario = case_when(
    scenario == "Liana" ~ "Liana-rich patches",
    scenario == "Tree" ~ "Liana-free patches"
  )) %>% mutate(param = case_when(
    param == "soil_moisture" ~ "θsoil",
    param == "orient_factor" ~ "ω [-]",
    param == "carotenoids" ~ "Car",
    param == "Cab" ~ "Cab", 
    param == "Cm" ~ "Cm", 
    param == "Cw" ~ "Cw",
    param == "Nlayers" ~ "Nlayers")) %>% mutate(order = case_when(
      param == "θsoil" ~ 7,
      param == "ω [-]" ~ 3,
      param == "Car" ~ 2,
      param == "Cab" ~ 4, 
      param == "Cm" ~ 6, 
      param == "Cw" ~ 5,
      param == "Nlayers" ~ 1)) %>% mutate(OP_var = case_when(
        OP_var == "NIR" ~ "NIR [801-2500nm]",
        OP_var == "PAR" ~ "PAR [400-800nm]"))

subplot.b2$pft <- factor(subplot.b2$pft, levels = c("Liana", "Tree", "Soil"))

subplotB <- 
 ggplot(data = subplot.b2 %>% filter(pft != "Soil"),
         aes(x = order,
             y = par.var_m,
             ymin = 0.95*par.var_m,
             ymax = par.var_m + par.var_sd,
             color = as.factor(pft),
             fill = as.factor(pft))) +
    geom_errorbar(position=position_dodge(width=0.8),width = 0.5) +
    geom_bar(stat = "identity",
             position = position_dodge(), width = 0.7) +
    facet_grid(OP_var ~ scenario,scales = "fixed") +
    coord_flip()  +
    scale_fill_manual(values = Values) +
    scale_color_manual(values = Values) +
    scale_y_continuous(limits = c(0,ymax),
                       breaks = seq(0,100,20),
                       expand = c(0.001,0.001)) +
    scale_x_continuous(breaks = seq(1:7),
                       labels=c("Nlayers","Car","ω [-]","Cab","Cw","Cm","θsoil")) + 
    theme_bw() + 
    labs(y = "Partial variance [%]",
         x = "",
         fill = "Source",
         colour = "Source") +
    theme(axis.text = element_text(size=12),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          panel.spacing.x = unit(1.5, "lines"))

plot_grid(subplotA,subplotB,
          align = c("hv"),
          nrow = 1,rel_widths = c(1,1.5))

ggsave(plot = last_plot(),
       dpi = 300,
       width = 30,
       height = 12,
       units = "cm",
       file = "~/data/RTM/Figure4.png")

##################################################################################
# Legend

legend <- 
  ggplot(data = subplot.b2,
         aes(x = order,
             y = par.var_m,
             ymin = 0.95*par.var_m,
             ymax = par.var_m + par.var_sd,
             color = as.factor(pft),
             fill = as.factor(pft))) +
  geom_errorbar(position=position_dodge(width=0.8),width = 0.5) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7) +
  facet_grid(OP_var ~ scenario,scales = "fixed") +
  coord_flip()  +
  scale_fill_manual(values = Values) +
  scale_color_manual(values = Values) +
  scale_y_continuous(limits = c(0,ymax),
                     breaks = seq(0,100,20),
                     expand = c(0.001,0.001)) +
  scale_x_continuous(breaks = seq(1:7),
                     labels=c("Nlayers","Car","ω [-]","Cab","Cw","Cm","θsoil")) + 
  theme_bw() + 
  labs(y = "Partial variance (%)",
       x = "",
       fill = "Source",
       colour = "Source") +
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        panel.spacing.x = unit(1.5, "lines"),
        legend.position = "top")

ggsave(plot = legend,
       dpi = 300,
       width = 30,
       height = 12,
       units = "cm",
       file = "~/data/RTM/legend_Figure4.png")

