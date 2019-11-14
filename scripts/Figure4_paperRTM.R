rm(list = ls())

load(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")
UA.step1 <- rbind(sensitivities_all %>% select(ref,pft,param,par.var) %>% mutate(OP_variable = "PAR"),
                  sensitivities_all %>% select(ref,pft,param,par.var) %>% mutate(OP_variable = "NIR"))


load(file = "~/data/RTM/Inverse_canopy_spectrum.Rdata")
UA.step2 <- model_sensitivities_all %>% dplyr::select(ref,OP_variable,pft,par.var,variable) %>%
  rename(param = variable) %>% filter(OP_variable == "Spectrum")


# To do:
  # - step1: OP_variable
  # - step2: NIR_liana
#######################################################################################################
# Subplot A

subplot.a <-
  UA.step1 %>% group_by(OP_variable,param,pft) %>% mutate(par.var = 100*par.var) %>%
  summarise(par.var_m = mean(par.var),
            par.var_sd = sd(par.var)) %>% filter(param != "ssigma") %>% mutate(
              pft = case_when(
                pft == "Liana_optical" ~ "Liana",
                pft == "Tree_optical"  ~ "Tree"))

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
  facet_wrap(~ OP_variable,scales = "free",nrow = 2) +
  theme_bw() + 
  scale_fill_manual(values =  as.character(df_PFT$Col)) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  coord_flip() +
  labs(y = "Partial variance (%)",
       x = "",
       fill = "Growth form",
       colour = "Growth form") +
  theme(axis.text = element_text(size=12),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

#######################################################################################################
# Subplot B

