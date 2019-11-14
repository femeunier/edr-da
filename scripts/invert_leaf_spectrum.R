rm(list=ls())

library(cowplot)
library(BayesianTools)
library(dplyr)
library(redr)
library(PEcAn.all)
library(pracma)
library(reshape2)
library(ggplot2)
library(ggridges)

# load(file="~/data/RTM/invert_leaf_liana.RData")

All_leaf_spectra <- readRDS(file= "~/data/RTM/All_leaf_spectra.rds")

iter  <- 25000
nrChains  <- 2
wl.max <- 2500
wl.min  <- 350
alpha  <- 0.05

Nrun_prospect <- 1000

files <- c("Figure1_Guzman.rds","Figure1_kalacska.rds",
           "Figure6_sanchez2009_PNM.rds","Figure6_sanchez2009_FTS.rds",
           "Figures4and5_castro_PNM.rds","Figures4and5_castro_FTS.rds")
Names <- c("Guzman","Kalacska",
           "Sanchez_PNM","Sanchez_FTS","Castro_PNM","Castro_FTS")

Colors <- c("#137300","#1E64C8")

use_meta.analysis <- FALSE
direct <- NULL
df_PFT <- data.frame(names = c("Liana_optical","Tree_optical"),PFTnum = c(17,3),Col = Colors)
df_PFT <- df_PFT %>% arrange(PFTnum)%>% mutate(Col = Colors)
df_PFT <- df_PFT %>% arrange(names)

npft <- nrow(df_PFT)

settings_MCMC <- list(iterations = iter, nrChains = nrChains)

PDA <- list()
best_set <- best_run <- list()

N=10000

# Define likelihood
create_likelihood <- function(observed, waves) {
  function(params) {
    
    # PROSPECT parameters
    Cab <- params[1]
    Car <- params[2]
    Cw <- params[3]
    Cm <- params[4]
    Nlayers <- params[5]
    
    ssigma <- params[6]
    
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
    # Calculate likelihood
    ll <- sum(dnorm(reflectance_waves, observed, ssigma, log = TRUE))
    return(ll)
  }
}

############################################################################################
# Loop over papers

marginalPlots <- spectra_post <- performance_prospect_plot <- PDA_all <- list()

performance_all <- best_set_all <- posterior_all <- ensemble_posterior_all <-
  ensemble_posterior_all2 <- data.frame()

model_performance <- sensitivities_all <- p.values_all <- data.frame()

for (ifile in seq(files)){
  data <- load_rds(file.path("~/data/RTM/",files[ifile]))
  
  #############################################################################################
  # Loop over PFTs
  
  # Uniform prior
  Prospect_param_names <- c("Ca","Cb","Car","Cw","Cm","Nlayers","ssigma")
  Prospect_param_names_short <- c("Cab","Car","Cw","Cm","Nlayers","ssigma")
  pft_lowers <- c(chlorophyll_a = 0,chlorophyll_b = 0, carotenoids = 0,Cw = 0, SLA = 1, Nlayers = 1, ssigma = 0)
  pft_uppers <-  c(chlorophyll_a = 100,chlorophyll_b = 50, carotenoids = 50,Cw = 0.1, SLA = 100, Nlayers = 5, ssigma = 1)
  
  dis2find <- c('chlorophyll_a','chlorophyll_b','carotenoids','Cw','SLA','Nlayers',"ssigma")
  
  for (ipft in seq(npft)){
    current_pft <- as.character(df_PFT$names[ipft])
    
    # Define priors
    
    # From posterior
    if ((!is.null(direct))) {
      postfile <- file.path(direct,"pft",current_pft,"post.distns.Rdata")
      if (use_meta.analysis){
        load(postfile)
        distns <- post.distns
        Distributions <- rownames(post.distns)
      } else {
        priorfile <- file.path(direct,"pft",current_pft,"prior.distns.Rdata")
        load(priorfile)
        distns <- prior.distns
        Distributions <- rownames(prior.distns)
      }
    } else {
      Distributions <- NULL
    }
    
    sampler <- matrix(NA,N,length(dis2find))
    lower <- upper <- best <- rep(NA,length(dis2find))
    
    for (idis in seq(dis2find)){
      current_dis = dis2find[idis]
      pos <- which(Distributions == current_dis)
      
      if (isempty(pos)){
        unif_distribution <- list(distn="unif",parama=pft_lowers[current_dis],paramb=pft_uppers[current_dis])
        sampling <- get.sample(unif_distribution,n=N)
      } else {
        sampling <- get.sample(distns[pos,],n=N)
      }
      
      if (current_dis == 'SLA'){
        sampling = 1/sampling/10 # unit conversions
      }
      sampler[,idis] <- sampling
      lower[idis] <- min(sampler[,idis])
      upper[idis] <- max(sampler[,idis])
      best[idis]  <- median(sampler[,idis])
    }
    
    colnames(sampler) <- names(lower) <- names(upper) <- names(best) <- Prospect_param_names
    
    # combine Ca and Cb
    sampler[,"Ca"] <- sampler[,"Ca"] + sampler[,"Cb"]
    colnames(sampler)[colnames(sampler)=="Ca"] = "Cab"
    sampler <- sampler[,!(colnames(sampler) %in% "Cb")]
    
    lower["Ca"] <- lower["Ca"] + lower["Cb"]
    lower <- lower[!(names(lower) %in% c("Cb"))]
    names(lower)[(names(lower) %in% c("Ca"))] <- "Cab"
    
    upper["Ca"] <- upper["Ca"] + upper["Cb"]
    upper <- upper[!(names(upper) %in% c("Cb"))]
    names(upper)[(names(upper) %in% c("Ca"))] <- "Cab"
    
    best["Ca"] <- best["Ca"] + best["Cb"]
    best <- best[!(names(best) %in% c("Cb"))]
    names(best)[(names(best) %in% c("Ca"))] <- "Cab"
    
    # 1/as.numeric(config_temp$SLA)/10*2
    prior <- createPriorDensity(sampler, method = "multivariate", eps = 1e-10,
                                lower = lower, upper = upper, best = best)
    
    # Read Leaf spectra 
    temp <- data %>% filter(pft==current_pft & wavelength>=wl.min & wavelength<=wl.max) %>% dplyr::select(c('wavelength','Reflectance'))
    
    observation <- temp[["Reflectance"]] 
    waves <- temp[["wavelength"]]
    
    likelihood <- create_likelihood(observation, waves)
    
    # Run inversion
    setup <- BayesianTools::createBayesianSetup(likelihood, prior, parallel = FALSE)
    samples <- BayesianTools::runMCMC(setup,settings = settings_MCMC)
    samples <- BayesianTools::runMCMC(samples, settings = settings_MCMC)
    
    for (ichain in seq(nrChains)){
      samples[[ichain]][["setup"]][["names"]] <- Prospect_param_names_short
    }
    
    PDA[[current_pft]] <- samples
    
  }
  
  ##############################################################
  # Plot
  PDA_results <- PDA
  PDA_all[[Names[ipft]]] <- PDA
  
  best_set <- best_run <- list()
  
  rsquare <- rmse <- slopes <- rep(NA,npft)
  prospect_performance <- ensemble_posterior <- data.frame()
  posterior_dis <- posterior_dis_small <- data.frame()
  
  for (ipft in seq(npft)){
    
    current_pft <- as.character(df_PFT$names[ipft])
    
    temp <- data %>% filter(pft==current_pft & wavelength>=wl.min & wavelength<=wl.max) %>% dplyr::select(c('wavelength','Reflectance'))
    
    observation <- temp %>% group_by(wavelength) %>% dplyr::summarize(R_m = mean(Reflectance),R_sd = sd(Reflectance))
    waves <- observation$wavelength
    
    samples <- PDA_results[[current_pft]]
    
    MAP_samples <- MAP(samples)$parametersMAP
    names(MAP_samples) <- Prospect_param_names_short
    best_set[[current_pft]]<- c(MAP_samples['Nlayers'],MAP_samples['Cab'],MAP_samples['Car'],MAP_samples['Cw'],MAP_samples['Cm'])
    best_run[[current_pft]] <- PEcAnRTM::prospect(best_set[[current_pft]], version = "5")
    
    
    best_set_all <- rbind(best_set_all,
                          as.data.frame(melt(MAP_samples)) %>% mutate(pft = current_pft,param = names(MAP_samples),data=Names[ifile]))
    
    best_run_interp <- interp1(x=PEcAnRTM::wavelengths(best_run[[current_pft]]),
                               y = as.vector(best_run[[current_pft]][,1]),
                               xi = waves)
    
    C <- Colors[ipft]
    
    current_model <- data.frame(obs = observation$R_m, mod = best_run_interp , pft = current_pft, data = Names[ifile])
    
    lm_pft <- lm(data = current_model,formula = mod ~ obs)
    rsquare[ipft] <- summary(lm_pft)$adj.r.squared
    rmse[ipft] <- sqrt(c(crossprod(lm_pft$residuals))/(length(lm_pft$residuals)-1))
    coef_pft <- coef(lm_pft)
    slopes[ipft] <- coef_pft[2]
    
    prospect_performance <- rbind(prospect_performance,current_model)
    
    
    temp_data <- All_leaf_spectra %>% filter(ref == Names[ifile]) %>% ungroup()
    temp_data_pft <- temp_data %>% filter(pft == current_pft)
    current_waves <- temp_data_pft %>% select(wavelength) %>% pull()
      
    posteriorMat <- getSample(samples, numSamples = Nrun_prospect,
                              parametersOnly = TRUE)
    
    colnames(posteriorMat) <- Prospect_param_names_short
    rownames(posteriorMat) <- 1:nrow(posteriorMat)
    
    posteriorMat_small <- getSample(samples, numSamples = 30,
                              parametersOnly = TRUE)
    
    colnames(posteriorMat_small) <- Prospect_param_names_short
    rownames(posteriorMat_small) <- 1:nrow(posteriorMat_small)
    
    posterior_all <- rbind(posterior_all,
                           as.data.frame(melt(posteriorMat) %>% rename(Param = Var2,sample = Var1) %>% mutate(pft = current_pft,
                                                                         site = Names[ifile])))
    
    posterior_dis <- rbind(posterior_dis, 
                           (melt(posteriorMat) %>% rename(Param = Var2) %>% dplyr::select(c(Param,value)) %>% mutate(pft = current_pft)))
    
    posterior_dis_small <- rbind(posterior_dis_small, 
                           (melt(posteriorMat_small) %>% rename(Param = Var2) %>% dplyr::select(c(Param,value)) %>% mutate(pft = current_pft)))
    
    
    temp_ensemble_prospect <- matrix(NA,nrow(posteriorMat),nrow(best_run[[current_pft]]))
    colnames(temp_ensemble_prospect) <- seq(400,2500)
    
    model_OP_ensemble_spectrum <- model_OP_ensemble_PAR <- model_OP_ensemble_NIR <- c()
    
    for (irun in seq(nrow(posteriorMat))){
      current_parameter_set <- c(posteriorMat[irun,'Nlayers'],posteriorMat[irun,'Cab'],posteriorMat[irun,'Car'],posteriorMat[irun,'Cw'],posteriorMat[irun,'Cm'])
      current_model_output <- PEcAnRTM::prospect(current_parameter_set, version = "5")
      temp_ensemble_prospect[irun,] <- current_model_output[,1]
      
      model_OP_ensemble_spectrum[irun] <- mean(interp1(x = 400:2500,
                                         y =  as.vector(current_model_output[,1]),
                                         xi = current_waves))
      
      model_OP_ensemble_PAR[irun] <- mean(interp1(x =  400:2500,
                                                       y =  as.vector(current_model_output[,1]),
                                                       xi = 400:800))
      
      model_OP_ensemble_NIR[irun] <- mean(interp1(x = 400:2500,
                                                  y =  as.vector(current_model_output[,1]),
                                                  xi = 801:2500))
    }
    
    lmdf <- as.data.frame(cbind(posteriorMat,model_OP_ensemble_spectrum,model_OP_ensemble_PAR,model_OP_ensemble_NIR))
    
    # Spectrum
    anovobj<-aov(lm(as.formula(paste0("model_OP_ensemble_spectrum"," ~ .^2")),data=lmdf))
    allssq<-summary(anovobj)[[1]][,2]
    varn<-names((anovobj)$coefficients)[-1]
    vars <- varn[1:ncol(posteriorMat)]
    
    ##########
    sensivities<-sapply(vars,function(nn){
      
      SSQ <- sum(allssq,na.rm = T)
      SSQ_d <- sum(allssq[which(nn == varn)],na.rm = T)
      SSQ_i <- sum(allssq[which(!(nn == varn) & grepl(nn, varn, fixed=T))],na.rm = T)
      return((SSQ_d+(SSQ_i/2))/SSQ)},USE.NAMES = T) 
    
    sensitivities_all <- rbind(sensitivities_all,
                               melt(sensivities) %>% rename(par.var = value) %>% 
                                 mutate(pft = current_pft,
                                        ref = Names[ifile],
                                        param = names(sensivities),
                                        OP_variable = "Spectrum"))
    
    # PAR
    anovobj<-aov(lm(as.formula(paste0("model_OP_ensemble_PAR"," ~ .^2")),data=lmdf))
    allssq<-summary(anovobj)[[1]][,2]
    varn<-names((anovobj)$coefficients)[-1]
    vars <- varn[1:ncol(posteriorMat)]
    
    sensivities<-sapply(vars,function(nn){
      
      SSQ <- sum(allssq,na.rm = T)
      SSQ_d <- sum(allssq[which(nn == varn)],na.rm = T)
      SSQ_i <- sum(allssq[which(!(nn == varn) & grepl(nn, varn, fixed=T))],na.rm = T)
      return((SSQ_d+(SSQ_i/2))/SSQ)},USE.NAMES = T) 
    
    sensitivities_all <- rbind(sensitivities_all,
                               melt(sensivities) %>% rename(par.var = value) %>% 
                                 mutate(pft = current_pft,
                                        ref = Names[ifile],
                                        param = names(sensivities),
                                        OP_variable = "PAR"))
    
    # NIR
    anovobj<-aov(lm(as.formula(paste0("model_OP_ensemble_NIR"," ~ .^2")),data=lmdf))
    allssq<-summary(anovobj)[[1]][,2]
    varn<-names((anovobj)$coefficients)[-1]
    vars <- varn[1:ncol(posteriorMat)]
    
    sensivities<-sapply(vars,function(nn){
      
      SSQ <- sum(allssq,na.rm = T)
      SSQ_d <- sum(allssq[which(nn == varn)],na.rm = T)
      SSQ_i <- sum(allssq[which(!(nn == varn) & grepl(nn, varn, fixed=T))],na.rm = T)
      return((SSQ_d+(SSQ_i/2))/SSQ)},USE.NAMES = T) 
    
    sensitivities_all <- rbind(sensitivities_all,
                               melt(sensivities) %>% rename(par.var = value) %>% 
                                 mutate(pft = current_pft,
                                        ref = Names[ifile],
                                        param = names(sensivities),
                                        OP_variable = "NIR"))
    
    ensemble_posterior <- rbind(ensemble_posterior,
                                as.data.frame(
                                  melt(temp_ensemble_prospect) %>% dplyr::select(Var2, value) %>% rename(wavelength = Var2, reflectance = value) %>%
                                    group_by(wavelength) %>% summarise(
                                      rmin = min(reflectance, na.rm = TRUE),
                                      rmax = max(reflectance, na.rm =
                                                   TRUE),
                                      alphamin = quantile(reflectance, alpha /
                                                            2, na.rm = TRUE),
                                      alphamax = quantile(reflectance, 1 -
                                                            alpha / 2, na.rm = TRUE),
                                      median = median(reflectance, na.rm =
                                                        TRUE)
                                    )
                                ) %>% mutate(pft = current_pft) #mutate
    ) #rbind
    
  }
  
  model_performance <- rbind(model_performance,
                             data.frame(rmse = rmse,r2 = rsquare, slope = slopes, pft = as.character(df_PFT$names)) %>%
                               mutate(ref = Names[ipft]))
  
  
  ensemble_posterior_all <- rbind(ensemble_posterior_all,
                                  ensemble_posterior %>% mutate(ref = Names[ifile]) %>%
                                    filter(wavelength >= min(temp_data %>% select(wavelength) %>% pull()) & 
                                             wavelength <= max(temp_data %>% select(wavelength) %>% pull())))
  
  ensemble_posterior_all2 <- rbind(ensemble_posterior_all,
                                  ensemble_posterior %>% mutate(ref = Names[ifile]))
  
  performance_all <- rbind(performance_all,prospect_performance)
  
  #################################################################################################
  # Compare mean 
  
  posterior_dis <- posterior_dis %>% group_by(Param) %>% mutate(rel_value = (value - min(value))/(max(value)-min(value))) %>% ungroup()
  posterior_dis_small <- posterior_dis_small %>% group_by(Param) %>% mutate(rel_value = (value - min(value))/(max(value)-min(value))) %>% ungroup()
  
  
  data_traits <- data.frame()
  p_values <- rep(NA,length(Prospect_param_names_short))
  names(p_values) <- Prospect_param_names_short
  
  for (trait in Prospect_param_names_short){
  
    data_aov <- posterior_dis_small %>% filter(Param == trait)
    
    Test <- kruskal.test(formula = value ~ as.factor(pft), data = data_aov)
    p_values[trait] <- Test$p.value
  }
  
  p.values_all <- rbind(p.values_all,
                        data.frame(param = names(p_values),
                                   p.values = p_values,
                                   ref = Names[ifile]))
  
}

######################################################################################
# Plots
# OP vs data

all_OP_plots <- 
  ggplot(data = ensemble_posterior_all,
       aes(x = wavelength,
           y = median,
           ymin = alphamin,
           ymax = alphamax,
           fill = pft,
           color = pft)) +
  geom_ribbon(alpha = 0.5,linetype = 0) +
  geom_line() +
  geom_line(data = All_leaf_spectra,
             aes(x = wavelength,
                 y = Reflectance_median,
                 ymin = Reflectance_min,
                 ymax = Reflectance_max,
                 fill = pft,
                 color = pft),linetype = 2) + 
  facet_wrap(ref ~ .,scales = "free_x") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","PDA_prospect_outputs.png"),dpi = 300,
       plot = all_OP_plots,
       width = 10, height = 5)

performance_all_plot <-
  ggplot(performance_all,aes(x=obs,y=mod,colour = as.factor(pft),fill = as.factor(pft),shape = as.factor(data))) +
  # geom_smooth(method = "lm",level = 1-alpha) + 
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  scale_shape_manual(values=seq(0,length(unique(performance_all$data)))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,colour = 'black',linetype=2) + 
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","performance_all.png"),dpi = 300,
       plot = performance_all_plot,
       width = 10, height = 5)

Pall <- 
  posterior_all %>% group_by(Param,pft,site) %>% mutate(rel_value = (value - min(value))/(max(value)-min(value))) %>%
  rename(ref = site)

p.values_all <- 
  p.values_all %>% mutate(signif = paste("Kruskal-Wallis test:",case_when(
    p.values < 0.001 ~ '***',
    p.values < 0.01 ~ '**',
    p.values < 0.05 ~ '*',
    TRUE ~ 'N.S.'))) %>% mutate(signif2 = case_when(
      p.values < 0.001 ~ '***',
      p.values < 0.01 ~ '**',
      p.values < 0.05 ~ '*',
      TRUE ~ 'N.S.'))

all_marg_plots <- 
  ggplot(data = Pall,
         aes(x = rel_value, y = Param, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="") +
  scale_x_continuous(limits = c(-0.25,1.5)) +
  theme_ridges(font_size = 13) + 
  geom_text(data = p.values_all %>% mutate(pft = "Liana_optical"),
            aes(x = 1.25, y = param, label = signif2),color = 'black',
            nudge_y = 0.5,
            size = 2) + 
  facet_wrap(ref ~ .) +
  theme_bw()
  
ggsave(filename = file.path("~/data/RTM/","PDA_prospect_marginalDis.png"),dpi = 300,
       plot = all_marg_plots,
       width = 20, height = 10, units = "cm")


posterior_all_m <- posterior_all %>% group_by(Param,pft,site) %>% summarise(value_m = mean(value),
                                                                            value_sd = sd(value),
                                                                            value_max = max(value)) %>%
  rename(ref = site,param = Param)

p.values_all <- left_join(p.values_all,posterior_all %>% group_by(Param,site) %>% filter(pft == "Tree_optical") %>%
                            summarise(value_max = 0.95*(mean(value) + sd(value))) %>%
                            rename(param = Param, ref = site))
temp_value_max <- posterior_all_m %>% ungroup() %>% group_by(param) %>% summarise(value_max_all = max(value_max))
p.values_all2 <- p.values_all %>% left_join(temp_value_max) %>% rename(Param = param)

all_marg_plots2 <- 
  ggplot(data = Pall %>% filter(Param!= "ssigma"),
         aes(x = value, y = ref, fill = pft)) +
  geom_density_ridges(alpha= 0.5) +
  scale_color_manual(values = as.character(df_PFT$Col)) +
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  labs(y="") +
  theme_ridges(font_size = 13) + 
  geom_text(data = p.values_all2 %>% mutate(pft = "Liana_optical") %>% filter(Param!= "ssigma"),
            aes(x = 1.15*value_max_all, y = ref, label = signif2),color = 'black',
            nudge_y = 0.5,
            size = 2) +
  facet_wrap(Param ~ .,scales = "free_x") +
  theme_bw()

ggsave(filename = file.path("~/data/RTM/","marginalPlot2.png"),dpi = 300,
       plot = all_marg_plots2,
       width = 10, height = 5)


barplot_parameters <-
  ggplot() + 
  geom_errorbar(data = posterior_all_m %>% group_by(param) %>% mutate(rel_value = value_m/max(value_max)),
                aes(x = pft,ymin = 0.95*value_m,ymax=value_m+value_sd),width=.2) +
  geom_bar(data = posterior_all_m %>% group_by(param) %>% mutate(rel_value = value_m/max(value_max)),
           aes(x = pft, y = value_m,fill = pft),stat = "identity",
           position = position_dodge(width=0.2), width = 0.7) + 
  scale_fill_manual(values = as.character(df_PFT$Col)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = p.values_all %>% mutate(pft = "Tree_optical"),
            aes(x = pft, y = 0.95*value_max, label = signif2), color = 'black',
            size = 3,
            nudge_x = 0.35) + 
  ylab("") + xlab("") +
  facet_grid(param ~ ref,scales = "free_y")

ggsave(filename = file.path("~/data/RTM/","barplot_param.png"),dpi = 300,
       plot = barplot_parameters,
       width = 10, height = 5)


sensitivity_plot <- ggplot(data = sensitivities_all %>% filter(param != "ssigma"),
       aes(x = param,y = 100*par.var,color = pft, fill =  pft)) +
  geom_bar(stat = "identity",
           position = position_dodge(), width = 0.7) +
  facet_wrap(. ~ ref,scales = "free_y",nrow = 2) +
  theme_bw() + 
  scale_fill_manual(values = as.character(df_PFT$Col)) +  
  scale_color_manual(values = as.character(df_PFT$Col)) +
  coord_flip() +
  ylab("Partial variance (%)") + 
  xlab("") 

ggsave(filename = file.path("~/data/RTM/","sensitivity_plot.png"),dpi = 300,
       plot = sensitivity_plot,
       width = 10, height = 5)

save.image(file = "~/data/RTM/Inverse_leaf_spectrum.Rdata")