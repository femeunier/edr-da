# Kalacska et al.

rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure1_kalacska.csv"
data <- read.csv(data2read)
PFTs <- c("Tree_optical","Liana_optical")
C <- c("#137300","#1E64C8")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
df <- data.frame()
for (ipft in seq(PFTs)){
  pft <- PFTs[ipft]
  columns <- ((ipft-1)*2+1):((ipft-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[pft]] <- temp[pos$ix,]
  lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
  names(data_R[[pft]]) <- c("wavelength","reflectance")
  
  df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                            Reflectance = data_R[[pft]] [,2],
                            pft = pft))
}

saveRDS(df,file = "~/data/RTM/Figure1_kalacska.rds")

##################################################################################
# Castro et al.

rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure6_castro.csv"

N <- 3
wl.min <- 450
wl.max <- 900
wl.all <- seq(wl.min,wl.max,1)

data <- read.csv(data2read)
PFTs <- c(rep("Tree_optical",N),rep("Liana_optical",N))
C <- c(rep("#137300",N),rep("#1E64C8",N))
plot(NA,NA,xlim=c(wl.min-50,wl.max+50),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

data_R <- data_raw <- data.frame()

for (ipft in seq(PFTs)){
  pft <- PFTs[ipft]
  columns <- ((ipft-1)*2+1):((ipft-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_temp <- data.frame(wavelength = wl.all,
                          R = (approxExtrap(x = temp[pos$ix,1],
                                           y = temp[pos$ix,2],
                                           xout = wl.all))$y,
                          pft = pft)
  lines(data_temp[["wavelength"]],data_temp[["R"]],col = C[ipft])
  data_R <- rbind(data_R,data_temp)
  data_raw <- rbind(data_raw,data.frame(wavelength = temp[pos$ix,1],
                                        Reflectance = temp[pos$ix,2],
                                        pft = pft))
}

saveRDS(data_raw,file = "~/data/RTM/Figure6_castro.rds")

##################################################################################
# de Guzman et al.

rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure1_Guzman.csv"
data <- read.csv(data2read)
WL <- data[,1]
Reflectance <- data[,-1]

PFTs <- c("Tree_optical","Liana_optical")
patterns <- c("Tree*","Liana*")

C <- c("#137300","#1E64C8")

plot(NA,NA,xlim=c(400,950),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

data_R <- data.frame()
for (ipft in seq(PFTs)){
  
  pft <- PFTs[ipft]
  is_pft <- grepl(patterns[ipft],colnames(Reflectance))
  
  columns <- which(is_pft)
  temp <- melt((Reflectance[,columns])) %>% rename(Reflectance = value) %>% select(Reflectance)
  matlines(WL,data[,columns],col = C[ipft])
  
  data_R <- rbind(data_R,data.frame(wavelength = rep(WL,length(columns)),
                                    Reflectance = temp,
                                    pft = pft))
  
}

saveRDS(data_R,file = "~/data/RTM/Figure1_Guzman.rds")

##################################################################################
# Castro et al. (4 and 5)

rm(list=ls())

sites <- c("FTS","PNM")
data2read <- "/home/carya/data/RTM/data/Figures4and5_castro.csv"
data <- read.csv(data2read)
PFTs <- c("Liana_optical","Tree_optical")
C <- c("#1E64C8","#137300")
plot(NA,NA,xlim=c(400,1000),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

for (isite in seq(sites)){
  data_R <- list()
  df <- data.frame()
  for (ipft in seq(PFTs)){
    pft <- PFTs[ipft]
    columns <- (isite -1)*4 + (((ipft-1)*2+1):((ipft-1)*2+2))
    temp <- data[!is.na(data[,columns[1]]),c(columns)]
    pos <- sort(temp[,1],index.return = TRUE)
    data_R[[pft]] <- temp[pos$ix,]
    lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
    names(data_R[[pft]]) <- c("wavelength","reflectance")
    
    df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                              Reflectance = data_R[[pft]] [,2],
                              pft = pft))
  }
  saveRDS(df,file = paste0("~/data/RTM/Figures4and5_castro_",sites[isite],".rds"))
}

##################################################################################
# Castro et al. (4 and 5)

rm(list=ls())

sites <- c("FTS","PNM")
data2read <- "/home/carya/data/RTM/data/Figure6_sanchez2009.csv"
data <- read.csv(data2read)
PFTs <- c("Tree_optical","Liana_optical")
C <- c("#137300","#1E64C8")
plot(NA,NA,xlim=c(400,1000),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

for (isite in seq(sites)){
  data_R <- list()
  df <- data.frame()
  for (ipft in seq(PFTs)){
    pft <- PFTs[ipft]
    columns <- (isite -1)*4 + (((ipft-1)*2+1):((ipft-1)*2+2))
    temp <- data[!is.na(data[,columns[1]]),c(columns)]
    pos <- sort(temp[,1],index.return = TRUE)
    data_R[[pft]] <- temp[pos$ix,]
    lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
    names(data_R[[pft]]) <- c("wavelength","reflectance")
    
    df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                              Reflectance = data_R[[pft]] [,2],
                              pft = pft))
  }
  saveRDS(df,file = paste0("~/data/RTM/Figure6_sanchez2009_",sites[isite],".rds"))
}

###########################################################################################
# All together
# readRDS(file = "~/data/RTM/Figure6_castro.rds")

rm(list=ls())

alpha = 0.05

All_leaf_spectra <-
  rbind(readRDS(file = "~/data/RTM/Figure1_Guzman.rds") %>% mutate(ref = "Guzman"),
      readRDS(file = "~/data/RTM/Figure1_kalacska.rds") %>% mutate(ref = "Kalacska"),
      readRDS(file = "~/data/RTM/Figures4and5_castro_FTS.rds") %>% mutate(ref = "Castro_FTS"),
      readRDS(file = "~/data/RTM/Figures4and5_castro_PNM.rds") %>% mutate(ref = "Castro_PNM"),
      readRDS(file = "~/data/RTM/Figure6_sanchez2009_PNM.rds") %>% mutate(ref = "Sanchez_PNM"),
      readRDS(file = "~/data/RTM/Figure6_sanchez2009_FTS.rds") %>% mutate(ref = "Sanchez_FTS"))

ref_alls <- unique(All_leaf_spectra$ref)

curves_all <- c()
for (ref in (ref_alls)){
  pos_ref <- which(All_leaf_spectra$ref == ref)
  pfts <- as.character(unique(All_leaf_spectra[pos_ref,"pft"]))
  for (pft in pfts){
    
    pos <- which(All_leaf_spectra$ref == ref &
                   All_leaf_spectra$pft == pft)
    compt <- 1
  
    curves <- rep(1,length(pos))
    
    for (i in seq(1,length(pos)-1)){
      if (All_leaf_spectra$wavelength[pos[i]] > All_leaf_spectra$wavelength[pos[i+1]]){
        compt = compt +1
      }
      curves[i] <- compt
    }
  
    curves_all <- c(curves_all,curves)
  }  
}

All_leaf_spectra$curve <- curves_all


count <- All_leaf_spectra %>% group_by(ref) %>% summarise(curve_max = max(curve))

# All_leaf_spectra_single <- All_leaf_spectra %>% filter(ref %in% (count %>% filter(curve_max ==1) %>% select(ref) %>%pull()))
# All_leaf_spectra_multiple <- All_leaf_spectra %>% filter(ref %in% (count %>% filter(curve_max >1) %>% select(ref) %>%pull()))

All_leaf_spectra <- All_leaf_spectra %>% group_by(ref,pft,wavelength) %>% summarise(Reflectance_min = min(Reflectance),
                                                                         Reflectance_max = max(Reflectance),
                                                                         Reflectance_median = median(Reflectance),
                                                                         Reflectance_alphamin = quantile(Reflectance,alpha/2),
                                                                         Reflectance_alphamax = quantile(Reflectance,1-alpha/2))


ggplot(data=All_leaf_spectra,
       aes(x = wavelength,
           y = Reflectance_median,
           ymin = Reflectance_min,
           ymax = Reflectance_max,
           fill = pft,
           color = pft)) +
  geom_ribbon(alpha = 0.5,linetype = 0) +
  geom_line() +
  facet_wrap(ref ~ .,scales = "free_x") +
  scale_color_manual(values = c("#137300","#1E64C8")) +
  scale_fill_manual(values = c("#137300","#1E64C8")) +
  theme_bw()

ggsave(filename = "~/data/RTM/All_leaf_spectra.png",
       plot = last_plot(),
       width = 20, height = 10,
       dpi = 300,units = "cm")

saveRDS(All_leaf_spectra,file= "~/data/RTM/All_leaf_spectra.rds")