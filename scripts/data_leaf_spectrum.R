# Kalacska et al.

rm(list=ls())

library(dplyr)
library(Hmisc)
library(reshape2)
library(pracma)

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
  temp <- melt((Reflectance[,columns])) %>% rename(Reflectance = value) %>% dplyr::select(Reflectance)
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
# Sanchez 2009
# 
# rm(list=ls())
# 
# sites <- c("FTS","PNM")
# data2read <- c("/home/carya/data/RTM/data/Fig5_castro2009_FTS2.csv","/home/carya/data/RTM/data/Fig5_castro2009_PNM.csv")
# 
# 
# PFTs <- c("Tree_optical","Liana_optical")
# C <- c("#137300","#1E64C8")
# plot(NA,NA,xlim=c(400,1000),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")
# 
# for (isite in seq(sites)){
#   data_R <- list()
#   df <- data.frame()
#   
#   datatemp <- read.csv(data2read[isite]) %>% mutate(GF = case_when(GF == "T" ~ "Tree_optical",
#                                                                    GF == "L" ~ "Liana_optical"))
# 
#   for (ipft in seq(PFTs)){
# 
#     pft <- PFTs[ipft]
#    
#     temp <- datatemp %>% filter(GF == pft) %>% dplyr::select(wl,R)
#     pos <- sort(temp[,1],index.return = TRUE)
#     
#     if (isite ==1){
#       temp[,1] <- temp[,1] -50
#     }
#     R <- rbind(c(450,temp[1,2]),
#                temp[pos$ix,],
#                c(1000,temp[length(pos$x),2]))
#     R_interp <- interp1(R[,1],R[,2],450:1000)
# 
# 
#     data_R[[pft]] <- cbind(450:1000,R_interp)
# 
#     lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
#     names(data_R[[pft]]) <- c("wavelength","reflectance")
# 
#     df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
#                               Reflectance = data_R[[pft]] [,2],
#                               pft = pft))
#   }
#   saveRDS(df,file = paste0("~/data/RTM/Figure6_sanchez2009_",sites[isite],".rds"))
# }


rm(list=ls())

sites <- c("FTS","PNM")
data2read <- "/home/carya/data/RTM/data/Figure6_sanchez2009.csv"
data <- read.csv(data2read)

dataT <- rbind(read.csv("/home/carya/data/RTM/data/Figure6_sanchez2009T_PNM.csv") %>% mutate(ref = "PNM"),
               read.csv("/home/carya/data/RTM/data/Figure6_sanchez2009T_FTS.csv") %>% mutate(ref = "FTS")) %>% mutate(GF = case_when(
                 GF == "L" ~ 'Liana_optical',
                 GF == "T" ~ 'Tree_optical'))

PFTs <- c("Tree_optical","Liana_optical")
C <- c("#137300","#1E64C8")


for (isite in seq(sites)){
  data_R <- list()
  df <- data.frame()
  plot(NA,NA,xlim=c(400,1000),ylim=c(0,1),xlab="Wavelength",ylab="Reflectance")
  for (ipft in seq(PFTs)){

    pft <- PFTs[ipft]
    columns <- (isite -1)*4 + (((ipft-1)*2+1):((ipft-1)*2+2))
    temp <- data[!is.na(data[,columns[1]]),c(columns)]
    pos <- sort(temp[,1],index.return = TRUE)
    A <- rbind(c(450,0),
               temp[pos$ix,],
               c(900,temp[length(pos$x),2]))
    A_interp <- interp1(A[,1],A[,2],450:900)

    temp <- dataT %>% filter(ref == sites[isite] & GF == pft) %>% dplyr::select(wl,L)
    dataT_temp <- rbind(c(450,temp[1,"L"]),
                        temp,
                        c(900,temp[nrow(temp),"L"]))
    T_interp <- interp1(dataT_temp[,1],dataT_temp[,2],450:900)

    data_R[[pft]] <- cbind(450:900,1-A_interp-T_interp)
    if (isite == 2 & ipft == 1) data_R[[pft]] [,1] = data_R[[pft]] [,1]
    lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
    names(data_R[[pft]]) <- c("wavelength","reflectance")
    # lines(data_R[[pft]] [,1],T_interp,col = C[ipft])
    
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

# 
# All_leaf_spectra <- All_leaf_spectra %>% filter(!(ref == "Sanchez_FTS" & wavelength < 500 | ref == "Sanchez_FTS" & wavelength > 650 & wavelength <700 |
#                                                   ref == "Sanchez_PNM" & wavelength < 500 | ref == "Sanchez_PNM" & wavelength > 650 & wavelength <700 ))


All_leaf_spectra2 <- All_leaf_spectra %>% ungroup() %>% mutate(pft = case_when(
  pft == "Tree_optical" ~ 'Tree',
  pft == "Liana_optical" ~ 'Liana'
))  %>% mutate(ref = case_when(
  ref == "Castro_FTS" ~ 'Castro (FTS)',
  ref == "Castro_PNM" ~ 'Castro (PNM)',
  ref == "Sanchez_FTS" ~ 'Sanchez (FTS)',
  ref == "Sanchez_PNM" ~ 'Sanchez (PNM)',
  ref == "Guzman" ~ 'Guzmán',
  TRUE ~ ref
)) %>% group_by(ref,pft)

ggplot(data=All_leaf_spectra2,
       aes(x = wavelength,
           y = Reflectance_median,
           ymin = Reflectance_min,
           ymax = Reflectance_max,
           group = pft,
           fill = pft,
           color = pft)) +
  # geom_ribbon(alpha = 0.5,linetype = 0) +
  geom_line() +
  facet_wrap(ref ~ .,scales = "free_x") +
  labs(y = "Reflectance [-]",
       x = "Wavelength [nm]",
       colour = "Growth form") +
  scale_color_manual(values = c("#1E64C8","#137300")) +
  scale_fill_manual(values = c("#1E64C8","#137300")) +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


ggsave(filename = "~/data/RTM/All_leaf_spectra.png",
       plot = last_plot(),
       width = 20, height = 10,
       dpi = 300,units = "cm")

saveRDS(All_leaf_spectra,file= "~/data/RTM/All_leaf_spectra.rds")

#############################################################
# Table data

bands_m <- c(550,600,400,800,1500,2000,680)
bands_M <- c(600,800,700,1400,1800,2500,700)
band_names <- c('green','Red edge','visible','NIR','SWIR','SWIR2','Rededge2')
studies <- unique(All_leaf_spectra$ref)

differences <- matrix(NA,length(studies),length(band_names))
rownames(differences) <- studies
colnames(differences) <- band_names
  
for (istudy in seq(studies)){
  for (iband in seq(band_names)){
    
    Lianas_temp <- All_leaf_spectra %>% filter(ref == studies[istudy] & pft == "Liana_optical" & wavelength >= bands_m[iband] &  wavelength <= bands_M[iband]) %>% ungroup() %>% dplyr::select(wavelength,Reflectance_median)
    
    if (nrow(Lianas_temp)>0){
      Lianas <- rbind(c(bands_m[iband],Lianas_temp[1,2] %>% pull()),
                      Lianas_temp,
                      c(bands_M[iband],Lianas_temp[nrow(Lianas_temp),2]%>% pull()))
      L_wl <- Lianas[,1]%>%pull()
      L_R <- Lianas[,2]%>%pull()
      Liana_interp <- interp1(L_wl,L_R,bands_m[iband]:bands_M[iband])
      
      Trees_temp <- All_leaf_spectra %>% filter(ref == studies[istudy] & pft == "Tree_optical" & wavelength >= bands_m[iband] &  wavelength <= bands_M[iband]) %>% ungroup() %>% dplyr::select(wavelength,Reflectance_median)
      Trees <- rbind(c(bands_m[iband],Trees_temp[1,2] %>% pull()),
                      Trees_temp,
                      c(bands_M[iband],Trees_temp[nrow(Trees_temp),2]%>% pull()))
      T_wl <- Trees[,1]%>%pull()
      T_R <- Trees[,2]%>%pull()
      Trees_interp <- interp1(T_wl,T_R,bands_m[iband]:bands_M[iband])
      
      if (band_names[iband] == "Red edge"){
        if (nrow(Lianas_temp)>0 & nrow(Trees_temp)>0 
            # & min(max(L_wl[2:(length(L_wl)-1)]),max(T_wl[2:(length(T_wl)-1)])) - bands_M[iband] > -100 &
            # max(min(L_wl[2:(length(L_wl)-1)]),min(T_wl[2:(length(T_wl)-1)])) - bands_m[iband] < 100
        ){
          red_edge_L <- mean(Liana_interp[(length(Liana_interp)-10):length(Liana_interp)]) - mean(Liana_interp[1:10])
          red_edge_T <- mean(Trees_interp[(length(Trees_interp)-10):length(Trees_interp)]) - mean(Trees_interp[1:10])
          
          if (mean(red_edge_T)>mean(red_edge_L)){
            differences[istudy,iband] <- "T"
          } else{
            differences[istudy,iband] <- "L"      
          }
        } 
        
      } else {
        if (nrow(Lianas_temp)>0 & nrow(Trees_temp)>0 
            # & min(max(L_wl[2:(length(L_wl)-1)]),max(T_wl[2:(length(T_wl)-1)])) - bands_M[iband] > -100 &
            # max(min(L_wl[2:(length(L_wl)-1)]),min(T_wl[2:(length(T_wl)-1)])) - bands_m[iband] < 100
        ){
          if (mean(Trees_interp)>mean(Liana_interp)){
            differences[istudy,iband] <- "T"
          } else{
            differences[istudy,iband] <- "L"      
          }
        } 
      }
    }
  }
}

