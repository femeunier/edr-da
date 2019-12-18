rm(list = ls())

library(dplyr)
library(hdf5r)
library(redr)
library(pracma)

crown_mod = 0
PFTselect = 1
alpha_frac = 0.5

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


##########################################################################


##########################################################################
best_set <- readRDS("~/data/RTM/current_samples_Kalacska_b.rds")
# ensemble <- getSample(best_set,numSamples = Nensemble)
params <- MAP(best_set)$parametersMAP
  
par.wl <- 400:2499 
nir.wl <- 2500 

ed2in <- PEcAn.ED2::read_ed2in(file.path(rundir,"ED2IN"))



outputs <- run_ED_RTM(rundir,modeloutdir,params,crown_mod,inventory,par.wl,nir.wl,
                      temp = FALSE,outputdir = OP_dir)

# remove weird patches
patches2remove <- c(15,42)
outputs$output_RTM[(1:Npatches %in% patches2remove) ,] <- NA

RTM <- outputs$output_RTM
Npatches <- nrow(RTM)
RTM_PAR <- matrix(RTM[,1:length(par.wl)],Npatches,length(par.wl))
RTM_NIR <- matrix(RTM[,(length(par.wl)+1):ncol(RTM)],Npatches,length(nir.wl))

par_filter <- t(matrix(rep(par.wl>=400 & par.wl <=800,Npatches),length(par.wl),Npatches))
nir_filter <- t(matrix(rep(par.wl>=801 & par.wl <=2500,Npatches),length(par.wl),Npatches))
red_filter <- t(matrix(rep(par.wl>=601 & par.wl <=700,Npatches),length(par.wl),Npatches))

PAR <- apply(RTM_PAR*par_filter,1,mean,na.rm=TRUE)
NIR <- apply(RTM_PAR*nir_filter,1,mean,na.rm=TRUE)
RED <- apply(RTM_PAR*red_filter,1,mean,na.rm=TRUE)
NDVI = (NIR-RED)/(NIR + RED)

plot(NDVI,LAItot)

par(mfrow=c(1,2))
plot(Lindex$fraction_LAI*Lindex$LAItot,(NDVI),pch=1,col="black")
plot(Lindex$LAItot,(NDVI),pch=1,col="black")

plot_patch_vertical_distribution(h5file_paste,ispatch = TRUE)

matplot(par.wl,t(RTM_PAR),type='l')