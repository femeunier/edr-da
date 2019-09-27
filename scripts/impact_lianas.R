#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#

# Library
library(redr)
library(PEcAn.all)
library(rlist)
library(hdf5r)
library(ggplot2)
library(pracma)
library(dplyr)

# Inputs
pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou.xml"
# pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_BCI.xml"

edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt") 
names <- c('Lianas','No_lianas')
list_these_pfts <- list(c(2,3,4,17),c(2,3,4,3))
rerun <- TRUE                                              # Rerun PEcAn workflow?
Lambdas <- 400:2500

# Output
output_RTM <- list ()

############################################################################################################
settings <- PEcAn.settings::read.settings(pecanxmlfile)

# PFT table
df_PFT <- data.frame(PFTnum = as.numeric(unlist(lapply(lapply(settings$pfts, "[[","constants"),"[[","num"))),names = as.character(unlist(lapply(settings$pfts, "[[","name"))))
df_PFT <- df_PFT %>% arrange(PFTnum)
############################################################################################################
# Run ED2 simulation
if (rerun){
  # Clean output folders
  system(paste("rm -rf",paste0(settings$outdir,"/*")))
  redr::run_workflow_pecan(pecanxmlfile)
}

############################################################################################################
# Meta-analysis

# Meta-analysis optical traits
# Prospect version 5
# Version 5 N, Cab, Car, Cw, Cm

default <- c(1.4,73.2,15.6,0.016,0.0105)
spectra_list_all <- optical_param <- list()

for (i in seq(df_PFT$PFTnum)){
  
  PFT_name <- as.character(df_PFT[i,"names"])
  postdist <- file.path(settings$outdir,"pft",PFT_name,"post.distns.MA.Rdata")
  
  optical_param[[PFT_name]] <- default
  
  if (file.exists(postdist)){
    load(postdist)
    variables <- rownames(post.distns)

    if ("Nlayers" %in% variables){
      mean_N <- distn.stats(post.distns["Nlayers","distn"], post.distns["Nlayers","parama"], post.distns["Nlayers","paramb"])[1]
      optical_param[[PFT_name]][1] <-mean_N
    }
    
   if (all(c("chlorophyll_a","chlorophyll_b") %in% variables)){
      mean_chla <- distn.stats(post.distns["chlorophyll_a","distn"], post.distns["chlorophyll_a","parama"], post.distns["chlorophyll_a","paramb"])[1]
      mean_chlb <- distn.stats(post.distns["chlorophyll_b","distn"], post.distns["chlorophyll_b","parama"], post.distns["chlorophyll_b","paramb"])[1]
      optical_param[[PFT_name]][2] <- mean_chla + mean_chlb
    }
    
    if ("carotenoids" %in% variables){
      mean_car <- distn.stats(post.distns["carotenoids","distn"], post.distns["carotenoids","parama"], post.distns["carotenoids","paramb"])[1]
     optical_param[[PFT_name]][3] <- mean_car
    }
    
    if ("carotenoids" %in% variables){
      mean_car <- distn.stats(post.distns["carotenoids","distn"], post.distns["carotenoids","parama"], post.distns["carotenoids","paramb"])[1]
      optical_param[[PFT_name]][3] <- mean_car
    }
    
    if ("Cw" %in% variables){
      mean_cw <- distn.stats(post.distns["Cw","distn"], post.distns["Cw","parama"], post.distns["Cw","paramb"])[1]
      optical_param[[PFT_name]][4] <- mean_cw
    }
    
    if ("SLA" %in% variables){
      mean_car <- 1/distn.stats(post.distns["SLA","distn"], post.distns["SLA","parama"], post.distns["SLA","paramb"])[1]/10
      optical_param[[PFT_name]][5] <- mean_car
    }
  }

  spectra_list_all[[PFT_name]] <-  PEcAnRTM::prospect(optical_param[[PFT_name]], version = "5")
}

# Meta-analysis Structural traits
trait_values_all <- list(
  Early_Hannes = list(clumping_factor = 0.8, b1Bl_large = 0.0096),
  Mid_Hannes = list(),
  Late_Hannes = list(clumping_factor = 0.8),
  Liana_optical = list())

############################################################################################################
# Plot MA optical values 

par(mfrow=c(1,2))
plot(Lambdas,spectra_list_all[[4]][,1],type='l',col='darkblue',ylim = c(0,0.5),ylab = "Reflectance")
lines(Lambdas,spectra_list_all[[3]][,1],col='darkgreen')

plot(Lambdas,spectra_list_all[[4]][,2],type='l',col='darkblue',ylim = c(0,0.5),ylab = "Transmittance")
lines(Lambdas,spectra_list_all[[3]][,2],col='darkgreen')

############################################################################################################
# RTM runs

data <- data.frame()

# Inputs
start.date <- settings$run$start.date
date <- ISOdatetime(
  lubridate::year(start.date),
  lubridate::month(start.date),
  lubridate::mday(start.date),
  12, 00, 00, tz = "UTC")

runfile <- file.path(settings$rundir,"runs.txt")
runs <- readLines(runfile)
rundir <- file.path(settings$rundir,runs[1])
outdir <- file.path(settings$outdir,"out",runs[1])

output_dir <- file.path(outdir,"RTM")
dir.create(output_dir)

ED2IN_file <- file.path(rundir,"ED2IN")
ed2in <- read_ed2in(ED2IN_file)

history.path <- dirname(ed2in$SFILOUT)
ed2in$MAXPATCH <- 0
ed2in$MAXCOHORT <- 0
# crown mod config !
ed2in$CROWN_MOD <- 0
PFT_ref <- ed2in$INCLUDE_THESE_PFT

historyfile <- file.path(outdir,paste('history','S',toString(lubridate::year(start.date)),
                                      sprintf("%02d", lubridate::month(start.date)),
                                      sprintf("%02d", lubridate::day(start.date)),
                                      '000000-g01.h5',sep='-'))

for (ifile in seq(names)){
  
  these_pfts <- list_these_pfts[[ifile]]
  pftsnames2keep <- as.character(df_PFT[match(these_pfts,df_PFT$PFTnum),'names'])
  
  trait_values <- trait_values_all[pftsnames2keep]
  names(trait_values) <- names(trait_values_all)
  
  spectra_list <- spectra_list_all[pftsnames2keep]
  names(spectra_list) <- names(spectra_list_all)
  
  output_dir <- file.path(outdir,"RTM",names[ifile])
  dir.create(output_dir)
  
  edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, date,patches=TRUE)

  
  output_RTM[[names[ifile]]] <- 
    PEcAnRTM::EDR(img_path = NULL,
                  ed2in_path = edr_ed2in,
                  spectra_list = spectra_list,
                  trait.values = trait_values,
                  edr_exe_path =  edr_exe_path,
                  patches = TRUE)

  
  if (ifile == 1){
    Npatches <- nrow(output_RTM[[ifile]])
    data_temp <- data.frame(lambda = rep(Lambdas,Npatches),
                            Reflectance = as.vector(t(output_RTM[[ifile]])),
                            patches = sort(rep(1:Npatches,length(Lambdas))),
                            scenario = names[ifile])
    data <- data_temp
  } else {
    data_temp <- data.frame(Reflectance = as.vector(t(output_RTM[[ifile]])))
    names(data_temp) <- c(paste0(c('Reflectance'),ifile))
    data <- cbind(data,data_temp)
  }
}

historyfile_out <- file.path(outdir,paste('history','S',toString(lubridate::year(start.date)),
                                              sprintf("%02d", lubridate::month(start.date)),
                                              sprintf("%02d", lubridate::day(start.date)),
                                              '000000-g01.h5',sep='-'))

hfile <- hdf5r::H5File$new(historyfile_out)

PACO_N=(readDataSet(hfile[['PACO_N']]))
LAI=readDataSet(hfile[['LAI_CO']])
HITE=readDataSet(hfile[['HITE']])
PFT=readDataSet(hfile[['PFT']])
DBH=readDataSet(hfile[['DBH']])
AREA=readDataSet(hfile[['AREA']])
NPATCHES=readDataSet(hfile[['NPATCHES_GLOBAL']])
nplant=readDataSet(hfile[['NPLANT']])

hfile$close_all()

ipaconow  <- rep(sequence(NPATCHES),times=PACO_N)
nplant_pa <- tapply(X=nplant, INDEX=ipaconow, FUN=sum)

LAI_temp <- LAI
LAI_tot <- tapply( X = LAI_temp, INDEX = ipaconow, FUN   = sum)#end tapply
filter <- rep(0,length(LAI))
filter[PFT == 17] <- 1
LAI_liana <- tapply( X = LAI_temp*filter, INDEX = ipaconow, FUN   = sum)#end tapply
order <- sort(LAI_liana, index.return=TRUE)$ix

rbPal <- colorRampPalette(c('Blue','red'))
C <- rbPal(Npatches)
Colors <- C
Colors[order] <- C

par(mfrow=c(1,1))
matplot(Lambdas,t(output_RTM[[1]]),type='l',col=Colors,lty=1)
matlines(Lambdas,t(output_RTM[[2]]),type='l',lty=2,col=Colors)

data <- data %>% mutate(error = Reflectance - Reflectance2,
                        relerror = (Reflectance - Reflectance2)/Reflectance)

ggplot(data = data, aes(Reflectance,Reflectance2, color = lambda)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  geom_abline(slope = 1,col = 'black',linetype = 2)

library(ggalt)
library(viridis)

ggplot(data %>% filter(abs(relerror)>0), aes(lambda, relerror)) + 
  geom_point(shape=16, size=0.25, show.legend = FALSE) +
  stat_bkde2d(aes(fill=..level..), geom="polygon") +
  theme_minimal() +
  scale_fill_viridis() 

ggplot(data, aes(lambda, error)) + 
  geom_point(shape=16, size=0.25) +
  stat_bkde2d(grid_size=c(128, 128), geom="polygon", aes(fill=..level..)) +
  scale_fill_viridis() +
  theme_minimal() +
  theme(panel.grid=element_blank())