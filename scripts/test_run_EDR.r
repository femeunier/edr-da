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

# Inputs
pecanxmlfile <- "/home/carya/R/inputs/RTM_workflow_Paracou.xml"
edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt") 
names <- c('Lianas','No_lianas')
keep_these_pfts <- list(c(2,3,4,17),c(2,3,4))
rerun <- TRUE # PEcAn workflow
Lambdas <- 400:2500

# Output
output_RTM <- list ()

############################################################################################################
settings <- PEcAn.settings::read.settings(pecanxmlfile)

# PFT table
df_PFT <- data.frame(PFTnum = as.numeric(unlist(lapply(lapply(settings$pfts, "[[","constants"),"[[","num"))),names = unlist(lapply(settings$pfts, "[[","name")))
############################################################################################################
# Meta-analysis

spectra_list_all <- list(Early_Hannes = PEcAnRTM::prospect(c(1.4, 73.2, 15.6,  0.016, 0.0105), version = "5"),
                     Mid_Hannes   = PEcAnRTM::prospect(c(1.4, 73.2, 15.6,  0.016, 0.0105), version = "5"),
                     Late_Hannes  = PEcAnRTM::prospect(c(1.4, 73.2, 15.6,  0.016, 0.0105), version = "5"),
                     Liana_BCI    = PEcAnRTM::prospect(c(1.4, 77, 15.8,  0.018, 0.0093), version = "5")) # Version 5 N, Cab, Car, Cw, Cm 

trait_values_all <- list(
  Early_Hannes = list(clumping_factor = 0.8, b1Bl_large = 0.0096),
  Mid_Hannes = list(),
  Late_Hannes = list(clumping_factor = 0.8),
  Liana_BCI = list())

############################################################################################################
# Plot MA optical values 

par(mfrow=c(1,2))
plot(Lambdas,spectra_list[[4]][,1],type='l',col='darkblue',ylim = c(0,0.75))
lines(Lambdas,spectra_list[[3]][,1],col='darkgreen')

plot(Lambdas,spectra_list[[4]][,2],type='l',col='darkblue',ylim = c(0,0.75))
lines(Lambdas,spectra_list[[3]][,2],col='darkgreen')

############################################################################################################
# RTM runs

data <- data.frame()

# Run ED2 simulation
if (rerun){
  # Clean output folders
  system(paste("rm -rf",paste0(settings$outdir,"/*")))
  redr::run_workflow_pecan(pecanxmlfile)
}

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
  
  these_pfts <- keep_these_pfts[[ifile]]
  pftsnames2keep <- as.character(df_PFT[df_PFT$PFTnum %in% these_pfts,'names'])
  remove_pfts <- df_PFT$PFTnum[!(df_PFT$PFTnum %in% these_pfts)]
  
  trait_values <- trait_values_all[pftsnames2keep]
  spectra_list <- spectra_list_all[pftsnames2keep]
  
  output_dir <- file.path(outdir,"RTM",names[ifile])
  dir.create(output_dir)

  historyfile_out <- file.path(output_dir,paste('history','S',toString(lubridate::year(start.date)),
                                                 sprintf("%02d", lubridate::month(start.date)),
                                                 sprintf("%02d", lubridate::day(start.date)),
                                                 '120000-g01.h5',sep='-'))
  
  ed2in$INCLUDE_THESE_PFT <- these_pfts   
  edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, date,patches=TRUE)
  
  remove_pft_from_history_file(file = historyfile,
                               file_mod = historyfile_out,
                               PFT2remove = remove_pfts)

  output_RTM[[names[ifile]]] <- 
    PEcAnRTM::EDR(img_path = NULL,
                  ed2in_path = edr_ed2in,
                  spectra_list = spectra_list,
                  trait.values = trait_values,
                  edr_exe_path =  edr_exe_path,
                  patches = TRUE)

  Npatches <- nrow(output_RTM[[ifile]])
  
  if (ifile == 1){
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

par(mfrow=c(1,1))
matplot(Lambdas,t(output_RTM[[1]]),type='l',col='black',lty=1)
matlines(Lambdas,t(output_RTM[[2]]),type='l',lty=2,col='red')


ggplot(data = data, aes(Reflectance,Reflectance2, color = lambda)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  geom_abline(slope = 1,col = 'black',linetype = 2)


#################################################################################################
# Analyze outputs

Npatches <- as.numeric(readLines(con = file.path(output_dir,
                                                 "patches/Npatches.dat")))


file <- "/home/carya/output/PEcAn_99000000002/out/SA-median/RTM/history-S-2004-01-01-120000-g01_keep.h5"
hfile <- hdf5r::H5File$new(file)

PACO_N_noL=(readDataSet(hfile[['PACO_N']]))
LAI=readDataSet(hfile[['LAI_CO']])
HITE=readDataSet(hfile[['HITE']])
PFT=readDataSet(hfile[['PFT']])
DBH=readDataSet(hfile[['DBH']])
LAI_PY=readDataSet(hfile[['LAI_PY']])
AREA_noL=readDataSet(hfile[['AREA']])
NPATCHES_GLOBAL_noL=readDataSet(hfile[['NPATCHES_GLOBAL']])
nplant=readDataSet(hfile[['NPLANT']])
Bleaf=readDataSet(hfile[['BLEAF']])
SLA=readDataSet(hfile[['SLA']])
CA=readDataSet(hfile[['CROWN_AREA_CO']])

hfile$close_all()

ipaconow  <- rep(sequence(NPATCHES_GLOBAL_noL),times=PACO_N_noL)
nplant_pa <- tapply(X=nplant, INDEX=ipaconow, FUN=sum)

CA_adult <- CA
CA_adult[DBH>10] <- 0.
CA_m=tapply(X=CA_adult*LAI,INDEX=ipaconow,FUN=sum)

LAI_temp <- LAI
LAI_tot <- tapply( X = LAI_temp, INDEX = ipaconow, FUN   = sum)#end tapply
filter <- rep(0,length(LAI))
filter[PFT == 17] <- 1
LAI_liana <- tapply( X = LAI_temp*filter, INDEX = ipaconow, FUN   = sum)#end tapply
order <- sort(LAI_tot, index.return=TRUE)$ix

rbPal <- colorRampPalette(c('grey','black'))
C <- rbPal(NPATCHES_GLOBAL_noL)
Colors_noL <- C
Colors_noL[order] <- C

pos_red <- 290
pos_nir <- 400

Red <- OP_noL[,pos_red]
NIR <- OP_noL[,pos_nir]
NDVI <- (NIR-Red)/(NIR+Red)

plot.new()
par(mfrow=c(1,1))
matplot(400:2500,t(OP_noL),col=Colors_noL,type='l',xlim=c(400,2500),ylim=c(0,0.75),lty=1)
abline(v=pos_red+400,col='red')
abline(v=pos_nir+400,col='red',lty=2)

plot.new()
par(mfrow=c(1,2))
plot(LAI_tot,Red,xlab='LAI')
plot(LAI_tot,NIR,xlab='LAI')

par(mfrow = c(1,1))
NDVI_sort <- sort(NDVI,index.return = TRUE)
plot(NDVI_sort$x,LAI_tot[NDVI_sort$ix],type='p')

LAI_liana_min <- which.min(LAI_liana)
OP_noL_minL <- t(array(rep(OP_noL[LAI_liana_min,],Npatches),c(2101,10)))
matplot(400:2500,t(OP_noL-OP_noL_minL),col=Colors_noL,type='l',xlim=c(400,2500),ylim=c(-0.25,0.25),lty=1)
abline(v=pos_red+400,col='red')
abline(v=pos_nir+400,col='red',lty=2)
