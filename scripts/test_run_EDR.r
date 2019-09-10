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

pecanxmlfile="/home/carya/R/inputs/RTM_workflow_Paracou.xml"
settings <- PEcAn.settings::read.settings(pecanxmlfile)
redr::run_workflow_pecan(pecanxmlfile)

# Inputs
edr_exe_path <- file.path("/home/carya/ED2/EDR/run","ed_2.1-opt") 

start.date <- settings$run$start.date
date <- ISOdatetime(
        lubridate::year(start.date),
        lubridate::month(start.date),
        lubridate::mday(start.date),
        12, 00, 00, tz = "UTC")

# Inputs
runfile <- file.path(settings$rundir,"runs.txt")
runs <- readLines(runfile)
rundir <- file.path(settings$rundir,runs[1])
outdir <- file.path(settings$outdir,runs[1])
output_dir <- '/home/carya/temp/nolianas'

ED2IN_file <- file.path(rundir,"ED2IN")
ed2in <- read_ed2in(ED2IN_file)

history.path <- dirname(ed2in$SFILOUT)
ed2in$MAXPATCH <- 0
ed2in$MAXCOHORT <- 0

# crown mod config !
ed2in$CROWN_MOD <- 0


edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, date,patches=TRUE)

trait_values <- list(
  Early_Hannes = list(clumping_factor = 0.8, b1Bl_large = 0.0096),
  Mid_Hannes = list(),
  Late_Hannes = list(clumping_factor = 0.8),
  Liana_BCI = list())

spectra_list <- list(Early_Hannes = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Mid_Hannes   = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Late_Hannes  = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Liana_BCI    = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"))

OP_noL <- PEcAnRTM::EDR(img_path = NULL,
                        ed2in_path = edr_ed2in,
                        spectra_list = spectra_list,
                        trait.values = trait_values,
                        edr_exe_path =  edr_exe_path,
                        patches = TRUE)

Npatches <- as.numeric(readLines(con = file.path(output_dir,
                                                 "patches/Npatches.dat")))

LAI_sim <- rep(0.,Npatches)
CA_sim  <- rep(0.,Npatches)
for (i in sequence(Npatches)){
  line <- readLines(con = file.path(output_dir,paste("patches/LAI",i,'.dat',sep='')))
  temp=as.numeric(strsplit(line,'\\s+')[[1]][-1])
  
  line <- readLines(con = file.path(output_dir,paste("patches/CA",i,'.dat',sep='')))
  temp_CA=as.numeric(strsplit(line,'\\s+')[[1]][-1])
  
  LAI_sim[i] <- sum(temp)
  CA_sim[i] <- mean(temp_CA)
}


file<-file.path(output_dir,paste('history','S',toString(lubridate::year(start.date)),
                                 sprintf("%02d", lubridate::month(start.date)),
                                 sprintf("%02d", lubridate::day(start.date)),
                                 '120000-g01.h5',sep='-'))

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
WAI=readDataSet(hfile[['WAI_CO']])
CA=readDataSet(hfile[['CROWN_AREA_CO']])

hfile$close_all()

ipaconow  <- rep(sequence(NPATCHES_GLOBAL_noL),times=PACO_N_noL)
nplant_pa <- tapply(X=nplant, INDEX=ipaconow, FUN=sum)

CA_adult <- CA
CA_adult[DBH>10] <- 0.
CA_m=tapply(X=CA_adult*LAI,INDEX=ipaconow,FUN=sum)

LAI_temp <- LAI
LAI_total_noL = tapply( X = LAI_temp, INDEX = ipaconow, FUN   = sum)#end tapply
order=sort(LAI_sim, index.return=TRUE)$ix

rbPal <- colorRampPalette(c('grey','black'))
C <- rbPal(NPATCHES_GLOBAL_noL)
Colors_noL <- C
Colors_noL[order] <- C

plot.new()
par(mfrow=c(1,1))
matplot(c(401:2500,2501),t(OP_noL),col=Colors_noL,type='l',xlim=c(400,2500),ylim=c(0,0.6),lty=1)

Red=OP_noL[,300]
NIR=OP_noL[,500]
NDVI=(NIR-Red)/(NIR+Red)
plot.new()
par(mfrow=c(1,2))
plot(LAI_total_noL,Red,xlab='LAI')
plot(LAI_total_noL,NIR,xlab='LAI')

par(mfrow = c(1,1))
NDVI <- sort(NDVI,index.return = TRUE)
plot(NDVI$x,LAI_total_noL[NDVI$ix],type='p')
