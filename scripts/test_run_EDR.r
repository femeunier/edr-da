#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#

# Library
library(rlang)
library(PEcAn.ED2)
library(PEcAnRTM)
library("hdf5r")
library("ggplot2")
library('pracma')

pecanxmlfile="/home/carya/output/PEcAn_99000000007/pecan.CONFIGS.xml"
settings <- PEcAn.settings::read.settings(pecanxmlfile)

# Inputs
edr_exe_path = file.path("/home/carya/ED2/EDR/build","ed_2.1-opt-master-5953f04a") 
pb = NULL
img_path = NULL
year_s=2004
month=1
all_WL=400:2500

date <- ISOdatetime(
        lubridate::year(paste(year_s,month,1,sep='-')),
        lubridate::month(paste(year_s,month,1,sep='-')),
        lubridate::mday(paste(year_s,month,1,sep='-')),
        12, 00, 00, tz = "UTC"
      )

################################################################################################
# Load properties

csvfile='/home/carya/data/reflectance.csv'
data=read.table(csvfile,sep=',',header=TRUE)

WL_T=data$Wavelength[!is.na(data$Wavelength)]
WL_L=data$WL[!is.na(data$WL)]
R_T=data$Rtree[!is.na(data$Rtree)]
R_L=data$Rliana[!is.na(data$Rliana)]

WL_T[length(WL_T)]<-WL_L[length(WL_L)]<-max(all_WL)

WL_interp=all_WL
Rall_T <- interp1(x=WL_T,y=R_T,xi=WL_interp)
Rall_L <- interp1(x=WL_L,y=R_L,xi=WL_interp)

#################################################################################################
# No liana simulation

# Inputs
dir='/home/carya/output/PEcAn_99000000007/run/BCI_nolianas'
odir= '/home/carya/output/PEcAn_99000000007/out/BCI_nolianas'
output_dir='/home/carya/temp/nolianas'

ED2IN_file <- file.path(dir,"ED2IN")
ed2in <- read_ed2in(ED2IN_file)
ed2in$SFILOUT <- file.path(odir,'history')
ed2in$POI_LAT <- 9
ed2in$POI_LON <- -79
ed2in$ED_MET_DRIVER_DB <- file.path("/home/carya/raw_inputED/met/BCI","BCI_DRIVER_HEADER2")
ed2in$VEG_DATABASE <- file.path("/home/carya/raw_inputED/OGE2",basename(ed2in$VEG_DATABASE))
ed2in$SOIL_DATABASE <- file.path("/home/carya/raw_inputED/FAO",basename(ed2in$SOIL_DATABASE))
ed2in$THSUMS_DATABASE <- file.path("/home/carya/raw_inputED/thermal_sums",basename(ed2in$THSUMS_DATABASE))

if (!lubridate::is.Date(date)) date <- lubridate::as_date(date)
dtime <- ISOdatetime(
  lubridate::year(date),
  lubridate::month(date),
  lubridate::mday(date),
  12, 00, 00, tz = "UTC"
)

history.path = dirname(ed2in$SFILOUT)
datetime = dtime

ed2in$INCLUDE_THESE_PFT <- c(2,3,4)
ed2in$MAXPATCH <- 0
ed2in$MAXCOHORT <- 0

# crown mod config !
ed2in$CROWN_MOD <- 0

edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, dtime,patches=TRUE)

trait_values=list(Early = list(clumping_factor=0.5,b1Bl_large=0.0096),
                  Mid = list(clumping_factor=0.5,b1Bl_large=0.0132),
                  Late = list(clumping_factor=0.5))

spectra_list <- list(Early = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Mid   = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Late  = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"))

spectra_list$Early[,1]=Rall_T
spectra_list$Mid[,1]=Rall_T
spectra_list$Late[,1]=Rall_T

OP_noL <- PEcAnRTM::EDR(img_path = img_path,
                        ed2in_path = edr_ed2in,
                        spectra_list = spectra_list,
                        trait.values = trait_values,
                        edr_exe_path =  edr_exe_path,
                        patches = TRUE)

# Change input file
# file<-file.path(output_dir,paste('history','-S-',toString(year_s), '-',sprintf("%02d", month),'-01-120000-g01.h5',sep=''))
# hfile <- hdf5r::H5File$new(file)
# 
# LeafT=(readDataSet(hfile[['LEAF_TEMP']]))
# LeafT=LeafT[1]*LeafT**0
# hfile$link_delete("LEAF_TEMP")
# hfile[["LEAF_TEMP"]] <- LeafT
# 
# WoodT=(readDataSet(hfile[['WOOD_TEMP']]))
# WoodT=WoodT[1]*WoodT**0
# hfile$link_delete("WOOD_TEMP")
# hfile[["WOOD_TEMP"]] <- WoodT
# 
# nplant=(readDataSet(hfile[['NPLANT']]))
# nplant=mean(nplant)*nplant**0
# hfile$link_delete("NPLANT")
# hfile[["NPLANT"]] <- nplant
# 
# #soilT=(readDataSet(hfile[['SOIL_TEMPK_PA']]))
# #soilT[,]=soilT[1,1]
# #hfile$link_delete("SOIL_TEMPK_PA")
# #hfile[["SOIL_TEMPK_PA"]] <- soilT
# 
# hfile$close_all()

# run model
#setwd(output_dir)
#ex <- system2(edr_exe_path, args = c("-f", edr_ed2in))

# par.wl = 400:2499 
# nir.wl = 2500
# # Read outputs
# output.path_patch = paste(output_dir, "patches", sep = "/")
Npatches <- as.numeric(readLines(con = file.path(output_dir,
                                                 "patches/Npatches.dat")))[[1]]
# OP_noL = matrix(NA, Npatches, length(par.wl) + length(nir.wl))
# for (ipatch in seq(1, Npatches)) {
#   OP_noL[ipatch, ] <- get.EDR.output(output.path_patch,
#                                      ipatch)
# }

LAI_sim=matrix(0.,Npatches)
CA_sim=matrix(0.,Npatches)
for (i in sequence(Npatches)){
  line <- readLines(con = file.path(output_dir,paste("patches/LAI",i,'.dat',sep='')))
  temp=as.numeric(strsplit(line,'\\s+')[[1]][-1])
  
  line <- readLines(con = file.path(output_dir,paste("patches/CA",i,'.dat',sep='')))
  temp_CA=as.numeric(strsplit(line,'\\s+')[[1]][-1])
  
  LAI_sim[i] <- sum(temp)
  CA_sim[i] <- mean(temp_CA)
}


file<-file.path(output_dir,paste('history','-S-',toString(year_s), '-',sprintf("%02d", month),'-01-120000-g01.h5',sep=''))

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


# sla=16.017942428
# CA_p= mapply(dbh2ca,DBH,HITE,sla,PFT)
# crown_area = pmin(1.0, nplant * CA_p)
# 
# dbh2ca <-function(dbh,hite,SLA,PFT){
#   
#   b1Ca<-1.1257275343
#   b2Ca<-1.0521197319
#   
#   dbh_crit=96.2577896118
#   ldbh = dbh_crit
#   
#   dbh2ca = b1Ca * min(dbh, ldbh ) ** b2Ca
#   return(dbh2ca)
# }

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
matplot(all_WL,t(OP_noL),col=Colors_noL,type='l',xlim=c(400,2500),ylim=c(0,0.75),lty=1)

Red=OP_noL[,300]
NIR=OP_noL[,500]
NDVI=(NIR-Red)/(NIR+Red)
plot.new()
par(mfrow=c(1,2))
plot(LAI_total_noL,Red,xlab='LAI')
plot(LAI_total_noL,NIR,xlab='LAI')


data<-data.frame(x=NDVI,y=LAI_total_noL)
exponential.model <- lm(log(y)~ x,data = data)
y_sort<-sort(as.vector(data$y),index.return=TRUE)
x<-y_sort$x
y<-as.vector(data$x[y_sort$ix])
#m <- nls(y~a*x/(b+x),start = list(a=1,b=1))
m <- nls(y~a/(1+exp(-b*x)),start = list(a=1,b=1))
m <- nls(x~a + b*exp(c*y) ,start = list(a=0.8,b=0.001,c=5))
residuals<-exp(exponential.model$residuals)

plot.new()
par(mfrow=c(1,1))
plot(x,y,type='p')
lines(sort(x),predict(m),lty=2,col="red",lwd=3)

var=nplant_pa
#var[28]<-var[44]<-NA
plot.new()
plot(var,residuals,xlab='Plant density',ylab='residuals')
lin.model <- lm(residuals ~ var)
lines(var[is.finite(var)],lin.model$fitted.values,col='black',lty=1)

summary(lin.model)

#################################################################################################
# Liana simulation
dir='/home/carya/output/PEcAn_99000000007/run/99000001833'
odir= '/home/carya/output/PEcAn_99000000007/out/99000001833'
output_dir='/home/carya/temp'

file<-file.path(output_dir,paste('history','-S-',toString(year_s), '-',sprintf("%02d", month),'-01-120000-g01.h5',sep=''))

ED2IN_file <- file.path(dir,"ED2IN")
ed2in <- read_ed2in(ED2IN_file)
ed2in$SFILOUT <- file.path(odir,'history')
ed2in$POI_LAT <- 9
ed2in$POI_LON <- -79
ed2in$ED_MET_DRIVER_DB <- file.path("/home/carya/raw_inputED/met/BCI","BCI_DRIVER_HEADER2")
ed2in$VEG_DATABASE <- file.path("/home/carya/raw_inputED/OGE2",basename(ed2in$VEG_DATABASE))
ed2in$SOIL_DATABASE <- file.path("/home/carya/raw_inputED/FAO",basename(ed2in$SOIL_DATABASE))
ed2in$THSUMS_DATABASE <- file.path("/home/carya/raw_inputED/thermal_sums",basename(ed2in$THSUMS_DATABASE))


history.path = dirname(ed2in$SFILOUT)
datetime = dtime

# The model
ed2in$INCLUDE_THESE_PFT <- c(2,3,4,17)
ed2in$MAXPATCH <- 0
ed2in$MAXCOHORT <- 0

edr_ed2in <- PEcAnRTM::setup_edr(ed2in, output_dir, dtime,patches=TRUE)

trait_values=list(Early = list(),
                  Mid = list(),
                  Late = list(),
                  Liana = list())

spectra_list <- list(Early = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Mid   = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Late  = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"),
                     Liana = PEcAnRTM::prospect(c(1.4, 40, 0.01, 0.010), version = "4"))

spectra_list$Early[,1]=Rall_T
spectra_list$Mid[,1]=Rall_T
spectra_list$Late[,1]=Rall_T
spectra_list$Liana[,1]=Rall_T

plot.new()
par(mfrow=c(1,2))
plot(all_WL,spectra_list$Early[,1],lty=1,col='red',type='l')
lines(all_WL,spectra_list$Mid[,1],lty=2,col='black',type='l')
lines(all_WL,spectra_list$Late[,1],lty=3,col='green',type='l')
lines(all_WL,spectra_list$Liana[,1],lty=3,col='blue',type='l')

plot(all_WL,spectra_list$Early[,2],lty=1,col='red',type='l')
lines(all_WL,spectra_list$Mid[,2],lty=2,col='black',type='l')
lines(all_WL,spectra_list$Late[,2],lty=3,col='green',type='l')
lines(all_WL,spectra_list$Liana[,2],lty=3,col='blue',type='l')

OP <- PEcAnRTM::EDR(
    img_path = img_path,
    ed2in_path = edr_ed2in,
    spectra_list = spectra_list,
    trait.values = trait_values,
    edr_exe_path =  edr_exe_path,
    patches = TRUE
    )

hfile <- hdf5r::H5File$new(file)

PACO_N=readDataSet(hfile[['PACO_N']])
DBH=readDataSet(hfile[['DBH']])
HITE=readDataSet(hfile[['HITE']])
PFT=readDataSet(hfile[['PFT']])
LAI=readDataSet(hfile[['LAI_CO']])
AREA=readDataSet(hfile[['AREA']])
BLEAF=readDataSet(hfile[['BLEAF']])
NPATCHES_GLOBAL=readDataSet(hfile[['NPATCHES_GLOBAL']])

hfile$close_all()

ipaconow  = rep(sequence(NPATCHES_GLOBAL),times=PACO_N)

LAI_L=LAI
LAI_L[PFT!=17]=0
LAI_L[DBH<2]=0

LAI_total = tapply( X = LAI, INDEX = ipaconow, FUN   = sum)#end tapply
LAI_liana = tapply( X = LAI_L, INDEX = ipaconow, FUN   = sum)#end tapply
Liana_contribution = LAI_liana/LAI_total

order=sort(LAI_total, index.return=TRUE)$ix
order_LianaLAI=sort(LAI_liana, index.return=TRUE)$ix

rbPal <- colorRampPalette(c('grey','darkgreen'))
C <- rbPal(NPATCHES_GLOBAL)
Colors <- Colors_liana <- C
Colors[order] <- C
Colors_liana[order_LianaLAI] <-C

plot.new()
par(mfrow=c(1,1))
matplot(all_WL,t(OP), type = 'l',lty=1,col = Colors,xlab='Wavelength',ylab='Reflectance') #plot

#################################################################################################
# Liana simulation, Liana properties

spectra_list$Liana[,1]=Rall_L

plot.new()
par(mfrow=c(1,2))
plot(all_WL,spectra_list$Early[,1],lty=1,col='red',type='l')
lines(all_WL,spectra_list$Liana[,1],lty=3,col='blue',type='l')
lines(all_WL,spectra_list$Mid[,1],lty=2,col='black',type='l')
lines(all_WL,spectra_list$Late[,1],lty=3,col='green',type='l')

plot(all_WL,spectra_list$Early[,2],lty=1,col='red',type='l')
lines(all_WL,spectra_list$Mid[,2],lty=2,col='black',type='l')
lines(all_WL,spectra_list$Late[,2],lty=3,col='green',type='l')
lines(all_WL,spectra_list$Liana[,2],lty=3,col='blue',type='l')

OP_Liana <- PEcAnRTM::EDR(
  img_path = img_path,
  ed2in_path = edr_ed2in,
  spectra_list = spectra_list,
  trait.values = trait_values,
  edr_exe_path =  edr_exe_path,
  patches = TRUE
  )

plot.new()
par(mfrow=c(1,2))
matplot(all_WL,t(OP_Liana), type = 'l',lty=1,col = Colors,xlab='Wavelength',ylab='Reflectance',ylim=c(0,0.7)) #plot
title('Including Liana properties')
legend('topright',legend=c('Increasing patch LAI'),col='darkgreen',lty=1)

Diff_spectrum<-(t(OP_Liana)-t(OP))
matplot(all_WL,Diff_spectrum, type = 'l',lty=1,col = Colors_liana,xlab='Wavelength',ylab='Reflectance',ylim=c(-0.1,0.1)) #plot
title('Liana param. - No liana param.')
abline(h=0,lty=2,col='black')
legend('topright',legend=c('Increasing patch Liana LAI'),col='darkgreen',lty=1)

plot.new()
par(mfrow=c(1,1))
plot(LAI_liana,Diff_spectrum[300,],type='p',col='red',ylim=c(-0.8,0.1)*0.05)
lines(LAI_liana,Diff_spectrum[100,],type='p',col='blue')
lines(LAI_liana,Diff_spectrum[200,],type='p',col='green')
lines(LAI_liana,Diff_spectrum[600,],type='p',col='black')
df=data.frame(x=LAI_liana,y=Diff_spectrum[600,])

LM <- lm(y~x,data = df)

mat_L=repmat(LAI_total,1,NPATCHES_GLOBAL_noL)
mat_noL=t(repmat(LAI_total_noL,1,NPATCHES_GLOBAL))
diff_sq=(mat_L-mat_noL)^2
pos_min=which(diff_sq==min(diff_sq))
pos_noL=which(LAI_total_noL==mat_noL[pos_min])
pos_L=which(LAI_total==mat_L[pos_min])

LAIm=sum(LAI_total*AREA)
LAIm_noL=sum(LAI_total_noL*AREA_noL)
simL=which(abs(LAI_total-LAIm)==min(abs(LAI_total-LAIm)))
simnoL=which(abs(LAI_total_noL-LAIm)==min(abs(LAI_total_noL-LAIm)))

plot.new()
par(mfrow=c(1,2))
plot(all_WL,OP_Liana[simL,], type = 'l',lty=1,col = 'blue',lwd=2,ylim=c(0,0.7))
lines(all_WL,OP_noL[simnoL,], type = 'l',lty=1,col = 'red',lwd=2)
title('Mean LAI')
legend('topright',legend=round(c(LAI_total[simL],LAI_total_noL[simnoL]),digits = 2),col=c('red','blue'),lty=1)

plot(all_WL,OP_Liana[pos_L,], type = 'l',lty=1,col = 'blue',lwd=2,ylim=c(0,0.7))
lines(all_WL,OP_noL[pos_noL,], type = 'l',lty=1,col = 'red',lwd=2)
title('Similar LAI')
legend('topright',legend=round(c(LAI_total[pos_L],LAI_total_noL[pos_noL]),digits = 2),col=c('red','blue'),lty=1)

plot.new()
par(mfrow=c(1,2))
matplot(all_WL,t(OP_Liana), type = 'l',lty=1,col = Colors,xlab='Wavelength',ylab='Reflectance',ylim=c(0,0.8)) #plot
title('Liana')

matplot(all_WL,t(OP_noL), type = 'l',lty=1,col = Colors_noL,xlab='Wavelength',ylab='Reflectance',ylim=c(0,0.8)) #plot
title('No liana')

