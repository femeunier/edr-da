
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

