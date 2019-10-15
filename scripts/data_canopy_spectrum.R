rm(list=ls())

data2read <- "~/data/RTM/data/figure3_marvin.csv"# rm(list=ls())

data <- read.csv(data2read)
scenarios <- c("0.75","0.25")
C <- c("#137300","#1E64C8")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
}

saveRDS(data_R,file = "~/data/RTM/Spectrum_canopy_data.R")



##################################################################

rm(list=ls())

data2read <- "~/data/RTM/data/Figure1_kalacslab.csv"

data <- read.csv(data2read)
scenarios <- c("0.4","0.")
C <- c("#1E64C8","#137300")

plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance",xaxt="n")
axis(side=1, at=c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244), labels = c(450,518,614,757,950,1172,1390,1593,1777,1946,2101,2244))

data_R <- list()

for (iscenario in seq(scenarios)){
  scenario <- scenarios[iscenario]
  columns <- ((iscenario-1)*2+1):((iscenario-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[scenario]] <- temp[pos$ix,]
  names(data_R[[scenario]]) <- c("wavelength","reflectance")
  
  th <- 1390
  fac <- 1.3
  data_R[[scenario]]["wavelength"][data_R[[scenario]]["wavelength"]<th ] <-  data_R[[scenario]]["wavelength"][data_R[[scenario]]["wavelength"]< th]/fac
  
  th1 <- 757
  th2 <- 1390
  fac <- 1.2
  pos <- which(data_R[[scenario]]["wavelength"]> th1 &
    data_R[[scenario]]["wavelength"]<th2)
  data_R[[scenario]][pos,"wavelength"] <-  
    data_R[[scenario]][pos,"wavelength"]*seq(1,fac,length.out = length(pos))
  
  th1 <- 0
  th2 <- 600
  fac <- 1.25
  pos <- which(data_R[[scenario]]["wavelength"]> th1 &
                 data_R[[scenario]]["wavelength"]<th2)
  data_R[[scenario]][pos,"wavelength"] <-  
    data_R[[scenario]][pos,"wavelength"]*seq(fac,1,length.out = length(pos))

  lines(data_R[[scenario]] [,1],data_R[[scenario]] [,2],col = C[iscenario],type='p',pch=19,lwd=0.01)

}

saveRDS(data_R,file = "~/data/RTM/Spectrum_canopy_data.R")

