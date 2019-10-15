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

