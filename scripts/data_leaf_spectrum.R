rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure1_kalacska.csv"
data <- read.csv(data2read)
PFTs <- c("Liana_optical","Tree_optical")
C <- c("#1E64C8","#137300")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
for (ipft in seq(PFTs)){
  pft <- PFTs[ipft]
  columns <- ((ipft-1)*2+1):((ipft-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[pft]] <- temp[pos$ix,]
  lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
  names(data_R[[pft]]) <- c("wavelength","reflectance")
}

saveRDS(data_R,file = "~/data/RTM/Spectrum_liana_data.R")

