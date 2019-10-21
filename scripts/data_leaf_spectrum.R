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

##################################################################################

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
                                        reflectance = temp[pos$ix,2],
                                        pft = pft))
}

saveRDS(data_R,file = "~/data/RTM/Figure6_castro.rds")
saveRDS(data_raw,file = "~/data/RTM/Figure6_castro_raw.rds")
