# Kalacska et al.

rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure1_kalacska.csv"
data <- read.csv(data2read)
PFTs <- c("Tree_optical","Liana_optical")
C <- c("#137300","#1E64C8")
plot(NA,NA,xlim=c(400,2500),ylim=c(0,0.5),xlab="Wavelength",ylab="Reflectance")
data_R <- list()
df <- data.frame()
for (ipft in seq(PFTs)){
  pft <- PFTs[ipft]
  columns <- ((ipft-1)*2+1):((ipft-1)*2+2)
  temp <- data[!is.na(data[,columns[1]]),c(columns)]
  pos <- sort(temp[,1],index.return = TRUE)
  data_R[[pft]] <- temp[pos$ix,]
  lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
  names(data_R[[pft]]) <- c("wavelength","reflectance")
  
  df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                            Reflectance = data_R[[pft]] [,2],
                            pft = pft))
}

saveRDS(df,file = "~/data/RTM/Figure1_kalacska.rds")

##################################################################################
# Castro et al.

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
                                        Reflectance = temp[pos$ix,2],
                                        pft = pft))
}

saveRDS(data_raw,file = "~/data/RTM/Figure6_castro.rds")

##################################################################################
# de Guzman et al.

rm(list=ls())

data2read <- "/home/carya/data/RTM/data/Figure1_Guzman.csv"
data <- read.csv(data2read)
WL <- data[,1]
Reflectance <- data[,-1]

PFTs <- c("Tree_optical","Liana_optical")
patterns <- c("Tree*","Liana*")

C <- c("#137300","#1E64C8")

plot(NA,NA,xlim=c(400,950),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

data_R <- data.frame()
for (ipft in seq(PFTs)){
  
  pft <- PFTs[ipft]
  is_pft <- grepl(patterns[ipft],colnames(Reflectance))
  
  columns <- which(is_pft)
  temp <- melt((Reflectance[,columns])) %>% rename(Reflectance = value) %>% select(Reflectance)
  matlines(WL,data[,columns],col = C[ipft])
  
  data_R <- rbind(data_R,data.frame(wavelength = rep(WL,length(columns)),
                                    Reflectance = temp,
                                    pft = pft))
  
}

saveRDS(data_R,file = "~/data/RTM/Figure1_Guzman.rds")

##################################################################################
# Castro et al. (4 and 5)

rm(list=ls())

sites <- c("FTS","PNM")
data2read <- "/home/carya/data/RTM/data/Figures4and5_castro.csv"
data <- read.csv(data2read)
PFTs <- c("Liana_optical","Tree_optical")
C <- c("#1E64C8","#137300")
plot(NA,NA,xlim=c(400,1000),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

for (isite in seq(sites)){
  data_R <- list()
  df <- data.frame()
  for (ipft in seq(PFTs)){
    pft <- PFTs[ipft]
    columns <- (isite -1)*4 + (((ipft-1)*2+1):((ipft-1)*2+2))
    temp <- data[!is.na(data[,columns[1]]),c(columns)]
    pos <- sort(temp[,1],index.return = TRUE)
    data_R[[pft]] <- temp[pos$ix,]
    lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
    names(data_R[[pft]]) <- c("wavelength","reflectance")
    
    df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                              Reflectance = data_R[[pft]] [,2],
                              pft = pft))
  }
  saveRDS(df,file = paste0("~/data/RTM/Figures4and5_castro_",sites[isite],".rds"))
}

##################################################################################
# Castro et al. (4 and 5)

rm(list=ls())

sites <- c("FTS","PNM")
data2read <- "/home/carya/data/RTM/data/Figure6_sanchez2009.csv"
data <- read.csv(data2read)
PFTs <- c("Tree_optical","Liana_optical")
C <- c("#137300","#1E64C8")
plot(NA,NA,xlim=c(400,1000),ylim=c(0,0.6),xlab="Wavelength",ylab="Reflectance")

for (isite in seq(sites)){
  data_R <- list()
  df <- data.frame()
  for (ipft in seq(PFTs)){
    pft <- PFTs[ipft]
    columns <- (isite -1)*4 + (((ipft-1)*2+1):((ipft-1)*2+2))
    temp <- data[!is.na(data[,columns[1]]),c(columns)]
    pos <- sort(temp[,1],index.return = TRUE)
    data_R[[pft]] <- temp[pos$ix,]
    lines(data_R[[pft]] [,1],data_R[[pft]] [,2],col = C[ipft])
    names(data_R[[pft]]) <- c("wavelength","reflectance")
    
    df <- rbind(df,data.frame(wavelength =  data_R[[pft]] [,1],
                              Reflectance = data_R[[pft]] [,2],
                              pft = pft))
  }
  saveRDS(df,file = paste0("~/data/RTM/Figure6_sanchez2009_",sites[isite],".rds"))
}

