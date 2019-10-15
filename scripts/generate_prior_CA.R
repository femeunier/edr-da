rm(list=ls())

b1Ca <- 1.1257275343 # m²
b2Ca <- 1.0521197319 # unitless

dbhs <- seq(0.0001,100,length.out = 100)

CA <- b1Ca*(dbhs**b2Ca)
CR <- sqrt(CA/pi)
plot(dbhs,CR,type='l',ylim=c(0,10))

N <- 10000
b1Ca_dist <-data.frame(distn = "lnorm",parama = log(b1Ca),paramb = b1Ca/10)
b2Ca_dist <-data.frame(distn = "lnorm", parama = log(b2Ca),paramb = b2Ca/3)

distn.stats("lnorm", b1Ca,0.1)

b1Ca_all <-get.sample(b1Ca_dist, n = N)
b2Ca_all <-get.sample(b2Ca_dist, n = N)
b1CR <- sqrt(b1Ca_all/pi)

CR_Sample_all <- matrix(NA,N,length(dbhs))
for (i in seq(N)){
  CA <- b1Ca_all[i]*(dbhs**b1Ca_all[i])
  CR_sample <- sqrt(CA/pi)
  CR_Sample_all[i,] <- CR_sample
  lines(dbhs,CR_sample,type='l',col='grey')
}

lines(dbhs,CR,type='l',lwd=2)
max(CR_Sample_all)