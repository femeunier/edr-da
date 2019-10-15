# rm(list = ls())


# Paper Contrasting leaf chemical traits in tropical lianas and trees:
# implications for future forest composition
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1461-0248.2012.01821.x
doi <- '10.1111/j.1461-0248.2012.01821.x'

# Table S2
Sites <- list(places = c('Africa','Australasia','Neotropics'),
              Nlianas = c(14+19+28,1+5+1+1+2+1+18+2,64+73+60+5+4+12+8+15+3+64+39+6+3+2+6+18+3+13+30+20+27+61+5+6+6+4+6+1+12+24+15+5+6+8+10+16+16),
              Ntrees = c(135+134+297,22+19+22+20+10+14+378+18,332+335 +65 +29+50+205+175+120+44+302+333+429+81+61+47+215+10+119+183+49+150+205+44+57+142+127+52+54+334+332+192+60+71+20+85+206+104))

data <- data.frame(LMA = c(117.23,	155.18,103.87,120.04,90.21,104.35),
                   LMA_sd = c(36.01,48.44,37.4,40.24,35.87,34.61),
                   Ca = c(3.97,3.28,4.93,4.26,5.97,5.04),
                   Ca_sd = c(1.31,1.23,1.79,1.62,2.45,1.83),
                   Cb = c(1.5,1.22,1.84,1.6,2.26,1.89),
                   Cb_sd = c(0.54,0.49,0.68,0.65,0.97,0.73),
                   Car = c(1.24,1.02,1.5,1.29,1.71,1.48),
                   Car_sd = c(0.32,0.35,0.5,0.44,0.65,0.5),
                   Water = c(62.61,58.76,60.94,56.04,62.38,58.87)/100,
                   Water_sd = c(6.9,6.77,7.65,6.77,9.05,7.63)/100,
                   N = c(Sites$Nlianas[1],Sites$Ntrees[1],Sites$Nlianas[2],Sites$Ntrees[2],Sites$Nlianas[3],Sites$Ntrees[3]),
                   Site = c('Madagascar','Madagascar','Australasia','Australasia','Neotropics','Neotropics'),
                   GF = c('Liana','Tree','Liana','Tree','Liana','Tree'))

# units (LMA (g/m2), Ca (mg/g), Cb (mg/g), Car (mg/g), Water (%)) --> paper
# EWT = (FM - DM)/(rho*A) = theta/(1-theta)*(SLA/rho)
# units prospect (Cab (µg/cm2), Car (µg/cm2), Cw (cm), Cm (g/cm2)))

Ca <- Cb <- Car <- Cw <- Cm <- Cab <- 
  Ca_sd <- Cb_sd <- Car_sd <- Cw_sd <- Cm_sd <- Cab_sd <- rep(NA,nrow(data))

Nsample = 10000
for (i in seq(nrow(data))){
  Ca_sample <- with(data,rnorm(n = N[i]*Nsample,mean = Ca[i], sd = Ca_sd[i]))
  Cb_sample <- with(data,rnorm(n = N[i]*Nsample,mean = Cb[i], sd = Cb_sd[i]))
  Car_sample <- with(data,rnorm(n = N[i]*Nsample,mean = Car[i], sd = Car_sd[i]))
  Water_sample <- with(data,rnorm(n = N[i]*Nsample,mean = Water[i], sd = Water_sd[i]))
  LMA_sample <- with(data,rnorm(n = N[i]*Nsample,mean = LMA[i], sd = LMA_sd[i]))
  
  pos <- which(Ca_sample      > 0.01*mean(Ca_sample) & 
                 Cb_sample    > 0.01*mean(Cb_sample) & 
                 Car_sample   > 0.01*mean(Car_sample) &
                 Water_sample > 0.01*mean(Water_sample) & Water_sample < 0.99*mean(Water_sample) & 
                 LMA_sample   > 0.01*mean(LMA_sample))
  
  Ca_sample <- Ca_sample[pos]
  Cb_sample <- Cb_sample[pos]
  Car_sample <- Car_sample[pos]
  Water_sample <- Water_sample[pos]
  LMA_sample <-  LMA_sample[pos]
  
  # mean
  Ca[i] <- mean(Ca_sample*LMA_sample/10)
  Cb[i] <- mean(Cb_sample*LMA_sample/10)
  Cab[i] <- Ca[i] + Cb[i]
  Car[i] <- mean(Car_sample*LMA_sample/10)
  Cw[i] <- mean(Water_sample/(1 - Water_sample)*LMA_sample/10000/1)
  Cm[i] <- mean(LMA_sample/10000)
  
  # sd
  Ca_sd[i] <- sd(Ca_sample*LMA_sample/10)
  Cb_sd[i] <- sd(Cb_sample*LMA_sample/10)
  Cab_sd[i] <- sd(Ca_sample + Cb_sample)
  Car_sd[i] <- sd(Car_sample*LMA_sample/10)
  Cw_sd[i] <- sd(Water_sample/(1 - Water_sample)*LMA_sample/10000/1)
  Cm_sd[i] <- sd(LMA_sample/10000)
}

data_units <- data.frame(Ca = Ca, 
                         Cb = Cb,
                         Cab = Cab,
                         Car = Car,
                         Cw = Cw,
                         Cm = Cm,
                         Ca_sd = Ca_sd, 
                         Cb_sd = Cb_sd,
                         Cab_sd = Cab_sd,
                         Car_sd = Car_sd,
                         Cw_sd = Cw_sd,
                         Cm_sd = Cm_sd,
                         N = c(Sites$Nlianas[1],Sites$Ntrees[1],Sites$Nlianas[2],Sites$Ntrees[2],Sites$Nlianas[3],Sites$Ntrees[3]),
                         Site = c('Madagascar','Madagascar','Australasia','Australasia','Neotropics','Neotropics'),
                         GF = c('Liana','Tree','Liana','Tree','Liana','Tree'))

data_units