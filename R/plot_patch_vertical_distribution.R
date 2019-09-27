#' Plot patch vertical distribution
#'
#' @export

plot_patch_vertical_distribution <- function(h5file,ispatch = TRUE){

hfile <- hdf5r::H5File$new(h5file)
RTM_path <-file.path(dirname(h5file))

if (!dir.exists(RTM_path)){
  dir.create(RTM_path)
}

if (ispatch){
  figure_path <- file.path(RTM_path,"patches")
} else{
  figure_path <- RTM_path
}

PACO_N=(readDataSet(hfile[['PACO_N']]))
patch_number <- length(PACO_N)
cohort_number <- sum(PACO_N)
LAI=readDataSet(hfile[['LAI_CO']])
height=readDataSet(hfile[['HITE']])
PFT=readDataSet(hfile[['PFT']])
DBH=readDataSet(hfile[['DBH']])
patch_area=readDataSet(hfile[['AREA']])
NPATCHES=readDataSet(hfile[['NPATCHES_GLOBAL']])
nplant=readDataSet(hfile[['NPLANT']])
CoArea <- rep(patch_area,PACO_N)
npft <- max(unique(PFT))
PAID <- rep(1:patch_number,PACO_N)

hfile$close_all()

Col_pft <- data.frame(col = c('darkgreen','darkblue'), num=c(3,17), names = c('Tree','Liana'))

init <- ifelse(ispatch,1,patch_number)
for (ipa in seq(init,patch_number)){
  
  pos <- (PAID == ipa)
  
  hiteconow   <- height[pos]
  dbhconow    <- DBH[pos]
  pftconow    <- PFT[pos]
  laiconow    <- LAI[pos]
  nplantconow <- nplant[pos]
  Nco <- length(hiteconow)
  
  pft_uni <- unique(pftconow)
  
  LAIs <- laiconow
  cumLAI <- cumsum(LAIs)
  hites <- hiteconow
  colors <- as.character(Col_pft$col[match(pftconow,Col_pft$num)])
  
  
  png(filename = file.path(figure_path,paste0("patch_",ipa,'.png')),
      width = 15, height = 10, units = "cm", pointsize = 12,
      bg = "white",  res = 300)
  
  
  layout(mat=rbind(2,1),heights=c(4.6,1.4))
  #---------------------------------------------------------------------------#
  
  #------ Legend. ------------------------------------------------------------#
  par(mar=c(0.1,4.6,.001,2.1))
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  legend( x      = "bottom"
          , inset  = 0.0
          , legend = c(as.character(Col_pft$names),'Total')
          , col   = c(as.character(Col_pft$col),'black')
          , lty   = 1
          , lwd   = 2 
          , ncol   = 3
          , title  = expression(bold("Plant functional type"))
          , cex    = 1
          , bg     = 0
          , xpd    = TRUE
          , box.lty= 0
  )#end legend
  #---------------------------------------------------------------------------#
  
  
  #----- Plot all monthly means together. ------------------------------------#
  par(mar=c(3.5,4.6,2.1,2.1), mgp=c(2.5,1,0))
  plot.window(xlim=c(0,1),ylim=c(0,1))
  
  
  letitre=paste("Patch area =",signif(patch_area[ipa],digits=2))
  
  where = pretty(cumLAI)
  
  xpos <- sample(seq(where[1],where[length(where)],length.out = Nco))
  
  plot(-1,0,lwd=0,ylim=c(0,1.02*max(hites)),xlim=c(0,where[length(where)]),ylab = 'Height [m]',xlab="Cumulate LAI [-]",main=letitre)
  for (ico in seq(1,Nco)){
    lines(c(0,where[length(where)]),c(hites[ico],hites[ico]),col=colors[ico],lwd=max(0.5,(LAIs[ico]/1))*2)
    lines(c(xpos[ico],xpos[ico]),c(0,hites[ico]),col='chocolate4',lwd=1)
  }
  
  
  for (ipft in c(pft_uni,npft+1)){
    
    if (ipft < npft + 1){
      hites_temp <- hiteconow[pftconow == ipft]
      LAI_temp <- LAIs[pftconow == ipft]    
      cumLAI <- cumsum(LAI_temp)
      C <- as.character(Col_pft$col[match(ipft,Col_pft$num)])
      lines(c(cumLAI,cumLAI[length(cumLAI)]),c(hites_temp,0),col=C,lwd=2,type='l')
    } else {
      hites_temp <- hiteconow
      LAI_temp <- LAIs
      cumLAI <- cumsum(LAI_temp)
      lines(c(cumLAI,cumLAI[length(cumLAI)]),c(hites_temp,0),col='black',lwd=2,type='l')
    }
  }
  dev.off()
}

}