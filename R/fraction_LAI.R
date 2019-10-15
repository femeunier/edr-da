#' Relative distribution of LAI for a specific PFT
#'
#' @export

fraction_LAI <- function(h5file,ispatch = TRUE,PFTselect = 17,alpha_frac = 0.5){
  
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
  CA =readDataSet(hfile[['CROWN_AREA_CO']])
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
  
  init <- ifelse(ispatch,1,patch_number)
  
  fraction_LAI <- COI <- COI2 <- ifelse(ispatch,rep(NA,patch_number),c())
  
  for (ipa in seq(init,patch_number)){
    
    pos <- (PAID == ipa)
    
    hiteconow   <- height[pos]
    dbhconow    <- DBH[pos]
    pftconow    <- PFT[pos]
    laiconow    <- LAI[pos]
    nplantconow <- nplant[pos]
    CAconow <- CA[pos]
    Nco <- length(hiteconow)
    
    if (ispatch) {
      fraction_LAI[ipa] <- sum(laiconow[pftconow %in% PFTselect])/sum(laiconow)
      
      # Compute COI
      LAItot <- sum(laiconow)
      LAI_alpha <- (1-alpha_frac)*LAItot
      h_alpha <- interp1(x = as.vector(cumsum(laiconow)), 
                         y = as.vector(hiteconow), 
                         xi = ifelse(LAI_alpha<laiconow[1],laiconow[1],LAI_alpha))
      
      COI[ipa] <- sum(laiconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/LAI_alpha
      COI2[ipa] <- sum(CAconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/sum(CAconow[hiteconow >= h_alpha])
        
    } else{
      fraction_LAI <- sum(laiconow[pftconow %in% PFTselect])/sum(laiconow) 
    }
    


  }
  writeLines(as.character(fraction_LAI),con = file.path(RTM_path,"LAI_fraction.dat"))
  writeLines(as.character(COI),con = file.path(RTM_path,"COI.dat"))
  writeLines(as.character(COI2),con = file.path(RTM_path,"COI2.dat"))
}