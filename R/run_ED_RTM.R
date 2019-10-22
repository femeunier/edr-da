#' Run ED RTM
#'
#' @export

run_ED_RTM <- function(rundir,outdir,params,crown_mod,inventory,par.wl,nir.wl){
  
 params_EDR <- params[-1]  
 
 dbh <- inventory$dbh
 pft <- inventory$pft
 hite <- inventory$hite
 nplant <- inventory$nplant
 Npatch <- inventory$Npatch
 Ncohort <- inventory$Ncohort
 PAN <- inventory$PAN
 PACO_N <- inventory$PACO_N
 PFTselect <- inventory$PFTselect

  soil_moisture <- params_EDR[1]
  
  pft_params <- matrix(params_EDR[-(1)], ncol = npft)
  
  # Calculate allometries
  b1leaf <- pft_params[1,]
  b2leaf <- pft_params[2,]
  sla <-  pft_params[3,]
  
  orient_factor <- pft_params[4,]
  clumping_factor <- pft_params[5,]
  
  # Extract Prospect parameters
  N <- pft_params[6,]
  Cab <- pft_params[7,]
  
  Car <- pft_params[8,]
  Cw <- pft_params[9,]
  Cm <-  1/pft_params[3,]/10
  
  if (crown_mod == 0){
    cai <- rep(1,Ncohort)
    b1Ca <- NULL
    b2Ca <- NULL
  } else {
    b1Ca <- pft_params[10,]
    b2Ca <- pft_params[11,]  
  }
  
  soil_reflect <- hapke_soil(rep(soil_moisture,2101))
  soil_file <- file.path("/tmp","soil_reflect_par.dat")
  writeLines(text=paste((soil_reflect),collapse="\t"),con = soil_file)
  
  # Prospect parameters
  spectra_list <- list()
  for (ipft in seq(npft)){
    current_pft <- as.character(df_PFT$names[ipft])
    spectra_list[[current_pft]] <-  PEcAnRTM::prospect(c(N[ipft],Cab[ipft],Car[ipft],Cw[ipft],Cm[ipft]), version = "5")   
    spectra_list[[current_pft]][,] <- na.approx(spectra_list[[current_pft]][,],rule =2)
  }
  
  # Rest of config
  config_file <- file.path(rundir,"config.xml")
  xml <- XML::xmlParse(config_file)
  config <- XML::xmlToList(xml)
  
  config_temp <- config[["pft"]]
  ipft <- 1
  trait_values <- list()
  
  while (!is.null(config_temp)){
    pft_num <- config_temp$num
    pft_name <- as.character(df_PFT$names[df_PFT$PFTnum==pft_num])
    
    params2remove <- c("chlorophyll_a","chlorophyll_b","carotenoids","Cw")
    for (iparam in seq(params2remove)){
      config_temp[[params2remove[iparam]]] <- NULL
    }
    
    config_temp$b1Bl_large<-b1leaf[ipft]
    config_temp$b2Bl_large<-b2leaf[ipft]
    
    config_temp$b1Bl_small<-b1leaf[ipft]
    config_temp$b2Bl_small<-b2leaf[ipft]
    
    if (!(is.null(b1Ca) & !is.null(b1Ca))){
      config_temp$b1Ca<-b1Ca[ipft]
      config_temp$b2Ca<-b2Ca[ipft]
    }
    
    config_temp$orient_factor<-orient_factor[ipft]
    config_temp$clumping_factor<-clumping_factor[ipft]
    
    if ("SLA" %in% names(config_temp)){
      config_temp$SLA <- sla[ipft]*0.48
    }
    if ("Vcmax" %in% names(config_temp)){
      config_temp$Vcmax <- as.character(as.numeric(config_temp$Vcmax)*2.4) # To go back to PEcAn values
    }
    
    trait_values[[pft_name]] <- config_temp
    
    config[["pft"]] <- NULL
    config_temp <- config[["pft"]]
    ipft <- ipft + 1
  }
  
  # create  temp directory to run
  # temp_dir <- outdir
  temp_dir <- system2("mktemp",list("-d","-p",outdir),stdout=TRUE)
  system2("mkdir",file.path(temp_dir,"patches"),stdout=NULL)
  dummy <- file.copy(file.path(outdir,"ED2IN"),temp_dir)
  dummy <- file.copy(list.files(outdir,pattern = "*.h5",full.names = TRUE),temp_dir)
  
  write.table(t(params),file = file.path(temp_dir,"params.txt"),sep = ' ',row.names = FALSE,col.names = FALSE)
  
  output_RTM <- tryCatch(
    PEcAnRTM::EDR(img_path = NULL,
                  ed2in_path = file.path(temp_dir,"ED2IN"),
                  spectra_list = spectra_list,
                  trait.values = trait_values,
                  soil_reflect_path = soil_file,
                  edr_exe_path =  edr_exe_path,
                  par.wl = par.wl, 
                  nir.wl = nir.wl,
                  patches = TRUE),
    error = function(e) NULL)
  
  if (output_RTM) return(list(output_RTM = NULL,COI = NULL))
  # read LAI
  lai <- rep(NA,Ncohort)
  PACO_ID <- cumsum(c(1,PACO_N))
  for (ipa in seq(Npatch)){
    lai[seq(PACO_ID[ipa],PACO_ID[ipa+1]-1)]<-scan(file.path(file.path(temp_dir,"patches"),paste0("LAI",ipa,".dat")),quiet=TRUE)
  }

  if (crown_mod !=0){
    cai <- cai_allometry(dbh,nplant,b1leaf[pft],b2leaf[pft],sla[pft],b1Ca[pft],b2Ca[pft])
  } else {
    cai <- 1+0*lai
  }
  
  # Compute Crown occupancy index
  COI <-rep(NA,Npatch)
  
  for (ipa in seq(1,Npatch)){
    
    pos <- (PAN == ipa)
    
    hiteconow   <- hite[pos]
    dbhconow    <- dbh[pos]
    pftconow    <- pft[pos]
    laiconow    <- lai[pos]
    nplantconow <- nplant[pos]
    CAconow <- cai[pos]
    
    # Compute COI
    LAItot <- sum(laiconow)
    LAI_alpha <- (1-alpha_frac)*LAItot
    h_alpha <- interp1(x = cumsum(laiconow), y = hiteconow, xi = ifelse(LAI_alpha<laiconow[1],laiconow[1],LAI_alpha))
    
    COI[ipa] <- max(0,min(1,sum(laiconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/LAI_alpha))
  }
  # To check
  # system2("cp",list("-r",temp_dir,"~/data/RTM/"))
  
  
  # remove temporary
  system2("rm",list("-rf",temp_dir),stdout=NULL)
  
  
  return(list(output_RTM = output_RTM,COI = COI))

  
}