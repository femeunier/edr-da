#' Calculate leaf biomass given DBH
#'
#' Allometry equation is of the form `Y = b1 * DBH ^ b2`.
#'
#' Original function in ED src is in `utils/allometry.f90`. See also
#' the `area_indices` subroutine in the same file.
#'
#' LAI is calculated as `bleaf * nplant * SLA`.
#'
#' @param dbh Cohort diameter at breast height
#' @param b1 Allometry equation base
#' @param b2 Allometry equation exponent
#' @return Leaf biomass
#' @author Alexey Shiklomanov
#' @export
size2bl <- function(dbh, b1, b2) {
  # Carbon to biomass ratio of plant tissues. Defined in
  # `memory/pft_coms`; initialized in `init/ed_params`.
  C2B <- 2.0
  b1 / C2B * dbh ^ b2
}

#' Wood allometry
#' 
#' Calculate wood area index (WAI) given DBH
#'
#' Two possibilities here, depending on `iallom`.
#' The default value is: `WAI = nplant * b1 * DBH ^ b2
#'
#' An alternative formulation is that it is always 11% of LAI: `WAI =
#' 0.11 * LAI`. But we don't use that here because it's trivial.
#' @param dbh Cohort diameter at breast height
#' @param nplant Stem density (stems m-2)
#' @param b1 Allometry equation base
#' @param b2 Allometry equation exponent
#' @return Wood area index
#' @author Alexey Shiklomanov
#' @export
wai_allometry <- function(dbh, nplant, b1, b2) {
  nplant * b1 * dbh ^ b2 
}

#' crown area allometry
#'
#' @export
#' 
dbh2ca <- function(dbh, b1, b2) {
  b1 * dbh ^ b2
}
#' cai allometry
#'
#' @export
cai_allometry <- function(dbh,nplant,b1Bl,b2Bl,sla,b1Ca,b2Ca) {
  Ca <- dbh2ca(dbh, b1Ca, b2Ca)
  Bl <- size2bl(dbh, b1Bl, b2Bl)
  loclai <- sla *Bl
  dbh2ca <- pmin(loclai, Ca)
  pmin(1.0, nplant * dbh2ca)
}

#' COI from h5 file
#'
#' @export

calc_COI <- function(h5file,PFTselect) {
  
  hfile <- hdf5r::H5File$new(h5file)
  
  dbh <- readDataSet(hfile[["DBH"]])
  nplant <- readDataSet(hfile[["NPLANT"]])
  Npatch <- readDataSet(hfile[["NPATCHES_GLOBAL"]])
  hite <- readDataSet(hfile[["HITE"]])
  lai <- readDataSet(hfile[["LAI_CO"]])
  cai <- readDataSet(hfile[["CROWN_AREA_CO"]])
  pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
  PACO_N <- readDataSet(hfile[["PACO_N"]])
  PAN <- rep(1:Npatch,PACO_N)
  
  hfile$close_all()
  
  COI <-rep(NA,Npatch)
  
  for (ipa in seq(1,Npatch)){
    
    pos <- (PAN == ipa)
    
    hiteconow   <- hite[pos]
    dbhconow    <- dbh[pos]
    pftconow    <- pft[pos]
    laiconow    <- lai[pos]
    nplantconow <- nplant[pos]
    CAconow <- cai[pos]
    
    COI[ipa] <- ifelse(any(PFTselect %in% pftconow),1,0)
  }

  return(COI)
}


#' Liana indexes from h5 file
#'
#' @export

calc_liana_index <- function(h5file,PFTselect,alpha_frac=0.8) {
  
  hfile <- hdf5r::H5File$new(h5file)
  
  dbh <- readDataSet(hfile[["DBH"]])
  nplant <- readDataSet(hfile[["NPLANT"]])
  Npatch <- readDataSet(hfile[["NPATCHES_GLOBAL"]])
  hite <- readDataSet(hfile[["HITE"]])
  lai <- readDataSet(hfile[["LAI_CO"]])
  cai <- readDataSet(hfile[["CROWN_AREA_CO"]])
  pft <- match(readDataSet(hfile[["PFT"]]),df_PFT$PFTnum) # Liana = 1, Tree = 2
  PACO_N <- readDataSet(hfile[["PACO_N"]])
  PAN <- rep(1:Npatch,PACO_N)
  
  hfile$close_all()
  
  is.liana <- COI  <- COI2 <- fraction_LAI <- LAI_total_pa <- rep(NA,Npatch)
  
  for (ipa in seq(1,Npatch)){
    
    pos <- (PAN == ipa)
    
    hiteconow   <- hite[pos]
    dbhconow    <- dbh[pos]
    pftconow    <- pft[pos]
    laiconow    <- lai[pos]
    nplantconow <- nplant[pos]
    CAconow <- cai[pos]
    Nco <- length(hiteconow)
    
    is.liana[ipa] <- ifelse(any(PFTselect %in% pftconow),1,0)
    
    ##########################################################################
    
    fraction_LAI[ipa] <- sum(laiconow[pftconow %in% PFTselect])/sum(laiconow)
    
    # Compute COI
    LAItot <- sum(laiconow)
    LAI_total_pa[ipa] <- LAItot
    LAI_alpha <- (1-alpha_frac)*LAItot
    h_alpha <- interp1(x = as.vector(cumsum(laiconow)), 
                       y = as.vector(hiteconow), 
                       xi = ifelse(LAI_alpha<laiconow[1],laiconow[1],LAI_alpha))
    
    COI[ipa] <- sum(laiconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/sum(laiconow[hiteconow >= h_alpha])
    COI2[ipa] <- sum(CAconow[pftconow %in% PFTselect & hiteconow >= h_alpha])/sum(CAconow[hiteconow >= h_alpha])
    
  }
  
  return(list(COI = COI, is.liana = is.liana, COI2 = COI2, fraction_LAI = fraction_LAI,LAItot = LAI_total_pa))
}
