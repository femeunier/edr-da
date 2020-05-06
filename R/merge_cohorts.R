#' @param patch patch
#' @return 
#' @author Félicien Meunier
#' @export

merge_cohorts <- function(patch,dbh_diff=1,nplant.col = "nplant"){
  
  patch_in <- patch
  patch_out <- patch_in[1,]
  patch_out <- patch_out[-c(1),] # empty
  
  nrow_beg <- nrow(patch_in)
  nrow_end <- nrow(patch_out)
  iter = 1
  while((nrow_beg != nrow_end) & iter < 50){
    for (i in seq(2,nrow_beg)){
      coh1 <- patch_in[i-1,]
      coh2 <- patch_in[i,]
      if ((coh1[["pft"]] == coh2[["pft"]]) & abs(coh1[["dbh"]]-coh2[["dbh"]]) < dbh_diff){ # merge
        
        coh2["dbh"] <-  (coh2[nplant.col]* coh2["dbh"] +  coh1[nplant.col]* coh1["dbh"])/( coh1[nplant.col] +  coh2[nplant.col])
        coh2[nplant.col] <-  coh1[nplant.col] +  coh2[nplant.col]
        coh1[nplant.col] <- 0
        
        patch_in[i-1,] <- coh1
        patch_in[i,] <- coh2
      }
    }
    
    patch_out <- patch_in[patch_in[,nplant.col]>0,]
    
    nrow_beg <- nrow(patch_in)
    nrow_end <- nrow(patch_out)
    
    patch_in <- patch_out
    nrow_beg <- nrow(patch_in)
    iter <- iter +1
  }
  return(patch_in)
}