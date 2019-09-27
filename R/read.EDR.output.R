#' Read all EDR outputs
#'
#' @export
read.EDR.output <- function(output.path,par.wl = 400:2499, 
                            nir.wl = 2500, patches = FALSE, ...){
  
  if (!patches) {
    albedo <- PEcAnRTM::get.EDR.output(output.path)
  }
  else {
    output.path_patch = paste(output.path, "patches", sep = "/")
    Npatches <- as.numeric(readLines(con = file.path(output.path, 
                                                     "patches/Npatches.dat")))[[1]]
    albedo = matrix(NA, Npatches, length(par.wl) + length(nir.wl))
    for (ipatch in seq(1, Npatches)) {
      albedo[ipatch, ] <- PEcAnRTM::get.EDR.output(output.path_patch, 
                                                   ipatch)
    }
  }
  return(albedo)
}

