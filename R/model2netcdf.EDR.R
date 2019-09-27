#' Write PEcAn standard output files
#'
#' @export
model2netcdf.EDR <- function(outdir,sitelat,sitelon,start_date,par.wl = 400:2499, 
                             nir.wl = 2500,patches = FALSE){
  
  start_year <- lubridate::year(start_date)
  
  output.path <- file.path(outdir,"RTM")
  if (patches){  
    flist <- list()
    
    flist[["albedo"]] <- dir(file.path(output.path,"patches"),"albedo_nir") #  Albedo NIR file
  } else {
    flist <- list()
    flist[["albedo"]] <- dir(output.path,"albedo_nir")           #  Albedo NIR file 
  }
  
  # check if there are files
  file.check <- sapply(flist, function (f) length(f) != 0)
  
  if (!any(file.check)) {
    # no output files
    PEcAn.logger::logger.warn("WARNING: No output files found for :", outdir)
    return(NULL)
  }
  
  EDR_outputs <- read.EDR.output(output.path,patches = patches)
  out <- list(PAR = EDR_outputs[,1:length(par.wl)],
              NIR = EDR_outputs[,(length(par.wl)+1):(length(par.wl)+length(nir.wl))])
  
  Npatches <- ifelse(patches,as.numeric(readLines(con = file.path(output.path, 
                                                                  "patches/Npatches.dat")))[[1]],1)
  # create lat/long/patches/lambda nc variables
  lat <- ncdf4::ncdim_def("lat", "degrees_north",
                          vals = as.numeric(sitelat),
                          longname = "station_latitude")
  lon <- ncdf4::ncdim_def("lon", "degrees_east",
                          vals = as.numeric(sitelon),
                          longname = "station_longitude")
  
  patch <- ncdf4::ncdim_def("Patches", "-",
                            vals = (1:Npatches),
                            longname = "Number of patches")
  
  WL_par <- ncdf4::ncdim_def("par.wl", "nm",
                             vals = par.wl,
                             longname = "PAR wavelength")
  
  WL_nir <- ncdf4::ncdim_def("nir.wl", "nm",
                             vals = nir.wl,
                             longname = "NIR wavelength")
  
  WL_all <- ncdf4::ncdim_def("Full", "nm",
                             vals = c(par.wl,nir.wl),
                             longname = "All wavelengths")
  
  
  # Filters 
  vis <- red <- rep(NA,length(par.wl))
  
  vis[par.wl>=380 & par.wl<= 740] <- 1
  vis_filter <- t(matrix(rep(vis,Npatches),ncol = Npatches))
  
  red[par.wl>=625 & par.wl<= 740] <- 1
  red_filter <- t(matrix(rep(red,Npatches),ncol = Npatches))
  
  out[['vis']] <- vis_filter*out[['PAR']]
  out[['red']] <- red_filter*out[['PAR']]
  out[['All']] <- EDR_outputs
    
  # Create netcdf files
  nc_var <- list()
  nc_var[[1]] <- ncdf4::ncvar_def("PAR", units = "-", dim = list(lon, lat, patch, WL_par), missval = -999, 
                                  longname = "PAR ecosystem reflectance")
  
  nc_var[[2]] <- ncdf4::ncvar_def("NIR", units = "-", dim = list(lon, lat, patch, WL_nir), missval = -999, 
                                  longname = "NIR ecosystem reflectance")
  
  nc_var[[3]] <- ncdf4::ncvar_def("Vis", units = "-", dim = list(lon, lat, patch, WL_par), missval = -999, 
                                  longname = "Vis ecosystem reflectance")
  
  nc_var[[4]] <- ncdf4::ncvar_def("Red", units = "-", dim = list(lon, lat, patch, WL_par), missval = -999, 
                                  longname = "Red ecosystem reflectance")
  
  nc_var[[5]] <- ncdf4::ncvar_def("Spectrum", units = "-", dim = list(lon, lat, patch, WL_all), missval = -999, 
                                  longname = "All ecosystem reflectance")
  
  nc <- ncdf4::nc_create(file.path(outdir, paste(start_year, "nc", sep = ".")), nc_var)
  
  for (i in seq(nc_var)){
    ncdf4::ncvar_put(nc, nc_var[[i]],out[[i]])
  }
  ncdf4::nc_close(nc)
}  