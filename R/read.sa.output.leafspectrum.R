#' read.sa.output.leafspectrum
#'
#' @export

read.sa.output.leafspectrum <- function(traits, quantiles, pecandir, outdir, variable, sa.run.ids = NULL,
                                    par.wl,nir.wl,spectra_list_all,pft.name){
  
  Lambdas <- c(par.wl,nir.wl)
  if (is.null(sa.run.ids)) {
    samples.file <- file.path(pecandir, "samples.Rdata")
    if (file.exists(samples.file)) {
      load(samples.file)
      sa.run.ids <- runs.samples$sa
    }
    else {
      PEcAn.logger::logger.error(samples.file, "not found, this file is required by the read.sa.output function")
    }
  }
  sa.output <- data.frame()
  
  for (trait in traits) {
    for (quantile in quantiles) {
      run.id <- sa.run.ids[[pft.name]][quantile, trait]
      
      out.tmp <- spectra_list_all[[run.id]]
      data_temp <- out.tmp[[pft.name]]
   
      rownames(data_temp) <- Lambdas
      
      df_temp <- as.data.frame(melt(t(data_temp))) %>% rename(var = Var1, lambda = Var2, Reflectance = value) %>%
        mutate(traits=trait,quantiles = quantile)
      sa.output <- rbind(sa.output,df_temp)
    }
    PEcAn.logger::logger.info("reading sensitivity analysis output for model run at ", 
                              quantiles, "quantiles of trait", trait)
  }
  # sa.output <- as.data.frame(sa.output)
  return(sa.output)
}



