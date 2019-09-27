#' read.sa.output.spectrum
#'
#' @export

read.sa.output.spectrum <- function(traits, quantiles, pecandir, outdir, variable, sa.run.ids = NULL,
                                    par.wl,nir.wl){
  
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
      
      out.tmp <- read.output(runid = run.id, outdir = file.path(outdir, 
                                                                run.id), start.year = start.year, end.year = end.year, 
                             variables = variable, pft.name = NULL)
      data_temp <- out.tmp[[variable]]
      Npatches <- ifelse(!is.null(dim(data_temp)),dim(data_temp)[1],1)
      rownames(data_temp) <- 1:Npatches
      colnames(data_temp) <- Lambdas
      
      df_temp <- as.data.frame(melt(t(data_temp))) %>% rename(patch = Var2, lambda = Var1, Reflectance = value) %>%
        mutate(traits=trait,quantiles = quantile)
      sa.output <- rbind(sa.output,df_temp)
    }
    PEcAn.logger::logger.info("reading sensitivity analysis output for model run at ", 
                              quantiles, "quantiles of trait", trait)
  }
  # sa.output <- as.data.frame(sa.output)
  return(sa.output)
}



