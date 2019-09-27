#' Plot variance decomposition plot for all PFTs
#'
#' @export
VDP_allPFTs <- function(variable,mapdf){
  
  settings$sensitivity.analysis$variable=variable
  
  if (is.null(settings$sensitivity.analysis$ensemble.id)){
    settings$sensitivity.analysis$ensemble.id <- "NOENSEMBLEID"
    settings$ensemble$ensemble.id <-"NOENSEMBLEID"
  }
  
  ensemble_id <- settings$sensitivity.analysis$ensemble.id
  
  st_year=substr(settings$run$start.date,1,4)
  end_year=substr(settings$run$end.date,1,4)
  
  runModule.get.results(settings)
  runModule.run.sensitivity.analysis(settings,plot = FALSE)
  
  results_file=paste("sensitivity.results",ensemble_id,variable,st_year,end_year,"Rdata",sep=".")
  load(file.path(settings$outdir,results_file))
  PFTs=ls(sensitivity.results)
  npft=length(PFTs)
  
  for (i in seq(npft)){
    if (i==1){
      Variances_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances
      sensitivities_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$sensitivities
      partial_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$partial.variances
      elasticity_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$elasticities
      coef.vars_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$coef.vars
      
      names_all=paste(names(sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances),PFTs[i],sep="_")
      names_all=names(sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances)
      
      names(Variances_temp) <- names_all
      names(sensitivities_temp) <- names_all
      names(partial_temp) <- names_all
      names(sensitivities_temp) <- names_all
      names(coef.vars_temp) <- names_all
      
      Variances_all <- Variances_temp
      sensitivities_all <- sensitivities_temp
      partial_all <- partial_temp
      elasticity_all <- elasticity_temp
      coef.vars_all <- coef.vars_temp
      PFT_all <- rep(PFTs[i],length(names_all))
      
    } else {
      Variances_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances
      sensitivities_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$sensitivities
      partial_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$partial.variances
      elasticity_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$elasticities
      coef.vars_temp=sensitivity.results[[PFTs[i]]]$variance.decomposition.output$coef.vars
      
      names_all=paste(names(sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances),PFTs[i],sep="_")
      names_all=names(sensitivity.results[[PFTs[i]]]$variance.decomposition.output$variances)
      
      names(Variances_temp) <- names_all
      names(sensitivities_temp) <- names_all
      names(partial_temp) <- names_all
      names(elasticity_temp) <- names_all
      names(coef.vars_temp) <- names_all
      
      Variances_all <- c(Variances_all,Variances_temp)
      sensitivities_all <- c(sensitivities_all,sensitivities_temp)
      partial_all <- c(partial_all,partial_temp)
      elasticity_all <- c(elasticity_all,elasticity_temp)
      coef.vars_all <- c(coef.vars_all,coef.vars_temp)
      PFT_all <- c(PFT_all,rep(PFTs[i],length(names_all)))
    }
  }
  
  Ordering=sort(Variances_all,index.return=TRUE)
  
  OP_all=ls()
  OP_all$coef.vars=coef.vars_all[Ordering$ix]
  OP_all$elasticities=elasticity_all[Ordering$ix]
  OP_all$sensitivities=sensitivities_all[Ordering$ix]
  OP_all$variances=Variances_all[Ordering$ix]
  OP_all$partial.variances=OP_all$variances/sum(OP_all$variances)
  OP_all$PFT_all=PFT_all[Ordering$ix]
  
  fname <- sensitivity.filename(settings, "variance.decomposition",
                                "pdf", all.var.yr = FALSE, pft = "all",
                                ensemble.id = ensemble_id, variable = variable,
                                start.year = st_year, end.year = end_year)
  C=list()
  all_col=mapdf$col[match(sort(OP_all$PFT_all),mapdf$pfts)]
  C$col<- as.vector(all_col[!duplicated(all_col)])
  C$PFT=OP_all$PFT_all
  
  vd.plots <- PEcAn.uncertainty::plot_variance_decomposition(plot.inputs=OP_all,Col=C$col)
  pdf(fname, width = 11, height = 8)
  do.call(grid.arrange, c(vd.plots, ncol = 4))
  dev.off()
}  
