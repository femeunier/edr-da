#' Run PEcAn Workflow
#'
#' @export
run_workflow_pecan <- function(pecanxmlfile,run_all = TRUE){
  
  source('~/R/scripts/add.constants.to.config.r')
  
  # make sure always to call status.end
  options(warn=1)
  
  options(error=quote({
    PEcAn.utils::status.end("ERROR")
    PEcAn.remote::kill.tunnel(settings)
    if (!interactive()) {
      q(status = 1)
    }
  }))
  
  
  # ----------------------------------------------------------------------
  # PEcAn Workflow
  # ----------------------------------------------------------------------
  # Open and read in settings file for PEcAn run.
  # args <- commandArgs(trailingOnly = TRUE)
  # if (is.na(args[1])){
  #   settings <- PEcAn.settings::read.settings(pecanxmlfile)
  # } else {
  #   settings.file = args[1]
  #   settings <- PEcAn.settings::read.settings(settings.file)
  # }
   
  # If remote, open tunnel
  if (!is.localhost(settings$host)){
    # Close it first just in case
    system(paste("ssh -S ",settings$host$tunnel,settings$host$name,"-O exit",sep=" "))
    system(paste("ssh -nNf -o ControlMaster=yes -S",settings$host$tunnel,settings$host$name,sep=" "))
  }
  
  # Check for additional modules that will require adding settings
  if("benchmarking" %in% names(settings)){
    library(PEcAn.benchmark)
    settings <- papply(settings, read_settings_BRR)
  }
  
  if("sitegroup" %in% names(settings)){
    if(is.null(settings$sitegroup$nSite)){
      settings <- PEcAn.settings::createSitegroupMultiSettings(settings, sitegroupId = settings$sitegroup$id)
    } else {
      settings <- PEcAn.settings::createSitegroupMultiSettings(settings, sitegroupId = settings$sitegroup$id,nSite = settings$sitegroup$nSite)
    }
    settings$sitegroup <- NULL ## zero out so don't expand a second time if re-reading
  }
  
  # Update/fix/check settings. Will only run the first time it's called, unless force=TRUE
  settings <- PEcAn.settings::prepare.settings(settings, force=FALSE)
  
  # Write pecan.CHECKED.xml
  PEcAn.settings::write.settings(settings, outputfile = "pecan.CHECKED.xml")
  
  # start from scratch if no continue is passed in
  statusFile <- file.path(settings$outdir, "STATUS")
  if (length(which(commandArgs() == "--continue")) == 0 && file.exists(statusFile)) {
    file.remove(statusFile)
  }
  
  # Do conversions
  settings <- PEcAn.workflow::do_conversions(settings)
  
  # Query the trait database for data and priors
  if (PEcAn.utils::status.check("TRAIT") == 0){
    PEcAn.utils::status.start("TRAIT")
    settings <- PEcAn.workflow::runModule.get.trait.data(settings)
    PEcAn.settings::write.settings(settings, outputfile='pecan.TRAIT.xml')
    PEcAn.utils::status.end()
  } else if (file.exists(file.path(settings$outdir, 'pecan.TRAIT.xml'))) {
    settings <- PEcAn.settings::read.settings(file.path(settings$outdir, 'pecan.TRAIT.xml'))
  }
  
  
  # Run the PEcAn meta.analysis
  if(!is.null(settings$meta.analysis)) {
    if (PEcAn.utils::status.check("META") == 0){
      PEcAn.utils::status.start("META")
      PEcAn.MA::runModule.run.meta.analysis(settings)
      PEcAn.utils::status.end()
    }
  }
  
  # Write model specific configs
  if (PEcAn.utils::status.check("CONFIG") == 0){
    PEcAn.utils::status.start("CONFIG")
    settings <- PEcAn.workflow::runModule.run.write.configs(settings)
    PEcAn.settings::write.settings(settings, outputfile='pecan.CONFIGS.xml')
    PEcAn.utils::status.end()
  } else if (file.exists(file.path(settings$outdir, 'pecan.CONFIGS.xml'))) {
    settings <- PEcAn.settings::read.settings(file.path(settings$outdir, 'pecan.CONFIGS.xml'))
  }
  
  # Add constants to empty pfts
  add.constants.to.config(settings)
  
  
  if ((length(which(commandArgs() == "--advanced")) != 0) && (PEcAn.utils::status.check("ADVANCED") == 0)) {
    PEcAn.utils::status.start("ADVANCED")
    q();
  }
  
  if (!run_all){
    run_file <- file.path(settings$rundir, "runs.txt")
    run2_file <- file.path(settings$rundir, "runs_all.txt")
    dummy <- file.copy(run_file,run2_file)
    
    run_list <- readLines(con = run_file)
    singlerun_list <- run_list[1]
    writeLines(text = singlerun_list,con = run_file)
  }
  
  # Start ecosystem model runs
  if (PEcAn.utils::status.check("MODEL") == 0) {
    PEcAn.utils::status.start("MODEL")
    PEcAn.remote::runModule.start.model.runs(settings,stop.on.error=FALSE)
    PEcAn.utils::status.end()
  }
  
  if (!run_all){
    from_dir <- file.path(settings$modeloutdir,singlerun_list)
    for (irun in seq(2,length(run_list))){
      to_dir <- file.path(settings$modeloutdir,run_list[irun])
      args <- paste(file.path(from_dir,'*'),to_dir,sep=' ')
      # system2('cp',args)
    }
    system2('rm',run_file)
    dummy <- file.copy(run2_file,run_file)
    system2('rm',run2_file)
    
  }
  
  
}
