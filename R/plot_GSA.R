#' plot_GSA
#'
#' @export

plot_GSA <- function(tmppath,var,mapdf){
  
  pecanxmlfile <- file.path(tmppath,"pecan.CONFIGS.xml")
  set <- PEcAn.settings::read.settings(pecanxmlfile) 
  all_runs <- readLines(file.path(tmppath, "run", "runs.txt"))
  
  ens.runs <- all_runs[unlist(map(.x = all_runs,.f = function(run){
    runtype <- readLines(file.path(tmppath,"run",run, "README.txt"),n = 1)
    return((grepl("ensemble", runtype)) )
  }),use.names = FALSE)]
  writeLines(ens.runs,con = file.path(tmppath, "run", "runs_ens.txt"))
  
  
  # Read outputs of ensemble files
  tryCatch({
    samples.env <- new.env()
    load(file.path(tmppath,'samples.Rdata'),samples.env)
    samples.env <- map(.x=PFT,.f = function(pft){as.list.environment(samples.env)$ensemble.samples[[pft]]})
    
    df <- data.frame(matrix(unlist(samples.env), nrow=nrow(samples.env[[1]]), byrow=F))
    names(df) <- unlist(mapply(function(x,pft) names(x)<-paste(names(x),pft,sep='_'),samples.env,PFT))
    
    samples.env <- df
    
    var %>%
      walk(function(var.tmp){
        readLines(file.path(tmppath, "run", "runs_ens.txt")) %>%
          map_dfr(function(run.id) {
            sapply(
              PEcAn.utils::read.output(
                run.id,
                file.path(set$outdir, "out", run.id),
                as.numeric(lubridate::year(set$run$start.date)),
                as.numeric(lubridate::year(set$run$end.date)),
                var.tmp
              ),
              mean,
              na.rm = TRUE
            ) %>%
              t %>%
              as.data.frame()
          }) %>%
          cbind(samples.env) %>%
          saveRDS(file.path(set$outdir,"out", paste0(var.tmp, "_GSA.rds")))
      })
    
  },
  error = function(e) {
    PEcAn.logger::logger.warn("Probably runs were faulty. Can't read the outputs")
  }
  )
  
  
  # GSA
  ready.for.plot <- list.files(file.path(set$outdir,"out"),"*_GSA.rds", recursive = T, full.names = T) %>%
    map_dfr(function(GSA_file){
      PFT <- strsplit(GSA_file,"/")[[1]][2]
      tmp.data <-readRDS(GSA_file)
      target.var <-names(tmp.data)[1]
      vars <-names(tmp.data)[2:length(names(tmp.data))]
      
      lmdf <- tmp.data
      anovobj<-aov(lm(as.formula(paste0(target.var," ~ .^2")),data=lmdf))
      allssq<-summary(anovobj)[[1]][,2]
      varn<-names((anovobj)$coefficients)[-1]
      
      ##########
      sensivities<-sapply(vars,function(nn){
        
        SSQ <- sum(allssq,na.rm = T)
        SSQ_d <- sum(allssq[which(nn == varn)],na.rm = T)
        SSQ_i <- sum(allssq[which(!(nn == varn) & grepl(nn, varn, fixed=T))],na.rm = T)
        return((SSQ_d+(SSQ_i/2))/SSQ)
        
      },USE.NAMES = T) %>%
        t %>%
        as.data.frame() %>%
        mutate(Var=target.var, 
               PFT=PFT)
    })
  
  pfts <- paste(paste0('_',as.character(df_PFT$names)),collapse ='|')
  
  ready.data.plot<- ready.for.plot %>%
    gather(Param,Value, -c(Var,PFT)) %>% select(Var,Param,Value) %>% 
    mutate(Param_name = gsub(pfts,'',Param)) %>%
    mutate(PFT =  gsub(paste(paste0(unique(Param_name),'_'),collapse ='|'),'',Param)) %>%
    group_by(Var, PFT, Param_name)%>%
    summarise(MeanV=mean(Value, na.rm=T)) %>%
    filter(MeanV>0.0000001, !is.nan(MeanV), !is.na(MeanV)) %>%
    arrange(desc(MeanV)) %>% rename(Param = Param_name) %>% ungroup() %>%
    mutate(Var=as.factor(Var),PFT=as.factor(PFT),Param=as.factor(Param)) 
  
  # complete
  ready.data.plot_all <- ready.data.plot %>% complete(Var,PFT,Param)
  
  wdth=0.5
  plot <- 
    ready.data.plot_all %>%
    ungroup() %>%
    ggplot(aes(Param, MeanV*100,group = PFT))+
    geom_bar(aes(fill=PFT), stat="identity", width = wdth,
             position=position_dodge(width = 0.5)) +
    facet_wrap(~Var)+
    scale_y_continuous(limits = c(0,max(ready.data.plot_all$MeanV*100))) +
    coord_flip()+
    scale_x_discrete() +
    scale_fill_manual(values = as.character(mapdf$col),limits = as.character(mapdf$pfts)) +
    theme_light(base_size = 13) +
    labs(y = "Variance explained by each parameter (%)", x = "", title = "",
         subtitle = "")
  
  ggsave(filename = file.path(tmppath,"pfts","all",paste0("GSA_all",".png")),dpi = 300,
         width = 10, height = 5, plot = plot)
  
}