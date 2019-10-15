# Merge gigante plots

library(plyr)

site_codes <- c("Gigante_control","Gigante_removal")
site_dir <- "/home/carya/data/Gigante"
final_path <- file.path(site_dir,"Gigante_all.lat9.000lon-79.000")

density <- 1/(60*60)

css_all <- c()
pss_all <- c()

Npatch_tot <- 0
Ncohort_tot <- 0

for (isite in seq(site_codes)){
  css_file <- tail(list.files(site_dir, pattern=paste0("^",site_codes[isite],".*","css"),full.names=TRUE), 1)
  pss_file <- tail(list.files(site_dir, pattern=paste0("^",site_codes[isite],".*","pss"),full.names=TRUE), 1)
  
  css_dat <- PEcAn.ED2::read_css(css_file)
  pss_dat <- read.table(pss_file, header = TRUE)
  
  # Patch
  old_patches <- pss_dat[["patch"]]
  Npatch <- length(unique(old_patches))
  new_patches <- Npatch_tot + seq(1,Npatch)
  pss_dat[["patch"]] <- new_patches
  css_dat[["patch"]] <- mapvalues(css_dat[["patch"]],old_patches,new_patches)
  
  # Cohort
  old_cohorts <- css_dat[["cohort"]]
  Ncohort <- length(unique(old_cohorts))
  new_cohorts <- Ncohort_tot + seq(1,Ncohort)
  css_dat[["cohort"]] <- mapvalues(css_dat[["cohort"]],old_cohorts,new_cohorts)
  
  css_all <- rbind(css_all,css_dat)
  pss_all <- rbind(pss_all,pss_dat)
  
  Npatch_tot <- Npatch_tot + Npatch
  Ncohort_tot <- Ncohort_tot + Ncohort
}

site_file <- tail(list.files(site_dir, pattern=paste0("^",".*","site"),full.names=TRUE), 1)

css_all[["n"]] <- density
# PEcAn.ED2::write_css(css_all,final_path,latitude = -9,longitude = -79)
# PEcAn.ED2::write_pss(pss_all,final_path,latitude = -9,longitude = -79)

file.copy(site_file,paste(final_path,"site",sep="."))
dummy   = write.table(css_all,
                      paste0(final_path,".css")
                      , append    = FALSE
                      , quote     = FALSE
                      , sep       = " "
                      , row.names = FALSE
                      , col.names = TRUE)
dummy   = write.table(pss_all,paste0(final_path,".pss")
                      , append    = FALSE
                      , quote     = FALSE
                      , sep       = " "
                      , row.names = FALSE
                      , col.names = TRUE)

css<-read.table("/home/carya/data/Gigante/Gigante_all.lat9.000lon-79.000.css", header = TRUE)
# pss<-read.table("/home/carya/data/Gigante/Gigante_all.lat9.000lon-79.000.lat-9lon-79.pss", header = TRUE)
