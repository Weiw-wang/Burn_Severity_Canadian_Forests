

####################### Add additional factors & save the data for each year & each ecozone ################################
years <- seq(2001, 2020, 1)
folder_years <- "DataPreparation/ValidSamples_Years"
files_years <- list.files(folder_years, pattern="*.rds$", full.names = TRUE)
PSD_thresholds <- read.csv("DataPreparation/FFMC_thresholds_p50_Wang2023.csv", head=T)
UniqZones <- c(4, 9, 11, 12, 14, 15, 16, 17, 18, 19)
for (zi in UniqZones) {
  for (yi in years){
    print(paste("Ecozone", zi, "; Year", yi))
    Samples_yi_zi <- c()
    fi_id <- grep(yi, files_years)
    for (fi in fi_id) {
      Samples_yifi <- readRDS(files_years[fi])
      Samples_yifi <- subset(Samples_yifi, Ecozone==zi)
      Samples_yi_zi <- rbind(Samples_yi_zi, Samples_yifi)
    }
    if (nrow(Samples_yi_zi)==0) {
      next
    }
    Samples_yi_zi_sub <- subset(Samples_yi_zi, (Pre_nfiLandCover %in% c(5,6,7)) & (dNBR1000>26 & dNBR1000<1301), select=-c(FireYear_growth))
    if (nrow(Samples_yi_zi_sub)==0) {
      next
    }
    ## derive seasons
    date_Spr <- as.Date(paste0(yi,"-03-01")) #spring (March-May)
    Jday_Spr <- as.numeric(date_Spr - as.Date(paste0(yi,"-01-01"))) + 1
    date_Sum <- as.Date(paste0(yi,"-06-01")) #summer (June-August)
    Jday_Sum <- as.numeric(date_Sum - as.Date(paste0(yi,"-01-01"))) + 1
    date_Aut <- as.Date(paste0(yi,"-09-01")) #autumn (September-November)
    Jday_Aut <- as.numeric(date_Aut - as.Date(paste0(yi,"-01-01"))) + 1
    date_Win <- as.Date(paste0(yi,"-12-01")) #winter (December-February)
    Jday_Win <- as.numeric(date_Win - as.Date(paste0(yi,"-01-01"))) + 1
    Samples_yi_zi_sub$Season <- cut(Samples_yi_zi_sub$GrowthDay,
                                    breaks = c(1, Jday_Spr, Jday_Sum, Jday_Aut, Jday_Win, 367), 
                                    labels = c(0,1,2,3,0), #c("Winter","Spring", "Summer", "Autumn","Winter"),
                                    include.lowest = TRUE, right = FALSE)
    ## derive potential spread days
    PSD_thres_zi <- PSD_thresholds$p50_FFMC[which(PSD_thresholds$Ecozone==zi)]
    Samples_yi_zi_sub$PSD <- cut(Samples_yi_zi_sub$FFMC,
                                 breaks = c(min(Samples_yi_zi_sub$FFMC)-1, PSD_thres_zi, max(Samples_yi_zi_sub$FFMC)),
                                 labels = c(0,1), ## Not OR Is a potential fire spread day
                                 include.lowest = FALSE, right = TRUE)
    saveRDS(Samples_yi_zi_sub, paste0("ValidSamples(Seasons3_Trees3_dNBR271300)_Zone",zi,"_Year",yi,".rds"))
  }
}



