
Categorize_dNBR1000 <- function(AllSamples, method_interest, nb_thresholds) {
  if (method_interest=="nb") { ### Categorize with natural breaks (Jenks)
    AllSamples$dNBRclass <- cut(AllSamples$dNBR1000,
                                breaks = c(26, round(nb_thresholds[1]), round(nb_thresholds[2]), 1301), 
                                labels = c(1,2,3), #Three classes
                                include.lowest = FALSE, right = TRUE)
  } else if (method_interest=="em") { ### Categorize with empirical breaks
    AllSamples$dNBRclass <- cut(AllSamples$dNBR1000,
                                breaks = c(26,276,517,1301),
                                labels = c(1,2,3), 
                                include.lowest = FALSE, right = TRUE)
  }
  return(AllSamples)
}


############## MAIN ##################
folder_ZonesYears <- "DataPreparation/ValidSamples_Years_Ecozones"
files_ZonesYears <- list.files(folder_ZonesYears, pattern="*.rds$", full.names = TRUE)
Jenks <- read.csv("DataPreparation/Jenks/ValidSamples(Seasons3_Trees3_dNBR271300)_ZonesFactors_Jenks_Ecozones.csv", head=T)
Categorization_Methods <- c("em", "nb")
UniqZones <- c(18, 4, 9, 11, 12, 14, 15, 16, 17, 19)
for (mi in Categorization_Methods) {
  for(zi in UniqZones) {
    Samples_MajorClass <- c()
    file_id <- grep(paste0("Zone",zi), files_ZonesYears)
    Jenks_zi <- subset(Jenks, Ecozone==zi & Repeats=="Average")
    nb_thresholds_zi <- Jenks_zi[c(grep("Low", names(Jenks_zi)), grep("Moderate", names(Jenks_zi)))]
    for (fi in file_id) {
      print(files_ZonesYears[fi])
      Samples_zi_fi <- readRDS(files_ZonesYears[fi])
      Samples_zi_fi$YearRow <- seq(1, nrow(Samples_zi_fi),1)
      Samples_zi_fi <- subset(Samples_zi_fi, select=c(FireYear_dNBR, YearRow, GrowthDay, Season, dNBR1000, FireID, DEM30_Aspect, Pre_nfiLandCover, PSD))
      Samples_zi_fi <- Categorize_dNBR1000(Samples_zi_fi, mi, nb_thresholds_zi)
      Fires <- unique(Samples_zi_fi$FireID)
      for (firei in Fires) {
        Samples_firei <- Samples_zi_fi[which(Samples_zi_fi$FireID == firei),]
        freq_firei <- table(Samples_firei$dNBRclass)
        majorClass_firei <- names(freq_firei)[which.max(as.numeric(freq_firei))]
        Samples_MajorClass_firei <- Samples_firei[which(Samples_firei$dNBRclass==majorClass_firei), ]
        Samples_MajorClass <- rbind(Samples_MajorClass, Samples_MajorClass_firei)
      }
    }
    saveRDS(Samples_MajorClass, paste0("ValidSamples(Seasons3_Trees3_dNBR271300)_", mi,"_Zone",zi,"_Factors_MajorClass.rds"))
  }
}


