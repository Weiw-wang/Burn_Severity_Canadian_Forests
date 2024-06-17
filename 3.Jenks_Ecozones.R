library(classInt)


folder_Zones <- "DataPreparation/ValidSamples_Years_Ecozones_Factors/"
UniqZones <- c(18, 4, 9, 11, 12, 14, 15, 16, 17, 19)
Jenks <- data.frame()
for(zi in UniqZones) {
  Samples_zi <- readRDS(paste0(folder_Zones, "ValidSamples(Seasons3_Trees3_dNBR271300)_Zone",zi,"_Factors.rds"))
  Jenks_zi <- data.frame()
  for (ti in 1:100) {
    print(paste("Ecozone:", zi, "; Repeat:", ti))
    data_subset <- sample(Samples_zi$dNBR1000, 5000)
    nb_classes <- classIntervals(data_subset, n = 3, style = "jenks", factor = FALSE)
    Jenks_zi_ti <- nb_classes$brks
    names(Jenks_zi_ti) <- c("[Min","Low]","Moderate]","Max]")
    Jenks_zi_ti <- t(data.frame(Jenks_zi_ti))
    Jenks_zi <- rbind(Jenks_zi, Jenks_zi_ti)
    Jenks <- rbind(Jenks, cbind(data.frame(Ecozone=zi, Repeats=ti), Jenks_zi_ti) )
  }
  Jenks_zi_ave <- colMeans(Jenks_zi);  Jenks_zi_ave <- t(data.frame(Jenks_zi_ave))
  Jenks <- rbind(Jenks, cbind(data.frame(Ecozone=zi, Repeats="Average"), Jenks_zi_ave) )
}
write.table(Jenks, file="ValidSamples(Seasons3_Trees3_dNBR271300)_ZonesFactors_Jenks_Ecozones.csv", sep=",", append=T, col.names=T, row.names=F)


