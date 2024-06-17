library(dplyr)
library(terra)
library(tidyr)
library(EnvStats)



############### Summary of annual number of days conducive to low/moderate/high burn severity ###################
Ecozone <- terra::rast("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones(ESRI_102002)_resampleFWI2020.tif")
Pre_nfiLandCover <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_att_nfiLandCover_SW_2020_v1_MODEresampleFWI2020(ESRI_102002)_3Trees.tif")
VegeType <- terra::resample(Pre_nfiLandCover, Ecozone, method="near")
crop_extent <- terra::ext(Ecozone); VegeType_crop <- terra::crop(VegeType, crop_extent)
VegeType_mask <- terra::mask(VegeType_crop, Ecozone); names(VegeType_mask) <- "Pre_nfiLandCover"
ValidCells <- terra::cells(VegeType_mask) 
Pixel_summary <- terra::extract(VegeType_mask, ValidCells, xy=TRUE)
Pixel_summary$Ecozone <- terra::extract(Ecozone, Pixel_summary[,c("x","y")], xy=FALSE)[, names(Ecozone)]
Pixel_summary$Pixel_ID <- seq(1, nrow(Pixel_summary), 1)
folder <- "Predictions_Classification/Predictions_Classification_AllDays/"
Years <- seq(1981,2020,1)
mi <- "em"
Pixel_summary_mi_low <- Pixel_summary; Pixel_summary_mi_low$dNBR_Class <- "Low"
Pixel_summary_mi_mode <- Pixel_summary; Pixel_summary_mi_mode$dNBR_Class <- "Mode"
Pixel_summary_mi_high <- Pixel_summary; Pixel_summary_mi_high$dNBR_Class <- "High"
for (yi in Years) {
  print(paste("Categorization:",mi,";Year:",yi))
  rast_fyi <- terra::rast(paste(folder, "/Predictions_",mi,"_Year",yi,"_March-November.tif", sep=""))
  rast_fyi_lowdays <- terra::app(rast_fyi, function(x) sum(x == 1, na.rm = TRUE))
  rast_fyi_modedays <- terra::app(rast_fyi, function(x) sum(x == 2, na.rm = TRUE))
  rast_fyi_highdays <- terra::app(rast_fyi, function(x) sum(x == 3, na.rm = TRUE))
  Pixel_summary_mi_low[, paste0("Year_",yi)] <- terra::extract(rast_fyi_lowdays, Pixel_summary[,c("x","y")], xy=FALSE)[, names(rast_fyi_lowdays)]
  Pixel_summary_mi_mode[, paste0("Year_",yi)] <- terra::extract(rast_fyi_modedays, Pixel_summary[,c("x","y")], xy=FALSE)[, names(rast_fyi_modedays)]
  Pixel_summary_mi_high[, paste0("Year_",yi)] <- terra::extract(rast_fyi_highdays, Pixel_summary[,c("x","y")], xy=FALSE)[, names(rast_fyi_highdays)]
}
Pixel_summary_mi <- rbind(Pixel_summary_mi_low, Pixel_summary_mi_mode, Pixel_summary_mi_high)
Years_40 <- seq(1981,2020,1)
Years_20a <- seq(1981,2000,1)
Years_20b <- seq(2001,2020,1)
Pixel_summary_mi$AverageDays <- apply(Pixel_summary_mi[, paste0("Year_", Years_40)], 1, function(x) mean(x, na.rm=TRUE)) 




#################### Trend test #######################
MK_FDR <- function(time_series_1pixel) {
  time_series_1pixel <- as.numeric(time_series_1pixel)
  time_series_1pixel <- time_series_1pixel[!is.na(time_series_1pixel)]
  if ( (sd(time_series_1pixel) == 0) | length(time_series_1pixel)<8 ) {
    trend_pi <- NA
  } else {
    ts_result <- kendallTrendTest(time_series_1pixel)
    adjusted_pvalues <- p.adjust(ts_result$p.value, method = "fdr")
    alpha_FDR <- 0.1
    alpha_global <- 0.5 * alpha_FDR
    if (adjusted_pvalues <= alpha_global) {
      trend_pi <- as.numeric(ts_result$estimate["slope"])
    } else {
      trend_pi <- NA
    }
  }
  return(trend_pi)
}

Pixel_summary_mi$Trend1981_2020 <- apply(Pixel_summary_mi[, paste0("Year_", Years_40)], 1, MK_FDR)
Pixel_summary_mi$Trend1981_2000 <- apply(Pixel_summary_mi[, paste0("Year_", Years_20a)], 1, MK_FDR)
Pixel_summary_mi$Trend2001_2020 <- apply(Pixel_summary_mi[, paste0("Year_", Years_20b)], 1, MK_FDR)
saveRDS(Pixel_summary_mi, paste0("PredictClassiAlldays_",mi,"_AnnualSeverityClassDays_Trends.rds"))
for (ci in unique(Pixel_summary_mi$dNBR_Class)) {
  data_mi <- subset(Pixel_summary_mi, dNBR_Class == ci)
  raster_mi_average <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$AverageDays))
  names(raster_mi_average) <- paste0("AverageAnnualDays_", ci)
  writeRaster(raster_mi_average, paste0("PredictClassiAlldays_",mi,"_Annual",ci,"SeverityDays_Average.tif"), overwrite=TRUE)
  
  raster_mi_Trend1981_2020 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend1981_2020))
  names(raster_mi_Trend1981_2020) <- paste0("Trend1981_2020_", ci)
  writeRaster(raster_mi_Trend1981_2020, paste0("PredictClassiAlldays_",mi,"_Annual",ci,"SeverityDays_Trend1981-2020.tif"), overwrite=TRUE)
  
  raster_mi_Trend1981_2000 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend1981_2000))
  names(raster_mi_Trend1981_2000) <- paste0("Trend1981_2000_", ci)
  writeRaster(raster_mi_Trend1981_2000, paste0("PredictClassiAlldays_",mi,"_Annual",ci,"SeverityDays_Trend1981-2000.tif"), overwrite=TRUE)
  
  raster_mi_Trend2001_2020 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend2001_2020))
  names(raster_mi_Trend2001_2020) <- paste0("Trend2001_2020_", ci)
  writeRaster(raster_mi_Trend2001_2020, paste0("PredictClassiAlldays_",mi,"_Annual",ci,"SeverityDays_Trend2001-2020.tif"), overwrite=TRUE)
}




################# t-test for comparing the two 20-year periods ###################
Years_40 <- seq(1981,2020,1)
Years_20a <- seq(1981,2000,1)
Years_20b <- seq(2001,2020,1)
mi <- "em"
Pixel_summary_mi <- readRDS(paste0("Predictions_Classification/Predictions_Classification_AnnualSeverityDays/PredictClassiAlldays_",mi,"_AnnualSeverityClassDays_Trends.rds"))
Ttest_annualDays <- data.frame()
for (ci in unique(Pixel_summary_mi$dNBR_Class)) {
  data_mi <- subset(Pixel_summary_mi, dNBR_Class == ci)
  for (zi in unique(data_mi$Ecozone)) {
    data_zi <- subset(data_mi, Ecozone==zi)
    data_zi$MeanDays_1981_2020 <- apply(data_zi[, paste0("Year_", Years_40)], 1, function(x) mean(x, na.rm=TRUE)) 
    data_zi$MeanDays_1981_2000 <- apply(data_zi[, paste0("Year_", Years_20a)], 1, function(x) mean(x, na.rm=TRUE)) 
    data_zi$MeanDays_2001_2020 <- apply(data_zi[, paste0("Year_", Years_20b)], 1, function(x) mean(x, na.rm=TRUE)) 
    Ttest_Paired <- t.test(data_zi$MeanDays_1981_2000, data_zi$MeanDays_2001_2020, paired=TRUE)
    Ttest_annualDays <- rbind(Ttest_annualDays, data.frame(dNBR_Class=ci, Ecozone=zi, MeanDays_1981_2020=mean(data_zi$MeanDays_1981_2020, na.rm=TRUE), 
                                                           MeanDays_1981_2000 = mean(data_zi$MeanDays_1981_2000, na.rm=TRUE), 
                                                           MeanDays_2001_2020 = mean(data_zi$MeanDays_2001_2020, na.rm=TRUE),
                                                           Difference_0120_8100 = mean(data_zi$MeanDays_2001_2020 - data_zi$MeanDays_1981_2000, na.rm=TRUE),
                                                           Ttest_p_Mean20yrs = as.numeric(Ttest_Paired$p.value)))
  }
}
write.table(Ttest_annualDays, file=paste0("PredictClassiAlldays_",mi,"_AnnualSeverityClassDays_MeanTtest_Summary.csv"), sep=",", append=T, col.names=T, row.names=F)



