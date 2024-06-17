library(dplyr)
library(terra)
library(tidyr)
library(EnvStats)



##################### Compute monthly high-severity potential for each year #############################
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
Pixel_potentials_mi <- c()
for (yi in Years) {
  print(paste("Categorization:",mi,";Year:",yi))
  raster_highsevP <- c()
  rast_fyi <- terra::rast(paste(folder, "/Predictions_",mi,"_Year",yi,"_March-November.tif", sep=""))
  days_fyi <- as.numeric(names(rast_fyi))
  date_objs <- as.Date(days_fyi - 1, origin = paste0(yi, "-01-01"))
  months_fyi <- as.numeric(format(date_objs, "%m"))
  for (mo in unique(months_fyi)) {
    data_yi_mi <- rast_fyi[[which(months_fyi==mo)]]
    alldays_mo <- as.numeric(difftime(as.Date(paste(yi, mo+1, "01", sep = "-")), as.Date(paste(yi, mo, "01", sep = "-")), units = "days"))
    if (length(names(data_yi_mi)) < alldays_mo) {
      next
    }
    Potential_calculate <- function(x, alldays, dNBR_Class) {
      valid_count <- sum(!is.na(x))
      if (valid_count < alldays) {
        return(NA)
      } else {
        class_count <- sum(x == dNBR_Class, na.rm=TRUE)
        percentage <- (class_count / valid_count) * 100
        return(percentage)
      }
    }
    Percent_classH <- terra::app(data_yi_mi, fun=Potential_calculate, alldays=alldays_mo, dNBR_Class=3)
    names(Percent_classH) <- paste0("Month_",mo) # each month one layer
    raster_highsevP <- c(raster_highsevP, Percent_classH)
  }
  raster_highsevP <- terra::rast(raster_highsevP)
  Pixel_summary_months <- terra::extract(raster_highsevP, Pixel_summary[,c("x","y")], xy=FALSE)[, names(raster_highsevP)]
  Pixel_summary_miyi <- Pixel_summary
  Pixel_summary_miyi$Year <- yi
  Pixel_summary_miyi[, c("Month_3", "Month_4", "Month_5", "Month_6", "Month_7", "Month_8", "Month_9", "Month_10", "Month_11")] <- NA
  Pixel_summary_miyi[, match(names(Pixel_summary_months), names(Pixel_summary_miyi))] <- Pixel_summary_months
  Pixel_potentials_mi <- rbind(Pixel_potentials_mi, Pixel_summary_miyi)
}
saveRDS(Pixel_potentials_mi, paste0("PredictClassiAlldays_",mi,"_40YrMonthlySeverePotential.rds"))




##################### Compute the most severe month & the most severe potential for each year #############################
SeveMonth_eachYear <- function(x) {
  if(all(is.na(x))) {  #all(is.na(x) | x == 0)
    MostSevereMonth <- NA # due to limited valid values for the potential calculation
  } else if(all(x == 0, na.rm = TRUE)) {
    MostSevereMonth <- 0 # due to all months' potential == 0 & there's no severe month in this year for this pixel
  } else {
    names_valid <- names(x)[!is.na(x)]
    MostSevereMonth <- names_valid[which.max(x[!is.na(x)])] # Otherwise, return the index of the max value (ignoring NA)
    MostSevereMonth <- strsplit(MostSevereMonth, "_")[[1]][2] ## only retain month number instead of "Month_"
  }
  return(MostSevereMonth)
}

SeveMonthPotential_eachYear <- function(x, monthtype) {
  monthi <- x[monthtype]
  if (is.na(monthi)){
    SeveMonthPotential <- NA
  } else if(monthi==0) {
    SeveMonthPotential <- 0
  } else {
    monthID <- grep(monthi, names(x), fixed=TRUE, value=FALSE)
    SeveMonthPotential <- as.numeric(x[monthID])
  }
  return(SeveMonthPotential)
}

HighSevePotent_em <- readRDS("PredictClassiAlldays_em_40YrMonthlySeverePotential.rds")
MonthColumns <- HighSevePotent_em[, paste0("Month_", seq(3,11,1))]
HighSevePotent_em$MostSevereMonth_annual <- as.numeric(apply(MonthColumns, 1, SeveMonth_eachYear))
HighSevePotent_em$MostSevereMonth_annualPotential <- as.numeric(apply(HighSevePotent_em, 1, SeveMonthPotential_eachYear, monthtype="MostSevereMonth_annual"))
HighSevePotent_em$MostSevereMonth_overall <- NA
for (i in unique(HighSevePotent_em$Pixel_ID)) {
  pixeli_rows <- which(HighSevePotent_em$Pixel_ID == i)
  data_pixeli <- HighSevePotent_em[pixeli_rows,]
  data_pixeli <- subset(data_pixeli, !is.na(MostSevereMonth_annual))
  if (nrow(data_pixeli)==0) {
    HighSevePotent_em$MostSevereMonth_overall[pixeli_rows] <- NA
    next
  }
  data_pixeli <- subset(data_pixeli, MostSevereMonth_annual!=0)
  if (nrow(data_pixeli)<5){
    HighSevePotent_em$MostSevereMonth_overall[pixeli_rows] <- 0
    next
  }
  majority_value <- data_pixeli %>% group_by(MostSevereMonth_annual) %>% summarise(freq = n(), avg_numeric = mean(MostSevereMonth_annualPotential)) %>%
    filter(freq == max(freq)) %>% arrange(desc(avg_numeric)) %>% slice(1) %>% pull(MostSevereMonth_annual)
  HighSevePotent_em$MostSevereMonth_overall[pixeli_rows] <- majority_value
}
HighSevePotent_em$MostSevereMonth_overallPotential <- as.numeric(apply(HighSevePotent_em, 1, SeveMonthPotential_eachYear, monthtype="MostSevereMonth_overall"))
saveRDS(HighSevePotent_em, "PredictClassiAlldays_em_40YrMonthlySeverePotential_SevereMonthPotential.rds")

overall_em <- HighSevePotent_em %>% group_by(x, y, Pre_nfiLandCover, Ecozone, Pixel_ID, MostSevereMonth_overall) %>% summarise(OverallMeanPotential=mean(MostSevereMonth_overallPotential, na.rm=TRUE))
saveRDS(overall_em, "PredictClassiAlldays_em_40YrMonthlySeverePotential_OverallSevereMonth&MeanPotentialSummary.rds")
overallMonth_em_raster <- terra::rasterize(x=as.matrix(overall_em[, c("x","y")]), y=VegeType_mask, values=as.numeric(overall_em$MostSevereMonth_overall))
names(overallMonth_em_raster) <- "MostSevereMonth_overall"
writeRaster(overallMonth_em_raster, "PredictClassiAlldays_em_40YrMonthlySeverePotential_OverallSevereMonthSummary.tif", overwrite=TRUE)
overallPotential_em_raster <- terra::rasterize(x=as.matrix(overall_em[, c("x","y")]), y=VegeType_mask, values=as.numeric(overall_em$OverallMeanPotential))
names(overallPotential_em_raster) <- "MostSevereMonth_overallMeanPotential"
writeRaster(overallPotential_em_raster, "PredictClassiAlldays_em_40YrMonthlySeverePotential_OverallSevereMonthMeanPotentialSummary.tif", overwrite=TRUE)





########################### Trend analysis for annual most severe month potential ###################################
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

Years_40 <- seq(1981,2020,1)
Years_20a <- seq(1981,2000,1)
Years_20b <- seq(2001,2020,1)
HighSevePotent_em_subAnnual <- subset(HighSevePotent_em, select=c(x,y,Pre_nfiLandCover,Ecozone,Pixel_ID,Year,MostSevereMonth_annualPotential))
HighSevePotent_em_subAnnual$Year <- paste0("Year_", HighSevePotent_em_subAnnual$Year)
HighSevePotentAnnual_em_reshape <- spread(HighSevePotent_em_subAnnual, Year, MostSevereMonth_annualPotential)
HighSevePotentAnnual_em_reshape$Trend1981_2020 <- apply(HighSevePotentAnnual_em_reshape[, paste0("Year_", Years_40)], 1, MK_FDR)
HighSevePotentAnnual_em_reshape$Trend1981_2000 <- apply(HighSevePotentAnnual_em_reshape[, paste0("Year_", Years_20a)], 1, MK_FDR)
HighSevePotentAnnual_em_reshape$Trend2001_2020 <- apply(HighSevePotentAnnual_em_reshape[, paste0("Year_", Years_20b)], 1, MK_FDR)
saveRDS(HighSevePotentAnnual_em_reshape, paste0("PredictClassiAlldays_em_40YrMonthlySeverePotential_AnnualSevereMonthPotential_Trends.rds"))

data_mi <- HighSevePotentAnnual_em_reshape
raster_pmi_Trend1981_2020 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend1981_2020))
names(raster_pmi_Trend1981_2020) <- "Trend1981_2020_AnnualPotential"
writeRaster(raster_pmi_Trend1981_2020, paste0("PredictClassiAlldays_em_40YrMonthlySeverePotential_AnnualSevereMonthPotential_Trend1981-2020.tif"), overwrite=TRUE)
raster_pmi_Trend1981_2000 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend1981_2000))
names(raster_pmi_Trend1981_2000) <- "Trend1981_2000_AnnualPotential"
writeRaster(raster_pmi_Trend1981_2000, paste0("PredictClassiAlldays_em_40YrMonthlySeverePotential_AnnualSevereMonthPotential_Trend1981-2000.tif"), overwrite=TRUE)
raster_pmi_Trend2001_2020 <- terra::rasterize(x=as.matrix(data_mi[, c("x","y")]), y=VegeType_mask, values=as.numeric(data_mi$Trend2001_2020))
names(raster_pmi_Trend2001_2020) <- "Trend2001_2020_AnnualPotential"
writeRaster(raster_pmi_Trend2001_2020, paste0("PredictClassiAlldays_em_40YrMonthlySeverePotential_AnnualSevereMonthPotential_Trend2001-2020.tif"), overwrite=TRUE)





########################### t-test for comparing the two 20-year periods ##########################
Years_40 <- seq(1981,2020,1)
Years_20a <- seq(1981,2000,1)
Years_20b <- seq(2001,2020,1)
mi <- "em"
HighSevePotentAnnual_em_reshape <- readRDS(paste0("Predictions_Classification/Predictions_Classification_SevereMonth/PredictClassiAlldays_",mi,"_40YrMonthlySeverePotential_AnnualSevereMonthPotential_Trends.rds"))
Ttest_mostsevere <- data.frame()
for (zi in unique(HighSevePotentAnnual_em_reshape$Ecozone)) {
  data_zi <- subset(HighSevePotentAnnual_em_reshape, Ecozone==zi)
  data_zi$MeanPotential_1981_2020 <- apply(data_zi[, paste0("Year_", Years_40)], 1, function(x) mean(x, na.rm=TRUE))  #mean(as.vector(as.matrix(data_zi[, paste0("Year_", Years_40)])), na.rm=TRUE)
  data_zi$MeanPotential_1981_2000 <- apply(data_zi[, paste0("Year_", Years_20a)], 1, function(x) mean(x, na.rm=TRUE)) 
  data_zi$MeanPotential_2001_2020 <- apply(data_zi[, paste0("Year_", Years_20b)], 1, function(x) mean(x, na.rm=TRUE))
  Ttest_Paired <- t.test(data_zi$MeanPotential_1981_2000, data_zi$MeanPotential_2001_2020, paired=TRUE)
  Ttest_mostsevere <- rbind(Ttest_mostsevere, data.frame(Ecozone=zi, MeanPotential_1981_2020=mean(data_zi$MeanPotential_1981_2020, na.rm=TRUE), 
                                                         MeanPotential_1981_2000 = mean(data_zi$MeanPotential_1981_2000, na.rm=TRUE), 
                                                         MeanPotential_2001_2020 = mean(data_zi$MeanPotential_2001_2020, na.rm=TRUE),
                                                         Difference_0120_8100 = mean(data_zi$MeanPotential_2001_2020 - data_zi$MeanPotential_1981_2000, na.rm=TRUE), #== mean(data_zi$MeanPotential_2001_2020, na.rm=TRUE) - mean(data_zi$MeanPotential_1981_2000, na.rm=TRUE),
                                                         Ttest_p_Mean20yrs = as.numeric(Ttest_Paired$p.value)))
}
write.table(Ttest_mostsevere, file=paste0("PredictClassiAlldays_",mi,"_40YrMonthlySeverePotential_AnnualSevereMonthPotential_MeanTtest_Summary.csv"), sep=",", append=T, col.names=T, row.names=F)




