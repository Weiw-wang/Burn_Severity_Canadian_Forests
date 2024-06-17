library(dplyr)
library(terra)
library(tidyr)
library(EnvStats)


############################### daily estimates during 1981-2020 for all ecozones #######################################################
Elevation025 <- terra::rast("Data/Fire_Growth/DEM_process/EarthEnv_dem100m_buffer_resampleFWI2020.tif")
DEM025_Aspect <- terra::rast("Data/Fire_Growth/DEM_process/EarthEnv_dem100m_buffer_resampleFWI2020_Aspect8Direction.tif")
DEM025_Slop <- terra::rast("Data/Fire_Growth/DEM_process/EarthEnv_dem100m_buffer_resampleFWI2020_Slope.tif")
DEM025_TPI <- terra::rast("Data/Fire_Growth/DEM_process/EarthEnv_dem100m_buffer_resampleFWI2020_TPI.tif")
Pre_nfiLandCover <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_att_nfiLandCover_SW_2020_v1_MODEresampleFWI2020(ESRI_102002)_3Trees.tif")
Pre_Biomass <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_att_biomass_SW_2020_v1_resampleFWI2020(ESRI_102002).tif")
Pre_Closure <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_att_closure_SW_2020_v1_resampleFWI2020(ESRI_102002).tif")
Pre_Height <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_att_height_SW_2020_v1_resampleFWI2020(ESRI_102002).tif")
Pre_prcC <- terra::rast("Data/Forest_SCANFI/National_ForestAttri/Resample_FWI2020/SCANFI_sps_prcC_ALL_SW_2020_v1_resampleFWI2020(ESRI_102002).tif")
Ecozone <- terra::rast("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones(ESRI_102002)_resampleFWI2020.tif")
VegeType <- terra::resample(Pre_nfiLandCover, Ecozone, method="near")
crop_extent <- terra::ext(Ecozone)
VegeType_crop <- terra::crop(VegeType, crop_extent)
VegeType_mask <- terra::mask(VegeType_crop, Ecozone); names(VegeType_mask) <- "Pre_nfiLandCover"
ValidCells <- terra::cells(VegeType_mask) 

Samples_Static <- terra::extract(VegeType_mask, ValidCells, xy=TRUE)
Samples_Static$Elevation30 <- terra::extract(Elevation025, Samples_Static[,c("x","y")], xy=FALSE)[, names(Elevation025)]
Samples_Static$DEM30_Aspect <- terra::extract(DEM025_Aspect, Samples_Static[,c("x","y")], xy=FALSE)[, names(DEM025_Aspect)]
Samples_Static$DEM30_Slop <- terra::extract(DEM025_Slop, Samples_Static[,c("x","y")], xy=FALSE)[, names(DEM025_Slop)]
Samples_Static$DEM30_TPI <- terra::extract(DEM025_TPI, Samples_Static[,c("x","y")], xy=FALSE)[, names(DEM025_TPI)]
Samples_Static$Pre_Biomass <- terra::extract(Pre_Biomass, Samples_Static[,c("x","y")], xy=FALSE)[, names(Pre_Biomass)]
Samples_Static$Pre_Closure <- terra::extract(Pre_Closure, Samples_Static[,c("x","y")], xy=FALSE)[, names(Pre_Closure)]
Samples_Static$Pre_Height <- terra::extract(Pre_Height, Samples_Static[,c("x","y")], xy=FALSE)[, names(Pre_Height)]
Samples_Static$Pre_prcC <- terra::extract(Pre_prcC, Samples_Static[,c("x","y")], xy=FALSE)[, names(Pre_prcC)]
Samples_Static$Ecozone <- terra::extract(Ecozone, Samples_Static[,c("x","y")], xy=FALSE)[, names(Ecozone)]

Xs_Elev <- c("Elevation30", "DEM30_Slop", "DEM30_TPI") 
Xs_Fuel <- c("Pre_Biomass", "Pre_Closure", "Pre_Height", "Pre_prcC")
Xs_FWeather <- c("BUI", "DMC", "FFMC", "FWI")  
Xs_numeric_All <- c(Xs_Elev, Xs_Fuel, Xs_FWeather) 
years <- seq(1981,2020,1)
Categorization_Methods <- c("em")
TrainSets <- seq(1,30,1)
for (mi in Categorization_Methods) {
  for (yi in years) {
    print(paste("Category:",mi, "; Year:",yi))
    FWI_outputs_folder <- "Data/ERA5_Weather/FWI_1979_2021/Outputs_CA(ESRI_102002)"
    BUI_yi <- terra::rast(paste(FWI_outputs_folder,"/build_up_index_",yi,".tif", sep=""))
    DMC_yi <- terra::rast(paste(FWI_outputs_folder,"/duff_moisture_code_",yi,".tif", sep=""))
    FFMC_yi <- terra::rast(paste(FWI_outputs_folder,"/fine_fuel_moisture_code_",yi,".tif", sep=""))
    FWI_yi <- terra::rast(paste(FWI_outputs_folder,"/fire_weather_index_",yi,".tif", sep=""))
    date_Spr <- as.Date(paste0(yi,"-03-01")) #spring (March-May)
    Jday_Spr <- as.numeric(date_Spr - as.Date(paste0(yi,"-01-01"))) + 1 
    date_Win <- as.Date(paste0(yi,"-12-01")) #winter (December-February)
    Jday_Win <- as.numeric(date_Win - as.Date(paste0(yi,"-01-01"))) + 1
    days <- seq(Jday_Spr, Jday_Win-1, 1)
    Year_maps <- c()
    for (di in days) {
      BUI <- BUI_yi[[di]]
      DMC <- DMC_yi[[di]]
      FFMC <- FFMC_yi[[di]]
      FWI <- FWI_yi[[di]]
      Samples_All <- Samples_Static
      Samples_All$BUI <- terra::extract(BUI, Samples_All[,c("x","y")], xy=FALSE)[, names(BUI)]
      Samples_All$DMC <- terra::extract(DMC, Samples_All[,c("x","y")], xy=FALSE)[, names(DMC)]
      Samples_All$FFMC <- terra::extract(FFMC, Samples_All[,c("x","y")], xy=FALSE)[, names(FFMC)]
      Samples_All$FWI <- terra::extract(FWI, Samples_All[,c("x","y")], xy=FALSE)[, names(FWI)]
      Samples_All <- na.omit(Samples_All)
      if(nrow(Samples_All)==0) {
        next
      }
      Samples_All$Pre_nfiLandCover <- as.factor(Samples_All$Pre_nfiLandCover)
      Samples_All$Ecozone <- as.factor(Samples_All$Ecozone)
      Samples_All$DEM30_Aspect <- as.factor(Samples_All$DEM30_Aspect)
      
      Predictions_si <- data.frame()
      Unq_zones <- unique(Samples_All$Ecozone)
      for (zi in Unq_zones) {
        Samples_zi <- subset(Samples_All, Ecozone==zi)
        Samples_zi[, Xs_numeric_All] <- scale(Samples_zi[, Xs_numeric_All], center=TRUE, scale=TRUE)
        Samples_zi <- na.omit(Samples_zi) 
        if(nrow(Samples_zi)==0) {
          next
        }
        models <- readRDS(paste("ModelBuilding_Classification/multinom_Models/Model_Classification_MajorClassSamples_Zone",zi,"_",mi,"_Set", min(TrainSets), "-", max(TrainSets),".rds", sep=""))
        predicts <- lapply(models, function(x) predict(x, newdata = Samples_zi, type = "raw"))
        predicts <- data.frame(predicts); names(predicts) <- paste("model",seq(1:ncol(predicts)),sep="")
        getmode <- function(v) {
          uniqv <- sort(unique(v), decreasing = TRUE)
          uniqv[which.max(tabulate(match(v, uniqv)))]
        }
        MajorVotes <- as.factor(apply(predicts, 1, getmode))
        Samples_zi$Prediction <- MajorVotes
        Predictions_si <- rbind(Predictions_si, Samples_zi)
      }
      if (nrow(Predictions_si)==0) {
        next
      }
      Predictions_si_raster <- terra::rasterize(x=as.matrix(Predictions_si[, c("x","y")]), y=VegeType_mask, values=as.numeric(Predictions_si$Prediction))
      names(Predictions_si_raster) <- di
      Year_maps <- c(Year_maps, Predictions_si_raster)
    }
    if(length(Year_maps)>0) {
      Year_maps <- terra::rast(Year_maps)
      writeRaster(Year_maps, paste("Predictions_", mi, "_Year", yi, "_March-November.tif", sep=""), overwrite=TRUE)
    }
  }
}




