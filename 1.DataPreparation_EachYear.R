library(dplyr)
library(terra)
library(automap)
library(sp)
library(gstat)
library(caret)



FW_interpolate <- function(fws, DEM30, firerast_fi_resa, di, ValidValues_fi_di) {
  proj_crs_string <- crs(DEM30, proj=TRUE) 
  firepoly_fi <- as.polygons(firerast_fi_resa, values = T, dissolve = T)
  names(firepoly_fi) <- "GrowthDay"
  firepoly_di <- firepoly_fi[firepoly_fi$GrowthDay==di]
  firepoly_di_buff <- buffer(firepoly_di, width = 50000) 
  fwi <- fws[[di]]
  fwi_crop <- terra::crop(fwi, firepoly_di_buff) 
  fwi_di <- terra::mask(fwi_crop, firepoly_di_buff)
  ValidCell_fwi <- terra::cells(fwi_di)
  ValidValues_fwi <- terra::extract(fwi_di, ValidCell_fwi, xy=TRUE)
  names(ValidValues_fwi)[which(names(ValidValues_fwi)==names(fwi_di))] <- "fweath"
  DemValues_fwi <- terra::extract(DEM30, ValidValues_fwi[,c("x","y")], xy=TRUE)
  names(DemValues_fwi)[which(names(DemValues_fwi)==names(DEM30))] <- "elev" 
  data_fit <- data.frame(fweath=ValidValues_fwi$fweath, elev=DemValues_fwi$elev)
  ID_valid <- which(complete.cases(data_fit))
  data_fit_valid <- data_fit[ID_valid, ]
  if(nrow(data_fit_valid)<30) {
    highres_fwi <- NA
  } else {
    low.res <- SpatialPointsDataFrame(coords = ValidValues_fwi[ID_valid, c("x","y")], 
                                      data = data_fit_valid, 
                                      proj4string = sp::CRS(proj_crs_string))
    high.res <- SpatialPointsDataFrame(coords = ValidValues_fi_di[,c("x","y")], 
                                       data = data.frame(elev=ValidValues_fi_di$Elevation),
                                       proj4string = sp::CRS(proj_crs_string))
    variogram <- automap::autofitVariogram(fweath ~ elev, input_data = low.res, model = "Sph")
    krige <- gstat::krige(fweath~elev, low.res, high.res, model = variogram$var_model)
    highres_fwi <- krige@data[,1]
  }
  return(highres_fwi)
}




########################################## Main ###################################################
### All other data except fire_growth are all projected by crs(DEM)=(ESRI:102002)
dNBR1000_proj <- terra::rast("Data/CanlaBS2/CanlaBS_dNBR1000_v1_Sal0Gap1Nbf1(ESRI_102002-MethodNear).tif")
dNBRYear_proj <- terra::rast("Data/CanlaBS2/CanLaBS_FireYear_v1(ESRI_102002).tif")
Ecozones_proj_rast <- terra::rast("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones(ESRI_102002)(DEM30).tif")
Elefolder <- "Data/Fire_Growth/DEM_process/"
DEM30 <- terra::rast(paste0(Elefolder,"EarthEnv_dem100m_buffer_Zones_resampledNBR1000.tif"))
DEM30_slop <- terra::rast(paste0(Elefolder,"EarthEnv_dem100m_buffer_Zones_resampledNBR1000_Slope.tif"))
DEM30_aspect_8direction <- terra::rast(paste0(Elefolder,"EarthEnv_dem100m_buffer_Zones_resampledNBR1000_Aspect8Direction.tif"))
DEM30_TPI <- terra::rast(paste0(Elefolder,"EarthEnv_dem100m_buffer_Zones_resampledNBR1000_TPI.tif"))
Pre_Forest_proj_folder <- "Data/Forest_SCANFI/Prefire_ForestAttri/Project_DEM100/"
Pre_Biomass_proj <- terra::rast(paste(Pre_Forest_proj_folder,"SCANFI_att_biomass_S_prefire_v0(ESRI_102002).tif", sep=""))
Pre_Closure_proj <- terra::rast(paste(Pre_Forest_proj_folder,"SCANFI_att_closure_S_prefire_v0(ESRI_102002).tif", sep=""))
Pre_Height_proj <- terra::rast(paste(Pre_Forest_proj_folder,"SCANFI_att_height_S_prefire_v0(ESRI_102002).tif", sep=""))
Pre_nfiLandCover_proj <- terra::rast(paste(Pre_Forest_proj_folder,"SCANFI_att_nfiLandCover_S_prefire_v0(ESRI_102002).tif", sep=""))
Pre_prcC_All_proj <- terra::rast(paste(Pre_Forest_proj_folder,"SCANFI_sps_prcC_ALL_S_prefire_v0(ESRI_102002).tif", sep=""))
proj_crs <- crs(DEM30)
FWI_inputs_folder <- "Data/ERA5_Weather/FWI_1979_2021/Inputs_CA(ESRI_102002)"
FWI_outputs_folder <- "Data/ERA5_Weather/FWI_1979_2021/Outputs_CA(ESRI_102002)"
years <- seq(2001, 2020, 1)
for (yi in years) {
  print(paste("Year:",yi,"/ Whole range:", min(years), "-", max(years)))
  Samples_outputs <- data.frame()
  folder_yi <- paste("Data/Fire_Growth/fire_growth_nbac_asc/", yi, sep="")
  files_yi <- list.files(folder_yi, pattern="*.asc$", full.names = TRUE)
  firepoly_yi <- terra::vect(paste("Data/Fire_Growth/fire_growth_nbac_poly/nbac_",yi,".shp", sep=""))
  
  BUI_yi <- terra::rast(paste(FWI_outputs_folder,"/build_up_index_",yi,".tif", sep=""))
  DMC_yi <- terra::rast(paste(FWI_outputs_folder,"/duff_moisture_code_",yi,".tif", sep=""))
  FFMC_yi <- terra::rast(paste(FWI_outputs_folder,"/fine_fuel_moisture_code_",yi,".tif", sep=""))
  FWI_yi <- terra::rast(paste(FWI_outputs_folder,"/fire_weather_index_",yi,".tif", sep=""))
  for (fi in files_yi) {
    print(fi)
    fid <- strsplit(fi, "_")[[1]]
    fid <- fid[length(fid)-1]
    firerast_fi <- terra::rast(fi)
    crs(firerast_fi) <- crs(firepoly_yi)
    firerast_fi <- terra::project(firerast_fi, proj_crs, method="near")
    crop_extent <- terra::ext(firerast_fi)
    
    DEM30_fi <- terra::crop(DEM30, crop_extent)
    firerast_fi_resa <-  terra::resample(firerast_fi, DEM30_fi, method="near")
    dNBR1000_fi <- terra::crop(dNBR1000_proj, crop_extent)
    dNBRYear_fi <- terra::crop(dNBRYear_proj, crop_extent)
    ecozone_fi <- terra::crop(Ecozones_proj_rast, crop_extent)
    DEM30_slop_fi <- terra::crop(DEM30_slop, crop_extent)
    DEM30_aspect_fi <- terra::crop(DEM30_aspect_8direction, crop_extent)
    DEM30_TPI_fi <- terra::crop(DEM30_TPI, crop_extent)
    Pre_Biomass_fi <- terra::crop(Pre_Biomass_proj, crop_extent)
    Pre_Closure_fi <- terra::crop(Pre_Closure_proj, crop_extent)
    Pre_Height_fi <- terra::crop(Pre_Height_proj, crop_extent)
    Pre_nfiLandCover_fi <- terra::crop(Pre_nfiLandCover_proj, crop_extent)
    Pre_prcC_All_fi <- terra::crop(Pre_prcC_All_proj, crop_extent)
    
    raster_stack_fi <- c(dNBRYear_fi, firerast_fi_resa, dNBR1000_fi, ecozone_fi, 
                         DEM30_fi, DEM30_slop_fi, DEM30_aspect_fi, DEM30_TPI_fi, 
                         Pre_nfiLandCover_fi, Pre_Biomass_fi, Pre_Closure_fi, Pre_Height_fi, Pre_prcC_All_fi) 
    ValidCell_fi <- terra::cells(firerast_fi_resa)
    ValidValues_fi <- terra::extract(raster_stack_fi, ValidCell_fi, xy=TRUE)
    ValidValues_fi <- na.omit(ValidValues_fi)
    names(ValidValues_fi) <- c("x", "y", "FireYear_dNBR", "GrowthDay", "dNBR1000", "Ecozone", 
                               "Elevation30", "DEM30_Slop", "DEM30_Aspect", "DEM30_TPI", 
                               "Pre_nfiLandCover", "Pre_Biomass", "Pre_Closure", "Pre_Height", "Pre_prcC")
    ValidValues_fi <- subset(ValidValues_fi, FireYear_dNBR==yi)
    if(nrow(ValidValues_fi)==0) {
      next
    }
    ValidValues_fi$FireYear_growth <- yi
    ValidValues_fi$FireID <- fid
    jdays <- unique(ValidValues_fi$GrowthDay)
    ValidValues_fi[, c("BUI", "DMC","FFMC", "FWI")] <- NA
    for (di in jdays) {
      di_ID <- which(ValidValues_fi$GrowthDay==di)
      ValidValues_fi_di <- ValidValues_fi[di_ID,]
      ValidValues_fi$BUI[di_ID] <- FW_interpolate(BUI_yi, DEM30, firerast_fi_resa, di, ValidValues_fi_di)
      ValidValues_fi$DMC[di_ID] <- FW_interpolate(DMC_yi, DEM30, firerast_fi_resa, di, ValidValues_fi_di)
      ValidValues_fi$FFMC[di_ID] <- FW_interpolate(FFMC_yi, DEM30, firerast_fi_resa, di, ValidValues_fi_di)
      ValidValues_fi$FWI[di_ID]  <- FW_interpolate(FWI_yi, DEM30, firerast_fi_resa, di, ValidValues_fi_di)
    }
    ValidValues_fi <- na.omit(ValidValues_fi)
    Samples_outputs <- rbind(Samples_outputs, ValidValues_fi)
  }
  saveRDS(Samples_outputs, paste("DataPreparation_ValidSamples_Year",yi,".rds", sep=""))
}





