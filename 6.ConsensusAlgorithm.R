library(dplyr)
library(caret)
library(rpart)
library(corrplot)
library(tidyr)
library(nnet)



############ Stratified Subsampling ##################
Startified_sampling <- function(Nsamples, AllSamples){
  stra_size <- floor( Nsamples/(length(unique(AllSamples$Pre_nfiLandCover)) * length(unique(AllSamples$DEM30_Aspect)) * length(unique(AllSamples$dNBRclass))) )
  Subsamples <- AllSamples %>%
    group_by(Pre_nfiLandCover, DEM30_Aspect, dNBRclass) %>% 
    sample_n(size=stra_size, replace=TRUE)
  Subsamples <- unique(Subsamples)
  
  if (nrow(Subsamples)<Nsamples) {
    add_size <- Nsamples - nrow(Subsamples)
    add_ID <- sample(AllSamples$ID[-Subsamples$ID], add_size, replace = FALSE)
    Subsamples <- rbind(Subsamples, AllSamples[add_ID,])
  } else if (nrow(Subsamples)>Nsamples){
    minus_size <- nrow(Subsamples)- Nsamples
    minus_ID <- sample(1:nrow(Subsamples), minus_size, replace = FALSE)
    Subsamples <- Subsamples[-minus_ID,]
  }
  return(Subsamples)
}


################## Categorize dNBR1000 ############## 
Categorize_dNBR1000 <- function(AllSamples, method_interest, nb_thresholds) {
  if (method_interest=="nb") { ### Categorize with natural breaks (Jenks)
    AllSamples$dNBRclass <- cut(AllSamples$dNBR1000,
                                breaks = c(26, round(nb_thresholds[1]), round(nb_thresholds[2]), 1301),
                                labels = c(1,2,3),
                                include.lowest = FALSE, right = TRUE)
  } else if (method_interest=="em") { ### Categorize with empirical breaks
    AllSamples$dNBRclass <- cut(AllSamples$dNBR1000,
                                breaks = c(26,276,517,1301),
                                labels = c(1,2,3), 
                                include.lowest = FALSE, right = TRUE)
  }
  return(AllSamples)
}



############################################# MAIN #################################################
Xs_factor <- c("Pre_nfiLandCover",  "DEM30_Aspect")
Xs_Elev <- c("Elevation30", "DEM30_Slop", "DEM30_TPI") 
Xs_Fuel <- c("Pre_Biomass", "Pre_Closure", "Pre_Height", "Pre_prcC") 
Xs_FWeather <- c("BUI", "DMC", "FFMC", "FWI") 
Xs_numeric_All <- c(Xs_Elev, Xs_Fuel, Xs_FWeather) 
Xs_All <- c(Xs_factor, Xs_numeric_All)
Categorization_Methods <- c("nb", "em")
Jenks <- read.csv("DataPreparation/Jenks/ValidSamples(Seasons3_Trees3_dNBR271300)_ZonesFactors_Jenks_Ecozones.csv", head=T)
ConSets <- seq(1,30,1) 
TrainSets <- seq(1,30,1) 
UniqZones <- c(4, 9, 11, 12, 14, 15, 16, 17, 18, 19)

folder_ZonesYears <- "DataPreparation/ValidSamples_Years_Ecozones"
files_ZonesYears <- list.files(folder_ZonesYears, pattern="*.rds$", full.names = TRUE)
folder_ZonesFactors <- "DataPreparation/ValidSamples_Years_Ecozones_Factors_MajorClass" 
files_ZonesFactors <- list.files(folder_ZonesFactors, pattern="*.rds$", full.names = TRUE)
Accuracy_overall <- data.frame()
Accuracy_byClass <- data.frame()
for (zi in UniqZones) {
  print(zi)
  Nsample_zi <- 5000 
  Jenks_zi <- subset(Jenks, Ecozone==zi & Repeats=="Average")
  nb_thresholds_zi <- Jenks_zi[c(grep("Low", names(Jenks_zi)), grep("Moderate", names(Jenks_zi)))] 
  for (mi in Categorization_Methods) {
    fi_zfactor <- grep(paste0(mi,"_Zone",zi), files_ZonesFactors)
    SamplesFactor <- readRDS(files_ZonesFactors[fi_zfactor])
    SamplesFactor$ID <- seq(1, nrow(SamplesFactor),1) 
    SamplesFactor$Pre_nfiLandCover <- as.factor(SamplesFactor$Pre_nfiLandCover)
    SamplesFactor$DEM30_Aspect <- as.factor(SamplesFactor$DEM30_Aspect)
    SamplesFactor <- Categorize_dNBR1000(SamplesFactor, mi, nb_thresholds_zi)
    for (ti in ConSets) {
      print(paste("Ecozone", zi, "; Category", mi, "; Set", ti))
      Samples_zi_fact <- Startified_sampling(Nsample_zi, SamplesFactor)
      Samples_zi <- c()
      for (yi in unique(Samples_zi_fact$FireYear_dNBR)) {
        file_id <- grep(paste0("Zone",zi,"_Year",yi), files_ZonesYears)
        FullSet <- readRDS(files_ZonesYears[file_id])
        yid <- which(Samples_zi_fact$FireYear_dNBR==yi)
        FullSet <- FullSet[Samples_zi_fact$YearRow[yid], Xs_numeric_All]
        Samples_zi_yi <- cbind(Samples_zi_fact[yid, ], FullSet)
        Samples_zi <- rbind(Samples_zi, Samples_zi_yi)
      }
      Xs_numeric <- Xs_numeric_All
      Samples_zi[, Xs_numeric] <- scale(Samples_zi[, Xs_numeric], center=TRUE, scale=TRUE)
      models <- readRDS(paste("ModelBuilding_Classification/multinom_Models/Model_Classification_MajorClassSamples_Zone",zi,"_",mi,"_Set", min(TrainSets), "-", max(TrainSets),".rds", sep=""))
      cons_set <- Samples_zi 
      predictions <- lapply(models, function(x) predict(x, newdata = cons_set, type = "raw"))
      predictions <- data.frame(predictions); names(predictions) <- paste("model",seq(1:ncol(predictions)),sep="")
      getmode <- function(v) {
        uniqv <- sort(unique(v), decreasing = TRUE)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      MajorVotes <- as.factor(apply(predictions, 1, getmode))
      ConMat <- confusionMatrix(MajorVotes, cons_set$dNBRclass)
      Accuracy_overall <- rbind(Accuracy_overall, cbind(t(data.frame(ConMat$overall)), Ecozone=zi, Categorization_Method=mi, ConsenSet=ti, Nsamples=Nsample_zi))
      Accuracy_byClass <- rbind(Accuracy_byClass, cbind(data.frame(ConMat$byClass), Class=rownames(ConMat$byClass), Ecozone=zi, Categorization_Method=mi, ConsenSet=ti, Nsamples=Nsample_zi))
    }
  }
}
write.table(Accuracy_overall, file="ClassificationConsensusModel_AccuracyOverallSummary.csv", sep=",", append=T, col.names=T, row.names=F)
write.table(Accuracy_byClass, file="ClassificationConsensusModel_AccuracyByClassSummary.csv", sep=",", append=T, col.names=T, row.names=F)

