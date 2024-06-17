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
Xs_factor <- c("Pre_nfiLandCover","DEM30_Aspect") 
Xs_Elev <- c("Elevation30", "DEM30_Slop", "DEM30_TPI") 
Xs_Fuel <- c("Pre_Biomass", "Pre_Closure", "Pre_Height", "Pre_prcC")  
Xs_FWeather <- c("BUI", "DMC", "FFMC", "FWI") 
Xs_numeric_All <- c(Xs_Elev, Xs_Fuel, Xs_FWeather) 
Xs_All <- c(Xs_factor, Xs_numeric_All)
Categorization_Methods <- c("nb", "em")
Jenks <- read.csv("DataPreparation/Jenks/ValidSamples(Seasons3_Trees3_dNBR271300)_ZonesFactors_Jenks_Ecozones.csv", head=T)
TrainSets <- seq(1,30,1)
UniqZones <- c(4, 9, 11, 12, 14, 15, 16, 17, 18, 19)

folder_ZonesYears <- "DataPreparation/ValidSamples_Years_Ecozones" # Full data for each zone (1 file per year)
files_ZonesYears <- list.files(folder_ZonesYears, pattern="*.rds$", full.names = TRUE)
folder_ZonesFactors <- "DataPreparation/ValidSamples_Years_Ecozones_Factors_MajorClass" # Only factors BUT 1 file for each zone
files_ZonesFactors <- list.files(folder_ZonesFactors, pattern="*.rds$", full.names = TRUE)
Variable_Selection <- data.frame()
Variable_Importance <- data.frame()
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
    models <- list()
    for (ti in TrainSets) {
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
      nzv <- caret::nearZeroVar(Samples_zi[,Xs_numeric], freqCut = 95/5, uniqueCut = 10)
      if(length(nzv)>0) {Xs_numeric <- Xs_numeric[-nzv]}
      nzv_results <- data.frame(matrix(1, ncol = length(Xs_All), nrow = 1)); colnames(nzv_results) <- Xs_All; nzv_results[match(Xs_numeric[nzv], Xs_All)] <- NA
      Variable_Selection <- rbind(Variable_Selection, data.frame(Ecozone=zi, Categorization_Method=mi,  TrainSet=ti, Process="near0var_NumericXs", nzv_results))
      
      cor <- caret::findCorrelation(cor(Samples_zi[,Xs_numeric]), cutoff = 0.7, names = TRUE)
      if(length(cor)>0) {Xs_numeric <- setdiff(Xs_numeric, cor)}
      cor_results <- data.frame(matrix(1, ncol = length(Xs_All), nrow = 1)); colnames(cor_results) <- Xs_All; cor_results[match(cor, Xs_All)] <- NA
      Variable_Selection <- rbind(Variable_Selection, data.frame(Ecozone=zi, Categorization_Method=mi,  TrainSet=ti, Process="HighCorr_NumericXs", cor_results))
      
      Xs_numeric <- unique(c(Xs_numeric, "FWI", "BUI", "Pre_prcC", "Pre_Biomass"))
      Samples_zi[, Xs_numeric] <- scale(Samples_zi[, Xs_numeric], center=TRUE, scale=TRUE)
      
      ##### Partition into calibration (70%) and validation (30%) sets
      train_size <- floor(0.7 * nrow(Samples_zi))
      train_indices <- sample(1:nrow(Samples_zi), train_size) 
      train_set <- Samples_zi[train_indices, ]
      test_set <- Samples_zi[-train_indices, ]
      
      ###### backward selection from all variables ######
      Xs_selected <- c(Xs_factor, Xs_numeric)
      Sele_results <- data.frame(matrix(NA, ncol = length(Xs_All), nrow = 1)); colnames(Sele_results) <- Xs_All; Sele_results[match(Xs_selected, Xs_All)] <- 1
      Variable_Selection <- rbind(Variable_Selection, data.frame(Ecozone=zi, Categorization_Method=mi,  TrainSet=ti, Process="Before_backward", Sele_results))
      formula_old <- as.formula(paste("dNBRclass ~", paste(Xs_selected, collapse="+")))
      model_mlr_old <- train(formula_old, data = train_set, method = "multinom") 
      pred_test_old <- predict(model_mlr_old, test_set, type = "raw") 
      ConMat_old <- confusionMatrix(pred_test_old, test_set$dNBRclass)
      Accuracy_old <- mean(data.frame(ConMat_old$byClass)[, "Balanced.Accuracy"])
      while (length(Xs_selected)>3) { ## at least remain 3 variables
        Accuracy_delete <- c()
        for (xi in Xs_selected) {
          Xs_delete <- setdiff(Xs_selected, xi)
          formula_delete <- as.formula(paste("dNBRclass ~", paste(Xs_delete, collapse="+")))
          model_mlr_delete <- train(formula_delete, data = train_set, method = "multinom") 
          pred_test_delete <- predict(model_mlr_delete, test_set, type = "raw")
          ConMat_delete <- confusionMatrix(pred_test_delete, test_set$dNBRclass)
          AccuClasses_delete <- mean(data.frame(ConMat_delete$byClass)[, "Balanced.Accuracy"])
          Accuracy_delete <- c(Accuracy_delete, AccuClasses_delete)
        }
        xi_ID <- which.max(Accuracy_delete)
        Accuracy_new <- Accuracy_delete[xi_ID]
        Xs_remove <- Xs_selected[xi_ID]
        if (Accuracy_new > Accuracy_old) {
          Xs_selected <- setdiff(Xs_selected, Xs_remove)
          Accuracy_old <- Accuracy_new
        } else {
          break
        }
      }
      Back_results <- data.frame(matrix(NA, ncol = length(Xs_All), nrow = 1)); colnames(Back_results) <- Xs_All; Back_results[match(Xs_selected, Xs_All)] <- 1
      Variable_Selection <- rbind(Variable_Selection, data.frame(Ecozone=zi, Categorization_Method=mi,  TrainSet=ti, Process="After_backward", Back_results))
      
      ######## Outputs by selected variables ########
      Xs <- Xs_selected
      formula <- as.formula(paste("dNBRclass ~", paste(Xs, collapse="+")))
      model_mlr <- train(formula, data = train_set, method = "multinom")
      models <- append(models, list(model_mlr)) 
      pred_test <- predict(model_mlr, test_set, type = "raw") 
      ConMat <- confusionMatrix(pred_test, test_set$dNBRclass)
      Accuracy_overall <- rbind(Accuracy_overall, cbind(t(data.frame(ConMat$overall)), Ecozone=zi, Categorization_Method=mi, TrainSet=ti, Nsamples=Nsample_zi))
      Accuracy_byClass <- rbind(Accuracy_byClass, cbind(data.frame(ConMat$byClass), Class=rownames(ConMat$byClass), Ecozone=zi, Categorization_Method=mi, TrainSet=ti, Nsamples=Nsample_zi))
      
      varImp_zi <- data.frame(varImp(model_mlr, scale = TRUE)$importance)
      rescale_to_1_100 <- function(x) { 1 + 99 * ((x - min(x)) / (max(x) - min(x))) }
      varImp_zi <- apply(varImp_zi, 2, rescale_to_1_100) 
      VegeType_ID <- grep("Pre_nfiLandCover", rownames(varImp_zi))
      Aspect_ID <- grep("DEM30_Aspect", rownames(varImp_zi))
      if (length(c(VegeType_ID, Aspect_ID))>0) { 
        varImp_zi_agg <- c(ifelse(length(VegeType_ID) >0, mean(varImp_zi[VegeType_ID, ]), 0), ifelse(length(Aspect_ID)>0, mean(varImp_zi[Aspect_ID, ]), 0), varImp_zi[-c(VegeType_ID, Aspect_ID), ])
        names(varImp_zi_agg) <- c("Pre_nfiLandCover","DEM30_Aspect", rownames(varImp_zi)[-c(VegeType_ID, Aspect_ID)]) 
      } else {
        varImp_zi_agg <- data.frame(t(varImp_zi))
      }
      varImportance_zi <- as.data.frame(matrix(0, nrow = 1, ncol = length(Xs_All)))
      colnames(varImportance_zi) <- Xs_All
      varImportance_zi[match(names(varImp_zi_agg), Xs_All)] <- varImp_zi_agg
      varImportance_zi$Ecozone <- zi; varImportance_zi$Categorization_Method <- mi; varImportance_zi$TrainSet <- ti
      Variable_Importance <- rbind(Variable_Importance, varImportance_zi)
    }
    saveRDS(models, paste("Model_Classification_MajorClassSamples_Zone",zi,"_", mi,"_Set", min(TrainSets), "-", max(TrainSets),".rds", sep="")) #Set1-10
  }
} 

write.table(Variable_Selection, file="ClassificationModel_VariableSelectionSummary.csv", sep=",", append=T, col.names=T, row.names=F)
write.table(Variable_Importance, file="ClassificationModel_VariableImportanceSummary.csv", sep=",", append=T, col.names=T, row.names=F)
write.table(Accuracy_overall, file="ClassificationModel_AccuracyOverallSummary.csv", sep=",", append=T, col.names=T, row.names=F)
write.table(Accuracy_byClass, file="ClassificationModel_AccuracyByClassSummary.csv", sep=",", append=T, col.names=T, row.names=F)

