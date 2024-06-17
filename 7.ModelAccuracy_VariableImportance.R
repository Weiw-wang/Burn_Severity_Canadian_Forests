library(dplyr)
library(tidyr)
library(terra)
library(cowplot)
library(ggplotify)


############## Summarize byClass accuracy: averaging over TRAIN/TEST subsets ############
outfile <- "Classification(Consensus)Model_AccuracyByClassSummary_Ave4Plot.csv"
Accuracy_byClass <- read.csv("ModelBuilding_Classification/Summary/ClassificationModel_AccuracyByClassSummary.csv", head=T)
Ave_Accuracy_byClass <- data.frame()
for (zi in unique(Accuracy_byClass$Ecozone)) {
  for (mi in unique(Accuracy_byClass$Categorization_Method)) {
    for (ci in unique(Accuracy_byClass$Class)) {
      sub_results <- subset(Accuracy_byClass, Class==ci & Ecozone==zi & Categorization_Method==mi)
      col_means <- colMeans(subset(sub_results, select=-c(Class, Ecozone, Categorization_Method, TrainSet, Nsamples)), na.rm = TRUE)
      col_means <- cbind(t(col_means), Class=ci, Ecozone=zi, Categorization_Method=mi, TrainSet="Average")
      Ave_Accuracy_byClass <- rbind(Ave_Accuracy_byClass, col_means)
    }
  }
}
F1_summary <- spread(Ave_Accuracy_byClass[,c("F1","Class","Ecozone","Categorization_Method")], Class, F1)
F1_summary$Measure <- "F1"; F1_summary$Sets <- "Train7/Test3"
BaAc_summary <- spread(Ave_Accuracy_byClass[,c("Balanced.Accuracy","Class","Ecozone","Categorization_Method")], Class, Balanced.Accuracy)
BaAc_summary$Measure <- "Balanced.Accuracy"; BaAc_summary$Sets <- "Train7/Test3"
write.table(cbind(F1_summary, "", BaAc_summary), file=outfile, sep=",", append=T, col.names=T, row.names=F)


############## Summarize byClass accuracy: averaging over CONSENSUS subsets #############
AccuracyCon_byClass <- read.csv("ModelBuilding_Classification/Summary/ClassificationConsensusModel_AccuracyByClassSummary.csv", head=T)
Ave_AccuracyCon_byClass <- data.frame()
for (zi in unique(AccuracyCon_byClass$Ecozone)) {
  for (mi in unique(AccuracyCon_byClass$Categorization_Method)) {
    for (ci in unique(AccuracyCon_byClass$Class)) {
      sub_results <- subset(AccuracyCon_byClass, Class==ci & Ecozone==zi & Categorization_Method==mi)
      col_means <- colMeans(subset(sub_results, select=-c(Class, Ecozone, Categorization_Method, ConsenSet, Nsamples)), na.rm = TRUE)
      col_means <- cbind(t(col_means), Class=ci, Ecozone=zi, Categorization_Method=mi, ConsenSet="Average")
      Ave_AccuracyCon_byClass <- rbind(Ave_AccuracyCon_byClass, col_means)
    }
  }
}
F1_summary <- spread(Ave_AccuracyCon_byClass[,c("F1","Class","Ecozone","Categorization_Method")], Class, F1)
F1_summary$Measure <- "F1"; F1_summary$Sets <- "Consensus"
BaAc_summary <- spread(Ave_AccuracyCon_byClass[,c("Balanced.Accuracy","Class","Ecozone","Categorization_Method")], Class, Balanced.Accuracy)
BaAc_summary$Measure <- "Balanced.Accuracy"; BaAc_summary$Sets <- "Consensus"
write.table(cbind(F1_summary, "", BaAc_summary), file=outfile, sep=",", append=T, col.names=T, row.names=F)





####################### Summarize Results of variable importance ##########################
Variable_Importance <- read.csv("ModelBuilding_Classification/Summary/ClassificationModel_VariableImportanceSummary.csv", head=T)
Ave_Variable_Importance <- data.frame()
for (zi in unique(Variable_Importance$Ecozone)) {
  for (mi in unique(Variable_Importance$Categorization_Method)) {
    sub_results <- subset(Variable_Importance, Ecozone==zi & Categorization_Method==mi)
    col_means <- colMeans(subset(sub_results, select=-c(Ecozone, Categorization_Method, TrainSet)), na.rm = TRUE)
    col_means <- cbind(t(col_means), Ecozone=zi, Categorization_Method=mi, TrainSet="Average")
    Ave_Variable_Importance <- rbind(Ave_Variable_Importance, col_means)
  }
}
write.table(Ave_Variable_Importance, file="ClassificationModel_VariableImportanceSummary_Ave4Plot.csv", sep=",", append=T, col.names=T, row.names=F)

################## To plot variable importance ###################
Ave_Variable_Importance <- read.csv("ClassificationModel_VariableImportanceSummary_Ave4Plot.csv", head=T)
Ave_VarImp <- gather(Ave_Variable_Importance, Variable, Importance, Pre_nfiLandCover:FWI)
Xs_weather <- c("FFMC","DMC","BUI","FWI"); Xs_topo <- c("Elevation30", "DEM30_Aspect", "DEM30_Slop", "DEM30_TPI"); 
Xs_vege <- c("Pre_nfiLandCover","Pre_Biomass","Pre_Closure","Pre_Height","Pre_prcC")
Varitypes <- c("Weather",  "Topography", "Fuel")
Ave_VarImp$VariType[which(Ave_VarImp$Variable %in% Xs_weather)] <- Varitypes[1]
Ave_VarImp$VariType[which(Ave_VarImp$Variable %in% Xs_topo)] <- Varitypes[2]
Ave_VarImp$VariType[which(Ave_VarImp$Variable %in% Xs_vege)] <- Varitypes[3]
Ave_VarImp$VariType <- factor(Ave_VarImp$VariType, levels = Varitypes, labels = Varitypes)

Variable_IDs <- c(Xs_weather, Xs_topo, Xs_vege)
Variable_names <- c("FFMC","DMC","BUI","FWI", 
                    "ELEV","ASPECT", "SLOP", "TPI", 
                    "VEGETYPE","BIOM","CRCL","HEIGHT","CONIC") 
Ave_VarImp$Variable <- factor(Ave_VarImp$Variable, levels = Variable_IDs, labels = Variable_names)
Ecozone_IDs <- c("4", "9", "11", "12", "14", "15", "16", "17", "18", "19")
Ecozone_fullnames <- c("Taiga Plains (TP)", "Boreal Plains (BP)", "Taiga Cordillera (TC)", "Boreal Cordillera (BC)", "Montane Cordillera (MC)", 
                       "Hudson Plains (HP)", "Taiga Shield West (TSW)", "Taiga Shield East (TSE)", "Boreal Shield West (BSW)", "Boreal Shield East (BSE)")
Ave_VarImp$Ecozone_full <- factor(Ave_VarImp$Ecozone, levels = Ecozone_IDs, labels = Ecozone_fullnames)
Ecozone_abbnames <- c("TP", "BP", "TC", "BC", "MC", "HP", "TSW", "TSE", "BSW", "BSE")
Ave_VarImp$Ecozone_abb <- factor(Ave_VarImp$Ecozone, levels = Ecozone_IDs, labels = Ecozone_abbnames)

ci <- "em"
Ave_VarImp_ci <- subset(Ave_VarImp, Categorization_Method==ci)
Ave_VarImp_ci_summ0 <- Ave_VarImp_ci %>% group_by(Ecozone_full, Ecozone_abb, Categorization_Method, VariType) %>% summarise(Mean_Importance=mean(Importance, na.rm=TRUE))
Ave_VarImp_ci_summ1 <- Ave_VarImp_ci_summ0 %>% group_by(Ecozone_full, Ecozone_abb) %>% mutate(Percent_Importance = Mean_Importance*100/sum(Mean_Importance))

Ecozone_poly <- terra::vect("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones.shp")
Ecozone_rast <- terra::rast("Data/Ecozones_4Fire/ShpTiff/Ecozones4Fire_10Zones(ESRI_102002)_resampleFWI2020.tif")
Ecozone_trans <- terra::project(Ecozone_poly, terra::crs(Ecozone_rast))
centers_points <- terra::centroids(Ecozone_trans, inside=TRUE)
centers <- data.frame(centers_points, terra::crds(centers_points))
centers <- centers[-c(9),] 
extent_zone <- terra::ext(Ecozone_trans)

g1 <- ggplot(data = Ave_VarImp_ci, aes(x = as.factor(Variable), y = as.factor(Ecozone_full), fill = as.numeric(Importance))) +
  geom_tile() + 
  xlab("Variable") + ylab("Ecozone") + labs(title="(a)", fill = "Importance (%)") + 
  scale_fill_gradient(low = "#F2FDFB", high = "#04C6A9" , na.value = "grey") +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust=1, vjust=1),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0.5, vjust=0.3)) +
  facet_grid(. ~ VariType, scales = "free_x", space = "free_x")

g2 <- ggplot(Ave_VarImp_ci_summ1, aes(x = Percent_Importance, y = Ecozone_abb, fill = VariType)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Ecozone") +  xlab("Average importance of each group (%)") + labs(title="(b)", fill = "Variable group") +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust=0, vjust=1),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        strip.text = element_text(face = "bold", hjust = 0.5, vjust=0.3)) +
  geom_text(aes(label = round(Percent_Importance)), size=3,
            position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = c("Weather" = "#F3CDAA", "Topography" = "#F7F7F8", "Fuel" = "#A4D9B7"))

g3 <- ggplotify::as.ggplot(function(){
  #### Assign category by checking values in g2
  color_mapping <- c("BSE"="#A4D9B7", "BSW"="#E5F5EA", "TSE"="#F3CDAA", "TSW"="#E5F5EA",
                     "HP"="#FAECDE", "MC"="#A4D9B7", "BC"="#FAECDE",
                     "TC"="#F3CDAA", "BP"="#F3CDAA", "TP"="#E5F5EA")
  color_to_use <- color_mapping[as.character(Ecozone_trans$ZONE_NOM)]
  par(bg = "white")
  terra::plot(Ecozone_trans, mar=c(0,1,0,0), col=color_to_use, border="grey35", lwd=0.6, axes=FALSE) # 
  text(centers$x, centers$y, centers$ZONE_NOM, cex=.7, col="grey35")
  mtext("(c)", side=3, adj=0.9, line=-1, cex=1.2, font=2)
  legend(extent_zone[2]*0.35, extent_zone[4]*1.1,  title="Most influential group",  box.lwd=0,  title.cex=0.9, cex=0.7, y.intersp=1.6, x.intersp=2, text.width=c(380000, 380000, 380000), pch=15, pt.cex = 2,
         legend=c("Weather dominant", "Weather foremost", "Fuel dominant", "Fuel foremost"), col=c("#F3CDAA","#FAECDE", "#A4D9B7", "#E5F5EA"))
})

combined_p2p3 <- cowplot::plot_grid(g2, g3, ncol=2, rel_widths=c(1,1))
combined_plot <- cowplot::plot_grid(g1, combined_p2p3, ncol=1, rel_heights = c(1,1))
ggsave(paste0("Variable importance_",ci,".tiff"), combined_plot, width = 8, height = 8)


