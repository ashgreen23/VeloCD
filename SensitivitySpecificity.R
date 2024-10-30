SensitivitySpecificity <- function(ResultFileLocation, GeneSignature, MyMetadata, MyMethod, Output, Class1, Class2, InterName, PredGroup, thresholdProb)  {
  ResFiles <- list.files(path=ResultFileLocation, pattern="SummaryTable_Sum_SelfNotRemoved_")
  ResFiles <- ResFiles[grepl(InterName, ResFiles) ]
  ResFiles <- ResFiles[grepl(PredGroup, ResFiles) ]
  Metadata <- read.csv(MyMetadata)
  setwd(ResultFileLocation)
  Sensitivity_set <- c()
  Specificity_set <- c()
  Severe_Perc <- c()
  Mild_Perc <- c()
  EmbedNNs <- c()
  Methods <- c()
  TPNNs <- c()
  for (file in ResFiles) {
    MyFile <- read.csv(file)
    NN_info <- gsub("SummaryTable_Sum_SelfNotRemoved_", "", file)
    NN_info <- gsub(".csv", "", NN_info)
    NN_info <- gsub(paste0("_", MyMethod), "", NN_info)
    NN_info <- gsub(paste0(MyMethod, "_"), "", NN_info)
    NN_info <- gsub("2D_", "", NN_info)
    NN_info <- gsub("3D_", "", NN_info)
    NN_info <- gsub("Sample", "", NN_info)
    NN_info <- gsub(paste0("_", InterName), "", NN_info)
    NN_info <- gsub("Group1", "", NN_info)
    NN_info <- gsub("Group2", "", NN_info)
    NN_info <- gsub("TPNN", "", NN_info)
    Methods <- append(Methods, MyMethod)
    if (MyMethod == "PCA") {
      EmbedNNs <- append(EmbedNNs, "Not Applicable to PCA")
      NN_info <- gsub("TPNN", "", NN_info)
      TPNNs <- append(TPNNs, as.numeric(NN_info))
    }
    else {
      NN_info <- strsplit(NN_info, "_", fixed=TRUE)
      TPNNs <- append(TPNNs, as.numeric(NN_info[[1]][2]))
      EmbedNNs <- append(EmbedNNs, as.numeric(gsub("NN", "", NN_info[[1]][1])))
    }
    Metadata2 <- Metadata[match(MyFile$Sample, Metadata$Sample),]
    if (identical(Metadata2$Sample, MyFile$Sample) == TRUE) {
      MyFile$Type <- Metadata2$Group
      Group1Severe <- MyFile[MyFile$Type == Class1,]
      Group2Milder <- MyFile[MyFile$Type == Class2,]
      Sensitivity <- nrow(Group1Severe[Group1Severe[Class1] >= thresholdProb,])/nrow(Group1Severe)
      Sensitivity_set <- append(Sensitivity_set, Sensitivity)
      Specificity <- nrow(Group2Milder[Group2Milder[Class1] < thresholdProb,])/nrow(Group2Milder)
      Specificity_set <- append(Specificity_set, Specificity)
      Percent_Severe <- nrow(Group1Severe[Group1Severe[Class1] >= thresholdProb,])/nrow(Group1Severe)*100
      Severe_Perc <- append(Severe_Perc, Percent_Severe)
      Percent_Mild <- nrow(Group2Milder[Group2Milder[Class1] < thresholdProb,])/nrow(Group2Milder)*100
      Mild_Perc <- append(Mild_Perc, Percent_Mild)
    }
  }
  setwd(Output)
  Results <- data.frame(EmbedNNs, TPNNs, Methods, Severe_Perc, Mild_Perc, Sensitivity_set, Specificity_set)
  colnames(Results) <- c("Embedding_NN", "Transition_Probability_NN", "Method", "Group1_Correctly_Predicted", "Group2_Correctly_Predicted", "Sensitivity", "Specificity")
  write.csv(Results, paste0(GeneSignature, "_", InterName, "_", PredGroup, "_SensitivitySpecificities.csv"), row.names=FALSE)
}


