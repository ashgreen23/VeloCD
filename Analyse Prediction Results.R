library("stringr"); library("pROC"); library("ROCit")
#on individual person level
AnalyseResultsIndividualLevel <- function(MyLocation, embedding, dimensions, outdir) {
  MyMetaFiles <- list.files(path = MyLocation, pattern = "Metadata_")
  MyMetaFiles <- gsub("Metadata_", "", MyMetaFiles)
  MyMetaFiles <- gsub(".csv", "", MyMetaFiles)
  for (Run in MyMetaFiles) {
    MySplit <- str_split(Run, "_")
    Person <- MySplit[[1]][1] #get info about the run
    if (embedding == "UMAP") {
      ResultsLocation <- paste0(MyLocation, "Results/Files/", embedding, "_", dimensions, "/Transition_Probability/Group_Level/", Run)
    }
    else {
      ResultsLocation <- paste0(MyLocation, "Results/Files/", embedding, "/Transition_Probability/Group_Level/", Run)
    }
    MyResultFiles <- list.files(path = ResultsLocation, pattern = "SummaryTable_Sum_SelfNotRemoved_") #get the files with the transition probabilities
    MyResultFiles <- MyResultFiles[grepl("Group2", MyResultFiles)]
    if (embedding == "PCA") {
        MyResultFiles <- MyResultFiles[grepl(dimensions, MyResultFiles)]
    }
    Within14Result <- c()
    NeverIRISResult <- c()
    Filename <- c()
    setwd(ResultsLocation)
    for (result in MyResultFiles) { #for each of these results files (different TPNN)
      result2 <- read.csv(result)
      Within14Result <- append(Within14Result, result2[result2$Sample == Person,]$Within14)
      NeverIRISResult <- append(NeverIRISResult, result2[result2$Sample == Person,]$NeverIRIS) #get the probability values
      if (embedding == "PCA") {
        Filename <- append(Filename, gsub(paste0("PCA_", dimensions, "_"), "", gsub("SummaryTable_Sum_SelfNotRemoved_", "", gsub("SampleGroup2.csv", "", result))))
      }
      else if (embedding == "UMAP" & dimensions == "3D") {
        Filename <- append(Filename, gsub("3D",  "", gsub("SummaryTable_Sum_SelfNotRemoved_", "", gsub("SampleGroup2.csv", "", gsub(paste0("_", embedding), "", result)))))
      }
      else {
        Filename <- append(Filename, gsub("SummaryTable_Sum_SelfNotRemoved_", "", gsub("SampleGroup2.csv", "", gsub(paste0("_", embedding), "", result))))
      }
    }
    ResultsSummary <- data.frame(Filename, Within14Result, NeverIRISResult)
    colnames(ResultsSummary) <- c("NN", "Within14_TP", "NeverIRIS_TP") #save the results
    setwd(outdir)
    write.csv(ResultsSummary, paste0(Run, "_", embedding, "_", dimensions, "_Results.csv"), row.names=FALSE)
    ResCol <- ResultsSummary[, colnames(ResultsSummary) == "Within14_TP"]
    Thresholds <- seq(0, 1, 0.05)
    NumGreater <- c()
    for (threshold in Thresholds) {
      NumGreater <- append(NumGreater, length(ResCol[ResCol >= threshold])) #number of runs at each threshold value where the subjects probability of transitioning to its actual outcome group is greater then the threshold
    }
    Summary <- data.frame(Thresholds, NumGreater)
    write.csv(Summary, paste0(Run, "_", embedding, "_", dimensions, "Thresholds.csv"), row.names=FALSE)
  }
}

AnalyseResultsNNLevel <- function(PCAYes, PCAYes3D, tSNEYes, UMAPYes, UMAPYes3D, MyLocation, outdir, ExampleOutputFilePre) { #
  setwd(outdir)
  if (UMAPYes3D == "yes") {
    UMAP_3D <- read.csv(paste0(ExampleOutputFilePre, "UMAP_3D_Results.csv")) #change based on user analysis
    NNList_UMAP <- UMAP_3D$NN   #get list of tpnn values
    MyMetaFiles <- list.files(path = MyLocation, pattern = "Metadata_") #
    MyMetaFiles <- gsub("Metadata_", "", MyMetaFiles)
    MyMetaFiles <- gsub(".csv", "", MyMetaFiles)
    for (TPNNvalue in NNList_UMAP) { #For each diff TPNN res:
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      for (subject in MyMetaFiles) {
        ResFile <- list.files(path = paste0(MyLocation, "Results/Files/UMAP_3D/Transition_Probability/Group_Level/", subject), pattern = paste0(TPNNvalue, "_"))
        ResFile <- ResFile[grepl("SummaryTable", ResFile)]
        ResFile <- ResFile[grepl("Group1", ResFile)]
        MyResFile <- read.csv(paste0(MyLocation, "Results/Files/UMAP_3D/Transition_Probability/Group_Level/", subject, "/", ResFile))
        MySplit <- str_split(subject, "_")
        Person <- MySplit[[1]][1]
        MyRow <- MyResFile[MyResFile$Sample == Person,] #get the row for this subject and save it to a dataframe
        df <- rbind(df, MyRow)
      }
      setwd(outdir)
      write.csv(df, paste0(TPNNvalue, "_Results_UMAP_3D.csv"), row.names=FALSE) #save the results
      Thresholds <- seq(0, 1, 0.05)
      Sensitivity <- c()
      Specificity <- c()
      for (threshold in Thresholds) { #for each probability threshold, get the sensitivity and specificity of the run
        Sensitivity <- append(Sensitivity, (nrow(df[df$Within14 >= threshold,])/nrow(df[df$Group == "Greater14",]))*100)
      }
      SensSpec <- data.frame(Thresholds, Sensitivity)
      write.csv(SensSpec, paste0(TPNNvalue, "_Sensitivity_UMAP_3D.csv"), row.names=FALSE)
    }
  }
  if (UMAPYes == "yes") {
    UMAP_2D <- read.csv(paste0(ExampleOutputFilePre, "UMAP_2D_Results.csv")) #change based on your analysis
    NNList_UMAP <- UMAP_2D$NN   #get list of tpnn values
    MyMetaFiles <- list.files(path = MyLocation, pattern = "Metadata_") #
    MyMetaFiles <- gsub("Metadata_", "", MyMetaFiles)
    MyMetaFiles <- gsub(".csv", "", MyMetaFiles)
    for (TPNNvalue in NNList_UMAP) { #For each diff TPNN res:
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      for (subject in MyMetaFiles) {
        ResFile <- list.files(path = paste0(MyLocation, "Results/Files/UMAP_2D/Transition_Probability/Group_Level/", subject), pattern = paste0(TPNNvalue, "_"))
        ResFile <- ResFile[grepl("SummaryTable", ResFile)]
        ResFile <- ResFile[grepl("Group1", ResFile)]
        MyResFile <- read.csv(paste0(MyLocation, "Results/Files/UMAP_2D/Transition_Probability/Group_Level/", subject, "/", ResFile))
        MySplit <- str_split(subject, "_")
        Person <- MySplit[[1]][1]
        MyRow <- MyResFile[MyResFile$Sample == Person,] #get the row for this subject and save it to a dataframe
        df <- rbind(df, MyRow)
      }
      setwd(outdir)
      write.csv(df, paste0(TPNNvalue, "_Results_UMAP_2D.csv"), row.names=FALSE) #save the results
      Thresholds <- seq(0, 1, 0.05)
      Sensitivity <- c()
      Specificity <- c()
      for (threshold in Thresholds) { #for each probability threshold, get the sensitivity and specificity of the run
        Sensitivity <- append(Sensitivity, (nrow(df[df$Within14 >= threshold,])/nrow(df[df$Group == "Greater14",]))*100)
      }
      SensSpec <- data.frame(Thresholds, Sensitivity)
      write.csv(SensSpec, paste0(TPNNvalue, "_Sensitivity_UMAP_2D.csv"), row.names=FALSE)
    }
  }
  if (tSNEYes == "yes") {
    tSNE_2D <- read.csv(paste0(ExampleOutputFilePre, "tSNE_2D_Results.csv"))
    NNList_tSNE <- tSNE_2D$NN   #get list of tpnn values
    MyMetaFiles <- list.files(path = MyLocation, pattern = "Metadata_") #
    MyMetaFiles <- gsub("Metadata_", "", MyMetaFiles)
    MyMetaFiles <- gsub(".csv", "", MyMetaFiles)
    for (TPNNvalue in NNList_tSNE) { #For each diff TPNN res:
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      for (subject in MyMetaFiles) {
        ResFile <- list.files(path = paste0(MyLocation, "Results/Files/tSNE/Transition_Probability/Group_Level/", subject), pattern = paste0(TPNNvalue, "_"))
        ResFile <- ResFile[grepl("SummaryTable", ResFile)]
        ResFile <- ResFile[grepl("Group1", ResFile)] # or 2?
        MyResFile <- read.csv(paste0(MyLocation, "Results/Files/tSNE/Transition_Probability/Group_Level/", subject, "/", ResFile))
        MySplit <- str_split(subject, "_")
        Person <- MySplit[[1]][1]
        MyRow <- MyResFile[MyResFile$Sample == Person,] #get the row for this subject and save it to a dataframe
        df <- rbind(df, MyRow)
      }
      setwd(outdir)
      write.csv(df, paste0(TPNNvalue, "_Results_tSNE_2D.csv"), row.names=FALSE) #save the results
      Thresholds <- seq(0, 1, 0.05)
      Sensitivity <- c()
      Specificity <- c()
      for (threshold in Thresholds) { #for each probability threshold, get the sensitivity and specificity of the run
        Sensitivity <- append(Sensitivity, (nrow(df[df$Within14 >= threshold,])/nrow(df[df$Group == "Greater14",]))*100)
      }
      SensSpec <- data.frame(Thresholds, Sensitivity)
      write.csv(SensSpec, paste0(TPNNvalue, "_Sensitivity_tSNE_2D.csv"), row.names=FALSE)
    }
  }
  if (PCAYes == "yes") {
    PCA_2D <- read.csv(paste0(ExampleOutputFilePre, "PCA_2D_Results.csv"))
    NNList_PCA <- PCA_2D$NN   #get list of tpnn values
    MyMetaFiles <- list.files(path = MyLocation, pattern = "Metadata_") #
    MyMetaFiles <- gsub("Metadata_", "", MyMetaFiles)
    MyMetaFiles <- gsub(".csv", "", MyMetaFiles)
    for (TPNNvalue in NNList_PCA) { #For each diff TPNN res:
      df <- data.frame(matrix(ncol = 3, nrow = 0))
      for (subject in MyMetaFiles) {
        ResFile <- list.files(path = paste0(MyLocation, "Results/Files/PCA/Transition_Probability/Group_Level/", subject), pattern = paste0(TPNNvalue, "Sample"))
        ResFile <- ResFile[grepl("SummaryTable", ResFile)]
        ResFile <- ResFile[grepl("2D", ResFile)]
        ResFile <- ResFile[grepl("Group1", ResFile)] # or 2?
        MyResFile <- read.csv(paste0(MyLocation, "Results/Files/PCA/Transition_Probability/Group_Level/", subject, "/", ResFile))
        MySplit <- str_split(subject, "_")
        Person <- MySplit[[1]][1]
        MyRow <- MyResFile[MyResFile$Sample == Person,] #get the row for this subject and save it to a dataframe
        df <- rbind(df, MyRow)
      }
      setwd(outdir)
      write.csv(df, paste0(TPNNvalue, "_Results_PCA_2D.csv"), row.names=FALSE) #save the results
      Thresholds <- seq(0, 1, 0.05)
      Sensitivity <- c()
      Specificity <- c()
      for (threshold in Thresholds) { #for each probability threshold, get the sensitivity and specificity of the run
        Sensitivity <- append(Sensitivity, (nrow(df[df$Within14 >= threshold,])/nrow(df[df$Group == "Greater14",]))*100)
      }
      SensSpec <- data.frame(Thresholds, Sensitivity)
      write.csv(SensSpec, paste0(TPNNvalue, "_Sensitivity_PCA_2D.csv"), row.names=FALSE)
    }
    setwd(outdir)
    if (PCAYes3D == "yes") {
      NNList_PCA_3D <- gsub("2D", "3D", NNList_PCA)
      for (TPNNvalue in NNList_PCA_3D) { #equivalent for the 3D PCAs
        df <- data.frame(matrix(ncol = 3, nrow = 0))
        for (subject in MyMetaFiles) {
          ResFile <- list.files(path = paste0(MyLocation, "Results/Files/PCA/Transition_Probability/Group_Level/", subject), pattern = paste0(TPNNvalue, "Sample"))
          ResFile <- ResFile[grepl("SummaryTable", ResFile)]
          ResFile <- ResFile[grepl("3D", ResFile)]
          ResFile <- ResFile[grepl("Group1", ResFile)]
          MyResFile <- read.csv(paste0(MyLocation, "Results/Files/PCA/Transition_Probability/Group_Level/", subject, "/", ResFile))
          MySplit <- str_split(subject, "_")
          Person <- MySplit[[1]][1]
          MyRow <- MyResFile[MyResFile$Sample == Person,]
          df <- rbind(df, MyRow)
        }
        setwd(outdir)
        write.csv(df, paste0(TPNNvalue, "_Results_PCA_3D.csv"), row.names=FALSE)
        Thresholds <- seq(0, 1, 0.05)
        Sensitivity <- c()
        Specificity <- c()
        for (threshold in Thresholds) {
          Sensitivity <- append(Sensitivity, (nrow(df[df$Within14 >= threshold,])/nrow(df[df$Group == "Greater14",]))*100)
        }
        SensSpec <- data.frame(Thresholds, Sensitivity)
        write.csv(SensSpec, paste0(TPNNvalue, "_Sensitivity_PCA_3D.csv"), row.names=FALSE)
      }
    }    
  }
}

AnalyseResultsIndividualLevel(MyLocation="/Velo_DISCO_analysis/", embedding="PCA", dimensions="2D", outdir="/Velo_DISCO_analysis/Summary/")
AnalyseResultsIndividualLevel(MyLocation="/Velo_DISCO_analysis/", embedding="PCA", dimensions="3D", outdir="/Velo_DISCO_analysis/Summary/")
AnalyseResultsIndividualLevel(MyLocation="/Velo_DISCO_analysis/", embedding="tSNE", dimensions="2D", outdir="/Velo_DISCO_analysis/Summary/")
AnalyseResultsIndividualLevel(MyLocation="/Velo_DISCO_analysis/", embedding="UMAP", dimensions="2D", outdir="/Velo_DISCO_analysis/Summary/")
AnalyseResultsIndividualLevel(MyLocation="/Velo_DISCO_analysis/", embedding="UMAP", dimensions="3D", outdir="/Velo_DISCO_analysis/Summary/")
AnalyseResultsNNLevel(PCAYes="yes", PCAYes3D="yes", tSNEYes="yes", UMAPYes="yes", UMAPYes3D="yes", MyLocation="/Velo_DISCO_analysis/", outdir="/Velo_DISCO_analysis/Summary/", ExampleOutputFilePre="LAI550A24_")

TPNNNfiles <- list.files(path="/Summary/", pattern="_Sensitivity_")
TPNNNfilesDay1PCA <- TPNNNfiles[grepl("PCA", TPNNNfiles)]
for (file in TPNNNfilesDay1PCA) {
  MyFile <- read.csv(paste0("", file)) #
  if (nrow(MyFile[MyFile$Sensitivity >= 80 & MyFile$Thresholds > 0.25,]) >= 1) {
    print(file)
    print(MyFile[MyFile$Sensitivity >= 80 & MyFile$Thresholds > 0.25,]) #
  } 
}


library("ggplot2")

TPNNNfilesDay1PCA <- TPNNNfiles[grepl("PCA", TPNNNfiles)]
TPNNNfilesDay1PCA <- TPNNNfilesDay1PCA[grepl("2D", TPNNNfilesDay1PCA)]

Vals <- c()
TPNNVal <- c()
for (file in TPNNNfilesDay1PCA) {
  MyFile <- read.csv(paste0("/Summary/", file))
  Vals <- append(Vals, MyFile[MyFile$Thresholds == 0.3,]$Sensitivity)
  TPNNVal <- append(TPNNVal, as.numeric(as.character(gsub("_Sensitivity_PCA_2D.csv", "", gsub("TPNN", "", file)))))
}

Resys <- data.frame(TPNNVal, Vals)
setwd("/Summary/")
png("TPNN_Range.png", height=500, width=700)
p <- ggplot(Resys, aes(x=TPNNVal, y=Vals)) + geom_point(size=3) + theme_bw() + theme(text=element_text(size=20)) + xlab("\nTPNN") + ylab("Correctly Preicted (%)\n")
print(p)
dev.off()