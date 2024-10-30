TransitionProbability <- function(filename, selfremoved, location, metaGroup, outputdir, grouptype, version, Randomsampling, UserThreshold) {
  `%ni%` = Negate(`%in%`)
  library("gtools"); library("reshape2"); library("ggplot2"); library("grid"); library("tidyverse"); library("scatterplot3d")
  metadata <- read.csv(metaGroup)
  locfilename <- paste0(location, filename, sep="")
  TPNN <- read.csv(locfilename)
  setwd(outputdir)
  rownames(TPNN) <- TPNN$X 
  TPNN <- TPNN[,-1]
  #rownames(TPNN) <- gsub("\\.", "", rownames(TPNN))
  rownames(TPNN) <- gsub('"', '', rownames(TPNN))
  #colnames(TPNN) <- gsub("\\.", "", colnames(TPNN))
  colnames(TPNN) <- gsub('"', '', colnames(TPNN))
  metadata <- metadata[metadata$Sample %in% rownames(TPNN),]
  metadata <- metadata[match(rownames(TPNN), metadata$Sample),]
  if (grouptype == "Strings") {
    metadata$Group <- unlist(lapply(metadata$Group, toString))
  }
  samples <- c()
  groups2 <- c()
  Meanvalues <- data.frame()
  Medianvalues <- data.frame()
  SumVal <- data.frame()
  samplesRS <- c()
  groups2RS <- c()
  MeanvaluesRS <- data.frame()
  MedianvaluesRS <- data.frame()
  SumValRS <- data.frame()
  randomsamplingnum <- c()
  minGroupSize <- min(as.data.frame(table(metadata$Group))$Freq) #will be 2 when self removed and your at the timepoint with both samples, otherwise 1 
  if (selfremoved == "TRUE" & minGroupSize > 1) {
    minGroupSize <- minGroupSize - 1
  }
  for (value in (1:nrow(TPNN))) { #for each sample:
    sampleName <- rownames(TPNN)[value]
    TPNNnome <- TPNN[value, ] #get self row
    if (selfremoved == "TRUE") {
      TPNNnome <- TPNNnome[, -value] #remove self column
    }
    groups <- metadata[metadata$Sample %in% names(TPNNnome),] #
    groups <- metadata[match(names(TPNNnome), metadata$Sample),]$Group #get groups
    TPNNnew <- data.frame(t(TPNNnome), names(TPNNnome), groups)
    colnames(TPNNnew) <- c("TransitionProbability", "Sample", "Group")
    ### Calculate using X number of randomly chosen samples
    randomsamplingnum <- append(minGroupSize, randomsamplingnum)
    randomDF <- data.frame()
    for (groupy in unique(TPNNnew$Group)) {
      SamplesInGroupy <- TPNNnew[TPNNnew$Group == groupy, ]
      chosensample <- sample_n(SamplesInGroupy, minGroupSize, replace=FALSE)
      randomDF <- rbind(randomDF, chosensample)
    }
    randomDF_means <- aggregate( TransitionProbability ~ Group, randomDF, mean )
    randomDF_medians <- aggregate( TransitionProbability ~ Group, randomDF, median )
    randomDF_sums <- aggregate( TransitionProbability ~ Group, randomDF, sum )
    randomDF_sumsPercent <- randomDF_sums[,2]*100
    RandomSamplingDF <- data.frame(randomDF_means$Group, randomDF_means$TransitionProbability, randomDF_medians$TransitionProbability, randomDF_sums$TransitionProbability, randomDF_sumsPercent)
    colnames(RandomSamplingDF) <- c("Group",  "MeanRS", "MedianRS",  "SumRS", "SumPercentageRS")
    sortedgroups <- sort(RandomSamplingDF$Group)
    RandomSamplingDF <- RandomSamplingDF[match(sortedgroups, RandomSamplingDF$Group),] #order alphanumerically so later on it matches the sort unique of the matadata
    #RandomSamplingDF$"MeanRank" <- rank(RandomSamplingDF[,2])
    #RandomSamplingDF$"MedianRank" <- rank(RandomSamplingDF[,3])
    #RandomSamplingDF$"SumRank" <- rank(RandomSamplingDF[,4])
    means <- aggregate( TransitionProbability ~ Group, TPNNnew, mean )     #get group mean and medians and sum for each grpups TP to the sample
    medians <- aggregate( TransitionProbability ~ Group, TPNNnew, median )
    sums <- aggregate( TransitionProbability ~ Group, TPNNnew, sum )
    sumsAsPercent <- sums[,2]*100
    if (identical(means$Group, medians$Group) == TRUE) {
      values <- data.frame(means$Group, means$TransitionProbability, medians$TransitionProbability, sums$TransitionProbability, sumsAsPercent)
      colnames(values) <- c("Group", "Mean", "Median", "Sum", "Sum (%)") #put together
      valsvalues <- sort(values$Group) #I believe this is what allows the colnames addiiton downstream
      values <- values[match(valsvalues, values$Group),] #order alphanumerically by Group
      #values$MeanRank <- rank(values$Mean)       #rank groups based on mean and median and add these to the df, if two values are the same they get the same rank
      #values$MedianRank <- rank(values$Median)
      #values$SumRank <- rank(values$Sum)
      #if (selfremoved == "TRUE") {
      #  filename2 <- paste(sampleName, "_SelfRemoved", filename, sep="")
      #  filename2RS <- paste(sampleName, "_SelfRemoved_Random_Sampling", filename, sep="")
      #}
      #else {
      #  filename2 <- paste(sampleName, "_SelfNotRemoved", filename, sep="")
      #  filename2RS <- paste(sampleName, "_SelfNotRemoved_Random_Sampling", filename, sep="")
      #}
      #write.csv(values, filename2, row.names=FALSE) #return as file
      #write.csv(RandomSamplingDF, filename2RS, row.names=FALSE) #return as file
      if (any(sort(unique(metadata$Group)) %ni% sort(values$Group))) {
        absent <- setdiff(metadata$Group, values$Group)
        positionbefore <- match(absent, sort(unique(metadata$Group)))-1
        positionafter <- match(absent, sort(unique(metadata$Group)))
        if (positionbefore == 0 & positionafter != nrow(metadata)) {
          myMeans <- c(NA, values$Mean[positionafter:length(values$Mean)])
          myMedians <- c(NA, values$Median[positionafter:length(values$Median)])
          mySum <- c(NA, values$Sum[positionafter:length(values$Sum)])
        }
        else if (positionbefore != 0  & positionafter == nrow(metadata)) {
          #all samples:
          myMeans <- c(values$Mean[1:positionbefore], NA)
          myMedians <- c(values$Median[1:positionbefore], NA)
          mySum <- c(values$Sum[1:positionbefore], NA)
        }
        else if (positionbefore != 0 & positionafter != nrow(metadata)) {
          #all samples:
          myMeans <- c(values$Mean[1:positionbefore], NA, values$Mean[positionafter:length(values$Mean)])
          myMedians <- c(values$Median[1:positionbefore], NA, values$Median[positionafter:length(values$Median)])
          mySum <- c(values$Sum[1:positionbefore], NA, values$Sum[positionafter:length(values$Sum)])
        }
        samples <- append(samples, sampleName)
        groups2 <- append(groups2, unique(metadata[metadata$Sample == sampleName,]$Group))
        Meanvalues <- rbind(Meanvalues, myMeans)
        Medianvalues <- rbind(Medianvalues, myMedians)
        SumVal <- rbind(SumVal, mySum)
      }
      else {
        samples <- append(samples, sampleName)
        groups2 <- append(groups2, unique(metadata[metadata$Sample == sampleName,]$Group))
        Meanvalues <- rbind(Meanvalues, values$Mean)
        Medianvalues <- rbind(Medianvalues, values$Median)
        SumVal <- rbind(SumVal, values$Sum)
      }
      if (any(sort(unique(metadata$Group)) %ni% sort(RandomSamplingDF$Group))) {
        absent <- setdiff(metadata$Group, RandomSamplingDF$Group)
        positionbefore <- match(absent, sort(unique(metadata$Group)))-1
        positionafter <- match(absent, sort(unique(metadata$Group)))
        if (positionbefore == 0 & positionafter != nrow(metadata)) {
          myMeansRS <- c(NA, RandomSamplingDF$MeanRS[positionafter:length(RandomSamplingDF$MeanRS)])
          myMediansRS <- c(NA, RandomSamplingDF$MedianRS[positionafter:length(RandomSamplingDF$MedianRS)])
          mySumRS <- c(NA, RandomSamplingDF$SumRS[positionafter:length(RandomSamplingDF$SumRS)])
        }
        else if (positionbefore != 0  & positionafter == nrow(metadata)) {
          myMeansRS <- c(RandomSamplingDF$MeanRS[1:positionbefore], NA)
          myMediansRS <- c(RandomSamplingDF$MedianRS[1:positionbefore], NA)
          mySumRS <- c(RandomSamplingDF$SumRS[1:positionbefore], NA)
        }
        else if (positionbefore != 0 & positionafter != nrow(metadata)) {
          myMeansRS <- c(RandomSamplingDF$MeanRS[1:positionbefore], NA, RandomSamplingDF$MeanRS[positionafter:length(RandomSamplingDF$MeanRS)])
          myMediansRS <- c(RandomSamplingDF$MedianRS[1:positionbefore], NA, RandomSamplingDF$MedianRS[positionafter:length(RandomSamplingDF$MedianRS)])
          mySumRS <- c(RandomSamplingDF$SumRS[1:positionbefore], NA, RandomSamplingDF$SumRS[positionafter:length(RandomSamplingDF$SumRS)])
        }
        samplesRS <- append(samplesRS, sampleName)
        groups2RS <- append(groups2RS, unique(metadata[metadata$Sample == sampleName,]$Group))
        MeanvaluesRS <- rbind(MeanvaluesRS, myMeansRS)
        MedianvaluesRS <- rbind(MedianvaluesRS, myMediansRS)
        SumValRS <- rbind(SumValRS, mySumRS)
      }
      else {
        samplesRS <- append(samplesRS, sampleName)
        groups2RS <- append(groups2RS, unique(metadata[metadata$Sample == sampleName,]$Group))
        MeanvaluesRS <- rbind(MeanvaluesRS, RandomSamplingDF$MeanRS)
        MedianvaluesRS <- rbind(MedianvaluesRS, RandomSamplingDF$MedianRS)
        SumValRS <- rbind(SumValRS, RandomSamplingDF$SumRS)
      }
    }
  }
  MeansTable <- data.frame(samples, groups2, Meanvalues)
  SumsTable <- data.frame(samples, groups2, SumVal)
  MediansTable <- data.frame(samples, groups2, Medianvalues)
  colnames(MeansTable) <- c("Sample", "Group", sort(unique(metadata$Group)))
  colnames(MediansTable) <- c("Sample", "Group", sort(unique(metadata$Group)))
  colnames(SumsTable) <- c("Sample", "Group", sort(unique(metadata$Group)))
  #Random Sampling:
  MeansTableRS <- data.frame(samplesRS, groups2RS, MeanvaluesRS)
  SumsTableRS <- data.frame(samplesRS, groups2RS, SumValRS)
  MediansTableRS <- data.frame(samplesRS, groups2RS, MedianvaluesRS)
  colnames(MeansTableRS) <- c("Sample", "Group", sort(unique(metadata$Group)))
  colnames(MediansTableRS) <- c("Sample", "Group", sort(unique(metadata$Group)))
  colnames(SumsTableRS) <- c("Sample", "Group", sort(unique(metadata$Group)))
  NumGenesRS <- data.frame(samplesRS, randomsamplingnum)
  colnames(NumGenesRS) <- c("Sample", "Number of genes used for random sampling")
  if (selfremoved == "TRUE") {
    RSGeneNumName <- paste0("Random_Sampling_Gene_Numbers_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MeansTabRS <- paste0("SummaryTable_Mean_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MedTabRS <- paste0("SummaryTable_Median_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    SumTabRS <- paste0("SummaryTable_Sum_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MeansTab <- paste0("SummaryTable_Mean_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MedTab <- paste0("SummaryTable_Median_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    SumTab <- paste0("SummaryTable_Sum_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
  }
  else {
    RSGeneNumName <- paste0("Random_Sampling_Gene_Numbers_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MeansTabRS <- paste0("SummaryTable_Mean_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MedTabRS <- paste0("SummaryTable_Median_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    SumTabRS <- paste0("SummaryTable_Sum_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MeansTab <- paste0("SummaryTable_Mean_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    MedTab <- paste0("SummaryTable_Median_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".csv", sep="")
    SumTab <- paste0("SummaryTable_Sum_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)),  version, ".csv", sep="")
  }
  #write.csv(MeansTable, MeansTab, row.names=FALSE)
  #write.csv(MediansTable, MedTab, row.names=FALSE)
  write.csv(SumsTable, SumTab, row.names=FALSE)
  if (Randomsampling == TRUE) {
    write.csv(NumGenesRS, RSGeneNumName, row.names=FALSE)
    #write.csv(MeansTableRS, MeansTabRS, row.names=FALSE)
    #write.csv(MediansTableRS, MedTabRS, row.names=FALSE)
    write.csv(SumsTableRS, SumTabRS, row.names=FALSE)
  }
  #heatmaps for mean, sums and medians values:
  MeansTableMelt <- melt(MeansTable[,c(1,3:ncol(MeansTable))], id="Sample")
  colnames(MeansTableMelt) <- c("Sample", "variable", "MeanTP")
  if (selfremoved == "TRUE") {
    fignameTimeRS <- paste0("TransitionProbabilityHeatMapMeans_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime3RS <- paste0("TransitionProbabilityHeatMapSums_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTimeMedRS <- paste0("TransitionProbabilityHeatMapMedians_Random_Sampling_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime <- paste0("TransitionProbabilityHeatMapMeans_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime3 <- paste0("TransitionProbabilityHeatMapSums_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTimeMed <- paste0("TransitionProbabilityHeatMapMedians_SelfRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
  }
  else {
    fignameTimeRS <- paste0("TransitionProbabilityHeatMapMeans_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime3RS <- paste0("TransitionProbabilityHeatMapSums_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTimeMedRS <- paste0("TransitionProbabilityHeatMapMedians_Random_Sampling_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime <- paste0("TransitionProbabilityHeatMapMeans_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTime3 <- paste0("TransitionProbabilityHeatMapSums_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
    fignameTimeMed <- paste0("TransitionProbabilityHeatMapMedians_SelfNotRemoved", gsub("Transition_Probability", "", gsub(".csv", "", filename)), version, ".png", sep="")
  }
  #png(fignameTime, width = 880, height = 780)
  #corrhm <- ggplot(data = MeansTableMelt, aes(x=variable, y=Sample, fill=MeanTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
  #print(corrhm)
  #dev.off()
  SumsTableMelt <- melt(SumsTable[,c(1,3:ncol(SumsTable))], id="Sample")
  colnames(SumsTableMelt) <- c("Sample", "variable", "SummedTP")
  #png(fignameTime3, width = 880, height = 780)
  #corrhm3 <- ggplot(data = SumsTableMelt, aes(x=variable, y=Sample, fill=SummedTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
  #print(corrhm3)
  #dev.off()
  MediansTableMelt <- melt(MediansTable[,c(1,3:ncol(MediansTable))], id="Sample")
  colnames(MediansTableMelt) <- c("Sample", "variable", "MedianTP")
  #png(fignameTimeMed, width = 880, height = 780)
  #corrhmMed <- ggplot(data = MediansTableMelt, aes(x=variable, y=Sample, fill=MedianTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
  #print(corrhmMed)
  #dev.off()
  #Random Sampling:
  MeansTableMeltRS <- melt(MeansTableRS[,c(1,3:ncol(MeansTableRS))], id="Sample")
  colnames(MeansTableMeltRS) <- c("Sample", "variable", "MeanTP")
  if (Randomsampling == TRUE) {
    #png(fignameTimeRS, width = 880, height = 780)
    #corrhmRS <- ggplot(data = MeansTableMeltRS, aes(x=variable, y=Sample, fill=MeanTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
    #print(corrhmRS)
    #dev.off()
    SumsTableMeltRS <- melt(SumsTableRS[,c(1,3:ncol(SumsTableRS))], id="Sample")
    colnames(SumsTableMeltRS) <- c("Sample", "variable", "SummedTP")
    #png(fignameTime3RS, width = 880, height = 780)
    #corrhm3RS <- ggplot(data = SumsTableMeltRS, aes(x=variable, y=Sample, fill=SummedTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
    #print(corrhm3RS)
    #dev.off()
    MediansTableMeltRS <- melt(MediansTableRS[,c(1,3:ncol(MediansTableRS))], id="Sample")
    colnames(MediansTableMeltRS) <- c("Sample", "variable", "MedianTP")
    #png(fignameTimeMedRS, width = 880, height = 780)
    #corrhmMedRS <- ggplot(data = MediansTableMeltRS, aes(x=variable, y=Sample, fill=MedianTP)) + geom_tile() + theme_minimal() + theme_bw() + xlab("\n\n...to these groups") + ylab("Transition probability from these samples...\n\n") + theme_grey(base_size =12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1), axis.text.y = element_text(size = 12), axis.text=element_text(size=12)) + scale_fill_distiller(palette = "Reds", direction = 1) #+ scale_fill_gradient(low = "red", high = "green") 
    #print(corrhmMedRS)
    #dev.off()
  }
  MeansTable2 <- MeansTable[, c(3:ncol(MeansTable))]
  MeansTable2[is.na(MeansTable2)] <- 0
  res <- apply(MeansTable2,1,function(x) which(x==max(x)))
  if (typeof(res) == "list") {
    my_groups <- c()
    valuelist <- c()
    meetThreshold <- c()
    for (table in 1:length(res)) {
      mytab <- res[table]
      groups <- rownames(data.frame(mytab))
      my_groups <- append(my_groups, toString(groups))
      valuelist <- append(valuelist, toString(MeansTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholds <- c()
      for (tpval in MeansTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholds <- append(setThresholds, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholds <- append(setThresholds, "No")
        } 
      }
      meetThreshold <- append(meetThreshold, toString(setThresholds))
    }
    Predictions <- data.frame(MeansTable[,1], valuelist, my_groups, meetThreshold)
  }
  else if (typeof(res) == "integer") {
    meetThreshold <- c()
    valuelist <- c()
    for (row in 1:nrow(MeansTable2)) {
      value <- MeansTable2[row, res[row]]
      valuelist <- append(valuelist, toString(value))
      if (value >= UserThreshold) {
        meetThreshold<- append(meetThreshold, "Yes")
      }
      else if (value < UserThreshold) {
        meetThreshold<- append(meetThreshold, "No")
      }
    }
    my_groups <- colnames(MeansTable2)[res]
    Predictions <- data.frame(MeansTable[,1], valuelist, my_groups, meetThreshold)
  }
  thresholdcolname <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilename <- paste0("PredictedGroupsMeanSelfRemoved_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilename <- paste0("PredictedGroupsMeanSelfNotRemoved_", version, "_",  gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(Predictions) <- c("Sample", "Highest Mean Transition Probability", "Highest Mean Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolname)
  #write.csv(Predictions, predsnfilename, row.names=FALSE)
  MediansTable2 <- MediansTable[, c(3:ncol(MediansTable))]
  MediansTable2[is.na(MediansTable2)] <- 0
  res <- apply(MediansTable2,1,function(x) which(x==max(x)))
  #MediansTable2 <- data.frame(c(2,8,1), c(7,3,5),c(7,6,4))   #fake example used to test it
  #colnames(MediansTable2) <- c("6", "12", "18")
  if (typeof(res) == "list") {
    my_groups <- c()
    valuelist <- c()
    meetThreshold <- c()
    for (table in 1:length(res)) {
      mytab <- res[table]
      groups <- rownames(data.frame(mytab))
      my_groups <- append(my_groups, toString(groups))
      valuelist <- append(valuelist, toString(MediansTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholds <- c()
      for (tpval in MediansTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholds <- append(setThresholds, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholds <- append(setThresholds, "No")
        }
      }
      meetThreshold <- append(meetThreshold, toString(setThresholds))
    }
    Predictions <- data.frame(MediansTable[,1], valuelist, my_groups, meetThreshold)
  }
  else if (typeof(res) == "integer") {
    meetThreshold <- c()
    valuelist <- c()
    for (row in 1:nrow(MediansTable2)) {
      value <- MediansTable2[row, res[row]]
      valuelist <- append(valuelist, toString(value))
      if (value >= UserThreshold) {
        meetThreshold <- append(meetThreshold, "Yes")
      }
      else if (value < UserThreshold) {
        meetThreshold<- append(meetThreshold, "No")
      }
    }
    my_groups <- colnames(MediansTable2)[res]
    Predictions <- data.frame(MediansTable[,1], valuelist, my_groups, meetThreshold)
  }
  thresholdcolname <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilename <- paste0("PredictedGroupsMedianSelfRemoved_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilename <- paste0("PredictedGroupsMedianSelfNotRemoved_", version, "_",  gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(Predictions) <- c("Sample", "Highest Median Transition Probability", "Highest Median Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolname)
  #write.csv(Predictions, predsnfilename, row.names=FALSE)
  SumsTable2 <- SumsTable[, c(3:ncol(SumsTable))]
  SumsTable2[is.na(SumsTable2)] <- 0
  res <- apply(SumsTable2,1,function(x) which(x==max(x)))
  if (typeof(res) == "list") {
    my_groups <- c()
    valuelist <- c()
    meetThreshold <- c()
    for (table in 1:length(res)) {
      mytab <- res[table]
      groups <- rownames(data.frame(mytab))
      my_groups <- append(my_groups, toString(groups))
      valuelist <- append(valuelist, toString(SumsTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholds <- c()
      for (tpval in SumsTable2[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholds <- append(setThresholds, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholds <- append(setThresholds, "No")
        }
      }
      meetThreshold <- append(meetThreshold, toString(setThresholds))
    }
    Predictions <- data.frame(SumsTable[,1], valuelist, my_groups, meetThreshold)
  }
  else if (typeof(res) == "integer") {
    meetThreshold <- c()
    valuelist <- c()
    for (row in 1:nrow(SumsTable2)) {
      value <- SumsTable2[row, res[row]]
      valuelist <- append(valuelist, toString(value))
      if (value >= UserThreshold) {
        meetThreshold<- append(meetThreshold, "Yes")
      }
      else if (value < UserThreshold) {
        meetThreshold<- append(meetThreshold, "No")
      }
    }
    my_groups <- colnames(SumsTable2)[res]
    Predictions <- data.frame(SumsTable[,1], valuelist, my_groups, meetThreshold)
  }
  thresholdcolname <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilename <- paste0("PredictedGroupsSumSelfRemoved_", version, "_",  gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilename <- paste0("PredictedGroupsSumSelfNotRemoved_", version, "_",  gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(Predictions) <- c("Sample", "Highest Sum Transition Probability", "Highest Sum Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolname)
  write.csv(Predictions, predsnfilename, row.names=FALSE)
  #Random Sampling:
  MeansTable2RS <- MeansTableRS[, c(3:ncol(MeansTableRS))]
  MeansTable2RS[is.na(MeansTable2RS)] <- 0
  resRS <- apply(MeansTable2RS,1,function(x) which(x==max(x)))
  if (typeof(resRS) == "list") {
    my_groupsRS <- c()
    valuelistRS <- c()
    meetThresholdRS <- c()
    for (table in 1:length(resRS)) {
      mytab <- resRS[table]
      groups <- rownames(data.frame(mytab))
      my_groupsRS <- append(my_groupsRS, toString(groups))
      valuelistRS <- append(valuelistRS, toString(MeansTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholdsRS <- c()
      for (tpval in MeansTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "No")
        }
      }
      meetThresholdRS <- append(meetThresholdRS, toString(setThresholdsRS))
    }
    PredictionsRS <- data.frame(MeansTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  else if (typeof(resRS) == "integer") {
    meetThresholdRS <- c()
    valuelistRS <- c()
    for (row in 1:nrow(MeansTable2RS)) {
      value <- MeansTable2RS[row, resRS[row]]
      valuelistRS <- append(valuelistRS, toString(value))
      if (value >= UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "Yes")
      }
      else if (value < UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "No")
      }
    }
    my_groupsRS <- colnames(MeansTable2RS)[resRS]
    PredictionsRS <- data.frame(MeansTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  thresholdcolnameRS <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilenameRS <- paste0("PredictedGroupsMeanSelfRemoved_RandomSampling_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilenameRS <- paste0("PredictedGroupsMeanSelfNotRemoved_RandomSampling_", version, "_",  gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(PredictionsRS) <- c("Sample", "Highest Mean Transition Probability", "Highest Mean Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolnameRS)
  if (Randomsampling == TRUE) {
    #write.csv(PredictionsRS, predsnfilenameRS, row.names=FALSE)
  }
  MediansTable2RS <- MediansTableRS[, c(3:ncol(MediansTableRS))]
  MediansTable2RS[is.na(MediansTable2RS)] <- 0
  resRS <- apply(MediansTable2RS,1,function(x) which(x==max(x)))
  #MediansTable2 <- data.frame(c(2,8,1), c(7,3,5),c(7,6,4))   #fake example used to test it
  #colnames(MediansTable2) <- c("6", "12", "18")
  if (typeof(resRS) == "list") {
    my_groupsRS <- c()
    valuelistRS <- c()
    meetThresholdRS <- c()
    for (table in 1:length(resRS)) {
      mytab <- resRS[table]
      groups <- rownames(data.frame(mytab))
      my_groupsRS <- append(my_groupsRS, toString(groups))
      valuelistRS <- append(valuelistRS, toString(MediansTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholdsRS <- c()
      for (tpval in MediansTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "No")
        }
      }
      meetThresholdRS <- append(meetThresholdRS, toString(setThresholdsRS))
    }
    PredictionsRS <- data.frame(MediansTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  else if (typeof(resRS) == "integer") {
    meetThresholdRS <- c()
    valuelistRS <- c()
    for (row in 1:nrow(MediansTable2RS)) {
      value <- MediansTable2RS[row, resRS[row]]
      valuelistRS <- append(valuelistRS, toString(value))
      if (value >= UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "Yes")
      }
      else if (value < UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "No")
      }
    }
    my_groupsRS <- colnames(MediansTable2RS)[resRS]
    PredictionsRS <- data.frame(MediansTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  thresholdcolnameRS <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilenameRS <- paste0("PredictedGroupsMedianSelfRemoved_Random_Sampling_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilenameRS <- paste0("PredictedGroupsMedianSelfNotRemoved_Random_Sampling_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(PredictionsRS) <- c("Sample", "Highest Median Transition Probability", "Highest Median Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolnameRS)
  if (Randomsampling == TRUE) {
    #write.csv(PredictionsRS, predsnfilenameRS, row.names=FALSE)
  }
  SumsTable2RS <- SumsTableRS[, c(3:ncol(SumsTableRS))]
  SumsTable2RS[is.na(SumsTable2RS)] <- 0
  resRS <- apply(SumsTable2RS,1,function(x) which(x==max(x)))
  if (typeof(resRS) == "list") {
    my_groupsRS <- c()
    valuelistRS <- c()
    meetThresholdRS <- c()
    for (table in 1:length(resRS)) {
      mytab <- resRS[table]
      groups <- rownames(data.frame(mytab))
      my_groupsRS <- append(my_groupsRS, toString(groups))
      valuelistRS <- append(valuelistRS, toString(SumsTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])]))
      setThresholdsRS <- c()
      for (tpval in SumsTable2RS[table, c(data.frame(mytab)[1:nrow(data.frame(mytab)),1:ncol(data.frame(mytab))])] ) {
        if (tpval >= UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "Yes")
        }
        else if (tpval < UserThreshold) {
          setThresholdsRS <- append(setThresholdsRS, "No")
        }
      }
      meetThresholdRS <- append(meetThresholdRS, toString(setThresholdsRS))
    }
    PredictionsRS <- data.frame(SumsTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  else if (typeof(resRS) == "integer") {
    meetThresholdRS <- c()
    valuelistRS <- c()
    for (row in 1:nrow(SumsTable2RS)) {
      value <- SumsTable2RS[row, resRS[row]]
      valuelistRS <- append(valuelistRS, toString(value))
      if (value >= UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "Yes")
      }
      else if (value < UserThreshold) {
        meetThresholdRS <- append(meetThresholdRS, "No")
      }
    }
    my_groupsRS <- colnames(SumsTable2RS)[resRS]
    PredictionsRS <- data.frame(SumsTableRS[,1], valuelistRS, my_groupsRS, meetThresholdRS)
  }
  thresholdcolnameRS <- paste0("Greater than ", toString(UserThreshold), "?", sep="")
  if (selfremoved == "TRUE") {
    predsnfilenameRS <- paste0("PredictedGroupsSumSelfRemoved_Random_Sampling_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  else {
    predsnfilenameRS <- paste0("PredictedGroupsSumSelfNotRemoved_Random_Sampling_", version, "_", gsub("Transition_Probability_", "", filename), sep="")
  }
  colnames(PredictionsRS) <- c("Sample", "Highest Sum Transition Probability", "Highest Sum Transition Probability (if two groups have the same value, both groups are returned as equally likely)", thresholdcolnameRS)
  if (Randomsampling == TRUE) {
    write.csv(PredictionsRS, predsnfilenameRS, row.names=FALSE)
  }
}
