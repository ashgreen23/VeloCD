InternalBootStrappingTransitonProbabilities <- function(filename, selfremoved, location, metaGroup, outputdir, grouptype, numberofIterations) {
  `%ni%` = Negate(`%in%`)
  library(boot); library(tidyverse); library(Rmisc); library("reshape2"); library(gridExtra)
  metadata <- read.csv(metaGroup)
  locfilename <- paste0(location, filename, sep="")
  TPNN <- read.csv(locfilename)
  setwd(outputdir)
  rownames(TPNN) <- TPNN$X
  TPNN <- TPNN[,-1]
  rownames(TPNN) <- gsub("X.", "", rownames(TPNN))
  rownames(TPNN) <- gsub("\\.", "", rownames(TPNN))
  colnames(TPNN) <- gsub("X.", "", colnames(TPNN))
  colnames(TPNN) <- gsub("\\.", "", colnames(TPNN))
  #rownames(TPNN) <- gsub("[^A-Za-z0-9]", "", rownames(TPNN))
  #colnames(TPNN) <- gsub("[^A-Za-z0-9]", "", colnames(TPNN))
  #metadata$Sample <- gsub("[^A-Za-z0-9,-]", "", metadata$Sample)
  metadata <- metadata[metadata$Sample %in% rownames(TPNN),]
  metadata <- metadata[match(rownames(TPNN), metadata$Sample),]
  if (grouptype == "Strings") {
    metadata$Group <- unlist(lapply(metadata$Group, toString))
  }
  #print(identical(rownames(TPNN), unlist(lapply(metadata$Sample, toString))))
  #print(identical(colnames(TPNN), unlist(lapply(metadata$Sample, toString))))
  sampleNameSet <- c()
  CIDF <- data.frame()
  samplestograb <- min(as.data.frame(table(metadata$Group))$Freq) #will be 2 when self removed and your at the timepoint with both samples, otherwise 1 
  if (selfremoved == "TRUE" & samplestograb > 1) {
    samplestograb <- samplestograb - 1
  }
  for (sample in 1:nrow(TPNN)) {
    sampleRow <- TPNN[sample,]
    sampleName <- rownames(TPNN)[sample]
    sampleNameSet <- append(sampleNameSet, sampleName)
    groups <- metadata[metadata$Sample %in% names(sampleRow),]
    groups <- metadata[match(names(sampleRow), metadata$Sample),]$Group #get groups
    TPNNRow <- data.frame(t(sampleRow), groups)
    rownames(TPNNRow) <- names(sampleRow)
    if (selfremoved == "TRUE") {
      TPNNRow <- TPNNRow[rownames(TPNNRow) != sampleName,] #remove self column
    }
    statsDF <- c()
    for (MyGroup in sort(unique(groups))) {
      GroupSub <- TPNNRow[TPNNRow$groups == MyGroup,]
      means <- c()
      sums <- c()
      for (i in 1:numberofIterations) {
        mybit <- GroupSub[,1]
        if (length(GroupSub[,1]) == 0 ) {
          mybit <- 0
        }
        RS <- sample(mybit, size=samplestograb, replace=FALSE)
        if (length(RS) == 0) {
          meanRS <- 0
          sumRS <- 0
          means <- append(means, meanRS)
          sums <- append(sums, sumRS)
        }
        else if (length(RS) == 1) {
          meanRS <- RS
          sumRS <- RS
          means <- append(means, meanRS)
          sums <- append(sums, sumRS)
        }
        else {
          meanRS <- mean(RS)
          sumRS <- sum(RS)
          means <- append(means, meanRS)
          sums <- append(sums, sumRS)
        }
      }
      MeanCI <- CI(means, ci=0.95)
      SumCI <- CI(sums, ci=0.95)
      LowerMeanCI <- MeanCI[3]
      upperMeanCI <- MeanCI[1]
      MeanMeanCI <- MeanCI[2]
      LowerSumCI <- SumCI[3]
      upperSumCI <- SumCI[1]
      MeanSumCI <- SumCI[2]
      statsDF <- append(statsDF, c(LowerMeanCI, upperMeanCI, MeanMeanCI, LowerSumCI, upperSumCI, MeanSumCI))
    }
    CIDF <- rbind(CIDF, statsDF)
  }
  CIDF <- data.frame(sampleNameSet, CIDF)
  names <- c("Sample")
  for (groupy in sort(unique(groups))) {
    names <- append(names, c(paste0(groupy, "_Mean_CI_Lower", sep=""), paste0(groupy, "_Mean_CI_Higher", sep=""), paste0(groupy, "_Mean_CI_Mean", sep=""), paste0(groupy, "_Sum_CI_Lower", sep=""), paste0(groupy, "_Sum_CI_Higher", sep=""), paste0(groupy, "_Sum_CI_Mean", sep="")))
  }
  colnames(CIDF) <- names
  if (selfremoved == "TRUE") {
    Newfilename <- paste0("SelfRemoved_Transition_Probabilities_Confidence_Intervals", gsub("Transition_Probability_", "", filename), sep="") #remove self column
    meansfigname <- paste0("SelfRemoved_Transition_Probabilities_Confidence_Intervals_Means", gsub("Transition_Probability", "", gsub(".csv", "", filename)), ".png", sep="")
    sumsfigname <- paste0("SelfRemoved_Transition_Probabilities_Confidence_Intervals_Sums", gsub("Transition_Probability", "", gsub(".csv", "", filename)), ".png", sep="")
  }
  else {
    Newfilename <- paste0("SelfNotRemoved_Transition_Probabilities_Confidence_Intervals", gsub("Transition_Probability_", "", filename), sep="") #remove self column
    meansfigname <- paste0("SelfNotRemoved_Transition_Probabilities_Confidence_Intervals_Means", gsub("Transition_Probability", "", gsub(".csv", "", filename)), ".png", sep="")
    sumsfigname <- paste0("SelfNotRemoved_Transition_Probabilities_Confidence_Intervals_Sums", gsub("Transition_Probability", "", gsub(".csv", "", filename)), ".png", sep="")
  }
  write.csv(CIDF, Newfilename, row.names=FALSE) #return as dataframe of sample rows, group cols with CIs in cells, save as csv
  #make a nice summary plot:
  metadata <- metadata[match(CIDF$Sample, metadata$Sample),]
  CIDF$CurrentGroup <- metadata$Group
  meanplots <- list()
  sumplots <- list()
  mylabs <- c()
  count = 1
  cols <- c( "#0072B2", "#E69F00", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#56B4E9")
  for (MYgroup in sort(unique(groups))) {
    mylabs <- append(mylabs, MYgroup)
    groupnamecols <- CIDF[,grepl(paste0(MYgroup, "_", sep=""), colnames(CIDF))]
    means <- groupnamecols[,1:3]
    means$CurrentGroup <- as.factor(metadata$Group)
    means$Sample <- CIDF$Sample
    colnames(means) <- c("LowerCI", "UpperCI", "meanCI", "CurrentGroup", "Sample")
    if (length(sort(unique(groups))) <= 8) {
      meanplots[[count]] <- ggplot(means, aes(x=Sample, y=meanCI, color=CurrentGroup)) +  geom_pointrange(aes(ymin=LowerCI, ymax=UpperCI)) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=12))  +  xlab("Sample") + ylab(paste0("Mean Transition Probability (Bootstrapping) to ", MYgroup, "\n", sep="")) + labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    }
    else {
      meanplots[[count]] <- ggplot(means, aes(x=Sample, y=meanCI, color=CurrentGroup)) +  geom_pointrange(aes(ymin=LowerCI, ymax=UpperCI)) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=12))  +  xlab("Sample") + ylab(paste0("Mean Transition Probability (Bootstrapping) to ", MYgroup, "\n", sep="")) + labs(color='Current Group\n')
    }
    Mysums <- groupnamecols[,4:6]
    Mysums$CurrentGroup <- as.factor(metadata$Group)
    Mysums$Sample <- CIDF$Sample
    colnames(Mysums) <- c("LowerCI", "UpperCI", "meanCI", "CurrentGroup", "Sample")
    if (length(sort(unique(groups))) <= 8) {
      sumplots[[count]] <- ggplot(Mysums, aes(x=Sample, y=meanCI, color=CurrentGroup)) + geom_point(size = 3) +  geom_pointrange(aes(ymin=LowerCI, ymax=UpperCI)) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=17)) + ylab(paste0("Mean Transition Probability (Bootstrapping) to ", MYgroup, "\n", sep="")) + labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    }
    else {
      sumplots[[count]] <- ggplot(Mysums, aes(x=Sample, y=meanCI, color=CurrentGroup)) + geom_point(size = 3) +  geom_pointrange(aes(ymin=LowerCI, ymax=UpperCI)) + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), text = element_text(size=17)) + ylab(paste0("Mean Transition Probability (Bootstrapping) to ", MYgroup, "\n", sep="")) + labs(color='Current Group\n')
    }
    count = count + 1
  }
  if (length(sort(unique(groups))) <= 3) {
    png(meansfigname, width = 1080, height = 980)
    meansfigys <-  do.call("grid.arrange", c(meanplots, ncol=1))
    #print(meansfigys)
    dev.off()
    png(sumsfigname, width = 1080, height = 980)
    sumsfigys <- do.call("grid.arrange", c(sumplots, ncol=1))
    #print(sumsfigys)
    dev.off()
  }
  else {
    for (i in 1:length(sort(unique(groups)))) {
      meansfignametemp <- gsub("Confidence_Intervals", paste("Confidence_Intervals_",toString(sort(unique(groups))[i]), sep=""), meansfigname)
      sumsfignametemp <- gsub("Confidence_Intervals", paste("Confidence_Intervals_",toString(sort(unique(groups))[i]), sep=""), sumsfigname)
      png(meansfignametemp, width = 780, height = 580)
      meansfigys <-  do.call("grid.arrange", c(meanplots[i], ncol=1))
      #print(meansfigys)
      dev.off()
      png(sumsfignametemp, width = 780, height = 580)
      sumsfigys <- do.call("grid.arrange", c(sumplots[i], ncol=1))
      #print(sumsfigys)
      dev.off()
    }
  }
}

#InternalBootStrappingTransitonProbabilities(filename="Transition_Probability_NN5_UMAP.csv", selfremoved=FALSE,
                                           # location="/media/claired/DATADRIVE1/AsySymHc/TransitionProbabilityTest/images/testy/",
                                          #  metaGroup="/media/claired/DATADRIVE1/AsySymHc/TransitionProbabilityTest/metadata.csv",
                                          #  outputdir="/media/claired/DATADRIVE1/AsySymHc/TransitionProbabilityTest/images/testy/", grouptype="Strings", numberofIterations=100000)