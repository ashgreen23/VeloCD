ConfidenceIntervalsAcrossTPNN <- function(wd, MainPattern, NN, embedType, outdir) {
  library(Rmisc); library(ggplot2) #; library("reshape2"); library(gridExtra); library(boot); library(tidyverse); 
  #get files
  setwd(wd)
  TP_files <- list.files(path = wd, pattern = MainPattern)
  TP_files <- grep(embedType, TP_files, value=TRUE)
  embedType <- gsub("Sample", "", embedType)
  TP_files <- grep(NN, TP_files, value=TRUE)
  #set up the table list
  TPNN1 <- read.csv(TP_files[1])
  rownames(TPNN1) <- TPNN1$Sample
  groups <- TPNN1[,2]
  grsam <- TPNN1[,c(1:2)]
  TPNN1 <- TPNN1[,-c(1:2)]
  my.names <- paste0(rownames(TPNN1), "_Table", sep="")
  mytables <- setNames(replicate(length(my.names), data.frame()), my.names)
  
  #open all files
  for (file in TP_files) {
    TPNN <- read.csv(file)
    rownames(TPNN) <- TPNN$Sample
    colnames(TPNN) <- gsub("X", "", colnames(TPNN))
    TPNN <- TPNN[,-c(1:2)]
    TPnnvalue <- gsub(MainPattern, "", file)
    TPnnvalue <- gsub(NN, "", TPnnvalue)
    TPnnvalue <- gsub(embedType, "", TPnnvalue)
    TPnnvalue <- gsub("_SampleGroup1.csv", "", TPnnvalue)
    TPnnvalue <- paste("TP_NN_", TPnnvalue, sep="")
    for (sampley in 1:nrow(TPNN)) {
      mytables[[sampley]] <- rbind(mytables[[sampley]], TPNN[sampley,])
      rownames(mytables[[sampley]]) <- c(rownames(mytables[[sampley]])[1:nrow(mytables[[sampley]])-1], TPnnvalue)
    }
  }
  setwd(outdir)
  #CI calculates column wise for each table
  CITab <- data.frame()
  for (table in 1:length(mytables)) {
    myTable <- mytables[[table]]
    CIs <- apply(as.matrix(myTable), 2, function(x) CI(x, ci=0.95)) #this is under a t-distribution
    #below agrees with online calculators:
    #CIs <- rbind(apply(as.matrix(myTable), 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[2,], apply(as.matrix(myTable), 2, function(x) mean(x)), apply(as.matrix(myTable), 2, function(x){mean(x)+c(-1.96,1.96)*sd(x)/sqrt(length(x))})[1,])
    #rownames(CIs) <- c("upper", "mean", "lower") #this is under a normal distribution, gives a less wide intrerval because the tails of a normal distribution are less wide than under t so 95% will be closer to the mean than under the t-distribution
    #smaller sample sizes so maybe the t-distribution is better?
    #the data may not even be anywhere normally distributed as for far away neighbours most will be zero then non-zero so left skew with one right hand tail no left hand tail as they can not be negative therefore might be useful to use an alternative method of measuring variance
    variance <- apply(as.matrix(myTable), 2, function(x) var(x))
    means <- apply(as.matrix(myTable), 2, function(x) mean(x))
    ranges <- apply(as.matrix(myTable), 2, function(x) range(x))
    sds <- apply(as.matrix(myTable), 2, function(x) sd(x))
    #Visualisation (CIs):
    mygroup <- grsam[grsam$Sample == gsub("_Table", "", names(mytables)[table]), 2]
    CIdf <- data.frame(colnames(CIs), CIs[1,], CIs[2,], CIs[3,])
    colnames(CIdf) <- c("Group", "Upper", "Mean", "Lower")
    CIdf$Group <- factor(CIdf$Group, levels=unique(CIdf$Group))
    #figname <- paste(gsub("SummaryTable", "CI_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", gsub("_Table", "", names(mytables)[table]),  ".png", sep="")
    #png(figname, width = 880, height = 680)
    #p <- ggplot(CIdf, aes(x=Group, y=Mean)) +  geom_pointrange(aes(ymin=Lower, ymax=Upper)) + theme_bw() + theme(text = element_text(size=17))  +  xlab("\nGroup") + ylab(paste0("Transition Probability of ", gsub("_Table", "", names(mytables)[table]), " (", toString(mygroup), ")", " to:", "\n", sep="")) #+ labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    #print(p)
    #dev.off()
    #Visualisation (Range)
    mygroup <- grsam[grsam$Sample == gsub("_Table", "", names(mytables)[table]), 2]
    Rangesdf <- data.frame(colnames(ranges), ranges[1,], means, ranges[2,])
    colnames(Rangesdf) <- c("Group", "Lower", "Mean", "Upper")
    Rangesdf$Group <- factor(Rangesdf$Group, levels=unique(Rangesdf$Group))
    #figname <- paste(gsub("SummaryTable", "Ranges_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", gsub("_Table", "", names(mytables)[table]),  ".png", sep="")
    #png(figname, width = 880, height = 680)
    #p <- ggplot(Rangesdf, aes(x=Group, y=Mean)) +  geom_pointrange(aes(ymin=Lower, ymax=Upper)) + theme_bw() + theme(text = element_text(size=17))  +  xlab("\nGroup") + ylab(paste0("Transition Probability of ", gsub("_Table", "", names(mytables)[table]), " (", toString(mygroup), ")", " to:", "\n", sep="")) #+ labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    #print(p)
    #dev.off()
    #Visualisation (variance)
    mygroup <- grsam[grsam$Sample == gsub("_Table", "", names(mytables)[table]), 2]
    Variancesdf <- data.frame(names(variance), variance)
    colnames(Variancesdf) <- c("Group", "Variance")
    Variancesdf$Group <- factor(Variancesdf$Group, levels=unique(Variancesdf$Group))
    #figname <- paste(gsub("SummaryTable", "Variance_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", gsub("_Table", "", names(mytables)[table]),  ".png", sep="")
    #png(figname, width = 880, height = 680)
    #p <- ggplot(Variancesdf, aes(x=Group, y=Variance)) + geom_point() + theme_bw() + theme(text = element_text(size=17))  +  xlab("\nGroup") + ylab(paste0("Transition Probability (Variance) of ", gsub("_Table", "", names(mytables)[table]), " (", toString(mygroup), ")", " to each group", "\n", sep="")) #+ labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    #print(p)
    #dev.off()
    #Visualisation (standard deviations)
    sdupper <- means + sds
    sdlower <- means - sds
    mygroup <- grsam[grsam$Sample == gsub("_Table", "", names(mytables)[table]), 2]
    StandardDevsdf <- data.frame(names(sds), sdupper, sdlower, means)
    colnames(StandardDevsdf) <- c("Group", "Upper", "Lower", "Means")
    StandardDevsdf$Group <- factor(StandardDevsdf$Group, levels=unique(StandardDevsdf$Group))
    #figname <- paste(gsub("SummaryTable", "StandardDeviations_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", gsub("_Table", "", names(mytables)[table]),  ".png", sep="")
    #png(figname, width = 880, height = 680)
    #p <- ggplot(StandardDevsdf, aes(x=Group, y=Means)) +  geom_pointrange(aes(ymin=Lower, ymax=Upper)) + theme_bw() + theme(text = element_text(size=17))  +  xlab("\nGroup") + ylab(paste0("Transition Probability of ", gsub("_Table", "", names(mytables)[table]), " (", toString(mygroup), ")", " to each group", "\n", sep="")) #+ labs(color='Current Group\n') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
    #print(p)
    #dev.off()
    
    #restructure this to put all the results on onw row then bind this to table with sample name as rownames
    Together <- c()
    for (row in 1:nrow(CIs)) {
      Together <- c(Together, CIs[row,])
    }
    names(Together) <- c(paste(colnames(CIs), "_upper", sep=""), paste(colnames(CIs), "_mean", sep=""), paste(colnames(CIs), "_lower", sep=""))
    CITab <- rbind(CITab, Together)
  }
  colnames(CITab) <- names(Together)
  rownames(CITab) <- rownames(TPNN1)
  Newfilename <- paste(gsub("SummaryTable", "CI_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, ".csv", sep="")
  write.csv(CITab, Newfilename, row.names=TRUE) #return as dataframe of sample rows, group cols with CIs in cells, save as csv
  
  #Visualisation
  cols <- c( "#0072B2", "#E69F00", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#000000", "#56B4E9")
  for (MyGroup in sort(unique(groups))) {
    mygroup2 <- c(paste(MyGroup, "_mean", sep=""), paste(MyGroup, "_upper", sep=""), paste(MyGroup, "_lower", sep=""))
    redTab  <- data.frame(rownames(CITab), groups, CITab[, colnames(CITab) %in% mygroup2])
    colnames(redTab)[1:2] <- c("Sample", "Group")
    colnames(redTab)[3:ncol(redTab)] <- c("Upper", "Mean", "Lower")
    redTab$Group <- factor(redTab$Group, levels=sort(unique(redTab$Group)))
    if (length(unique(groups)) <= 8) {
      figname <- paste(gsub("SummaryTable", "CI_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", MyGroup,  ".png", sep="")
      png(figname, width = 980, height = 680)
      p <- ggplot(redTab, aes(x=Sample, y=Mean, color=Group)) +  geom_pointrange(aes(ymin=Lower, ymax=Upper)) + theme_bw() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), text = element_text(size=17))  +  xlab("\nSample") + ylab(paste0("Transition Probability to Group: ", MyGroup, "\n", sep="")) + labs(color='Current Group') + scale_colour_manual(values = cols, aesthetics = c("colour", "fill"))
      print(p) 
      dev.off()
    }
    else {
      figname <- paste(gsub("SummaryTable", "CI_AcrossTPNN", MainPattern), "_", gsub("_", "", NN), "_", embedType, "_", MyGroup,  ".png", sep="")
      png(figname, width = 980, height = 680)
      p <- ggplot(redTab, aes(x=Sample, y=Mean, color=Group)) +  geom_pointrange(aes(ymin=Lower, ymax=Upper)) + theme_bw() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), text = element_text(size=17))  +  xlab("\nSample") + ylab(paste0("Transition Probability to Group: ", MyGroup, "\n", sep="")) + labs(color='Current Group')
      print(p)
      dev.off()
    }
  }
}
#ConfidenceIntervalsAcrossTPNN(wd="/media/claired/DATADRIVE1/Pvivax_New/TransitionProbWork/images/Pfal_test", MainPattern="SummaryTable_Mean_SelfNotRemoved", NN="_NN4_", embedType="tSNE")
#ConfidenceIntervalsAcrossTPNN(wd="/media/claired/DATADRIVE1/Pvivax_New/TransitionProbWork/images/Pfal_test", MainPattern="SummaryTable_Mean_SelfNotRemoved", NN="_NN4_", embedType="tSNESample")
#ConfidenceIntervalsAcrossTPNN(wd="/media/claired/DATADRIVE1/Pvivax_New/TransitionProbWork/images/NNchangetest_forCCVals", MainPattern="SummaryTable_Sum_SelfNotRemoved", NN="_NN4_", embedType="tSNE")
#ConfidenceIntervalsAcrossTPNN(wd="/media/claired/DATADRIVE1/Pvivax_New/TransitionProbWork/images/Pfal_Mezo_PfNF", MainPattern="SummaryTable_Sum_SelfNotRemoved", NN="_NN4_", embedType="tSNE")
#ConfidenceIntervalsAcrossTPNN(wd="/media/claired/DATADRIVE1/Pvivax_New/TransitionProbWork/images/NNchangetest_testCIacrossNN", MainPattern="SummaryTable_Sum_SelfNotRemoved", NN="_NN4_", embedType="tSNE")