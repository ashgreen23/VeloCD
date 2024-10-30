TransitionTypeEstimates <- function(filename, location) {
  setwd(location)
  TPGroupTable <- read.csv(filename, check.names=FALSE)
  filename2 <- gsub(".*/", "", filename)
  #colnames(TPGroupTable)[3:ncol(TPGroupTable)] <- gsub("X", "", colnames(TPGroupTable)[3:ncol(TPGroupTable)]) #because X is added to the colnames of columns that are numeric in R
  TPGroupTable[is.na(TPGroupTable)] <- 0
  res <- apply(TPGroupTable[,3:ncol(TPGroupTable)],1,function(x) which(x==max(x)))
  SameGroup <- 0
  HigherGroup <- 0
  SmallerGroup <- 0
  if (typeof(res) == "integer") {
    HighestGroups <- colnames(TPGroupTable[,3:ncol(TPGroupTable)])[res]
    restable <- data.frame(TPGroupTable$Group, HighestGroups)
    colnames(restable) <- c("Group", "TP_HighestGroup")
    for (sample in 1:nrow(restable)) {
      myrow <- restable[sample,]
      if (as.integer(as.character(myrow$Group)) == as.integer(as.character(myrow$TP_HighestGroup))) {
        SameGroup <- SameGroup + 1
      }
      else if (as.integer(as.character(myrow$Group)) > as.integer(as.character(myrow$TP_HighestGroup))) {
        SmallerGroup <- SmallerGroup + 1
      }
      else if (as.integer(as.character(myrow$Group)) < as.integer(as.character(myrow$TP_HighestGroup))) {
        HigherGroup <- HigherGroup + 1
      }
    }
  }
  if (typeof(res) == "list") {
    selfGroups <- TPGroupTable$Group
    for (table in 1:length(res)) {
      mytab <- res[table]
      groups <- rownames(data.frame(mytab))
      selfGroupy <- selfGroups[table]
      if (length(groups) > 1) {
        for (group in groups) {
          if (as.integer(as.character(group)) == as.integer(as.character(selfGroupy))) {
            SameGroup <- SameGroup + 1/length(groups)
          }
          else if (as.integer(as.character(group)) < as.integer(as.character(selfGroupy))) {
            SmallerGroup <- SmallerGroup + 1/length(groups)
          }
          else if (as.integer(as.character(group)) > as.integer(as.character(selfGroupy))) {
            HigherGroup <- HigherGroup + 1/length(groups)
          }
        }
      }
      else {
        if (as.integer(as.character(groups)) == as.integer(as.character(selfGroupy))) {
          SameGroup <- SameGroup + 1
        }
        else if (as.integer(as.character(groups)) > as.integer(as.character(selfGroupy))) {
          HigherGroup <- HigherGroup + 1
        }
        else if (as.integer(as.character(groups)) < as.integer(as.character(selfGroupy))) {
          SmallerGroup <- SmallerGroup + 1
        }
      }
    }
  }
  SameGroupPer <- SameGroup/nrow(TPGroupTable)*100
  HigherGroupPer <- HigherGroup/nrow(TPGroupTable)*100
  SmallerGroupPer <- SmallerGroup/nrow(TPGroupTable)*100
  groupstable <- data.frame(c("Same Group", "Later Group", "Earlier Group"), c(SameGroup, HigherGroup, SmallerGroup), c(SameGroupPer, HigherGroupPer, SmallerGroupPer))
  colnames(groupstable) <- c("Probability of transition to:", "Number of samples", "% of samples")
  write.csv(groupstable, gsub("SummaryTable_", "", gsub(".csv", "_PercentagePredictionsTypes.csv", filename2)), row.names=FALSE)
}
#TransitionTypeEstimates(filename="SummaryTable_Mean_SelfRemoved_NN5_UMAP.csv", location="/media/claired/DATADRIVE1/CircadianMouseLiver/RF/Trimmed/Counts/cat/FateMappingWT/images/NoFS/")