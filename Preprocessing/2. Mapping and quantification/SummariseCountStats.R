library(reshape2); library(ggplot2); library(data.table)

#################Gene-level expression################

genesum <- list.files(path = "", pattern = "_Gene.txt.summary") #add directory location
sampleName <- list()
assigned_reads <- list()
Unassigned_NoFeatures <- list()
Unassigned_Ambiguity <- list()
for (featurecount_out_summ in genesum) {
  nodups <- read.csv(featurecount_out_summ, sep="\t")
  sampleName <- append(sampleName, sub('_Gene.txt.summary', "", featurecount_out_summ))
  assigned_reads <- append(assigned_reads, nodups[1,2])
  Unassigned_NoFeatures <- append(Unassigned_NoFeatures, nodups[4,2])
  Unassigned_Ambiguity <- append(Unassigned_Ambiguity, nodups[2,2])
}
counts <- data.table(unlist(sampleName), unlist(assigned_reads), unlist(Unassigned_NoFeatures), unlist(Unassigned_Ambiguity))
colnames(counts) <- c("ID", "Assigned", "Unassigned: No Feature", "Unassigned: Ambiguous")
melted_counts <- melt(counts, id.vars="ID", factorsAsStrings=F)
ggplot(melted_counts, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + ylab("Read Pairs\n") + xlab("\nSample ID") + theme_bw() + labs(fill="Quantified") + theme(axis.text.x=element_text(size=4)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
counts$total <- rowSums(counts[,2:4])
counts$assigned_percent <- counts$Assigned/counts$total*100
counts$no_feature_percent <- counts$"Unassigned: No Feature"/counts$total*100
counts$ambiguous_percent <- counts$"Unassigned: Ambiguous"/counts$total*100
percentages <- data.table(counts$ID, counts$assigned_percent, counts$no_feature_percent, counts$ambiguous_percent)
colnames(percentages) <- c("ID", "Asigned", "No Feature", "Ambiguous")
melted_percentages <- melt(percentages, id.vars="ID", factorsAsStrings=F)
ggplot(melted_percentages, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + ylab("Read Pairs (%)\n") + xlab("\nSample ID") + theme_bw() + labs(fill="Quantified") + theme(axis.text.x=element_text(size=4)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
range(percentages$Asigned) #  
mean(percentages$Asigned) #
range(percentages$`No Feature`) # 
mean(percentages$`No Feature`) #
range(percentages$Ambiguous) # 
mean(percentages$Ambiguous) #
write.csv(percentages, "GeneSummaryStats_Percentages.csv", row.names=FALSE)
write.csv(counts, "GeneSummaryStats.csv", row.names=FALSE)

##############Exon (for each gene: "meta")-level expression################

genesum <- list.files(path = "", pattern = "_ExonMeta.txt.summary") #add directory location
sampleName <- list()
assigned_reads <- list()
Unassigned_NoFeatures <- list()
Unassigned_Ambiguity <- list()
for (featurecount_out_summ in genesum) {
  nodups <- read.csv(featurecount_out_summ, sep="\t")
  sampleName <- append(sampleName, sub('_ExonMeta.txt.summary', "", featurecount_out_summ))
  assigned_reads <- append(assigned_reads, nodups[1,2])
  Unassigned_NoFeatures <- append(Unassigned_NoFeatures, nodups[4,2])
  Unassigned_Ambiguity <- append(Unassigned_Ambiguity, nodups[2,2])
}
counts <- data.table(unlist(sampleName), unlist(assigned_reads), unlist(Unassigned_NoFeatures), unlist(Unassigned_Ambiguity))
colnames(counts) <- c("ID", "Assigned", "Unassigned: No Feature", "Unassigned: Ambiguous")
melted_counts <- melt(counts, id.vars="ID", factorsAsStrings=F)
ggplot(melted_counts, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + ylab("Read Pairs\n") + theme_bw() + xlab("\nSample ID") + labs(fill="Quantified") + theme(axis.text.x=element_text(size=2)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
counts$total <- rowSums(counts[,2:4])
counts$assigned_percent <- counts$Assigned/counts$total*100
counts$no_feature_percent <- counts$"Unassigned: No Feature"/counts$total*100
counts$ambiguous_percent <- counts$"Unassigned: Ambiguous"/counts$total*100
percentages <- data.table(counts$ID, counts$assigned_percent, counts$no_feature_percent, counts$ambiguous_percent)
colnames(percentages) <- c("ID", "Asigned", "No Feature", "Ambiguous")
melted_percentages <- melt(percentages, id.vars="ID", factorsAsStrings=F)
ggplot(melted_percentages, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + ylab("Read Pairs (%)\n") + theme_bw() + xlab("\nSample ID") + labs(fill="Quantified") + theme(axis.text.x=element_text(size=2)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
range(percentages$Asigned) #
mean(percentages$Asigned) #
range(percentages$`No Feature`) # 
mean(percentages$`No Feature`) #
range(percentages$Ambiguous) # 
mean(percentages$Ambiguous) #
write.csv(percentages, "ExonSummaryStats_Percentages.csv", row.names=FALSE)
write.csv(counts, "ExonSummaryStats.csv", row.names=FALSE)

#################Intron-(for each gene: "meta")-level expression################

genesum <- list.files(path = "", pattern = "_IntronMeta.txt.summary") #add directory location
sampleName <- list()
assigned_reads <- list()
Unassigned_NoFeatures <- list()
Unassigned_Ambiguity <- list()
for (featurecount_out_summ in genesum) {
  nodups <- read.csv(featurecount_out_summ, sep="\t")
  sampleName <- append(sampleName, sub('_IntronMeta.txt.summary', "", featurecount_out_summ))
  assigned_reads <- append(assigned_reads, nodups[1,2])
  Unassigned_NoFeatures <- append(Unassigned_NoFeatures, nodups[4,2])
  Unassigned_Ambiguity <- append(Unassigned_Ambiguity, nodups[2,2])
}
counts <- data.table(unlist(sampleName), unlist(assigned_reads), unlist(Unassigned_NoFeatures), unlist(Unassigned_Ambiguity))
colnames(counts) <- c("ID", "Assigned", "Unassigned: No Feature", "Unassigned: Ambiguous")
melted_counts <- melt(counts, id.vars="ID", factorsAsStrings=F)
ggplot(melted_counts, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + theme_bw() + ylab("Read Pairs\n") + xlab("\nSample ID") + labs(fill="Quantified") + theme(axis.text.x=element_text(size=2)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
counts$total <- rowSums(counts[,2:4])
counts$assigned_percent <- counts$Assigned/counts$total*100
counts$no_feature_percent <- counts$"Unassigned: No Feature"/counts$total*100
counts$ambiguous_percent <- counts$"Unassigned: Ambiguous"/counts$total*100
percentages <- data.table(counts$ID, counts$assigned_percent, counts$no_feature_percent, counts$ambiguous_percent)
colnames(percentages) <- c("ID", "Asigned", "No Feature", "Ambiguous")
melted_percentages <- melt(percentages, id.vars="ID", factorsAsStrings=F)
ggplot(melted_percentages, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + theme_bw() + ylab("Read Pairs (%)\n") + xlab("\nSample ID") + labs(fill="Quantified") + theme(axis.text.x=element_text(size=2)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
range(percentages$Asigned) # 
mean(percentages$Asigned) #
range(percentages$`No Feature`) #
mean(percentages$`No Feature`) #
range(percentages$Ambiguous) #  
mean(percentages$Ambiguous) #
write.csv(percentages, "IntronSummaryStats_Percentages.csv", row.names=FALSE)
write.csv(counts, "IntronSummaryStats.csv", row.names=FALSE)
