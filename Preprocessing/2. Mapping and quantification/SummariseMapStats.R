library(reshape2); library(ggplot2); library(data.table)

setwd("") #location of the .final.out files from STAR
star_output_files <- list.files(path = "", pattern = "final.out") #same for path
print(length(star_output_files)) #check number of files
samples <- list()
input_reads <- list()
average_read_length <- list()
uniquely_mapped_reads <- list()
uniquely_mapped_reads_percentage <- list()
av_mapped_length <- list()
too_many_loci <- list()
too_many_loci_percent <- list()
unmapped_mismatched <- list()
unmapped_mismatched_percent <- list()
unmapped_short <- list()
unmapped_short_percent <- list()
unmapped_other <- list()
unmapped_other_percent <- list()
chimeric <- list()
chimeric_percent <- list()
for (star_out in star_output_files) {
  test <- read.csv(star_out, sep="\t")
  samples <- append(samples, sub("_R1Log.final.out", "", toString(star_out)))
  input_reads <- append(input_reads, as.numeric(toString(test[4,2])))
  average_read_length <- append(average_read_length, as.numeric(toString(test[5,2])))
  uniquely_mapped_reads_percentage <- append(uniquely_mapped_reads_percentage, as.numeric(sub('%', '', toString(test[8,2]))))
  uniquely_mapped_reads <- append(uniquely_mapped_reads, as.numeric(toString(test[7,2])))
  av_mapped_length <- append(av_mapped_length, as.numeric(toString(test[9,2])))
  too_many_loci <- append(too_many_loci, as.numeric(toString(test[24,2])))
  too_many_loci_percent <- append(too_many_loci_percent, as.numeric(sub('%', '', toString(test[25,2]))))
  unmapped_mismatched <- append(unmapped_mismatched, as.numeric(toString(test[27,2])))
  unmapped_mismatched_percent <- append(unmapped_mismatched_percent, as.numeric(sub('%', '', toString(test[28,2])))) 
  unmapped_short <- append(unmapped_short, as.numeric(toString(test[29,2])))
  unmapped_short_percent <- append(unmapped_short_percent, as.numeric(sub('%', '', toString(test[30,2]))))
  unmapped_other <- append(unmapped_other, as.numeric(toString(test[31,2])))
  unmapped_other_percent <- append(unmapped_other_percent, as.numeric(sub('%', '', toString(test[32,2])))) 
  chimeric <- append(chimeric, as.numeric(toString(test[34,2])))
  chimeric_percent <- append(chimeric_percent, as.numeric(sub('%', '', toString(test[35,2])))) ###
}
trimmed_stats <- cbind(samples, input_reads, average_read_length, uniquely_mapped_reads, uniquely_mapped_reads_percentage, av_mapped_length, too_many_loci, too_many_loci_percent, unmapped_mismatched, unmapped_mismatched_percent, unmapped_short, unmapped_short_percent, unmapped_other, unmapped_other_percent, chimeric, chimeric_percent)
rownames(trimmed_stats) <- samples

ID <- samples
needed_trimmed_stats <-  data.table(unlist(ID), unlist(uniquely_mapped_reads), unlist(too_many_loci),  unlist(unmapped_mismatched), unlist(unmapped_short), unlist(chimeric), unlist(unmapped_other))
rownames(needed_trimmed_stats) <- samples
colnames(needed_trimmed_stats) <- c("ID", "Uniquely Mapped", "Unmapped: Multimapped", "Unmapped: Mismatched",  "Unmapped: Too short", "Unmapped: Chimeric", "Unmapped: Other")
melted_trimmed_stats <- melt(needed_trimmed_stats, id.vars="ID", factorsAsStrings=F)
write.csv(needed_trimmed_stats, "MappedStats.csv", row.names=FALSE)

mean(needed_trimmed_stats$`Uniquely Mapped`) #
range(needed_trimmed_stats$`Uniquely Mapped`) #

# Stacked
ggplot(melted_trimmed_stats, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + theme_bw() + ylab("Number of read pairs\n") + xlab("\nSample ID") + labs(fill="Mapping State") + theme(axis.text.x=element_text(size=5)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))

ID <- samples
percent_trimmed_stats <-  data.table(unlist(ID), unlist(uniquely_mapped_reads_percentage), unlist(too_many_loci_percent),  unlist(unmapped_mismatched_percent), unlist(unmapped_short_percent), unlist(chimeric_percent), unlist(unmapped_other_percent))
rownames(percent_trimmed_stats) <- samples
colnames(percent_trimmed_stats) <- c("ID", "Uniquely Mapped", "Unmapped: Multimapped", "Unmapped: Mismatched",  "Unmapped: Too short", "Unmapped: Chimeric", "Unmapped: Other")
melted_percent_trimmed_stats <- melt(percent_trimmed_stats, id.vars="ID", factorsAsStrings=F)

# Stacked
ggplot(melted_percent_trimmed_stats, aes(fill=variable, y=value, x=ID)) + geom_bar(position="stack", stat="identity") + ylab("Percentage of read pairs\n") + xlab("\nSample ID") + labs(fill="Mapping State") + theme(axis.text.x=element_text(size=5)) + scale_y_continuous(expand = c(0,0)) + scale_fill_manual(values=c("#CC6677","#882255","#117733", "#44AA99", "#DDCC77", "#88CCEE"))
write.csv(percent_trimmed_stats, "PercentageMapped.csv", row.names=FALSE)

percent_trimmed_stats <- percent_trimmed_stats[percent_trimmed_stats$ID %in% MyMeta$lane_id,]
mean(percent_trimmed_stats$`Uniquely Mapped`) #
range(percent_trimmed_stats$`Uniquely Mapped`) #
