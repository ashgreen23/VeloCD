Gene_Count_files <- list.files(path = "", pattern = "_Gene.txt") #path is the location of the count files

`%ni%` = Negate(`%in%`)

summary_files <- grep('summary', Gene_Count_files, value=TRUE)

 Gene_Count_files <- Gene_Count_files[Gene_Count_files %ni% summary_files]

length(Gene_Count_files)

setwd("") #location of the count files

GeneFile <- read.delim('', comment.char = '#') #the first sample in Gene_Count_files
GeneFileUseful <- GeneFile[,c(1,6:7)]
colnames(GeneFileUseful) <- c(colnames(GeneFileUseful)[1:2], "Count")
Gene_count_DataFrame <- GeneFileUseful
write.csv(GeneFileUseful[,1:2], "Lengths_Gene.csv", row.names=FALSE)
colnames(Gene_count_DataFrame) <- c("ID", "Length", "") #the final argument "" should be the same ID of the first sample file
Gene_count_DataFrame <- Gene_count_DataFrame[,-2]
GeneFileUseful <- NULL
for (Gene_Count_file in Gene_Count_files[2:length(Gene_Count_files)]) {
  GeneFile <- read.delim(Gene_Count_file, comment.char = '#')
  GeneFileUseful <- GeneFile[,c(1,7)]
  print(Gene_Count_file)
  print(dim(GeneFileUseful))
  print(identical(GeneFileUseful$Geneid, Gene_count_DataFrame$ID)) #check the Genes are in the same order
  colnames(GeneFileUseful) <- c("ID", "Count")
  Gene_count_DataFrame <- data.frame(Gene_count_DataFrame, GeneFileUseful$Count)
  new_col_name <- toString(Gene_Count_file)
  new_col_name <- paste0("Sample_", sub('_Gene.txt', "", new_col_name))
  print(new_col_name)
  colnames(Gene_count_DataFrame) <- c(colnames(Gene_count_DataFrame)[1:ncol(Gene_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
} 

write.csv(Gene_count_DataFrame, "FeatureCounts_Gene.csv", row.names=FALSE)

dim(Gene_count_DataFrame) #
