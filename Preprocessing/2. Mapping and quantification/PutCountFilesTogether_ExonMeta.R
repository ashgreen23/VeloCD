ExonMeta_Count_files <- list.files(path = "", pattern = "_ExonMeta.txt") #add directory location

`%ni%` = Negate(`%in%`)

summary_files <- grep('summary', ExonMeta_Count_files, value=TRUE)

ExonMeta_Count_files <- ExonMeta_Count_files[ExonMeta_Count_files %ni% summary_files]

length(ExonMeta_Count_files) #

setwd("") #add directory location

ExonMetaFile <- read.delim('', comment.char = '#') #add name of the first file in this vector
ExonMetaFileUseful <- ExonMetaFile[,c(1,6:7)]
colnames(ExonMetaFileUseful) <- c(colnames(ExonMetaFileUseful)[1:2], "Count")
ExonMeta_count_DataFrame <- ExonMetaFileUseful
write.csv(ExonMetaFileUseful[,1:2], "Lengths_ExonMeta.csv", row.names=FALSE)
colnames(ExonMeta_count_DataFrame) <- c("ID", "Length", "") #add sample ID of the sample in this first file
ExonMeta_count_DataFrame <- ExonMeta_count_DataFrame[,-2]
ExonMetaFileUseful <- NULL
for (ExonMeta_Count_file in ExonMeta_Count_files[2:length(ExonMeta_Count_files)]) {
  ExonMetaFile <- read.delim(ExonMeta_Count_file, comment.char = '#')
  ExonMetaFileUseful <- ExonMetaFile[,c(1,7)]
  print(ExonMeta_Count_file)
  print(dim(ExonMetaFileUseful))
  print(identical(ExonMetaFileUseful$Geneid, ExonMeta_count_DataFrame$ID)) #check the ExonMetas are in the same order
  colnames(ExonMetaFileUseful) <- c("ID", "Count")
  ExonMeta_count_DataFrame <- data.frame(ExonMeta_count_DataFrame, ExonMetaFileUseful$Count)
  new_col_name <- toString(ExonMeta_Count_file)
  new_col_name <- paste0("Sample_", sub('_ExonMeta.txt', "", new_col_name))
  print(new_col_name)
  colnames(ExonMeta_count_DataFrame) <- c(colnames(ExonMeta_count_DataFrame)[1:ncol(ExonMeta_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
}

write.csv(ExonMeta_count_DataFrame, "FeatureCounts_ExonMeta.csv", row.names=FALSE)

dim(ExonMeta_count_DataFrame) 
