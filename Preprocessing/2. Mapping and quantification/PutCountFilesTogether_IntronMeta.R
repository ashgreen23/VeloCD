IntronMeta_Count_files <- list.files(path = "", pattern = "_IntronMeta.txt") #add directory location

`%ni%` = Negate(`%in%`)

summary_files <- grep('summary', IntronMeta_Count_files, value=TRUE)

IntronMeta_Count_files <- IntronMeta_Count_files[IntronMeta_Count_files %ni% summary_files]

length(IntronMeta_Count_files) #

setwd("") #add directory location

IntronMetaFile <- read.delim('', comment.char = '#') #first file in the vector to initialize the data frame and get the  lengths
IntronMetaFileUseful <- IntronMetaFile[,c(1,6:7)]
colnames(IntronMetaFileUseful) <- c(colnames(IntronMetaFileUseful)[1:2], "Count")
IntronMeta_count_DataFrame <- IntronMetaFileUseful
write.csv(IntronMetaFileUseful[,1:2], "Lengths_IntronMeta.csv", row.names=FALSE)
colnames(IntronMeta_count_DataFrame) <- c("ID", "Length", "") #missing argument is the sample id of the first sample
IntronMeta_count_DataFrame <- IntronMeta_count_DataFrame[,-2]
IntronMetaFileUseful <- NULL
for (IntronMeta_Count_file in IntronMeta_Count_files[2:length(IntronMeta_Count_files)]) {
  IntronMetaFile <- read.delim(IntronMeta_Count_file, comment.char = '#')
  IntronMetaFileUseful <- IntronMetaFile[,c(1,7)]
  print(IntronMeta_Count_file)
  print(dim(IntronMetaFileUseful))
  print(identical(IntronMetaFileUseful$Geneid, IntronMeta_count_DataFrame$ID)) #check the IntronMetas are in the same order
  colnames(IntronMetaFileUseful) <- c("ID", "Count")
  IntronMeta_count_DataFrame <- data.frame(IntronMeta_count_DataFrame, IntronMetaFileUseful$Count)
  new_col_name <- toString(IntronMeta_Count_file)
  new_col_name <- paste0("Sample_", sub('_IntronMeta.txt', "", new_col_name))
  print(new_col_name)
  colnames(IntronMeta_count_DataFrame) <- c(colnames(IntronMeta_count_DataFrame)[1:ncol(IntronMeta_count_DataFrame)-1], new_col_name)
  print("-------------------------------------")
}

write.csv(IntronMeta_count_DataFrame, "FeatureCounts_IntronMeta.csv", row.names=FALSE)

dim(IntronMeta_count_DataFrame)
