library("ape")

options(scipen=999)

gtfHs <- rtracklayer::import('Homo_sapiens.GRCh38.104_IntronsFinalSorted.gff3')
gtf_HsDf <- as.data.frame(gtfHs)
TranscriptInfo <- gtf_HsDf[gtf_HsDf$type == "transcript",]
TranscriptInfoRequired <- data.frame(TranscriptInfo$ID, TranscriptInfo$gene_id)
colnames(TranscriptInfoRequired) <- c("TranscriptID", "GeneID")

write.csv(TranscriptInfoRequired, "Gene_Transcript_Information.csv", row.names=FALSE)
