setwd("K:/projects/ChIP_seq/RNA_DATA/")
sample.name = "Sample_RH4Seq_T_D21KAACXX"

###load protein coding list
coding <- read.table("RSEM/HGNCprotein-coding_gene_Refseq.txt", sep="\t", header=T)

###load cufflinks (.TPM) output
EXP <- read.table(paste(sample.name,"/",sample.name,".gene.TPM.txt",sep=""), sep="\t", header=T)

###remove non-coding RNA entries
EXP$coding = EXP$GeneID %in% coding$symbol
EXP.coding = subset(EXP, EXP$coding %in% c("TRUE"))

##write files
write.table(EXP.coding, file=paste(sample.name,"/",sample.name,".gene.coding.TPM.txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
