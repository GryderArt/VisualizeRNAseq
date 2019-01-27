### Prep Matrix for EDEN


### 1. Set directories, define samples to work on.  getwd()
        setwd("K:/projects/ChIP_seq/RNA_DATA/")
        project.folder = "RNA_projects"
        sample.set = "shSNAI2"

### 2. Read in matrix from buildTPM_Matrix output
        TPMmatrix = read.table(file = "RNA_projects/ExpressionMatrices/shSNAI2.coding.TPM.matrix.txt",header = T, sep="\t")
        TPMmatrix$EDENstyle = paste(TPMmatrix$GeneID,"_",TPMmatrix$Chr,":",TPMmatrix$Start,"-",TPMmatrix$Stop,sep = "")
        
### 3. Transform into EDEN input "MATRIX" style
        EDENmatrix = TPMmatrix[,c(ncol(TPMmatrix),5:(ncol(TPMmatrix)-1))]
        write.table(EDENmatrix,file = "RNA_projects/ExpressionMatrices/shSNAI2.coding.EDEN.matrix.txt", col.names = T, row.names = F, sep = "\t", quote = F)

