### Correlation plot for matrix data


### 1. Set directories, define samples to work on.  getwd()
        setwd("K:/projects/ChIP_seq/RNA_DATA/RNA_projects/")
        project.folder = "SNAI2"

### 2. Read in matrix
        FPKMmatrix = read.table("SNAI2/SNAI2_MYOD1_correlation.txt", sep="\t",header = T)
        FPKMmatrix$NewGroups = FPKMmatrix$Group; FPKMmatrix$NewGroups <- sub("FP-", "", FPKMmatrix$NewGroups); FPKMmatrix$NewGroups <- sub("FN-", "", FPKMmatrix$NewGroups)
        FPKMmatrix$NewGroups <- sub("Normal Tissue - muscle", "RMS", FPKMmatrix$NewGroups)
        
### 3. Graph data
        library(LSD)
        sample.name = "SNAI2_v_GAPDH"
        output = paste(project.folder,sample.name,".pdf",sep = '')
        pdf(output) 
        comparisonplot(log2(FPKMmatrix$SNAI2+1),log2(FPKMmatrix$GAPDH+1), colpal="ylgnbu", main = sample.name, add.density = TRUE)#, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
        dev.off()
        
        library(ggplot2)
        ggplot(FPKMmatrix,aes(x=log2(SNAI2+1),y=log2(MYOD1+1),colour=factor(Group))) + geom_point() +
          theme_bw() + geom_smooth(method="lm")
        
### 4. Stats for correlations
        FPKMmatrix.sub = subset(FPKMmatrix, FPKMmatrix$NewGroups %in% c("RMS"))
        cor.test(log2(FPKMmatrix.sub$SNAI2+1), log2(FPKMmatrix.sub$GAPDH+1), method = "pearson", conf.level = 0.95)
        
        require(plyr)
        func <- function(FPKMmatrix)
        {
          return(data.frame(COR = cor(FPKMmatrix$SNAI2, FPKMmatrix$MYOD1)))
        }
        
        ddply(FPKMmatrix, .(NewGroups), func)
        
        
        
        
        