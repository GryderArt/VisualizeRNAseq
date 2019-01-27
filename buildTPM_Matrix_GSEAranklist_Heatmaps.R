#build expression matrix from ChIP_seq/RNA_DATA/

### 1. Set directories, define samples to work on.  getwd()
    setwd("K:/projects/ChIP_seq/RNA_DATA/")
    sample.list.all = list.dirs(full.names = F, recursive = F)
    sample.list = sample.list.all[grep("shCtr|shScr|shSNAI2", sample.list.all)]
    #sample.list = sample.list[grep("MYOD_rep1|rep2_RNA_SRR5151763", sample.list)]
    project.folder = "RNA_projects"
    sample.set = "shSNAI2"
    
### 2. Make Protein coding TPM file

    lapply(sample.list, function(x) {
        
        ##load protein coding list
        coding <- read.table("RSEM/HGNCprotein-coding_gene_Refseq.txt", sep="\t", header=T)
      
        ##load RSEM (.TPM) output
        EXP <- read.table(paste(x,"/",x,".gene.TPM.txt",sep=""), sep="\t", header=T)
        
        ##remove non-coding RNA entries
        EXP$coding = EXP$GeneID %in% coding$symbol
        EXP.coding <<- subset(EXP, EXP$coding %in% c("TRUE"))
        
        ##write files
        write.table(EXP.coding[,1:5], file=paste(x,"/",x,".gene.coding.TPM.txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    })
    
### 3. Build matrix
    
    ###initiate matrix
    EXP.coding.matrix = EXP.coding[,1:4]
        
    ###loop through samples again, extract TPM
    lapply(sample.list, function(x) {
      
        ##read data, name column after "x"
        EXP.sample <- read.table(paste(x,"/",x,".gene.coding.TPM.txt",sep=""), sep="\t", header=T)
        
            EXP.sample = as.data.frame(EXP.sample[,5])
            removable.string = "Sample_"
            sample.name = gsub(removable.string,"",x)
            colnames(EXP.sample) = c(sample.name)
            
        #cbind writing outside the loop
        EXP.coding.matrix <<- cbind(EXP.coding.matrix, EXP.sample)
    })
    
    write.table(EXP.coding.matrix, file=paste(project.folder,"/ExpressionMatrices/",sample.set,".coding.TPM.matrix.txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
### 4. Make some comparison plots
    
    ###heatmaps
    cutoff.expression.min = 10
    EXP.coding.matrix$maxTPM = apply(EXP.coding.matrix[,c(5:ncol(EXP.coding.matrix))], 1, FUN=max)
    EXP.expressed.matrix = subset(EXP.coding.matrix, EXP.coding.matrix$maxTPM > cutoff.expression.min)
    library(pheatmap)
    pheatmap(log2(as.matrix(EXP.expressed.matrix[,5:(ncol(EXP.coding.matrix)-1)])+1))
    
    ###scatter
    library(LSD)
    comparisonplot(log2(EXP.expressed.matrix$J0621_T1R_T_HMFT5BGX7+1),log2(EXP.expressed.matrix$J0341_T1R_T_HMFT5BGX7+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    comparisonplot(log2(EXP.expressed.matrix$Fibroblast_MYOD_rep1_RNA_SRR5151765_T_GSE93263+1),log2(EXP.expressed.matrix$Fibroblast_rep2_RNA_SRR5151763_T_GSE93263+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    comparisonplot(log2(EXP.expressed.matrix$RH4_DMEM_T_HH3JVBGX7+1),log2(EXP.expressed.matrix$RH4_LW3_10_T_HH3JVBGX7+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    
### 5. GSEA ranklist maker
    
    ##Create a space for GSEA ranklists
    dir.create(file.path(project.folder, "GSEA_ranklist", sample.set))
    
    ##Add columnes for log2FC comparisons
    EXP.GSEA = EXP.coding.matrix
    
      EXP.GSEA$CTRshSNAI2vsScr = log2(EXP.GSEA$CTR_shSNAI2_10d_T_UTHSCSA+1)-log2(EXP.GSEA$CTR_shScr_10d_T_UTHSCSA+1)
      EXP.GSEA$RDshSNAI2vsScr = log2(EXP.GSEA$RD_shSNAI2_10d_T_UTHSCSA+1)-log2(EXP.GSEA$RD_shScr_10d_T_UTHSCSA+1)
      
    EXP.GSEA[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)]=round(EXP.GSEA[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)], digits = 5)
    
    for (i in (ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)) {
      Ranklist <- data.frame(EXP.GSEA[, 4])
      Ranklist$DeltaFPKM <- EXP.GSEA[, i]
      Ranklist = Ranklist[rev(order(Ranklist$DeltaFPKM)),]
      SampleName = colnames(EXP.GSEA)[i]
      mytime <- format(Sys.time(), "%b_%d_%Y")
      myfile <- file.path(project.folder, "GSEA_ranklist", sample.set,paste0(SampleName,"_",mytime,".rnk"))
      write.table(Ranklist, file = myfile, sep = "\t", row.names = FALSE, col.names = FALSE,
                  quote = FALSE, append = FALSE)
    }
    
### 6. Heatmap the data
 ## PREP DATA   
    GeneSet1.name = "BALL_GOI"
    #GeneSet1 = read.table(file = "Q:/RNAseq/Probe Lists/Qlucore_Format/SE_RA_Wave_Correlated_TFs.genelist.txt", sep="\t", header=F)
    GeneSet1 = as.data.frame(c("CREBBP","KMT2A","PAX5","MYC","KLF4","KLF2","LMNA"))  #custom cut
    EXP.GSEA.GeneSet1 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet1[,1])
    EXP.GSEA.GeneSet1.log2 = EXP.GSEA.GeneSet1; EXP.GSEA.GeneSet1.log2[,c(5,6,7,8,9,10,11,12,13)] = log2(EXP.GSEA.GeneSet1.log2[,c(5,6,7,8,9,10,11,12,13)]+1) #manual column selection 
      dir.create(file.path(project.folder,sample.set))
      myfile <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.txt"))
        write.table(EXP.GSEA.GeneSet1, file = myfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
      myfile.log2 <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.log2.txt"))
        write.table(EXP.GSEA.GeneSet1.log2, file = myfile.log2, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
      allfile <- file.path(project.folder,sample.set, paste0(sample.set,".allgenes.GSEA_matrix.txt"))
        write.table(EXP.GSEA, file = allfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)

 ## PLOT
    EXP.GSEA.GeneSet1.plot = EXP.GSEA.GeneSet1[,c(1,2,3,4,15,16)] #manual column selection  
    matrix.GeneSet1.plot = as.matrix(EXP.GSEA.GeneSet1.plot[,5:(ncol(EXP.GSEA.GeneSet1.plot))])
    rownames(matrix.GeneSet1.plot) = EXP.GSEA.GeneSet1.plot$GeneID
    #matrix.log.bounds = 2
    #matrix.GeneSet1.plot[matrix.GeneSet1.plot>matrix.log.bounds] <- matrix.log.bounds; matrix.GeneSet1.plot[matrix.GeneSet1.plot<(-matrix.log.bounds)] <- -matrix.log.bounds
    library(pheatmap)
    #breaklist = seq(-matrix.log.bounds, matrix.log.bounds, by = 0.1)
    pheatmap(log2(matrix.GeneSet1.plot+1),scale = "none", cluster_rows = T, cluster_cols = T)#, color = colorRampPalette(c("blue", "white", "orange")))
    
    
### 7. Rank scatter?
    
    library(ggplot2)
    ggplot(EXP.GSEA, aes(x=sgMYOD_cnt_CEM, y=sgNT_cnt_CEM))+geom_point(alpha=0.05,color="grey50")+geom_point(data=EXP.GSEA.GeneSet1, color="red")+theme_bw()
    ggplot(EXP.GSEA, aes(x=sgMYOD_cnt_CEM)) + geom_histogram()+theme_bw()
    
    
    library(LSD)
    sample.name = "sgMYOD_CEM" 
    comparisonplot(EXP.GSEA$sgMYOD_cnt_CEM,EXP.GSEA$sgNT_cnt_CEM, colpal="ylgnbu", main = sample.name, add.density = TRUE, ... = abline(h = 0.5),xlim = c(-4,4),ylim = c(-4,4))
    comparisonplot(EXP.GSEA.GeneSet1$sgMYOD_cnt_CEM,EXP.GSEA.GeneSet1$sgNT_cnt_CEM, colpal="ylgnbu", main = sample.name, add.density = TRUE, ... = abline(h = 0.5),xlim = c(-4,4),ylim = c(-4,4))
    
    comparisonplot(log2(EXP.GSEA$HEK_sgMYOD_cntrl_RNA_038_T_HMNTWBGX7+0.001),EXP.GSEA$sgMYOD_cnt_CEM, colpal=c("grey95","dodgerblue","red"), main = sample.name, add.density = TRUE, xlim = c(-2,12),ylim = c(-4,4))
    comparisonplot(log2(EXP.GSEA.GeneSet1$HEK_sgMYOD_cntrl_RNA_038_T_HMNTWBGX7+0.001),EXP.GSEA.GeneSet1$sgMYOD_cnt_CEM, colpal=c("grey10","red"), main = sample.name, add.density = TRUE, ... = abline(h = 0.5),xlim = c(-2,12),ylim = c(-4,4))
    comparisonplot(log2(EXP.GSEA.GeneSet1$HEK_sgNT_cntrl_RNA_038_T_HMNTWBGX7+0.001),EXP.GSEA.GeneSet1$sgNT_cnt_CEM, colpal=c("grey10","red"), main = sample.name, add.density = TRUE, ... = abline(h = 0.5),xlim = c(-2,12),ylim = c(-4,4))
    comparisonplot(EXP.GSEA$sgNT_cnt_CEM,EXP.GSEA$sgMYOD_cnt_CEM, colpal=c("grey80","red"), main = sample.name, add.density = TRUE, xlim = c(-4,4),ylim = c(-4,4))
    comparisonplot(EXP.GSEA.GeneSet1$sgNT_cnt_CEM,EXP.GSEA.GeneSet1$sgMYOD_cnt_CEM, colpal=c("grey10","red"), main = sample.name, add.density = TRUE, xlim = c(-4,4),ylim = c(-4,4))
    
    
    ranklists <- list.files(path=(file.path(project.folder, "GSEA_ranklist", sample.set)), full.names=T, recursive=F)  #DONT USE "." in GSEA Analysis Names, will cause incorrect parsing

    lapply(ranklists, function(x) {
          rankscat = read.table(x, sep="\t", header = F)
          #ranklistname = base(x)
          library(dplyr)
          rankscat = dplyr::mutate(rankscat, rank = row_number())
          rankscat.GeneSet1 = subset(rankscat, rankscat$V1 %in% GeneSet1[,1])
          library(ggplot2)
          rankplot = ggplot(rankscat, aes(x=rank, y=V2)) + geom_line(color="black") +
            geom_point(data=rankscat.GeneSet1,color = "red") +  
            scale_x_reverse() + theme_bw() + #scale_y_continuous(labels = comma,limits=c(0,70000)) +
            labs(y = "delta log2(TPM)") + ggtitle(x) +  
            geom_text(data=rankscat.GeneSet1,label=rankscat.GeneSet1$V1, hjust=1.2, size=5) + 
            theme(axis.text = element_text(size = 12), #aspect.ratio=1,
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()
            )
          
          output = paste(project.folder,basename(x),".ranklist.pdf",sep = '')
          ggsave(output,plot = rankplot,width = 7, height = 5)
    })
    
    
    
    