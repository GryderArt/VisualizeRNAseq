##################################################
###  build transcripts per million data matrix ###
###  then perform GSEA, heatmap comparisons,   ###
###  PCA plots, ranked scatter plots, and more ###
##################################################
###  by Berkley Gryder, gryderart@gmail.com    ###
##################################################

#build expression matrix from ChIP_seq/RNA_DATA/

### 1. Set directories, define samples to work on.  getwd()
    setwd("K:/projects/ChIP_seq/RNA_DATA/")
    sample.list.all = list.dirs(full.names = F, recursive = F)
    sample.list = sample.list.all[grep("LAmp|GSE43785", sample.list.all)]
    #sample.list = sample.list[grep("_sh", sample.list, invert = T)]
    project.folder = "RNA_projects"
    sample.set = "NB_RA"
    
### 2. Make Protein coding TPM file

    file.exists(paste(sample.list,"/",sample.list,".gene.TPM.txt",sep="")) #check to see if TPM.txt outputs exist
    
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
    EXP.expressed.matrix = EXP.coding.matrix
    EXP.expressed.matrix$maxTPM = apply(EXP.expressed.matrix[,c(5:ncol(EXP.expressed.matrix))], 1, FUN=max)
    EXP.expressed.matrix = subset(EXP.expressed.matrix, EXP.expressed.matrix$maxTPM > cutoff.expression.min)
    AllHumanTFs = read.table(file = "K:/projects/ChIP_seq/RNA_DATA/RNA_projects/Genesets/Qlucore_format/TranscriptionFactorsGryder.txt", sep="\t", header=F)
    EXP.expressed.matrix.TFs = subset(EXP.expressed.matrix, EXP.expressed.matrix$GeneID %in% AllHumanTFs$V1)
    library(pheatmap)
    pheatmap(log2(as.matrix(EXP.expressed.matrix[,5:(ncol(EXP.expressed.matrix.TFs)-1)])+1))
    
    ###matrix of scatters and Pearson Corr.
    library(ggplot2)
    library(GGally)
    GGally::ggpairs(log2((EXP.expressed.matrix[,5:(ncol(EXP.expressed.matrix.TFs)-1)])+1))
    
    ###scatter
    library(LSD)
    comparisonplot(log2(EXP.expressed.matrix$LAmp122_BG15n30uM_48h_RNA_T_H37TWBGXC+1),log2(EXP.expressed.matrix$LAmp122_DMSOhi_48h_RNA_T_H37TWBGXC+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    comparisonplot(log2(EXP.expressed.matrix$Fibroblast_MYOD_rep1_RNA_SRR5151765_T_GSE93263+1),log2(EXP.expressed.matrix$Fibroblast_rep2_RNA_SRR5151763_T_GSE93263+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    comparisonplot(log2(EXP.expressed.matrix$RH4_DMEM_T_HH3JVBGX7+1),log2(EXP.expressed.matrix$RH4_LW3_10_T_HH3JVBGX7+1), colpal="ylgnbu", main = sample.set, add.density = TRUE) #, ... = abline(h = 0.5),xlim = c(0,12),ylim = c(-5,5))
    
    ###PCA plot.  Metadata needs to be curated to meet your particular sample set

      EXP.pca = as.matrix(EXP.coding.matrix[,c(5:ncol(EXP.coding.matrix))])
      rownames(EXP.pca) = EXP.coding.matrix$GeneID
      EXP.pca.log2 = log2(EXP.pca+1)
      sample.name.list = c(colnames(EXP.pca))
      treatment.list = c("DMSO","DMSO","DMSO","DMSO","Rigo","Rigo","Rigo","Rigo")
      time.list = c("24 hrs","24 hrs","6 hrs","6 hrs","24 hrs","24 hrs","6 hrs","6 hrs")
      EXP.coldata = data.frame(sample.name.list,treatment.list,time.list)#;rownames(EXP.coldata) = sample.name.list
      
      # Generate PCA Data & Proportion of variability
      library(tidyverse) #CRAN - install.packages("tidyverse")
      library(ggrepel)   #CRAN - install.packages("ggrepel")
      pca        <- EXP.pca.log2 %>% t %>% prcomp
      EXP.d      <- pca$x %>% as.data.frame; EXP.d$sample.name.list = EXP.coldata$sample.name.list
      library(plyr)
      EXP.dmeta  <- join(EXP.d, EXP.coldata, by = "sample.name.list")
      pcv        <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

      # Graph the PCR plot
      plot.pca   <- ggplot(EXP.dmeta, aes(PC1,PC2,colour = treatment.list)) +
        geom_point() +
        xlab(label=paste0("PC1 (", pcv[1], "%)")) +
        ylab(label=paste0("PC2 (", pcv[2], "%)")) +
        theme_bw() +
        geom_label_repel(aes(label = sample.name.list), show.legend = F) +
        theme(axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15)) +
        labs(title    = "PCA",
             subtitle = "From RMS cell RNA-seq")
      print(plot.pca)
      
      
### 5. GSEA ranklist maker
    
    ##Create a space for GSEA ranklists
    dir.create(file.path(project.folder, "GSEA_ranklist", sample.set))
    
    ##Add columnes for log2FC comparisons
    EXP.GSEA = EXP.coding.matrix
    
      EXP.GSEA$L2FC_LAmp_LNCaP = log2(EXP.GSEA$LAmp122_DMSOhi_48h_RNA_T_H37TWBGXC+1)-log2(EXP.GSEA$LNCaP_ctrl_RNA_SRR653217_T_GSE43785+1)
      EXP.GSEA$L2FC_LAmp_LNCaPDHT = log2(EXP.GSEA$LAmp122_DMSOhi_48h_RNA_T_H37TWBGXC+1)-log2(EXP.GSEA$LNCaP_DHT_RNA_SRR653221_T_GSE43785+1)
      EXP.GSEA$L2FC_LAmp_BG15vEnza = log2(EXP.GSEA$LAmp122_BG15n30uM_48h_RNA_T_H37TWBGXC+1)-log2(EXP.GSEA$LAmp122_Enza30uM_48h_RNA_T_H37TWBGXC+1)
      
    EXP.GSEA[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)]=round(EXP.GSEA[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)], digits = 5)
    
    for (i in (ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)) {
      Ranklist <- data.frame(EXP.GSEA[, 4])
      Ranklist$DeltaTPM <- EXP.GSEA[, i]
      Ranklist = Ranklist[rev(order(Ranklist$DeltaTPM)),]
      SampleName = colnames(EXP.GSEA)[i]
      mytime <- format(Sys.time(), "%b_%d_%Y")
      myfile <- file.path(project.folder, "GSEA_ranklist", sample.set,paste0(SampleName,"_",mytime,".rnk"))
      write.table(Ranklist, file = myfile, sep = "\t", row.names = FALSE, col.names = FALSE,
                  quote = FALSE, append = FALSE)
    }
    
### 6. Heatmap the data
 ## PREP DATA   
    GeneSet1.name = "BG15n_v_Enza48_MTFs"
    GeneSet1 = read.table(file = "K:/projects/ChIP_seq/RNA_DATA/RNA_projects/Genesets/Qlucore_format/AnyCluster_MTF_top30.geneset.txt", sep="\t", header=F)
    #GeneSet1 = as.data.frame(c("MYC","MYCN","MYCL","MYOD1","SOX8","MYOG","PAX3","FOXM1","RARA","GAPDH","RPL11","POLR2L","HMGB1","CLIC1","ACTB","B2M"))  #custom cut
    EXP.GSEA.GeneSet1 = subset(EXP.GSEA, EXP.GSEA$GeneID %in% GeneSet1[,1])
    EXP.GSEA.GeneSet1.log2 = EXP.GSEA.GeneSet1;# EXP.GSEA.GeneSet1.log2[,c(5,6,7,8)] = log2(EXP.GSEA.GeneSet1.log2[,c(5,6,7,8)]+1) #manual column selection 
      dir.create(file.path(project.folder,sample.set))
      myfile <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.txt"))
        write.table(EXP.GSEA.GeneSet1, file = myfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
      myfile.log2 <- file.path(project.folder,sample.set, paste0(GeneSet1.name,".GSEA_matrix.log2.txt"))
        write.table(EXP.GSEA.GeneSet1.log2, file = myfile.log2, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)
      allfile <- file.path(project.folder,sample.set, paste0(sample.set,".allgenes.GSEA_matrix.txt"))
        write.table(EXP.GSEA, file = allfile, sep = "\t", row.names = F, col.names = T, quote = FALSE, append = FALSE)

 ## PLOT DATA
    #plot log2 scale, not delta values
    EXP.GSEA.GeneSet1.plot = EXP.GSEA.GeneSet1#[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20)] #manual column selection  
    matrix.GeneSet1.plot = as.matrix(EXP.GSEA.GeneSet1.plot[,5:(ncol(EXP.GSEA.GeneSet1.plot))])
    rownames(matrix.GeneSet1.plot) = EXP.GSEA.GeneSet1.plot$GeneID
    matrix.GeneSet1.plot.log2 = log2(matrix.GeneSet1.plot+1)
    matrix.log.bounds = 11
    matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot.log2>matrix.log.bounds] <- matrix.log.bounds; #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot<(-matrix.log.bounds)] <- -matrix.log.bounds
    library(pheatmap)
    #breaklist = seq(-matrix.log.bounds, matrix.log.bounds, by = 0.1)
    pheatmap(matrix.GeneSet1.plot.log2,scale = "none", cluster_rows = T, cluster_cols = F,main=paste(sample.set," ",GeneSet1.name," log2 TPM heatmap",sep=""))#, color = colorRampPalette(c("blue", "white", "orange")))
    
    #plot log 2 fold change values
    EXP.GSEA.GeneSet1.plot.log2fc = EXP.GSEA.GeneSet1[,(ncol(EXP.coding.matrix)+1):ncol(EXP.GSEA)] #manual column selection  
    matrix.GeneSet1.plot.log2fc = as.matrix(EXP.GSEA.GeneSet1.plot.log2fc)
    rownames(matrix.GeneSet1.plot.log2fc) = EXP.GSEA.GeneSet1$GeneID
    #matrix.log.bounds = 11
    #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot.log2>matrix.log.bounds] <- matrix.log.bounds; #matrix.GeneSet1.plot.log2[matrix.GeneSet1.plot<(-matrix.log.bounds)] <- -matrix.log.bounds
    library(pheatmap)
    #breaklist = seq(-matrix.log.bounds, matrix.log.bounds, by = 0.1)
    pheatmap(matrix.GeneSet1.plot.log2fc,scale = "none", cluster_rows = T, cluster_cols = T,main=paste(sample.set," ",GeneSet1.name," log2 TPM heatmap",sep=""))#, color = colorRampPalette(c("blue", "white", "orange")))
    
    
 ## Prettier heatmaps
    ### with improved clustering
    #sort it out a plot again. from http://slowkow.com/notes/heatmap-tutorial/
    
    mat_cluster_cols <- hclust(dist(t(matrix.GeneSet1.plot.log2)))
    plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
    
    library(dendsort)
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    mat_cluster_cols <- sort_hclust(mat_cluster_cols)
    plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
    mat_cluster_rows <- sort_hclust(hclust(dist(matrix.GeneSet1.plot)))
    
    pheatmap(matrix.GeneSet1.plot.log2,scale='none',cluster_cols=mat_cluster_cols,cluster_rows=T,  main=paste(sample.set," ",GeneSet1.name," log2 TPM heatmap",sep=""))
    
### 7. Rank scatter (used for Nature Biotechnology paper https://www.ncbi.nlm.nih.gov/pubmed/31712774)
    
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
    
    
    
    