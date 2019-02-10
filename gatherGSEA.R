#########################################################################################################
##  Collect GSEA results.  Works with a limited set of gene sets all in the same GSEA output folder.
##
##  Berkley Gryder, 2018-2019 (gryderart@gmail.com): https://github.com/GryderArt/VisualizeRNAseq
#########################################################################################################

### 1. Find all the folder names
    setwd("K:/projects/ChIP_seq/RNA_DATA/RNA_projects/GSEA_ranklist/CRTF_HDACi/")
    GSEA.folder = "NatComm"
    samples <- list.dirs(path=GSEA.folder, full.names=F, recursive=F)  #DONT USE "." in GSEA Analysis Names, will cause incorrect parsing
    samples = samples[grepl("somestring", samples)]

### 2. Edit string from folder names to get report names (xls file name is the same as the sheet name, minus a few digits)
    sample.df = read.table(text = samples, sep=".")
    sample.df$V3 = as.character(sample.df$V3)
    colnames(sample.df)=c("Drug","GSEA_Style","Digits")
    sample.df$neg.xls.path = paste(GSEA.folder,"/",paste(sample.df$Drug,sample.df$GSEA_Style,sample.df$Digits,sep="."),"/gsea_report_for_na_neg_",sample.df$Digits,".xls",sep="")
    sample.df$pos.xls.path = paste(GSEA.folder,"/",paste(sample.df$Drug,sample.df$GSEA_Style,sample.df$Digits,sep="."),"/gsea_report_for_na_pos_",sample.df$Digits,".xls",sep="")
    
    #subset
    #sample.df = sample.df[grep("weight", sample.df$Drug), ]

### 3. Import GSEA data 
    #initiate dataframe
    allsamples = t(data.frame(sample = c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val","Drug")))
    colnames(allsamples) = c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val","Drug")
    allsamples=allsamples[-1,]
    
    lapply(sample.df$Drug, function(x) {
      #read in fake xls file
      temp.df <<- subset(sample.df,sample.df$Drug %in% x)
      df.pos = read.table(temp.df$neg.xls.path,sep="\t",header=T)
      df.neg = read.table(temp.df$pos.xls.path,sep="\t",header=T)
      df = rbind(df.pos,df.neg)
      #make it smaller, add column for identification
      df = df[,c("NAME","ES","NES","NOM.p.val","FDR.q.val","FWER.p.val")]
      df$Drug = x
      allsamples <<- rbind(allsamples,df)
    })
    
    #calculate some statistics
    allsamples$NOM.p.val[is.na(allsamples$NOM.p.val)] <- 1
    allsamples$NOM.p.val = allsamples$NOM.p.val + 0.00001
    allsamples$log10.NOM.p.val= log10(allsamples$NOM.p.val)
    
    #reorder by Enrichment Score (ES)
    #df$Species <- factor(df$Species, levels=df[order(df$x,decreasing=T),]$Species)  #not written yet
    ES.min = min(abs(allsamples$ES))
    NES.min = min(abs(allsamples$NES))
    
    #allsamples$NAME <- factor(allsamples$NAME, levels = c("HALLMARK_APOPTOSIS","GRYDER_HOUSEKEEPING","TFS_NO_EPIMACHINES","GRYDER_RH4_SE_GENES","GRYDER_RH4_TOP_SE_TFS")) #order genesets

### 4. Plots
  # a. Bubble chart plots
        library(ggplot2)
        ggplot(allsamples, aes(x = Drug, y = NAME, size = abs(NES)-NES.min,colour = NES)) + geom_point() +scale_colour_gradient2(low="tomato",mid="white", high="lightskyblue") +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ggplot(allsamples, aes(x = Drug, y = NAME, size = -log10.NOM.p.val,colour = NES)) + geom_point() +scale_colour_gradient2(low="tomato",mid="white", high="lightskyblue") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        ggplot(subset(allsamples,(allsamples$NAME %in% c("GRYDER_RH4_SELECTIVECRTFS","HALLMARK_APOPTOSIS"))), aes(x = Drug, y = NAME, size = -log10.NOM.p.val,colour,colour = NES)) + geom_point() +scale_colour_gradient2(low="tomato",mid="white", high="lightskyblue") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
  # b. GSEA ranklist plot for 1 gene set across many samples
        #extract out RANK and RUNNING ES data for a given gene set
        Geneset = "GRYDER_RH4_SELECTIVECRTFS"
        Geneset.allsamples = data.frame(RANK.IN.GENE.LIST=numeric(),RUNNING.ES=numeric(),Drug=character())
        
        lapply(sample.df$Drug, function(x) {
          #read in fake xls file
          ESplot.sample <- subset(sample.df,sample.df$Drug %in% x)
          
          ESplot.sample$geneset.path = paste(GSEA.folder,"/",paste(ESplot.sample$Drug,ESplot.sample$GSEA_Style,ESplot.sample$Digits,sep="."),"/",Geneset,".xls",sep="")
          
          Geneset.df = read.table(ESplot.sample$geneset.path,sep="\t",header=T)
         
          #make it smaller, add column for identification
          Geneset.df = Geneset.df[,c("RANK.IN.GENE.LIST","RUNNING.ES")]
          
          Geneset.df$Drug = x
          
          Geneset.allsamples <<- rbind(Geneset.allsamples,Geneset.df)
        })
        
        library(grid)
        library(gridExtra)
        library(wesanderson)
        Drugcolors = c("red","tomato","lightskyblue","goldenrod1","blue3")
        
        p1 = ggplot(Geneset.allsamples,aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES, group = Drug, color = Drug)) + geom_line(size=1.2) + geom_hline(yintercept=0,lty='dashed')+ scale_color_manual(values=Drugcolors)+
          theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = c(0.85, 0.85))+
          ylab("Enrichment Score") + xlab(paste("Genes Ranked by log2 fold change")) +ggtitle(paste("Geneset:", Geneset,sep=" "))
        
        p2 = ggplot(Geneset.allsamples,group = Drug) + stat_density(aes(x=RANK.IN.GENE.LIST,y=0.5,fill=..density..),geom="tile",position="identity")+facet_wrap(~Drug, ncol=1)+ scale_fill_gradient(low="white",high="darkgrey")+
          theme(axis.line=element_blank(),
                  axis.title.x=element_blank(), legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank())
          
        p3 = p2+ geom_linerange(aes(x=RANK.IN.GENE.LIST,ymin=0,ymax=1,color = Drug))+facet_wrap(~Drug, ncol=1)+theme_bw() + scale_color_manual(values=Drugcolors)+
          theme(axis.line=element_blank(),
                axis.title.x=element_blank(), legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank())
        
        grid.arrange(arrangeGrob(p1,p3,ncol=1))
