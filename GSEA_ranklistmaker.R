###GSEA rank order list prep code
#Prepare expression matrix of delta DMSO
require(XLConnect)
wb = loadWorkbook("Expression_Matrix/ALL_RMS_CL/Treatments/EpiDrugs.all.PC_UCSC_Summary_deltaLog2.xlsx")
df = readWorksheet(wb, sheet = "deltaLog2", header = TRUE)

#or make a delta log2 matrix from cuff.frame
df = cuff.frame.coding[,1:(ncol(cuff.frame.coding)-1)]
df$RH4_DvEnt_1h_log2deltaFPKM = log(df$RH4_ENT1_T_HVNVFBGX2+1,2)-log(df$RH4_D1_T_HVNVFBGX2+1,2)
df$RH4_DvEnt_6h_log2deltaFPKM = log(df$RH4_ENT6_T_HVNVFBGX2+1,2)-log(df$RH4_D6_T_HVNVFBGX2+1,2)
df$RH4_DvEnt_24h_log2deltaFPKM = log(df$RH4_ENT24_T_HVNVFBGX2+1,2)-log(df$RH4_D24_T_HVNVFBGX2+1,2)
df = df[,c("gene_id","RH4_DvEnt_1h_log2deltaFPKM","RH4_DvEnt_6h_log2deltaFPKM","RH4_DvEnt_24h_log2deltaFPKM")]

df.coltron = subset(df,df$gene_id %in% c("MYOD1","SOX8","MYCN","MYC","PAX3","FOXO1","MYOG"))

df[,2:ncol(df)]=round(df[,2:ncol(df)], digits = 5)

for (i in 2:ncol(df)) {
  Ranklist <- data.frame(df[, 1])
  Ranklist$DeltaFPKM <- df[, i]
  Ranklist = Ranklist[rev(order(Ranklist$DeltaFPKM)),]
  SampleName = colnames(df)[i]
  mytime <- format(Sys.time(), "%b_%d_%H_%M_%Y")
  myfile <- file.path(paste0("GSEA/ranklists/",mytime, "_", SampleName, ".rnk"))
  write.table(Ranklist, file = myfile, sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE, append = FALSE)
}

