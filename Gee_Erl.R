# Gee_Erl.R
# Process Erlotinib clinical drug data from Geeleher 

library(data.table)
library(reshape2)
library(affy)
library(GEOquery)
library(ggplot2)
library(hugene10sthsentrezg.db)
library(pd.hugene.1.0.st.v1)

Function_cel_process <- function(geo.folder){
  geo.files<-paste(geo.folder, list.files(geo.folder), sep="/")
  
  raw.data=ReadAffy(filenames=geo.files)
  #raw.data=ReadAffy(verbose=TRUE, filenames=geo.files)
  raw.data@cdfName <- "hugene10st_Hs_ENTREZG"
  data.rma.norm=affy::rma(raw.data)
  rma=exprs(data.rma.norm)
  
  #Extract probe ids, entrez symbols, and entrez ids
  probes=row.names(rma)
  # Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
  Symbols = unlist(mget(probes, hugene10sthsentrezgSYMBOL, ifnotfound=NA))
  # Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))
  
  #Assign hugo symbols to matrix
  rma<-rma[!is.na(Symbols),]
  rownames(rma)<-Symbols[!is.na(Symbols)]
  
  #Obtain max of gene for more than one probe per gene
  rma<-data.table(rma,keep.rownames = T  )
  rma<-rma[,lapply(.SD, max), by=rn]
  rma<-data.frame(rma,row.names = 1)
  rma<-data.matrix(rma)
  
  #Clean headers
  colnames(rma)<-substr(colnames(rma), 1, 9)  
  
  #Return
  return(rma)
}

# Load data
date       <- Sys.Date()
nsclc      <- getGEO("GSE33072", getGPL = F)

# Process
exp_table  <- Function_cel_process("/Users/jzamalloa/Documents/Rotation/DATABASES/GEELEHER/ERLOTINIB_SORAFENIB/GSE33072_RAW")

pDataAll   <- data.table(pData(phenoData(nsclc[[1]])), stringsAsFactors = F)
p_erlo     <- pDataAll[,c("geo_accession", "characteristics_ch1.3", "characteristics_ch1.7"),with=F]
p_erlo     <- p_erlo[characteristics_ch1.3=="treatment: erlotinib"][grepl("progression-free",characteristics_ch1.7)]
p_erlo$Compound <- "Erlotinib"
p_erlo$target   <- as.numeric(sapply(as.character(p_erlo$characteristics_ch1.7), function(x) strsplit(x, ": ")[[1]][2]))
p_erlo$cell_name <- p_erlo$geo_accession
p_erlo     <- p_erlo[,c("cell_name", "Compound", "target"),with=F]

print(dim(exp_table))
print(exp_table[1:5,1:5])
print(p_erlo)

# Save files
saveRDS(list(exp_table = exp_table, feat_table = p_erlo), file="/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030417_GEE_ERLOTINIB.rds")

# Make some nice plots
pdf(file=paste0("~/Documents/FOLDER/LAB/GM/033017/",date,"_Byers_Kim_set_1.pdf"), height=8, width=8)
ggplot(p_erlo,
       aes(target)) +
  geom_histogram() +
  theme_bw() + ylab("Sample frequency") + xlab("Months to Progression") +
  ggtitle("Non-small cell lung cancer samples\nclassified by Progression-free survival time", subtitle = "Byers et al., 2103; Kim et al., 2011") +
  theme(axis.text  = element_text(size=12, colour="grey20"),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        strip.text.x = element_text(size = 12, colour = "black")) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

print("Done")