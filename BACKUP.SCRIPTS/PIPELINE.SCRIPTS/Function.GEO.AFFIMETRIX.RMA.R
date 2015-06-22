library(affy)
library("hugene10stv1cdf")
library(hugene10stprobeset.db)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(colorRamps)
library(oligo)

MATRIX.AFFY<-function(geo.folder){
  geo.files<-paste(geo.folder, list.files(geo.folder), sep="/")
  
  raw.data=read.celfiles(verbose=TRUE, filenames=geo.files)
  #raw.data=ReadAffy(verbose=TRUE, filenames=geo.files)
  data.rma.norm=rma(raw.data)
  rma=exprs(data.rma.norm)
  
  #Extract probe ids, entrez symbols, and entrez ids
  probes=row.names(rma)
  Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
  #Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
  Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))
  
  #Assign hugo symbols to matrix
  rma<-rma[!is.na(Symbols),]
  rownames(rma)<-Symbols[!is.na(Symbols)]
  
  #Obtain max of gene for more than one probe per gene
  rma<-data.table(rma,keep.rownames = T  )
  rma<-rma[,lapply(.SD, max), by=rn]
  rma<-data.frame(rma,row.names = 1)
  
  #Clean headers
  colnames(rma)<-substr(colnames(rma), 1, 9)  
  
  #Return
  return(rma)
}

#WORKING WITH HUANG GEO
geo.folder<-"DATABASES/METABOLOMICS/HUANG.2013/GSE48390_RAW"

anno.hs<-as.data.frame(org.Hs.egALIAS2EG)
setnames(anno.hs, c("EntrezGene.ID", "Hugo_Symbol"))

huang.subtypes<-intrinsic.cluster.predict(pam50.robust, t(scale(rma-apply(rma,1,median))), anno.hs, std="robust")
huang.subtypes<-as.data.table(huang.subtypes$subtype, keep.rownames=T)
setnames(huang.subtypes, c("SAMPLE", "SUBTYPE"))
huang.subtypes$SAMPLE<-substr(huang.subtypes$SAMPLE,1,10)

huang.true<-data.table(t(fread("DATABASES/METABOLOMICS/HUANG.2013/SUBTYPES.csv", header=F)))
setnames(huang.true, c("REAL.SUBTYPE", "SAMPLE"))
huang.true$REAL.SUBTYPE<-sapply(huang.true$REAL.SUBTYPE, function(x) unlist(strsplit(x, ": "))[2])
huang.true$REAL.SUBTYPE<-ifelse(huang.true$REAL.SUBTYPE=="Luminal A", "LumA",
                                ifelse(huang.true$REAL.SUBTYPE=="Luminal B", "LumB",
                                       ifelse(huang.true$REAL.SUBTYPE=="Normal breast-like", "Normal",
                                              ifelse(huang.true$REAL.SUBTYPE=="HER2-enriched", "Her2",
                                                     ifelse(huang.true$REAL.SUBTYPE=="Basal-like subtype", "Basal", huang.true$REAL.SUBTYPE)))))

huang.true<-merge(huang.true, huang.subtypes, by="SAMPLE")
sum(huang.true$REAL.SUBTYPE==huang.true$SUBTYPE)
huang.true[REAL.SUBTYPE!=SUBTYPE,]
#none=52
#median=52
#median+scale=53

#WORKING WITH TERUNUMA GEO
geo.folder<-"DATABASES/METABOLOMICS/TERUNUMA.2014/GSE37751_RAW"

#Assign normals and not
gsm.class<-data.table(t(fread("DATABASES/METABOLOMICS/TERUNUMA.2014/NORMAL.VS.TUMOR.GSM.csv", header=F, stringsAsFactors=F)))
setnames(gsm.class, c("CLASS", "SAMPLE"))
gsm.class$CLASS<-sapply(gsm.class$CLASS, function(x) unlist(strsplit(x, split=" "))[1])

#Store
saveRDS(list(MATRIX=data.matrix(rma), CLASS=gsm.class), "DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds")

#Plot
pheatmap(rma[c("PHGDH", "L2HGDH", "ADHFE1"),  ],scale="none", annotation_col=data.frame(CLASS=gsm.class$CLASS, row.names=gsm.class$SAMPLE), color=heat.colors(400))
pheatmap(log2(rma[c("PHGDH", "L2HGDH", "ADHFE1"), gsm.class[CLASS=="Tumor",]$SAMPLE] / apply(rma[c("PHGDH", "L2HGDH", "ADHFE1"), gsm.class[CLASS=="Normal",]$SAMPLE],1, median)),
         scale="none", color=heat.colors(400))

#WORKING WITH ROHLE GEO
geo.folder<-"DATABASES/METABOLOMICS/ROHLE.2013/GSE45197_RAW/"
saveRDS(rma, "DATABASES/METABOLOMICS/ROHLE.2013/051515.ROHLE.RMA.rds")
