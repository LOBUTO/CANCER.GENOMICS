#NEW CLAUDIN LOW DATASETS

#BUFFA - AFFIMETRIX


#ENERLY 
enerly.mrna<-fread("COLLABORATION/KANG/CL/ENERLY/mrna.txt", header=T)
setnames(enerly.mrna, c("ID", colnames(enerly.mrna)[2:ncol(enerly.mrna)]))
enerly.labels<-fread("COLLABORATION/KANG/CL/ENERLY/GPL6480-9577.txt", header=T,na.strings="")
enerly.labels<-enerly.labels[!is.na(GENE_SYMBOL),]
enerly.labels<-enerly.labels[!duplicated(GENE_SYMBOL),]
  
enerly.mrna<-merge(enerly.mrna, enerly.labels, by="ID")

enerly.mrna<-as.data.frame(enerly.mrna)
enerly.mrna$ID<-NULL
rownames(enerly.mrna)<-enerly.mrna$GENE_SYMBOL
enerly.mrna$GENE_SYMBOL<-NULL

enerly.mrna<-enerly.mrna[complete.cases(enerly.mrna),]
enerly.mrna<-data.matrix(enerly.mrna)
enerly.mrna<-normalizeBetweenArrays(enerly.mrna, method="quantile")

ENERLY.PRED<-Function.CL.PRED(CL.CENTROIDS, enerly.mrna, scale=F)
table(ENERLY.PRED$PRED)
