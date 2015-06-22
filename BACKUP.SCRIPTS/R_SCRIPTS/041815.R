#041815.R
#Analyze expression sources in differential expression
library(edgeR)

Function.brca.clinical<-function(clinical.file){
  #Process clinical file
  clinical<-fread(clinical.file, skip=3, header=F, stringsAsFactors=F, drop=c(1,3:43,45:109))
  setnames(clinical, c("SAMPLE", "ER"))
  clinical<-clinical[ER %in% c("Negative", "Positive"),]
  clinical$SAMPLE<-gsub("-",".",clinical$SAMPLE)
  clinical<-rbind(clinical[,list(SAMPLE=paste(SAMPLE, "01A", sep=".")  ), by="ER"], clinical[,list(SAMPLE=paste(SAMPLE, "01B", sep=".")  ), by="ER"])
  return(clinical)
}

Function.scale.exp.er<-function(brca.exp, brca.clinical){
  
  #Assign clinical info to samples
  clin<-brca.clinical[SAMPLE %in% colnames(brca.exp),]
  
  #Calculate imbalance
  POS<-table(clin$ER)[["Positive"]]
  NEG<-table(clin$ER)[["Negative"]]
  
  if (POS<NEG){
    #Then pos is limiting, take negative limited sample
    POS.SAMPLES<-clin[ER=="Positive",]$SAMPLE
    print (length(POS.SAMPLES))
    NEG.SAMPLES.LIM<-sample(clin[ER=="Negative",]$SAMPLE, POS)
    print (length(NEG.SAMPLES.LIM))
    
    #Calculate median based on mixture of balanced population
    gene.medians<-apply(brca.exp[,c(POS.SAMPLES, NEG.SAMPLES.LIM)], 1, median)
    
  } else{
    #Otherwise neg is limiting
    NEG.SAMPLES<-clin[ER=="Negative",]$SAMPLE
    print (length(NEG.SAMPLES))
    POS.SAMPLES.LIM<-sample(clin[ER=="Positive",]$SAMPLE, NEG)
    print (length(POS.SAMPLES.LIM))
    
    gene.medians<-apply(brca.exp[,c(NEG.SAMPLES, POS.SAMPLES.LIM)], 1, median)
  }
  
  #Scale based on corrected subtype mixture
  brca.exp<-scale(brca.exp - gene.medians)
  
  #Return
  return(brca.exp)
}

brca.clinical<-Function.brca.clinical("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab//nationwidechildrens.org_clinical_patient_brca.txt")
######Expression of Microarray data in breast cancer from agilent platform#####
brca.agilent<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041815.BRCA.MATRICES.AGILENT.rds")
data.matrix(brca.agilent$tumor)[1:3,1:3]

require(limma)
combined.matrices<-cbind(brca.agilent$normal, brca.agilent$tumor[rownames(brca.agilent$normal),])
combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
combined.matrices<-data.matrix(combined.matrices)
combined.matrices[1:3,1:3]

design.matrix = data.frame(type=factor(c(rep("normal", ncol(brca.agilent$normal)), rep("cancer", ncol(brca.agilent$tumor)))) )
design.matrix = model.matrix(~type, design.matrix)
colnames(design.matrix)<- c("normal", "tumor")

combined.matrices = normalizeBetweenArrays(combined.matrices, method="quantile")

fit<-lmFit(combined.matrices, design.matrix)
eb = eBayes(fit)
volcanoplot(eb, highlight=20)
exp.table<-data.table(topTable(eb, coef=2, n=Inf),keep.rownames=T)
setnames(exp.table, c("Hugo_Symbol", setdiff(colnames(exp.table), "rn")))
ma<-cbind(exp.table$logFC, exp.table$AveExpr)
plotMA(ma)

exp.table[adj.P.Val<0.05,]
exp.table[adj.P.Val<0.05 & abs(logFC)>1,]

BRCA.MA.RAW.TCGA.PRED<-Function.CL.PRED(CL.CENTROIDS, Function.scale.exp.er(data.matrix(brca.agilent$tumor), brca.clinical), scale=F)
PREDICTED.CL.SAMPLES<-BRCA.MA.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE
BRCA.MA.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]
#####Expression of raw counts from RNASEQ V1 in breast cancer######
brca.V1.RAW<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041815.BRCA.RNASEQ.MATRICES.V1.RAW.rds")
dim(brca.V1.RAW$tumor)

require(DESeq2)
library("BiocParallel")

combined.matrices<-cbind(brca.V1.RAW$normal, brca.V1.RAW$tumor[rownames(brca.V1.RAW$normal),])
combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
combined.matrices<-data.matrix(combined.matrices)

design.matrix = data.frame(type=as.factor(c(rep("normal", ncol(brca.V1.RAW$normal)), rep("cancer", ncol(brca.V1.RAW$tumor)))) )
design.matrix = model.matrix(~type, design.matrix)

#diff.exp
G.all<-copy(combined.matrices)
isexpr<-rowSums(G.all)>=20
G.all<-G.all[isexpr,]
G.all<-DGEList(counts=G.all) #For scale normalization
G.all<-calcNormFactors(G.all)
G.all<-voom(G.all,design.matrix, normalize.method = "quantile")

BRCA.V1.RAW.TO.MA.FIT<-lmFit(G.all, design.matrix)
BRCA.V1.EB<-eBayes(BRCA.V1.RAW.TO.MA.FIT)

hist(BRCA.V1.EB$p.value[, 2])
volcanoplot(BRCA.V1.EB, highlight=20)
brca.top.table<-topTable(BRCA.V1.EB, coef=2, number=25000)
ma<-cbind(brca.top.table$logFC, brca.top.table$AveExpr)
plot(ma)

data.table(brca.top.table, keep.rownames = T)[abs(logFC)>1 & adj.P.Val<0.05,]

#####Expression of RPKM from V1#####
brca.RPKM<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/050415.BRCA.RNASEQ.MATRICES.V1.rds")

BRCA.RPKM.TCGA.PRED<-Function.CL.PRED(CL.CENTROIDS, Function.scale.exp.er(data.matrix(brca.RPKM$tumor), brca.clinical), scale=F)
PREDICTED.CL.SAMPLES<-BRCA.RPKM.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE

#####Expression of "raw counts" (RSEM non-upper quartile normalized) from RNASEQ V2 in breast cancer######
brca.V2.RAW<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041715.BRCA.RNASEQ.MATRICES.V2.RAW.rds")

"28: TCGA.A8.A08H.01A 157.9823    158.6941 CLAUDIN.LOW
29: TCGA.A2.A25F.01A 151.7387    152.7063 CLAUDIN.LOW
30: TCGA.AC.A2QH.01A 163.0151    166.2887 CLAUDIN.LOW"

BRCA.V2.RAW.PREDICT<-Function.CL.PRED(CL.CENTROIDS, Function.scale.exp.er(data.matrix(brca.V2.RAW$tumor), brca.clinical), scale=F)
PREDICTED.CL.SAMPLES<-BRCA.V2.RAW.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE

#####Expression of RSEM V2 from Breast cancer#####
BRCA.V2.RSEM<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041715.BRCA.RNASEQ.MATRICES.V2.RSEM.UQ.rds")

BRCA.V2.RSEM.PREDICT<-Function.CL.PRED(CL.CENTROIDS, Function.scale.exp.er(data.matrix(BRCA.V2.RSEM$tumor), brca.clinical), scale=F)
PREDICTED.CL.SAMPLES<-BRCA.V2.RSEM.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE
write.table(BRCA.V2.RSEM.PREDICT, "COLLABORATION/KANG/CL/052015.TCGA.RSEM.SUBTYPES.txt", quote = F,sep = "\t",col.names = T,row.names = F)

intersect(BRCA.V2.RSEM.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE, BRCA.RPKM.TCGA.PRED[PRED=="CLAUDIN.LOW",]$SAMPLE)
intersect(BRCA.V2.RSEM.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE,BRCA.MA.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE)
intersect(BRCA.V2.RSEM.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE,BRCA.V1.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE)
intersect(BRCA.V2.RSEM.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE,BRCA.V2.RAW.PREDICT[PRED=="CLAUDIN.LOW",]$SAMPLE)

##################################################

#####Expression of raw counts from miRNASEQ in breast cancer######
brca.mirnaseq.iga<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041915.BRCA.MIRNASEQ.ILLUMINAGA.MATRICES.RAW.rds")
brca.mirnaseq.iga$tumor[1:3,1:3]

brca.mirnaseq.ihiseq<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041915.BRCA.MIRNASEQ.ILLUMINAHISEQ.MATRICES.RAW.rds")
dim(brca.mirnaseq.ihiseq$tumor)

BRCA.V2.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE %in% colnames(brca.mirnaseq.iga$tumor)
BRCA.V2.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE %in% colnames(brca.mirnaseq.ihiseq$tumor)

#Try with hiseq first
mirna.rnaseq.to.voom<-function(mirna.objs){
  
  main.list<-lapply(mirna.objs, function(i) {
    i[["normal"]]<-i[["normal"]][,unique(colnames(i[["normal"]]))]
    i[["tumor"]]<-i[["tumor"]][,unique(colnames(i[["tumor"]]))]
    combined.matrices<-cbind(i[["normal"]], i[["tumor"]][rownames(i[["normal"]]),])
    print (dim(combined.matrices))
    combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
    combined.matrices<-data.matrix(combined.matrices)
    
    design.matrix = data.frame(type=c(rep("normal", ncol(i[["normal"]])), rep("cancer", ncol(i[["tumor"]]))) )
    design.matrix = model.matrix(~type, design.matrix)
    
    G.isexpr<- rowSums(cpm(combined.matrices)>1) >= 20 #Keep genes with at least n count-per-million reads (cpm) in at least 20 samples
    combined.matrices<-combined.matrices[G.isexpr,]
    
    mirna.ma.matrix<-voom(combined.matrices, design.matrix, normalize.method="quantile")    
    
    return(mirna.ma.matrix$E)
  })
  
  #Find common mirna in all lists
  common.mirna<-intersect(rownames(main.list[[1]]), rownames(main.list[[2]]))
  
  #Filter for unique samples
  first.samples<-colnames(main.list[[1]])
  second.samples<-setdiff(colnames(main.list[[2]]), colnames(main.list[[1]]))
  
  #Build final matrix
  main.matrix<-cbind(main.list[[1]][common.mirna,], main.list[[2]][common.mirna, second.samples])
  return(main.matrix)
}

normalized.mirna.matrix<-mirna.rnaseq.to.voom( list(brca.mirnaseq.ihiseq, brca.mirnaseq.iga))

#PREDICTED.CL.SAMPLES<-c(PREDICTED.CL.SAMPLES, paste(substr(PREDICTED.CL.SAMPLES, 1,12), ".01B", sep="") )
normalized.mirna.matrix<-normalized.mirna.matrix[,intersect(colnames(BRCA.V2.RSEM$tumor), colnames(normalized.mirna.matrix))]

mirna.hiseq.cl.samples<-intersect(PREDICTED.CL.SAMPLES, colnames(normalized.mirna.matrix))
mirna.cl.exp.melt<-data.table(melt(normalized.mirna.matrix[c("hsa-mir-199a-1", "hsa-mir-199a-2"),!(grepl("11A", colnames(normalized.mirna.matrix)))]))
setnames(mirna.cl.exp.melt, c("miRNA", "SAMPLE", "EXP"))
mirna.cl.exp.melt<-mirna.cl.exp.melt[!(grepl("11B", SAMPLE)),]

mirna.cl.exp.melt$TYPE<-mirna.cl.exp.melt$SAMPLE %in% mirna.hiseq.cl.samples
mirna.cl.exp.melt[order(EXP, decreasing=T),][miRNA=="hsa-mir-199a-1",]

ggplot(mirna.cl.exp.melt, aes(miRNA, EXP, colour=TYPE)) + geom_boxplot() + theme.format +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

wilcox.test(mirna.cl.exp.melt[TYPE==TRUE & miRNA=="hsa-mir-199a-1",]$EXP,
            mirna.cl.exp.melt[TYPE==FALSE & miRNA=="hsa-mir-199a-1",]$EXP,
            alternative="greater", paired=F, var.equal=F)
wilcox.test(mirna.cl.exp.melt[TYPE==TRUE & miRNA=="hsa-mir-199a-2",]$EXP,
            mirna.cl.exp.melt[TYPE==FALSE & miRNA=="hsa-mir-199a-2",]$EXP,
            alternative="greater", paired=F)
#

#Obtain p-values for all miRNA
mirna.cl.exp.melt<-data.table(melt(normalized.mirna.matrix))
setnames(mirna.cl.exp.melt, c("miRNA", "SAMPLE", "EXP"))
mirna.cl.exp.melt$TYPE<-mirna.cl.exp.melt$SAMPLE %in% mirna.hiseq.cl.samples
mirna.cl.exp.melt<-mirna.cl.exp.melt[!(grepl("11A", SAMPLE)),][!(grepl("11B", SAMPLE)),]

mirna.pvalue<-sapply(as.vector(unique(mirna.cl.exp.melt$miRNA)), function(x) {
  main.table<-mirna.cl.exp.melt[miRNA==x,]
  p.value<-wilcox.test(main.table[TYPE=="TRUE",]$EXP, main.table[TYPE=="FALSE"]$EXP, alternative="greater")$p.value
  return(p.value)
})
mirna.pvalue<-data.table(miRNA=as.vector(unique(mirna.cl.exp.melt$miRNA)), P.VAL=mirna.pvalue)
mirna.pvalue$P.VAL.ADJ<-p.adjust(mirna.pvalue$P.VAL, method="fdr")
mirna.pvalue<-mirna.pvalue[order(P.VAL.ADJ),]
write.table(mirna.pvalue, "COLLABORATION/KANG/CL/052015.BRCA.TCGA.PEROU.PVALS.txt", quote=F,sep="\t",row.names=F, col.names=T)
mirna.pvalue[P.VAL.ADJ<0.05,][1:50,]
mirna.pvalue[miRNA=="hsa-mir-199a-1",]

###Classify subtypes through genefu and PAM50###
library(genefu)
library(org.Hs.eg.db)
data(pam50)
anno.hs<-as.data.frame(org.Hs.egALIAS2EG)
setnames(anno.hs, c("EntrezGene.ID", "Hugo_Symbol"))

pam50$centroids

brca.subtypes<-intrinsic.cluster.predict(pam50, t(Function.scale.exp.er(data.matrix(BRCA.V2.RSEM$tumor), brca.clinical)), anno.hs)
brca.subtypes<-as.data.table(brca.subtypes$subtype,keep.rownames=T)
setnames(brca.subtypes, c("SAMPLE", "SUBTYPE"))

###Classify all into PAM50 subtypes and CL subtypes####
#brca.subtypes<-brca.subtypes[grepl("01A", SAMPLE),]
brca.subtypes<-brca.subtypes[!(grepl("11A", SAMPLE)),][!(grepl("11B", SAMPLE)),]
brca.subtypes$SUBTYPE<-ifelse(brca.subtypes$SAMPLE %in% PREDICTED.CL.SAMPLES, "CL", brca.subtypes$SUBTYPE)
write.table(brca.subtypes, "COLLABORATION/KANG/CL/052015.RSEM.INTRINSIC.SUBTYPES.txt", quote=F, sep="\t", row.names=F)
ggplot(brca.subtypes, aes(SUBTYPE)) + geom_histogram() + theme.format + ylab("Samples") + xlab("Subtypes") + 
  ggtitle("TCGA Samples - Breast Cancer Intrinsic Subtypes")

brca.subtypes<-merge(brca.subtypes, mirna.cl.exp.melt, by="SAMPLE")
ggplot(brca.subtypes[miRNA %in% c("hsa-mir-199a-2", "hsa-mir-199a-1"),], aes(miRNA, EXP, colour=SUBTYPE)) + geom_boxplot() + theme.format

##LCOR analysis
LCOR.RSEM<-data.table(t(data.frame(scale(BRCA.V2.RSEM$tumor- apply(BRCA.V2.RSEM$tumor, 1, median)))["LCOR",]), keep.rownames = T)
setnames(LCOR.RSEM, c("SAMPLE", "LCOR"))
LCOR.RSEM<-merge(LCOR.RSEM, brca.subtypes, by="SAMPLE")
ggplot(LCOR.RSEM, aes(SUBTYPE, LCOR)) + geom_boxplot() + theme.format

###Try Enerly again
enerly.mrna<-fread("COLLABORATION/KANG/CL/ENERLY/mrna.txt", header=T)
enerly.probes<-fread("COLLABORATION/KANG/CL/ENERLY/GPL6480-9577.txt", header=T, na.strings = "")
enerly.probes<-enerly.probes[!is.na(GENE_SYMBOL),]
setnames(enerly.probes, c("ID_REF", "GENE"))
enerly.mrna<-merge(enerly.mrna, enerly.probes, by="ID_REF")
enerly.mrna$ID_REF<-NULL
enerly.mrna<-enerly.mrna[,lapply(.SD, max), by="GENE"]
enerly.mrna<-data.frame(enerly.mrna,row.names = 1)
enerly.mrna<-data.matrix(enerly.mrna)
dim(enerly.mrna)
enerly.mrna[1:3,1:3]

enerly.clinical<-fread("COLLABORATION/KANG/CL/ENERLY/enerly.er.csv", header=F,na.strings = "")
enerly.clinical<-enerly.clinical[!is.na(V2),]
table(enerly.clinical$V2) #more positive than negative
sam.enerly<-c(sample(enerly.clinical[V2=="estrogen receptor status: Positive",]$V1, 40), enerly.clinical[V2=="estrogen receptor status: Negative",]$V1)
med.genes<-apply(enerly.mrna[,sam.enerly], 1, median)
PRED.CL.ENERLY<-Function.CL.PRED(CL.CENTROIDS, scale(enerly.mrna - med.genes) , scale=F)
PRED.ENERLY.SAMPLES<-PRED.CL.ENERLY[PRED=="CLAUDIN.LOW",]$SAMPLE

enerly.ids<-fread("COLLABORATION/KANG/CL/ENERLY/enerly.ids.csv", header=T)
enerly.ids$PRED<-enerly.ids$mrna.ID %in% PRED.ENERLY.SAMPLES

enerly.mirna<-fread("COLLABORATION/KANG/CL/ENERLY/mirna.clean.txt", header = T, sep="\t")
enerly.mirna<-data.matrix(data.frame(enerly.mirna, row.names = 1))
enerly.mirna<-enerly.mirna[complete.cases(enerly.mirna),]

enerly.mirna<-scale(enerly.mirna- apply(enerly.mirna, 1, median))
enerly.199<-data.table(t(enerly.mirna[c("hsa-miR-199a-5p", "hsa-miR-199b-3p", "hsa-miR-199b-5p"),]), keep.rownames = T)
enerly.199<-melt(enerly.199, id.vars = "rn")
enerly.199<-enerly.199[rn %in% enerly.ids$mirna.ID,]
enerly.199$PRED<-enerly.199$rn %in% enerly.ids[PRED==TRUE,]$mirna.ID

ggplot(enerly.199, aes(variable, value, colour=PRED)) + geom_boxplot() + theme.format
wilcox.test(enerly.199[variable=="hsa-miR-199a-5p" & PRED==TRUE,]$value,
            enerly.199[variable=="hsa-miR-199a-5p" & PRED==FALSE,]$value,
            alternative="greater", paired=F)

###Differential expression between CL and rest of samples
diff.exp.cl<-brca.V2.RAW.TO.MA$E[,grepl("01A", colnames(brca.V2.RAW.TO.MA$E))]
cl.samples<-colnames(diff.exp.cl)[colnames(diff.exp.cl) %in% BRCA.V2.RAW.TCGA.PRED[PRED=="CLAUDIN.LOW" & grepl("01A", SAMPLE),]$SAMPLE]
non.cl.samples<-setdiff(colnames(diff.exp.cl), cl.samples)

cl.design.matrix<-data.frame(TYPE=c(rep("NON.CL", length(non.cl.samples)), rep("CL", length(cl.samples))))
cl.design.matrix<-model.matrix(~TYPE, cl.design.matrix)

diff.exp.cl<-diff.exp.cl[,c(cl.samples, non.cl.samples)]
cl.fit = lmFit(diff.exp.cl, cl.design.matrix)
cl.eb<-eBayes(cl.fit)

cl.diff.res<-as.data.table(topTable(cl.eb,coef=2,number=25000),keep.rownames=T)
cl.diff.res[adj.P.Val<0.05,]
cl.diff.res[rn=="LCOR",]

brca.MIRNA.HISEQ.TO.MA$E[1:3,1:3]
write.table(brca.MIRNA.HISEQ.TO.MA$E, "DATABASES/CANCER_DATA/PEROU/BRCA.mirna.expression", quote=F, sep="\t",row.names=T, col.names=T)
brca.V2.RAW.TO.MA$E[1:3,1:3]
write.table(brca.V2.RAW.TO.MA$E, "DATABASES/CANCER_DATA/PEROU/BRCA.gene.expression", quote=F, sep="\t",row.names=T, col.names=T)
dim(brca.V2.RAW.TO.MA$E)
dim(normalized.mirna.matrix)
length(unique(colnames(normalized.mirna.matrix)))