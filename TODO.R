cl.subtypes<-fread("~/Dropbox/BC datasets collaboration/results/052015.RSEM.INTRINSIC.SUBTYPES.txt", header=T, sep="\t")
table(cl.subtypes$SUBTYPE)

wilcox.test(as.numeric(brca.exp$tumor["LCOR",cl.subtypes[SUBTYPE=="CL"]$SAMPLE]), as.numeric(brca.exp$normal["LCOR",]))
log2(median(as.numeric(brca.exp$tumor["LCOR",cl.subtypes[SUBTYPE=="CL"]$SAMPLE]))/ median(as.numeric(brca.exp$normal["LCOR",])))

cl.lcor<-data.table(LCOR.EXP=c(as.numeric(brca.exp$tumor["LCOR",cl.subtypes[SUBTYPE=="CL"]$SAMPLE]), as.numeric(brca.exp$normal["LCOR",]) ),
           TYPE=c(rep("CL", length(as.numeric(brca.exp$tumor["LCOR",cl.subtypes[SUBTYPE=="CL"]$SAMPLE]))),
                  rep("Normal Tissue", length(as.numeric(brca.exp$normal["LCOR",])))))
write.table(cl.lcor, file = "~/Dropbox/BC datasets collaboration/results/081915.CL.VS.NORMAL.LCOR.txt", quote = F,sep = "\t", row.names = F, col.names = F)

emf("COLLABORATION/KANG/FIGURES/CL.VS.NORMAL.LCOR.emf")
ggplot(cl.lcor, aes(TYPE, LCOR.EXP)) + geom_boxplot()
dev.off()

x<-merge(cl.subtypes[SUBTYPE=="CL",], TCGA.BRCA.CLINICAL[,c("ER.STATUS", "SAMPLE"), with=F], by="SAMPLE")
table(x$ER.STATUS)
lcor.cancer<-data.table(t(brca.exp$tumor["LCOR", ]), keep.rownames = T)
setnames(lcor.cancer, c("SAMPLE", "EXP"))
lcor.normal<-data.table(t(brca.exp$normal["LCOR", ]), keep.rownames = T)
setnames(lcor.normal, c("SAMPLE", "EXP"))
lcor.normal$SUBTYPE<-"NORMAL"
lcor.normal$ER.STATUS<-"NORMAL"

x<-merge(x, lcor.cancer, by="SAMPLE")
y<-rbind(x, lcor.normal)
y$ER.STATUS<-as.factor(y$ER.STATUS)
write.table(y, "~/Dropbox/BC datasets collaboration/results/LCOR.txt", quote = F, sep = "\t", row.names = F)

emf("COLLABORATION/KANG/FIGURES/CL.ER.VS.NORMAL.LCOR.EXP.emf")
ggplot(y, aes(SUBTYPE, EXP, colour=ER.STATUS)) + geom_boxplot()  + ylab("LCOR.EXP")
dev.off()

table(y$ER.STATUS)

table(TCGA.BRCA.CLINICAL$ER.STATUS)/2

#Notes:
#1. En TCGA CL vs Normal (non-tumor) en 199a (er- vs er+)
#2. En Enerly CL vs Normal (non-tumor) en 199a y LCOR (todos y er- vs er+)
#3. En TCGA, CL vs Normal (non-tumor) 199a y 199b boxplots
#4. Save 052015.BRCA.CL.9.RSEM.199.ALL.SUBTYPES.tiff as enhanced metafile (and others)
#5. Farazi all types (inc. medullary and metaplastic) and compare 199a/b - GEO GSE28884

####### 1 #######
cl.subtypes
TCGA.BRCA.CLINICAL<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/081415.BRCA.CLINICAL.rds")

normalized.mirna.matrix

normal.samples<-colnames(normalized.mirna.matrix)[grepl("11A",colnames(normalized.mirna.matrix))]

table.1<-melt(normalized.mirna.matrix[c("hsa-mir-199a-1","hsa-mir-199a-2"), union(cl.subtypes[SUBTYPE=="CL",]$SAMPLE, normal.samples)])
setnames(table.1, c("miRNA", "SAMPLE", "Expression"))
table.1$TYPE<-ifelse(table.1$SAMPLE %in% normal.samples, "NORMAL", "CL")
table.1
table.1$MOL.SUBTYPE<-ifelse(table.1$SAMPLE %in% TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, "ER.POS", 
                            ifelse(table.1$SAMPLE %in% TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, "ER.NEG", "NORMAL"))
table.1<-as.data.table(table.1)
write.table(table.1, "~/Dropbox/BC datasets collaboration/results/table.1.txt", quote = F, sep = "\t", row.names = F)

emf("COLLABORATION/KANG/FIGURES/CL.VS.NORMAL.199A.emf" )
ggplot(table.1, aes(miRNA, Expression, colour=TYPE)) + geom_boxplot()
dev.off()

emf("COLLABORATION/KANG/FIGURES/CL.VS.NORMAL.199A.ER.emf" )
ggplot(table.1[!(TYPE=="CL" & MOL.SUBTYPE=="NORMAL"),], aes(miRNA, Expression, colour=MOL.SUBTYPE)) + geom_boxplot()
dev.off()

######## 2 ########
enerly

####### 3 #######
table.3<-melt(normalized.mirna.matrix[c("hsa-mir-199a-1","hsa-mir-199a-2", "hsa-mir-199b"), union(cl.subtypes[SUBTYPE=="CL",]$SAMPLE,normal.samples)])
setnames(table.3, c("miRNA", "SAMPLE", "Expression"))
table.3<-as.data.table(table.3)
table.3$TYPE<-ifelse(table.3$SAMPLE %in% normal.samples, "NORMAL", "CL")
write.table(table.3, "~/Dropbox/BC datasets collaboration/results/table.3.txt", quote = F, sep = "\t", row.names = F)

emf("COLLABORATION/KANG/FIGURES/CL.VS.NORMAL.199AB.emf")
ggplot(table.3, aes(miRNA, Expression, colour=TYPE)) + geom_boxplot()
dev.off()

####### 4 #######
table.4<-as.data.table(melt(normalized.mirna.matrix[c("hsa-mir-199a-1","hsa-mir-199a-2"),]))
setnames(table.4, c("miRNA", "SAMPLE", "Expression"))
table.4<-merge(table.4 , cl.subtypes, by="SAMPLE")
write.table(table.4, "~/Dropbox/BC datasets collaboration/results/table.4.txt", quote = F, sep = "\t", row.names = F)

emf("COLLABORATION/KANG/FIGURES/052015.BRCA.CL.9.RSEM.199.ALL.SUBTYPES.emf")
ggplot(table.4, aes(miRNA, Expression, colour=SUBTYPE)) + geom_boxplot()
dev.off()

####### 5 ########
library(Biobase)
library(GEOquery)
gse29174<-getGEO("GSE28884", GSEMatrix = T)

gse.test<-getGEO(filename = system.file("DATABASES/CANCER_DATA/FARAZI/GSE29173_family.soft.gz", package = "GEOquery"))
