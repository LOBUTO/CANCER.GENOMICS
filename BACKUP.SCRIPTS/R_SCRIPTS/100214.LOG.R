library(ggplot2)
library(data.table)

setwd("Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS")

BRCA.EXP.OBJ<-readRDS("BRCA/091514.CANCER.MATRICES.NORMALIZED.OBJ.rds")
#table.1.plus<-readRDS("BRCA/082714.BRCA.Table.1.PLUS.rds") - NOT ACCURATE
table.1<-readRDS("BRCA/091714.BRCA.Table.1.rds")
CLINICAL<-readRDS("BRCA/093014.BRCA.CLINICAL.rds")
PAIRED.IDS<-substr(BRCA.EXP.OBJ$normal.patients,1,12)

########101214########
#Do differential expression and cluster patients based on top n genes
#   Aim is to have all normal cluster together away from cancer so we can have a good distance separation

#Use Function.RNAseq.Differential.Expression.V2() from Functions.R()
BRCA.DIFF.EXP<-Function.RNAseq.Differential.Expression.V2(BRCA.EXP.OBJ, BRCA.EXP.OBJ$cancer.patients)

ggplot(BRCA.DIFF.EXP, aes(logFC)) + geom_histogram() + theme.format

BRCA.EXP.OBJ.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)),
        rep("CANCER", length(BRCA.EXP.OBJ$cancer.patients))
    ) )
rownames(BRCA.EXP.OBJ.ANN)<-c(BRCA.EXP.OBJ$normal.patients, BRCA.EXP.OBJ$cancer.patients)

pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>1.5,]$ID),
    c(BRCA.EXP.OBJ$normal.patients, BRCA.EXP.OBJ$cancer.patients)],
    scale="none", annotation=BRCA.EXP.OBJ.ANN,
    clustering_distance_cols="correlation")

#Ordered by distance
pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>1.5,]$ID),
    c(BRCA.EXP.OBJ$normal.patients, as.vector(BRCA.CANCER.DIST.1.5[order(DIST.TO.NORMAL),]$PATIENT))],
    scale="none", annotation=BRCA.EXP.OBJ.ANN, cluster_cols=FALSE)

#####TEST DIFFERENT THRESHOLDS#######

#@1.5

#Get correlation matrix of threshold at abs(logFC)>1.5
BRCA.CORR.1.5<-cor(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>1.5,]$ID),c(BRCA.EXP.OBJ$normal.patients,
    BRCA.EXP.OBJ$cancer.patients)], method="spearman")
pheatmap(BRCA.CORR.1.5, scale="none", annotation=BRCA.EXP.OBJ.ANN, clustering_distance_cols="correlation")

#Get distance measure from normal for each cancer patient
BRCA.DISS.1.5<-1-BRCA.CORR.1.5

BRCA.CANCER.DIST.1.5<-sapply(BRCA.EXP.OBJ$cancer.patients, function(x) median(BRCA.DISS.1.5[x, BRCA.EXP.OBJ$normal.patients]))
BRCA.CANCER.DIST.1.5<-as.data.table(BRCA.CANCER.DIST.1.5, keep.rownames=T)
setnames(BRCA.CANCER.DIST.1.5, c("PATIENT", "DIST.TO.NORMAL"))

ggplot(BRCA.CANCER.DIST.1.5, aes(DIST.TO.NORMAL)) + geom_histogram() + theme.format

#Test if correlation to number of mutations still exist...
table.1.pval.count
BRCA.CANCER.DIST.PVAL.1.5<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST.1.5), as.data.frame(table.1.pval.count), by="PATIENT"))
BRCA.CANCER.DIST.PVAL.1.5$PAIRED<-substr(BRCA.CANCER.DIST.PVAL.1.5$PATIENT,1,12) %in% PAIRED.IDS

ggplot(BRCA.CANCER.DIST.PVAL.1.5, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + facet_wrap(~PAIRED) + scale_y_log10()
cor.test(BRCA.CANCER.DIST.PVAL.1.5$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.1.5$N.MUTATIONS, method="spearman")
cor.test(BRCA.CANCER.DIST.PVAL.1.5[PAIRED==TRUE,]$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.1.5[PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Test clinical with new distance
test.1.5<-copy(BRCA.CANCER.DIST.PVAL.1.5)
test.1.5<-as.data.table(merge(as.data.frame(test.1.5), as.data.frame(CLINICAL), by="PATIENT"))

#SURVIVAL/PROCUREMENT
ggplot(test.1.5[GENDER=="FEMALE",], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL)) + geom_point() + scale_y_log10()
cor.test(test.1.5[GENDER=="FEMALE",]$DIST.TO.NORMAL, test.1.5[GENDER=="FEMALE",]$OVERALL.SURVIVAL, method="spearman")

ggplot(test.1.5[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=VITAL.STATUS)) + geom_point() + scale_y_log10() + facet_wrap(~VITAL.STATUS)
cor.test(test.1.5[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$DIST.TO.NORMAL,
    test.1.5[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$OVERALL.SURVIVAL,
    method="spearman")
cor.test(as.vector(test.1.5[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$DIST.TO.NORMAL),
    as.vector(test.1.5[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$OVERALL.SURVIVAL),
    method="spearman")

ggplot(test.1.5[DAYS.TO.PROCUREMENT>=0,], aes(DIST.TO.NORMAL, DAYS.TO.PROCUREMENT)) + geom_point() + scale_y_log10()

#@2.0

#Get correlation matrix of threshold at abs(logFC)>1.5
BRCA.CORR.2.0<-cor(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>2.0,]$ID),c(BRCA.EXP.OBJ$normal.patients,
    BRCA.EXP.OBJ$cancer.patients)], method="spearman")
pheatmap(BRCA.CORR.2.0, scale="none", annotation=BRCA.EXP.OBJ.ANN, clustering_distance_cols="correlation")

#Get distance measure from normal for each cancer patient
BRCA.DISS.2.0<-1-BRCA.CORR.2.0

BRCA.CANCER.DIST.2.0<-sapply(BRCA.EXP.OBJ$cancer.patients, function(x) median(BRCA.DISS.2.0[x, BRCA.EXP.OBJ$normal.patients]))
BRCA.CANCER.DIST.2.0<-as.data.table(BRCA.CANCER.DIST.2.0, keep.rownames=T)
setnames(BRCA.CANCER.DIST.2.0, c("PATIENT", "DIST.TO.NORMAL"))

ggplot(BRCA.CANCER.DIST.2.0, aes(DIST.TO.NORMAL)) + geom_histogram() + theme.format

#Test if correlation to number of mutations still exist...
table.1.pval.count
BRCA.CANCER.DIST.PVAL.2.0<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST.2.0), as.data.frame(table.1.pval.count), by="PATIENT"))
BRCA.CANCER.DIST.PVAL.2.0$PAIRED<-substr(BRCA.CANCER.DIST.PVAL.2.0$PATIENT,1,12) %in% PAIRED.IDS

ggplot(BRCA.CANCER.DIST.PVAL.2.0, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + facet_wrap(~PAIRED) + scale_y_log10()
cor.test(BRCA.CANCER.DIST.PVAL.2.0$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.2.0$N.MUTATIONS, method="spearman")
cor.test(BRCA.CANCER.DIST.PVAL.2.0[PAIRED==TRUE,]$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.2.0[PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Test clinical with new distance
test.2.0<-copy(BRCA.CANCER.DIST.PVAL.2.0)
test.2.0<-as.data.table(merge(as.data.frame(test.2.0), as.data.frame(CLINICAL), by="PATIENT"))

#SURVIVAL/PROCUREMENT
ggplot(test.2.0[GENDER=="FEMALE",], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL)) + geom_point() + scale_y_log10()
cor.test(test.2.0[GENDER=="FEMALE",]$DIST.TO.NORMAL, test.2.0[GENDER=="FEMALE",]$OVERALL.SURVIVAL, method="spearman")

ggplot(test.2.0[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=VITAL.STATUS)) + geom_point() + scale_y_log10() + facet_wrap(~VITAL.STATUS)
cor.test(test.2.0[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$DIST.TO.NORMAL,
    test.2.0[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$OVERALL.SURVIVAL,
    method="spearman")
cor.test(as.vector(test.2.0[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$DIST.TO.NORMAL),
    as.vector(test.2.0[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$OVERALL.SURVIVAL),
    method="spearman")

ggplot(test.2.0[DAYS.TO.PROCUREMENT>=0,], aes(DIST.TO.NORMAL, DAYS.TO.PROCUREMENT)) + geom_point() + scale_y_log10()

#@3.0

#Get correlation matrix of threshold at abs(logFC)>1.5
BRCA.CORR.3.0<-cor(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>3.0,]$ID),c(BRCA.EXP.OBJ$normal.patients,
    BRCA.EXP.OBJ$cancer.patients)], method="spearman")
pheatmap(BRCA.CORR.3.0, scale="none", annotation=BRCA.EXP.OBJ.ANN, clustering_distance_cols="correlation")

#Get distance measure from normal for each cancer patient
BRCA.DISS.3.0<-1-BRCA.CORR.3.0

BRCA.CANCER.DIST.3.0<-sapply(BRCA.EXP.OBJ$cancer.patients, function(x) median(BRCA.DISS.3.0[x, BRCA.EXP.OBJ$normal.patients]))
BRCA.CANCER.DIST.3.0<-as.data.table(BRCA.CANCER.DIST.3.0, keep.rownames=T)
setnames(BRCA.CANCER.DIST.3.0, c("PATIENT", "DIST.TO.NORMAL"))

ggplot(BRCA.CANCER.DIST.3.0, aes(DIST.TO.NORMAL)) + geom_histogram() + theme.format

#Test if correlation to number of mutations still exist...
table.1.pval.count
BRCA.CANCER.DIST.PVAL.3.0<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST.3.0), as.data.frame(table.1.pval.count), by="PATIENT"))
BRCA.CANCER.DIST.PVAL.3.0$PAIRED<-substr(BRCA.CANCER.DIST.PVAL.3.0$PATIENT,1,12) %in% PAIRED.IDS

ggplot(BRCA.CANCER.DIST.PVAL.3.0, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + facet_wrap(~PAIRED) + scale_y_log10()
cor.test(BRCA.CANCER.DIST.PVAL.3.0$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.3.0$N.MUTATIONS, method="spearman")
cor.test(BRCA.CANCER.DIST.PVAL.3.0[PAIRED==TRUE,]$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL.3.0[PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Test clinical with new distance
test.3.0<-copy(BRCA.CANCER.DIST.PVAL.3.0)
test.3.0<-as.data.table(merge(as.data.frame(test.3.0), as.data.frame(CLINICAL), by="PATIENT"))

#SURVIVAL/PROCUREMENT
ggplot(test.3.0[GENDER=="FEMALE",], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL)) + geom_point() + scale_y_log10()
cor.test(test.3.0[GENDER=="FEMALE",]$DIST.TO.NORMAL, test.3.0[GENDER=="FEMALE",]$OVERALL.SURVIVAL, method="spearman")

ggplot(test.3.0[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,][ER.STATUS %in% c("Positive", "Negative"),], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=ER.STATUS)) + geom_point() + scale_y_log10() + facet_wrap(~VITAL.STATUS) + stat_smooth(method="lm", formula = y ~ log(x))

ggplot(test.3.0[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL/30, colour=VITAL.STATUS)) + geom_point() + facet_wrap(~VITAL.STATUS)
cor.test(test.3.0[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$DIST.TO.NORMAL,
    test.3.0[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS=="Alive",]$OVERALL.SURVIVAL,
    method="spearman")
cor.test(as.vector(test.3.0[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$DIST.TO.NORMAL),
    as.vector(test.3.0[OVERALL.SURVIVAL>0,][VITAL.STATUS=="Dead",]$OVERALL.SURVIVAL),
    method="spearman")

ggplot(test.3.0[DAYS.TO.PROCUREMENT>=0,], aes(DIST.TO.NORMAL, DAYS.TO.PROCUREMENT)) + geom_point() + scale_y_log10()

####Influence of receptor on Expression Patterns####
CLINICAL.TYPE<-CLINICAL[, c(1,2,10:12,15), with=F]

#@ER
CLINICAL.ER<-CLINICAL.TYPE[ER.STATUS %in% c("Positive", "Negative") & GENDER=="FEMALE",]
ER.POSITIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.ER[ER.STATUS=="Positive",]$PATIENT)
ER.NEGATIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.ER[ER.STATUS=="Negative",]$PATIENT)

ER.ANN<-data.frame(ER.TYPE=c(rep("Normal", length(BRCA.EXP.OBJ$normal.patients)),
    rep("Negative", length(ER.NEGATIVE.PATIENTS)),
    rep("Positive", length(ER.POSITIVE.PATIENTS))))
rownames(ER.ANN)<-c(BRCA.EXP.OBJ$normal.patients,ER.NEGATIVE.PATIENTS, ER.POSITIVE.PATIENTS)

pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>3.0,]$ID),
    c(BRCA.EXP.OBJ$normal.patients,ER.NEGATIVE.PATIENTS, ER.POSITIVE.PATIENTS)],
    scale="none", annotation=ER.ANN,
    clustering_distance_cols="euclidean")

#@PR
CLINICAL.PR<-CLINICAL.TYPE[PR.STATUS %in% c("Positive", "Negative"),]
PR.POSITIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.PR[PR.STATUS=="Positive",]$PATIENT)
PR.NEGATIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.PR[PR.STATUS=="Negative",]$PATIENT)

PR.ANN<-data.frame(PR.TYPE=c(rep("Normal", length(BRCA.EXP.OBJ$normal.patients)),
    rep("Negative", length(PR.NEGATIVE.PATIENTS)),
    rep("Positive", length(PR.POSITIVE.PATIENTS))))
rownames(PR.ANN)<-c(BRCA.EXP.OBJ$normal.patients,PR.NEGATIVE.PATIENTS, PR.POSITIVE.PATIENTS)

pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>3.0,]$ID),
    c(BRCA.EXP.OBJ$normal.patients,PR.NEGATIVE.PATIENTS, PR.POSITIVE.PATIENTS)],
    scale="none", annotation=PR.ANN,
    clustering_distance_cols="euclidean")

#@HER2
CLINICAL.HER2<-CLINICAL.TYPE[HER2.STATUS %in% c("Positive", "Negative"),]
HER2.POSITIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.HER2[HER2.STATUS=="Positive",]$PATIENT)
HER2.NEGATIVE.PATIENTS<-intersect(BRCA.EXP.OBJ$cancer.patients, 
    CLINICAL.HER2[HER2.STATUS=="Negative",]$PATIENT)

HER2.ANN<-data.frame(HER2.TYPE=c(rep("Normal", length(BRCA.EXP.OBJ$normal.patients)),
    rep("Negative", length(HER2.NEGATIVE.PATIENTS)),
    rep("Positive", length(HER2.POSITIVE.PATIENTS))))
rownames(HER2.ANN)<-c(BRCA.EXP.OBJ$normal.patients,HER2.NEGATIVE.PATIENTS, HER2.POSITIVE.PATIENTS)

pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>3.0,]$ID),
    c(BRCA.EXP.OBJ$normal.patients,HER2.NEGATIVE.PATIENTS, HER2.POSITIVE.PATIENTS)],
    scale="none", annotation=HER2.ANN,
    clustering_distance_cols="correlation")