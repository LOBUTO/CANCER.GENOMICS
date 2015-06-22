library(ggplot2)
library(data.table)

setwd("Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS")

BRCA.EXP.OBJ<-readRDS("BRCA/091514.CANCER.MATRICES.NORMALIZED.OBJ.rds")
table.1.plus<-readRDS("BRCA/082714.BRCA.Table.1.PLUS.rds")
table.1<-readRDS("BRCA/091714.BRCA.Table.1.rds")
CLINICAL<-readRDS("BRCA/093014.BRCA.CLINICAL.rds")
PAIRED.IDS<-substr(BRCA.EXP.OBJ$normal.patients,1,12)


#Get cancer patient distances in term of expression to normal
#This distance will be calculated based on the median distance of a cancer patient to all normal

BRCA.CORR<-cor(BRCA.EXP.OBJ$combined.matrices[,c(BRCA.EXP.OBJ$normal.patients,
    BRCA.EXP.OBJ$cancer.patients)], method="spearman")
BRCA.DISS<-1-BRCA.CORR

BRCA.CANCER.DIST<-sapply(BRCA.EXP.OBJ$cancer.patients, function(x) median(BRCA.DISS[x, BRCA.EXP.OBJ$normal.patients]))
BRCA.CANCER.DIST<-as.data.table(BRCA.CANCER.DIST, keep.rownames=T)
setnames(BRCA.CANCER.DIST, c("PATIENT", "DIST.TO.NORMAL"))

#Including CNV
table.1.plus.count<-table.1.plus[,list(N.MUTATIONS=length(Hugo_Symbol)), by="PATIENT"]
ggplot(table.1.plus.count, aes(N.MUTATIONS))+geom_histogram()

BRCA.CANCER.DIST.PLUS<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST), 
    as.data.frame(table.1.plus.count)))
BRCA.CANCER.DIST.PLUS$PAIRED<-substr(BRCA.CANCER.DIST.PLUS$PATIENT,1,12) %in% PAIRED.IDS

ggplot(BRCA.CANCER.DIST.PLUS[PATIENT %in% BRCA.CANCER.DIST.NOCNV$PATIENT,], aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + scale_y_log10()
cor.test(BRCA.CANCER.DIST.PLUS[PATIENT %in% BRCA.CANCER.DIST.NOCNV$PATIENT,][PAIRED==TRUE,]$DIST.TO.NORMAL, BRCA.CANCER.DIST.PLUS[PATIENT %in% BRCA.CANCER.DIST.NOCNV$PATIENT,][PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Without CNV
table.1.count<-table.1$table.1[,list(N.MUTATIONS=sum(Missense)), by="PATIENT"]
ggplot(table.1.count, aes(N.MUTATIONS)) + geom_histogram()

BRCA.CANCER.DIST.NOCNV<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST), 
    as.data.frame(table.1.count), by="PATIENT"))
BRCA.CANCER.DIST.NOCNV$PAIRED<-substr(BRCA.CANCER.DIST.NOCNV$PATIENT,1,12) %in% PAIRED.IDS

ggplot(BRCA.CANCER.DIST.NOCNV, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + scale_y_log10()
cor.test(BRCA.CANCER.DIST.NOCNV$DIST.TO.NORMAL, BRCA.CANCER.DISTw.NOCNV$N.MUTATIONS, method="spearman")
cor.test(BRCA.CANCER.DIST.NOCNV[PAIRED==TRUE,]$DIST.TO.NORMAL, 
    BRCA.CANCER.DIST.NOCNV[PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Get significance level per gene without CNV - Done with Function.table.1.bmr.R
ggplot(table.1$table.1, aes(Hugo_Symbol, Missense)) + geom_point()

table.1.pval<-readRDS("BRCA/100414.BRCA.Table.1.bmr.rds")
table.1.pval<-table.1.pval[P.VAL.ADJ<0.05,] #Still covering all patients
PATIENT.WITH.COSMIC<-unique(as.vector(table.1.pval[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]$PATIENT))

PATIENTS.WITH.TP53<-unique(as.vector(table.1.pval[Hugo_Symbol=="TP53",]$PATIENT))
PATIENTS.WITH.TTN<-unique(as.vector(table.1.pval[Hugo_Symbol=="TTN",]$PATIENT))
PATIENTS.WITH.KLK15<-unique(as.vector(table.1.pval[Hugo_Symbol=="KLK15",]$PATIENT))

table.1.pval.count<-table.1.pval[,list(N.MUTATIONS=sum(Missense)), by="PATIENT"]

BRCA.CANCER.DIST.PVAL<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST), 
    as.data.frame(table.1.pval.count), by="PATIENT"))

BRCA.CANCER.DIST.PVAL$PAIRED<-substr(BRCA.CANCER.DIST.PVAL$PATIENT,1,12) %in% PAIRED.IDS
BRCA.CANCER.DIST.PVAL$COSMIC<-BRCA.CANCER.DIST.PVAL$PATIENT %in% PATIENT.WITH.COSMIC
BRCA.CANCER.DIST.PVAL$VITAL.STATUS<-ifelse(as.vector(substr(BRCA.CANCER.DIST.PVAL$PATIENT, 1,12)) %in% as.vector(CLINICAL[VITAL.STATUS=="Alive",]$PATIENT), "Alive",
    ifelse(as.vector(substr(BRCA.CANCER.DIST.PVAL$PATIENT, 1,12)) %in% as.vector(CLINICAL[VITAL.STATUS=="Dead",]$PATIENT), "Dead", "NO.DATA"))

BRCA.CANCER.DIST.PVAL$TP53<-BRCA.CANCER.DIST.PVAL$PATIENT %in% PATIENTS.WITH.TP53
BRCA.CANCER.DIST.PVAL$TTN<-BRCA.CANCER.DIST.PVAL$PATIENT %in% PATIENTS.WITH.TTN
BRCA.CANCER.DIST.PVAL$KLK15<-BRCA.CANCER.DIST.PVAL$PATIENT %in% PATIENTS.WITH.KLK15

ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=PAIRED)) + geom_point() + scale_y_log10()
ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=COSMIC)) + geom_point() + scale_y_log10()
ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=VITAL.STATUS)) + geom_point() + scale_y_log10() + facet_wrap(~VITAL.STATUS)

ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=TP53)) + geom_point() + scale_y_log10()
ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=TTN)) + geom_point() + scale_y_log10()
ggplot(BRCA.CANCER.DIST.PVAL, aes(DIST.TO.NORMAL, N.MUTATIONS, colour=KLK15)) + geom_point() + scale_y_log10()

cor.test(BRCA.CANCER.DIST.PVAL$DIST.TO.NORMAL, BRCA.CANCER.DIST.PVAL$N.MUTATIONS, method="spearman")
cor.test(BRCA.CANCER.DIST.PVAL[PAIRED==TRUE,]$DIST.TO.NORMAL, 
    BRCA.CANCER.DIST.PVAL[PAIRED==TRUE,]$N.MUTATIONS, method="spearman")

#Same significant, however, compare most frequent genes to before - TTN IS NOT HIGH ANYMORE
table.1.NOCNV.FREQ<-as.data.table(table(table.1$table.1$Hugo_Symbol), keep.rownames=T)
setnames(table.1.NOCNV.FREQ, c("Hugo_Symbol", "Freq"))
table.1.NOCNV.FREQ<-table.1.NOCNV.FREQ[order(Freq, decreasing=T),]
head(table.1.NOCNV.FREQ,10)

table.1.PVAL.FREQ<-as.data.table(table(table.1.pval$Hugo_Symbol), keep.rownames=T)
setnames(table.1.PVAL.FREQ, c("Hugo_Symbol", "Freq"))
table.1.PVAL.FREQ<-table.1.PVAL.FREQ[order(Freq, decreasing=T),]
head(table.1.PVAL.FREQ,10)

COSMIC<-as.data.table(read.csv("~/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/COSMIC/cancer_gene_census.csv", header=T,sep=","))
COSMIC.BRCA<-COSMIC[apply(COSMIC, 1, function(x) any(grepl("breast",x,ignore.case=T))),]
saveRDS(as.vector(COSMIC.BRCA$Symbol), file="BRCA/091814.COSMIC.BRCA.rds")
COSMIC.GBM<-COSMIC[apply(COSMIC, 1, function(x) any(grepl("glioblastoma",x,ignore.case=T))),]
COSMIC.OV<-COSMIC[apply(COSMIC, 1, function(x) any(grepl("ovar",x,ignore.case=T))),]
COSMIC.READ<-COSMIC[apply(COSMIC, 1, function(x) any(grepl("rectal",x,ignore.case=T))),]

table.1.PVAL.FREQ[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

length(unique(as.vector(table.1$table.1$Hugo_Symbol))) #20143
length(unique(as.vector(table.1.pval$Hugo_Symbol))) #16001

#Develop algorithm to choose significant genes from table.1.pval starting to those patients that are closer to normal - DEVELOPED
BRCA.CANCER.DIST
table.1.pval

#####ANALYZE####
CNV<-readRDS("BRCA/092214.BRCA.GISTIC.TH.2.rds")
TH.1<-readRDS("BRCA/092314.BRCA.GISTIC.FILTERED.TH.1.rds")
TH.2<-readRDS("BRCA/092314.BRCA.GISTIC.FILTERED.TH.2.rds")

CNV.COUNT<-CNV[,list(N.CNV=sum(abs(CNV.TH))), by="PATIENT"]
CNV.COUNT<-as.data.table(merge(as.data.frame(CNV.COUNT), as.data.frame(BRCA.CANCER.DIST), by="PATIENT"))
CNV.COUNT$PAIRED<-substr(CNV.COUNT$PATIENT,1,12) %in% PAIRED.IDS
ggplot(CNV.COUNT, aes(DIST.TO.NORMAL, N.CNV, colour=PAIRED)) + geom_point() 
cor.test(CNV.COUNT$DIST.TO.NORMAL, CNV.COUNT$N.CNV, method="spearman")
cor.test(CNV.COUNT[PAIRED==TRUE,]$DIST.TO.NORMAL, CNV.COUNT[PAIRED==TRUE,]$N.CNV, method="spearman")

TH.1.COUNT<-TH.1[,list(N.CNV=sum(abs(CNV.TH))), by="PATIENT"]
TH.1.COUNT<-as.data.table(merge(as.data.frame(TH.1.COUNT), as.data.frame(BRCA.CANCER.DIST), by="PATIENT"))
TH.1.COUNT$PAIRED<-substr(TH.1.COUNT$PATIENT,1,12) %in% PAIRED.IDS
ggplot(TH.1.COUNT, aes(DIST.TO.NORMAL, N.CNV, colour=PAIRED)) + geom_point()
cor.test(TH.1.COUNT$DIST.TO.NORMAL, TH.1.COUNT$N.CNV, method="spearman")
cor.test(TH.1.COUNT[PAIRED==TRUE,]$DIST.TO.NORMAL, TH.1.COUNT[PAIRED==TRUE,]$N.CNV, method="spearman")

TH.2.COUNT<-TH.2[,list(N.CNV=sum(abs(CNV.TH))), by="PATIENT"]
TH.2.COUNT<-as.data.table(merge(as.data.frame(TH.2.COUNT), as.data.frame(BRCA.CANCER.DIST), by="PATIENT"))
TH.2.COUNT$PAIRED<-substr(TH.2.COUNT$PATIENT,1,12) %in% PAIRED.IDS
ggplot(TH.2.COUNT, aes(DIST.TO.NORMAL, N.CNV,colour=PAIRED)) + geom_point() 
cor.test(TH.2.COUNT$DIST.TO.NORMAL, TH.2.COUNT$N.CNV, method="spearman")
cor.test(TH.2.COUNT[PAIRED==TRUE,]$DIST.TO.NORMAL, TH.2.COUNT[PAIRED==TRUE,]$N.CNV, method="spearman")

CNV.WITH.DCAF13<-unique(as.vector(TH.1[Hugo_Symbol=="DCAF13",]$PATIENT))
CNV.WITH.ARF1<-unique(as.vector(TH.1[Hugo_Symbol=="ARF1",]$PATIENT))
CNV.WITH.KLK13<-unique(as.vector(TH.1[Hugo_Symbol=="KLK13",]$PATIENT))
TH.1.COUNT$DCAF13<-TH.1.COUNT$PATIENT %in% CNV.WITH.DCAF13
TH.1.COUNT$ARF1<-TH.1.COUNT$PATIENT %in% CNV.WITH.ARF1
TH.1.COUNT$KLK13<-TH.1.COUNT$PATIENT %in% CNV.WITH.KLK13
ggplot(TH.1.COUNT, aes(DIST.TO.NORMAL, N.CNV,colour=KLK13)) + geom_point() 

#For GM
MUT.METHOD.1<-as.data.table(data.frame(SEED=c(0,10,15,20,22:30,35,40), N.GENES=c(124,115,112,112,111,109,109,109,108,108,108,109,110,112,112)))
ggplot(MUT.METHOD.1, aes(SEED, N.GENES)) + geom_point() + geom_line()

MUT.METHOD.2<-as.data.table(data.frame(SEED=c(0,5,10,15,20,25,30,35), N.GENES=c(116,112,105,102,100,97,99,101)))
ggplot(MUT.METHOD.2, aes(SEED, N.GENES)) + geom_point() + geom_line()

main.table[Hugo_Symbol =="TP53",][TYPE=="MUTATION",]
main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]
main.table[Hugo_Symbol =="PIK3CA",][TYPE=="MUTATION",]

TP53.ARF1<-intersect(main.table[Hugo_Symbol =="TP53",]$PATIENT, main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT )
TP53.NO.ARF1<-setdiff(main.table[Hugo_Symbol =="TP53",]$PATIENT, main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT )
ARF1.NO.TP53<-setdiff(main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT, main.table[Hugo_Symbol =="TP53",]$PATIENT )
ARF1.TTN<-intersect(main.table[Hugo_Symbol =="TTN",]$PATIENT, main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT )
ARF1.NO.TTN<-setdiff(main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT, main.table[Hugo_Symbol =="TTN",]$PATIENT )
ARF1.PIK3CA<-intersect(main.table[Hugo_Symbol =="PIK3CA",]$PATIENT, main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT )
ARF1.NO.PIK3CA<-setdiff(main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT, main.table[Hugo_Symbol =="PIK3CA",]$PATIENT )
ARF1.CCND1<-intersect(main.table[Hugo_Symbol =="CCND1",]$PATIENT, main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT )
ARF1.NO.CCND1<-setdiff(main.table[Hugo_Symbol =="ARF1",][TYPE=="CNV",]$PATIENT, main.table[Hugo_Symbol =="CCND1",]$PATIENT )
TTN.TP53<-intersect(main.table[Hugo_Symbol =="TP53",]$PATIENT, main.table[Hugo_Symbol =="TTN",]$PATIENT )
TTN.NO.TP53<-setdiff(main.table[Hugo_Symbol =="TTN",]$PATIENT, main.table[Hugo_Symbol =="TP53",]$PATIENT )

TEST<-main.table[TYPE=="MUTATION",][,list(N.GENES=length(Hugo_Symbol)), by=c("PATIENT", "DIST.TO.NORMAL")]
TEST$TYPE<-ifelse(TEST$PATIENT %in% TP53.ARF1, "TP53 AND ARF1", 
    ifelse(TEST$PATIENT %in% ARF1.NO.TP53, "ARF1 ONLY", "NONE"))
ggplot(TEST[TYPE!="NONE",], aes(N.GENES,DIST.TO.NORMAL, colour=TYPE)) + geom_point() + scale_x_log10()+  geom_smooth(method = lm, se=FALSE)

TEST$TYPE<-ifelse(TEST$PATIENT %in% ARF1.TTN, "ARF1 AND TTN", 
    ifelse(TEST$PATIENT %in% ARF1.NO.TTN, "ARF1 ONLY", "NONE"))
ggplot(TEST[TYPE!="NONE",], aes(N.GENES,DIST.TO.NORMAL, colour=TYPE)) + geom_point() + scale_x_log10()+  geom_smooth(method = lm, se=FALSE)

TEST$TYPE<-ifelse(TEST$PATIENT %in% ARF1.PIK3CA, "ARF1 AND PIK3CA", 
    ifelse(TEST$PATIENT %in% ARF1.NO.PIK3CA, "ARF1 ONLY", "NONE"))
ggplot(TEST[TYPE!="NONE",], aes(N.GENES,DIST.TO.NORMAL, colour=TYPE)) + geom_point() + scale_x_log10()+  geom_smooth(method = lm, se=FALSE)

TEST$TYPE<-ifelse(TEST$PATIENT %in% ARF1.CCND1, "ARF1 AND CCND1", 
    ifelse(TEST$PATIENT %in% ARF1.NO.CCND1, "ARF1 ONLY", "NONE"))
ggplot(TEST[TYPE!="NONE",], aes(N.GENES,DIST.TO.NORMAL, colour=TYPE)) + geom_point() + scale_x_log10()+  geom_smooth(method = lm, se=FALSE)

TEST$TYPE<-ifelse(TEST$PATIENT %in% TTN.TP53, "TTN AND TP53", 
    ifelse(TEST$PATIENT %in% TTN.NO.TP53, "TTN ONLY", "NONE"))
ggplot(TEST[TYPE!="NONE",], aes(N.GENES,DIST.TO.NORMAL, colour=TYPE)) + geom_point() + scale_x_log10()+  geom_smooth(method = lm, se=FALSE) 

#kernel density based global two-sample comparison test
TEST[TYPE=="ARF1 ONLY",]
TEST[TYPE=="TP53 AND ARF1",]
KDE<-kde.test(as.matrix(TEST[TYPE=="ARF1 ONLY",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]),
    as.matrix(TEST[TYPE=="TP53 AND ARF1",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]))
KDE
kde.test(as.matrix(TEST[TYPE=="TP53 AND ARF1",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]),
    as.matrix(TEST[TYPE=="ARF1 ONLY",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]))

kde.test(as.matrix(TEST[TYPE=="ARF1 AND CCND1",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]),
    as.matrix(TEST[TYPE=="ARF1 ONLY",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]))

kde.test(as.matrix(TEST[TYPE=="TTN AND TP53",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]),
    as.matrix(TEST[TYPE=="TTN ONLY",][,c("DIST.TO.NORMAL", "N.GENES"),with=F]))

#Do KL Divergence Test on Hidden distribution
TEST[TYPE!="NONE",][TYPE=="TTN ONLY",]
H.TTN.NO.TP53<-lapply(1:nrow(TEST[TYPE!="NONE",][TYPE=="TTN ONLY",]) , function(x) 
    rep(TEST[TYPE!="NONE",][TYPE=="TTN ONLY",]$DIST.TO.NORMAL, TEST[TYPE!="NONE",][TYPE=="TTN ONLY",]$N.GENES) )
H.TTN.NO.TP53<-do.call(cbind, H.TTN.NO.TP53)
hist(H.TTN.NO.TP53)
H.TTN.TP53<-lapply(1:nrow(TEST[TYPE!="NONE",][TYPE=="TTN AND TP53",]) , function(x) 
    rep(TEST[TYPE!="NONE",][TYPE=="TTN AND TP53",]$DIST.TO.NORMAL, TEST[TYPE!="NONE",][TYPE=="TTN AND TP53",]$N.GENES) )
H.TTN.TP53<-do.call(cbind, H.TTN.TP53)
hist(H.TTN.TP53)

library(FNN)
KL.divergence(H.TTN.NO.TP53, H.TTN.TP53, k=10)

#####TESTING CLINICAL#######
test<-copy(BRCA.CANCER.DIST.PVAL)
CLINICAL[is.na(CLINICAL$OVERALL.SURVIVAL),]
test<-as.data.table(merge(as.data.frame(test), as.data.frame(CLINICAL), by="PATIENT"))

#SURVIVAL/PROCUREMENT
ggplot(test[GENDER=="FEMALE",], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL)) + geom_point() + scale_y_log10()
cor.test(test[GENDER=="FEMALE",]$DIST.TO.NORMAL, test[GENDER=="FEMALE",]$OVERALL.SURVIVAL, method="spearman")

ggplot(test[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=VITAL.STATUS.x)) + geom_point() + scale_y_log10() + facet_wrap(~VITAL.STATUS.x)
cor.test(test[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS.x=="Alive",]$DIST.TO.NORMAL,
    test[GENDER=="FEMALE",][OVERALL.SURVIVAL>0,][VITAL.STATUS.x=="Alive",]$OVERALL.SURVIVAL,
    method="spearman")
cor.test(as.vector(test[OVERALL.SURVIVAL>0,][VITAL.STATUS.x=="Dead",]$DIST.TO.NORMAL),
    as.vector(test[OVERALL.SURVIVAL>0,][VITAL.STATUS.x=="Dead",]$OVERALL.SURVIVAL),
    method="spearman")
ggplot(test[GENDER=="FEMALE",][OVERALL.SURVIVAL>200,], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=VITAL.STATUS.x)) + geom_point() + facet_wrap(~VITAL.STATUS.x)

ggplot(test[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,][!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),], aes(DIST.TO.NORMAL, OVERALL.SURVIVAL, colour=VITAL.STATUS.x)) + geom_point() + scale_y_log10()  + facet_grid(VITAL.STATUS.x~TUMOR.STAGE.SET)
ggplot(test[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,][!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),], aes(N.MUTATIONS, OVERALL.SURVIVAL, colour=VITAL.STATUS.x)) + geom_point() + scale_y_log10() + facet_grid(VITAL.STATUS.x~TUMOR.STAGE.SET)


ggplot(test[GENDER=="FEMALE",][OVERALL.SURVIVAL!=0,], aes(N.MUTATIONS, OVERALL.SURVIVAL, colour=VITAL.STATUS.x)) + geom_point() + scale_y_log10() +scale_x_log10()+ facet_wrap(~VITAL.STATUS.x)
ggplot(test, aes(N.MUTATIONS, OVERALL.SURVIVAL)) + geom_point() + scale_x_log10()

ggplot(test[DAYS.TO.PROCUREMENT!=0,], aes(DIST.TO.NORMAL, DAYS.TO.PROCUREMENT,colour=VITAL.STATUS.x)) + geom_point() + scale_y_log10()+ facet_wrap(~VITAL.STATUS.x)
cor.test(test[DAYS.TO.PROCUREMENT!=0,][VITAL.STATUS.x=="Alive",]$DAYS.TO.PROCUREMENT,
    test[DAYS.TO.PROCUREMENT!=0,][VITAL.STATUS.x=="Alive",]$DIST.TO.NORMAL, 
    method="spearman")
cor.test(test[DAYS.TO.PROCUREMENT!=0,][VITAL.STATUS.x=="Dead",]$DAYS.TO.PROCUREMENT,
    test[DAYS.TO.PROCUREMENT!=0,][VITAL.STATUS.x=="Dead",]$DIST.TO.NORMAL, 
    method="spearman")
cor.test(test[DAYS.TO.PROCUREMENT!=0,]$DAYS.TO.PROCUREMENT,
    test[DAYS.TO.PROCUREMENT!=0,]$DIST.TO.NORMAL, 
    method="spearman")

ggplot(test, aes(N.MUTATIONS, DAYS.TO.PROCUREMENT)) + geom_point() 
ggplot(test, aes(DIST.TO.NORMAL, TCGA.SURVIVAL)) + geom_point() + scale_y_log10()

#STAGES
ggplot(test[!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),], aes(TUMOR.STAGE.SET, DIST.TO.NORMAL)) + geom_boxplot()
ggplot(test[!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),],
    aes(DIST.TO.NORMAL, N.MUTATIONS, colour=TUMOR.STAGE.SET)) + geom_point() + scale_y_log10() + facet_wrap(~TUMOR.STAGE.SET)
cor.test(test[!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),][TUMOR.STAGE.SET=="Stage IV",]$DIST.TO.NORMAL,
    test[!(TUMOR.STAGE.SET %in% c("[Not Available]", "Stage Tis", "Stage X")),][TUMOR.STAGE.SET=="Stage IV",]$N.MUTATIONS, method="spearman")
data.table(STAGE=c("I", "II", "III","IV"), P.VAL=c(3.684e-05,0, 7.883e-09, 0.3264),
    RHO=c(0.3163255,0.3698251, 0.3800364, 0.2721585))

#LYMPH NODE PERCENTAGE
ggplot(test[!(is.na(PERCENT.LYMPH.NODES)),], aes(DIST.TO.NORMAL, PERCENT.LYMPH.NODES)) + geom_point() 

#100114
##REMOVE MALE PATIENTS FROM EXP MATRIX (NORMAL AND CANCER) PRIOR TO DIST TO MEDIAN NORMAL CALCULATION!!
BRCA.CORR
heatmap.2(BRCA.CORR, trace="none", scale="none", 
    ColSideColors=c(rep("green", length(BRCA.EXP.OBJ$normal.patients)),
        rep("red", length(BRCA.EXP.OBJ$cancer.patients))))

BRCA.CORR.FEMALE<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[GENDER=="FEMALE",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[GENDER=="FEMALE",]$PATIENT)]]
heatmap.2(BRCA.CORR.FEMALE, trace="none", scale="none", 
    ColSideColors=c(rep("green", 96),
        rep("red", 1080-96)))

#Do correlation heatmap per race to see if we can get a common clustered normal
CLINICAL.RACE<-CLINICAL[!(RACE %in% c("[Not Available]", "[Not Evaluated]")),]

#For White
BRCA.CORR.WHITE<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="WHITE",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="WHITE",]$PATIENT)]]
heatmap.2(BRCA.CORR.WHITE, trace="none", scale="none", 
    ColSideColors=c(rep("green", 103),
        rep("red", 816-103)))

#For White and Female only
BRCA.CORR.WHITE.FEMALE<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="WHITE" & GENDER=="FEMALE",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="WHITE" & GENDER=="FEMALE",]$PATIENT)]]
heatmap.2(BRCA.CORR.WHITE.FEMALE, trace="none", scale="none", 
    ColSideColors=c(rep("green", 102),
        rep("red", 807-102)))

#For Asian - Unfortunately only 1 normal
BRCA.CORR.ASIAN<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="ASIAN",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="ASIAN",]$PATIENT)]]
heatmap.2(BRCA.CORR.ASIAN, trace="none", scale="none", 
    ColSideColors=c(rep("green", 1),
        rep("red", 58-1)))

#For Black
BRCA.CORR.BLACK<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="BLACK OR AFRICAN AMERICAN",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="BLACK OR AFRICAN AMERICAN",]$PATIENT)]]
heatmap.2(BRCA.CORR.BLACK, trace="none", scale="none", 
    ColSideColors=c(rep("green", 6),
        rep("red", 134-6)))

#For Black and Female only
BRCA.CORR.BLACK.FEMALE<-BRCA.CORR[rownames(BRCA.CORR)[rownames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="BLACK OR AFRICAN AMERICAN" & GENDER=="FEMALE",]$PATIENT)],
colnames(BRCA.CORR)[colnames(BRCA.CORR) %in% as.vector(CLINICAL[RACE=="BLACK OR AFRICAN AMERICAN" & GENDER=="FEMALE",]$PATIENT)]]
heatmap.2(BRCA.CORR.BLACK.FEMALE, trace="none", scale="none", 
    ColSideColors=c(rep("green", 6),
        rep("red", 133-6)))