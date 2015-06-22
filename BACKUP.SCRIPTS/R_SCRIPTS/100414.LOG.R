#100414.LOG.R

#####100414####

#LOAD PACKAGES
library(ggplot2)
library(data.table)

#SETUP WD
setwd("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS")

#OBJECTS
BRCA.EXP.OBJ<-readRDS("BRCA/091514.CANCER.MATRICES.NORMALIZED.OBJ.rds")
table.1<-readRDS("BRCA/100514.BRCA.Table.1.rds")
CLINICAL<-readRDS("BRCA/093014.BRCA.CLINICAL.rds")
PAIRED.IDS<-substr(BRCA.EXP.OBJ$normal.patients,1,12)
BRCA.DIFF.EXP<-Function.RNAseq.Differential.Expression.V2(BRCA.EXP.OBJ, BRCA.EXP.OBJ$cancer.patients)
table.1.pval<-readRDS("BRCA/100514.BRCA.Table.1.bmr.rds")
brca.cnv.table<-readRDS("BRCA/092314.BRCA.GISTIC.FILTERED.TH.1.rds")

###### WORK #######

#@IDEA 1
#Use hclust to cluster patients based on expression patters pass a logFC threshold and use it to apply wilcoxon on significant genes versus rest to find encirhd mutated genes
table.1.pval

BRCA.DIST.1.5<-dist(t(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[logFC>=2.0,]$ID),BRCA.EXP.OBJ$cancer.patients]), method="euclidean")
BRCA.HCLUST.1.5<-hclust(BRCA.DIST.1.5, method="ward")
BRCA.HCLUST.1.5.GROUPS<-as.data.table(cutree(BRCA.HCLUST.1.5, k=5), keep.rownames=T)
BRCA.HCLUST.1.5.GROUPS<-BRCA.HCLUST.1.5.GROUPS[order(V2),]

HCLUST.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)) ,as.vector(BRCA.HCLUST.1.5.GROUPS$V2)))
rownames(HCLUST.ANN)<-c(BRCA.EXP.OBJ$normal.patients, as.vector(BRCA.HCLUST.1.5.GROUPS$V1))

pheatmap(BRCA.EXP.OBJ$combined.matrices[as.vector(BRCA.DIFF.EXP[abs(logFC)>2.0,]$ID),
     c(BRCA.EXP.OBJ$normal.patients, as.vector(BRCA.HCLUST.1.5.GROUPS$V1) )],
      scale="none", annotation=HCLUST.ANN, cluster_cols=FALSE
      ,
      filename="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/FIGURES/100414.BRCA.HCLUST.10.DIFF.EXP.2.0.jpeg", width=10, height=10) 

#TEST survival of samples with and without TTN
TTN.PATIENTS<-as.vector(table.1.pval[Hugo_Symbol=="TTN",]$PATIENT)
TP53.PATIENTS<-as.vector(table.1.pval[Hugo_Symbol=="TP53",]$PATIENT)
PIK3CA.PATIENTS<-as.vector(table.1.pval[Hugo_Symbol=="PIK3CA",]$PATIENT)
MUC4.PATIENTS<-as.vector(table.1.pval[Hugo_Symbol=="MUC4",]$PATIENT)
test<-table.1.pval[PATIENT %in% setdiff(union(TP53.PATIENTS, MUC4.PATIENTS), intersect(TP53.PATIENTS, MUC4.PATIENTS)),][Hugo_Symbol %in% c("MUC4", "TP53"),]

test<-test[Missense!=0,][,c("PATIENT","Hugo_Symbol","P.VAL.ADJ"),with=F]
test.clinical<-CLINICAL[,c("PATIENT", "GENDER","VITAL.STATUS","DAYS.TO.LAST.CONTACT", "OVERALL.SURVIVAL"), with=F]
test.clinical$status<-ifelse(test.clinical$DAYS.TO.LAST.CONTACT=="[Not Available]", 1, 2)
test<-as.data.table(merge(as.data.frame(test), as.data.frame(test.clinical)))

library(survival)
library(GGally)
test.survival<-survfit(Surv(OVERALL.SURVIVAL, status)~Hugo_Symbol, data=test[P.VAL.ADJ<1.05,])
ggsurv(test.survival) + theme(legend.position="bottom") 
survdiff(Surv(OVERALL.SURVIVAL, status)~Hugo_Symbol, data=test[P.VAL.ADJ<1.05,])

########WORK ON PROPOSAL FOR GRANT - METABOLITES########
table.2<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
setnames(table.2, c("METABOLITE", "KEGG_ID", "Hugo_Symbol"))
table.3<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.3.rds")
setnames(table.3, c("Hugo_Symbol", "KEGG_ID"))

table.v.pval<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071814.BRCA.Table.v.pval.rds")
table.u.pval<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071814.BRCA.Table.u.rds")

#Load BRCA step.01 file
BRCA.MUT<-as.data.table(read.csv("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/SOFTWARE/pipeline/step.01/100614.BRCA.mut", header=F, sep="\t", stringsAsFactors=F))
setnames(BRCA.MUT, c("ISOFORM", "PATIENT", "AA.POSITION", "WT", "MT"))
BRCA.MUT$PATIENT<-sapply(BRCA.MUT$PATIENT, function(x) paste0(strsplit(x,"-")[[1]][1:4], collapse="."))
BRCA.MUT$Hugo_Symbol<-sapply(BRCA.MUT$ISOFORM, function(x) strsplit(x, ".0")[[1]][1])

#Get significant metabolites
V.KEGG.SIG<-as.vector(table.v.pval[ADJ.P.VAL<0.05,]$KEGG_ID)
U.KEGG.SIG<-as.vector(table.u.pval$STATS[WILCOXON.P.ADJ<0.05,]$KEGG_ID)

###Work on Clinical per Pathway### - RESULTS
PATHWAY.TABLE.1<-copy(table.1.pval)

#STEROID HORMONE BIOSYNTHESIS - 75 PATIENTS AFFECTED
STEREOID.KEGG<-table.2[KEGG_ID %in% c("C01780", "C02140","C05284", "C00735"),]

STEREOID.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C01780", "C02140","C05284", "C00735"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

PATHWAY.TABLE.1$STEROID<-PATHWAY.TABLE.1$PATIENT %in% STEREOID.PATIENTS
STEROID.TABLE<-unique(PATHWAY.TABLE.1[,c("PATIENT", "STEROID"), with=F])
STEROID.TABLE<-as.data.table(merge(as.data.frame(STEROID.TABLE), as.data.frame(test.clinical)))

STEROID.SURVIVAL<-survfit(Surv(OVERALL.SURVIVAL,status)~STEROID, data=STEROID.TABLE)
ggsurv(STEROID.SURVIVAL) + theme(legend.position="bottom")
plot(STEROID.SURVIVAL, col= c("red", "black"), main="Steroid hormone biosynthesis - Survival Curves")
title(xlab="Days", ylab="Overall Survival (Probability)") 
legend(locator(1),legend=c("Steroid Mutations","Rest"), fill=c("red","black"))
a<-survdiff(Surv(OVERALL.SURVIVAL, status)~STEROID, data=STEROID.TABLE)

#PHOSPHOTIDYLINOSITOL SIGNALING SYSTEM - 697 PATIENTS AFFECTED
PS.KEGG<-table.2[KEGG_ID %in% c("C01194","C05981","C00165"),]

PS.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C01194","C05981","C00165"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

PATHWAY.TABLE.1$PS<-PATHWAY.TABLE.1$PATIENT %in% PS.PATIENTS
PS.TABLE<-unique(PATHWAY.TABLE.1[,c("PATIENT", "PS"), with=F])
PS.TABLE<-as.data.table(merge(as.data.frame(PS.TABLE), as.data.frame(test.clinical)))

PS.SURVIVAL<-survfit(Surv(OVERALL.SURVIVAL,status)~PS, data=PS.TABLE)
ggsurv(PS.SURVIVAL) + theme(legend.position="bottom")
plot(PS.SURVIVAL, col= c("red", "black"), main="Phosphatidylinositol signaling system - Survival Curves")
title(xlab="Days", ylab="Overall Survival (Probability)") 
legend(locator(1),legend=c("Phosphatidylinositol Mutations","Rest"), fill=c("red","black"))
survdiff(Surv(OVERALL.SURVIVAL, status)~PS, data=PS.TABLE)

#PLATELET ACTIVATION (NO ATP - C00002) - 684 PATIENTS AFFECTED
PLATELET.KEGG<-table.2[KEGG_ID %in% c("C00165","C05981"),]

PLATELET.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C00165","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

PATHWAY.TABLE.1$PLATELET<-PATHWAY.TABLE.1$PATIENT %in% PLATELET.PATIENTS
PLATELET.TABLE<-unique(PATHWAY.TABLE.1[,c("PATIENT", "PLATELET"), with=F])
PLATELET.TABLE<-as.data.table(merge(as.data.frame(PLATELET.TABLE), as.data.frame(test.clinical)))

PLATELET.SURVIVAL<-survfit(Surv(OVERALL.SURVIVAL,status)~PLATELET, data=PLATELET.TABLE)
ggsurv(PLATELET.SURVIVAL) + theme(legend.position="bottom")
plot(PLATELET.SURVIVAL, col= c("red", "black"), main="Platelet activation - Survival Curves")
title(xlab="Days", ylab="Overall Survival (Probability)") 
legend(locator(1),legend=c("Platelet Mutations","Rest"), fill=c("red","black"))
survdiff(Surv(OVERALL.SURVIVAL, status)~PLATELET, data=PLATELET.TABLE)

#ALDOSTERONE-REGULATED SODIUM REABSORPTION - 536 PATIENTS AFFECTED
ALDOSTERONE.KEGG<-table.2[KEGG_ID %in% c("C00735","C01780","C05981"),]

ALDOSTERONE.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C00735","C01780","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

PATHWAY.TABLE.1$ALDOSTERONE<-PATHWAY.TABLE.1$PATIENT %in% ALDOSTERONE.PATIENTS
ALDOSTERONE.TABLE<-unique(PATHWAY.TABLE.1[,c("PATIENT", "ALDOSTERONE"), with=F])
ALDOSTERONE.TABLE<-as.data.table(merge(as.data.frame(ALDOSTERONE.TABLE), as.data.frame(test.clinical)))

ALDOSTERONE.SURVIVAL<-survfit(Surv(OVERALL.SURVIVAL,status)~ALDOSTERONE, data=ALDOSTERONE.TABLE)
ggsurv(ALDOSTERONE.SURVIVAL) + theme(legend.position="bottom")
plot(ALDOSTERONE.SURVIVAL, col= c("red", "black"), main="Aldosterone-regulated sodium reabsorption - Survival Curves")
title(xlab="Days", ylab="Overall Survival (Probability)") 
legend(locator(1),legend=c("Aldosterone Mutations","Rest"), fill=c("red","black"))
survdiff(Surv(OVERALL.SURVIVAL, status)~ALDOSTERONE, data=ALDOSTERONE.TABLE)

###Work on Expression per Pathway### - NO RESULTS
PATHWAY.TABLE.1
table.3
BRCA.EXP

#STEROID HORMONE BIOSYNTHESIS - 75 PATIENTS AFFECTED
STEREOID.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C01780", "C02140","C05284", "C00735"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))
STEREOID.HUGOS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C01780", "C02140","C05284", "C00735"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$Hugo_Symbol))

EXP.STEREOID.PATIENTS<-intersect(STEREOID.PATIENTS, BRCA.EXP.OBJ$cancer.patients)
EXP.STEREOID.HUGOS<-intersect(STEREOID.HUGOS, rownames(BRCA.EXP.OBJ$combined.matrices))

STEREOID.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)) ,
    rep("STEROID", length(EXP.STEREOID.PATIENTS)),
    rep("NON.STEREOID", length(setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.STEREOID.PATIENTS)))   ))
rownames(STEREOID.ANN)<-c(BRCA.EXP.OBJ$normal.patients, 
    EXP.STEREOID.PATIENTS,
    setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.STEREOID.PATIENTS))

pheatmap(BRCA.EXP.OBJ$combined.matrices[EXP.STEREOID.HUGOS, c(BRCA.EXP.OBJ$normal.patients,
    EXP.STEREOID.PATIENTS, setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.STEREOID.PATIENTS))], scale="none", annotation=STEREOID.ANN,cluster_cols=FALSE)

#PHOSPHOTIDYLINOSITOL SIGNALING SYSTEM - 697 PATIENTS AFFECTED
PS.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C01194","C05981","C00165"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))
PS.HUGOS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C01194","C05981","C00165"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$Hugo_Symbol))

EXP.PS.PATIENTS<-intersect(PS.PATIENTS, BRCA.EXP.OBJ$cancer.patients)
EXP.PS.HUGOS<-intersect(PS.HUGOS, rownames(BRCA.EXP.OBJ$combined.matrices))

PS.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)) ,
    rep("PS", length(EXP.PS.PATIENTS)),
    rep("NON.PS", length(setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PS.PATIENTS)))   ))

rownames(PS.ANN)<-c(BRCA.EXP.OBJ$normal.patients, EXP.PS.PATIENTS,
    setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PS.PATIENTS))

pheatmap(BRCA.EXP.OBJ$combined.matrices[EXP.PS.HUGOS, c(BRCA.EXP.OBJ$normal.patients,
    EXP.PS.PATIENTS, setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PS.PATIENTS))], scale="none", annotation=PS.ANN,cluster_cols=FALSE)

#PLATELET ACTIVATION (NO ATP - C00002) - 684 PATIENTS AFFECTED
PLATELET.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C00165","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))
PLATELET.HUGOS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.3[KEGG_ID %in% c("C00165","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$Hugo_Symbol))

EXP.PLATELET.PATIENTS<-intersect(PLATELET.PATIENTS, BRCA.EXP.OBJ$cancer.patients)
EXP.PLATELET.HUGOS<-intersect(PLATELET.HUGOS, rownames(BRCA.EXP.OBJ$combined.matrices))

PLATELET.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)) ,
    rep("PLATELET", length(EXP.PLATELET.PATIENTS)),
    rep("NON.PLATELET", length(setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PLATELET.PATIENTS)))   ))

rownames(PLATELET.ANN)<-c(BRCA.EXP.OBJ$normal.patients, EXP.PLATELET.PATIENTS,
    setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PLATELET.PATIENTS))

pheatmap(BRCA.EXP.OBJ$combined.matrices[EXP.PLATELET.HUGOS, c(BRCA.EXP.OBJ$normal.patients,
    EXP.PLATELET.PATIENTS, setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.PLATELET.PATIENTS))], scale="none", annotation=PLATELET.ANN,cluster_cols=FALSE)

#ALDOSTERONE-REGULATED SODIUM REABSORPTION - 536 PATIENTS AFFECTED
ALDOSTERONE.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C00735","C01780","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))
ALDOSTERONE.HUGOS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% c("C00735","C01780","C05981"),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$Hugo_Symbol))

EXP.ALDOSTERONE.PATIENTS<-intersect(ALDOSTERONE.PATIENTS, BRCA.EXP.OBJ$cancer.patients)
EXP.ALDOSTERONE.HUGOS<-intersect(ALDOSTERONE.HUGOS, rownames(BRCA.EXP.OBJ$combined.matrices))

ALDOSTERONE.ANN<-data.frame(TYPE=c(rep("NORMAL", length(BRCA.EXP.OBJ$normal.patients)) ,
    rep("ALDOSTERONE", length(EXP.ALDOSTERONE.PATIENTS)),
    rep("NON.ALDOSTERONE", length(setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.ALDOSTERONE.PATIENTS)))   ))

rownames(ALDOSTERONE.ANN)<-c(BRCA.EXP.OBJ$normal.patients, EXP.ALDOSTERONE.PATIENTS,
    setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.ALDOSTERONE.PATIENTS))

pheatmap(BRCA.EXP.OBJ$combined.matrices[EXP.ALDOSTERONE.HUGOS, c(BRCA.EXP.OBJ$normal.patients,
    EXP.ALDOSTERONE.PATIENTS, setdiff(BRCA.EXP.OBJ$cancer.patients , EXP.ALDOSTERONE.PATIENTS))], scale="none", annotation=ALDOSTERONE.ANN,cluster_cols=FALSE)

###Work on Significant KEGG Mutation Binding Scores### - RESULTS

#Load alignment Hugo to PDB file
ALN.PDB<-as.data.table(read.csv("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/SOFTWARE/pipeline/step.02/refSeqsAgainstBioLip.tab.parsed", head=T, sep="\t", stringsAsFactors=F))
ALN.PDB$Hugo_Symbol<-sapply(ALN.PDB$QUERY, function(x) strsplit(x, ".0")[[1]][1])

table.2
table.2$HAS.PDB<-table.2$Hugo_Symbol %in% unique(as.vector(ALN.PDB$Hugo_Symbol))
sum(unique(table.2[,3:4,with=F])$HAS.PDB) #2345 out of 3897 have PDB structures (60%)

#Look for coverage of significant KEGG - Remove ATP!!
U.KEGG.SIG
V.KEGG.SIG
KEGG.SIG<-c(U.KEGG.SIG, V.KEGG.SIG)
KEGG.SIG<-KEGG.SIG[KEGG.SIG!="C00002"]
table.2.count<-table.2[,list(COUNT=length(Hugo_Symbol)), by="KEGG_ID"]

KEGG.SIG.PDB.COV<-as.data.table(sapply(KEGG.SIG, function(x) sum(as.vector(table.2[KEGG_ID==x,]$HAS.PDB))), keep.rownames=T)
setnames(KEGG.SIG.PDB.COV, c("KEGG_ID", "PDB.COVERAGE"))
KEGG.SIG.PDB.COV<-join(KEGG.SIG.PDB.COV, table.2.count, by="KEGG_ID")

#Load Dario's aa mutation table
#GET P-VALUE PER PATIENT USING HYPERGEOMETRIC - AND PER METABOLITE???? 
BRCA.MUT
table.2.KEGG<-table.2[KEGG_ID %in% KEGG.SIG,]

#Load structural files into single table
dario.structural.file<-as.data.table(read.csv("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/SOFTWARE/pipeline/step.04/structFiles", header=F, stringsAsFactors=F))

dario.structural.table<-lapply(dario.structural.file$V1, function(x) {
    target.struct<-as.data.table(read.csv(x, header=F, sep="\t", stringsAsFactors=F))
    ISOFORM<-strsplit(strsplit(x, "BINDID90/")[[1]][2], ".struct")[[1]][1]
    target.struct$ISOFORM<-ISOFORM
    return(target.struct)
})
dario.structural.table<-do.call(rbind, dario.structural.table)
setnames(dario.structural.table, c("AA.POSITION", "WT", "TYPE","ISOFORM"))
dario.structural.count<-dario.structural.table[,list(STRUCT=sum(TYPE=="S"),
    UNKNOWN=sum(TYPE=="U"), BINDING=sum(TYPE=="B"), ALL=length(TYPE)), by="ISOFORM"]
dario.structural.count
test<-copy(dario.structural.count)
test$Hugo_Symbol<-sapply(test$ISOFORM, function(x) strsplit(x,".0")[[1]][1])
length(unique(as.vector(test$Hugo_Symbol))) #3032

#Obtain genes affected by sifinicant KEGG_IS and define their TYPE
BRCA.MUT.KEGG<-as.data.table(merge(as.data.frame(BRCA.MUT), as.data.frame(table.2.KEGG), by="Hugo_Symbol"))
BRCA.MUT.KEGG<-as.data.table(merge(as.data.frame(BRCA.MUT.KEGG),as.data.frame(dario.structural.table)))

ACTUAL.KEGG.SIG.COV<-unique(BRCA.MUT.KEGG[,c("Hugo_Symbol","KEGG_ID"), with=F])[,list(PDB.COVERAGE=length(Hugo_Symbol)), by="KEGG_ID"]
ACTUAL.KEGG.SIG.COV$COUNT<-c(272,80,14,9,9,7, 7, 10, 25)

#Test significant metabolites by hypergeometric (using isoforms)
BRCA.MUT.KEGG<-as.data.table(merge(as.data.frame(BRCA.MUT.KEGG), as.data.frame(dario.structural.count)))
BRCA.MUT.KEGG.PVAL<-BRCA.MUT.KEGG[,list(P.VAL=phyper(q=sum(TYPE=="B")-1, m=sum(BINDING), n=sum(ALL)-sum(BINDING), k=length(TYPE), lower.tail=F)), by="KEGG_ID"]

BRCA.MUT.KEGG.PVAL$P.VAL.ADJ<-p.adjust(BRCA.MUT.KEGG.PVAL$P.VAL, method="fdr")

#Test survival on those found by enrichment
ENRICHED.KEGG<-table.2[KEGG_ID %in% as.vector(BRCA.MUT.KEGG.PVAL[P.VAL.ADJ<0.05,]$KEGG_ID),]

ENRICHED.PATIENTS<-unique(as.vector(PATHWAY.TABLE.1[Hugo_Symbol %in% table.2[KEGG_ID %in% as.vector(BRCA.MUT.KEGG.PVAL[P.VAL.ADJ<0.05,]$KEGG_ID),]$Hugo_Symbol,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

PATHWAY.TABLE.1$ENRICHED<-PATHWAY.TABLE.1$PATIENT %in% ENRICHED.PATIENTS
ENRICHED.TABLE<-unique(PATHWAY.TABLE.1[,c("PATIENT", "ENRICHED"), with=F])
ENRICHED.TABLE<-as.data.table(merge(as.data.frame(ENRICHED.TABLE), as.data.frame(test.clinical)))

ENRICHED.SURVIVAL<-survfit(Surv(OVERALL.SURVIVAL,status)~ENRICHED, data=ENRICHED.TABLE)
ggsurv(ENRICHED.SURVIVAL) + theme(legend.position="bottom")
survdiff(Surv(OVERALL.SURVIVAL, status)~ENRICHED, data=ENRICHED.TABLE)

#Get P.VALUES for all KEGG, regardless of significant
#Obtain genes affected by sifinicant KEGG_IS and define their TYPE
BRCA.MUT.KEGG.ALL<-as.data.table(merge(as.data.frame(BRCA.MUT), as.data.frame(table.2), by="Hugo_Symbol"))
BRCA.MUT.KEGG.ALL<-as.data.table(merge(as.data.frame(BRCA.MUT.KEGG.ALL),as.data.frame(dario.structural.table)))

ACTUAL.KEGG.SIG.COV.ALL<-unique(BRCA.MUT.KEGG.ALL[,c("Hugo_Symbol","KEGG_ID"), with=F])[,list(PDB.COVERAGE=length(Hugo_Symbol)), by="KEGG_ID"]

#Test significant metabolites by hypergeometric (using isoforms)
BRCA.MUT.KEGG.ALL<-as.data.table(merge(as.data.frame(BRCA.MUT.KEGG.ALL), as.data.frame(dario.structural.count)))
BRCA.MUT.KEGG.ALL.PVAL<-BRCA.MUT.KEGG.ALL[,list(P.VAL=phyper(q=sum(TYPE=="B")-1, m=sum(BINDING), n=sum(ALL)-sum(BINDING), k=length(TYPE), lower.tail=F)), by="KEGG_ID"]

BRCA.MUT.KEGG.ALL.PVAL$P.VAL.ADJ<-p.adjust(BRCA.MUT.KEGG.ALL.PVAL$P.VAL, method="fdr")
BRCA.MUT.KEGG.ALL.PVAL[P.VAL.ADJ<0.05,]

#Develop function for pathway enrichment
PATH.FILE<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND"
KEGGS<-copy(KEGG.SIG)
KEGGS<-c(KEGGS, "C00002")
universe.keggs<-unique(c(as.vector(table.2$KEGG_ID), as.vector(table.3$KEGG_ID)))
KEGG.PATH.ENRICHMENT<-function(KEGGS, PATH.FILE, universe.keggs){
    #Look for enrichment of kegg IDs in each pathway
    #PATHWAY.FILE in the form "050314_PATHWAY_TO_COMPOUND"
    #universe.keggs are all the keggs that we originially started with
    #NOTE - Filtering for KEGGs found in 2 pathways or more

    require(data.table)

    #Load pathway table
    KEGG.PATHWAY<-as.data.table(read.csv(PATH.FILE, sep="\t", header=T, stringsAsFactors=F))

    #Filter pathway file for those kegg that actually exist in our universe
    KEGG.PATHWAY<-KEGG.PATHWAY[COMPOUND %in% universe.keggs,]

    #Calculate enrichment
    KEGG.P.VAL<-KEGG.PATHWAY[,list(P.VAL=
        phyper(q=sum(KEGGS %in% COMPOUND)-1,
            m=length(COMPOUND),
            n=length(universe.keggs)-length(COMPOUND),
            k=length(KEGGS), lower.tail=F)
        , KEGG.IN.PATH=sum(KEGGS %in% COMPOUND)),
         by=c("PATHWAY", "DESCRIPTION")]

    #Adjust for multiple hypothesis correction
    KEGG.P.VAL$P.VAL.ADJ<-p.adjust(KEGG.P.VAL$P.VAL, method="fdr")

    #Clean up and return
    KEGG.P.VAL<-KEGG.P.VAL[KEGG.IN.PATH>2,]
    KEGG.P.VAL<-KEGG.P.VAL[order(P.VAL.ADJ),]
    return(KEGG.P.VAL)
}

KEGG.P.VAL[P.VAL.ADJ<0.05,][,c(2,4,5),with=F]

#Check rest found by significant binding alone
KEGG.P.VAL.ALL<-KEGG.PATH.ENRICHMENT(as.vector(BRCA.MUT.KEGG.ALL.PVAL[P.VAL.ADJ<0.05,]$KEGG_ID),
    "/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND")
KEGG.P.VAL.ALL[P.VAL.ADJ<0.05,]

#####Check for Enrichment on Metabolic Pathways based on Wilcoxon tables#####
wilcoxon.table
HUGOS<-as.vector(wilcoxon.table[P.VAL.ADJ<0.05,]$Hugo_Symbol)
universe.hugos<-unique(as.vector(table.1$table.1$Hugo_Symbol))
enzyme.product.file<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/KEGG/061214_ENZYME_PRODUCT"

HUGO.PATH.ENRICHMENT<-function(HUGOS, PATH.FILE, enzyme.product.file, universe.hugos){
    #Look for enrichment of Hugos in each pathway
    #PATHWAY.FILE in the form "050314_PATHWAY_TO_COMPOUND"
    #enyzme.product.file in the form "061214_ENZYME_PRODUCT"

    require(data.table)

    #Load pathway table
    KEGG.PATHWAY<-as.data.table(read.csv(PATH.FILE, sep="\t", header=T, stringsAsFactors=F))

    #Load enzyme to KEGG table
    KEGG.ENZYME<-as.data.table(read.csv(enzyme.product.file, sep="\t", header=T, stringsAsFactors=F))
    KEGG.ENZYME<-unique(KEGG.ENZYME[,1:2,with=F])

    #Make Hugo to Pathway table [PATHWAY, DESCRIPTION, Hugo_Symbol]
    HUGO.PATHWAY<-as.data.table(merge(as.data.frame(KEGG.PATHWAY), as.data.frame(KEGG.ENZYME), by.x="COMPOUND", by.y="KEGG_ID"))
    HUGO.PATHWAY<-unique(HUGO.PATHWAY[,c("PATHWAY", "DESCRIPTION", "Enzyme"), with=F])
    setnames(HUGO.PATHWAY, c("PATHWAY", "DESCRIPTION","Hugo_Symbol"))

    #Filter pathway table for those that are actually found in our universe
    HUGO.PATHWAY<-HUGO.PATHWAY[Hugo_Symbol %in% universe.hugos,]

    #Calculate enrichment
    HUGO.P.VAL<-HUGO.PATHWAY[,list(P.VAL=
        phyper(q=sum(Hugo_Symbol %in% HUGOS)-1,
            m=length(Hugo_Symbol),
            n=length(universe.hugos)-length(Hugo_Symbol),
            k=length(HUGOS), lower.tail=F)
        , HUGOS.IN.PATH=sum(Hugo_Symbol %in% HUGOS)),
         by=c("PATHWAY", "DESCRIPTION")]

    #Adjust for multiple hypothesis correction
    HUGO.P.VAL$P.VAL.ADJ<-p.adjust(HUGO.P.VAL$P.VAL, method="fdr")

    #Clean up and return
    HUGO.P.VAL<-HUGO.P.VAL[HUGOS.IN.PATH>2,]
    HUGO.P.VAL<-HUGO.P.VAL[order(P.VAL.ADJ),]
    return(HUGO.P.VAL)
}   
HUGO.P.VAL[P.VAL.ADJ<0.05,]
HUGO.P.VAL[1:20,]

MUTSIG.KEGG.PATH<-HUGO.PATH.ENRICHMENT(as.vector(MUTSIG[q<0.05,]$gene), PATH.FILE, enzyme.product.file, as.vector(MUTSIG$gene) )
MUTSIG.KEGG.PATH[P.VAL.ADJ<0.05,]

length(unique(as.vector(table.1$table.1$PATIENT)))
length(unique(as.vector(table.1$table.1[Hugo_Symbol %in% as.vector(wilcoxon.table[P.VAL.ADJ<0.05 & Hugo_Symbol!="TTN",]$Hugo_Symbol),]$PATIENT)))

#FOR GM
wilcoxon.table #Results from Function.u.p.wilcoxon
library(ROCR)
WILCOXON.PRED<-prediction(as.vector(wilcoxon.table[P.VAL.ADJ<0.05,]$W), 
    as.vector(wilcoxon.table[P.VAL.ADJ<0.05,]$Hugo_Symbol) %in% as.vector(COSMIC.BRCA$Symbol)  )
WILCOXON.PERF<-performance(WILCOXON.PRED, "prec", "rec")
test.main<-data.table(x=WILCOXON.PERF@x.values[[1]],y=WILCOXON.PERF@y.values[[1]])
test.main.bmr<-data.table(x=WILCOXON.PERF@x.values[[1]],y=WILCOXON.PERF@y.values[[1]])
setnames(test.main, c("RECALL", "PRECISSION"))
setnames(test.main.bmr, c("RECALL", "PRECISSION"))
test.main$TYPE<-"WILCOX"
test.main.bmr$TYPE<-"BMR+WILCOX"
ggplot(rbind(test.main.bmr, test.main), aes(RECALL, PRECISSION, colour=TYPE))+geom_path()

#Load MUTSIG
MUTSIG<-as.data.table(read.csv("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/SOFTWARE/MUTSIG/MutSigCV_1.4/072014.BRCA.CURATED.RESULTS.sig_genes.txt", head=T, sep="\t", stringsAsFactors=F))
MUTSIG$SCORE<--log(MUTSIG$q)
MUTSIG$SCORE[is.infinite(MUTSIG$SCORE)]<-100
MUTSIG.PRED<-prediction( as.vector(MUTSIG$SCORE), 
    as.vector(MUTSIG$gene) %in% as.vector(COSMIC.BRCA$Symbol))
MUTSIG.PERF<-performance(MUTSIG.PRED, "prec","rec")
MUTSIG.MAIN<-data.table(x=MUTSIG.PERF@x.values[[1]],y=MUTSIG.PERF@y.values[[1]])
setnames(MUTSIG.MAIN, c("RECALL", "PRECISSION"))
MUTSIG.MAIN$TYPE<-"MUTSIG"

ggplot(rbind(test.main.bmr, test.main,MUTSIG.MAIN), aes(RECALL, PRECISSION, colour=TYPE))+geom_path()

#Load v.p
TABLE.V.P.PRED<-prediction( as.vector(table.v.p$v.PROTEIN), 
    as.vector(table.v.p$Hugo_Symbol) %in% as.vector(COSMIC.BRCA$Symbol))
TABLE.V.P.PERF<-performance(TABLE.V.P.PRED, "prec","rec")
TABLE.V.P.MAIN<-data.table(x=TABLE.V.P.PERF@x.values[[1]],y=TABLE.V.P.PERF@y.values[[1]])
setnames(TABLE.V.P.MAIN, c("RECALL", "PRECISSION"))
TABLE.V.P.MAIN$TYPE<-"EXPRESSION"

ggplot(rbind(test.main.bmr, test.main,MUTSIG.MAIN, TABLE.V.P.MAIN), aes(RECALL, PRECISSION, colour=TYPE))+geom_path()

########100814######
#Do enrichment of wilcoxon results on pathways
brca.bmr.w.path<-Function.path.enrichment( as.vector(test.2[P.VAL.ADJ<0.05,]$Hugo_Symbol) ,path, unique(as.vector(table.1$table.1$Hugo_Symbol)))
brca.bmr.w.path[P.VAL.ADJ<0.05,]

gbm.bmr.w.path<-Function.path.enrichment( as.vector(test.2.GBM[P.VAL.ADJ<0.05,]$Hugo_Symbol) ,path, unique(as.vector(GBM.table.1$table.1$Hugo_Symbol)))
gbm.bmr.w.path[P.VAL.ADJ<0.05,]

coad.bmr.w.path<-Function.path.enrichment( as.vector(test.2.COAD[P.VAL.ADJ<0.05,]$Hugo_Symbol) ,path, unique(as.vector(COAD.table.1$table.1$Hugo_Symbol)))
coad.bmr.w.path[P.VAL.ADJ<0.05,]

ov.bmr.w.path<-Function.path.enrichment( as.vector(test.2.OV[P.VAL.ADJ<0.05,]$Hugo_Symbol) ,path, unique(as.vector(OV.table.1$table.1$Hugo_Symbol)))
ov.bmr.w.path[P.VAL.ADJ<0.05,]

read.bmr.w.path<-Function.path.enrichment( as.vector(test.2.READ[P.VAL.ADJ<0.05,]$Hugo_Symbol) ,path, unique(as.vector(READ.table.1$table.1$Hugo_Symbol)))
read.bmr.w.path[P.VAL.ADJ<0.05,]

#Next, do patient precission recall based on wilcoxon[P.VAL.ADJ<0.05,]


####
KEGG.TARGET<-as.data.table(read.csv("/Users/jzamalloa/Desktop/FOLDER/GRANTS/NIH/CANCER/KEGG_TARGET.csv", header=T, sep=",", stringsAsFactors=F))
setnames(KEGG.TARGET, c("METABOLITE","KEGG_ID", "PATHWAY", "ID","TYPE"))

BRCA.DIFF.EXP[ID %in% KEGG.TARGET$ID, ][,c(1,5,7), with=F]
KEGG.TARGET<-as.data.table(merge(as.data.frame(KEGG.TARGET), as.data.frame(BRCA.DIFF.EXP[,c(1,5,7), with=F])))
unique(KEGG.TARGET[,c("ID", "logFC","adj.P.Val"), with=F])
KEGG.TARGET[PATHWAY=="T cell receptor signaling pathway",]

KEGG.PVAL.PATH<-KEGG.TARGET[,list(DIFF.EXPRESSED=sum(adj.P.Val<0.05),
    NON.DIFF.EXPRESSED=sum(adj.P.Val>=0.05)), by="PATHWAY"]
KEGG.PVAL.PATH<-melt(KEGG.PVAL.PATH, id.vars="PATHWAY")
setnames(KEGG.PVAL.PATH, c("PATHWAY", "DIFFERENTIAL.EXPRESSION", "GENES"))
ggplot(KEGG.PVAL.PATH, aes(PATHWAY, GENES, colour=DIFFERENTIAL.EXPRESSION)) + geom_histogram(stat="identity", position="dodge") +theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Function to assign p-value on survival of patients for kegg pathway related mutations
KEGG.PATHWAY
KEGG.P.VAL[P.VAL.ADJ<0.05,][,c(2,4,5),with=F]
table.1.pval
table.2
test.clinical
Function.KEGG.PATH.SURV<-function(KEGG.PATHWAY, INTEREST.PATHWAYS.IDS, table.2, table.1.pval, clinical.table, INTEREST.KEGG.IDS=c(), EXCLUDE.KEGG=c()){

    #Obtains survival fits for patients affected with mutations in each KEGG pathway containing target KEGG_IDs

    require(data.table)
    require(survival)
    require(reshape2)

    #Filter pathway table for target pathways
    KEGG.PATHWAY<-KEGG.PATHWAY[PATHWAY %in% INTEREST.PATHWAYS.IDS,]

    #Filter KEGG.PATHWAY for target KEGGS if desired, if not, use all
    if (length(INTEREST.KEGG.IDS)>0){
        KEGG.PATHWAY<-KEGG.PATHWAY[COMPOUND %in% INTEREST.KEGG.IDS,]
    }

    #Any desired KEGG.IDs to be excluded from dataset?
    if (length(EXCLUDE.KEGG)>0){
        KEGG.PATHWAY<-KEGG.PATHWAY[!(COMPOUND %in% EXCLUDE.KEGG),]
    }

    #Apply function per pathway
    internal.surv.function<-function(compounds, clinical.table){

        #Get genes related to compounds
        target.genes<-as.vector(table.2[KEGG_ID %in% compounds,]$Hugo_Symbol)

        #Get target patients - using bmr p.vals
        target.patients<-unique(as.vector(table.1.pval[Hugo_Symbol %in% target.genes,][Missense!=0 & P.VAL.ADJ<0.05,]$PATIENT))

        #Build table for survival fit
        surv.table<-unique(table.1.pval[,"PATIENT", with=F])
        surv.table$IN.PATHWAY<-surv.table$PATIENT %in% target.patients
        surv.table<-as.data.table(merge(as.data.frame(surv.table), as.data.frame(clinical.table)))

        #Perform fit
        surv.fit<-survdiff(Surv(OVERALL.SURVIVAL, status)~IN.PATHWAY, data=surv.table)

        #Get pvalue
        P.VAL<-1-pchisq(surv.fit$chisq, length(surv.fit$n)-1)

        return(list(P.VAL=P.VAL, N.PATIENTS=length(target.patients), N.COMPOUNDS=length(compounds)))
    }

    SURVIVAL.P.VAL<-KEGG.PATHWAY[,internal.surv.function(COMPOUND), by=c("PATHWAY", "DESCRIPTION")]

    #Clean up and return
    SURVIVAL.P.VAL$P.VAL.ADJ<-p.adjust(SURVIVAL.P.VAL$P.VAL, method="fdr")
    SURVIVAL.P.VAL<-SURVIVAL.P.VAL[order(P.VAL.ADJ),]
    return(SURVIVAL.P.VAL)

}

KEGG.SIG.TABLE<-data.table(KEGG_ID=KEGG.SIG, COUNT=1:13)
KEGG.SIG.TABLE<-unique(join(KEGG.SIG.TABLE, table.2[,1:2, with=F], by="KEGG_ID"))
KEGG.SIG.TABLE$COUNT<-NULL
KEGG.SIG.TABLE<-unique(KEGG.SIG.TABLE)
KEGG.SIG.TABLE<-KEGG.SIG.TABLE[,internal.surv.function(KEGG_ID, test.clinical), by=c("METABOLITE")]
KEGG.SIG.TABLE$P.VAL.ADJ<-p.adjust(KEGG.SIG.TABLE$P.VAL, method="fdr")
KEGG.SIG.TABLE[,c(1,2,5),with=F]
write.table(KEGG.SIG.TABLE[,c(1,2,5),with=F], file="/Users/jzamalloa/Desktop/FOLDER/GRANTS/NIH/CANCER/P.VALS.COMPOUND", quote=F, sep="\t", col.names=T, row.names=F)

b<-Function.KEGG.PATH.SURV(KEGG.PATHWAY, as.vector(KEGG.P.VAL[P.VAL.ADJ<0.05,]$PATHWAY), 
    table.2, table.1.pval, test.clinical, KEGG.SIG)
write.table(b[,c(2,3,5,6),with=F], file="/Users/jzamalloa/Desktop/FOLDER/GRANTS/NIH/CANCER/P.VALS", quote=F, sep="\t", col.names=T, row.names=F)

KEGG.PATHWAY[DESCRIPTION=="Carbohydrate digestion and absorption - Homo sapiens (human)",][COMPOUND %in% KEGG.SIG,]


a<-Function.KEGG.PATH.SURV(KEGG.PATHWAY, as.vector(KEGG.P.VAL[P.VAL.ADJ<0.05,]$PATHWAY), 
    table.2, table.1.pval, test.clinical)

a.b<-as.data.table(merge(as.data.frame(a[,c(2,5), with=F]), 
    as.data.frame(b[,c(2,5), with=F]), by=c("DESCRIPTION") ))
setnames(a.b, c("DESCRIPTION", "P.VAL.ADJ.ALL","P.VAL.ADJ.SIG"))

a.b$NEG.LOG.ALL<--log(a.b$P.VAL.ADJ.ALL)
a.b$NEG.LOG.SIG<--log(a.b$P.VAL.ADJ.SIG)
setnames(a.b, c("Pathway", "P.VAL.ADJ.ALL", "P.VAL.ADJ.SIG", "All Metabolites", "Significant Metabolites"))
a.b.melt<-melt(a.b[, c(1,4,5), with=F], id.vars="Pathway")

ggplot(a.b.melt, aes(Pathway, value, fill=variable)) + geom_histogram(stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_hline(aes(yintercept=-log(0.05), color="red")) + ylab("-logP")