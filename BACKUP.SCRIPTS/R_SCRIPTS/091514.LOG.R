#####TRY CLUSTERING BY PCAs######
library(pheatmap)
library(data.table)
library(ggplot2)
library(reshape2)
library(gplots)
setwd("Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS")

Function.heatmap.2.mut.color<-function(corr.matrix, table.1, target.gene){
    #Takes corr.matrix and produces a vector color for colnames based on given mutated gene in table.1
    #Patients with target gene will be colored in red while the rest will be in green

    GENE.COLOR<-as.character(colnames(corr.matrix) %in% as.vector(table.1[Hugo_Symbol==target.gene]$PATIENT))
    GENE.COLOR<-replace(GENE.COLOR, GENE.COLOR=="TRUE", "red")
    GENE.COLOR<-replace(GENE.COLOR, GENE.COLOR=="FALSE", "green")
    return(GENE.COLOR)
}

BRCA.EXP.OBJ<-readRDS("BRCA/082614.CANCER.MATRICES.NORMALIZED.OBJ.rds")
BRCA.EXP<-BRCA.EXP.OBJ$combined.matrices[,BRCA.EXP.OBJ$cancer.patients]
pheatmap(BRCA.EXP, scale="row")

#Do principal components on cancer matrix - TOO MUCH REDUCTION
BRCA.EXP.PCA<-prcomp(BRCA.EXP, scores=T, scale.=T)
names(BRCA.EXP.PCA)

plot(BRCA.EXP.PCA)
screeplot(BRCA.EXP.PCA, type="lines")
dim(BRCA.EXP.PCA$rotation)
head(BRCA.EXP.PCA$rotation[,1:2])
ggplot(as.data.frame(BRCA.EXP.PCA$rotation), aes(PC1, PC2)) + geom_point()

BRCA.EXP.PCA.COR<-cor(t(BRCA.EXP.PCA$rotation), method="spearman")
head(BRCA.EXP.PCA.COR[,1:3])
pheatmap(BRCA.EXP.PCA.COR, scale="none")

#######GO BACK TO MORE NOISY DATA, NON-BATCHED EFFECT REDUCED########
BRCA.CANCER.MATRICES<-readRDS("BRCA/061314.BRCA.CANCER.MATRICES")
names(BRCA.CANCER.MATRICES)
dim(BRCA.CANCER.MATRICES$tumor)
head(BRCA.CANCER.MATRICES$tumor[,1:3])

#First normalized RNAseq counts without batch effect reduction
BRCA.CANCER.MATRICES.NORM<-Function.RNAseq.Matrices.Normalization(BRCA.CANCER.MATRICES$normal,
    BRCA.CANCER.MATRICES$tumor, rm.batch.effect=F)
names(BRCA.CANCER.MATRICES.NORM)
dim(BRCA.CANCER.MATRICES.NORM$combined.matrices)
saveRDS(object=BRCA.CANCER.MATRICES.NORM, file="BRCA/091514.CANCER.MATRICES.NORMALIZED.OBJ.rds")

#Correlation on cancer patients
BRCA.CORR<-cor(BRCA.CANCER.MATRICES.NORM$combined.matrices[,BRCA.CANCER.MATRICES.NORM$cancer.patients
    ], method="spearman")
pheatmap(BRCA.CORR, scale="none")   

TP53.COLOR<-as.character(colnames(BRCA.CORR) %in% as.vector(table.1[Hugo_Symbol=="TP53"]$PATIENT))
TP53.COLOR<-replace(TP53.COLOR, TP53.COLOR=="TRUE", "red")
TP53.COLOR<-replace(TP53.COLOR, TP53.COLOR=="FALSE", "green")

TTN.COLOR<-as.character(colnames(BRCA.CORR) %in% as.vector(table.1[Hugo_Symbol=="TTN"]$PATIENT))
TTN.COLOR<-replace(TTN.COLOR, TTN.COLOR=="TRUE", "red")
TTN.COLOR<-replace(TTN.COLOR, TTN.COLOR=="FALSE", "green")

#test.patients<-sample(colnames(BRCA.CORR), 50)
heatmap.2(BRCA.CORR, trace="none", scale="none", 
    ColSideColors=Function.heatmap.2.mut.color(BRCA.CORR, table.1, "MYC"))

BRCA.DISS<-1-BRCA.CORR
BRCA.DISS<-sqrt(1-BRCA.CORR^2)
pheatmap(BRCA.DISS, scale="none")

#Apply diss.matrix to silhouette(p) and silohuette (j) functions - SOME FOUND IN S(j)
#FOUND 4 CANCER CENSUS
BRCA.S.j<-readRDS("BRCA/091514.BRCA.j.silhouette.median.rds")
names(BRCA.S.j)
BRCA.S.j$S.I.C.TABLE[order(S.I.CLUSTER,decreasing=T),]

BRCA.S.j.ai.count<-BRCA.S.j$A.I.TABLE[,list(CLUSTER.A.I=median(A.I)), by="KEGG_ID"]
BRCA.S.j.ai.count<-as.data.table(merge(as.data.frame(BRCA.S.j.ai.count), as.data.frame(patient.count)))
BRCA.S.j.NULL<-readRDS("BRCA/091514.BRCA.j.silhouette.ai.NULL.rds")
BRCA.S.j.ai.pval<-BRCA.S.j.ai.count[,
list(P.VAL=mean(as.vector(BRCA.S.j.NULL[SIZE==N.PATIENTS,]$RANDOM.CLUSTER.A.I)<CLUSTER.A.I)), by="KEGG_ID"]
BRCA.S.j.ai.pval$P.VAL.ADJ<-p.adjust(BRCA.S.j.ai.pval$P.VAL, method="fdr")
head(BRCA.S.j.ai.pval[order(P.VAL.ADJ),],10)

BRCA.S.p<-readRDS("BRCA/091514.BRCA.p.Silhouette.median.rds")
BRCA.S.p$S.I.C.TABLE[order(S.I.CLUSTER,decreasing=T),]
BRCA.S.p.ai.count<-BRCA.S.p$A.I.TABLE[,list(CLUSTER.A.I=median(A.I)), by="Hugo_Symbol"]
BRCA.S.p.ai.count<-as.data.table(merge(as.data.frame(BRCA.S.p.ai.count), as.data.frame(table.1.coverage)))
BRCA.S.p.NULL<-readRDS("BRCA/091514.BRCA.p.Silhouette.ai.NULL.rds")
BRCA.S.p.ai.pval<-BRCA.S.p.ai.count[,
list(P.VAL=mean(as.vector(BRCA.S.p.NULL[SIZE==N.PATIENTS,]$RANDOM.CLUSTER.A.I)<CLUSTER.A.I)), by="Hugo_Symbol"]
BRCA.S.p.ai.pval$P.VAL.ADJ<-p.adjust(BRCA.S.p.ai.pval$P.VAL, method="fdr")
head(BRCA.S.p.ai.pval[order(P.VAL.ADJ),],10)
BRCA.S.p.ai.pval[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

#Correlation cluster in matrix for production of only (Table 3)??
table.3<-readRDS("061214_Table.3.rds")
BRCA.CORR.KEGG<-cor(BRCA.CANCER.MATRICES.NORM$combined.matrices[
    intersect(rownames(BRCA.CANCER.MATRICES.NORM$combined.matrices),unique(as.vector(table.3$Enzyme))),
    BRCA.CANCER.MATRICES.NORM$cancer.patients], method="spearman")

pheatmap(BRCA.CORR.KEGG, scale="none")
heatmap.2(BRCA.CORR.KEGG, trace="none", scale="none",
    ColSideColors=c(rep(red,524), rep(blue,524)))

#Do driver genes cluster better cancer patients from normal tissue than passanger genes
BRCA.CORR.ALL<-cor(BRCA.CANCER.MATRICES.NORM$combined.matrices[,c(BRCA.CANCER.MATRICES.NORM$normal.patients,BRCA.CANCER.MATRICES.NORM$cancer.patients)], method="spearman")
length(BRCA.CANCER.MATRICES.NORM$normal)
colnames(BRCA.CORR.ALL)[1:112]
heatmap.2(BRCA.CORR.ALL, trace="none", scale="none", 
    ColSideColors=c(rep("green", length(BRCA.CANCER.MATRICES.NORM$normal)),
        rep("red", length(BRCA.CANCER.MATRICES$tumor))))

Function.Vis.heatmap.2.normal.cancer<-function(corr.matrix, table.1, normal.patients, target.gene){
    require(gplots)

    #Draws heatmap of normal vs cancer patients with gene of interest
    
    target.patients<-intersect(colnames(corr.matrix), as.vector(table.1[Hugo_Symbol==target.gene]$PATIENT))
    target.matrix<-corr.matrix[c(normal.patients, target.patients), c(normal.patients, target.patients)]

    heatmap.2(target.matrix, trace="none", scale="none",
        ColSideColors=c(rep("green", length(normal.patients)), 
            rep("red", length(target.patients))))
}
Function.Vis.heatmap.2.normal.cancer(BRCA.CORR.ALL, table.1, BRCA.CANCER.MATRICES.NORM$normal.patients, "AA06")

#####Analyze expression of genes known to be involved in cancer only and then cluster######
gene.path<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/NCI.PATHWAYS/091614.path.gene.rds")
gene.path
gene.path.death.apop<-gene.path[grepl(Path, ignore.case=T, pattern="apoptosis") | grepl(Path, ignore.case=T, pattern="death"),]

BRCA.CORR.DEATH.APOP<-cor(BRCA.CANCER.MATRICES.NORM$combined.matrices[
    intersect(as.vector(gene.path.death.apop$Hugo_Symbol), rownames(BRCA.CANCER.MATRICES.NORM$combined.matrices) ),BRCA.CANCER.MATRICES.NORM$cancer.patients], 
    method="spearman")

heatmap.2(BRCA.CORR.DEATH.APOP, trace="none", scale="none", 
    ColSideColors=Function.heatmap.2.mut.color(BRCA.CORR.DEATH.APOP, table.1, "TP53"))

#score.paths
gene.path.test<-gene.path[Path=="Regulation of ornithine decarboxylase (ODC)",]
BRCA.CORR.PATH.TEST<-cor(BRCA.CANCER.MATRICES.NORM$combined.matrices[
    intersect(as.vector(gene.path.test$Hugo_Symbol), rownames(BRCA.CANCER.MATRICES.NORM$combined.matrices) ),BRCA.CANCER.MATRICES.NORM$cancer.patients], 
    method="spearman")
dim(BRCA.CORR.PATH.TEST)
heatmap.2(BRCA.CORR.PATH.TEST, trace="none", scale="none", 
    ColSideColors=Function.heatmap.2.mut.color(BRCA.CORR.PATH.TEST, table.1, "TP53"))
sum(BRCA.CORR.PATH.TEST)

#hclust distances?
dim(BRCA.DISS)
BRCA.HCLUST<-hclust(as.dist(BRCA.DISS[1:10,1:10]), method="complete")
names(BRCA.HCLUST)
plot(BRCA.HCLUST)
BRCA.HCLUST$labels[BRCA.HCLUST$order]
as.data.table(cutree(BRCA.HCLUST, k=2), keep.rownames=T)


#Set up function to minimize pairwise distances over all 

#######DO RANDOM WALK FROM PATIENT THAT IS FARTHEST FROM NORMAL AND DO +/- DEPENDING ON CANCER OR NORMAL, THIS IS FOR EACH CLUSTER-GENE TO NORMAL, SCORE IS HIGHEST PEAK DIVIDED BY POSSIBLE HIGHEST PEAK - ANALOGOUS TO GSEA
