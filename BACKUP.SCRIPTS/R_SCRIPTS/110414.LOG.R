####Dividing breast cancer patients into subtypes
library(data.table)
library(genefu)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

theme.format<-theme(axis.text.y=element_text(size=rel(2.5)), axis.text.x=element_text(size=rel(2.5)),
                    axis.title.y=element_text(size=22), axis.title.x=element_text(size=22),
                    legend.text = element_text(size = 22))

Function.RNAseq.Differential.Expression.V2<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  
  require(limma)
  require(edgeR)
  require(data.table)
  
  #Get target combined matrix
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  target.combined.matrix<-normalized.matrices.object$combined.matrices[,c(cancer.samples.in.matrix,normal.samples)]
  
  #Get design matrix
  G.design.matrix<-data.frame(G=c(rep("G1", length(cancer.samples.in.matrix)), rep("G0", length(normal.samples))))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Perform differential expression
  G.fit = lmFit(target.combined.matrix, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Convert to data.table
  if ("ID" %in% colnames(all.G.fit)) {
    all.G.fit<-as.data.table(all.G.fit)
  } else {
    all.G.fit$ID<-rownames(all.G.fit)
    all.G.fit<-as.data.table(all.G.fit)
  }
  
  #Return topTable
  return(all.G.fit)
}

BRCA.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/110414.CANCER.MATRICES.NORMALIZED.OBJ.NB.rds")
BRCA.EXP$combined.matrices

BRCA.DIFF<-Function.RNAseq.Differential.Expression.V2(BRCA.EXP, BRCA.EXP$cancer.patients)
saveRDS(BRCA.DIFF, file="PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/111114.BRCA.DIFF.EXP.rds")

#Using top 1000 differentially expressed genes on batch normalized exp matrix
pheatmap(BRCA.EXP$combined.matrices[as.vector(BRCA.DIFF$ID)[1:1000],BRCA.EXP$cancer.patients], scale="none")

#Using PAM50
data(pam50)
data(pam50.scale)
data(pam50.robust)
head(pam50$centroids)
head(pam50.robust$centroids)
pam50.robust$std

pam50$centroids.map
pheatmap(BRCA.EXP$combined.matrices[intersect(pam50$centroids.map$probe, rownames(BRCA.EXP$combined.matrices)),BRCA.EXP$cancer.patients], scale="none")
heatmap(BRCA.EXP$combined.matrices[intersect(pam50$centroids.map$probe, rownames(BRCA.EXP$combined.matrices)),BRCA.EXP$cancer.patients], scale="none")

#Using PAM50 centroids
Function.BRCA.SUBTYPE<-function(brca.normalized.obj, version=1){
  #Takes normalized expression object from breast cancer and returns subtype score and assigns subtype to each patient based on maximum score
  
  require(genefu)
  require(data.table)
  
  #Choose pam50 method
  if (version==1){
    data(pam50)
    pam<-copy(pam50)
  } else if(version==2){
    data(pam50.scale)
    pam<-copy(pam50.scale)
  } else if(version==3){
    data(pam50.robust)
    pam<-copy(pam50.robust)
  }
  
  #Obtain target pam gene in expression matrix
  target.genes<-intersect(rownames(pam$centroids), rownames(brca.normalized.obj$combined.matrices))
  
  #Simplify pam and expression matrices
  simplified.pam<-pam$centroids[target.genes,]
  simplified.exp<-brca.normalized.obj$combined.matrices[target.genes, brca.normalized.obj$cancer.patients]
  
  #Obtain predicted type based on correlation
  BRCA.SCORES<-cor(simplified.exp, simplified.pam, method="spearman")
  BRCA.SCORES<-as.data.frame(BRCA.SCORES)
  BRCA.SCORES$PATIENT<-rownames(BRCA.SCORES)
  BRCA.SCORES$TYPE<-colnames(BRCA.SCORES)[max.col(BRCA.SCORES[,1:5])]
  
  #Clean up and return
  BRCA.SCORES<-as.data.table(BRCA.SCORES)
  BRCA.SCORES<-BRCA.SCORES[order(TYPE),]
  return(BRCA.SCORES)
}

test<-Function.BRCA.SUBTYPE(BRCA.EXP,version=3)
saveRDS(test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/100414.BRCA.PATIENT.SUBTYPES.rds")
table(test$TYPE)

pheatmap(BRCA.EXP$combined.matrices[intersect(pam50$centroids.map$probe, rownames(BRCA.EXP$combined.matrices)),
                                    as.vector(test$PATIENT)], scale="none", 
         annotation=data.frame(TYPE=as.vector(test$TYPE), row.names=as.vector(test$PATIENT)),
         cluster_cols=FALSE, fontsize=16)

#Look at mutations patterns per molecular subtype
table.1.obj<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/100514.BRCA.Table.1.rds")
table.1<-merge(table.1.obj$table.1[Missense!=0,], test[,c("PATIENT", "TYPE"), with=F], by="PATIENT")
table.1.count<-table.1[,list(N.PATIENTS=length(Missense)), by=c("Hugo_Symbol", "TYPE")]
table.1.count<-merge(table.1.count, table.1[,list(TOTAL.TYPE.PATIENTS=length(unique(PATIENT))), by="TYPE"], by="TYPE")
table.1.count$PATIENT.PROPORTION<-table.1.count$N.PATIENTS/table.1.count$TOTAL.TYPE.PATIENTS

ggplot(unique(table.1.count[,c("TYPE", "TOTAL.TYPE.PATIENTS"), with=F]), aes(TYPE, TOTAL.TYPE.PATIENTS)) + geom_histogram(stat="identity") + theme.format
ggplot(table.1.count[Hugo_Symbol %in% c("TP53", "PIK3CA","BRCA1","TTN", "CDH1"),], aes(TYPE, PATIENT.PROPORTION, fill=Hugo_Symbol)) + 
  geom_histogram(stat="identity", position="dodge") + theme.format
ggplot(table.1.count[Hugo_Symbol %in% c("TP53", "PIK3CA","BRCA1","TTN", "CDH1"),], aes(TYPE, N.PATIENTS, fill=Hugo_Symbol)) + 
  geom_histogram(stat="identity", position="dodge") + theme.format

#Look at staging within each type
brca.clinical<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/110414.BRCA.CLINICAL.rds")
brca.clinical<-brca.clinical[TUMOR.STAGING!="Stage X"]
brca.clinical<-merge(test[,6:7,with=F], brca.clinical[,c("PATIENT", "TUMOR.STAGING"), with=F], by="PATIENT")
ggplot(brca.clinical, aes(TYPE, fill=TUMOR.STAGING)) + geom_histogram(position="dodge") + theme.format

#Look at mutation count per patient at different stages in each subtype
table.1.gene.count<-table.1[,list(N.MUTATIONS=sum(Missense)), by=c("PATIENT","TYPE")]
table.1.gene.count<-merge(table.1.gene.count, brca.clinical, by=c("PATIENT","TYPE"))

table.1.gene.count.wilcox<-table.1.gene.count[,{
  target.frame<-data.frame(N.MUT=N.MUTATIONS, STAGE=TUMOR.STAGING)
  early.mut<-as.vector(target.frame[TUMOR.STAGING=="Early",]$N.MUT)
  late.mut<-as.vector(target.frame[TUMOR.STAGING=="Late",]$N.MUT)
  
  w.test<-wilcox.test(early.mut, late.mut, paired=F)
  list(P.VAL=w.test$p.value, W=w.test$statistic)
}, by="TYPE"]
table.1.gene.count.wilcox

ggplot(table.1.gene.count, aes(TYPE, N.MUTATIONS, fill=TUMOR.STAGING)) + geom_boxplot() + theme.format + scale_y_log10() +
  annotation_custom(tableGrob(table.1.gene.count.wilcox[,1:2,with=F]),xmin=4.5, ymin=3)

#Perform glm(binomial) on Luminal A patients from wilcoxon output
brca.wilcoxon<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/110514.BRCA.Table.1.bmr.subtype.rds")
brca.wilcoxon$COVERAGE.TABLE
ggplot(melt(brca.wilcoxon$COVERAGE.TABLE[,c(1,2,4)]), aes(SUBTYPE, value, fill=variable)) + geom_histogram(stat="identity", position="dodge")+theme.format

brca.wilcoxon$SUBTYPES$LumA$wilcoxon.map

brca.model<-acast(brca.wilcoxon$SUBTYPES$LumA$wilcoxon.map[,c("Hugo_Symbol","Missense", "PATIENT"),with=F], PATIENT~Hugo_Symbol,fill=0, value.var="Missense")
brca.model<-as.data.frame(brca.model)
brca.model$PATIENT<-rownames(brca.model)
brca.model<-as.data.table(brca.model)
brca.model<-merge(brca.model, brca.clinical[,c(1,3),with=F, ], by="PATIENT")
brca.model$TUMOR.STAGING<-as.factor(brca.model$TUMOR.STAGING)

library(exactRankTests)
brca.model.binary<-brca.model[,2:(ncol(brca.model)-1),with=F]
brca.model.binary<-as.data.table(ifelse(brca.model.binary>=1,1,0))
brca.model.binary$PATIENT<-brca.model$PATIENT
brca.model.binary$TUMOR.STAGING<-brca.model$TUMOR.STAGING

brca.lumA.glm<-glm(TUMOR.STAGING~., data=brca.model.binary[,c(1:16,18),with=F], family=binomial())
brca.lumA.glm
summary(brca.lumA.glm)

#Perform decission tree
library(C50)
brca.lumA.dec.tree<-C5.0(brca.model.binary[,1:16, with=F], brca.model.binary$TUMOR.STAGING,10)
summary(brca.lumA.dec.tree)

#Is there anything enriched
brca.model<-brca.model[order(TUMOR.STAGING),]
brca.model.matrix<-as.matrix(brca.model[,c(2:17), with=F])
brca.model.matrix<-ifelse(brca.model.matrix>0, 1, 0) #Replace with 0-1
rownames(brca.model.matrix)<-brca.model$PATIENT
pheatmap(t(brca.model.matrix), scale="none",
         annotation=data.frame(TYPE=as.vector(brca.model$TUMOR.STAGING),row.names=brca.model$PATIENT),
         cluster_cols=FALSE, fontsize=16)