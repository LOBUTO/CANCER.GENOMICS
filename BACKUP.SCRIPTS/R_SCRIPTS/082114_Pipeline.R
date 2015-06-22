#####Function.v.p.NULL
######082114
#Takes Table.1 and CNV Table to get null v(p) of number of samples per gene p
#Sampling is done 100x
#Needs cancer matrices and v(p) table as well to provided corrected p.values
#cnv.table obtained from Function.process.GISTIC.TH.R

#####################################################################################################################################################################
Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  
  require(plyr)
  require(limma)
  require(edgeR)
  
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 20 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data 
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Return topTable
  return(all.G.fit)
}

Function.v.p.NULL<-function(table.1.rds, table.cnv.rds, cancer.matrices.rds, v.p.table.object.rds) {
  
  require(parallel)
  require(reshape2)
  require(data.table)
  require(base)
  require(plyr)
  require(limma)
  require(edgeR)
  
  internal.function.100<-function(cancer.matrices,n){
    CANCER.PATIENTS<-colnames(cancer.matrices$tumor)
    REPLICATIONS<-replicate(100, Function.RNAseq.Differential.Expression(cancer.matrices$normal, cancer.matrices$tumor[, sample(CANCER.PATIENTS,n)]) ,simplify=F)
    DISTRIBUTION<-sapply(REPLICATIONS, function(x) nrow(x[x$adj.P.Val<0.05,])/nrow(x))
    return(DISTRIBUTION)
  }
  
  #Load tables and matrices
  #table.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds")
  table.1<-readRDS(table.1.rds)
  table.1<-table.1$table.1
  
  cnv.table<-readRDS(table.cnv.rds)
  
  cancer.matrices<-readRDS(cancer.matrices.rds)
  
  v.p.table.object<-readRDS(v.p.table.object.rds)
  v.p.table<-v.p.table.object$v.table[, c("Hugo_Symbol", "v.PROTEIN"), with =F]
  
  table.1<-copy(BRCA.Table.1$table.1)
  cnv.table<-copy(BRCA.CNV.Table.1)
  cancer.matrices<-copy(BRCA.CANCER.MATRICES)
  v.p.table<-copy(BRCA.v.test)
  
  #Make master Table.1
  Table.1.PLUS<-table.1
  Table.1.PLUS$Tumor_Sample_Barcode<-as.character(Table.1.PLUS$Tumor_Sample_Barcode)
  Table.1.PLUS$PATIENT<-sapply(Table.1.PLUS$Tumor_Sample_Barcode, function(x) paste0(strsplit(x, "-")[[1]][1:4], collapse="."))
  Table.1.PLUS$Tumor_Sample_Barcode<-NULL
  Table.1.PLUS<-unique(rbind(Table.1.PLUS[,c(1,3), with=F], cnv.table))
  
  #Get Patient coverage per gene to obtain "n" counts per gene for random sampling
  Table.1.PATIENT.COVERAGE<-Table.1.PLUS[,list(PATIENT.COVERAGE=length(PATIENT)), by="Hugo_Symbol"]
  Table.1.PATIENT.COVERAGE<-Table.1.PATIENT.COVERAGE[PATIENT.COVERAGE>1,]
  GENE.COVERAGE.LIST<-unique(as.vector(Table.1.PATIENT.COVERAGE$PATIENT.COVERAGE)) #Vector of NULL distributions to build
  
  #Set up Parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "cancer.matrices", "internal.function.100"),envir=environment())  
  
  #Run random sampling for NULLs
  RAN.NULL<-parLapply(cl, GENE.COVERAGE.LIST, function(x) internal.function.100(cancer.matrices, x))
  names(RAN.NULL)<-GENE.COVERAGE.LIST
  stopCluster(cl)
  
  RAN.NULL.50<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082314.BRCA.v.NULL.50")
  RAN.NULL.100<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082314.BRCA.v.NULL.100")
  RAN.NULL.150<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082314.BRCA.v.NULL.150")
  RAN.NULL.220<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082314.BRCA.v.NULL.220")
  RAN.NULL<-c(RAN.NULL.50,RAN.NULL.100, RAN.NULL.150,RAN.NULL.220)
  RAN.NULL$`216`
  RAN.NULL[["216"]]
  
  #Get p-values per gene
  v.p.table<-as.data.table(merge(as.data.frame(v.p.table), as.data.frame(Table.1.PATIENT.COVERAGE), by="Hugo_Symbol"))
  v.p.table<-v.p.table[,list(P.VAL=sum(RAN.NULL[[as.character(PATIENT.COVERAGE)]]>=v.PROTEIN)/ length(RAN.NULL[[as.character(PATIENT.COVERAGE)]])),
                       by=c("Hugo_Symbol", "v.PROTEIN", "SAMPLE.COVERAGE")]
  v.p.table$P.VAL.ADJ<-p.adjust(v.p.table$P.VAL, method="fdr")
  hist(v.p.table$P.VAL.ADJ)
  
  #Return
  return(v.p.table)
  
}
test<-copy(v.p.table)
test$SIG<-test$P.VAL.ADJ<0.05
ggplot(test, aes(v.PROTEIN, SAMPLE.COVERAGE, colour=SIG)) + geom_point() + theme.format
ggplot(test, aes(x=SIG, y=SAMPLE.COVERAGE)) + geom_boxplot() +theme.format 
wilcox.test(as.vector(test[SIG==TRUE,]$SAMPLE.COVERAGE) , as.vector(test[SIG==FALSE,]$SAMPLE.COVERAGE)  )

v.p.table[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

#Try with re-sampling instead
dim(cancer.matrices$normal)

AKT1.BOOT<-boot(as.vector(Table.1.PLUS[Hugo_Symbol=="AKT1",]$PATIENT),function(x,d){
  DIFF=Function.RNAseq.Differential.Expression(cancer.matrices$normal, cancer.matrices$tumor[,x[d]])
  v.P= nrow(DIFF[DIFF$adj.P.Val<0.05,])/nrow(DIFF)
  return(v.P)  
} ,100)
AKT1.BOOT
quantile(AKT1.BOOT$t,c(0.025,0.975))
plot(AKT1.BOOT)

TP53.BOOT<-boot(intersect(as.vector(Table.1.PLUS[Hugo_Symbol=="TP53",]$PATIENT), colnames(cancer.matrices$tumor)) ,function(x,d){
  DIFF=Function.RNAseq.Differential.Expression(cancer.matrices$normal, cancer.matrices$tumor[,x[d]])
  v.P= nrow(DIFF[DIFF$adj.P.Val<0.05,])/nrow(DIFF)
  return(v.P)  
} ,100)
TP53.BOOT
quantile(TP53.BOOT$t,c(0.025,0.975))
plot(TP53.BOOT)

A1BG.BOOT<-boot(intersect(as.vector(Table.1.PLUS[Hugo_Symbol=="A1BG",]$PATIENT), colnames(cancer.matrices$tumor)) ,function(x,d){
  DIFF=Function.RNAseq.Differential.Expression(cancer.matrices$normal, cancer.matrices$tumor[,x[d]])
  v.P= nrow(DIFF[DIFF$adj.P.Val<0.05,])/nrow(DIFF)
  return(v.P)  
} ,100)
A1BG.BOOT
quantile(A1BG.BOOT$t,c(0.025,0.975))
plot(A1BG.BOOT)

library(boot)
test.boot<-boot(1:10,function(x,d){
  print (d)
  return(mean(x[d]))  } ,R=100)
test.boot
mean(test.boot$t)
test.boot$t0
test.boot$weights
plot(test.boot)

sd(test.boot$t)/sqrt(100)

sample(1:10, replace=T)

#####BOOTSRAP MAY NOT BE USEFUL

######082514#######
#####Classify Normals in BRCA########
head(cancer.matrices$normal[,1:3])
dummy.matrix<-DGEList(counts=cancer.matrices$normal)
dummy.matrix<-calcNormFactors(dummy.matrix)
dummy.matrix<-voom(dummy.matrix, normalize.method="quantile")
pheatmap(dummy.matrix$E, scale="none")

head(dummy.matrix$E[,1:3])

#Cluster normals
BRCA.NORMAL.CLUSTERS<-kmeans(t(dummy.matrix$E), 2,nstart=100)
head(BRCA.NORMAL.CLUSTERS)
BRCA.NORMAL.CLUSTERS$size

BRCA.NORMAL.CLUSTERS<-as.data.table(BRCA.NORMAL.CLUSTERS$cluster,keep.rownames=T)

CLINICAL.RACE<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/99d90670-467c-4c6b-b870-027712c19014/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt", header=T, sep="\t", stringsAsFactors=F,skip=2))
CLINICAL.RACE<-CLINICAL.RACE[,c(1,7,9,10), with=F]
CLINICAL.RACE$