
Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  
  require(plyr)
  require(limma)
  require(edgeR)
  require(sva)
  
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
  
  #Quantile normalization
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data 
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Return topTable
  return(all.G.fit)
}

Function.RNAseq.Matrices.Normalization<-function(normal.matrix, cancer.matrix) {
  #Normalize RNAseq matrices for posterior differential expression analysis
  #May take a while depending on size of matrices
  
  require(plyr)
  require(limma)
  require(edgeR)
  require(sva)
  
  normal.matrix<-copy(cancer.matrices$normal)
  cancer.matrix<-copy(cancer.matrices$tumor)
  
  #Build design matrix
  G1.patients<-colnames(cancer.matrix)
  G0.patients<-colnames(normal.matrix)
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
  
  #Obtain batch normalization parameters
  batch.design=data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  batch.mod = model.matrix(~ G, batch.design) #What we have
  batch.mod0 = model.matrix(~ 1, batch.design) #Null
  batch.ss = svaseq(G.all$counts, batch.mod, batch.mod0) #Apply function to get batch effect parameters
  
  #Quantile normalization
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Obtain matrices corrected from batch effect
  f.ss = lmFit(G.all, model.matrix(~G + batch.ss$sv, batch.design)) #Just like regular fit but added batch effect normalization
  ss.removed = G.all$E - f.ss$coefficients[, -c(1, 2)] %*% t(f.ss$design[, -c(1, 2)]) #Get fully normalized matrices
  
  #Return combined matrix fully normalized (quantile, log cpm converated + batch effect) GENESxPATIENTS
  return(list(combined.matrices=ss.removed, cancer.patients=G1.patients, normal.patients=G0.patients))
}

batch.test<-Function.RNAseq.Matrices.Normalization(cancer.matrices$normal, cancer.matrices$tumor)
batch.test<-list(combined.matrices=ss.removed, cancer.patients=G1.patients, normal.patients=G0.patients)
batch.test$cancer.patients
saveRDS(batch.test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082614.CANCER.MATRICES.NORMALIZED.OBJ.rds")
saveRDS(Table.1.PLUS, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PLUS.rds")

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

test.diff.exp<-Function.RNAseq.Differential.Expression.V2(batch.test, as.vector(Table.1.PLUS[Hugo_Symbol=="TP53",]$PATIENT))
Table.1.PATIENT.COVERAGE
test.diff.exp
hist(test.diff.exp$adj.P.Val)
hist(abs(test.diff.exp$logFC))

Function.v.p.Version.3<-function(normalized.matrices.object, table1plus) {
  #Takes object from Function.RNAseq.Differential.Expression.V2 and Table.1.Plus to obtain v(p) object
  #This will produce a diff.expression list per gene - Not to be confused with single table object per v(p) as past versions
  #Requires:
  #   Function.RNAseq.Differential.Expression.V2()
  
  require(parallel)
  require(reshape2)
  require(data.table)
  require(limma)
  require(edgeR)
  
  #Filter table 1 plus for genes that have greater than 5 patients
  patient.coverage<-table1plus[,list(size=length(PATIENT)), by="Hugo_Symbol"]
  patient.coverage<-patient.coverage[size>5,]
  table1plus<-table1plus[Hugo_Symbol %in% as.vector(unique(patient.coverage$Hugo_Symbol)),]
  
  #Break table.1 into lists per gene
  table.1.lists<-split(table1plus, table1plus$Hugo_Symbol, drop=T)
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression.V2", "table.1.lists", "normalized.matrices.object"),envir=environment())
  
  #v(p) list object
  v.p.object<-parLapply(cl, names(table.1.lists), 
                        function(x) Function.RNAseq.Differential.Expression.V2(normalized.matrices.object, as.vector(table.1.lists[[x]]$PATIENT) ) )
  names(v.p.object)<-names(table.1.lists)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(v.p.object)  
}

v.p.object.test<-Function.v.p.Version.3(batch.test, Table.1.PLUS[Hugo_Symbol %in% unique(as.vector(Table.1.PLUS$Hugo_Symbol))[15001:20000],])
saveRDS(v.p.object.test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.15001.20000.rds")
Table.1.PLUS[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
length(unique(as.vector(Table.1.PLUS$Hugo_Symbol)))

#####082714#######
v.p.object.test<-list()
for (record in list.files("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS",pattern="082714")) {
  print (record)
  v.p.object.test<-c(v.p.object.test, readRDS(paste0("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/",record)) )
}
length(names(v.p.object.test))
saveRDS(v.p.object.test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")
rm(v.p.object.test)
Function.v.p.Version.3.NULL<-function(normalized.matrices.object, null.patient.count.vector) {
  #Obtain Null lists for 100x random sampling from cancer patients in normalized matrices object
  
  require(limma)
  require(data.table)
  require(data.table)
  require(reshape2)
  require(parallel)
  
  #Make sure patient null vector includes only samples greater than 5
  null.vector<-null.patient.count.vector[null.patient.count.vector>5]
  
  #All cancer patients to sample from
  cancer.patients<-normalized.matrices.object$cancer.patients
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression.V2", "null.vector", "normalized.matrices.object", "cancer.patients"),envir=environment())
  
  #Build NULL distributions
  Null.dist<-parLapply(cl, null.vector, function(x) replicate(100, Function.RNAseq.Differential.Expression.V2(normalized.matrices.object, sample(cancer.patients,x)),simplify=F))
  names(Null.dist)<-null.vector
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(Null.dist)
}

#Run above function in cluster
BRCA.Table.1.PATIENT.COVERAGE<-Table.1.PLUS[,list(PATIENT.COVERAGE=length(PATIENT)), by="Hugo_Symbol"]
saveRDS(BRCA.Table.1.PATIENT.COVERAGE, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PATIENT.COVERAGE.rds")
v.p.object.NULL<-Function.v.p.Version.3.NULL(batch.test, unique(as.vector(BRCA.Table.1.PATIENT.COVERAGE$PATIENT.COVERAGE))[1:5]) #Ran in cluster

#######082814#####
#Got v.p.object.test, which is a differential expression table for all v(p) genes
BRCA.Table.1.PATIENT.COVERAGE

v.p.table.mean.abs.logFC<-as.data.table(sapply(v.p.object.test, function(x) mean(abs(x$logFC))  ,USE.NAMES=T),keep.rownames=T)
setnames(v.p.table.mean.abs.logFC, c("Hugo_Symbol", "mean.abs.logFC"))

v.p.table.mean.abs.logFC<-as.data.table(merge(as.data.frame(v.p.table.mean.abs.logFC), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE), by="Hugo_Symbol"))
ggplot(v.p.table.mean.abs.logFC, aes(mean.abs.logFC, PATIENT.COVERAGE)) + geom_point() + theme.format

#Counting number of genes at different logFC thresholds
v.p.table.logFC.TH.count<-as.data.table(t(sapply(v.p.object.test, function(x) {
  y1.5<-sum(abs(x$logFC)>1.5)
  y2<-sum(abs(x$logFC)>2)
  y2.5<-sum(abs(x$logFC)>2.5)
  y3<-sum(abs(x$logFC)>3)
  y3.5<-sum(abs(x$logFC)>3.5)
  y4<-sum(abs(x$logFC)>4)
  return(c(y1.5,y2, y2.5, y3, y3.5, y4))
} , USE.NAMES=T)), keep.rownames=T) 

setnames(v.p.table.logFC.TH.count, c("Hugo_Symbol", "1.5", "2.0","2.5","3.0", "3.5","4.0"))
saveRDS(v.p.table.logFC.TH.count, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082814.BRCA.v.p.logFC.abs.TH.count")

v.p.table.logFC.TH.count<-melt(v.p.table.logFC.TH.count, id="Hugo_Symbol")
v.p.table.logFC.TH.count<-as.data.table(merge(as.data.frame(v.p.table.logFC.TH.count), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE),by="Hugo_Symbol"))
v.p.table.logFC.TH.count
ggplot(v.p.table.logFC.TH.count, aes(PATIENT.COVERAGE, value, colour=variable)) + geom_point() + theme.format
ggplot(v.p.table.logFC.TH.count, aes(as.factor(PATIENT.COVERAGE), value, colour=variable)) + geom_boxplot() + theme.format+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Need to look at NULL distribution
v.p.object.NULL<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v.p.object.NULL.rds")
v.p.object.NULL[["15"]][2]

v.p.table.null.logFC.TH.count<-lapply(names(v.p.object.NULL), function(x) {
  z=t(sapply(v.p.object.NULL[[x]], function(y) {
  y.all=abs(y$logFC)
  y.1.5<-sum(y.all>1.5)
  y.2<-sum(y.all>2)
  y.2.5<-sum(y.all>2.5)
  y.3<-sum(y.all>3)
  y.3.5<-sum(y.all>3.5)
  y.4<-sum(y.all>4)
  return(c(y.1.5, y.2, y.2.5, y.3, y.3.5, y.4))
}))
  z<-as.data.table(z)
  setnames(z, c("1.5", "2.0", "2.5", "3.0", "3.5", "4.0"))
  z$PATIENT.COVERAGE<-x
  z<-melt(z, id="PATIENT.COVERAGE")
  return(z)
})
names(v.p.table.null.logFC.TH.count)<-names(v.p.object.NULL)
v.p.table.null.logFC.TH.count.table<-do.call(rbind, v.p.table.null.logFC.TH.count)
v.p.table.null.logFC.TH.count.table$PATIENT.COVERAGE<-as.numeric(v.p.table.null.logFC.TH.count.table$PATIENT.COVERAGE)

ggplot(v.p.table.null.logFC.TH.count.table, aes(as.factor(PATIENT.COVERAGE), value, colour=variable)) + geom_boxplot() + theme.format +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Variance of Null distribution appears to decrease with increasing Sample coverage for all thersholds
ggplot(v.p.table.null.logFC.TH.count.table[, list(variance=var(value)), by=c("PATIENT.COVERAGE", "variable")], aes(PATIENT.COVERAGE, variance, colour=variable)) +
  geom_point()+ theme.format
ggplot(v.p.table.null.logFC.TH.count.table[, list(variance=var(value)), by=c("PATIENT.COVERAGE", "variable")], aes(PATIENT.COVERAGE, variance)) +
  geom_point()+ theme.format + facet_wrap(~variable,ncol=1,scales="free_y") + theme(strip.text.x=element_text(size=20))

######082914######
#Apply NULL logFC TH to obtained logFC counts per gene per TH 
v.p.table.logFC.TH.count
v.p.table.null.logFC.TH.count.table

ggplot(v.p.table.logFC.TH.count, aes(factor(PATIENT.COVERAGE), value)) + geom_point(shape=8) + geom_line(colour="red")+ theme.format + 
  geom_boxplot(data=v.p.table.null.logFC.TH.count.table, aes(factor(PATIENT.COVERAGE), value)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~variable,scales="free_y") + theme(strip.text.x=element_text(size=20)) +
  ylab("Number of Differentially Expressed Genes Beyond Threshold")

#Get p-value per gene by applying null across all thresholds
v.p.table.logFC.TH.count[as.numeric(variable)==2 & as.numeric(PATIENT.COVERAGE)==21,]
setnames(v.p.table.null.logFC.TH.count.table, c("NULL.PATIENT.COVERAGE", "null.variable", "null.value"))

v.p.table.logFC.TH.count.p.val<-v.p.table.logFC.TH.count[,list(p.val= mean(as.vector(v.p.table.null.logFC.TH.count.table[as.numeric(null.variable)==as.numeric(variable) &
                                                                              as.numeric(NULL.PATIENT.COVERAGE)==PATIENT.COVERAGE,]$null.value)>=value)),
                         by=c("Hugo_Symbol","variable")]

v.p.table.logFC.TH.count.p.val$p.val.adj<-p.adjust(v.p.table.logFC.TH.count.p.val$p.val, method="fdr")
v.p.table.logFC.TH.count.p.val[order(p.val.adj),]
hist(v.p.table.logFC.TH.count.p.val$p.val.adj)
ggplot(v.p.table.logFC.TH.count.p.val[Hugo_Symbol %in% COSMIC.BRCA$Symbol,], aes(Hugo_Symbol, p.val.adj, colour=variable)) + geom_bar(stat="identity", position="dodge") + 
  theme.format +theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_hline(yintercept=0.05, color="red")

#######083114##########
#Look at house-keeping genes logFC variation with respect to rest of genes
HK.genes<-as.data.table(read.csv("DATABASES/HK_genes.EISENBERG.txt", header=F, sep="\t", stringsAsFactors=F))
setnames(HK.genes, c("Hugo_Symbol", "Refseq"))
HK.genes$Hugo_Symbol<-sapply(HK.genes$Hugo_Symbol, function(x) strsplit(x," ")[[1]][1])

BRCA.diff.exp.table<-Function.RNAseq.Differential.Expression.V2(batch.test, batch.test$cancer.patients)
BRCA.diff.exp.table$House.keeping<-as.character(BRCA.diff.exp.table$ID) %in% as.character(as.vector(HK.genes$Hugo_Symbol))

ggplot(BRCA.diff.exp.table, aes(House.keeping, logFC)) + geom_boxplot() + theme.format
ggplot(BRCA.diff.exp.table, aes(logFC, colour=House.keeping)) + geom_histogram() + theme.format
ggplot(BRCA.diff.exp.table, aes(logFC)) + geom_histogram(, binwidth=0.01) + theme.format + facet_wrap(~House.keeping,ncol=1)
var.test(as.vector(BRCA.diff.exp.table[House.keeping==FALSE,]$logFC), as.vector(BRCA.diff.exp.table[House.keeping==TRUE,]$logFC)) 
#Look at logFC expressions at different thresholds of house-keeping genes

#Correlation within matrix of column patients for each gene v(p)
library(Hmisc)
a<-rcorr(as.matrix(mtcars), type="spearman")
a$r
batch.test
head(batch.test$combined.matrices[,1:4])
Table.1.PLUS

Function.v.p.corr.matrix.obj<-function(normalized.matrices.object, table1plus, hk.genes) {
  #Produces correlation matrix per causal gene in table1plus after removal of house-keeping genes (low noise) from expression matrix
  #Takes object from Function.RNAseq.Differential.Expression.V2 and Table.1.Plus to obtain v(p) correlation object
  #This will produce a correlation matrix per v(p) gene
  
  require(parallel)
  require(reshape2)
  require(data.table)
  require(Hmisc)
  
  #Filter table 1 plus for genes that have greater than 5 patients and are in our gene expression matrix only
  table1plus<-table1plus[PATIENT %in% normalized.matrices.object$cancer.patients,]
  patient.coverage<-table1plus[,list(size=length(PATIENT)), by="Hugo_Symbol"]
  patient.coverage<-patient.coverage[size>5,]
  table1plus<-table1plus[Hugo_Symbol %in% as.vector(unique(patient.coverage$Hugo_Symbol)),]
  
  #Break table.1 into lists per gene
  table.1.lists<-split(table1plus, table1plus$Hugo_Symbol, drop=T)
  
  #Get normalized matrix for cancer patient only
  cancer.normalized.matrix<-normalized.matrices.object$combined.matrices[, normalized.matrices.object$cancer.patients]
  
  #Remove house keeping genes from cancer expression matrix
  matrix.genes<-rownames(cancer.normalized.matrix)
  non.hk.genes<-matrix.genes[!(matrix.genes %in% hk.genes)]
  cancer.normalized.matrix<-cancer.normalized.matrix[non.hk.genes, ]
  print(dim(cancer.normalized.matrix))
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("table.1.lists", "cancer.normalized.matrix", "rcorr"),envir=environment())
  
  #Obtain correlation matrix per gene
  #cancer.corr.matrix.object<-parLapply(cl, table.1.lists, function(x) rcorr(cancer.normalized.matrix[, as.vector(x$PATIENT)] , type="spearman")$r )
  #names(cancer.corr.matrix.object)<-names(table.1.lists)
  cancer.corr.matrix.object<-parSapply(cl, table.1.lists, function(x) mean(apply(cancer.normalized.matrix[,as.vector(x$PATIENT)], 1, var)) ,USE.NAMES=T)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(cancer.corr.matrix.object)  
}

v.p.corr.matrix.obj<-Function.v.p.corr.matrix.obj(batch.test, 
                                                  Table.1.PLUS[Hugo_Symbol %in% c("TTN",
                                                                                  sample(unique(as.vector(Table.1.PLUS$Hugo_Symbol)),50),
                                                                                  as.vector(COSMIC.BRCA$Symbol)),], as.vector(HK.genes$Hugo_Symbol) )
sum(v.p.corr.matrix.obj$A2M)/ncol(v.p.corr.matrix.obj$A2M)^2
sum(v.p.corr.matrix.obj$TP53)/ncol(v.p.corr.matrix.obj$TP53)^2
sum(v.p.corr.matrix.obj$A1CF)/ncol(v.p.corr.matrix.obj$A1CF)^2
sum(v.p.corr.matrix.obj$TTN)/ncol(v.p.corr.matrix.obj$TTN)^2
as.vector(COSMIC.BRCA$Symbol)

as.data.table(v.p.corr.matrix.obj, keep.rownames=T)[rn %in% COSMIC.BRCA$Symbol,]
hist(as.data.table(v.p.corr.matrix.obj, keep.rownames=T)$v.p.corr.matrix.obj)

length(as.vector(Table.1.PLUS[Hugo_Symbol=="TP53",]$PATIENT))
length(as.vector(Table.1.PLUS[Hugo_Symbol=="TTN",]$PATIENT))
length(as.vector(Table.1.PLUS[Hugo_Symbol=="MYC",]$PATIENT))

length(intersect(as.vector(Table.1.PLUS[Hugo_Symbol=="PIK3CA",]$PATIENT), as.vector(Table.1.PLUS[Hugo_Symbol=="TTN",]$PATIENT)))

Table.1.PLUS[,list(SIZE=length(PATIENT)), by="Hugo_Symbol"][order(SIZE, decreasing=T),]

####Two questions
#   How likely is it to hit all genes found in the true targets? (Since it seems that at any point we have the capacity to hit the true number of targets anyways)
#   How likely is it that the patients uncovered by gene v(p) clustered together? - from Function.v.p.corr.matrix.obj()

BRCA.diff.exp.table
v.p.object.test<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")

Function.v.p.prediction<-function(v.p.object, cancer.diff.exp.table) {
  #Uses object from Function.v.p.Version.3() to obtain a coefficient per gene to determine how close it is to predicted diff gene for all thresholds in main cancer table
  
  require(data.table)
  require(reshape2)
  require(plyr)
  require(parallel)
  
  #Break cancer.diff.table into vectors of genes 
  diff.1.5<-as.vector(cancer.diff.exp.table[abs(logFC)>1.5,]$ID)
  diff.2.0<-as.vector(cancer.diff.exp.table[abs(logFC)>2.0,]$ID)
  diff.2.5<-as.vector(cancer.diff.exp.table[abs(logFC)>2.5,]$ID)
  diff.3.0<-as.vector(cancer.diff.exp.table[abs(logFC)>3.0,]$ID)
  diff.3.5<-as.vector(cancer.diff.exp.table[abs(logFC)>3.5,]$ID)
  diff.4.0<-as.vector(cancer.diff.exp.table[abs(logFC)>4.0,]$ID)
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("v.p.object", "diff.1.5", "diff.2.0", "diff.2.5", "diff.3.0", "diff.3.5", "diff.4.0"),envir=environment())
  
  #Run
  v.p.pred.table<-parSapply(cl, v.p.object, function(x) {
    y.1.5<-as.vector(x[abs(x$logFC)>1.5,]$ID)
    y.2.0<-as.vector(x[abs(x$logFC)>2.0,]$ID)
    y.2.5<-as.vector(x[abs(x$logFC)>2.5,]$ID)
    y.3.0<-as.vector(x[abs(x$logFC)>3.0,]$ID)
    y.3.5<-as.vector(x[abs(x$logFC)>3.5,]$ID)
    y.4.0<-as.vector(x[abs(x$logFC)>4.0,]$ID)
    
    #How many of the ones we predict are right (prediction)
    p.1.5<-length(intersect(y.1.5, diff.1.5))/length(y.1.5)
    p.2.0<-length(intersect(y.2.0, diff.2.0))/length(y.2.0)
    p.2.5<-length(intersect(y.2.5, diff.2.5))/length(y.2.5)
    p.3.0<-length(intersect(y.3.0, diff.3.0))/length(y.3.0)
    p.3.5<-length(intersect(y.3.5, diff.3.5))/length(y.3.5)
    p.4.0<-length(intersect(y.4.0, diff.4.0))/length(y.4.0)
    
    #How many of the target we get right (coverage)
    c.1.5<-length(intersect(y.1.5, diff.1.5))/length(diff.1.5)
    c.2.0<-length(intersect(y.2.0, diff.2.0))/length(diff.2.0)
    c.2.5<-length(intersect(y.2.5, diff.2.5))/length(diff.2.5)
    c.3.0<-length(intersect(y.3.0, diff.3.0))/length(diff.3.0)
    c.3.5<-length(intersect(y.3.5, diff.3.5))/length(diff.3.5)
    c.4.0<-length(intersect(y.4.0, diff.4.0))/length(diff.4.0)
    
    return(c(p.1.5, p.2.0, p.2.5, p.3.0, p.3.5, p.4.0, c.1.5, c.2.0, c.2.5, c.3.0,c.3.5,c.4.0))
    } ,USE.NAMES=T)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Clean up
  v.p.pred.table<-t(v.p.pred.table)
  v.p.pred.table<-as.data.table(v.p.pred.table, keep.rownames=T)
  setnames(v.p.pred.table, c("Hugo_Symbol", "p1.5","p2.0","p2.5","p3.0","p3.5","p4.0", "c1.5","c2.0","c2.5","c3.0","c3.5","c4.0"))
  
  #Return
  return(v.p.pred.table)
  
}
closeAllConnections()
v.p.prediction.table<-Function.v.p.prediction(v.p.object.test, BRCA.diff.exp.table)