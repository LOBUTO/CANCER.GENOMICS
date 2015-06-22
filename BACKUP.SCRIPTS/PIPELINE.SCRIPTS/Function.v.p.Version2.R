#Function.v.p.R - VERSION 2 
#081014
#Calculates the differential expression influence of protein p
#It compares cancer(or disease) phenotype vs normal phenotype and returns values based on all, pos or neg FC values for differentially expressed genes
#It returns:
# Hugo_Symbol, v.PROTEIN (pos and neg), v.PROTEIN.POS (pos only), MeanLogFC(pos and neg), SUM.LogFC(pos+neg), POS.SUM.LogFC(pos only), SAMPLE.COVERAGE
#PREFERABLY TO BE USED IN A COMPUTATIONALLY INTENSIVE MACHINE (32+ cores) WITH PARALLELIZATION!!!

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

Function.v.p<-function(cancer.matrices, processed.table.1.rds, processed.th.cnv) {
  #Calculates v(j) based on degree of influence of protein on cancer phenotype when compared to normal
  #Need functions pre-loaded:
  #   Function.read.RNAseq.files()
  #   Function.process.RNAseq.map.files()
  
  require(base)
  require(data.table)
  require(parallel)
  internal.function<-function(patients, rnaseq.matrices){
    
    #Get expression matrix
    dummy.G1.expression<-rnaseq.matrices$tumor[,patients]
    
    #Do differential expression against normal
    dummy.diff.exp<-Function.RNAseq.Differential.Expression(rnaseq.matrices$normal, dummy.G1.expression)
    print (dummy.diff.exp)
    
    #Get v(p) by dividing that were differentially expressed over all genes tested
    dummy.G.diff.exp.genes<-dummy.diff.exp[dummy.diff.exp$adj.P.Val<0.05,] #Value set at p<0.05 for bonferroni corrected p-values
    dummy.v.protein<-nrow(dummy.G.diff.exp.genes)/nrow(dummy.diff.exp)
    dummy.v.protein.pos<-nrow(dummy.G.diff.exp.genes[dummy.G.diff.exp.genes$logFC>0,])/nrow(dummy.diff.exp)
    
    if ("ID" %in% colnames(dummy.G.diff.exp.genes)){
      dummy.v.proteins.neg<-as.vector(dummy.G.diff.exp.genes[dummy.G.diff.exp.genes$logFC<0,]$ID) #Negative logFC diff exp genes
      dummy.v.proteins<-as.vector(dummy.G.diff.exp.genes$ID) #All logFC diff exp genes
    } else {
      dummy.v.proteins.neg<-as.vector(rownames(dummy.G.diff.exp.genes[dummy.G.diff.exp.genes$logFC<0,]))
      dummy.v.proteins<-as.vector(rownames(dummy.G.diff.exp.genes))
    }
    
    #Get mean absolute logFC, sum and POS.SUM
    dummy.v.logFC<-mean(abs(as.vector(dummy.G.diff.exp.genes$logFC)))
    dummy.v.SUM.logFC<-sum(as.vector(dummy.G.diff.exp.genes$logFC))
    dummy.v.POS.SUM.logFC<-sum(as.vector(dummy.G.diff.exp.genes[dummy.G.diff.exp.genes$logFC>0,]$logFC))
    
    #Return v(p)
    return(list(c(dummy.v.protein, dummy.v.protein.pos, dummy.v.logFC, dummy.v.SUM.logFC, dummy.v.POS.SUM.logFC), dummy.v.proteins, dummy.v.proteins.neg))
  }
  
  #####CALL FILES##### - Only used when called from within function
  #processed.table.1.rds<-....
  cancer.matrices<-"PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES"
  processed.th.cnv<-"PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.GISTIC.TH.2.rds"
  
  ###################
  
  #Load RDS files 
  dummy.table.1<-readRDS(processed.table.1.rds)
  dummy.table.1<-dummy.table.1$table.1[,1:2, with=F] #[Tumor_Sample_Barcode, Hugo_Symbol]
  
  #Create tumor and cancer matrices using map file on RNA-seq folder
  RNASEQ.MATRICES<-readRDS(cancer.matrices)
  
  #Get cancer samples in processed.table.1 that have expression information in cancer matrix [Hugo_Symbol, PATIENT]
  cancer.rnaseq.patients<-colnames(RNASEQ.MATRICES$tumor)
  dummy.table.1$PATIENT<-sapply(as.character(dummy.table.1$Tumor_Sample_Barcode), function(x) paste(strsplit(x, "-")[[1]][1:4] , collapse="."))
  dummy.table.1<-dummy.table.1[PATIENT %in% cancer.rnaseq.patients, ] #Patients in table.1 that have gene expresion information
  dummy.table.1$Tumor_Sample_Barcode<-NULL
  
  #Integrate CNV data and keep only those that have expression information
  cnv.table<-readRDS(processed.th.cnv)
  cnv.table<-cnv.table[PATIENT %in% cancer.rnaseq.patients, ] #Patients in CNV data that have expression information
  dummy.table.1<-unique(as.data.table(rbind(dummy.table.1, cnv.table)))
  
  #If a protein is not mutated in at least 10 patients remove from analysis - [Hugo_Symbol, PATIENT] - THRESHOLD TO 10 PATIENTS
  target.genes<-as.data.table(table(as.vector(dummy.table.1$Hugo_Symbol))) #Number of patients per gene
  target.genes<-unique(as.vector(target.genes[N>10,]$V1)) #Genes that occur in at least ten patients
  dummy.table.1<-dummy.table.1[Hugo_Symbol %in% target.genes,] 
  
  #Split into list of patients per gene
  dummy.table.1.split<-split(dummy.table.1, dummy.table.1$Hugo_Symbol, drop=T)
  
  #START parallelization 
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("internal.function", "Function.RNAseq.Differential.Expression", "RNASEQ.MATRICES"),envir=environment())
  clusterEvalQ(cl, library(limma))
  clusterEvalQ(cl, library(edgeR)) 
  clusterEvalQ(cl, library(plyr))  
  
  #Get differential expression versus normal
  diff.table.v<-parLapply(cl, dummy.table.1.split, function(x) internal.function(as.vector(x$PATIENT), RNASEQ.MATRICES))
  
  #Store values in 2 datasets
  diff.table.stat<-as.data.table(do.call(rbind, lapply(diff.table.v, function(x) x[[1]])), keep.rownames=T)
  setnames(diff.table.stat, colnames(diff.table.stat), c("Hugo_Symbol", "v.PROTEIN", "v.PROTEIN.POS","Mean.LogFC", "SUM.LogFC", "POS.SUM.LogFC"))
  diff.list.genes<-lapply(diff.table.v, function(x) x[[2]])
  diff.list.neg.genes<-lapply(diff.table.v, function(x) x[[3]])
  
  #Add Sample Column for later normalization by content
  N.SAMPLES<-length(unique(as.vector(dummy.table.1$PATIENT)))
  GENE.COVERAGE<-dummy.table.1[,list(SAMPLE.COVERAGE=length(PATIENT)/N.SAMPLES), by="Hugo_Symbol"]
  diff.table<-merge(diff.table.stat, GENE.COVERAGE, by="Hugo_Symbol")
  
  #Combine into list for both type of data
  dummy.table.v<-list(v.table=diff.table, diff.neg.list=diff.list.neg.genes, diff.all.list=diff.list.genes)
  
  #STOP parallelization
  stopCluster(cl)
  
  #Return v(p) table
  return(dummy.table.v)
}

#PROCESS ENTRIES
library(base)

args<-commandArgs(trailingOnly=T)

input.cancer.matrices<-args[1]
input.table.1.rds<-args[2]
input.cnv.rds<-args[3]
output.file<-args[4]

#EXECUTE
Run<-Function.v.p(input.cancer.matrices, input.table.1.rds, input.cnv.rds)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)
