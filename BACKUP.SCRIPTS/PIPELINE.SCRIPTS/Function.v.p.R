#Function.v.p.R
#072014
#Calculates the differential expression influence of protein p
#It compares cancer(or disease) phenotype vs normal phenotype and gives a ratio of the number of differentially expressed genes
#Returns a table per protein with its respective v(j) value
#PREFERABLY TO BE USED IN A COMPUTATIONALLY INTENSIVE MACHINE WITH PARALLELIZATION

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
    print (length(patients))
    
    #Get expression matrix
    dummy.G1.expression<-rnaseq.matrices$tumor[,patients]
    
    #Do differential expression against normal
    dummy.diff.exp<-Function.RNAseq.Differential.Expression(rnaseq.matrices$normal, dummy.G1.expression)
    
    #Get v(p) by dividing that were differentially expressed over all genes tested
    dummy.G.diff.exp.genes<-nrow(dummy.diff.exp[dummy.diff.exp$adj.P.Val<0.05,]) #Value set at p<0.05 for bonferroni corrected p-values
    dummy.v.protein<-dummy.G.diff.exp.genes/nrow(dummy.diff.exp)
    
    #Return v(p)
    return(dummy.v.protein)
  }
  
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
  print (dim(cnv.table))
  cnv.table<-cnv.table[PATIENT %in% cancer.rnaseq.patients, ] #Patients in CNV data that have expression information
  print (dim(cnv.table))
  dummy.table.1<-unique(as.data.table(rbind(dummy.table.1, cnv.table)))
  print(dummy.table.1)
  
  #If a protein is not mutated in at least 2 patients remove from analysis - [Hugo_Symbol, PATIENT]
  target.genes<-as.data.table(table(as.vector(dummy.table.1$Hugo_Symbol))) #Number of patients per gene
  target.genes<-unique(as.vector(target.genes[N>1,]$V1)) #Genes that occur in at least two patients
  dummy.table.1<-dummy.table.1[Hugo_Symbol %in% target.genes,] 
  
  print(length(unique(as.vector(dummy.table.1$Hugo_Symbol))))
  #Split into list of patients per gene
  dummy.table.1.split<-split(dummy.table.1, dummy.table.1$Hugo_Symbol, drop=T)
  
  #START parallelization 
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("internal.function", "Function.RNAseq.Differential.Expression", "RNASEQ.MATRICES"),envir=environment())
  
  #Get differential expression versus normal 
  diff.table.v<-parSapply(cl, dummy.table.1.split, function(x) internal.function(as.vector(x$PATIENT), RNASEQ.MATRICES) ,USE.NAMES=T)
  diff.table.v<-as.data.table(diff.table.v, keep.rownames=T)
  setnames(diff.table.v, colnames(diff.table.v), c("Hugo_Symbol", "v.PROTEIN"))
  
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
