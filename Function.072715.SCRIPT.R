library(Hmisc)

Function.cor.matrix<-function(exp.matrix){
  #Assumes expression matrix has rows as genes and samples as columns
  require (Hmisc)
  
  #Filter correlation matrix for values that have 0 sd
  exp.matrix<-exp.matrix[apply(exp.matrix, 1, sd)!=0,]
  
  #Calculate correlation and p-values
  print ("calculating correlation")
  cor.matrix<-rcorr(t(data.matrix(exp.matrix)), type="spearman")
  
  #Correct p-values for multiple hypothesis testing
  print ("correcting for fdr")
  cor.matrix$P<-matrix(p.adjust(cor.matrix$P, method = "fdr"), ncol=ncol(cor.matrix$P), dimnames = list(rownames(cor.matrix$P), colnames(cor.matrix$P)))
  cor.matrix$P<- cor.matrix$P < 0.1
  diag(cor.matrix$P)<-FALSE
  
  #Use correlation value only if passes multiple hypothesis test
  cor.matrix<-cor.matrix$r * cor.matrix$P
  
  #Return
  return(cor.matrix)
}

#Arguments
args<-commandArgs(trailingOnly=T)
exp.obj<-readRDS(args[1])
output.file<-args[2]
print ("done loading files")

MAIN.OBJ<-Function.cor.matrix(exp.obj$tumor)

#Save to output
saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")