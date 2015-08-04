library(data.table)
library(reshape2)

Function.tang.cor.matrix<-function(exp.matrix, tang.clean){
  
  require (Hmisc)
  
  #Load tang matrix
  tang<-read.csv(tang.clean, header=F, stringsAsFactors=F)
  setnames(tang, as.vector(as.matrix(tang[1,])) )
  tang<-tang[2:nrow(tang),]
  rownames(tang)<-tang$METABOLITE
  tang$METABOLITE<-NULL
  colnames(tang)<-sapply(colnames(tang), function(x) {
    if (x=="Normal Breast"){
      cn<-"NORMAL"
    } else {
      cn<-paste0(c("TCGA",strsplit(x, "-")[[1]], "01A"),collapse = ".")
    }
    return(cn)
  })
  
  #Construct tang expression matrix 
  exp.matrix<-exp.matrix[,intersect(colnames(tang), colnames(exp.matrix))]
  
  #Filter expression matrix for values that have 0 sd
  exp.matrix<-exp.matrix[apply(exp.matrix, 1, sd)!=0,]
  print ("Done cleaning up expression matrix")
  print (dim(exp.matrix))
  
  #Calculate correlation and p-values
  cor.matrix<-rcorr(t(data.matrix(exp.matrix)), type="spearman")
  print ("Done with correlation matrix")
  
  #Correct p-values for multiple hypothesis testing
  print ("Correcting for multiple hypothesis testing")
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
tang.clean<-args[2]

output.file<-args[3]
print ("done loading files")

MAIN.OBJ<-Function.tang.cor.matrix(exp.obj$tumor, tang.clean)

#Save to output
saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")
