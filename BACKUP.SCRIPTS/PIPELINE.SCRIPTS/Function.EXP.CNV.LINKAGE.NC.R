#Function.EXP.CNV.LINKAGE.NC.R
#021715
#Calculates the linkage between mutated gene (regardless of position) and expression between normal and cancer populations
#Will only take into account sites that are seen mutated at least 2 times in population

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.EXP<-function(exp.rds, paired=T){
  
  #Load expression matrix file
  exp.matrix<-readRDS(exp.rds)
  
  #Work with paired samples only
  if (paired==T){
    normal.patients<-substr(exp.matrix$normal.patients, 1, 12)
    
    #Get cancer patients with a matched normal
    cancer.patients<-data.table(samples=exp.matrix$cancer.patients, patients=substr(exp.matrix$cancer.patients, 1, 12))
    cancer.patients$MATCH<-cancer.patients$patients %in% normal.patients
    cancer.patients<-cancer.patients[MATCH==TRUE, ]$samples
    
    #Filter expression matrix for it
    exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, cancer.patients)]
    exp.matrix$cancer.patients<-cancer.patients
    
    #Return filtered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)
    
  } else {
    #Return unfiltered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)  
  }
  
}

Function.Main<-function(cancer.cnv, exp.matrix){
  
  #Load CNV
  cnv<-readRDS(cancer.cnv)
  cnv<-cnv$combined.matrices
  
  #Get set of total patients with expression and cnv data
  patients<-intersect(colnames(cnv), colnames(exp.matrix$combined.matrices))
  
  #Filter both tables for these patients
  cnv<-cnv[,patients]
  exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,patients]
  print(dim(cnv))
  print(dim(exp.matrix$combined.matrices))
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("cnv", "as.data.table","exp.matrix","data.table","patients") ,envir=environment())
  print ("Done exporting values")
  
  #Execute parallelization
  print ("Finding linkage")
  main.list<-parLapply(cl, rownames(cnv), function(x) {
    
    mut.gene<-x
    cnv.vector<-as.vector(cnv[x,patients])
    
    cnv.list<-apply(exp.matrix$combined.matrices, 1, function(y) {
      spearman<-cor.test(cnv.vector, y[patients])
      s.table<-data.table(RHO=spearman$estimate, P.VAL=spearman$p.value)
      return(s.table)
    })
    
    cnv.table<-do.call(rbind, cnv.list)
    cnv.table$Hugo_EXP=names(cnv.list)
    cnv.table$Hugo_MUT=mut.gene
    cnv.table$N.PATIENTS=length(patients)
    
    return(cnv.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with linkages")
  
  #Merge gene tables
  main.table<-do.call(rbind, main.list)
  
  #Correct for multiple hypothesis
  print ("Correcting for multiple hypothesis testing")
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Cleaup and return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cancer.cnv<-args[1] #Of the form 021815.BRCA.PROCESSED.CNV.rds
exp.rds<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
exp.matrix<-Function.Prep.EXP(exp.rds, paired=F)
print ("done prepping expression rds")

main.function<-Function.Main(cancer.cnv, exp.matrix)
print ("Done with functions")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)
print ("Done writing file")