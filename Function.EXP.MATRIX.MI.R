#######Function.EXP.MATRIX.MI.R######
#030315
#Calculate the pairwise degree of mutual information between genes in an expression matrix

##################FUNCTIONS################
library(data.table)

Function.PREP.EXP.MATRIX<-function(exp.rds, patients="NORMAL"){
  
  #Load matrix obj
  exp.obj<-readRDS(exp.rds)
  
  #Filter for desired subset
  if (patients=="NORMAL"){
    exp.matrix<-exp.obj$combined.matrices[,exp.obj$normal.patients]
  } else if (patients=="CANCER") {
    exp.matrix<-exp.obj$combined.matrices[,exp.obj$cancer.patients]
  }
  
  #Return
  return(exp.matrix)
}

Function.GENE.MATRIX.MI<-function(exp.matrix){
  
  require(entropy)
  require(parallel)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix", "discretize2d", "mi.empirical") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  genes<-rownames(exp.matrix)
  count<-0
  main.list<-lapply(genes, function(x) {
    
    #Get expression for looping gene
    gene.exp<-as.vector(exp.matrix[x,,drop=F]) 
    
    #Apply over all others pairwise
    mi.gene.vector<-parApply(cl, exp.matrix[,,drop=F], 1, function(y) {
      
      #Discretize pairs
      gene.disc<-discretize2d(gene.exp, as.vector(y), numBins1=5, numBins2=5)
      mi.pair<-mi.empirical(gene.disc)
      
      #Return pair mi
      return (mi.pair)
      
    })
    
    #Construct table for gene table
    gene.table<-data.table(Hugo.1=x, Hugo.2=rownames(exp.matrix[,,drop=F]), MI=mi.gene.vector)
    
    #Remove gene from matrix 
    exp.matrix<<-exp.matrix[setdiff(rownames(exp.matrix[,,drop=F]), x),]
    
    #Update count - Don't need to know matrix dimensions, it throws an error when it reaches zero, removed
    count<<-count+1
    print (count/length(genes))
    
    #Return mi table
    return(gene.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")  
  
  #Combine mi tables
  main.table<-do.call(rbind, main.list)
  
  #Clean up and return
  main.table<-main.table[Hugo.1!=Hugo.2, ]
  main.table<-main.table[order(MI, decreasing=T),]
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
exp.rds<-args[1]
output.file<-args[2]
print("opened files")
##########################################

##################EXECUTE#################
exp.matrix<-Function.PREP.EXP.MATRIX(exp.rds, patients="NORMAL")
print ("done prepping expression rds")

main.function<-Function.GENE.MATRIX.MI(exp.matrix)
print ("Done with main function")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)
print ("Done writing to file")