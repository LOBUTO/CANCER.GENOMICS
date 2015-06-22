library(parallel)
library(data.table)
library(reshape2)

Function.Inter.MH.Score<-function(cleaned.recon, cleaned.cor, beta){
  
  #Filter recon table and correlation matrices for common genes
  common.genes<-intersect(cleaned.recon$Hugo_Symbol, intersect(rownames(cleaned.cor$COR.NORMAL), rownames(cleaned.cor$COR.CANCER)))
  recon.table<-cleaned.recon[Hugo_Symbol %in% common.genes,]
  cor.normal<-cleaned.cor$COR.NORMAL[common.genes, common.genes]
  cor.tumor<-cleaned.cor$COR.CANCER[common.genes, common.genes]
  
  #Obtain recon counts
  recon.count<-recon.table[,list(N.GENES=length(Hugo_Symbol)), by ="ID"]
  recon.count<-recon.count[N.GENES>=3, ]
  recon.table<-recon.table[ID %in% recon.count$ID,]
  recon.count<-unique(recon.count$N.GENES)
  print (recon.count)
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "cor.normal", "cor.tumor",
                              "recon.count", "recon.table", "beta") ,envir=environment())
  
  #Obtain Inter-dysregulation score per metabolite
  print ("Obtaining permutations per count")
  
  recon.genes<-unique(recon.table$Hugo_Symbol)
  main.list<-lapply(recon.count, function(x) {
    print (x)
    
    #Sample recon.count n genes 100 times from recon.gene pool
    count.samples<-replicate(100, sample(recon.genes, x), simplify = F)
    
    #Obtain count scores 
    count.scores<-parSapply(cl, count.samples, function(y) {
      
      Hugo.target<-unique(y)
      other.hugos<-setdiff(recon.genes, Hugo.target)
      
      dys.diff<-abs(cor.tumor[Hugo.target, other.hugos] - cor.normal[Hugo.target, other.hugos])
      n.score<-apply(dys.diff, 1, median)
      n.score<-mean(n.score>beta)
      
      return(n.score)
    })
    
    #Return count iteration scores per count
    return(count.scores)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Clean up
  names(main.list)<-as.character(recon.count)
  
  #Return
  return(main.list)
}

#Arguments
args<-commandArgs(trailingOnly=T)
cleaned.recon<-readRDS(args[1])
cleaned.cor<-readRDS(args[2])
beta<-as.numeric(args[3])
output.file<-args[4]
print ("done loading files")

MAIN.OBJ<-Function.Inter.MH.Score(cleaned.recon, cleaned.cor, beta=beta)

saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")