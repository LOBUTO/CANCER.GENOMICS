#Function.MET.HUB.INTER.NULL.R
#052215
#Function to build distribution for MET.HUB.AGG.SCORE from Function.Teru.Met.Hub.Met.Impact() function

library(data.table)
library(reshape2)
library(parallel)

#Functions
Function.Main<-function(teru.exp.obj, teru.met.hub.obj){
  
  #Load met.hub.obj
  teru.met.hub<-readRDS(teru.met.hub.obj)
  
  #Load gene expression and scale expression to z-scores
  teru.exp<-readRDS(teru.exp.obj)
  teru.exp$MATRIX<-scale(teru.exp$MATRIX-apply(teru.exp$MATRIX,1,median))
  
  teru.normal<-teru.exp$CLASS[CLASS=="Normal",]$SAMPLE
  teru.cancer<-teru.exp$CLASS[CLASS=="Tumor",]$SAMPLE
  
  #Filter hub genes for those present in expression matrix
  teru.met.hub<-teru.met.hub[Hugo_Symbol %in% rownames(teru.exp$MATRIX),]
  
  #Get gene count to obtain null distribution for and gene pool from available metabolic genes
  gene.counts<-teru.met.hub[,list(N.GENES=length(Hugo_Symbol)), by="MET"]
  gene.counts<-unique(gene.counts$N.GENES)
  gene.pool<-unique(teru.met.hub$Hugo_Symbol)
  
  #Obtain hub dysregulation 100x per random sample of hub size
  mets<-unique(teru.met.hub$MET)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "mets", "teru.met.hub", "teru.exp", 
                              "teru.normal", "teru.cancer", "gene.counts", "gene.pool") ,envir=environment())
  print ("Done exporting values")
  
  #Loop through random sizes
  main.list<-lapply(gene.counts, function(x) {
    print (x)
    
    #Obtain 100 random hubs of size x
    random.hubs<-replicate(100, sample(gene.pool, x), simplify = F)
    
    #Calculate aggregated inter score per each random hub 
    size.random.list<-parSapply(cl, random.hubs, function(y) {
     
      random.hub.genes<-y
      
      #For random hub genes, calculate mean absolute pair correlation difference (MAPCD)
      hub.inter.MAPCD<-sapply(mets, function(z){
        
        non.hub.genes<-unique(teru.met.hub[MET==z,]$Hugo_Symbol)
        
        #normal cor
        normal.cor<-c(cor(t(teru.exp$MATRIX[random.hub.genes, teru.normal, drop=F]), t(teru.exp$MATRIX[non.hub.genes, teru.normal, drop=F]),method = "spearman"))
        
        #cancer cor
        cancer.cor<-c(cor(t(teru.exp$MATRIX[random.hub.genes, teru.cancer, drop=F]), t(teru.exp$MATRIX[non.hub.genes, teru.cancer, drop=F]),method = "spearman"))
        
        #Calculate mean absolute pair correlation difference (MAPCD)
        MAPCD<-mean(abs(normal.cor- cancer.cor))
        
        #Return
        return(MAPCD)
      })
      
      #Obtain aggregate MAPCD score for random iteration of size x
      AMAPCD=mean(hub.inter.MAPCD)
      
      #Return
      return(AMAPCD)
    })
    
    #Construct null table of AMPACD for size x
    random.table<-data.table(NULL.SIZE=x, AMAPCD=size.random.list)
    
    #Return
    return(random.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Reconstruct main table for all null distributions
  main.table<-do.call(rbind, main.list)
  
  #Return
  return(main.table)
}


#Arguments
args<-commandArgs(trailingOnly=T)
teru.exp.obj<-args[1]
teru.met.hub.obj<-args[2]
output.file<-args[3]
print("opened files")

main.obj<-Function.Main(teru.exp.obj, teru.met.hub.obj)

#Write out
saveRDS(main.obj, output.file)
print ("done writing to output")
