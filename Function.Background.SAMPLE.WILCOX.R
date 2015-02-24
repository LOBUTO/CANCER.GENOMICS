#Function.Background.SAMPLE.WILCOX.R
#022315
#Function to build background distributions for wilcoxon diff expression across number of samples

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Main<-function(linkage.nc.file, exp.rds){
  
  #Load linkage file
  main.table<-fread(linkage.nc.file, header=T, sep="\t", stringsAsFactors=F)
  print (main.table)
  print ("Done loading linkage file")
  
  #Obtain vector of population sizes
  population<-unique(main.table$N.PATIENTS)
  population<-population[1:2]
  print (population)
  
  #Load expression file
  exp.matrix<-readRDS(exp.rds)
  print ("Done loading expression file")
  
  #Subset samples
  normal<-unique(exp.matrix$normal.patients)
  cancer<-unique(exp.matrix$cancer.patients)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("population", "normal","cancer", "as.data.table","exp.matrix","data.table") ,envir=environment())
  print ("Done exporting values")
  
  #Executing
  print ("Finding background")
  main.list<-list()
  
  #Loop for each sample size in population of samples
  for (pop in population){
    
    #Create random populations of pop size
    random.sampling<-replicate(100, sample(cancer,pop), simplify=F)
    
    #Calculate proportion of diff expressed genes for each random population
    random.prop<-parSapply(cl, random.sampling, function(x){
      
      #Proportion is based on fdr corrected diff exp genes
      wilcox.pvals<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[x], y[normal], paired=F)$p.value)
      wilcox.vector<-p.adjust(as.vector(wilcox.pvals), method="fdr" )
      wilcox.ratio<-mean(wilcox.vector<0.05)
      
      #Return diff exp ratio
      return(wilcox.ratio)
      
    })
    
    #Append random diff exp ratio population to list
    main.list[[as.character(pop)]]<-random.prop
  }
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with background distributions")
  
  #Return
  return(main.list)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
linkage.nc.file<-args[1]
exp.rds<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
main.function<-Function.Main(linkage.nc.file, exp.rds)
print ("Done")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing file")