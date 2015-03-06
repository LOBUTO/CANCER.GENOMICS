#Function.TEST.SIZE.WILCOX.R
#031615
#Randomly sample one population and compare others at diferent sizes

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.PREP.EXP.MATRIX<-function(exp.rds){
  
  #Load matrix obj
  exp.obj<-readRDS(exp.rds)
  
  #Filter for desired subset - Since we are only interested in cancer, this would be cancer
  exp.matrix<-exp.obj$combined.matrices[,exp.obj$cancer.patients]
  
  #Return
  return(exp.matrix)
}

Function.Main<-function(exp.matrix){
  
  #Randomly sample one cohort of 100 individuals
  set.seed(43)
  test.cohort<-sample(colnames(exp.matrix),100)
  
  #Set rest of patients
  rest<-setdiff(colnames(exp.matrix), test.cohort)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix", "test.cohort", "rest") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-list()
  sizes<-c(2,5,10,20,50,75,100,150,200,400)
  for (pop in sizes){
    print (pop)
    
    #Create random populations lists per each pop size
    random.sampling<-replicate(100, sample(rest, pop), simplify=F) 
    
    #Calculate proportion of diff expressed genes for each random population
    random.prop<-parSapply(cl, random.sampling, function(x){
      
      #Proportion is based on fdr corrected diff exp genes
      wilcox.pvals<-apply(exp.matrix, 1, function(y) wilcox.test(y[x], y[test.cohort], paired=F)$p.value)
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
  print ("Done parallelizing")  
  
  #Return
  return(main.list)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
exp.rds<-args[1]
output.file<-args[2]
print("opened files")
##########################################

##################EXECUTE#################
exp.matrix<-Function.PREP.EXP.MATRIX(exp.rds)
print ("done prepping expression rds")

main.function<-Function.Main(exp.matrix)
print ("Done with main function")

###############WRITING OUTPUT############
saveRDS(object=main.function, output.file)
print ("Done writing to file")