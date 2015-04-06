######Function.METABOLIC.COR.BACKGROUND.R#######
#040615
#Calculate bacgkground deviation of correlation score for all sample sizes in results from Function.METABOLIC.COR.R

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.Cor<-function(cor.obj){
  
  #Load obj
  obj<-readRDS(cor.obj)
  
  #Extract unique sample sizes
  n.samples<-sapply(obj, function(x) x[["N.SAMPLES"]])
  n.samples<-unique(n.samples)
  
  #Return
  return(n.samples)
}

Function.exp<-function(exp.obj){
  
  #Load expresssion file
  BRCA.EXP<-readRDS(exp.obj)
  
  #Obtain normal and cancer correlations
  #BRCA.EXP.COR.NORMAL<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$normal.patients]), method="spearman")
  #BRCA.EXP.COR.CANCER<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients]), method="spearman")
  
  #Return
  return(list(NORMAL=BRCA.EXP$combined.matrices[,BRCA.EXP$normal.patients],
              CANCER=BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients]))
}

Function.Enzyme<-function(enzyme.file){
  
  #Read enzyme file
  ENZYME<-fread(enzyme.file, header=T, sep="\t", stringsAsFactors=F)
  
  #Return
  return(ENZYME)
  
}

Function.Main<-function(n.cor, cancer.exp, normal.exp, enzymes){
  
  #Filter exp matrix for target enzymes
  cancer.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  normal.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  
  #Obtain background enzymatic correlation
  BRCA.EXP.COR.NORMAL<-cor(t(normal.exp), method="spearman")
  genes.order<-rownames(BRCA.EXP.COR.NORMAL)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "n.cor",  "cancer.exp", "maf", "BRCA.EXP.COR.NORMAL",
                              "genes.order") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-lapply(n.cor, function(x) {
    print (x)
    
    #Obtain random samples
    ran.samples<-replicate(1000,sample(colnames(cancer.exp), x), simplify=F)
    
    #Ran correlation simulations
    simul.list<-parLapply(cl, ran.samples, function(y) {
      
      #Construct correlation matrix for random samples
      main.cor<-cor(t(cancer.exp[,y]), method="spearman")
      
      #Substract
      signal.cor<-main.cor[gene.order, gene.order]-back.cor
      
      #Get sum of absolute value of change of enzyme with respect to all - LOG-TRANSFORMED
      sum.delta.abs<-log(apply(signal.cor, 1, function(z) sum(abs(z))))
      
      #Return
      return(sum.delta.abs)
    })
    
    #Assign to genes
    main.table<-do.call(cbind, simul.list)
    rownames(main.table)<-gene.order
    
    #Clean up of NAs
    main.table<-t(main.table)
    main.table<-main.table[complete.cases(main.table),]
    main.table<-t(main.table)
    
    #Return
    return(main.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Assign cluster size names to list of random matrices produced during simulation
  names(main.list)<-as.character(n.cor)
  
  #Return obj
  return(main.list)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cor.obj<-args[1]
exp.rds<-args[2]
enzyme.file<-args[3]
output.file<-args[4]
print("opened files")
##########################################

##################EXECUTE#################
n.cor<-Function.Prep.Cor(cor.obj)
print ("done extracting correlation sample sizes")

exp<-Function.exp(exp.rds)
print ("done building cancer and normal correlations")

enzyme<-Function.Enzyme(enzyme.file)
print ("done reading enzyme file")

main.function<-Function.Main(n.cor, exp$CANCER, exp$NORMAL, unique(as.vector(enzyme$Enzyme)))
print ("done with main function")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing to file")