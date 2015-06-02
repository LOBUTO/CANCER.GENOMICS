
library(parallel)
library(data.table)

Function.predict.features.lm<-function(target.matrix, feat, target.feat){
  
  #1000 ITERATIONS FOR BEST FEATURE SELECTION
  
  feat.max<-nrow(target.matrix)-1
  main.list<-list()
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "target.matrix", "feat", "target.feat") ,envir=environment())
  print ("Done exporting values")
  
  main.list<-parLapply(cl, 1:1000, function(x){
    print (x)
    
    #Randomize feature samples after each iteration
    features<-sample(feat)
    
    #Reset AIC
    current.AIC<-Inf
    current.feat<-c("ER.STATUS")
    count=1
    f.adj.sq<-1
    
    #Continue step-wise AIC search till we reach maximum allowed features
    while(is.na(f.adj.sq)==FALSE | count<length(features)){
      testing.feat<-features[count]
      
      f.formula<-as.formula(paste(target.feat ," ~ ", paste(c(current.feat, testing.feat), collapse= "+")))
      print (f.formula)
      f.lm<-lm(f.formula, target.matrix)
      f.aic<-extractAIC(f.lm)[2]
      f.adj.sq<-summary(f.lm)$adj.r.squared
      
      #Keep feature if lower than current AIC and update AIC
      if ((f.aic<current.AIC) & (is.na(f.adj.sq)==FALSE)){
        current.feat<-c(current.feat, testing.feat)
        current.AIC<-f.aic
      }
      print (c(current.feat, current.AIC))
      
      #Remove features if we are below 0.5 adjusted r-square and we haven't improved the aic
      if ((f.adj.sq<0.5) & (f.aic<current.AIC) ){
        feat<<-setdiff(feat, testing.feat) #Notice we remove from global feat, and not from local features, so next count will still keep the flow of search for next upcoming feature
      }
      
      #Update count
      count<-count+1
    }
    return(current.feat)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  return(main.list)
}

#Arguments
args<-commandArgs(trailingOnly=T)
teru.obj<-args[1]
output.file<-args[2]

TERU.2HG.EXP<-readRDS(teru.obj)
TERU.2HG.SELECTED.FEATURES<-Function.predict.features.lm(data.frame(TERU.2HG.EXP), setdiff(colnames(data.frame(TERU.2HG.EXP)),c("METABOLITE","ER.STATUS")) ,"METABOLITE")

saveRDS(object = TERU.2HG.SELECTED.FEATURES, file = output.file)
print ("Done writing output")