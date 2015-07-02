library(parallel)
library(data.table)
library(reshape2)
library(ROCR)
library(pROC)

Function.met.scoring<-function(network.table, met.targets, d=0.85, d.filter=1000) {
  
  #Filter network table by degree (filter out currency metabolites)
  degree.filter<-unique(rbind(data.table(ID.1=network.table$SUBSTRATE, ID.2=network.table$PRODUCT),
                              data.table(ID.1=network.table$PRODUCT, ID.2=network.table$SUBSTRATE)))
  degree.filter<-degree.filter[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
  degree.filter<-degree.filter[DEGREE<d.filter,]
  
  network.table<-network.table[SUBSTRATE %in% unique(degree.filter$ID.1),][PRODUCT %in% unique(degree.filter$ID.1),]
  
  #Create transition matrix
  t.m<-network.table[,list(DEGREE=length(unique(PRODUCT))), by="SUBSTRATE"]
  t.m<-merge(t.m, network.table[,c("SUBSTRATE", "PRODUCT"), with=F], by = "SUBSTRATE")
  t.m$trans<-1/t.m$DEGREE #may modify depending on weighted jaccard edges
  degree.info<-unique(data.table(SUBSTRATE=t.m$SUBSTRATE, DEGREE=t.m$DEGREE))
  degree.info<-rbind(degree.info, data.table(SUBSTRATE=unique(t.m$PRODUCT)[!(unique(t.m$PRODUCT) %in% degree.info$SUBSTRATE)], DEGREE=0))
  t.m$DEGREE<-NULL
  t.m<-rbind(t.m, data.table(SUBSTRATE=unique(t.m$SUBSTRATE), PRODUCT=unique(t.m$SUBSTRATE), trans=0))
  
  all.mets<-unique(c(t.m$SUBSTRATE, t.m$PRODUCT))
  adj.t.m<-matrix(ncol=length(all.mets), nrow=length(all.mets), dimnames = list(all.mets, all.mets))
  adj.t.m[is.na(adj.t.m)]<-0
  apply(as.matrix(t.m), 1, function(x) {
    adj.t.m[x[2], x[1]]<<-as.numeric(x[3])
  })
  
  #Obtain personalized vector for path
  #personalized.vector<-unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  personalized.vector<-met.targets[met.targets %in% rownames(adj.t.m)]
  personalized.vector<-ifelse(rownames(adj.t.m) %in% personalized.vector, 1/length(personalized.vector), 0)
  
  #####Iterate to minimize mean change per row and converge to stable network####
  initial.rank<-rep(1/nrow(adj.t.m), nrow(adj.t.m))
  
  #Obtain first iteration page rank
  A.1<-(adj.t.m*d)%*%initial.rank + (1-d)*personalized.vector #Implementing personalized vector
  
  #Iterate to converge on optimal page rank
  delta<-Inf
  while(delta>0.0000000001){
    delta<-mean(abs(A.1 - ((adj.t.m*d)%*%A.1 + (1-d)*personalized.vector)))
    A.1<-((adj.t.m*d)%*%A.1 + (1-d)*personalized.vector)
  }
  
  #####Obtained the non-personalized version####
  personalized.vector<-rep(1/nrow(adj.t.m), nrow(adj.t.m)) #all equal for all
  
  #Obtain first iteration page rank
  A.2<-(adj.t.m*d)%*%initial.rank + (1-d)*personalized.vector #Implementing personalized vector
  
  #Iterate to converge on optimal page rank
  delta<-Inf
  while(delta>0.0000000001){
    delta<-mean(abs(A.2 - ((adj.t.m*d)%*%A.2 + (1-d)*personalized.vector)))
    A.2<-((adj.t.m*d)%*%A.2 + (1-d)*personalized.vector)
  }
  ##############################################
  A.1<-as.vector(A.1)
  A.2<-as.vector(A.2)
  PAGE.RANK<-A.1-A.2
  
  #Clean up and return
  main.table<-data.table(SUBSTRATE=rownames(adj.t.m), PAGE.RANK, A.1, A.2)
  main.table$PAGE.RANK<-main.table$PAGE.RANK / sum(main.table$PAGE.RANK) #Normalizing page rank to 1
  main.table<-merge(main.table, degree.info, by="SUBSTRATE")
  main.table$TARGET.PATH<-main.table$SUBSTRATE %in% met.targets
  #main.table$TARGET.PATH<-main.table$SUBSTRATE %in% unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  main.table<-main.table[order(PAGE.RANK, decreasing = T),]
  return(main.table)
}

Function.HOO.ROC<-function(HOO){
  require(ROCR)
  require(pROC)
  
  #Remove true targets that were not held out
  HOO<-HOO[TARGET.PATH!=T,]
  
  #For each hold out obtain normalized ranks and predictions
  hold.outs<-unique(HOO$HOLD.OUT)
  
  pred.vector<-c()
  pred.scores<-c()
  
  for (hold in hold.outs){
    hoo.target<-HOO[HOLD.OUT==hold,]
    hoo.scores<-rank(hoo.target$PAGE.RANK, ties.method="max")
    hoo.pred<-hoo.target$SUBSTRATE %in% hold.outs
    
    pred.vector<-c(pred.vector, hoo.pred)
    pred.scores<-c(pred.scores, hoo.scores)
  }
  
  #Obtain ROC
  roc.pred<-prediction(pred.scores, pred.vector)
  roc.perf<-performance(roc.pred, "tpr", "fpr")
  
  #Clean up and return
  main.table<-data.table(FPR=roc.perf@x.values[[1]], TPR=roc.perf@y.values[[1]])
  print (auc(pred.vector, pred.scores))
  return(list(ROC=main.table, AUC=auc(pred.vector, pred.scores)))
}

Function.HOO.PR<-function(target.path, uni.network, kegg.path, d=0.85){
  
  #Obtain kegg targets
  kegg.targets<-unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  kegg.targets<-kegg.targets[kegg.targets %in% unique(uni.network$SUBSTRATE)]
  
  #Get page ranks with hold one out if any target present
  main.table<-data.table()
  if (length(kegg.targets)>1){
    for (cpd in kegg.targets){
      hold.out<-cpd
      train<-setdiff(kegg.targets, hold.out)
      
      #Filter kegg.path
      kegg.path.filt<-kegg.path[DESCRIPTION==target.path,][COMPOUND %in% train,]
      
      #Obtain page ranks
      train.pr<-Function.met.scoring(uni.network, kegg.path.filt$COMPOUND, d)
      train.pr$HOLD.OUT<-hold.out
      train.pr$SUBSTRATE.HOO<-train.pr$SUBSTRATE==hold.out
      
      #Store
      main.table<-rbind(main.table, train.pr)
    }
  } else {
    print ("not enough metabolites for page rank hold out")
  }
  
  #Return
  return(main.table)
}

Function.Main.HOO.ROC<-function(kegg.path, recon.directed, path.filter.l=3, path.filter.r=50, d=0.85, degree.filter=1000){
  
  #####Prefiltering steps
  #Consider metabolites that pass the degree filter (This is to filter high-degree currency metabolites)
  recon.directed[,S.DEGREE:=length(unique(PRODUCT)), by="SUBSTRATE"]
  recon.directed<-recon.directed[S.DEGREE<degree.filter,]
  recon.directed$S.DEGREE<-NULL
  
  #Consider paths that have enough metabolites found in network
  kegg.path<-kegg.path[COMPOUND %in% unique(recon.directed$SUBSTRATE),]
  kegg.path[,N.MET:=length(unique(COMPOUND)),DESCRIPTION]
  kegg.path<-kegg.path[N.MET<path.filter.r & N.MET>path.filter.l,]
  kegg.path$N.MET<-NULL
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - "))[1])
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "kegg.path", "Function.HOO.PR",
                              "Function.HOO.ROC", "recon.directed", "kegg.path", "d", "Function.met.scoring") ,envir=environment())
  
  #Obtain AUC for all cancer paths that pass initial threshold
  #paths<-unique(kegg.path$DESCRIPTION)
  paths<-unique(kegg.path[grepl("cancer", DESCRIPTION, ignore.case = T),]$DESCRIPTION)
  print (length(paths))
  main.auc<-parSapply(cl, paths, function(x) {
    print (x)
    hoo.path<-Function.HOO.PR(x, recon.directed, kegg.path, d)
    hoo.roc<-Function.HOO.ROC(hoo.path)[["AUC"]][[1]]
    return(hoo.roc)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Construct main table
  main.table<-data.table(PATH=paths, AUC=main.auc)
  
  #Clean up and return
  main.table<-main.table[order(AUC, decreasing = T),]
  return(main.table)
}

#Arguments
args<-commandArgs(trailingOnly=T)
kegg.path<-fread(args[1], header=T, sep="\t")
recon.directed<-readRDS(args[2])
output.file<-args[3]
print ("done loading files")

#Execute
d.filters=seq(0.05,0.95,0.05)
degree.filters<-c(60,70,80,90,100,150,250,300)

MAIN.OBJ<-data.table()
for (d.f in d.filters){
  for (deg in degree.filters){
    print (c(d.f, deg))
    TEMP.ROC<-Function.Main.HOO.ROC(kegg.path, recon.directed, path.filter.l=3, path.filter.r=60, d=d.f, degree.filter=deg)
    TEMP.ROC$d<-d.f
    TEMP.ROC$DEGREE.FILTER<-deg
    MAIN.OBJ<-rbind(MAIN.OBJ, TEMP.ROC)
  }
}

saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")