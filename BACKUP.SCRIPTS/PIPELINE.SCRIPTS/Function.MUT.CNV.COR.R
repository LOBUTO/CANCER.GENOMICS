#Function.MUT.CNV.COR.R - FIX!!! TAKES TOO LONG TO RUN
#102714
#Finds correlation between non-silent mutations and structural variations across patients
#Need gistic file and table.1

cnv.correlation<-function(patients, gistic.table) {
  require(data.table)
  require(exactRankTests)
  
  internal.function<-function(target.patients,cnv){
    all<-data.table(TARGET=target.patients, CNV=cnv)
    TARGET.CNV<-as.vector(all[TARGET==TRUE,]$CNV)
    NON.TARGET.CNV<-as.vector(all[TARGET==FALSE,]$CNV)
    
    if( (length(TARGET.CNV)!=0) & (length(TARGET.CNV)!=0) ){
      test<-wilcox.exact(TARGET.CNV, NON.TARGET.CNV)
      P.VAL=test$p.value
    } else {
      P.VAL=1
    }
    
    return(list(P.VAL=P.VAL))    
  }
  
  gistic.table$TARGET.PATIENT<-gistic.table$PATIENT %in% patients
  gistic.table<-gistic.table[,internal.function(TARGET.PATIENT, CNV.TH), by="Hugo_Symbol"]
  
  count<<-count+1
  print (count)
  return(list(HUGO.CNV=as.vector(gistic.table$Hugo_Symbol), P.VAL=as.vector(gistic.table$P.VAL)))
}


Function.main.v2<-function(table.1,gistic){ #FIX!!!!
  require(data.table)
  require(parallel)
  require(reshape2)
  require(exactRankTests)
  require(plyr)
  
  #Only keep Missense mutations
  table.1<-copy(table.1$table.1)
  table.1<-table.1[Missense!=0,]
  
  #Create gistic list of gene cnv vectors
  gistic.cnv.vectors<-split(gistic[,2:3,with=F], gistic$Hugo_Symbol)
  gistic.cnv.vectors$A1BG
    
  #Merge table.1 and gistic to get pairwise data.table, rename gistic
  setnames(gistic, c("PATIENT", "CNV.HUGO", "CNV.TH"))
  mock.patients<-sample(unique(as.vector(table.1$PATIENT)),10)
  table.1.CNV<-merge(unique(table.1[PATIENT %in% mock.patients,c(1,2,4), with=F]), gistic[PATIENT %in% mock.patients], by="PATIENT",allow.cartesian=TRUE)
  
  #Filter for setdiff!=integer(0)
  cnv.filter<-table.1.CNV[,list(SET.DIFF.LENGTH=length(setdiff(as.vector(gistic.cnv.vectors[[CNV.HUGO]]$CNV.TH),CNV.TH))), by=c("Hugo_Symbol", "CNV.HUGO")]
  table.1.CNV<-join(table.1.CNV, cnv.filter, by=c("Hugo_Symbol", "CNV.HUGO"))
  table.1.CNV<-table.1.CNV[SET.DIFF.LENGTH!=0,]
  
  #Apply wilcoxon test
  table.1.wilcox<-table.1.CNV[,list(P.VAL=wilcox.test(setdiff(as.vector(gistic.cnv.vector[[CNV.HUGO]]$CNV.TH),CNV.TH), CNV.TH)), 
                              by=c("Hugo_Symbol", "CNV.HUGO")]
  table.1.wilcox$P.VAL.ADJ<-p.adjust(table.1.wilcox$P.VAL, method="fdr")
  table.1.wilcox$SIG<-table.1.wilcox$P.VAL.ADJ<0.05
  
  #Return
  return(table.1.wilcox)
  
}

Function.main<-function(table.1, gistic){
  require(data.table)
  require(parallel)
  require(reshape2)
  
  #Only keep Missense mutations
  table.1<-table.1[Missense!=0,]
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  
  #Split table.1 into number of nodes
  table.1.genes<-unique(as.vector(table.1$Hugo_Symbol))
  table.1.genes<-data.frame(Hugo_Symbol=table.1.genes, GROUP=sample(1:nodes, length(table.1.genes),replace=T))
  table.1<-as.data.table(merge(as.data.frame(table.1), table.1.genes, by="Hugo_Symbol"))
  table.1.split<-split(table.1, table.1$GROUP)
  
  #Introduce parallelizable variables
  clusterExport(cl, varlist=c("count","cnv.correlation", "table.1.split", "gistic","data.table"),envir=environment())
  print ("Done setting up parallelization")
  
  #Execute
  print ("Running parallelization")
  cnv.mut.cor<-parLapply(cl,table.1.split, function(x) {
    x<-as.data.table(x)
    cnv.mut.part<-x[,cnv.correlation(PATIENT, gistic), by="Hugo_Symbol"]
    return(cnv.mut.part)
  })
  
  print ("Done with parallelization")
  stopCluster(cl)
  
  #Merging
  cnv.mut.cor<-do.call(rbind, cnv.mut.cor)
  
  #Multiple hypothesis correction
  cnv.mut.cor$P.VAL.ADJ<-p.adjust(cnv.mut.cor$P.VAL, method="fdr")
  cnv.mut.cor$SIG<-cnv.mut.cor$P.VAL.ADJ<0.05
  
  #Return
  return(cnv.mut.cor)
}

count<<-0
args<-commandArgs(trailingOnly=T)

table.1.obj<-readRDS(args[1]) 
gistic<-readRDS(args[2])
output.file<-args[3]

main.object<-Function.main(table.1.obj$table.1, gistic)
print ("Done with main function")

saveRDS(object=main.object, file=output.file)
print ("done saving")
