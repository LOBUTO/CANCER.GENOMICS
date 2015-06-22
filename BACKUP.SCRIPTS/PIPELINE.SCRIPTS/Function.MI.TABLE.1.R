#Function to calculate mutual information for genes across patients

Function.main<-function(table.1){
  #Calculates the mutual information for genes across patients
  
  require(data.table)
  require(infotheo)
  require(reshape2)
  require(parallel)
  
  #Clean up table.1 - Keep only Missense
  table.1<-table.1[Missense!=0,]
  
  #
  
}

table.1.obj<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/100514.BRCA.Table.1.rds")
table.1.obj$table.1
dual.table.1<-acast(data=table.1.obj$table.1[Missense!=0,][Hugo_Symbol %in% unique(as.vector(table.1.obj$table.1$Hugo_Symbol)),], 
                    PATIENT~Hugo_Symbol, value.var="Missense", fill=0)
table.mi<-mutinformation(X=as.data.frame(dual.table.1))
pheatmap(table.mi, scale="none")
system.time(mutinformation(X=as.data.frame(dual.table.1)))

plot(c(10,100,200,500), c(0.001,0.367,2.894,42.733))

dim(dual.table.1)
dual.table.1[1:3,1:3]
dual.table.1.combn<-as.data.table(t(combn(colnames(dual.table.1),2)))
head(dual.table.1.combn)

dual.table.1.combn
table.mi<-sapply(1:100, function(x) mutinformation(dual.table.1[,as.character(dual.table.1.combn[x,1,with=F])],
                                                 dual.table.1[,as.character(dual.table.1.combn[x,2,with=F])]))

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("dual.table.1","dual.table.1.combn", "mutinformation"),envir=environment())

system.time(parSapply(cl,1:1000, function(x) mutinformation(dual.table.1[,as.character(dual.table.1.combn[x,1])],
                                                     dual.table.1[,as.character(dual.table.1.combn[x,2])])))

stopCluster(cl)