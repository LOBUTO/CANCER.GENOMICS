#######Function.MI.CANCER.ANALYSIS.R########
#030415


##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Process.Cluster<-function(cluster.file){
  #Processes clusters into lists of gene sets
  
  #Load file
  clusters<-scan(cluster.file, what="", sep="\n")
  
  #Process into list of cluster gene vectors
  cluster.list<-lapply(clusters, function(x) unlist(strsplit(x, "\t")))
  
  #Name clusters to keep track
  name(cluster.list)<-letters[1:length(cluster.list)]
  
  #Return
  return(cluster.list)
}

Function.SIMPLIFY.NORMAL.MI<-function(normal.mi.table, cluster.list){
  #Function to simplify normal.mi.table based on only genes found in cluster.list
  
  #Obtain total gene list
  cluster.genes<-unlist(cluster.list)
  
  #First filtering passage on normal.mi.table
  normal.mi<-normal.mi.table[(Hugo.1 %in% cluster.genes) | (Hugo.2 %in% cluster.genes),]
  print ("Done with first pass of mi filtering")
  
  #Construct tables Hugo.1 - Hugo.2 table according to order in function Function.CANCER.COHORT.CLUSTER.MI
  main.list<-lapply(cluster.list, function(x) {
    comb.table<-data.table(t(combn(x,2)))
    return(comb.table)
  })
  
  main.table<-do.call(rbind, main.list)
  setnames(main.table, c("Hugo.1", "Hugo.2"))
  
  #Double normal.mi structure
  normal.mi.comp<-copy(normal.mi[,c("Hugo.2","Hugo.1", "MI"), with=F])
  setnames(normal.mi.comp, c("Hugo.1", "Hugo.2", "MI"))
  normal.mi<-rbind(normal.mi, normal.mi.comp)
  
  #Second filtering of normal.mi by gene pairs of interest
  normal.mi<-merge(normal.mi, main.table, by=c("Hugo.1", "Hugo.2"))
  print("Done with second pass of mi filtering")
  
  #Return
  return (normal.mi)
  
}

Function.CANCER.COHORT.CLUSTER.MI<-function(cluster.list, samples, cancer.matrix) {
  #Calculates the mutual information for all pairs of genes in each cluster for the designated cancer samples
  require(entropy)
  require(parallel)
  
  #Filter cancer expression matrix for patients of interest
  cancer.matrix<-cancer.matrix[,samples]
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","cancer.matrix", "discretize2d", "mi.empirical", "cluster.list") ,envir=environment())
  print ("Done exporting values")
  
  #Calculate mutual information
  cluster.mi.list<-lapply(names(cluster.list), function(x) {
    
    #Obtain genes in cluster
    genes<-cluster.list[[x]]
    
    #Construct gene comnb table
    comb.table<-data.table(t(combn(genes, 2)))
    
    #Calculate MI per gene comb pair in table
    comb.mi<-parApply(cl,comb.table, 1, function(y){
      
      #Obtain expression matrix for genes in loop
      exp.1<-cancer.matrix[y[1],]
      exp.2<-cancer.matrix[y[2],]
      
      #Calculate MI for pair
      gene.disc<-discretize2d(exp.1, exp.2, numBins1=10, numBins2=10)
      mi.pair<-mi.empirical(gene.disc)
      
      #Return MI
      return(mi.pair)
    })
    
    #Assign MIs and cluster name to pairs
    comb.table$MI<-comb.mi
    comb.table$CLUSTER<-x
    
    #Clean up and return table
    setnames(comb.table, c("Hugo.1", "Hugo.2", "CANCER.MI", "CLUSTER"))
    return(comb.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")  
  
  #Combine cluster list into one
  cluster.mi.table<-do.call(rbind, cluster.mi.list)
  
  #Return
  return(cluster.mi.table)
}

Function.CLUSTERS.MI.WICOXON<-function(filtered.normal.mi, cancer.mi) {
  #Calculates the wilcoxon statistic between the cancer clusters and our clusters
  
  #Merge cancer and normal mi tables
  main.table<-merge(filtered.normal.mi, cancer.mi)
  
  #Calcualte wilcoxon per cluster
  #NOTE: This will be a paired test, so p-value when comparing cancer to normal per cluster across all mutated genes will be independent of sample size
  
  
  
}