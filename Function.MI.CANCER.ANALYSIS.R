#######Function.MI.CANCER.ANALYSIS.R########
#030415


##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)
library(ggplot2)

theme.format<-theme(axis.text.y=element_text(size=rel(2.5)), axis.text.x=element_text(size=rel(2.5)),
                    axis.title.y=element_text(size=22), axis.title.x=element_text(size=22),
                    legend.text = element_text(size = 22))

Function.Process.Cluster<-function(cluster.file){
  #Processes cluster files from SPICI results into lists of gene sets
  
  #Load file
  clusters<-scan(cluster.file, what="", sep="\n")
  
  #Process into list of cluster gene vectors
  cluster.list<-lapply(clusters, function(x) unlist(strsplit(x, "\t")))
  
  #Name clusters to keep track
  names(cluster.list)<-paste("C", 1:length(cluster.list) , sep=".")
  
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
  samples<-intersect(samples, colnames(cancer.matrix))
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

Function.CLUSTERS.MI.WICOXON<-function(filtered.normal.mi, cancer.mi, figures.folder, gene.pos, samples) {
  #Calculates the wilcoxon statistic between the cancer clusters and our clusters
  
  #Merge cancer and normal mi tables
  main.table<-merge(filtered.normal.mi, cancer.mi, by=c("Hugo.1", "Hugo.2"))
  
  #Number of samples covered by gene.pos mutation
  n.samples<-length(samples)
  
  #Calcualte wilcoxon per cluster
  #NOTE: This will be a paired test, so p-value when comparing cancer to normal per cluster across all mutated genes will be independent of sample size
  main.wilcox<-main.table[,list(P.VAL=wilcox.test(MI, MI.CANCER, paired=T)$p.value),
                         by="CLUSTER"]
  
  #Correct for multiple hypothesis testing
  main.wilcox$P.VAL.ADJ<-p.adjust(main.wilcox$P.VAL, method="fdr")
  
  #Split table to plot figures
  main.list<-split(main.table, main.table$CLUSTER)
  
  #Plot figures
  for (cluster in main.list) {
    cluster.name<-unique(cluster$CLUSTER)
    
    #Melt table to separate for plotting 
    cluster$PAIR.ID<-paste("G", c(1:nrow(cluster)), sep=".")
    cluster.melt<-melt(cluster[,c("MI","CANCER.MI","PAIR.ID"), with=F], id.vars=c("PAIR.ID"))
    setnames(cluster.melt, c("PAIR.ID", "TYPE", "MI"))
    
    #Plot
    my.plot<-ggplot(cluster.melt, aes(TYPE, MI, group=PAIR.ID, colour=PAIR.ID)) + geom_line() + theme.format+
      geom_point( size=4, shape=21, fill="white") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste(gene.pos, " - ","CLUSTER:",cluster.name," - ", n.samples, sep=""))
    
    filename=paste(gene.pos, cluster.name, sep=".")
    jpeg(filename=paste(figures.folder, filename, sep="/") ,width=1200,height=800, quality=100,type="quartz")
    print (my.plot)
    dev.off() 
  }
  
  #Clean up and Return 
  main.wilcox<-main.wilcox[order(P.VAL.ADJ),]
  return(main.wilcox)
}

###TESTING###
test.process.file<-Function.Process.Cluster("PIPELINES/METABOLIC.DRIVERS/NETWORKS/SPICI.RESULTS/BRCA.EXP.NORMAL.TH.0.75.cluster")

test.simplify<-Function.SIMPLIFY.NORMAL.MI(BRCA.NORMAL.MI, test.process.file)

PIK3CA<-unique(brca.maf[Hugo_Symbol=="PIK3CA",]$SAMPLE)
TTN<-unique(brca.maf[Hugo_Symbol=="TTN",]$SAMPLE)
PIK3CA.NO.TTN<-setdiff(PIK3CA, TTN)
PIK3CA.NO.TTN.cohort.cluster.mi<-Function.CANCER.COHORT.CLUSTER.MI(test.process.file, PIK3CA.NO.TTN, brca.exp.nb$combined.matrices)

