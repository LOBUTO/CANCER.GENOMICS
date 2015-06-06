library(parallel)
library(data.table)
library(reshape2)
library(Hmisc)

Function.Main<-function(brca.exp, recon.table, HMDB.table, alpha=0.3){
  #This will use pre-filtered recon.table that contains either a viable KEGG or HMDB identifier
  
  #Double check and clean up expression datasets
  brca.exp$normal<-brca.exp$normal[complete.cases(brca.exp$normal),]
  brca.exp$tumor<-brca.exp$tumor[complete.cases(brca.exp$tumor),]
  
  #Get unique identifiers for recon table (Use primarily kegg id, only use hmdb id if kegg id not present)
  kegg.ids<-recon.table[KEGG_ID!="NONE",]
  hmdb.ids<-recon.table[KEGG_ID=="NONE" & HMDB!="NONE",]
  
  kegg.ids<-kegg.ids[,c("Hugo_Symbol", "KEGG_ID"),with=F]
  setnames(kegg.ids, c("Hugo_Symbol", "ID"))
  hmdb.ids<-hmdb.ids[,c("Hugo_Symbol", "HMDB"), with=F]
  setnames(hmdb.ids, c("Hugo_Symbol", "ID"))
  
  recon.table<-unique(rbind(kegg.ids, hmdb.ids))
  print ("Processed recon table")
  
  #Add HMDB processed KEGG targeted table
  print (recon.table)
  hmdb.source<-fread(HMDB.table, header=T, sep="\t", stringsAsFactors = F)
  hmdb.source<-hmdb.source[,c("KEGG_ID", "Hugo_Symbol"),with=F]
  setnames(hmdb.source, c("ID", "Hugo_Symbol"))
  print (hmdb.source)
  recon.table<-unique(rbind(recon.table, hmdb.source[,c("Hugo_Symbol","ID"),with=F]))
  print (recon.table)
  
  #Filter recon.table and brca.exp for common genes
  common.genes<-intersect(unique(recon.table$Hugo_Symbol), intersect(rownames(brca.exp$tumor), rownames(brca.exp$normal)))
  recon.table<-recon.table[Hugo_Symbol %in% common.genes,]
  brca.exp$tumor<-brca.exp$tumor[common.genes,]
  brca.exp$normal<-brca.exp$normal[common.genes,]
  
  #Filter recon.table for mets that have at least 3 genes assigned to them
  print (length(unique(recon.table$ID)))
  recon.count<-recon.table[,list(N.GENES=length(Hugo_Symbol)), by ="ID"]
  recon.count<-recon.count[N.GENES>=3, ]
  recon.table<-recon.table[ID %in% recon.count$ID,]
  print (length(unique(recon.table$ID)))
  
  #Calculate correlation matrix in cancer and normal
  cor.normal<-rcorr(t(data.matrix(brca.exp$normal)), type = "spearman")
  cor.tumor<-rcorr(t(data.matrix(brca.exp$tumor)), type = "spearman")
  print ("Done with correlation tables")
  
  #Correct correlation p-values for multiple hypothesis testing and obtain their signficance status (P-val<0.05)
  cor.normal$P<-matrix(p.adjust(cor.normal$P, method = "fdr"), ncol=ncol(cor.normal$P), 
                       dimnames = list(rownames(cor.normal$P), colnames(cor.normal$P)))
  cor.tumor$P<-matrix(p.adjust(cor.tumor$P, method = "fdr"), ncol=ncol(cor.tumor$P), 
                      dimnames = list(rownames(cor.tumor$P), colnames(cor.tumor$P)))
  cor.normal$P<-cor.normal$P<0.05
  cor.tumor$P<-cor.tumor$P<0.05
  diag(cor.normal$P)<-FALSE 
  diag(cor.tumor$P)<-FALSE
  
  #Use corrrelation value only if it passes multiple hypothesis testing (Self correlation value will be zero)
  cor.normal$r<-cor.normal$r * cor.normal$P
  cor.tumor$r<-cor.tumor$r * cor.tumor$P
  print ("Done with fdr correcting correlation tables")
  
  #Obtain gene pool and counts for replication
  gene.pool<-unique(as.vector(recon.table$Hugo_Symbol))
  gene.counts<-unique(recon.table[,list(N.GENES=length(unique(Hugo_Symbol))), by="ID"]$N.GENES)
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "gene.counts", "gene.pool",
                              "cor.normal", "cor.tumor", "alpha") ,envir=environment())
  
  #Random 100x per gene count
  main.list<-lapply(gene.counts, function(x) {
    print (x)
    
    #Obtain 100 samples per gene count
    rep.sample.genes<-replicate(100, sample(gene.pool, x), simplify = F)
    
    #Get score per random sample
    rep.scores<-parSapply(cl, rep.sample.genes, function(y) {
      
      norm.vec<-cor.normal$r[y, y][upper.tri(cor.normal$r[y,y])]
      cancer.vec<-cor.tumor$r[y, y][upper.tri(cor.tumor$r[y,y])]
      SCORE=mean(abs(cancer.vec - norm.vec)>alpha)
      
      #Return score for randoms sample
      return(SCORE)
    })
    #Return vector of scores for replications per gene count
    return(rep.scores)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Make sure to name lists with the corresponding gene count
  names(main.list)<-as.character(gene.counts)
  
  #Clean up and return
  return(main.list)
}

#Arguments
args<-commandArgs(trailingOnly=T)
brca.exp<-readRDS(args[1])
recon.table<-readRDS(args[2])
HMDB.table<-args[3]
alpha=as.numeric(args[4])
output.file<-args[5]
print ("done loading files")

MAIN.OBJ<-Function.Main(brca.exp, recon.table, HMDB.table, alpha=alpha)

saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")