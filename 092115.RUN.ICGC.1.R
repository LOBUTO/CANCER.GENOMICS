library(data.table)
library(reshape2)
library(parallel)

#Functions
Function.icgc.enrich.1<-function(icgc.mut, brca.exp, table.2, edges=F){
  #Finds enriched driver metabolites based on mutation and expression analysis
  
  #Are we supplying kegg.edges or table.2?
  if (edges==T){
    table.2<-rbind(data.table(KEGG.ID=table.2$SUBSTRATE, Hugo_Symbol=table.2$Hugo_Symbol),
                   data.table(KEGG.ID=table.2$PRODUCT, Hugo_Symbol=table.2$Hugo_Symbol))
    table.2<-unique(table.2)
  } else {
    table.2<-unique(table.2[,c("KEGG_ID", "GENE"),with=F])
    setnames(table.2, c("KEGG.ID", "Hugo_Symbol")) 
  }
  
  #Filter expression and mutation info against each other for common samples
  common.samples<-intersect(unique(icgc.mut$SAMPLE), unique(colnames(brca.exp$tumor)))
  icgc.mut<-icgc.mut[SAMPLE %in% common.samples,]
  brca.exp$tumor<-brca.exp$tumor[,common.samples]
  
  #Filter expression matrix by non-zero expression
  brca.exp$tumor<-brca.exp$tumor[rowSums(brca.exp$tumor)!=0,]
  
  #Assign sample mutations to kegg identifiers based on associated genes
  icgc.mut<-icgc.mut[MUTATION=="MISSENSE",]
  mut.table<-merge(icgc.mut, table.2, by="Hugo_Symbol")
  
  #Filter mutation table by those metabolites that have at least 10 samples assocaited to them
  mut.table[,N.SAMPLES:=length(unique(SAMPLE)), by="KEGG.ID"]
  mut.table<-mut.table[N.SAMPLES>=10,]
  
  #Filter expression metabolites by those that have at least 4 genes associated to them
  exp.table<-table.2[Hugo_Symbol %in% rownames(brca.exp$tumor),]
  exp.table[,N.GENES:=length(unique(Hugo_Symbol)), by="KEGG.ID"]
  exp.table<-exp.table[N.GENES>=4,]
  exp.mets<-unique(exp.table$KEGG.ID)
  
  #Iteration through filtered mutation mets and obtain significantly differentiated values for expression mets
  count<-1
  internal.function<-function(samples, KEGG.ID){
    
    not.samples<-setdiff(colnames(brca.exp$tumor), samples)
    met.pvals<-parSapply(cl, exp.mets, function(x) {
      met.hugos<-unique(exp.table[KEGG.ID==x, ]$Hugo_Symbol)
      PVAL=wilcox.test(rowMeans(brca.exp$tumor[met.hugos, samples]), rowMeans(brca.exp$tumor[met.hugos, not.samples]), paired = T)$p.value
      return(PVAL)
    })
    main.mets<-data.table(EXP.METS=exp.mets, PVAL=met.pvals)
    main.mets$PVAL.ADJ<-p.adjust(main.mets$PVAL, method="fdr")
    
    count<<-count+1
    print (count)
    return(list(EXP.METS=main.mets$EXP.METS, PVAL=main.mets$PVAL, PVAL.ADJ=main.mets$PVAL.ADJ))
  }
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "mut.table", "internal.function", "brca.exp", "exp.table", "exp.mets",
                              "count") ,envir=environment())
  print ("Done exporting values")
  
  #Executing
  print (length(unique(mut.table$KEGG.ID)))
  main.table<-mut.table[,internal.function(unique(SAMPLE), unique(KEGG.ID)),by="KEGG.ID"]
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Return
  return(main.table)
}

#Arguments
args<-commandArgs(trailingOnly=T)
brca.icgc.mut<-readRDS(args[1])
brca.exp<-readRDS(args[2])
table.2<-readRDS(args[3])
output.file<-args[4]

#Execute
icgc.brca.enrich.1<-Function.icgc.enrich.1(brca.icgc.mut, brca.exp, table.2, edges=T)

#Save
saveRDS(object = icgc.brca.enrich.1 , file = output.file)
print ("Done writing output")