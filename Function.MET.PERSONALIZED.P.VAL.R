#Function MET.PERSONALIZED.P.VAL
#021015
#Calculate P.VAL per metabolite with associated mutations in each individual

##################FUNCTIONS################
library(data.table)
library(reshape2)

Function.Prep.MAF<-function(maf.file) {
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode"),with=F]
  setkey(maf)
  maf<-unique(maf)
  
  #For now classify as change and no-change, later we will also calculate separate probabilities for INS, DEL and non-missense SNPs
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "NO.CHANGE","CHANGE")
  
  #Filter for non-silent only
  maf<-maf[TYPE=="CHANGE",]
  maf$TYPE<-NULL
  
  #Return
  return(maf)
}

Function.Prep.EXON<-function(exon.file){
  
  #Prep length info
  exon<-fread(exon.file, header=T, sep="\t",stringsAsFactors=F)
  
  #Calculate exon length per gene
  exon<-exon[,list(exon.nt=sum(FEAT_END-FEAT_START)),by=Hugo_Symbol]
  
  #Return
  return(exon)
}

Function.Prep.TABLE.2<-function(table.2){
  
  #Introduce metabolic info
  table.2<-readRDS(table.2)
  setnames(table.2, c("METABOLITE","KEGG_ID","Hugo_Symbol"))
  
  #Return
  return(table.2)
}

Function.Empirical.PVALS<-function(maf, table.2, exon){
  
  #Internal functions
  internal.function<-function(sample.n.mut, met.rate, KEGG) {
    
    met.table<-copy(exon)
    
    #Get each gene probabilities for sampling
    gene.prob<-met.table$exon.nt/sum(exon$exon.nt)
    
    #Classify table into gene that are and aren't assocaited with metabolite
    met.table$MET<-met.table$Hugo_Symbol %in% unique(table.2[KEGG_ID==KEGG,]$Hugo_Symbol)
    
    #Simulate sampling 1000 times
    met.dist<-replicate(1000, {
      main<-sample(met.table$Hugo_Symbol, sample.n.mut, replace=T, prob=gene.prob)
      main<-met.table[Hugo_Symbol %in% main, ]
      met.test<-sum(main$MET)
      return(met.test)
    })
    
    #Calcualte p-val from sampling distribution for values that are as extreme or greater than our metabolic rate
    p.val<-mean(met.dist>=met.rate)
    
    #return p-value
    return(list(P.VAL=p.val))
  }
  
  #Obtain number of mutations per patient
  maf[,sample.n.mut:=length(Start_Position), by="Tumor_Sample_Barcode"]
  
  #Merge with metabolic info
  main.table<-merge(maf, table.2, by="Hugo_Symbol")
  
  #Calculate metabolic mutation count per metabolite in each patient
  main.table[,met.mut.count:=length(Start_Position), by=c("Tumor_Sample_Barcode", "METABOLITE","KEGG_ID")]
  
  #Split table to parallelize job
  main.split<-split(main.table, list(main.table$Tumor_Sample_Barcode, main.table$METABOLITE, main.table$KEGG_ID), drop=T)
  print ("Done spliting tables")
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  library(parallel)
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("main.split", "as.data.table", "internal.function", "data.table", "table.2","copy", "exon") ,envir=environment())
  
  #Calculate empirical P.VAL per metabolite per patient (use sample.n.mut, met.rate)
  print (main.table) #erase later
  print ("Simulations")
  
  main.table<-parLapply(cl, main.split, function(x){
    met.table<-x[,internal.function(unique(sample.n.mut), unique(met.mut.count), unique(KEGG_ID)),
                 by=c("Tumor_Sample_Barcode", "METABOLITE","KEGG_ID")]
    return(met.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print("done calculating p-vals")
  
  #Restoing table format
  main.table<-do.call(rbind, main.table)
  
  #Correct for multiple hypothesis testing
  print ("correcting for multiple hypothesis")
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Clean up and return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cancer.maf<-args[1]
exon.file<-args[2]
table.2<-args[3]
output.file<-args[4]
print("opened files")
##########################################

##################EXECUTE#################
maf<-Function.Prep.MAF(cancer.maf)
print ("done prepping maf")

exon<-Function.Prep.EXON(exon.file)
print ("done prepping exon")

table.2<-Function.Prep.TABLE.2(table.2)
print ("done prepping table.2")

main.function<-Function.Empirical.PVALS(maf, table.2, exon)

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)
