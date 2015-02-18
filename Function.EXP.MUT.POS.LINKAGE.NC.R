#Function.EXP.MUT.POS.LINKAGE.NC.R
#021515
#Calculates the linkage between mutation position and expression between normal and cancer populations
#Will initially only compare paired samples
#Will only take into account sites that are seen mutated at least 2 times in population

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.MAF<-function(maf.file) {
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #For now classify as change and no-change, later we will also calculate separate probabilities for INS, DEL and non-missense SNPs
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "NO.CHANGE","CHANGE")
  
  #Filter for non-silent only
  maf<-maf[TYPE=="CHANGE",]
  maf$TYPE<-NULL
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Return
  return(maf)
}

Function.Prep.EXP<-function(exp.rds, paired=T){
  
  #Load expression matrix file
  exp.matrix<-readRDS(exp.rds)
  
  #Work with paired samples only
  if (paired==T){
    normal.patients<-substr(exp.matrix$normal.patients, 1, 12)
    
    #Get cancer patients with a matched normal
    cancer.patients<-data.table(samples=exp.matrix$cancer.patients, patients=substr(exp.matrix$cancer.patients, 1, 12))
    cancer.patients$MATCH<-cancer.patients$patients %in% normal.patients
    cancer.patients<-cancer.patients[MATCH==TRUE, ]$samples
    
    #Filter expression matrix for it
    exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, cancer.patients)]
    exp.matrix$cancer.patients<-cancer.patients
    
    #Return filtered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)
    
  } else {
    #Return unfiltered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)  
  }
  
}

internal.function<-function(n, exp.matrix, patients, normal.patients){
  
  #Obtain count of genes that are significantly differentiated between sample and normal cohorts
  target.patients<-sample(patients, n)
  p.vals<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[target.patients], y[normal.patients], paired=F)$p.value)
  p.vals.adj<-p.adjust(p.vals, method="fdr")
  
  return(sum(p.vals.adj<0.05))
}

Function.Main<-function(maf, exp.matrix){
  
  #Get set of total patients with expression and mutation data
  patients<-intersect(unique(maf$SAMPLE), colnames(exp.matrix$combined.matrices))
  
  #Filter both tables for these patients
  maf<-maf[SAMPLE %in% patients,]
  exp.matrix$cancer.patients<-intersect(exp.matrix$cancer.patients, patients)
  exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, exp.matrix$cancer.patients)]
  
  #Filter for samples that are seeing mutated in at least more than 2 samples at the same position
  maf[,N.SAMPLES.POS:=length(unique(SAMPLE)), by=c("Start_Position","Hugo_Symbol")]
  maf<-maf[N.SAMPLES.POS>=2,]
  
  #Simplify maf table
  maf<-maf[,c("Hugo_Symbol","SAMPLE", "Start_Position", "N.SAMPLES.POS"),with=F]
  setkey(maf)
  maf<-unique(maf)
  print (maf[order(N.SAMPLES.POS, decreasing=T),])
  
  #####################Create null distribution for empirical p-val calculation##########################
  print ("Creating null sample distribution")
  sample.dist<-unique(as.vector(maf$N.SAMPLES.POS))
  print (length(sample.dist))
  
  #Prepping first parallelization
  print ("prepping for first parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("sample.dist", "patients","as.data.table","exp.matrix","data.table", "internal.function") ,envir=environment())
  print ("Done exporting values")
  
  normal.patients<-exp.matrix$normal.patients
  
  print ("calculating null distributions")
  main.dist<-parLapply(cl, sample.dist, function(x) {
    sample.pop<-replicate(100, internal.function(x, exp.matrix, patients, normal.patients )) ###TESTING####
    return(sample.pop)
  })
  
  names(main.dist)<-as.character(sample.dist)
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with sample distributions")
  #############################################################################################
  
  #Split into tables 
  main.list<-split(maf, list(maf$Hugo_Symbol, maf$Start_Position),drop=T)
  
  #Prepping parallelization
  print ("prepping for second parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("main.list", "maf","as.data.table","exp.matrix","data.table") ,envir=environment())
  print ("Done exporting values")
  
  #Execute parallelization
  print ("Finding linkage")
  main.table<-parLapply(cl, main.list, function(x) {
    
    mut.gene<-unique(x$Hugo_Symbol)
    target.patients<-unique(x$SAMPLE)
    non.target.patients<-exp.matrix$normal.patients
    
    wilcox.matrix<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[target.patients], y[non.target.patients], paired=F)$p.value)
    wilcox.matrix<-data.table(Hugo_MUT=mut.gene, Position.MUT=unique(x$Start_Position) , Hugo_EXP=names(wilcox.matrix), 
                              N.PATIENTS=length(target.patients), P.VAL=as.vector(wilcox.matrix))  
    
    #Correct for multiple hypothesis 
    wilcox.matrix$P.VAL.ADJ<-p.adjust(wilcox.matrix$P.VAL, method="fdr")
    
    return(wilcox.matrix)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with linkages")
  
  #Merge gene tables
  main.table<-do.call(rbind, main.table)
  
  #Calculate per site empirical p-value based on empirical distribution
  print ("Calculating empirical p-values")
  main.table<-main.table[,list(EXP.P.VAL=   mean(main.dist[[as.character(unique(N.PATIENTS))]]  >= sum(P.VAL.ADJ<0.05)) ), 
                         by=c("Hugo_MUT","Position.MUT","N.PATIENTS")]
  
  print (main.table)
  
  #Correct for multiple hypothesis
  print ("Correcting for multiple hypothesis testing")
  main.table$EXP.P.VAL.ADJ<-p.adjust(main.table$EXP.P.VAL, method="fdr")
  
  #Cleaup and return
  main.table<-main.table[order(EXP.P.VAL.ADJ),]
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cancer.maf<-args[1]
exp.rds<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
maf<-Function.Prep.MAF(cancer.maf)
print ("done prepping maf")

exp.matrix<-Function.Prep.EXP(exp.rds, paired=F)
print ("done prepping expression rds")

main.function<-Function.Main(maf, exp.matrix)
print ("done")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)