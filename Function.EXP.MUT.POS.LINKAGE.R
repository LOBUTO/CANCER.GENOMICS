#Function.EXP.MUT.POS.LINKAGE.R
#021215
#Calculates the linkage between mutation position and expression
#Will only take into account sites that are seen mutated at least 2 times in population

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

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
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Return
  return(maf)
}

Function.Prep.EXP<-function(exp.rds){
  
  exp.matrix<-readRDS(exp.rds)
  
  return(exp.matrix$combined.matrices)
}

Function.Main<-function(maf, exp.matrix){
  
  #Filter for samples that are seeing mutated in at least more than the same position
  maf[,N.SAMPLES.POS:=length(unique(SAMPLE)), by=c("Start_Position","Hugo_Symbol")]
  maf<-maf[N.SAMPLES.POS>=2,]
  
  #Get set of total patients with expression and mutation data
  patients<-intersect(unique(maf$SAMPLE), colnames(exp.matrix))
  
  #Filter both tables for these patients
  maf<-maf[SAMPLE %in% patients,]
  exp.matrix<-exp.matrix[,patients]
  
  #Simplify maf table
  maf<-maf[,c("Hugo_Symbol","SAMPLE", "Start_Position"),with=F]
  setkey(maf)
  maf<-unique(maf)
  
  #Split into tables 
  print ("Spliting maf")
  main.list<-split(maf, list(maf$Hugo_Symbol, maf$Start_Position))
  
  #Obtain mutated gene list
  #gene.list<-unique(maf$Hugo_Symbol)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("main.list", "maf","as.data.table","exp.matrix","data.table") ,envir=environment())
  print ("Done exporting values")
  
  #Execute parallelization
  print ("Finding linkage")
  main.table<-parLapply(cl, main.list, function(x) {
    
    #To perform wilcoxon need to have at least 2 samples with mutations
    target.patients<-unique(x$SAMPLE)
    
    if (length(target.patients)>1){
      non.target.patients<-setdiff(patients, target.patients)
      
      wilcox.matrix<-apply(exp.matrix, 1, function(y) wilcox.test(y[target.patients], y[non.target.patients], paired=F)$p.value)
      wilcox.matrix<-data.table(Hugo_MUT=x, Position.MUT=unique(x$Start_Position) , Hugo_EXP=names(wilcox.matrix), 
                                N.PATIENTS=length(target.patients), P.VAL=as.vector(wilcox.matrix))  
      
    } else {
      wilcox.matrix<-data.table(Hugo_MUT=x, Position.MUT=unique(x$Start_Position), Hugo_EXP=rownames(exp.matrix), 
                                N.PATIENTS=length(target.patients), P.VAL="NONE")  
    }
    
    return(wilcox.matrix)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with linkages")
  
  #Merge gene tables
  main.table<-do.call(rbind, main.table)
  
  #Correct for multiple hypothesis
  print ("Correcting for multiple hypothesis testing")
  main.table<-main.table[P.VAL!="NONE",]
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Cleaup and return
  main.table<-main.table[order(P.VAL.ADJ),]
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

exp.matrix<-Function.Prep.EXP(exp.rds)
print ("done prepping expression rds")

main.function<-Function.Main(maf, exp.matrix)
print ("done")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)