library(data.table)
library(reshape2)
library(h2o)
library(mlbench)
library(caret)

Function.PROCESS.ANNOVAR<-function(MAF.ANNOVAR){
  #Assuming columns 1 and 6 have been dropped from avinput_exonic_variant_function file
  #Assumes files of type 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST.avinput.exonic_variant_function as input
  
  require(Biostrings)
  require(car)
  
  #Parse into hugo_symbol
  MAF.ANNOVAR$V3<-sapply(MAF.ANNOVAR$V3, function(x) strsplit(x,":")[[1]][1]) 
  
  #Assign labels
  setnames(MAF.ANNOVAR,  c("TYPE", "Hugo_Symbol","Chrom", "Position", "REF","ALT","TRIMER","MT","EXON","AF", "PHAST"))
  
  #Filter out unknown types
  MAF.ANNOVAR<-MAF.ANNOVAR[TYPE!="unknown",]
  
  #Get MAF from AF
  MAF.ANNOVAR$MAF<-ifelse(MAF.ANNOVAR$AF>0.5, 1-MAF.ANNOVAR$AF, MAF.ANNOVAR$AF)
  
  #Clean up
  MAF.ANNOVAR<-MAF.ANNOVAR[!is.na(MT),]
  MAF.ANNOVAR<-MAF.ANNOVAR[EXON==TRUE,]
  
  #Equalize REF-ALT pairs
  MAF.ANNOVAR$PRE.REF.ALT<-paste(as.vector(MAF.ANNOVAR$REF), as.vector(MAF.ANNOVAR$ALT), sep="_")
  MAF.ANNOVAR$REF.ALT<-recode(MAF.ANNOVAR$PRE.REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  
  #Equalize trimers
  MAF.ANNOVAR$RECODE<-MAF.ANNOVAR$PRE.REF.ALT!=MAF.ANNOVAR$REF.ALT
  TEMP.RECODE<-MAF.ANNOVAR[RECODE==T,]
  TEMP.RECODE$TRIMER<-substring(complement(DNAStringSet(TEMP.RECODE$TRIMER)),1,3)
  MAF.ANNOVAR<-rbind(MAF.ANNOVAR[RECODE==F,],TEMP.RECODE)
  
  #Clean up and Return
  MAF.ANNOVAR$PRE.REF.ALT<-NULL
  MAF.ANNOVAR$RECODE<-NULL
  return(MAF.ANNOVAR)
}

Function.PREP.MAF.CLASS.ML<-function(MAF.ANNOVAR, REP.TIME, HUGO=F){
  #Added replication time and rep class
  
  require(data.table)
  
  #Enter rep info
  MAF.ANNOVAR<-merge(MAF.ANNOVAR, REP.TIME, by="Hugo_Symbol")
  
  #Classify MAF into categories
  MAF.ANNOVAR$MAF.CUT<-cut(MAF.ANNOVAR$MAF, quantile(MAF.ANNOVAR$MAF, c(0,0.33,0.66,1)), include.lowest=T, 
                           labels=as.factor(c("low", "medium", "high")))
  
  #Filter out not needed columns and place classificaiton column at the end
  if (HUGO==F){
    MAF.ANNOVAR<-MAF.ANNOVAR[,c("TYPE", "Chrom", "PHAST", "REF.ALT","REP.CLASS", "TRIMER", "MT", "MAF.CUT", "REP.TIME"), with=F]  
  } else{
    MAF.ANNOVAR<-MAF.ANNOVAR[,c("TYPE", "Hugo_Symbol", "Chrom", "PHAST", "REF.ALT","REP.CLASS", "TRIMER", "MT", "MAF.CUT", "REP.TIME"), with=F]  
  }
  
  #Return
  return(MAF.ANNOVAR)
}

args<-commandArgs(trailingOnly=T)

maf.annovar<-args[1] #Such as 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST
chen.rep<-args[2] #Such as 123114.BRCA.MAF.TRIMERS.csv
output.file<-args[3]

print ("Processing annovar")
MAF.ANNOVAR<-fread(maf.annovar, drop=c(1,6),sep="\t", header=F, stringsAsFactors=F)
MAF.ANNOVAR<-Function.PROCESS.ANNOVAR(MAF.ANNOVAR)

CHEN.REP<-fread(chen.rep, header=T, sep="\t",stringsAsFactors=F, drop=2)

MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR,  CHEN.REP)

######## Predicting with H2O on classifiable MAF variable#########################
print ("building annovar class")
MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR,  CHEN.REP)

# Start a local cluster with 8GB RAM
print ("setting up h2o")
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
h2o_maf <- as.h2o(localH2O, MAF.ANNOVAR.CLASS) 

#Break into train and test set
set.seed(1234)
ALL_ROWS<-1:nrow(h2o_maf)
RAND_FOLDS<-createFolds(ALL_ROWS,3)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:2])
TEST_ROWS<-unlist(RAND_FOLDS[3])
TRAIN_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TRAIN_ROWS]) #Train with 2/3
TEST_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TEST_ROWS]) #Test on 1/3

#Model without dropout
print ("building model")
MAF.MODEL<-h2o.deeplearning(x=c(1:4,6:7,9), y=8, data=h2o_maf[TRAIN_ROWS[1:10000],],
                            activation = "Tanh", balance_classes = TRUE, hidden = c(1000,500,100), epochs = 400)
print ("done with dnn")
#Save model
h2o.saveModel(MAF.MODEL, "RESULTS", "MAF.MODEL.10000.3.1000.500.100.E.400.h2o")
print ("done saving")
