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

Function.PREP.MAF.CLASS.ML<-function(MAF.ANNOVAR, REP.TIME, NT.LENGTH, EXTRA=c(), FILTER=F){
  #Added replication time and rep class
  
  require(data.table)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Filter out low counts?
  if (FILTER==T){
    MIN.MAF<-min(MAF.ANNOVAR$MAF[MAF.ANNOVAR$MAF!=0])
    MAF.ANNOVAR<-MAF.ANNOVAR[MAF>MIN.MAF,]
  }
  
  #Enter rep info
  MAF.ANNOVAR<-merge(MAF.ANNOVAR, REP.TIME, by="Hugo_Symbol")
  
  #Enter length info
  MAF.ANNOVAR<-merge(MAF.ANNOVAR, NT.LENGTH, by="Hugo_Symbol")
  MAF.ANNOVAR$NT<-normalize.vector(MAF.ANNOVAR$NT)
  
  #Classify MAF into categories
  MAF.ANNOVAR$MAF.CUT<-cut(MAF.ANNOVAR$MAF, quantile(MAF.ANNOVAR$MAF, c(0,0.33,0.66,1)), include.lowest=T, 
                           labels=as.factor(c("low", "medium", "high")))
  
  #Get cummulative MAF per hugo
  MAF.ANNOVAR[,HUGO.MAF.SUM:=sum(MAF), by="Hugo_Symbol"]
  
  #Filter out not needed columns and place classificaiton column at the end
  MAF.ANNOVAR<-MAF.ANNOVAR[,c("TYPE", "Chrom", "PHAST", "REF.ALT","REP.CLASS", "TRIMER", "MT", "MAF.CUT", "REP.TIME", "HUGO.MAF.SUM", "NT",
                              EXTRA), with=F]  
  
  #Return
  return(MAF.ANNOVAR)
}

args<-commandArgs(trailingOnly=T)

maf.annovar<-args[1] #Such as 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST
chen.rep<-args[2] #Such as 123114.BRCA.MAF.TRIMERS.csv
exon.nt<-args[3]
output.file<-args[4]

print ("Processing annovar")
MAF.ANNOVAR<-fread(maf.annovar, drop=c(1,6),sep="\t", header=F, stringsAsFactors=F)
MAF.ANNOVAR<-Function.PROCESS.ANNOVAR(MAF.ANNOVAR)

CHEN.REP<-fread(chen.rep, header=T, sep="\t",stringsAsFactors=F, drop=2)

EXON.NT<-fread(exon.nt, header=T, sep="\t", stringsAsFactors=F, drop=2:6)
EXON.NT<-EXON.NT[,list(NT=sum(FEAT_END-FEAT_START)), by=Hugo_Symbol]

######## Predicting with H2O on classifiable MAF variable#########################
print ("building annovar class")
MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR,  REP.TIME=CHEN.REP, NT.LENGTH=EXON.NT, EXTRA="MAF", FILTER=T)

# Start a local cluster with 30GB RAM
print ("setting up h2o")
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
key.v<-sample(letters,1)
h2o_maf <- as.h2o(localH2O, MAF.ANNOVAR.CLASS, key=key.v) 

#Break into train and test set
set.seed(1234)
ALL_ROWS<-1:nrow(h2o_maf)
RAND_FOLDS<-createFolds(ALL_ROWS,5)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:3]) #Train with 60%
TEST_ROWS<-unlist(RAND_FOLDS[4:5]) #Test on 40%
TRAIN_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TRAIN_ROWS]) 
TEST_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TEST_ROWS]) 

#Model without dropout
print ("building model")
MAF.MODEL<-h2o.deeplearning(x=c(1:7,11), y=8, data=h2o_maf[TRAIN_ROWS[1:100000],],
                            activation = "Tanh", balance_classes = TRUE, hidden = c(5000,1000,100), epochs = 500)
                            #input_dropout_ratio = 0.2, hidden_dropout_ratios = c(0.5,0.5, 0.5))
print ("done with dnn")
#Save model
h2o.saveModel(MAF.MODEL, "RESULTS", "MAF.MODEL.FILTERED.NT.REP.CLASS.100000.3.5000.1000.100.E.500.h2o")
print ("done saving")
