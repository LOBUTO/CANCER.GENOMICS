######Function.RNASEQ.CANCERONLY.R#####
#031315
#Function to process normalized matrix form text files directly downloaded from the broad institute repository:
#   Such as file GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt in the
#   site http://gdac.broadinstitute.org/runs/stddata__latest/data/GBM/20150204/ (for RNASEQ) - OR -
#   file GBM.transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt (for MA)
#Designed to filter for cancer samples only

#NOTE: Will also normalize MA data depending on the argument given ("RNASEQ" or "MA")

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(edgeR)
library(limma)

Function.Main<-function(rnaseq.file, type, rnaseq.file.2=c()){
  
  #Load file
  RNASEQ<-fread(rnaseq.file, header=T, sep="\t",stringsAsFactors=F)
  
  #Parse to obtain correct gene names and remove duplicates
  RNASEQ<-RNASEQ[2:nrow(RNASEQ),]
  setnames(RNASEQ, c("Hugo_Symbol", colnames(RNASEQ)[2:ncol(RNASEQ)]))
  RNASEQ$Hugo_Symbol<-sapply(RNASEQ$Hugo_Symbol, function(x)  unlist(strsplit(x, "[|]"))[1])
  RNASEQ<-RNASEQ[Hugo_Symbol!="?",]
  RNASEQ<-RNASEQ[!duplicated(RNASEQ$Hugo_Symbol),]
  
  #Add additional file if present
  if (length(rnaseq.file.2)>0){
    #Load file
    RNASEQ.2<-fread(rnaseq.file.2, header=T, sep="\t",stringsAsFactors=F)
    
    #Parse to obtain correct gene names and remove duplicates
    RNASEQ.2<-RNASEQ.2[2:nrow(RNASEQ.2),]
    setnames(RNASEQ.2, c("Hugo_Symbol", colnames(RNASEQ.2)[2:ncol(RNASEQ.2)]))
    RNASEQ.2$Hugo_Symbol<-sapply(RNASEQ.2$Hugo_Symbol, function(x)  unlist(strsplit(x, "[|]"))[1])
    RNASEQ.2<-RNASEQ.2[Hugo_Symbol!="?",]
    RNASEQ.2<-RNASEQ.2[!duplicated(RNASEQ.2$Hugo_Symbol),]
    
    #Merge with original table
    RNASEQ<-merge(RNASEQ, RNASEQ.2, by="Hugo_Symbol")
  }
  
  HUGO.NAMES<-copy(RNASEQ$Hugo_Symbol)
  RNASEQ$Hugo_Symbol<-NULL
  
  #Parse patient ID and remove duplicate samples
  RNASEQ<-as.data.frame(RNASEQ)
  setnames(RNASEQ, sapply(colnames(RNASEQ), function(x) paste(unlist(strsplit(x, "-"))[1:4], collapse=".")))   
  RNASEQ<-RNASEQ[,!duplicated(colnames(RNASEQ))]
  
  #Convert to matrix while keeping gene names
  RNASEQ<-as.matrix(RNASEQ)
  RNASEQ<-apply(RNASEQ, 2, function(x) as.numeric(x))
  rownames(RNASEQ)<-HUGO.NAMES
  
  #Remove NA
  RNASEQ<-RNASEQ[complete.cases(RNASEQ),]
  
  #NOTE: Filter for cancer samples only
  cancer.patients<-data.table(SAMPLE=colnames(RNASEQ), TYPE=substr(colnames(RNASEQ), 14, 16))
  cancer.patients<-cancer.patients[TYPE %in% c("01A", "01B","01C"),]
  RNASEQ<-RNASEQ[,cancer.patients$SAMPLE]
  
  #Choose normalization scheme depending on type of data (RNASEQ or MA)
  if (type=="RNASEQ") {
    
    #Normalize and filter low counts  
    EXPR<-rowSums(cpm(RNASEQ)>1) >= 3
    RNASEQ<-RNASEQ[EXPR,]
    RNASEQ<-DGEList(counts=RNASEQ)
    RNASEQ<-calcNormFactors(RNASEQ)
    
    #Convert to counts per million
    VOOM<-voom(RNASEQ, normalize.method="quantile")
    
    #Return
    return(list(combined.matrices=VOOM$E))
    
  } else if(type=="MA"){
    
    #Normalize arrays
    RNASEQ<-normalizeBetweenArrays(RNASEQ, method="quantile")
    
    #Return
    return(list(combined.matrices=RNASEQ))
  }
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
rnaseq.file<-args[1]
type<-args[2]
output.file<-args[3]
if (length(args)>3){
  rnaseq.file.2<-args[4]
} else {
  rnaseq.file.2<-c()
}

print("opened files")
##########################################

##################EXECUTE#################
main.function<-Function.Main(rnaseq.file, type, rnaseq.file.2)
print ("Done with main function")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing to file")