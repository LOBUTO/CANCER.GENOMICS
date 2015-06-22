#MAF PREDICTION TOOL
#010116 
#Based on trimer model with phastcon inclusion
#LACKS BASE ALGORITHM!!!

library(data.table)

args<-commandArgs(trailingOnly=T)

thousand.trimers.phast<-args[1] #Such as 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST
tcga.trimers.phast<-args[2] #Such as 123114.BRCA.MAF.TRIMERS.csv

thousand.clean<-function(thousand.table, filtering=0){
  
  #Load file
  MAF.TRIMERS.EXON.PHAST<-fread(thousand.table, sep="\t", stringsAsFactors=F)
  setnames(MAF.TRIMERS.EXON.PHAST, c("Chrom", "Position", "REF","ALT","TRIMER","MT","EXON","AF", "PHAST.100"))
  
  #Clean up
  MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[AF!=0,] 
  MAF.TRIMERS.EXON.PHAST$MAF<-ifelse(MAF.TRIMERS.EXON.PHAST$AF>0.5, 1-MAF.TRIMERS.EXON.PHAST$AF, MAF.TRIMERS.EXON.PHAST$AF)
  MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[!is.na(MT),]#clear of potential NXN trimers MT of NA!!!
  
  #Apply filtering if necessary (factor obtained from lowest count as of 2013 1000GENOME release)
  if (filtering!=0){
    filter=filtering*0.000199682
    MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[MAF>=filter,] #Filter smallest counts  
  }
  
  #Convert into factors
  MAF.TRIMERS.EXON.PHAST$Chrom<-as.factor(MAF.TRIMERS.EXON.PHAST$Chrom)
  MAF.TRIMERS.EXON.PHAST$TRIMER<-as.factor(MAF.TRIMERS.EXON.PHAST$TRIMER)
  
  return (MAF.TRIMERS.EXON.PHAST)

}

tcga.clean<-function(tcga.table){
  tcga<-fread(tcga.table, sep="\t", stringsAsFactors=F, header=F)
  setnames(tcga, c("Chrom", "Position", "Hugo_Symbol", "MUT.FREQ", "TYPE", "REF", "ALT", "TRIMER", "MT", "PHAST.100"))
  
}

maf.model<-function(thousand.cleaned){
  
  #Consider exon only instances since we are predicting exome data. NOTE: COULD CHANGE IN FUTURE!!!
  thousand.cleand<-thousand.cleaned[EXON==TRUE,]
  
  #Build linear models
  Chromosomes<-as.character(1:22,"X","Y")
  for (chrm in Chromosomes) {
    model<-lm(formula=MAF~exp(PHAST.100)+MT+TRIMER, data=thousand.cleaned[Chrom=="chrm",])  
    assign(paste0("mod.", chrm), )
    
  }
  
}



output.file<-
