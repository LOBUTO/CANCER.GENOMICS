#Function.MAF.BMR.R
#011915
#Obtain background mutation rates of silent and missense mutations 

library(data.table)

args<-commandArgs(trailingOnly=T) 

maf.in<-args[1] #TCGA maf file 
exon.in<-args[2] #This is of the form 011915.EXON.COORDIANTES.FIXED.OVERLAPS
output.file<-args[3]

######FUNCTIONS#######
Function.Process.Exon<-function(exon.file){
  #Function to obtain nt length per hugo symbol
  
  exon<-fread(exon.file, sep="\t", header=T, stringsAsFactors=F)
  exon.nt<-exon[,list(nt.length=sum(FEAT_END-FEAT_START)), by="Hugo_Symbol"]
  
  return (exon.nt)
}

Function.Process.GL<-function(gene.length.file){
  
  gl<-fread(gene.length.file,sep="\t",header=T,stringsAsFactors=F)
  setnames(gl, c("Hugo_Symbol","nt.length"))
  setkey(gl)
  gl<-unique(gl)
  return(gl)
}

Function.Process.MAF<-function(maf.file){
  #Function to clean maf files straight from TCGA
  
  maf.record<-fread(maf.file, sep="\t", header=T, stringsAsFactors=F, drop=c(2:4,7:8, 10,12,14:15,17:37))
  maf.record<-maf.record[Variant_Classification %in% c("Missense_Mutation","Silent"),]
  maf.record<-maf.record[nchar(Tumor_Seq_Allele2)==1 & nchar(Reference_Allele)==1,]
  setkey(maf.record)
  maf.record<-unique(maf.record)
  
  #Return
  return(maf.record)
}

Function.Integrate.BMR<-function(maf.proc, exon.proc){
  
  #Assign length to each gene and obtain mutation rate at the patient level per gene and type
  
  #Total nt length - May have to take into account that there are inverse coding regions
  #total.nt<-sum(exon.proc$nt.length)
  total.nt<-sum(exon.proc[Hugo_Symbol %in% as.vector(maf.proc$Hugo_Symbol),]$nt.length) #smaller number
  
  #Integrate datasets
  main.table<-merge(maf.proc, exon.proc, by="Hugo_Symbol")
  
  #Obtain background mutation rate first per patient including both silent and missense
  main.table[,BMR:=length(Reference_Allele)/total.nt, by=Tumor_Sample_Barcode]
  
  #Filter out silent mutations 
  main.table<-main.table[Variant_Classification=="Missense_Mutation",]
  
  #Filter out genes that do not occur in at least 5 patients
  gene.count<-main.table[,list(gene.count=length(unique(Tumor_Sample_Barcode))), by="Hugo_Symbol"]
  gene.count<-gene.count[gene.count>=5,]
  main.table<-main.table[Hugo_Symbol %in% gene.count$Hugo_Symbol,]
  
  #Obtain background mutation rate first per patient
  #main.table[,BMR:=length(Reference_Allele)/total.nt, by=Tumor_Sample_Barcode]
  
  #Integrate background mutation rate
  #main.table<-merge(main.table, bmr.table, by="Tumor_Sample_Barcode", all.x=T)
  #main.table[is.na(main.table),]<-0
  
  #Then mutation rate per gene per patient
  main.table[,GMR:=length(Reference_Allele)/unique(nt.length),by=c("Tumor_Sample_Barcode","Hugo_Symbol")]
  
  #Clean table
  main.table<-main.table[,c("Hugo_Symbol","Tumor_Sample_Barcode","BMR","GMR"),with=F]
  setkey(main.table)
  main.table<-unique(main.table)
  
  #Return
  return(main.table)
}

#######EXECUTION######
#exon.proc<-Function.Process.Exon(exon.in)
exon.proc<-Function.Process.GL(exon.in)
maf.proc<-Function.Process.MAF(maf.in)
main.table<-Function.Integrate.BMR(maf.proc, exon.proc)

######OUTPUT#########
write.table(file=output.file, main.table, sep="\t", quote=F, row.names=F, col.names=T)
cat ("done pre-processing maf file")

#maf.file<-"DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
