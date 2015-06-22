#######Function.MAF.METABOLITE.MR.R########
#012015
#Integrate maf information and find mutation rate per metabolite per cancer maf
#Called by MET.MR.CANCERS.ENRICH.sh

library(data.table)

######LOAD ARGUMENTS######
args<-commandArgs(trailingOnly=T) 

maf.in<-args[1] #TCGA maf file 
exon.in<-args[2] #This is of the form 011915.EXON.COORDIANTES.FIXED.OVERLAPS
table.2<-readRDS(args[3]) #This is of the form 061214_Table.2.rds
output.file<-args[4] #Output file

######FUNCTIONS#######
Function.Process.Exon<-function(exon.file){
  #Function to obtain nt length per hugo symbol
  
  exon<-fread(exon.file, sep="\t", header=T, stringsAsFactors=F)
  exon.nt<-exon[,list(nt.length=sum(FEAT_END-FEAT_START)), by="Hugo_Symbol"]
  
  return (exon.nt)
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

Function.MR.Discovery<-function(maf.proc, exon.proc){
  #Function to obtain sum of mutation rate across patients per gene
  
  #Obtain number of mutations per gene in cancer
  main.table<-maf.proc[,list(N.MUT=length(Reference_Allele)), by=c("Hugo_Symbol","Variant_Classification")]
  
  #Obtain mut freq by dividing over length
  main.table<-merge(main.table, exon.proc, by="Hugo_Symbol")
  main.table$HUGO.MUT.FREQ<-main.table$N.MUT/main.table$nt.length
  
  #Clean up and return [Hugo_Symbol, Variant_Classification, HUGO.MUT.FREQ]
  main.table<-main.table[,c("Hugo_Symbol","Variant_Classification","HUGO.MUT.FREQ"),with=F]
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.Meta.Enrichment<-function(main.table, table.2){
  #Obtain enrichment score (p.val) per metabolite 
  
  internal.function<-function(n.hugo, met.cum.mr){
    ran.test<-replicate(1000, sum(sample(mr.pool, n.hugo)))
    ran.test<-mean(ran.test>=met.cum.mr)
    return(list(P.VAL=ran.test))
  }
  
  #Deal with missense first
  main.table<-main.table[Variant_Classification=="Missense_Mutation",]
  
  #Integrate metabolite information [Hugo_Symbol, HUGO.MUT.FREQ, METABOLITE, KEGG_ID]
  setnames(table.2, c("METABOLITE","KEGG_ID","Hugo_Symbol"))
  main.table<-merge(main.table, table.2, by="Hugo_Symbol")
  
  #Obtain MR per metabolite - Investement of cell in mutations associated with metabolite
  met.mr<-main.table[,list(MET.CUM.MR=sum(HUGO.MUT.FREQ), N.HUGO=length(unique(Hugo_Symbol))),by=c("KEGG_ID","METABOLITE")]
  
  #Perform enrichment test - Drawing as many number of hugo, do we get equal or great MET.CUM.MR by chance
  mr.pool<-as.vector(main.table$HUGO.MUT.FREQ)
  met.test<-met.mr[,internal.function(N.HUGO,MET.CUM.MR),by=c("KEGG_ID", "METABOLITE")]
  
  #Correct for multiple hypothesis testing
  met.test$P.VAL.ADJ<-p.adjust(met.test$P.VAL, method="fdr")
  
  #Clean up and return
  met.test<-met.test[order(P.VAL.ADJ),]
  return(met.test)
}

#maf.file<-"DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"

######EXECTUTE#######
exon.proc<-Function.Process.Exon(exon.in)
cat ("Done processing exon\n")

maf.proc<-Function.Process.MAF(maf.in)
cat ("Done processing maf file\n")

mr.discovery<-Function.MR.Discovery(maf.proc, exon.proc)
cat ("Done calculating mutation rates\n")

cat ("Running tests\n")
command.run<-Function.Meta.Enrichment(mr.discovery, table.2)

######OUTPUT#########
cat("Writing to output\n")
write.table(file=output.file, command.run, sep="\t", quote=F, row.names=F, col.names=T)
cat("Done")

