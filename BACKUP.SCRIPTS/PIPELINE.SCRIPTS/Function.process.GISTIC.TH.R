#######FUNCTION TO PROCESS CNV THRESHOLDED GENES RESULT FROM GISTIC########
###080714####

#PROCESS ENTRIES
library(base)
library(data.table)
library(reshape2)

args<-commandArgs(trailingOnly=T)

input.threshold.file<-args[1] #e.g. "all_thresholded.by_genes.txt"
input.threshold<-args[2] #Option of 1 and 2
output.file<-args[3]

#EXECUTE
#Process file to melted table
BRCA.CNV.THRESHOLDED<-as.data.table(read.csv(input.threshold.file,header=T, sep="\t",stringsAsFactors=F))
BRCA.CNV.TABLE<-BRCA.CNV.THRESHOLDED[,c(1,4:ncol(BRCA.CNV.THRESHOLDED)),with=F]
BRCA.CNV.TABLE<-melt(BRCA.CNV.TABLE, id.vars="Gene.Symbol",variable.name="Sample")

#Remove genes depending on threshold
BRCA.CNV.TABLE<-BRCA.CNV.TABLE[abs(value)>=as.numeric(input.threshold),]

#Get per PATIENT sample
BRCA.CNV.TABLE$PATIENT<-substr(BRCA.CNV.TABLE$Sample,1,16)

#Clean up
BRCA.CNV.TABLE$Sample<-NULL
#BRCA.CNV.TABLE$value<-NULL
BRCA.CNV.TABLE<-BRCA.CNV.TABLE[,c("PATIENT", "Gene.Symbol","value"), with=F]
setnames(BRCA.CNV.TABLE, colnames(BRCA.CNV.TABLE), c("PATIENT", "Hugo_Symbol", "CNV.TH"))

#WRITE TO OUTPUT
saveRDS(BRCA.CNV.TABLE, file=output.file)