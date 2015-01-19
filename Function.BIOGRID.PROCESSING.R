#####Function.BIOGRID.PROCESSING.R
#011615

library(data.table)

args<-commandArgs(trailingOnly=T) 

file.in<-args[1] #This of the form BIOGRID-ORGANISM-Homo_sapiens-3.2.120.tab2.txt 
output.file<-args[2]

#Load file
file.in<-"DATABASES/BIOGRID/BIOGRID-ORGANISM-3.2.120.tab2/BIOGRID-ORGANISM-Homo_sapiens-3.2.120.tab2.txt"
biogrid<-fread(file.in, header=T, sep="\t", stringsAsFactors=T, drop=c(1:7,12:24))
setnames(biogrid, c("A","B", "SYNA","SYNB"))
setkey(biogrid)
biogrid<-unique(biogrid)

#Work synonyms
biogrid<-biogrid[,list(Hugo.A=unique(c(A, unlist(strsplit(SYNA,"|",fixed=T))))),by=c("A","B","SYNA", "SYNB")]
biogrid<-biogrid[,list(Hugo.B=unique(c(B, unlist(strsplit(SYNB,"|",fixed=T))))),by=c("A","B","SYNA", "SYNB","Hugo.A")]

#Work on complementary to eliminate doubles
biogrid.comp<-biogrid[,c("B","A"),with=F]
setnames(biogrid.comp,c("A","B"))

biogrid<-unique(rbind(biogrid,biogrid.comp))

length(unique(intersect(TCGA.MUT.45.MID.2$Hugo_Symbol, biogrid$Hugo.A)))

#Get degree
biogrid.degree<-biogrid[,list(degree=length(B)),by="A"]
biogrid.degree<-biogrid.degree[order(degree, decreasing=T),]
setnames(biogrid.degree, c("Hugo_Symbol", "degree"))
biogrid.degree$DEGREE.CUT<-cut(biogrid.degree$degree, quantile(biogrid.degree$degree, c(0,0.25,0.50,0.75,1)), include.lowest=T)
