#Function.hugo.pfam
#Processes uniprot and pfams tables to provide pfam coordinates per hugo_symbol plus amino acid lengths per protein
#070814

args<-commandArgs(trailingOnly=T)

input.pfam.file<- args[1]
input.uniprot.file<-args[2]
output.data.file<-args[3]

require(data.table)

#Process files
test.p<-fread(input.pfam.file, sep="\t",stringsAsFactors=F,drop=2:4) 
test.u<-fread(input.uniprot.file, header=T, sep="\t", drop=c(2:4, 6))
setnames(test.u, colnames(test.u), c("UNIPROT", "GENES", "LENGTH"))
test.u$Hugo_Symbol<-sapply(test.u$GENES, function(x) strsplit(x," ")[[1]][1])
test.u$GENES<-NULL

#Merge to match uniprot with pfam to hugo_symbol
dummy.table<-as.data.table(merge(as.data.frame(test.p), as.data.frame(test.u), by.x="V1", by.y="UNIPROT"))
dummy.table<-dummy.table[complete.cases(dummy.table),] #remove NAs

#Clean up and return hugo_symbol pfams + gene start and end
dummy.table$START<-0
dummy.table<-dummy.table[,c(6,2:4,7,5), with=F]
dummy.table<-unique(dummy.table)
setnames(dummy.table, colnames(dummy.table), c("Hugo_Symbol", "PFAM", "FEAT.START", "FEAT.END", "START", "END"))

#Save
saveRDS(dummy.table, file=output.data.file)