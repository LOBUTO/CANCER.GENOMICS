args<-commandArgs(trailingOnly=T)
output.file<-args[2]

input.pre.table.2.file<-args[1]
input.table.3.data<-args[2]
output.data.file<-args[3]

Function.post.process.table.2<-function(table.2, KEGG.IDS) {
  #Filters table.2 for KEGG.IDS, preferrably from table.3
  #Table.2 of the form BRCA.table.2 (data.table with KEGG_ID column)
  
  require (data.table)
  table.2<-table.2[KEGG_ID %in% KEGG.IDS,]
  return(table.2)
}

require(data.table)

#Load table.2
input.table.2<-as.data.table(read.csv(input.pre.table.2.file, header=T, sep="\t"))

#Get KEGG identifiers from table.3
table.3<-readRDS(input.table.3.data)
KEGGS<-unique(as.vector(table.3$KEGG_ID))

#Call function
Processed.table.2<-Function.post.process.table.2(input.table.2, KEGGS)

#Write
saveRDS(Processed.table.2, file=output.data.file)

