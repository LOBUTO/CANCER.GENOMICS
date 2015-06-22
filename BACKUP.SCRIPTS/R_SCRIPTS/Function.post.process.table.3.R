
#options(echo=T) #To see commands in output file
args<-commandArgs(trailingOnly=T)

input.table.3.filename<-args[1]
output.file<-args[2]

Function.post.process.table.3<-function(table.3) {
  #Small fine tuning
  #   unique() to filter for duplicates obtained from 
  #   Filters out "Product" column from table.3
  #   Filters out water and H2O
  #Table.3 of the form BRCA.Table.3 (data.table with 3 columns)
  
  require(data.table)
  
  #Open table 3 file
  table.3<-as.data.table(read.csv(table.3, header=T, sep="\t"))
  
  #Process
  table.3<-table.3[Product!="H2O",]
  table.3<-table.3[Product!="Water",]
  table.3<-unique(table.3[,c(1,2),with=F])
  
  #Return
  return (table.3)
}

#Call
Processed.table.3<-Function.post.process.table.3(input.table.3.filename)
saveRDS(Processed.table.3, file=output.file)

