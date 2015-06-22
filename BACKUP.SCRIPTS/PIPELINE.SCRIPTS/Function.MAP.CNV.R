#####Function.MAP.CNV.R
#021615
#Function to map cnv file regions to genes
#PART OF: TCGA.CNV.MATRIX.CONSTRUCT.sh

##################FUNCTIONS################
library(data.table)
library(reshape2)

Function.main<-function(cnv.folder) {
  
  #Filter for reduced noise files
  files<-list.files(cnv.folder, pattern="nocnv_hg19")
  
  #Combine files
  main.list<-lapply(files, function(x) {
    cnv.file<-paste0(cnv.folder, "/",x)
    cnv.file<-fread(cnv.file, sep="\t",header=T, stringsAsFactors=F, drop=c(1,5,6))
    return(cnv.file)
  })
  
  main.file<-do.call(rbind, main.list)
  
  #Filter for uniqueness of regions
  setkey(main.file)
  main.file<-unique(main.file)
  
  #Add ref and alt dummy columns needed for annovar
  main.file$REF<-0
  main.file$ALT<-"-"
  
  #Return
  return(main.file)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cnv.folder<-args[1]
output.file<-args[2]
print("Files loaded")
##########################################

##################EXECUTE#################
main.function<-Function.main(cnv.folder)
print ("done")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=F)