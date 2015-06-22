######Function.RNASEQ.DESEQ.R######
#Executes Differential expression analysis on raw count matrix using DESeq2
#041815

####FUNCTIONS#####
library(DESeq2)
library(BiocParallel)

Function.MAIN<-function(rnaseq.obj){
  
  brca.V1.RAW<-readRDS(rnaseq.obj)

  combined.matrices<-cbind(brca.V1.RAW$normal, brca.V1.RAW$tumor[rownames(brca.V1.RAW$normal),])
  combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
  combined.matrices<-data.matrix(combined.matrices)
  
  col.data<-data.frame(type=c(rep("normal", ncol(brca.V1.RAW$normal)), rep("cancer", ncol(brca.V1.RAW$tumor))) )
  rownames(col.data)<-c(colnames(brca.V1.RAW$normal), colnames(brca.V1.RAW$tumor))
  
  dds<-DESeqDataSetFromMatrix(countData=combined.matrices, colData=col.data, design=~type)
  
  dds$type<-relevel(dds$type, "normal")
  dds$type<-droplevels(dds$type)
  
  register(MulticoreParam(16))
  dds<-DESeq(dds, parallel=T)
  res<-results(dds)
  
  res<-res[order(res$padj),]
  return(res)
}

#################

#Arguments
args<-commandArgs(trailingOnly=T)
rnaseq.obj<-args[1]
output.file<-args[2]
print("opened files")

#Execute
main.obj<-Function.MAIN(rnaseq.obj)

#Write out
saveRDS(main.obj, output.file)
print ("done writing file")