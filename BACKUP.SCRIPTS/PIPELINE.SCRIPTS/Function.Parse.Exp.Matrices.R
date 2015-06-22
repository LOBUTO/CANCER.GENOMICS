######Function.Parse.Exp.Matrices.R
#041715
#Constructs expression matrix using TCGA data, DOES NOT normalize data, just constructs matrix

library(data.table)
Function.process.RNAseq.map.files<-function(map.file, folder, rna.version) {  
  #Processes map.file from RNAseq data
  #Filters for gene name containing files and shortens barcode to represent patient.sample
  #Produces a 2-column matrix of "filename" and "patient.sample" columns
  
  dummy.map<-read.csv(map.file, header=T, sep="\t")
  
  #Keep only files that are for processed genes - The type of files depend on the rna processing version we are using
  if (rna.version=="V2"){
    dummy.map<-dummy.map[grepl("rsem.genes.normalized_results", dummy.map$filename),]  
  } else if ((rna.version=="V1") | (rna.version=="V1.raw")){
    dummy.map<-dummy.map[grepl("gene.quantification.txt", dummy.map$filename),]  
  } else if (rna.version=="V2.raw"){
    dummy.map<-dummy.map[grepl("rsem.genes.results", dummy.map$filename),]  
  } else if (rna.version=="agilent"){
    dummy.map<-dummy.map
  } else if (rna.version=="miRNASEQ"){
    dummy.map<-dummy.map[grepl("mirna.quantification.txt", dummy.map$filename),]  
  }
  
  if (rna.version=="miRNASEQ"){
    dummy.map$patient.sample<-substr(dummy.map$filename,1,16)
  } else{
    #Process sample name from barcode
    dummy.map$patient.sample<-substr(dummy.map$barcode.s.,1,16) 
  }
  dummy.map$barcode.s.<-NULL
  
  #Filter for files that actually exist in the rnaseq.folder 
  dummy.map<-dummy.map[dummy.map$filename %in% list.files(folder),]
  
  #Return
  return(as.matrix(dummy.map))
}

Function.read.RNAseq.files<- function(folder, processed.map.matrix, cancer.sep=T, rna.version) {
  #Produces matrices of genesxpatient.samples based on folder locations of RNAseq files and processed.map.matrix obtained from
  # Function.process.RNAseq.map.files()
  # Could produce 2 matrices corresponding to normal 11A and cancer 01A or simply a single matrix 
  #processed.map.matrix has first column as filename and second column as patient.sample
  
  Internal.Function.1<-function(input.matrix, folder){
    #Merges RNAseq files in "folder" depending on instructions from input matrix
    
    #Read files into list depending on rna version processing platform
    print ("extracting gene expression data from each file")
    if (rna.version=="V2"){
      matrices<-lapply(1:nrow(input.matrix), function(f) read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"), 
                                                                  stringsAsFactors=F, col.names=c("gene_id", input.matrix[f,2])))  
    } else if (rna.version=="V1"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"),
                            stringsAsFactors=F, col.names=c("gene_id", "u.1", "u.2", input.matrix[f,2]))
        pre.table<-pre.table[,c(1,4)]
        return(pre.table)
      })
    } else if (rna.version=="V2.raw"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"),
                            stringsAsFactors=F, col.names=c("gene_id", input.matrix[f,2] , "u.1", "u.2"))
        pre.table<-pre.table[,c(1,2)]
        return(pre.table) 
      })
    } else if(rna.version=="agilent"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"),
                            stringsAsFactors=F, col.names=c("gene_id", input.matrix[f,2]))
        pre.table<-pre.table[2:nrow(pre.table),c(1,2)]
        return(pre.table) 
      })
    } else if (rna.version=="V1.raw"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"),
                            stringsAsFactors=F, col.names=c("gene_id", input.matrix[f,2], "u.1", "u.2"))
        pre.table<-pre.table[,c(1,2)]
        return(pre.table)
      })
    } else if (rna.version=="miRNASEQ"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", na.strings=c("null"),
                            stringsAsFactors=F, col.names=c("gene_id", input.matrix[f,2], "u.1", "u.2"))
        pre.table<-pre.table[,c(1,2)]
        return(pre.table)
      })
    }
    
    print ("joining expression data into matrix")
    #Merge vector matrices to create expression table
    dummy.expression.table<-join_all(matrices, by="gene_id", type="inner")
    dummy.expression.table$gene_id<-as.character(dummy.expression.table$gene_id)
    
    #Minor fix to account for ANNOTATION ERROR
    dummy.expression.table[dummy.expression.table=="SLC35E2|728661"]<-"SLC35E2B|728661"
    
    #Remove "?" genes
    dummy.expression.table$gene<-sapply(dummy.expression.table$gene_id, function(x) strsplit(x, "[|]")[[1]][1])
    dummy.expression.table<-dummy.expression.table[dummy.expression.table$gene!="?",]
    
    #Remove duplicated genes
    dummy.expression.table<-dummy.expression.table[!duplicated(dummy.expression.table$gene),]
    
    #Assign genes as rownames and clean up
    rownames(dummy.expression.table)<-dummy.expression.table$gene
    dummy.expression.table$gene<-NULL
    dummy.expression.table$gene_id<-NULL
    dummy.expression.table<-dummy.expression.table[complete.cases(dummy.expression.table),]
    
    #Return
    return(dummy.expression.table)
  }
  
  require(plyr)
  require(data.table)
  
  if (cancer.sep==T) {
    processed.map.matrix<-as.data.frame(processed.map.matrix, stringAsFactors=F)
    
    #Separate based on normal or cancer - Will not include metastatic samples (06)
    processed.map.matrix$type<-substr(processed.map.matrix$patient.sample, 14,15)
    processed.map.matrix.normal<-processed.map.matrix[processed.map.matrix$type=="11",]
    processed.map.matrix.cancer<-processed.map.matrix[processed.map.matrix$type=="01",]
    
    #Throw error if don't have any samples of normal expression to compare to, BUT CONTINUE
    if (nrow(processed.map.matrix.normal)==0){
      print ("No normal samples found (Tissue - 11)!!!!")
      
      dummy.cancer<-Internal.Function.1(as.matrix(processed.map.matrix.cancer), folder)
      dummy.return<-list(tumor=dummy.cancer)
      
    } else {
      #If no error continue with normal processing
      dummy.normal<-Internal.Function.1(as.matrix(processed.map.matrix.normal), folder)
      dummy.cancer<-Internal.Function.1(as.matrix(processed.map.matrix.cancer), folder)
      
      dummy.return<-list(dummy.normal, dummy.cancer)
      names(dummy.return)<-c("normal", "tumor")  
    }
    
    
  } else if (cancer.sep==F) {
    dummy.return<-Internal.Function.1(processed.map.matrix, folder)
  } else
    print ("Please choose correct cancer.sep")
  
  #Return
  return(dummy.return)
}

#Arguments
args<-commandArgs(trailingOnly=T)
map.file<-args[1]
folder<-args[2]
rna.version<-args[3]
output.file<-args[4]
print("opened files")

#Execute
processed.map.matrix<-Function.process.RNAseq.map.files(map.file, folder, rna.version)

main.obj<-Function.read.RNAseq.files(folder, processed.map.matrix, cancer.sep=T, rna.version)

#Write out
saveRDS(main.obj, output.file)
print ("done writing to output")
