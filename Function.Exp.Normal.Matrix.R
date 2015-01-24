#Function.Exp.Normal.Matrix.R
#102514
#Function to normalize RNA seq expression files obtained from map files and quantile normalized rsem.genes.normalized_results

Function.process.RNAseq.map.files<-function(map.file, folder, rna.version) {  
  #Processes map.file from RNAseq data
  #Filters for gene name containing files and shortens barcode to represent patient.sample
  #Produces a 2-column matrix of "filename" and "patient.sample" columns
  
  dummy.map<-read.csv(map.file, header=T, sep="\t")
 
  #Keep only files that are for processed genes - The type of files depend on the rna processing version we are using
  if (rna.version=="V2"){
    dummy.map<-dummy.map[grepl("rsem.genes.normalized_results", dummy.map$filename),]  
  } else if (rna.version=="V1"){
    dummy.map<-dummy.map[grepl("gene.quantification.txt", dummy.map$filename),]  
  }
  
  #Process sample name from barcode
  dummy.map$patient.sample<-substr(dummy.map$barcode.s.,1,16)
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
    if (rna.version=="V2"){
      matrices<-lapply(1:nrow(input.matrix), function(f) read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t",
                                                                  col.names=c("gene_id", input.matrix[f,2])))  
    } else if (rna.version=="V1"){
      matrices<-lapply(1:nrow(input.matrix), function(f) {
        pre.table<-read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t", 
                 col.names=c("gene_id", "u.1", "u.2", input.matrix[f,2]))
        pre.table<-pre.table[,c(1,4)]
        return(pre.table)
      })
    }
    
    #Merge vector matrices to create expression table
    dummy.expression.table<-join_all(matrices, by="gene_id", type="inner")
    dummy.expression.table$gene_id<-as.character(dummy.expression.table$gene_id)
    
    #Minor fix to account for ANNOTATION ERROR
    dummy.expression.table[dummy.expression.table=="SLC35E2|728661"]<-"SLC35E2B|728661"
    
    #Remove "?" genes
    dummy.expression.table$gene<-sapply(dummy.expression.table$gene_id, function(x) strsplit(x, "[|]")[[1]][1])
    dummy.expression.table<-dummy.expression.table[dummy.expression.table$gene!="?",]
    
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
    print (processed.map.matrix)
    
    #Separate based on normal or cancer
    processed.map.matrix$type<-substr(processed.map.matrix$patient.sample, 14,15)
    processed.map.matrix.normal<-processed.map.matrix[processed.map.matrix$type=="11",]
    processed.map.matrix.cancer<-processed.map.matrix[processed.map.matrix$type=="01",]
    
    dummy.normal<-Internal.Function.1(as.matrix(processed.map.matrix.normal), folder)
    dummy.cancer<-Internal.Function.1(as.matrix(processed.map.matrix.cancer), folder)
    
    dummy.return<-list(dummy.normal, dummy.cancer)
    names(dummy.return)<-c("normal", "tumor")
    
  } else if (cancer.sep==F) {
    dummy.return<-Internal.Function.1(processed.map.matrix, folder)
  } else
    print ("Please choose correct cancer.sep")
  
  #Return
  return(dummy.return)
}

Function.RNAseq.Matrices.Normalization<-function(normal.matrix, cancer.matrix, rm.batch.effect=T) {
  #Normalize RNAseq matrices for posterior differential expression analysis and returns fully normalized expression matrix 
  #+ list of normal and cancer patients
  #batch effect normalization may not necessarily advantageous and end up being more time consuming, apply with care!
  
  require(plyr)
  require(limma)
  require(edgeR)
  require(sva)
  
  #Get patient cohorts
  G1.patients<-colnames(cancer.matrix)
  G0.patients<-colnames(normal.matrix)
  
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 3 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  
  if (rm.batch.effect==T) {
    print("performing batch effect correction")

    #Obtain batch normalization parameters
    batch.design=data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
    batch.mod = model.matrix(~ G, batch.design) #What we have
    batch.mod0 = model.matrix(~ 1, batch.design) #Null
    batch.ss = svaseq(G.all$counts, batch.mod, batch.mod0) #Apply function to get batch effect parameters
    print("done with batch effect correction")
  }
  
  #Quantile normalization
  print("converting log log-count per million")
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Need to correct for batch effects?
  if (rm.batch.effect==T){
    #Obtain matrices corrected from batch effect
    f.ss = lmFit(G.all, model.matrix(~G + batch.ss$sv, batch.design)) #Just like regular fit but added batch effect normalization
    ss.removed = G.all$E - f.ss$coefficients[, -c(1, 2)] %*% t(f.ss$design[, -c(1, 2)])
  } else {
    ss.removed=G.all$E
  }
  
  #Return combined matrix fully normalized (quantile, log cpm converated + batch effect) GENESxPATIENTS + list of cancer and normal patients
  ss.removed<-list(combined.matrices=ss.removed, cancer.patients=G1.patients, normal.patients=G0.patients)
  return(ss.removed)
}

args<-commandArgs(trailingOnly=T)

map.file<-args[1] 
rnaseq.folder<-args[2]
rna.version<-args[3]
output.file<-args[4]

#map.file<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/102514/FILE_SAMPLE_MAP.txt"
#rnaseq.folder<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/102514/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3"
#output.file<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/102514.CANCER.MATRICES.NORMALIZED.OBJ.rds"

processed.map.matrix<-Function.process.RNAseq.map.files(map.file, rnaseq.folder, rna.version)
print ("map file constructed\n")

print ("Building main expression matrix")
cancer.matrices<-Function.read.RNAseq.files(rnaseq.folder, processed.map.matrix, cancer.sep=T, rna.version)
print ("main table built")

print ("Normalizing expression values")
output.obj<-Function.RNAseq.Matrices.Normalization(cancer.matrices$normal, cancer.matrices$tumor, rm.batch.effect=F)
print ("done normalizing")

#saveRDS(object=ss.removed, file=output.file)
saveRDS(object=output.obj, file=output.file)
print ("done saving")