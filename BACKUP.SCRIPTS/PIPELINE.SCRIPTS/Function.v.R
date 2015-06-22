#Function.v.R
#071714
#Calculates the differential expression influence of metabolite j
#It compares cancer(or disease) phenotype vs normal phenotype and gives a ratio of the number of differentially expressed genes
#Returns a table per metabolite with its respective v(j) value

#####################################################################################################################################################################

Function.read.RNAseq.files<- function(folder, processed.map.matrix, cancer.sep=T) {
  #Produces matrices of genesxpatient.samples based on folder locations of RNAseq files and processed.map.matrix obtained from
  # Function.process.RNAseq.map.files()
  # Could produce 2 matrices corresponding to normal 11A and cancer 01A or simply a single matrix 
  #processed.map.matrix has first column as filename and second column as patient.sample
  
  Internal.Function.1<-function(input.matrix, folder){
    #Merges RNAseq files in "folder" depending on instructions from input matrix
    
    #Read files into list
    matrices<-lapply(1:nrow(input.matrix), function(f) read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t",
                                                                col.names=c("gene_id", input.matrix[f,2])))
    
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
  
  if (cancer.sep==T) {
    processed.map.matrix<-as.data.frame(processed.map.matrix, stringAsFactors=F)
    
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

Function.process.RNAseq.map.files<-function(map.file) {  
  #Processes map.file from RNAseq data
  #Filters for gene name containing files and shortens barcode to represent patient.sample
  #Produces a 2-column matrix of "filename" and "patient.sample" columns
  
  dummy.map<-read.csv(map.file, header=T, sep="\t")
  
  #Keep only files that are for processed genes
  dummy.map<-dummy.map[grepl("rsem.genes.normalized.results", dummy.map$filename),]
  
  #Process sample name from barcode
  dummy.map$patient.sample<-substr(dummy.map$barcode.s.,1,16)
  dummy.map$barcode.s.<-NULL
  
  #Return
  return(as.matrix(dummy.map))
}

Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  
  require(plyr)
  require(limma)
  require(edgeR)
  
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
  G.isexpr<- rowSums(cpm(G.all)>1) >= 20 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data 
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Return topTable
  return(all.G.fit)
}

Function.v<-function(RNA.SEQ.FOLDER, cancer.map.file, processed.table.1.rds, processed.table.3.rds, Table.u.KEGG_IDs) {
  #Calculates v(j) based on degree of influence of KEGG(j) mutations on cancer phenotype when compared to normal
  #Need functions pre-loaded:
  #   Function.read.RNAseq.files()
  #   Function.process.RNAseq.map.files()
  require(base)
  require(data.table)
  
  #Load RDS files
  dummy.table.1<-readRDS(processed.table.1.rds)
  dummy.table.1<-dummy.table.1$table.1[,1:2, with=F]
  dummy.table.3<-readRDS(processed.table.3.rds)
  
  #Create tumor and cancer matrices using map file on RNA-seq folder
  RNASEQ.MATRICES<-Function.read.RNAseq.files(RNA.SEQ.FOLDER, cancer.sep=T,
                                              Function.process.RNAseq.map.files(cancer.map.file))
  
  #Get cancer samples in processed.table.1 that have expression information in cancer matrix [Hugo_Symbol, PATIENT]
  cancer.rnaseq.patients<-colnames(RNASEQ.MATRICES$tumor)
  dummy.table.1$PATIENT<-sapply(as.character(dummy.table.1$Tumor_Sample_Barcode), function(x) paste(strsplit(x, "-")[[1]][1:4] , collapse="."))
  dummy.table.1<-dummy.table.1[PATIENT %in% cancer.rnaseq.patients, ]
  dummy.table.1$Tumor_Sample_Barcode<-NULL
  
  #Prep storeage matrix for v(j)
  dummy.table.v<-data.frame(KEGG_ID=c(), v.METABOLITE=c())
  
  #To keep count
  x=0
  
  #Process each metabolite for v(j)
  for (metabolite in Table.u.KEGG_IDs) {
    
    #Get enzymes that produce metabolite
    dummy.metabolite.enzymes<-unique(as.vector(dummy.table.3[KEGG_ID==metabolite,]$Enzyme))
    
    #Obtain cancer samples in dummy.table.1 that contain mutations in metabolite's enzymes
    dummy.G1.samples<-unique(as.vector(dummy.table.1[Hugo_Symbol %in% dummy.metabolite.enzymes, ]$PATIENT))
    
    #Account for fact that there may not be more than 1 sample with mutation, cannot do diff expression in that case
    if (length(dummy.G1.samples)>1) {
      
      #Get expression matrix for patients affected with mutations associated with metabolic enzymes
      dummy.G1.expression<-RNASEQ.MATRICES$tumor[,dummy.G1.samples]
      
      #Do differential expression against normal
      dummy.diff.exp<-Function.RNAseq.Differential.Expression(RNASEQ.MATRICES$normal, dummy.G1.expression)
      
      #Get v(j) by dividing that were differentially expressed over all genes tested
      dummy.G.diff.exp.genes<-nrow(dummy.diff.exp[dummy.diff.exp$adj.P.Val<0.05,]) #Value set at p<0.05 for bonferroni corrected p-values
      dummy.v.metabolite<-dummy.G.diff.exp.genes/nrow(dummy.diff.exp)
      
    } else #If 1 or no sample have mutation, then its influence will be equal to zero
      dummy.v.metabolite<-0
    
    #Add v(metabolite) to table
    dummy.table.v<-rbind(dummy.table.v, data.frame(KEGG_ID=metabolite, v.METABOLITE=dummy.v.metabolite))
    
    #To count
    x=x+1
    print (x/length(Table.u.KEGG_IDs))
  }
  #Return v(j) table
  dummy.table.v<-as.data.table(dummy.table.v)
  return(dummy.table.v)
}

#PROCESS ENTRIES
require(data.table)
require(base)

args<-commandArgs(trailingOnly=T)

input.RNA.SEQ.FOLDER<-args[1]
input.cancer.map.file<-args[2]
input.table.1.rds<-args[3]
input.table.3.rds<-args[4]
input.table.u<-args[5]
output.file<-args[6]

u.KEGG.IDS<-unique(as.vector(input.table.u$STATS$KEGG_ID))

#EXECUTE
Run<-Function.v(input.RNA.SEQ.FOLDER, input.cancer.map.file, input.table.1.rds, input.table.3.rds, u.KEGG.IDS)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)

