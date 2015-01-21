#FUNCTION.DIFF.EXP.MR.R
#012115
#Check for differential expression between patients containing mutations in desired genes and those that do not

######INTRODUCE ARGUMENTS######
library(data.table)

args<-commandArgs(trailingOnly=T) 

maf.in<-args[1] #TCGA maf file 
mr.in<-readRDS(args[2]) #Target mr
exp.in<-readRDS(args[3]) #Expression object of the form 102514.CANCER.MATRICES.NORMALIZED.OBJ.rds
cancer<-args[4] #Type of cancer to be analyzed
output.file<-args[5]

#######FUNCTIONS#######
Function.Process.MAF<-function(maf.file){
  #Function to clean maf files straight from TCGA
  
  #Process for missense only
  maf.record<-fread(maf.file, sep="\t", header=T, stringsAsFactors=F, drop=c(2:4,7:8, 10,12,14:15,17:37))
  maf.record<-maf.record[Variant_Classification=="Missense_Mutation",]
  maf.record<-maf.record[nchar(Tumor_Seq_Allele2)==1 & nchar(Reference_Allele)==1,]
  
  #Convert to patient ID for later expression processing
  maf.record$PATIENT<-sapply(maf.record$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x,"-"))[1:4], collapse="."))
  
  #Clean up and return
  maf.record<-maf.record[,c("Hugo_Symbol","PATIENT"),with=F]
  setkey(maf.record)
  maf.record<-unique(maf.record)
  
  #Return
  return(maf.record)
}

Function.RNAseq.Differential.Expression.V2<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  
  require(limma)
  require(edgeR)
  require(data.table)
  
  #Get target combined matrix
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  target.combined.matrix<-normalized.matrices.object$combined.matrices[,c(cancer.samples.in.matrix,normal.samples)]
  
  #Get design matrix
  G.design.matrix<-data.frame(G=c(rep("G1", length(cancer.samples.in.matrix)), rep("G0", length(normal.samples))))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Perform differential expression
  G.fit = lmFit(target.combined.matrix, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Convert to data.table
  if ("ID" %in% colnames(all.G.fit)) {
    all.G.fit<-as.data.table(all.G.fit)
  } else {
    all.G.fit$ID<-rownames(all.G.fit)
    all.G.fit<-as.data.table(all.G.fit)
  }
  
  #Return topTable
  return(all.G.fit)
}

Function.MR.Process<-function(mr.in, cancer){
  #Filter mr table for target cancer
 
  main.table<-mr.in[CANCER==cancer,]
  
  #Return
  return(main.table)
}

Function.Main<-function(main.table, maf.table, exp.in){
  #Calculates number of differentially expressed genes per metabolite
  
  internal.function<-function(gene){
    samples<-as.vector(maf.table[Hugo_Symbol %in% gene,]$PATIENT)
    not.samples<-as.vector(maf.table[!(Hugo_Symbol %in% gene),]$PATIENT)
    
    #Check that we have at least one sample to do differential expression
    if length(samples)>1{
      
      #Do differential expression
      diff.exp<-Function.RNAseq.Differential.Expression.V2(exp.in, samples)
      diff.exp.not<-Function.RNAseq.Differential.Expression.V2(exp.in, not.samples)
      
      #Count differentially expressed genes
      n.diff<-nrow(diff.exp[adj.P.Val<0.05,])
      n.not.diff<-nrow(diff.exp.not[adj.P.Val<0.05,])
      
    } else {
      #If only one sample then experiment will be inconclusive
      n.diff<-NA
      n.not.diff<-NA
    }
    
    #Return number of differentially expressed genes
    return(list(N.DIFF=n.diff, N.NOT.DIFF=n.not.diff))
  }
  
  #Execute comparisson of metabolite related gene differential expression vs rest
  main.table<-main.table[,internal.function(GENE),by=c("METABOLITE","KEGG_ID")]
  
  #Return
  return(main.table)
}

#######EXECUTE########
