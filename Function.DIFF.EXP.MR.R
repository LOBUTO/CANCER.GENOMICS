#FUNCTION.DIFF.EXP.MR.R
#012115
#Check for differential expression between patients containing mutations in desired genes and those that do not

######INTRODUCE ARGUMENTS######
library(data.table)

args<-commandArgs(trailingOnly=T) 

maf.in<-args[1] #TCGA maf file used in mr.in build!!!
mr.in<-readRDS(args[2]) #Target mr as in 012115.MET.ENRICHED.FILTERED.rds
exp.in<-readRDS(args[3]) #Expression object of the form 102514.CANCER.MATRICES.NORMALIZED.OBJ.rds
cancer<-args[4] #Type of cancer to be analyzed (i.e. "BRCA")
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

Function.RNAseq.Differential.Expression.V3<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  #MODIFIED: To count number of differentially expressed gene at adj.p.val<0.05
  
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
  
  #Count differentially expressed genes
  n.diff<-nrow(all.G.fit[adj.P.Val<0.05,])
  
  #Return topTable
  return(n.diff)
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
    samples<-unique(as.vector(maf.table[Hugo_Symbol %in% gene,]$PATIENT))
    not.samples<-unique(as.vector(maf.table[!(PATIENT %in% samples),]$PATIENT))
    
    #Check that we have at least one sample to do differential expression
    if (length(samples)>1){
      
      #Do differential expression to obtain number of diff expression
      diff.exp<-Function.RNAseq.Differential.Expression.V3(exp.in, samples)
      diff.exp.not<-Function.RNAseq.Differential.Expression.V3(exp.in, not.samples)
      
      #Count samples
      n.samples<-length(samples)
      n.not.samples<-length(not.samples)
      
    } else {
      #If only one sample then experiment will be inconclusive
      diff.exp<-NA
      diff.exp.not<-NA
      n.samples<-NA
      n.not.samples<-NA
    }
    
    #Return number of differentially expressed genes
    return(list(N.DIFF=diff.exp, N.NOT.DIFF=diff.exp.not, N.SAMPLES=n.samples, N.NOT.SAMPLES=n.not.samples))
  }
  
  #Execute comparisson of metabolite related gene differential expression vs rest
  main.table<-main.table[,internal.function(GENE),by=c("METABOLITE","KEGG_ID")]
  
  #Return
  return(main.table)
}

Function.Main.PVAL<-function(main.command, maf.table,exp.in){
  #Calculates p-values per metabolite in main.command run
  
  #Function for simulation 
  internal.function<-function(n.samples,n.diff, n.not.samples, n.not.diff){
    
    n.diff.pval<-mean(replicate(100, Function.RNAseq.Differential.Expression.V3(exp.in, sample(all.samples, n.samples))) >=n.diff)
    n.not.diff.pval<-mean(replicate(100, Function.RNAseq.Differential.Expression.V3(exp.in, sample(all.samples, n.not.samples))) >=n.not.diff)
    
    return(list(N.DIFF.PVAL=n.diff.pval, N.NOT.DIFF.PVAL=n.not.diff.pval))
  }
  
  #Exectute
  all.samples<-unique(as.vector(maf.table$PATIENT))
  main.table<-main.command[,internal.function(N.SAMPLES, N.DIFF, N.NOT.SAMPLES, N.NOT.DIFF),by=c("METABOLITE","KEGG_ID")]
  
  #Clean up and return
  main.table<-main.table[order(N.DIFF.PVAL),]
  return(main.table)
}

#######EXECUTION########
maf.table<-Function.Process.MAF(maf.in)
cat ("Processed maf file\n")

main.table<-Function.MR.Process(mr.in, cancer)
cat ("Processed mr file\n")

cat ("Performing expression analysis\n")
main.command<-Function.Main(main.table, maf.table, exp.in)

cat ("Obtaining p values (Simmulations) \n")
main.command.pval<-Function.Main.PVAL(main.command, maf.table, exp.in)

######OUTPUT#########
write.table(file=output.file, main.command.pval, sep="\t", quote=F, row.names=F, col.names=T)
cat ("DONE")

