#######Function.MI.CANCER.ANALYSIS.R########
#030715


##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.MAF<-function(maf.file, mut.filter=T) {
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #Filter for Nonfunctional mutations types if needed
  if (mut.filter==T){
    non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
    non.func.type<-c("DEL","INS")
    maf<-maf[!(Variant_Classification %in% non.func.class),]
    maf<-maf[!(Variant_Type %in% non.func.type),]  
  }
  
  #For now classify as change and no-change, later we will also calculate separate probabilities for INS, DEL and non-missense SNPs
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "NO.CHANGE","CHANGE")
  
  #Filter for non-silent only
  maf<-maf[TYPE=="CHANGE",]
  maf$TYPE<-NULL
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Return
  setkey(maf)
  maf<-unique(maf)
  return(maf)
}

Function.Prep.EXP<-function(exp.rds, paired=T){
  
  #Load expression matrix file
  exp.matrix<-readRDS(exp.rds)
  
  #Work with paired samples only
  if (paired==T){
    normal.patients<-substr(exp.matrix$normal.patients, 1, 12)
    
    #Get cancer patients with a matched normal
    cancer.patients<-data.table(samples=exp.matrix$cancer.patients, patients=substr(exp.matrix$cancer.patients, 1, 12))
    cancer.patients$MATCH<-cancer.patients$patients %in% normal.patients
    cancer.patients<-cancer.patients[MATCH==TRUE, ]$samples
    
    #Filter expression matrix for it
    exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, cancer.patients)]
    exp.matrix$cancer.patients<-cancer.patients
    
    #Return filtered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)
    
  } else {
    #Return unfiltered matrix
    print (dim(exp.matrix$combined.matrices))
    return(exp.matrix)  
  }
  
}

Function.PREP.MAF<-function(maf, train.patients, POS=T){
  #Preparing the maf table for solving with expression data as a linear system of equations
  
  require(data.table)
  
  #Filter maf table for patients of interest
  maf<-maf[SAMPLE %in% train.patients,]
  
  #Merge hugo and position to give a unique ID to mutation -or- just MUT
  if (POS==T){
    maf$CHANGE<-paste(maf$Hugo_Symbol, maf$Start_Position, sep=".")  
  } else {
    maf<-maf[,c("SAMPLE", "Hugo_Symbol"), with=F]
    setnames(maf, c("SAMPLE", "CHANGE"))
  }

  #Simplify table
  maf<-maf[,c("SAMPLE", "CHANGE"), with=F]
  
  #Add count to cast table
  maf$COUNT<-1
  maf[,POS.SAMPLES:=sum(COUNT), by="CHANGE"]
  
  #Cast table to long format
  maf.cast<-acast(maf, SAMPLE~CHANGE ,fill=0,value.var="COUNT", fun.aggregate=sum)
  
  #Clean up and Return as list
  maf<-maf[,c("CHANGE","POS.SAMPLES"),with=F]
  setkey(maf)
  maf<-unique(maf)
  return(list(CAST=maf.cast, TABLE=maf)) 
}

Function.EXP.DIFF<-function(exp.obj, genes, target.patients){
  require(data.table)
  
  #Extract data
  if (genes=="ALL"){
    exp.matrix<-exp.obj$combined.matrices
  } else {
    exp.matrix<-exp.obj$combined.matrices[genes,,drop=F]  
  }
  normal<-exp.obj$normal.patients
  cancer<-target.patients
  
  #Get difference to mean normal for genes of interense
  main.table<-data.table(apply(exp.matrix, 1, function(x) x[cancer]-mean(x[normal])), keep.rownames=T)
  
  #Label data
  setnames(main.table, c("SAMPLE", colnames(main.table)[2: ncol(main.table)]))
  
  #Return
  return(main.table)
}

Function.SOLVE.EXP<-function(prepped.maf.cast, exp.diff.table, prepped.maf.table, train.patients){
  #Solve for coefficients
  
  require(limSolve)
  
  #Get expression of genes that we are trying to solve for
  exp.genes<-colnames(exp.diff.table)[2: ncol(exp.diff.table)]
  
  #Sample order that we are solving for
  target.samples<-exp.diff.table$SAMPLE
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "exp.genes", "prepped.maf.cast", "target.samples","exp.diff.table", 
                              "Solve","setnames") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-parLapply(cl, exp.genes, function(x) {
    
    #Solve equation
    solved.table<-data.table(as.matrix(Solve(prepped.maf.cast[target.samples,], exp.diff.table[[x]]) ), keep.rownames=T)
    
    #Clean table
    setnames(solved.table, c("CHANGE", x))
    
    #Return
    return(solved.table)
  } )
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Re-merge tables
  main.table<-Reduce(function(x,y) merge(x,y, by="CHANGE"), main.list)
  
  #Merge with position table to keep count
  main.table<-merge(main.table, prepped.maf.table[,c("CHANGE","POS.SAMPLES"),with=F], by="CHANGE")
  
  #Return
  return(list(TABLE=main.table, TRAIN.PATIENTS=train.patients))
  
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cancer.maf<-args[1]
exp.rds<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
maf<-Function.Prep.MAF(cancer.maf, mut.filter=F)
print ("done prepping maf")

exp.matrix<-Function.Prep.EXP(exp.rds, paired=F)
print ("done prepping expression rds")

train.patients<-sample(intersect(unique(maf$SAMPLE), exp.matrix$cancer.patients),750)

maf.obj<-Function.PREP.MAF(maf, train.patients, POS=F)
print ("done casting maf")

exp.diff<-Function.EXP.DIFF(exp.matrix, "ALL", train.patients)
print ("done prepping expression diff table")

main.function<-Function.SOLVE.EXP(maf.obj$CAST, exp.diff, maf.obj$TABLE, train.patients)
print ("Done with main function")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing to file")