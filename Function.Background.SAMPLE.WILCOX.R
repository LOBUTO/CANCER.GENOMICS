#Function.Background.SAMPLE.WILCOX.R
#022315
#Function to build background distributions for wilcoxon diff expression across number of samples

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.MAF<-function(maf.file) {
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #Filter for Nonfunctional mutations types
#   non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
#   non.func.type<-c("DEL","INS")
#   maf<-maf[!(Variant_Classification %in% non.func.class),]
#   maf<-maf[!(Variant_Type %in% non.func.type),]
  
  #For now classify as change and no-change, later we will also calculate separate probabilities for INS, DEL and non-missense SNPs
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "NO.CHANGE","CHANGE")
  
  #Filter for non-silent only
  maf<-maf[TYPE=="CHANGE",]
  maf$TYPE<-NULL
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Return
  return(maf)
}

Function.Main<-function(linkage.nc.file, exp.rds, maf){
  
  #Load linkage file
  main.table<-fread(linkage.nc.file, header=T, sep="\t", stringsAsFactors=F)
  print ("Done loading linkage file")
  
  #Obtain vector of population sizes
  population<-unique(main.table$N.PATIENTS)
  print (population)
  
  #Load expression file
  exp.matrix<-readRDS(exp.rds)
  print ("Done loading expression file")
  
  #Get set of total patients with expression and mutation data
  patients<-intersect(unique(maf$SAMPLE), colnames(exp.matrix$combined.matrices))
  
  #Filter both tables for these patients
  maf<-maf[SAMPLE %in% patients,]
  exp.matrix$cancer.patients<-intersect(exp.matrix$cancer.patients, patients)
  exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, exp.matrix$cancer.patients)]
  
  #Subset samples
  normal<-unique(exp.matrix$normal.patients)
  cancer<-unique(exp.matrix$cancer.patients)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("population", "normal","cancer", "as.data.table","exp.matrix","data.table") ,envir=environment())
  print ("Done exporting values")
  
  #Executing
  print ("Finding background")
  main.list<-list()
  
  #Loop for each sample size in population of samples
  for (pop in population){
    print (pop)
    
    #Create random populations of pop size
    random.sampling<-replicate(200, sample(cancer,pop), simplify=F)
    
    #Calculate proportion of diff expressed genes for each random population
    random.prop<-parSapply(cl, random.sampling, function(x){
      
      #Proportion is based on fdr corrected diff exp genes
      wilcox.pvals<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[x], y[normal], paired=F)$p.value)
      wilcox.vector<-p.adjust(as.vector(wilcox.pvals), method="fdr" )
      wilcox.ratio<-mean(wilcox.vector<0.05)
      
      #Return diff exp ratio
      return(wilcox.ratio)
      
    })
    
    #Append random diff exp ratio population to list
    main.list[[as.character(pop)]]<-random.prop
  }
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with background distributions")
  
  #Return
  return(main.list)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
linkage.nc.file<-args[1]
exp.rds<-args[2]
cancer.maf<-args[3]
output.file<-args[4]
print("opened files")
##########################################

##################EXECUTE#################
maf<-Function.Prep.MAF(cancer.maf)
print ("done prepping maf")

main.function<-Function.Main(linkage.nc.file, exp.rds, maf)
print ("Done")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing file")