#Function.EXP.MUT.POS.LINKAGE.NC.R
#021515
#Calculates the linkage between mutation position and expression between normal and cancer populations
#Will initially only compare paired samples
#Will only take into account sites that are seen mutated at least 2 times in population

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)
library(qvalue)

Function.BRCA.SUBTYPE<-function(brca.normalized.obj, version=1){
  #Takes normalized expression object from breast cancer and returns subtype score and assigns subtype to each patient based on maximum score
  
  require(genefu)
  require(data.table)
  
  #Choose pam50 method
  if (version==1){
    data(pam50)
    pam<-copy(pam50)
  } else if(version==2){
    data(pam50.scale)
    pam<-copy(pam50.scale)
  } else if(version==3){
    data(pam50.robust)
    pam<-copy(pam50.robust)
  }
  
  #Obtain target pam gene in expression matrix
  target.genes<-intersect(rownames(pam$centroids), rownames(brca.normalized.obj$combined.matrices))
  
  #Simplify pam and expression matrices
  simplified.pam<-pam$centroids[target.genes,]
  simplified.exp<-brca.normalized.obj$combined.matrices[target.genes, brca.normalized.obj$cancer.patients]
  
  #Obtain predicted type based on correlation
  BRCA.SCORES<-cor(simplified.exp, simplified.pam, method="spearman")
  BRCA.SCORES<-as.data.frame(BRCA.SCORES)
  BRCA.SCORES$PATIENT<-rownames(BRCA.SCORES)
  BRCA.SCORES$TYPE<-colnames(BRCA.SCORES)[max.col(BRCA.SCORES[,1:5])]
  
  #Clean up and return
  BRCA.SCORES<-as.data.table(BRCA.SCORES)
  BRCA.SCORES<-BRCA.SCORES[order(TYPE),]
  return(BRCA.SCORES)
}

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
  non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
  non.func.type<-c("DEL","INS")
  maf<-maf[!(Variant_Classification %in% non.func.class),]
  maf<-maf[!(Variant_Type %in% non.func.type),]
  
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

internal.function<-function(n, exp.matrix, patients, normal.patients){
  
  #Obtain count of genes that are significantly differentiated between sample and normal cohorts
  target.patients<-sample(patients, n)
  p.vals<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[target.patients], y[normal.patients], paired=F)$p.value)
  p.vals.adj<-p.adjust(p.vals, method="fdr")
  
  return(sum(p.vals.adj<0.05))
}

Function.Main<-function(maf, exp.matrix){
  
  #Get set of total patients with expression and mutation data
  patients<-intersect(unique(maf$SAMPLE), colnames(exp.matrix$combined.matrices))
  
  #Filter both tables for these patients
  maf<-maf[SAMPLE %in% patients,]
  exp.matrix$cancer.patients<-intersect(exp.matrix$cancer.patients, patients)
  exp.matrix$combined.matrices<-exp.matrix$combined.matrices[,c(exp.matrix$normal.patients, exp.matrix$cancer.patients)]
  
  #Filter for samples that are seeing mutated in at least more than 2 samples at the same position
  maf[,N.SAMPLES.POS:=length(unique(SAMPLE)), by=c("Start_Position","Hugo_Symbol")]
  maf[,N.MUT.HUGO.SAMPLE:=length(Start_Position), by=c("Hugo_Symbol","SAMPLE")] #Number of mutations per gene in each sample
  maf<-maf[N.SAMPLES.POS>=2,]
  
  #Check that we have enough data to pursue
  if(nrow(maf)>1){
    
    #Simplify maf table
    maf<-maf[,c("Hugo_Symbol","SAMPLE", "Start_Position", "N.SAMPLES.POS"),with=F]
    setkey(maf)
    maf<-unique(maf)
    print (maf[order(N.SAMPLES.POS, decreasing=T),])
    
    #Split into tables 
    main.list<-split(maf, list(maf$Hugo_Symbol, maf$Start_Position),drop=T)
    
    #Prepping parallelization
    print ("prepping for parallelization")
    
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    clusterExport(cl, varlist=c("main.list", "maf","as.data.table","exp.matrix","data.table") ,envir=environment())
    print ("Done exporting values")
    
    #Execute parallelization
    print ("Finding linkage")
    main.table<-parLapply(cl, main.list, function(x) {
      
      mut.gene<-unique(x$Hugo_Symbol)
      target.patients<-unique(x$SAMPLE)
      non.target.patients<-exp.matrix$normal.patients
      
      wilcox.matrix<-apply(exp.matrix$combined.matrices, 1, function(y) wilcox.test(y[target.patients], y[non.target.patients], paired=F)$p.value)
      wilcox.matrix<-data.table(Hugo_MUT=mut.gene, Position.MUT=unique(x$Start_Position) , Hugo_EXP=names(wilcox.matrix), 
                                N.PATIENTS=length(target.patients), P.VAL=as.vector(wilcox.matrix))  
      
      wilcox.matrix$P.VAL.ADJ<-p.adjust(wilcox.matrix$P.VAL, method="fdr" )
      
      return(wilcox.matrix)
    })
    
    #Stop parallelization
    stopCluster(cl)
    print ("Done with linkages")
    
    #Merge gene tables
    main.table<-do.call(rbind, main.table)
    
    #Cleaup and return
    main.table<-main.table[order(P.VAL.ADJ),]
    return(main.table)
    
    #If no maf data tto work with, return empty table
  } else
    return(data.table())
}

Function.Main.Subtype<-function(maf, exp.matrix){
  
  #Classify based on brca subtype
  subtype.table<-Function.BRCA.SUBTYPE(exp.matrix)
  subtypes<-unique(as.vector(subtype.table$TYPE))
  
  #Apply main function by subtype
  subtype.list<-list()
  print ("Calculating linkage per subtype")
  for (type in subtypes) {
    print (type)
    print (subtype.table[TYPE==type,])
    
    #Check that we have enough data in each subtype
    if (nrow(subtype.table[TYPE==type,])>1){
      
      #Obtain subtype samples
      type.samples<-unique(as.vector(subtype.table[TYPE==type,]$PATIENT))
      type.maf<-maf[SAMPLE %in% type.samples,]
      
      #Check that we have enough data in processed maf file
      print (type.maf)
      if (nrow(type.maf)>1){
        #Execute
        subtype.linkage<-Function.Main(type.maf, exp.matrix)
        
        #Check that we have linkage data per subtype
        if(nrow(subtype.linkage)>0){
          #Classify and append to main list
          subtype.linkage$SUBTYPE<-type
          subtype.list[[type]]<-subtype.linkage      
        }
      }
    }
  }
  
  #Combine all linkage tables
  subtype.main<-do.call(rbind, subtype.list)
  
  #Return
  return(subtype.main)
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
maf<-Function.Prep.MAF(cancer.maf)
print ("done prepping maf")

exp.matrix<-Function.Prep.EXP(exp.rds, paired=F)
print ("done prepping expression rds")

main.function<-Function.Main.Subtype(maf, exp.matrix)
print ("Done with main function")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)
print ("Done writing to file")