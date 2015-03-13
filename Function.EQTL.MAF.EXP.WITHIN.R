#####Function.EQTL.MAF.EXP.WITHIN.R####
#031215
#Function to describe expression cause by mutated class within cancer phenotype#

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.PREP.EXP.MATRIX<-function(exp.rds){
  
  #Load matrix obj
  exp.obj<-readRDS(exp.rds)
  
  #Return
  return(exp.obj)
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
  
  #Remove silent mutations
  maf<-maf[Variant_Classification!="Silent",]
  
  #Separate by type of mutation
  non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
  non.func.type<-c("DEL","INS")
  maf$TYPE<-ifelse(maf$Variant_Classification %in% non.func.class, "NON.FUNC", 
                   ifelse(maf$Variant_Type %in% non.func.type, "NON.FUNC", "MISS"))
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Get one sample representing all vials
  maf$VIAL<-substr(maf$SAMPLE, 14,16) #Get vial
  maf$SAMPLE<-substr(maf$SAMPLE,1,12)
  maf[,REP.VIAL:= names(sort(table(VIAL),decreasing=T))[1] , by=SAMPLE] #Get highest counting vial ID
  maf$SAMPLE<-paste( maf$SAMPLE, maf$REP.VIAL, sep="." ) #Rename sample with highest counting vial ID
  maf$VIAL<-NULL
  maf$REP.VIAL<-NULL
  
  #Count how many samples are covered by a type of mutation
  maf[,POP.MUT:=length(SAMPLE),by=c("Hugo_Symbol","Start_Position","TYPE")]
  
  #Cleand up and Return
  setkey(maf)
  maf<-unique(maf)
  maf<-maf[order(POP.MUT, decreasing=T),]
  return(maf)
}

Function.Main<-function(maf, exp.obj){
  
  #Filter maf and exp.obj by common patients
  samples<-intersect(unique(maf$SAMPLE), colnames(exp.obj$combined.matrices))
  maf<-maf[SAMPLE %in% samples,]
  exp.matrix<-exp.obj$combined.matrices[,samples]
  
  #Classify mutation by functional impact, that is MISS are position dependent and NON.FUNC
  maf$CLASS<-ifelse(maf$TYPE=="MISS", paste(maf$Hugo_Symbol, maf$Start_Position, maf$TYPE, sep="."), paste(maf$Hugo_Symbol, maf$TYPE))
  
  #Filter out those mutations classes that do not appear in at least 2 samples
  maf[,CLASS.POP:=length(SAMPLE), by="CLASS"]
  maf<-maf[CLASS.POP>=2,]
  
  #Obtain classes and all samples
  mut.classes<-unique(as.vector(maf$CLASS))
  all.samples<-unique(as.vector(maf$SAMPLE))
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix","data.table", "mut.classes", "maf","all.samples") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-parLapply(cl,mut.classes, function(x) {
    
    #Obtain population for mutation class
    subtable<-maf[CLASS==x,]
    pop<-unique(subtable$SAMPLE)
    non.pop<-setdiff(all.samples, pop)
    
    #Calculate p.values for both sides
    greater.pval<-apply(exp.matrix, 1, function(x) wilcox.test(x[pop], x[non.pop], paired=F, alternative="greater")$p.value)
    less.pval<-apply(exp.matrix, 1, function(x) wilcox.test(x[pop], x[non.pop], paired=F, alternative="less")$p.value)
    
    #Correct for multiple hypothesis with fdr
    greater.pval.adj<-p.adjust(greater.pval, method="fdr")
    less.pval.adj<-p.adjust(less.pval, method="fdr")
    
    #Build table per mutation class
    class.table<-data.table(Hugo_EXP=rownames(exp.matrix), PLUS.PVAL=greater.pval, MINUS.PVAL=less.pval, 
                            PLUS.PVAL.ADJ=greater.pval.adj, MINUS.PVAL.ADJ=less.pval.adj, N.SAMPLES=length(pop),
                            MUT.CLASS=x)
    
    #Return
    return(class.table)
    
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with parallelization")
  
  #Combine tables
  main.table<-do.call(rbind, main.list)
  
  #Return
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
maf.file<-args[1]
exp.rds<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
exp.matrix<-Function.PREP.EXP.MATRIX(exp.rds)
print ("done prepping expression rds")

maf<-Function.Prep.MAF(maf.file)
print ("done prepping maf")

main.function<-Function.Main(maf, exp.matrix)
print ("Done with main function")

###############WRITING OUTPUT############
write.table(file=output.file, main.function, sep="\t", quote=F, row.names=F, col.names=T)
print ("Done writing to file")