######Function.METABOLIC.COR.R#######
#040615
#Calculate deviation of correlation score for all enzymatic genes in cancer expression matrix

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.Prep.MAF<-function(maf.file) {
  
  require(car)
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode",
              "Reference_Allele", "Tumor_Seq_Allele2"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #Remove silent mutations
  #maf<-maf[Variant_Classification!="Silent",]
  
  #Equalize REF-ALT pairs
  maf$PRE.REF.ALT<-ifelse(maf$Variant_Type=="SNP", paste(as.vector(maf$Reference_Allele), as.vector(maf$Tumor_Seq_Allele2), sep="_"), "-")
  maf$REF.ALT<-ifelse(maf$PRE.REF.ALT!="-",recode(maf$PRE.REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" '),"-")
  maf$PRE.REF.ALT<-NULL
  
  #Separate by type of mutation
  non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
  non.func.type<-c("DEL","INS")
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "SILENT",
                   ifelse(maf$Variant_Classification %in% non.func.class, "NON.FUNC", 
                          ifelse(maf$Variant_Type %in% non.func.type, "NON.FUNC", "MISS")))
  
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
  
  #Classify mutations by classes
  maf$CLASS<-ifelse((maf$TYPE=="MISS" | maf$TYPE=="SILENT"), 
                    paste(maf$Hugo_Symbol, maf$Start_Position, maf$REF.ALT, maf$TYPE, sep="."), 
                    paste(maf$Hugo_Symbol, maf$TYPE, sep="."))
  
  #Count how many samples are covered by a type of mutation
  maf[,POP.CLASS:=length(SAMPLE),by="CLASS"]
  
  #Filter out silent
  maf<-maf[TYPE!="SILENT",]
  
  #Classify mut class into two categories non.functional and gain of function - CHANGE LATER FOR SPECIFIC GAIN OF FUNCTION!!
  maf$CLASS<-paste(maf$Hugo_Symbol, maf$TYPE, sep=".")
  
  #Cleand up and Return
  maf<-maf[,c("Hugo_Symbol", "SAMPLE", "CLASS"), with=F]
  setkey(maf)
  maf<-unique(maf)
  #maf<-maf[order(POP.CLASS, decreasing=T),]
  return(maf)
}

Function.exp<-function(exp.obj){
  
  #Load expresssion file
  BRCA.EXP<-readRDS(exp.obj)
  
  #Obtain normal and cancer correlations
  #BRCA.EXP.COR.NORMAL<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$normal.patients]), method="spearman")
  #BRCA.EXP.COR.CANCER<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients]), method="spearman")
  
  #Return
  return(list(NORMAL=BRCA.EXP$combined.matrices[,BRCA.EXP$normal.patients],
              CANCER=BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients]))
}

Function.Enzyme<-function(enzyme.file){
  
  #Read enzyme file
  ENZYME<-fread(enzyme.file, header=T, sep="\t", stringsAsFactors=F)
  
  #Return
  return(ENZYME)
  
}

Function.Main<-function(maf, cancer.exp, normal.exp, enzymes){
  
  ########Construct expression correlation matrix for each mut class using cancer expression data#######
  
  #Crossfilter maf and exp matrix
  all.patients<-intersect(unique(maf$SAMPLE), colnames(cancer.exp))
  maf<-maf[SAMPLE %in% all.patients,]
  cancer.exp<-cancer.exp[, all.patients]
  
  #Extract mutation classes
  mut.class<-unique(as.vector(maf$CLASS))
  
  #Filter exp matrix for target enzymes
  cancer.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  normal.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  
  #Obtain background enzymatic correlation
  BRCA.EXP.COR.NORMAL<-cor(t(normal.exp), method="spearman")
  genes.order<-rownames(BRCA.EXP.COR.NORMAL)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "mut.class", "cancer.exp", "maf", "BRCA.EXP.COR.NORMAL",
                              "gene.order") ,envir=environment())
  print ("Done exporting values")
  
  main.list<-parLapply(cl, mut.class, function(x) {
    
    #Obtain samples under class
    samples<-as.vector(maf[CLASS %in% x,]$SAMPLE)
    
    #Perform class correlation in cancer
    cancer.cor<-cor(t(cancer.exp[,samples]), method="spearman")
    
    #############Substract background correlation(NORMAL) from main correlation (CANCER)###############
    #Substract
    signal.cor<-cancer.cor[genes.order, genes.order]-BRCA.EXP.COR.NORMAL
    
    #Get sum of absolute value of change of enzyme with respect to all (LOG-CONVERTED)
    sum.delta.abs<-log(apply(signal.cor, 1, function(x) sum(abs(x))))
    
    #Organize to table
    main.table<-matrix(sum.delta.abs)
    rownames(main.table)<-rownames(signal.cor)
    
    #Clean up and return
    return(list(MAIN.TABLE=main.table, N.SAMPLES=length(samples)))
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Return
  return(main.list)
  
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
cancer.maf<-args[1]
exp.rds<-args[2]
enzyme.file<-args[3]
output.file<-args[4]
print("opened files")
##########################################

##################EXECUTE#################
maf<-Function.Prep.MAF(cancer.maf)
print ("done prepping maf")

exp<-Function.exp.cor(exp.rds)
print ("done building cancer and normal correlations")

enzyme<-Function.Enzyme(enzyme.file)
print ("done reading enzyme file")

main.function<-Function.Main(maf, exp$CANCER, exp$NORMAL, unique(as.vector(enzyme$Enzyme)))
print ("done with main function")

###############WRITING OUTPUT############
saveRDS(object=main.function, file=output.file)
print ("Done writing to file")