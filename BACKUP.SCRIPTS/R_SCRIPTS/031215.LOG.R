#031215
########FUNCTIONS#############
library(data.table)
library(reshape2)
library(ggplot2)
library(parallel)
library(gplots)
library(pheatmap)
library(ROCR)
library(h2o)
library(mlbench)
library(caret)

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
  
  #Classify mutations by classes
  maf$CLASS<-ifelse(maf$TYPE=="MISS", paste(maf$Hugo_Symbol, maf$Start_Position, maf$TYPE, sep="."), paste(maf$Hugo_Symbol, maf$TYPE, sep="."))
  
  #Count how many samples are covered by a type of mutation
  maf[,POP.CLASS:=length(SAMPLE),by="CLASS"]
  
  #Cleand up and Return
  setkey(maf)
  maf<-unique(maf)
  maf<-maf[order(POP.CLASS, decreasing=T),]
  return(maf)
}

Function.Main<-function(maf, exp.obj){
  
  #Filter maf and exp.obj by common patients
  samples<-intersect(unique(maf$SAMPLE), colnames(exp.obj$combined.matrices))
  maf<-maf[SAMPLES %in% samples,]
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
  main.list<-lapply(mut.classes, function(x) {
    
    #Obtain population for mutation class
    subtable<-maf[CLASS==x,]
    pop<-unique(subtable$SAMPLE)
    non.pop<-setdiff(all.samples, pop)
    
    #Calculate p.values for both sides
    greater.pval<-parApply(cl, exp.matrix, 1, function(x) wilcox.test(x[pop], x[non.pop], paired=F, alternative="greater")$p.value)
    less.pval<-parApply(cl, exp.matrix, 1, function(x) wilcox.test(x[pop], x[non.pop], paired=F, alternative="less")$p.value)
    
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

Function.PLOT.CLASS.COMP<-function(maf, exp.obj, mut.class, exp.gene){
  
  #Get all available samples and filter tables
  all.samples<-intersect(unique(maf$SAMPLE), colnames(exp.obj$combined.matrices))
  maf<-maf[SAMPLE %in% all.samples,]
  exp.obj$combined.matrices<-exp.obj$combined.matrices[,all.samples]
  
  #Get target samples
  samples<-as.vector(maf[CLASS==mut.class, ]$SAMPLE)
  non.samples<-setdiff(all.samples, samples)
  
  #Construct plot
  target.exp<-exp.obj$combined.matrices[exp.gene, samples]
  non.target.exp<-exp.obj$combined.matrices[exp.gene, non.samples]
  plot.table<-data.table(EXP=c(target.exp, non.target.exp), SOURCE=c(rep(mut.class, length(target.exp)), rep("REST", length(non.target.exp))))
  
  #Plot
  ggplot(plot.table, aes(SOURCE, EXP, colour=SOURCE)) + geom_boxplot() + geom_jitter() + theme.format
  
}

Function.ROC<-function(target.table, target.genes, first="gene", second="q", METHOD="CANCER", curve="PR", p.val=T){
  #Assumes that first column is gene symbol and second column is p-value column
  
  if (p.val==T){
    pred.scores<-ifelse(as.vector(target.table[[second]])==0, 100, -log(as.vector(target.table[[second]])))  
  } else {
    pred.scores<-as.vector(target.table[[second]])
  }
  
  pred.genes<-as.vector(target.table[[first]])
  target.pred<-prediction(pred.scores, pred.genes %in% target.genes)
  
  #Choose type of curve
  if (curve=="PR"){
    target.perf<-performance(target.pred, "prec", "rec")  
  } else if (curve=="ROC"){
    target.perf<-performance(target.pred, "tpr", "fpr")  
  }
  
  #Convert to ready to ggplot format
  pred.table<-data.table(RECALL=target.perf@x.values[[1]], PRECISSION=target.perf@y.values[[1]])
  pred.table$METHOD<-METHOD
  
  #Return
  return(pred.table) 
}

Function.EQTL.SIG<-function(eqtl.table, maf){
  EQTL.SIG<-eqtl.table[MINUS.PVAL.ADJ<0.05 | PLUS.PVAL.ADJ<0.05,]
  
  EQTL.SIG<-EQTL.SIG[, list(PLUS.COUNT=sum(PLUS.PVAL.ADJ<0.05), MINUS.COUNT=sum(MINUS.PVAL.ADJ<0.05)), by="MUT.CLASS"] #get influence counts per side
  EQTL.SIG$COUNT<-EQTL.SIG$PLUS.COUNT + EQTL.SIG$MINUS.COUNT #get total count
  EQTL.SIG$Hugo_Symbol<-sapply(EQTL.SIG$MUT.CLASS, function(x) unlist(strsplit(unlist(strsplit(x," "))[1], "[.]"))[1] )
  EQTL.SIG<-EQTL.SIG[,list(SCORE=mean(COUNT)), by="Hugo_Symbol"]
  EQTL.SIG<-EQTL.SIG[order(SCORE, decreasing=T),] 
  EQTL.SIG<-rbind(EQTL.SIG, 
                  data.table(Hugo_Symbol=setdiff(unique(maf$Hugo_Symbol),EQTL.SIG$Hugo_Symbol), SCORE=0 ))
  
  return (EQTL.SIG)
}

Function.PROCESS.ANNOVAR<-function(MAF.ANNOVAR){
  #Assuming columns 1 and 6 have been dropped from avinput_exonic_variant_function file
  #Assumes files of type 122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST.avinput.exonic_variant_function as input
  
  require(Biostrings)
  require(car)
  
  #Parse into hugo_symbol
  MAF.ANNOVAR$V3<-sapply(MAF.ANNOVAR$V3, function(x) strsplit(x,":")[[1]][1]) 
  
  #Assign labels
  setnames(MAF.ANNOVAR,  c("TYPE", "Hugo_Symbol","Chrom", "Position", "REF","ALT","TRIMER","MT","EXON","AF", "PHAST"))
  
  #Filter out unknown types
  MAF.ANNOVAR<-MAF.ANNOVAR[TYPE!="unknown",]
  
  #Get MAF from AF
  MAF.ANNOVAR$MAF<-ifelse(MAF.ANNOVAR$AF>0.5, 1-MAF.ANNOVAR$AF, MAF.ANNOVAR$AF)
  
  #Clean up
  MAF.ANNOVAR<-MAF.ANNOVAR[!is.na(MT),]
  MAF.ANNOVAR<-MAF.ANNOVAR[EXON==TRUE,]
  
  #Equalize REF-ALT pairs
  MAF.ANNOVAR$PRE.REF.ALT<-paste(as.vector(MAF.ANNOVAR$REF), as.vector(MAF.ANNOVAR$ALT), sep="_")
  MAF.ANNOVAR$REF.ALT<-recode(MAF.ANNOVAR$PRE.REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  
  #Equalize trimers
  MAF.ANNOVAR$RECODE<-MAF.ANNOVAR$PRE.REF.ALT!=MAF.ANNOVAR$REF.ALT
  TEMP.RECODE<-MAF.ANNOVAR[RECODE==T,]
  TEMP.RECODE$TRIMER<-substring(complement(DNAStringSet(TEMP.RECODE$TRIMER)),1,3)
  MAF.ANNOVAR<-rbind(MAF.ANNOVAR[RECODE==F,],TEMP.RECODE)
  
  #Clean up and Return
  MAF.ANNOVAR$PRE.REF.ALT<-NULL
  MAF.ANNOVAR$RECODE<-NULL
  return(MAF.ANNOVAR)
}

Function.PREP.MAF.CLASS.ML<-function(MAF.ANNOVAR, REP.TIME, NT.LENGTH, EXTRA=c(), FILTER=F){
  #Added replication time and rep class
  
  require(data.table)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Filter out low counts?
  if (FILTER==T){
    MIN.MAF<-min(MAF.ANNOVAR$MAF[MAF.ANNOVAR$MAF!=0])
    MAF.ANNOVAR<-MAF.ANNOVAR[MAF>MIN.MAF,]
  }
  
  #Enter rep info
  MAF.ANNOVAR<-merge(MAF.ANNOVAR, REP.TIME, by="Hugo_Symbol")
  
  #Enter length info
  MAF.ANNOVAR<-merge(MAF.ANNOVAR, NT.LENGTH, by="Hugo_Symbol")
  MAF.ANNOVAR$NT<-normalize.vector(MAF.ANNOVAR$NT)
  
  #Classify MAF into categories
  MAF.ANNOVAR$MAF.CUT<-cut(MAF.ANNOVAR$MAF, quantile(MAF.ANNOVAR$MAF, c(0,0.33,0.66,1)), include.lowest=T, 
                           labels=as.factor(c("low", "medium", "high")))
  
  #Get cummulative MAF per hugo
  MAF.ANNOVAR[,HUGO.MAF.SUM:=sum(MAF), by="Hugo_Symbol"]
  
  #Filter out not needed columns and place classificaiton column at the end
  MAF.ANNOVAR<-MAF.ANNOVAR[,c("TYPE", "Chrom", "PHAST", "REF.ALT","REP.CLASS", "TRIMER", "MT", "MAF.CUT", "REP.TIME", "HUGO.MAF.SUM", "NT",
                              EXTRA), with=F]  
  
  #Return
  return(MAF.ANNOVAR)
}
##############################

########TESTS################
BRCA.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf")

GBM.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/GBM/031315/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf")

KIRC.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/KIRC/031615/a7d9d1ef-89f3-49b0-8555-251df7709f88/Somatic_Mutations/BCM__Mixed_DNASeq/Level_2/hgsc.bcm.edu__Mixed_DNA_Sequencing_level2.maf")

OV.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/OV/031415/29b94157-be45-41e6-bcfd-de5f49fbe5b5/Somatic_Mutations/BCM__SOLiD_DNASeq_curated/Level_2/hgsc.bcm.edu__SOLiD_curated_DNA_sequencing_level2.maf")

BRCA.MAF[CLASS=="NOTCH3 NON.FUNC",]
BRCA.MAF[Hugo_Symbol=="TTN" & POP.CLASS>1,]
Function.PLOT.CLASS.COMP(BRCA.MAF, brca.exp.nb, "TTN NON.FUNC", "ELN")

#####Analyze results of within cancer class influence on phenotype###
#BRCA
BRCA.EQTL<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/031215.BRCA.EQTL.WITHIN", header=T, sep="\t", stringsAsFactors=F)

BRCA.MAF[CLASS %in% unique(c(unique(BRCA.EQTL[PLUS.PVAL.ADJ<0.05,]$MUT.CLASS), unique(BRCA.EQTL[MINUS.PVAL.ADJ<0.05,]$MUT.CLASS))),]
length(unique(BRCA.MAF[CLASS %in% unique(c(unique(BRCA.EQTL[PLUS.PVAL.ADJ<0.05,]$MUT.CLASS), unique(BRCA.EQTL[MINUS.PVAL.ADJ<0.05,]$MUT.CLASS))),]$SAMPLE))

BRCA.EQTL.SIG<-BRCA.EQTL[MINUS.PVAL.ADJ<0.05 | PLUS.PVAL.ADJ<0.05,] #get significant scores

BRCA.EQTL.SIG<-BRCA.EQTL.SIG[, list(PLUS.COUNT=sum(PLUS.PVAL.ADJ<0.05), MINUS.COUNT=sum(MINUS.PVAL.ADJ<0.05)), by="MUT.CLASS"] #get influence counts per side
BRCA.EQTL.SIG$COUNT<-BRCA.EQTL.SIG$PLUS.COUNT + BRCA.EQTL.SIG$MINUS.COUNT #get total count
BRCA.EQTL.SIG$Hugo_Symbol<-sapply(BRCA.EQTL.SIG$MUT.CLASS, function(x) unlist(strsplit(unlist(strsplit(x," "))[1], "[.]"))[1] )
BRCA.EQTL.SIG<-BRCA.EQTL.SIG[,list(SCORE=mean(COUNT)), by="Hugo_Symbol"]
BRCA.EQTL.SIG<-BRCA.EQTL.SIG[order(SCORE, decreasing=T),] 
BRCA.EQTL.SIG<-rbind(BRCA.EQTL.SIG, 
                     data.table(Hugo_Symbol=setdiff(unique(BRCA.MAF$Hugo_Symbol),BRCA.EQTL.SIG$Hugo_Symbol), SCORE=0 ))

BRCA.EQTL.SIG.PLOT<-Function.ROC(BRCA.EQTL.SIG, as.vector(COSMIC.BRCA$Hugo_Symbol), "Hugo_Symbol", "SCORE", "WITHIN.EQTL","ROC",p.val=F)
MUTSIG.BRCA.PLOT<-Function.ROC(BRCA.CURATED.MUTSIG[,c("gene","q"),with=F], as.vector(COSMIC.BRCA$Hugo_Symbol), "gene", "q", "MUTSIG", "ROC" )

ggplot(rbind(BRCA.EQTL.SIG.PLOT, MUTSIG.BRCA.PLOT), aes(RECALL, PRECISSION, colour=METHOD)) + geom_line() +theme.format +
  xlab("FPR") + ylab("TPR")

#GBM
GBM.EQTL<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/GBM/031415.GBM.EQTL.WITHIN.MA", header=T, sep="\t", stringsAsFactors=F)
GBM.EQTL.SIG<-Function.EQTL.SIG(GBM.EQTL, GBM.MAF)

GBM.CURATE.MUTSIG<-fread("DATABASES/CANCER_DATA/MUTSIG/GBM/GBM.SIG.txt", header=T, sep="\t", stringsAsFactors=F)

GBM.EQTL.SIG.PLOT<-Function.ROC(GBM.EQTL.SIG, COSMIC.GBM$Hugo_Symbol, "Hugo_Symbol","SCORE", "WITHIN.EQTL", "ROC", p.val=F)
MUTSIG.GBM.PLOT<-Function.ROC(GBM.CURATE.MUTSIG[,c("gene","q"),with=F], as.vector(COSMIC.GBM$Hugo_Symbol), "gene", "q", "MUTSIG", "ROC" )
ggplot(rbind(GBM.EQTL.SIG.PLOT, MUTSIG.GBM.PLOT), aes(RECALL, PRECISSION, colour=METHOD)) + geom_line() +theme.format +
  xlab("FPR") + ylab("TPR")

#######Create Function to normalize CNV data, into binary, high(2), low(1) or 0 as maf files##### - POSTPONED
length(unique(GBM.MAF$SAMPLE))
length(unique(GBM.MAF[CLASS %in% c("IDH1.209113112.MISS", "TP53.NON.FUNC","ATRX.NON.FUNC","NF1.NON.FUNC","RB1.NON.FUNC" ),]$SAMPLE))
GBM.MAF[Hugo_Symbol=="PIK3R1",]

########FOCUS ON ML APPROACH TO PREDICT BACKGROUND MUTATION RATE USING 1000G DATA####
MAF.ANNOVAR<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/ANNOVAR/122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST.avinput.exonic_variant_function", drop=c(1,6),
                   sep="\t", header=F, stringsAsFactors=F)

MAF.ANNOVAR<-Function.PROCESS.ANNOVAR(MAF.ANNOVAR)
table(MAF.ANNOVAR$TYPE)

ggplot(MAF.ANNOVAR[,list(COUNT=length(TYPE)), by="Hugo_Symbol"], aes(COUNT)) + geom_histogram() w +
  scale_x_log10()
ggplot(MAF.ANNOVAR, aes(MAF, PHAST, colour=REF.ALT)) + geom_point() + theme.format + facet_grid(REF.ALT~TYPE)

#Introduce replication time
CHEN.REP<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/010815.CHEN.REP.TIMES.MID.2", header=T, sep="\t",stringsAsFactors=F, drop=2)

#Introduce length
EXON.NT<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/011915.EXON.COORDIANTES.FIXED.OVERLAPS", header=T, sep="\t", stringsAsFactors=F, drop=2:6)
EXON.NT<-EXON.NT[,list(NT=sum(FEAT_END-FEAT_START)), by=Hugo_Symbol]

######## Predicting with H2O on continous MAF variable#########################

# Start a local cluster with 8GB RAM
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g') 

#Convert data into h2o
h2o_maf <- as.h2o(localH2O, MAF.ANNOVAR, key = 'MAF.ANNOVAR') 

#Break into train and test set
set.seed(1234)
ALL_ROWS<-1:nrow(h2o_maf)
RAND_FOLDS<-createFolds(ALL_ROWS,5)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:3])
TEST_ROWS<-unlist(RAND_FOLDS[4:5])
TRAIN_MAF<-as.factor(MAF.ANNOVAR$MAF[TRAIN_ROWS]) #Train with 60%
TEST_MAF<-as.factor(MAF.ANNOVAR$MAF[TEST_ROWS]) #Test on 40%

#Model without dropout
MAF.MODEL<-h2o.deeplearning(x=c(1,3,7,8,11,13), y=12, data=h2o_maf[TRAIN_ROWS,], classification=F,
                            activation = "Tanh", balance_classes = TRUE, hidden = c(100,100), epochs = 100)

#Evaluate 
maf_train <- h2o.predict(MAF.MODEL, h2o_maf[TRAIN_ROWS,])$predict
maf_train <- as.factor(as.matrix(maf_train))
maf_test <- h2o.predict(MAF.MODEL, h2o_maf[TEST_ROWS, ])$predict
maf_test <- as.factor(as.matrix(maf_test))

confusionMatrix(maf_train, TRAIN_MAF)$overall[1]

TRAIN_MAF$PREDICTED<-maf_train[,1] #Differences are obvious

#Shut down H2o
h2o.shutdown(localH2O)

#############################################################################################

######## Predicting with H2O on classifiable MAF variable#########################
MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR,  CHEN.REP, EXON.NT, "MAF",FILTER=T) #NOTE: Using filtered version!!!!

# Start a local cluster with 10GB RAM
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '10g', nthreads=-1) 

#Convert data into h2o
h2o_maf <- as.h2o(localH2O, MAF.ANNOVAR.CLASS) 

#Break into train and test set
set.seed(1234)
ALL_ROWS<-1:nrow(h2o_maf)
RAND_FOLDS<-createFolds(ALL_ROWS,5)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:3])
TEST_ROWS<-unlist(RAND_FOLDS[4:5])
TRAIN_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TRAIN_ROWS]) #Train with 60%
TEST_MAF<-as.factor(MAF.ANNOVAR.CLASS$MAF.CUT[TEST_ROWS]) #Test on 40%

#Model without dropout
MAF.MODEL<-h2o.deeplearning(x=c(1:7,11), y=8, data=h2o_maf[TRAIN_ROWS[1:100000],],
                            activation = "Tanh", balance_classes = TRUE, hidden = c(500,100,10), epochs = 300)
                            
MAF.MODEL.FILTERED.NT.1000.2.50.50.E.300<-copy(MAF.MODEL)
MAF.MODEL.FILTERED.1000.2.50.50.E.300<-copy(MAF.MODEL)
MAF.MODEL.FILTERED.NT.10000.2.100.100.E.300<-copy(MAF.MODEL)
MAF.MODEL.FILTERED.NT.10000.3.500.250.50.E.300<-copy(MAF.MODEL)
MAF.MODEL.FILTERED.NT.REP.TIME.10000.2.100.100.E.300<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O//MAF.MODEL.FILTERED.NT.REP.TIME.10000.2.100.100.E.300.h2o")
MAF.MODEL.LAST<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/MAF.MODEL.FILTERED.NT.REP.CLASS.100000.3.5000.1000.100.E.500.h2o/")

#Evaluate 
maf_train <- h2o.predict(MAF.MODEL.FILTERED.NT.REP.TIME.10000.2.100.100.E.300, h2o_maf[TRAIN_ROWS[1:10000],])$predict
maf_train <- as.factor(as.matrix(maf_train))
maf_test <- h2o.predict(MAF.MODEL.FILTERED.NT.REP.TIME.10000.2.100.100.E.300, h2o_maf[TEST_ROWS[1:10000], ])$predict
maf_test <- as.factor(as.matrix(maf_test))

x<-prop.table(table(maf_train, TRAIN_MAF[1:10000]))
mean(diag(x[colnames(x), colnames(x)])/rowSums(x[colnames(x), colnames(x)]))

prop.table(table(maf_train,TRAIN_MAF[1:10000]),1)
table(maf_test,TEST_MAF[1:10000])

Function.EVALUATE.DEEP<-function(learners.vector, h2o_maf, train_rows, test_rows, train_maf, test_maf){
  
  main.list<-lapply(learners.vector, function(x) {
    
    #Obtain model
    x.model<-get(x)
    
    #Obtain trainers length
    if (grepl("REP.TIME", x)){
      l<-unlist(strsplit(x, "[.]"))[7]
    } else if(grepl("NT", x)){
      l<-unlist(strsplit(x, "[.]"))[5]
    } else{
      l<-unlist(strsplit(x, "[.]"))[4]
    }
    train_rows<-train_rows[1:as.numeric(l)]
    train_maf<-train_maf[1:as.numeric(l)]
    
    #Evaluate 
    maf_train <- h2o.predict(x.model, h2o_maf[train_rows,])$predict
    maf_train <- as.factor(as.matrix(maf_train))
    maf_test <- h2o.predict(x.model, h2o_maf[test_rows, ])$predict
    maf_test <- as.factor(as.matrix(maf_test))
    
    #Table
    x.train<-prop.table(table(maf_train, train_maf))
    x.test<-prop.table(table(maf_test, test_maf))
    
    #Calculate
    train.score<-mean(diag(x.train[colnames(x.train), colnames(x.train)])/rowSums(x.train[colnames(x.train), colnames(x.train)]))
    test.score<-mean(diag(x.test[colnames(x.test), colnames(x.test)])/rowSums(x.test[colnames(x.test), colnames(x.test)]))
    
    #Return
    return(data.table(METHOD=x, TRAIN.SCORE=train.score, TEST.SCORE=test.score))
  })
  
  #Combined
  main.table<-do.call(rbind, main.list)
}

MAIN.EVAL<-Function.EVALUATE.DEEP(ls(pattern="MAF.MODEL.FILTERED"), h2o_maf, TRAIN_ROWS, TEST_ROWS[1:10000], TRAIN_MAF, TEST_MAF[1:10000])
MAIN.EVAL<-melt(MAIN.EVAL, id.vars="METHOD")
ggplot(MAIN.EVAL, aes(METHOD, value, colour=variable)) + geom_histogram(position="dodge",stat="identity") +theme.format + 
  theme(text = element_text(size=6),axis.text.x = element_text(angle = 65, hjust = 1))

#Shut down H2o
h2o.shutdown(localH2O)

###############TESTING RESULTS FROM CLUSTER RUN#########################
MAF.MODEL.10000.3.1000.500.100.E.400<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/MAF.MODEL.10000.3.1000.500.100.E.400.REP.CLASS.h2o/")
MAF.MODEL.10000.3.100.50.10.E.500<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/MAF.MODEL.10000.3.100.50.10.E.500.REP.CLASS.h2o/")
MAF.MODEL.10000.3.100.50.10.E.1000<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/MAF.MODEL.10000.3.100.50.10.E.1000.REP.CLASS.h2o/")
MAF.MODEL.10000.3.100.50.10.E.2000<-h2o.loadModel(localH2O, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/MAF.MODEL.10000.3.100.50.10.E.2000.REP.CLASS.h2o/")

MAF.TEST.PREDICT<-h2o.predict(MAF.MODEL.10000.3.100.50.10.E.2000, h2o_maf[TRAIN_ROWS[1:25000] ,])$predict
table(as.factor(as.matrix(MAF.TEST.PREDICT)), TRAIN_MAF[1:25000])

########################################################################

######## LINEAR PREDICTIONS OF MAF CLASS#########################
MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR, CHEN.REP, EXON.NT, "MAF")

ALL_ROWS<-1:nrow(MAF.ANNOVAR.CLASS)
RAND_FOLDS<-createFolds(ALL_ROWS,3)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:2]) #Train with 66%
TEST_ROWS<-unlist(RAND_FOLDS[3]) #Test on 33%
MAF.ANNOVAR.TEST<-MAF.ANNOVAR.CLASS[TEST_ROWS,]
MAF.ANNOVAR.TEST$SET<-"TEST"
MAF.ANNOVAR.TRAIN<-MAF.ANNOVAR.CLASS[TRAIN_ROWS,]
MAF.ANNOVAR.TRAIN$SET<-"TRAIN"

#Using binomial model first
GLM.MODEL.BINOMIAL.REP.CLASS<-glm(MAF.CUT~TYPE+PHAST+NT+REF.ALT+REP.CLASS+Chrom+TRIMER, 
                                  data=MAF.ANNOVAR.CLASS[TRAIN_ROWS,], family=binomial)
MAF.ANNOVAR.TRAIN$PREDICT<-predict.glm(GLM.MODEL.BINOMIAL.REP.CLASS, MAF.ANNOVAR.CLASS[TRAIN_ROWS, ], type="response")
MAF.ANNOVAR.TEST$PREDICT<-predict.glm(GLM.MODEL.BINOMIAL.REP.CLASS, MAF.ANNOVAR.CLASS[TEST_ROWS, ], type="response")

ggplot(rbind(MAF.ANNOVAR.TEST, MAF.ANNOVAR.TRAIN), aes(MAF.CUT, PREDICT, colour=MAF.CUT)) + geom_boxplot() + geom_jitter(size=0.3) + theme.format +
  facet_wrap(~SET)

GLM.MODEL.BINOMIAL.REP.TIME<-glm(MAF.CUT~, data=MAF.ANNOVAR.CLASS[TRAIN_ROWS, ], family=binomial) 
MAF.ANNOVAR.TRAIN$PREDICT<-predict.glm(GLM.MODEL.BINOMIAL.REP.TIME, MAF.ANNOVAR.CLASS[TRAIN_ROWS, ], type="response")
MAF.ANNOVAR.TEST$PREDICT<-predict.glm(GLM.MODEL.BINOMIAL.REP.TIME, MAF.ANNOVAR.CLASS[TEST_ROWS, ], type="response")

ggplot(rbind(MAF.ANNOVAR.TEST, MAF.ANNOVAR.TRAIN), aes(MAF.CUT, PREDICT, colour=MAF.CUT)) + geom_boxplot() + geom_jitter(size=0.3) + theme.format +
  facet_wrap(~SET)

#Using poison
MAF.ANNOVAR.CLASS.POISON<-copy(MAF.ANNOVAR.CLASS)
MAF.ANNOVAR.CLASS.POISON$MAF.CUT<-as.numeric(MAF.ANNOVAR.CLASS.POISON$MAF.CUT) #Converting data to poisson friendly

GLM.MODEL.POISSON.REP.CLASS<-glm(MAF.CUT~., data=MAF.ANNOVAR.CLASS.POISON[TRAIN_ROWS, c(1:8), with=F], family=poisson)
MAF.ANNOVAR.TRAIN$PREDICT<-predict.glm(GLM.MODEL.POISSON.REP.CLASS, MAF.ANNOVAR.CLASS.POISON[TRAIN_ROWS, ], type="response")
MAF.ANNOVAR.TEST$PREDICT<-predict.glm(GLM.MODEL.POISSON.REP.CLASS, MAF.ANNOVAR.CLASS.POISON[TEST_ROWS, ], type="response")

ggplot(rbind(MAF.ANNOVAR.TEST, MAF.ANNOVAR.TRAIN), aes(MAF.CUT, PREDICT, colour=MAF.CUT)) + geom_boxplot() + geom_jitter(size=0.3) + theme.format +
  facet_wrap(~SET)

GLM.MODEL.POISSON.REP.TIME<-glm(MAF.CUT~., data=MAF.ANNOVAR.CLASS.POISON[TRAIN_ROWS, c(1:4,6:9), with=F], family=poisson)
MAF.ANNOVAR.TRAIN$PREDICT<-predict.glm(GLM.MODEL.POISSON.REP.TIME, MAF.ANNOVAR.CLASS.POISON[TRAIN_ROWS, ], type="response")
MAF.ANNOVAR.TEST$PREDICT<-predict.glm(GLM.MODEL.POISSON.REP.TIME, MAF.ANNOVAR.CLASS.POISON[TEST_ROWS, ], type="response")

ggplot(rbind(MAF.ANNOVAR.TEST, MAF.ANNOVAR.TRAIN), aes(MAF.CUT, PREDICT, colour=MAF.CUT)) + geom_boxplot() + geom_jitter(size=0.3) + theme.format +
  facet_wrap(~SET)