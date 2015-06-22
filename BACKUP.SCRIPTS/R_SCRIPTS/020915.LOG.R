#020915
library(data.table)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gplots)

#Load cancer data
brca.maf<-fread("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t",stringsAsFactors=F)
brca.maf<-brca.maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode"),with=F]
setkey(brca.maf)
brca.maf<-unique(brca.maf)

#Analyze changes
table(brca.maf$Variant_Classification)
unique(brca.maf[,4:5, with=F])[order(Variant_Type),]

#Filter for Nonfunctional mutations types
non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
non.func.type<-c("DEL","INS")
brca.maf<-brca.maf[!(Variant_Classification %in% non.func.class),]
brca.maf<-brca.maf[!(Variant_Type %in% non.func.type),]

#For now classify as change and no-change, later we will also calculate separate probabilities for INS, DEL and non-missense SNPs
brca.maf$TYPE<-ifelse(brca.maf$Variant_Classification=="Silent", "NO.CHANGE","CHANGE")

#Filter for non-silent only
brca.maf<-brca.maf[TYPE=="CHANGE",]
brca.maf$TYPE<-NULL

brca.maf$SAMPLE<-sapply(brca.maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
brca.maf$Tumor_Sample_Barcode<-NULL

#Introduce length info
exon<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/011915.EXON.COORDIANTES.FIXED.OVERLAPS", header=T, sep="\t",stringsAsFactors=F)
exon<-exon[,list(exon.nt=sum(FEAT_END-FEAT_START)),by=Hugo_Symbol]

#Introduce metabolic info
table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
setnames(table.2, c("METABOLITE","KEGG_ID","Hugo_Symbol"))

Function.Empirical.PVALS<-function(maf, table.2, exon){
  
  #Internal functions
  internal.function<-function(sample.n.mut, met.rate, KEGG) {
    
    met.table<-copy(exon)
    
    #Get each gene probabilities for sampling
    gene.prob<-met.table$exon.nt/sum(exon$exon.nt)
    
    #Classify table into gene that are and aren't assocaited with metabolite
    met.table$MET<-met.table$Hugo_Symbol %in% unique(table.2[KEGG_ID==KEGG,]$Hugo_Symbol)
    
    #Simulate sampling 1000 times
    met.dist<-replicate(1000, {
      main<-sample(met.table$Hugo_Symbol, sample.n.mut, replace=T, prob=gene.prob)
      main<-met.table[Hugo_Symbol %in% main, ]
      met.test<-sum(main$MET)
      return(met.test)
    })
    
    #Calcualte p-val from sampling distribution for values that are as extreme or greater than our metabolic rate
    p.val<-mean(met.dist>=met.rate)
    
    #return p-value
    return(list(P.VAL=p.val))
  }
  
  #Obtain number of mutations per patient
  maf[,sample.n.mut:=length(Start_Position), by="Tumor_Sample_Barcode"]
  
  #Merge with metabolic info
  main.table<-merge(maf, table.2, by="Hugo_Symbol")
  
  #Add a mutation count per line
  #main.table$MUT.COUNT<-1
  
  #Calculate metabolic mutation count per metabolite in each patient
  main.table[,met.mut.count:=length(Start_Position), by=c("Tumor_Sample_Barcode", "METABOLITE","KEGG_ID")]
  
  #Split table to parallelize job
  main.split<-split(main.table, list(main.table$Tumor_Sample_Barcode, main.table$METABOLITE, main.table$KEGG_ID), drop=T)
  main.names<-names(main.split)
  print ("Done spliting tables")
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  library(parallel)
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("main.split", "as.data.table", "internal.function", "data.table", "table.2","copy", "exon") ,envir=environment())
  
  #Calculate empirical P.VAL per metabolite per patient (use sample.n.mut, met.rate)
  print (main.table)
  print ("Simulations")
  
  main.table<-parLapply(cl, main.split, function(x){
   met.table<-x[,internal.function(unique(sample.n.mut), unique(met.mut.count), unique(KEGG_ID)),
                         by=c("Tumor_Sample_Barcode", "METABOLITE","KEGG_ID")]
   return(met.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print("done calculating p-vals")
  
  #Restoing table format
  main.table<-do.call(rbind, main.table)
  
  #Correct for multiple hypothesis testing
  print ("correcting for multiple hypothesis")
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Clean up and return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}

#Test function
function.test<-Function.Empirical.PVALS(brca.maf[Tumor_Sample_Barcode %in% unique(brca.maf$Tumor_Sample_Barcode)[1:10],],
                                        table.2, exon)

###############021115#############

#Analysis
brca.met<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/021015.BRCA.PERSONALIZED.MET", head=T, sep="\t", stringsAsFactors=F)
brca.met$P.VAL.ADJ<-p.adjust(brca.met$P.VAL, method="fdr")
brca.met[P.VAL.ADJ<0.1,]
length(unique(brca.met[P.VAL.ADJ<0.1,]$Tumor_Sample_Barcode))

length(unique(brca.met[P.VAL.ADJ<0.05,]$KEGG_ID))

#Draw a heatmap at 
brca.met.heat<-brca.met[P.VAL.ADJ<0.05,]
brca.met.heat$UNIT<-1
brca.met.heat<-acast(brca.met.heat, Tumor_Sample_Barcode~METABOLITE, value.var="UNIT",fill=0)
dim(brca.met.heat)
pheatmap(brca.met.heat, scale="none",drop_levels=T)

#Introduce pathway info
KEGG.PATH<-fread("DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND", header=T, sep="\t", stringsAsFactors=F)
KEGG.PATH$PATH<-sapply(KEGG.PATH$DESCRIPTION, function(x) unlist(strsplit(x," - "))[1])
KEGG.PATH$DESCRIPTION<-NULL
KEGG.PATH$PATHWAY<-NULL
setnames(KEGG.PATH, c("KEGG_ID","PATH"))

KEGG.PATH<-merge(brca.met.heat[,c(1,3),with=F], KEGG.PATH, by="KEGG_ID")
KEGG.PATH$UNIT<-1
KEGG.PATH$KEGG_ID<-NULL
setkey(KEGG.PATH)
KEGG.PATH<-unique(KEGG.PATH)

KEGG.PATH[,SAMPLE.COVERAGE:=length(unique(Tumor_Sample_Barcode)),by="PATH"]
KEGG.PATH<-KEGG.PATH[!(PATH %in% c("Metabolic pathways", "Biosynthesis of secondary metabolites")),]
length(unique(KEGG.PATH$PATH))

KEGG.HEAT<-KEGG.PATH[SAMPLE.COVERAGE>=7,]
length(unique(KEGG.HEAT$Tumor_Sample_Barcode))


KEGG.HEAT<-acast(KEGG.HEAT, PATH~Tumor_Sample_Barcode, value.var="UNIT", fill=0)

pheatmap(KEGG.HEAT, scale="none",drop_levels=T)

heatmap.2(KEGG.HEAT, scale="none", trace="none")

############## 021314 ###############

#Load cancer census
COSMIC<-fread("DATABASES/CANCER_DATA/COSMIC/cancer_gene_census.csv", head=T, sep=",", stringsAsFactors=F, drop=c(2:5,7,9:12,14:17))
COSMIC<-COSMIC[Somatic=="yes",]
COSMIC$Somatic<-NULL

setnames(COSMIC, c("Hugo_Symbol","TUMOR","TYPE"))
COSMIC<-COSMIC[grepl("Mis",TYPE),]
COSMIC$TYPE<-NULL

COSMIC.BRCA<-rbind(COSMIC[grepl("breast",TUMOR,ignore.case=T),], COSMIC[grepl("brca",TUMOR,ignore.case=T),] )

#Plot site frequencies of cosmic genes
for (gene in c(COSMIC.BRCA$Hugo_Symbol, "BRCA1","GATA3")){
  file.name<-paste("FIGURES/COSMIC/BREAST",gene, "jpeg",sep=".")
  q<-ggplot(brca.maf[Hugo_Symbol==gene,][,list(N.SAMPLE=length(SAMPLE)),by=Start_Position], aes(Start_Position, N.SAMPLE)) + geom_point() + theme.format
  ggsave(filename=file.name, plot=q,scale=3)  
}

##Use BRCA.maf build linear model of mutation to phenotype trait
Function.MUT.MET<-function(mut.gene, exp.gene, position=F, paired=F){
  require(data.table)
  
  #Filter tables
  all.patients<-intersect(unique(brca.maf$SAMPLE), colnames(brca.exp.nb$combined.matrices))
  brca.maf<-brca.maf[SAMPLE %in% all.patients,]
  
  #Analysis based on paired only?
  if (paired==T){
    pair.ids<-sapply(brca.exp.nb$normal.patients, function(x) substr(x, 1, 12))
    paired.cancer<-brca.exp.nb$cancer.patients[substr(brca.exp.nb$cancer.patients,1,12) %in% pair.ids] 
    brca.maf<-brca.maf[SAMPLE %in% paired.cancer,]
    brca.exp.nb$cancer.patients<-paired.cancer
  }
  
  #Target patients
  if (position==F){
    
    if (mut.gene=="ALL"){
      target.patients<-brca.exp.nb$cancer.patients
    } else{
      target.patients<-unique(brca.maf[Hugo_Symbol %in% mut.gene,]$SAMPLE)  
    }
    non.target.patients<-brca.exp.nb$normal.patients 
    
  } else {
    
    if (mut.gene=="ALL"){
      target.patients<-brca.exp.nb$cancer.patients 
    } else{
      target.patients<-unique(brca.maf[Hugo_Symbol %in% mut.gene & Start_Position==position,]$SAMPLE)  
    }
    non.target.patients<-brca.exp.nb$normal.patients
  }
  
  #Apply
  gene.exp<-as.vector(brca.exp.nb$combined.matrices[exp.gene, target.patients])
  gene.not<-as.vector(brca.exp.nb$combined.matrices[exp.gene, non.target.patients])
  
  #Prepare for plot
  main.table<-data.table(EXP=c(gene.exp, gene.not), TYPE=c( rep("MUT",length(gene.exp)), rep("NON.MUT", length(gene.not)) ))
  
  #Plot
  ggplot(main.table, aes(TYPE, EXP, colour=TYPE)) + geom_boxplot() + geom_jitter() + theme.format +
    ggtitle(paste(mut.gene, exp.gene, sep=" - "))
  
}

Function.MUT.MET("TTN","A1BG", 179514309)
Function.MUT.MET("PIK3CA","A1BG",178916938)
brca.maf[Hugo_Symbol=="TTN" & Start_Position==179514309,]
brca.maf[Hugo_Symbol=="PIK3CA" & Start_Position==178916938,]

brca.maf[,list(N.SAMPLE.POS=length(unique(SAMPLE))),by=c("Start_Position","Hugo_Symbol")][order(N.SAMPLE.POS, decreasing=T),]
y<-brca.maf[,list(N.SAMPLE.POS=length(unique(SAMPLE))),by=c("Start_Position","Hugo_Symbol")][N.SAMPLE.POS>=2,]
y<-split(y, y$Hugo_Symbol)
y[[1]]

#######Test results using ROC curves

COSMIC<-fread("DATABASES/CANCER_DATA/COSMIC/021815.CENSUS.csv", head=T, sep=",", stringsAsFactors=F, drop=c(2:5,7,9:12,14:17))
COSMIC<-COSMIC[Somatic=="yes",]
COSMIC$Somatic<-NULL

setnames(COSMIC, c("Hugo_Symbol","TUMOR","TYPE"))
COSMIC<-COSMIC[grepl("Mis",TYPE),]
COSMIC$TYPE<-NULL

COSMIC.BRCA<-COSMIC[grepl("breast",TUMOR,ignore.case=T),]
COSMIC.GBM<-rbind(COSMIC[grepl(c("glioblastoma"), TUMOR, ignore.case=T),], COSMIC[grepl(c("gbm"), TUMOR, ignore.case=T),])

library(ROCR)
Function.ROC<-function(target.table, target.genes, first="gene", second="q", METHOD="CANCER", curve="PR", p.val=T){
  #Assumes that first column is gene symbol and second column is p-value column
  
  if (p.val==T){
    pred.scores<-ifelse(as.vector(target.table[[second]])==0, 100, -log(as.vector(target.table[[second]])))  
  } else {
    pred.scores<-as.vector(target.table[[second]])
  }
  
  pred.genes<-as.vector(target.table[[first]])
  target.pred<- (pred.scores, pred.genes %in% target.genes)
  
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

#Test mutsig on all cosmic
BRCA.CURATED.MUTSIG<-fread("DATABASES/CANCER_DATA/MUTSIG/sig_genes/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.1.1.0.curated.somatic-OUT.sig_genes.txt",
                           header=T, sep="\t", stringsAsFactors=F)

#Test on position only method
BRCA.POS.LINKAGE<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/021615.BRCA.MUT.EXP.LINKAGE.POS.NC.ALL", header=T, sep="\t", stringsAsFactors=F)
BRCA.POS.LINKAGE<-unique(BRCA.POS.LINKAGE[,c("Hugo_MUT", "P.VAL.ADJ"), with=F])
BRCA.POS.LINKAGE<-BRCA.POS.LINKAGE[,list(q=min(P.VAL.ADJ)), by="Hugo_MUT"]

#Test on position missense method
BRCA.POS.LINKAGE.MIS<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/022515.BRCA.EXP.POS.MUT.LINKAGE.NC.ALL.MISSENSE", header=T, sep="\t", stringsAsFactors=F)

#Test on position all method
BRCA.POS.LINKAGE.ALL<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/022515.BRCA.EXP.POS.MUT.LINKAGE.NC.ALL", header=T, sep="\t", stringsAsFactors=F)

#Test bayes method
brca.bayes<-test.bayes[,list(PROB=max(BAYES.PROB)),by="Hugo_Symbol"]
brca.bayes.degree<-test.bayes[,list(PROB=max(BAYES.PROB)),by="Hugo_Symbol"]

#Compare methods
MUTSIG.BRCA.ROC<-Function.ROC(BRCA.CURATED.MUTSIG[,c("gene","q"),with=F], as.vector(COSMIC.BRCA$Hugo_Symbol), "gene", "q", "MUTSIG", "PR" )
POS.BRCA.ROC<-Function.ROC(BRCA.POS.LINKAGE, as.vector(COSMIC.BRCA$Hugo_Symbol), "Hugo_MUT","q", "POS.LINKAGE", "ROC")
BAYES.BRCA.ROC<-Function.ROC(brca.bayes, as.vector(COSMIC.BRCA$Hugo_Symbol), "Hugo_Symbol", "PROB", "BAYES", "ROC",p.val=F)
BAYES.BRCA.DEGREE.ROC<-Function.ROC(brca.bayes.degree, as.vector(COSMIC.BRCA$Hugo_Symbol), "Hugo_Symbol", "PROB", "BAYES.DEGREE", "ROC",p.val=F)

#Compare methods
ggplot(rbind(POS.BRCA.ROC, MUTSIG.BRCA.ROC, BAYES.BRCA.ROC, BAYES.BRCA.DEGREE.ROC), aes(RECALL, PRECISSION, colour=METHOD)) + geom_line() + theme.format


########022315########
#Analysis for Missense
BRCA.POS.LINKAGE.BG.MIS<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/022515.BRCA.LINKAGE.BG.MISSENSE.rds")
BRCA.POS.LINKAGE.BG.MIS<-do.call(rbind, lapply(names(BRCA.POS.LINKAGE.BG.MIS), function(x) data.table(N.PATIENTS=as.numeric(x), DISTRIBUTION=BRCA.POS.LINKAGE.BG.MIS[[x]]) ))
setnames(BRCA.POS.LINKAGE.BG.MIS, c("N.PATIENTS.BG","DISTRIBUTION"))

BRCA.LINKAGE.RATIO<-BRCA.POS.LINKAGE.MIS[,list(EXP.RATIO=mean(P.VAL.ADJ<0.05), N.PATIENTS=unique(N.PATIENTS)), by=c("Hugo_MUT", "Position.MUT")]

ggplot(BRCA.POS.LINKAGE.BG.MIS, aes(as.factor(N.PATIENTS.BG), DISTRIBUTION)) + geom_boxplot() + 
  geom_jitter(data=BRCA.LINKAGE.RATIO, aes(as.factor(N.PATIENTS), EXP.RATIO, colour="EXP.RATIO")) + theme.format + 
  ylab("Differentially Expressed Ratio - Random Sampling") + xlab("Number of Sample on Test")

BRCA.LINKAGE.RATIO[,P.VAL:=mean(as.vector(BRCA.POS.LINKAGE.BG.MIS[N.PATIENTS.BG==N.PATIENTS,]$DISTRIBUTION)>=EXP.RATIO)  , 
                   by=c("Hugo_MUT", "Position.MUT")]
BRCA.LINKAGE.RATIO$P.VAL.ADJ<-p.adjust(BRCA.LINKAGE.RATIO$P.VAL, method="fdr")
BRCA.LINKAGE.RATIO<-BRCA.LINKAGE.RATIO[order(P.VAL.ADJ),]

BRCA.LINKAGE.RATIO[P.VAL.ADJ<0.05,]

#Analysis for All aberrations
BRCA.POS.LINKAGE.BG.ALL<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/022515.BRCA.LINKAGE.BG.ALL.rds")
BRCA.POS.LINKAGE.BG.ALL<-do.call(rbind, lapply(names(BRCA.POS.LINKAGE.BG.ALL), function(x) data.table(N.PATIENTS=as.numeric(x), DISTRIBUTION=BRCA.POS.LINKAGE.BG.ALL[[x]]) ))
setnames(BRCA.POS.LINKAGE.BG.ALL, c("N.PATIENTS.BG","DISTRIBUTION"))

BRCA.LINKAGE.RATIO<-BRCA.POS.LINKAGE.ALL[,list(EXP.RATIO=mean(P.VAL.ADJ<0.05), N.PATIENTS=unique(N.PATIENTS)), by=c("Hugo_MUT", "Position.MUT")]

ggplot(BRCA.POS.LINKAGE.BG.ALL, aes(as.factor(N.PATIENTS.BG), DISTRIBUTION)) + geom_boxplot() + 
  geom_jitter(data=BRCA.LINKAGE.RATIO, aes(as.factor(N.PATIENTS), EXP.RATIO, colour="EXP.RATIO")) + theme.format + 
  ylab("Differentially Expressed Ratio - Random Sampling") + xlab("Number of Sample on Test")

BRCA.LINKAGE.RATIO[,P.VAL:=mean(as.vector(BRCA.POS.LINKAGE.BG.ALL[N.PATIENTS.BG==N.PATIENTS,]$DISTRIBUTION)>=EXP.RATIO)  , 
                   by=c("Hugo_MUT", "Position.MUT")]
BRCA.LINKAGE.RATIO$P.VAL.ADJ<-p.adjust(BRCA.LINKAGE.RATIO$P.VAL, method="fdr")
BRCA.LINKAGE.RATIO<-BRCA.LINKAGE.RATIO[order(P.VAL.ADJ),]

BRCA.LINKAGE.RATIO[P.VAL.ADJ<0.05,]
BRCA.LINKAGE.RATIO[Hugo_MUT=="CDH1",]

#########022415##########
BRCA.POS.LINKAGE.MIS[Hugo_MUT=="CDH1",][order(P.VAL.ADJ),]
Function.MUT.MET("ALL","COL10A1", paired=T)

brca.maf[Hugo_Symbol=="CDH1" & Start_Position==68844139,]$SAMPLE

x<-brca.exp$combined.matrices["A1BG", brca.maf[Hugo_Symbol=="CDH1" & Start_Position==68844139,]$SAMPLE]
y<-brca.exp$combined.matrices["A1BG", normal.patients]
wilcox.test(x, y, paired=F, alternative="less")

brca.subtype<-Function.BRCA.SUBTYPE(brca.exp,version=1)
unique(brca.subtype$TYPE)

##########022515##########
#Disentangling expression using machine learning
brca.exp.nb<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/022415.BRCA.CANCER.MATRICES.NORMALIZED.OBJ.NB.rds")
brca.exp.nb.ann<-data.frame(TYPE=ifelse(colnames(brca.exp.nb$combined.matrices) %in% brca.exp.nb$cancer.patients, "cancer", "normal"), 
                            row.names=colnames(brca.exp.nb$combined.matrices))
pheatmap(brca.exp.nb$combined.matrices[c("A1BG","CDH1"),], scale="none", annotation=brca.exp.nb.ann,cluster_cols=F)

Function.Exp.Pairs<-function(gene.exp, mut.gene=F, position){
  
  #Get pair IDs
  sub.normal<-sapply(brca.exp.nb$normal.patients, function(x) substr(x, 1,12))
  sub.cancer<-sapply(brca.exp.nb$cancer.patients, function(x) substr(x, 1,12))
  
  #Normal table
  normal.table<-data.table(NORMAL=brca.exp.nb$normal.patients, PAIR.ID=sub.normal)
  cancer.table<-data.table(CANCER=brca.exp.nb$cancer.patients, PAIR.ID=sub.cancer)
  
  #Combined
  main.table<-merge(normal.table, cancer.table, by="PAIR.ID")
  
  #Introduce expression level
  main.table$NORMAL.EXP<-as.vector(brca.exp.nb$combined.matrices[gene.exp, main.table$NORMAL])
  main.table$CANCER.EXP<-as.vector(brca.exp.nb$combined.matrices[gene.exp, main.table$CANCER])
  
  if (mut.gene!=F){
    brca.maf<-brca.maf[Hugo_Symbol==mut.gene & Start_Position==position,]
    main.table$CLASSIFICATION<-ifelse(main.table$CANCER %in% as.vector(brca.maf$SAMPLE), paste(mut.gene, " - ", "Position ", position, sep=""), "REST")
    N.PATIENTS<-unique(as.vector(main.table[CLASSIFICATION!="REST",]$CANCER))
    
    #Save table for analysis
    ratio.table<-copy(main.table)
    
    #Melt table
    main.table<-melt(main.table[,c(1,4,5,6),with=F], id.vars=c("PAIR.ID","CLASSIFICATION"))
    setnames(main.table, c("PAIR.ID", "CLASSIFICATION", "TYPE","EXP"))
    
    #Plot object
    my.plot<-ggplot(main.table, aes(TYPE, EXP, group=PAIR.ID, colour=CLASSIFICATION)) + geom_line() + theme.format+
      geom_point( size=4, shape=21, fill="white") + facet_wrap(~CLASSIFICATION) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste(unique(main.table$CLASSIFICATION[main.table$CLASSIFICATION!="REST"]), " - ", "N.PATIENTS: ", length(N.PATIENTS) ,sep=""))
    
    #Classify data
    ratio.table<-ratio.table[CLASSIFICATION!="REST",]
    if(nrow(ratio.table)==0){
      info.vector<-"NO.DATA"
    } else{
      exp.diff<-median(ratio.table$CANCER.EXP-ratio.table$NORMAL.EXP)
      if(abs(exp.diff)<1){
        info.vector<-"NO.CHANGE"
      } else if( exp.diff<1){
        info.vector<-"DOWN"
      } else{
        info.vector<-"UP"
      }
    }
    
    #Return as list
    return(list(plot=my.plot, info=info.vector))
    
  } else{
    
    #Melt table
    main.table<-melt(main.table[,c(1,4,5),with=F], id.vars="PAIR.ID")
    setnames(main.table, c("PAIR.ID", "TYPE","EXP"))  
    
    #Plot
    ggplot(main.table, aes(TYPE, EXP, group=PAIR.ID, colour=PAIR.ID)) + geom_line() + theme.format+
      geom_point( size=4, shape=21, fill="white")
  } 
}

#Try DMD, LOXL4
Function.Exp.Pairs("COL10A1", "PIK3CA", 178952085)
Function.Exp.Pairs("COL10A1")

Function.Exp.Pairs("DMD", "PIK3CA", 178952085)
Function.Exp.Pairs("DMD")

Function.Exp.Pairs("A1BG", "PIK3CA", 178952085)
Function.Exp.Pairs("A1BG")

Function.Exp.Pairs("ABCA10","FAM21A", 51853633)
Function.Exp.Pairs("ABCA10")

Function.Exp.Pairs("SPRR4","PIK3CA", 178936082) 
Function.Exp.Pairs("SPRR4")

Function.Exp.Pairs("ZNF169","PIK3CA", 178952085)
Function.Exp.Pairs("ZNF169")

Function.Exp.Pairs("ZNF169","TP53", 7577539)
Function.Exp.Pairs("ZNF169")

Function.Exp.Pairs("ADAMTS5","PIK3CA", 178936082)
Function.Exp.Pairs("ADAMTS5")

Function.Exp.Pairs("A1BG","CDH1", 68772218)
Function.Exp.Pairs("A1BG")

#BUILD MODEL - Focus on DMD
BRCA.POS.LINKAGE.ALL[Hugo_EXP=="DMD",]
apply(BRCA.POS.LINKAGE.ALL[Hugo_EXP=="ZNF169",], 1, function(x){
  
  filename=paste("022715.", x[1],".",as.numeric(x[2]),".jpeg" , sep="")
  plot.obj<-Function.Exp.Pairs("ZNF169", x[1], as.numeric(x[2]) )
  
  #Store depending on change
  folder=paste("FIGURES/LINKAGE/ZNF169.TEST/", plot.obj$info, sep="")
  jpeg(filename=paste(folder, filename, sep="/") ,width=1200,height=800, quality=100,type="quartz")
  print (plot.obj$plot)
  dev.off()
} )

Function.EXP.UNEARTH<-function(brca.exp.nb){
  
  #Get pair IDs
  sub.normal<-sapply(brca.exp.nb$normal.patients, function(x) substr(x, 1,12))
  sub.cancer<-sapply(brca.exp.nb$cancer.patients, function(x) substr(x, 1,12))
  
  #Normal table
  normal.table<-data.table(NORMAL=brca.exp.nb$normal.patients, PAIR.ID=sub.normal)
  cancer.table<-data.table(CANCER=brca.exp.nb$cancer.patients, PAIR.ID=sub.cancer)
  
  #Combine
  main.table<-merge(normal.table, cancer.table, by="PAIR.ID")
  
  #And then trusting order of matrix combine
  paired.count<-nrow(main.table)
  paired.matrix<-brca.exp.nb$combined.matrices[,c(main.table$NORMAL, main.table$CANCER)]
  
  #And calculate
  diff.vector<-sapply(rownames(paired.matrix), function(x) {
    diff.exp<-as.vector(paired.matrix[x, main.table$CANCER]-paired.matrix[x, main.table$NORMAL])
    
    if ( mean(diff.exp>0.5)==1){
      return(median(diff.exp))
    } else if( mean(diff.exp< -0.5)==1){
      return(median(diff.exp))
    } else {
      return(0)
    }
    
  })
  
  #Annotate diff exp to genes
  diff.table<-data.table(Hugo_Symbol=rownames(paired.matrix), DIFF.EXP=diff.vector)
  
  #Clean up and Return
  diff.table<-diff.table[order(abs(DIFF.EXP), decreasing=T),]
  return(diff.table)
}

diff.table<-Function.EXP.UNEARTH(brca.exp.nb)
diff.table[DIFF.EXP!=0,]
Function.Exp.Pairs("MMP11")

Function.Exp.Pairs("COL10A1")
Function.Exp.Pairs("MIA")

####Cluster patients by expression (pairs only)###### 
pair.ids<-substr(brca.exp.nb$normal.patients, 1, 12)
cancer.pairs<-brca.exp.nb$cancer.patients[substr(brca.exp.nb$cancer.patients,1,12) %in% pair.ids]

pheatmap(brca.exp.nb$combined.matrices[c("COL10A1","MMP11"),c(brca.exp.nb$normal.patients,cancer.pairs)], trace="none", scale="none",cluster_cols=F)

##Do PCA on expression
PC.MATRIX<-t(brca.exp.nb$combined.matrices[, brca.exp.nb$cancer.patients])
PC.EXP<-prcomp(PC.MATRIX,scale.=T)
plot(PC.EXP, type="l") #First 4 PC seem to capture most of the data (elbow point)

PC.EXP.4<-data.frame(PC.EXP$x[,1:4])
plot(PC.EXP.4, pch=16,col=rgb(0,0,0,0.5) ) #we can see the separation between the 4

##Now that we've simplified expression data into 4 dimensions, let's so k-means clustering on it to cluster patients

Function.K.CLUSTERS.DET<-function(mydata){
  # Determine number of clusters
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers=i,nstart=100, iter.max=1000)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")  
}
Function.K.CLUSTERS.DET(PC.EXP.4) #Around 4 clusters?

PC.EXP.4.CLUSTERS<-kmeans(PC.EXP.4, 4, nstart=25, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(PC.EXP.4, col=PC.EXP.4.CLUSTERS$clust, pch=16)

#3D plot
library(rgl)
plot3d(PC.EXP.4$PC1, PC.EXP.4$PC2, PC.EXP.4$PC3,col=PC.EXP.4.CLUSTERS$cluster)
plot3d(PC.EXP.4$PC1, PC.EXP.4$PC3, PC.EXP.4$PC4,col=PC.EXP.4.CLUSTERS$cluster)

#Look at clustering on expression datahead(PC.EXP.4)
EXP.CLUSTERS<-as.data.table(PC.EXP.4.CLUSTERS$cluster,keep.rownames=T)
setnames(EXP.CLUSTERS, c("SAMPLE", "CLUSTER"))
EXP.CLUSTERS<-EXP.CLUSTERS[order(CLUSTER),]
pheatmap(brca.exp.nb$combined.matrices[sample(1:18185, 100), EXP.CLUSTERS$SAMPLE], scale="none", cluster_cols=F,
         annotation=data.frame(SAMPLE=as.factor(EXP.CLUSTERS$CLUSTER), row.names=EXP.CLUSTERS$SAMPLE))

####The test may seem valid, but perhaps is better to analyze on genes that vary the most from normal#####
BRCA.DIFF.EXP<-Function.RNAseq.Differential.Expression.V2(brca.exp.nb, brca.exp.nb$cancer.patients)

Function.DIFF.EXP.ANALYSIS<-function(exp.matrix, gene){
  
  require(data.table)
  
  normal.exp<-as.vector(exp.matrix$combined.matrices[gene, exp.matrix$normal.patients])
  cancer.exp<-as.vector(exp.matrix$combined.matrices[gene, exp.matrix$cancer.patients])
  
  main.table<-data.table(samples=c(normal.exp, cancer.exp),
                         type=c(rep("normal", length(normal.exp)), rep("cancer", length(cancer.exp))))
  
  ggplot(main.table, aes(type, samples, colour=type)) + geom_boxplot() + geom_jitter() + theme.format
  
}

Function.DIFF.EXP.ANALYSIS(brca.exp.nb, "DMD")

#Now repeat PCA with the top 1000 differentially expressed genes
top.1000.genes<-as.vector(BRCA.DIFF.EXP$ID)[1:14835]
PC.MATRIX<-t(brca.exp.nb$combined.matrices[top.1000.genes , brca.exp.nb$cancer.patients])
PC.EXP<-prcomp(PC.MATRIX,scale.=T,center=T)
plot(PC.EXP, type="l") #First 4 PC seem to capture most of the data (elbow point)

PC.EXP.4<-data.frame(PC.EXP$x[,1:4])
plot(PC.EXP.4, pch=16,col=rgb(0,0,0,0.5) ) #we can see the separation between the 4

Function.K.CLUSTERS.DET(PC.EXP.4) #Around 4 clusters?

PC.EXP.4.CLUSTERS<-kmeans(PC.EXP.4, 4, nstart=100, iter.max=1000)
library(RColorBrewer)
library(scales)
palette(alpha(brewer.pal(9,'Set1'), 0.5))
plot(PC.EXP.4, col=PC.EXP.4.CLUSTERS$clust, pch=16)

EXP.CLUSTERS<-as.data.table(PC.EXP.4.CLUSTERS$cluster,keep.rownames=T)
setnames(EXP.CLUSTERS, c("SAMPLE", "CLUSTER"))
EXP.CLUSTERS<-EXP.CLUSTERS[order(CLUSTER),]
EXP.CLUSTERS<-rbind(EXP.CLUSTERS, data.table(SAMPLE=brca.exp.nb$normal.patients, CLUSTER="NORMAL"))

pheatmap(brca.exp.nb$combined.matrices[top.1000.genes[1:1000], c(EXP.CLUSTERS$SAMPLE, brca.exp.nb$normal.patients)], scale="none", cluster_cols=F,
         annotation=data.frame(SAMPLE=as.factor(EXP.CLUSTERS$CLUSTER), row.names=EXP.CLUSTERS$SAMPLE))

####Test if paired cancers have representative distributions to non-paired#####
Function.cancer.wilcox<-function(exp.obj, pair.ids) {
    
  require(data.table)
  
  #Get both subsets
  paired.cancer<-exp.obj$cancer.patients[substr(exp.obj$cancer.patients, 1,12) %in% pair.ids]
  rest<-setdiff(exp.obj$cancer.patients, paired.cancer)
  
  #Test wilcoxon
  wilcox.vector<-apply(exp.obj$combined.matrices, 1, function(x) wilcox.test(as.vector(x[paired.cancer]), as.vector(x[rest]), paired=F)$p.value)
  
  #Build table
  wilcox.table<-data.table(Hugo_Symbol=rownames(exp.obj$combined.matrices), P.VAL=wilcox.vector)
  
  #Correct for multiple hypothesis
  wilcox.table$P.VAL.ADJ<-p.adjust(wilcox.table$P.VAL, method="fdr")
  
  #Clean up and return
  wilcox.table<-wilcox.table[order(P.VAL.ADJ),]
  return(wilcox.table)
}

repres.test<-Function.cancer.wilcox(brca.exp.nb, pair.ids)
repres.test[P.VAL.ADJ<0.05,]
hist(repres.test$P.VAL.ADJ)

Function.REPRES.TEST<-function(exp.obj, gene.vector, pair.ids){
  
  require(data.table)
  
  #Get both subsets
  paired.cancer<-exp.obj$cancer.patients[substr(exp.obj$cancer.patients, 1,12) %in% pair.ids]
  rest<-setdiff(exp.obj$cancer.patients, paired.cancer)
  
  #Get exp vectors for all genes being tested
  main.list<-list()
  for (gene in gene.vector){
    paired.exp<-as.vector(exp.obj$combined.matrices[gene, paired.cancer])
    rest.exp<-as.vector(exp.obj$combined.matrices[gene, rest])
    normal.exp<-as.vector(exp.obj$combined.matrices[gene, exp.obj$normal.patients])  
    
    main.table<-data.table(samples=c(paired.exp, rest.exp, normal.exp),
                           type=c(rep("paired", length(paired.exp)), rep("rest", length(rest.exp)), rep("normal", length(normal.exp)) ))
    main.table$Hugo_Symbol<-gene
    
    main.list[[gene]]<-main.table
  }
  
  #Combine tables
  main.table<-do.call(rbind, main.list)
  
  #Plot
  ggplot(main.table, aes(type, samples, colour=type)) + geom_boxplot() + geom_jitter(size=1.2) + theme.format + facet_wrap(~Hugo_Symbol)
  
}

Function.REPRES.TEST(brca.exp.nb, c("CAPN10", "WASH3P","WASH7P","ZGPAT", "SPATA2", "PPP4R1",
                                    "ZNF775", "TUBG2", "TGFBR3"), pair.ids)

repres.test[P.VAL.ADJ<0.05,]

Function.RANDOM.WILCOX.PAIRED<-function(exp.obj, target.genes, pair.ids){
  #Test if 10% of the paired-cancer population has a significantly 
  # different distribution than the rest of the population by chance
  
  require(data.table)
  require(parallel)
  
  #Get paired cancer samples
  paired.cancer<-exp.obj$cancer.patients[substr(exp.obj$cancer.patients, 1,12) %in% pair.ids]
  print (length(paired.cancer))
  
  #Reduce exp matrix
  exp.matrix<-exp.obj$combined.matrices[target.genes, paired.cancer]
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix","data.table", "paired.cancer") ,envir=environment())
  print ("Done exporting values")
  
  #Permute
  gene.prob<-parApply(cl,exp.matrix, 1, function(x) {
    
    permut.vector<-replicate(1000, {
      ten<-sample(paired.cancer, 10)
      rest<-setdiff(paired.cancer, ten)
      wilcoxon<-wilcox.test(rep(as.vector(x[ten]),10), rep(as.vector(x[rest]),10),paired=F)$p.value
    })
    permut.vector<-p.adjust(permut.vector, method="fdr")
    permut.result<-mean(permut.vector<0.05)
    return(permut.result)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done with background distributions")
  
  #Construct into table
  main.table<-data.table(Hugo_Symbol=target.genes, PROB=gene.prob)
  
  return(main.table)
}

paired.repres.random.test<-Function.RANDOM.WILCOX.PAIRED(brca.exp.nb, as.vector(repres.test[P.VAL.ADJ<0.05,]$Hugo_Symbol) , pair.ids)
paired.repres.random.test.b<-Function.RANDOM.WILCOX.PAIRED(brca.exp.nb, as.vector(repres.test[P.VAL.ADJ>0.90,]$Hugo_Symbol) , pair.ids)

ggplot(paired.repres.random.test, aes(PROB)) + geom_histogram() + theme.format + 
  ggtitle("Cancer pairs - Expression Permutations comparing 10% of Population by Unpaired Wilcoxon")
ggplot(paired.repres.random.test.b, aes(PROB)) + geom_histogram() + theme.format + 
  ggtitle("Cancer pairs - Expression Permutations comparing 10% of Population by Unpaired Wilcoxon")


###Now
BRCA.POS.LINKAGE.ALL
Function.Exp.Pairs("COL10A1","PIK3CA", 178952085) 
Function.Exp.Pairs("A1BG","PIK3CA", 178952085) 

Function.PAIR.EXP.LM.FIT<-function(exp.obj, target.gene, pair.ids, maf){
  
  require(data.table)
  
  #Get target patients
  paired.cancer<-exp.obj$cancer.patients[substr(exp.obj$cancer.patients, 1,12) %in% pair.ids]
  normal.cancer<-exp.obj$normal.patients
  
  #Simplify matrix
  exp.matrix<-exp.obj$combined.matrices[, c(normal.cancer, paired.cancer)]
  
  #Filter and simplify maf 
  maf<-maf[SAMPLE %in% paired.cancer,]
  maf<-maf[,c("Hugo_Symbol", "Start_Position", "SAMPLE"), with=F]
  setkey(maf)
  maf<-unique(maf)
  
  #Filter for mutations found in at least two samples
  maf[,N.POS.SAMPLES:=length(unique(SAMPLE)), by=c("Hugo_Symbol", "Start_Position")]
  maf<-maf[N.POS.SAMPLES>1,]
  
  #Assign pair ids
  maf$PAIR.IDS<-substr(maf$SAMPLE, 1, 12)
  
  #Split maf
  maf.split<-split(maf, list(maf$Hugo_Symbol, maf$Start_Position), drop=T)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix","data.table", "maf.split", "target.gene","normal.cancer", "paired.cancer") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  print ("Executing parallelization")
  main.list<-parLapply(cl,maf.split, function(x) {
    
    #Get paired samples
    samples<-unique(as.vector(x$SAMPLE))
    paired.normals<-normal.cancer[substr(normal.cancer, 1, 12) %in% unique(x$PAIR.IDS)]
    
    #Get rest
    non.samples<-setdiff(paired.cancer, samples)
    non.paired.normals<-setdiff(normal.cancer, paired.normals)
      
    #First wilcoxon
    wilcox.pval<-wilcox.test(as.vector(exp.matrix[target.gene, samples]), as.vector(exp.matrix[target.gene, paired.normals]), paired=F)$p.value
    wilcox.rest.pval<-wilcox.test(as.vector(exp.matrix[target.gene, non.samples]), as.vector(exp.matrix[target.gene, non.paired.normals]), paired=F)$p.value
    
    #Then linear model
    linear.data<-data.table(TYPE=c(rep(1,length(samples)), rep(0, length(paired.normals))),
                            EXP= c(as.vector(exp.matrix[target.gene, samples]), as.vector(exp.matrix[target.gene, paired.normals])) )
    l.model<-lm(EXP~TYPE, data=linear.data)
    l.ars<-summary(l.model)[["adj.r.squared"]]
    l.slope<-as.vector(l.model$coefficients[2])
    
    linear.data.rest<-data.table(TYPE=c(rep(1,length(non.samples)), rep(0, length(non.paired.normals))),
                            EXP= c(as.vector(exp.matrix[target.gene, non.samples]), as.vector(exp.matrix[target.gene, non.paired.normals])) )
    l.model.rest<-lm(EXP~TYPE, data=linear.data.rest)
    l.ars.rest<-summary(l.model.rest)[["adj.r.squared"]]
    l.slope.rest<-as.vector(l.model.rest$coefficients[2])
    
    #Return
    hugo.table<-data.table(Hugo_Symbol=unique(x$Hugo_Symbol), POS=unique(x$Start_Position), N.SAMPLES=length(samples), 
                           WILCOX.P.VAL=wilcox.pval, LINEAR.ARS=l.ars, LINEAR.SLOPE=l.slope,
                           WILCOX.REST.P.VAL=wilcox.rest.pval, LINEAR.REST.ARS=l.ars.rest, LINEAR.REST.SLOPE=l.slope.rest)
    return(hugo.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Join table and correct for fdr
  main.table<-do.call(rbind, main.list)
  main.table$P.VAL.ADJ<-p.adjust(main.table$WILCOX.P.VAL, method="fdr")
  
  #Clean up and return
  main.table<-main.table[order(main.table$P.VAL.ADJ),]
  return(main.table)
}

LM.WILCOX.COL10A1<-Function.PAIR.EXP.LM.FIT(brca.exp.nb, "COL10A1", pair.ids, brca.maf)
LM.WILCOX.COL10A1[order(LINEAR.ARS, decreasing=T),]

LM.WILCOX.SPRR4<-Function.PAIR.EXP.LM.FIT(brca.exp.nb, "SPRR4", pair.ids, brca.maf)
LM.WILCOX.SPRR4[order(abs(LINEAR.SLOPE), decreasing=T),]
ggplot(LM.WILCOX.SPRR4, aes(LINEAR.SLOPE, LINEAR.REST.ARS)) + geom_point()

Function.Exp.Pairs("SPRR4","MUC20", 195452814)

s<-sample(seq(2,3,0.1),20, replace=T)
x<-data.table(Y=c(s+3, s), TYPE=c(rep(1,length(s)), rep(0, length(s))))
ggplot(x, aes(TYPE, Y)) + geom_point()

y<-lm(Y~TYPE, data=x )
summary(y)[["adj.r.squared"]]
as.vector(y$coefficients[2])

#Get pairwise mutual information on normals
library(entropy)
brca.normal.matrix<-brca.exp.nb$combined.matrices[,brca.exp.nb$normal.patients]

Function.GENE.MATRIX.MI<-function(exp.matrix){
  
  require(entropy)
  require(parallel)
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","exp.matrix", "discretize2d", "mi.empirical") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  genes<-rownames(exp.matrix)
  main.list<-lapply(genes, function(x) {
    
    #Get expression for looping gene
    gene.exp<-as.vector(exp.matrix[x,])
    
    #Apply over all others pairwise
    mi.gene.vector<-parApply(cl, exp.matrix, 1, function(y) {
      
      #Discretize pairs
      gene.disc<-discretize2d(gene.exp, as.vector(y), numBins1=10, numBins2=10)
      mi.pair<-mi.empirical(gene.disc)
      
      #Return pair mi
      return (mi.pair)
      
    })
    
    #Construct table for gene table
    gene.table<-data.table(Hugo.1=x, Hugo.2=rownames(exp.matrix), MI=mi.gene.vector)
    
    #Remove gene from matrix
    exp.matrix<-exp.matrix[setdiff(rownames(exp.matrix), x),]
    
    #Return mi table
    return(gene.table)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")  
  
  #Combine mi tables
  main.table<-do.call(rbind, main.list)
  
  #Clean up and return
  main.table<-main.table[Hugo.1!=Hugo.2, ]
  main.table<-main.table[order(MI, decreasing=T),]
  return(main.table)
}

BRCA.NORMAL.MI<-Function.GENE.MATRIX.MI(brca.normal.matrix[1:1000,])
BRCA.NORMAL.MI[MI>0.8 & Hugo.1=="AQP7",]
cat(BRCA.NORMAL.MI[MI>0.8 & Hugo.1=="AQP7",]$Hugo.2)

pheatmap(brca.normal.matrix[c("APOC3","ADAM7"),], scale="none", trace="none", color=brewer.pal(9,"YlOrRd"))

##### Analyze MI breast cancer normal network ######
BRCA.NORMAL.MI<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/030415.BRCA.EXP.MI.NORMALS", header=T, sep="\t", stringsAsFactors=F)
BRCA.NORMAL.MI

ggplot(BRCA.NORMAL.MI, aes(MI)) + geom_histogram() + theme.format

#Check how many edges we keep at each MI cut
MI.NORMAL.EDGES<-data.table()
for (th in seq(0,max(BRCA.NORMAL.MI$MI) ,0.05)){
  
  mi.edges<-nrow(BRCA.NORMAL.MI[MI>th,])
  MI.NORMAL.EDGES<-rbind(MI.NORMAL.EDGES, data.table(TH=th, N.EDGES=mi.edges))
}

ggplot(MI.NORMAL.EDGES, aes(TH, N.EDGES)) + geom_histogram(stat="identity") + theme.format 

#Check how many genes we keep at each MI cut
MI.NORMAL.GENES<-data.table()
for (th in seq(0,max(BRCA.NORMAL.MI$MI) ,0.05)){
  
  mi.edges<-BRCA.NORMAL.MI[MI>th,]
  mi.genes<-length(unique(c(mi.edges$Hugo.1, mi.edges$Hugo.2)))
  MI.NORMAL.GENES<-rbind(MI.NORMAL.GENES, data.table(TH=th, N.GENES=mi.genes))
}

ggplot(MI.NORMAL.GENES, aes(TH, N.GENES)) + geom_histogram(stat="identity") + theme.format

#Check at the degree distribution at each MI cut
MI.NORMAL.DEGREE.DIST<-data.table()
for (th in seq(0,max(BRCA.NORMAL.MI$MI) ,0.05)){
  print (th)
  
  mi.edges<-BRCA.NORMAL.MI[MI>th,]
  
  mi.edges.1<-mi.edges[,list(DEGREE=length(MI)), by=Hugo.1]  
  setnames(mi.edges.1, c("Hugo", "DEGREE"))
  mi.edges.2<-mi.edges[,list(DEGREE=length(MI)), by=Hugo.2]  
  setnames(mi.edges.2, c("Hugo", "DEGREE"))
  
  mi.edges<-rbind(mi.edges.1, mi.edges.2)
  setkey(mi.edges)
  mi.edges<-unique(mi.edges)
  mi.edges<-mi.edges[,list(DEGREE=sum(DEGREE)), by="Hugo"]
  mi.edges$TH<-th
  
  MI.NORMAL.DEGREE.DIST<-rbind(MI.NORMAL.DEGREE.DIST, mi.edges)
}

MI.NORMAL.DEGREE.DIST.2<-MI.NORMAL.DEGREE.DIST[TH %in% c(seq(0,1.4, 0.1),1.45),] #Without 0.5 threshold points

ggplot(MI.NORMAL.DEGREE.DIST.2, aes(x=DEGREE, fill=factor(TH))) +   
  geom_density(aes(y=..count..), alpha=0.4) + theme.format + coord_trans(x="log1p") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_vline(xintercept=100, color="red", linetype="dashed") 

#Produce networks at each cut
for (th in seq(0.75,1,0.05)){
  print (th)
  
  mi.th<-BRCA.NORMAL.MI[MI>th,]
  filename<-paste("PIPELINES/METABOLIC.DRIVERS/NETWORKS/BRCA.EXP.NORMAL.TH.",th, sep="")
  write.table(file=filename, mi.th, sep="\t", quote=F, row.names=F, col.names=F)
  
}

###Analyzing clusters####
cluster.count<-data.table()
for (th in c("0.75", "0.8", "0.85", "0.9","0.95")){
  test<-scan(paste("PIPELINES/METABOLIC.DRIVERS/NETWORKS/SPICI.RESULTS/BRCA.EXP.NORMAL.TH.",th ,".cluster" ,sep=""), what="",sep="\n")
  test.sizes<-unlist(lapply(test, function(x) length(unlist(strsplit(x,"\t")))))
  cluster.count<-rbind(cluster.count, data.table(CLUSTER=as.numeric(th), SIZES=test.sizes))
}

ggplot(cluster.count, aes(as.factor(CLUSTER), SIZES, colour=CLUSTER)) + geom_boxplot() + geom_jitter() + theme.format +
  coord_trans(y="log1p")

########030615#######

#Test wilcoxon within cancer samples as proof of concept#
brca.within.wilcox<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/030615.BRCA.TEST.WITHIN.CANCER.WILCOX.rds")
brca.within.wilcox.table<-data.table()
for (wil in names(brca.within.wilcox)){
  brca.within.wilcox.table<-rbind(brca.within.wilcox.table, data.table(SIZE=as.numeric(wil), DISTRIBUTION=brca.within.wilcox[[wil]]))
}

ggplot(brca.within.wilcox.table, aes(as.factor(SIZE), DISTRIBUTION, colour=DISTRIBUTION)) + geom_boxplot() + geom_jitter() + theme.format +
  scale_y_sqrt()

####BUILDING FUNCTIONS TO TEST SLE####
library(limSolve)

Function.PREP.MAF<-function(maf, train.patients, POS=T, cnv){
  #Preparing the maf table for solving with expression data as a linear system of equations
  
  require(data.table)
  
  #Filter maf table for patients of interest
  maf<-maf[SAMPLE %in% train.patients,]
  
  #Separate by type of mutation
  non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
  non.func.type<-c("DEL","INS")
  maf$TYPE<-ifelse(maf$Variant_Classification %in% non.func.class, "NON.FUNC", 
                        ifelse(maf$Variant_Type %in% non.func.type, "NON.FUNC", "MISS"))
  
  #Merge hugo and position to give a unique ID to mutation -or- just MUT
  if (POS==T){
    maf$CHANGE<-paste(maf$Hugo_Symbol, maf$Start_Position, maf$TYPE, sep=".")  
  } else {
    maf$CHANGE<-paste(maf$Hugo_Symbol, maf$TYPE, sep=".")  
  }
  
  #Load CNV info and filter for patients of interest
  cnv$COUNT<- 2*(2^cnv$CNV)
  cnv<-cnv[COUNT>2.35 | COUNT<1.65,] #Filter for copy threshold
  cnv$CHANGE<-paste(cnv$Hugo_Symbol, "CNV", sep=".")
  cnv<-cnv[, c("SAMPLE","CHANGE", "COUNT"), with=F]
  cnv<-cnv[SAMPLE %in% train.patients,]
  
  #Simplify maf table
  maf<-maf[,c("SAMPLE", "CHANGE"), with=F]
  
  #Add count to cast maf table
  maf$COUNT<-1
  maf[,POS.SAMPLES:=sum(COUNT), by="CHANGE"]
  
  #Integrate maf and cnv tables to cast
  maf.cast<-rbind(maf[,c("SAMPLE","CHANGE","COUNT"), with=F], cnv)
  
  #Cast table to long format
  maf.cast<-acast(maf.cast, SAMPLE~CHANGE ,fill=0,value.var="COUNT", fun.aggregate=sum)
  
  #Clean up and Return as list
  maf<-maf[,c("CHANGE","POS.SAMPLES"),with=F]
  setkey(maf)
  maf<-unique(maf)
  return(list(CAST=maf.cast, TABLE=maf)) 
}

Function.EXP.RATIO<-function(exp.obj, genes, target.patients){
  require(data.table)
  
  #Extract data
  exp.matrix<-exp.obj$combined.matrices[genes,,drop=F]
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
  #main.table<-merge(main.table, prepped.maf.table[,c("CHANGE","POS.SAMPLES"),with=F], by="CHANGE")
  
  #Return
  return(list(TABLE=main.table, TRAIN.PATIENTS=train.patients))
  
}

##APPLY FUNCTIONS TO SOLVE TRAIN SET###
train.patients<-sample(intersect(intersect(unique(brca.maf$SAMPLE), brca.exp.nb$cancer.patient), unique(BRCA.CNV$main.table$SAMPLE))  ,750) #Take care of intersect when sampling
tested.genes<-c("C20orf46","CDKN3","MED21","MGC16142","OR7D2","PFDN2","PROZ","SILV","TMEM167A")

PREPPED.MAF<-Function.PREP.MAF(brca.maf, train.patients,POS=F, BRCA.CNV$main.table)

EXP.DIFF<-Function.EXP.RATIO(brca.exp.nb, tested.genes , train.patients)

SOLVE.TABLE<-Function.SOLVE.EXP(PREPPED.MAF$CAST, EXP.DIFF, PREPPED.MAF$TABLE, train.patients)
SOLVE.TABLE$TABLE

##ANALYZE FOR COHERENCE OF FUNCTIONS
SOLVE.TABLE.TEST<-melt(SOLVE.TABLE$TABLE, id.vars=c("CHANGE","POS.SAMPLES"))
setnames(SOLVE.TABLE.TEST, c("CHANGE", "POS.SAMPLES","EXP.GENE","COEF"))
ggplot(SOLVE.TABLE.TEST, aes(abs(COEF), POS.SAMPLES, colour=EXP.GENE)) + geom_point() + facet_wrap(~EXP.GENE) + theme.format

Function.TEST.SLE.RESULTS<-function(SOL.TABLE, EXP.DIFF, MAF.CAST) {

  require(data.table)
  
  #Get genes of interest
  genes<-colnames(EXP.DIFF)[2:length(colnames(EXP.DIFF))]
  
  #Order "A" mutation matrix (MAF.CAST) by order to "x' solution mutations in SOL.TABLE and "B" sample order in EXP.DIFF
  MAF.CAST<-MAF.CAST[EXP.DIFF$SAMPLE, SOL.TABLE$CHANGE]
  
  #Get euclidean distance of solution to actual expression values
  euc.vector<-sapply(genes, function (x)  {
    
    A.PRIME<-MAF.CAST%*%SOL.TABLE[[x]]
    A.EUC<-sqrt(sum((EXP.DIFF[[x]]-A.PRIME)^2))
    
    return(A.EUC)
  })
  
  #Build euclidean result table
  main.table<-data.table(Hugo_EXP=genes, EUC=euc.vector)
  
  #Clean up and return
  main.table<-main.table[order(EUC),]
  return(main.table)
}
TEST.SLE.RESULTS<-Function.TEST.SLE.RESULTS(SOLVE.TABLE$TABLE, EXP.DIFF, PREPPED.MAF$CAST)

####TEST ON TEST COHORT####
test.patients<-intersect(intersect(setdiff(unique(brca.maf$SAMPLE), train.patients), brca.exp.nb$cancer.patients), unique(BRCA.CNV$main.table$SAMPLE))

PREPPED.MAF.TEST<-Function.PREP.MAF(brca.maf, test.patients,POS=F,BRCA.CNV$main.table)
EXP.DIFF.TEST<-Function.EXP.RATIO(brca.exp.nb, tested.genes , test.patients)

Function.TEST.SLE.ON.TEST<-function(SOL.TABLE, EXP.DIFF, MAF.CAST) {
  
  require(data.table)
  
  #Get exp genes of interest
  genes<-colnames(EXP.DIFF)[2:length(colnames(EXP.DIFF))]
  
  #Get CHANGE (mutation/position) of interest
  change<-intersect(SOL.TABLE$CHANGE, colnames(MAF.CAST))
  print (length(SOL.TABLE$CHANGE))
  print (length(change))
  
  #Filter sol.table by available "CHANGE" features in solution from train cohort- FILTERING BY CHANGE
  SOL.TABLE<-SOL.TABLE[CHANGE %in% change,]
  
  #Order "A" mutation matrix (MAF.CAST) by order to "x' solution mutations in SOL.TABLE and "B" sample order in EXP.DIFF - FILTERING BY CHANGE
  MAF.CAST<-MAF.CAST[EXP.DIFF$SAMPLE, SOL.TABLE$CHANGE]
  print (sum(rowSums(MAF.CAST)>0))
  
  #Get euclidean distance of solution to actual expression values
  euc.vector<-sapply(genes, function (x)  {
    
    A.PRIME<-MAF.CAST%*%SOL.TABLE[[x]]
    A.EUC<-sqrt(sum((EXP.DIFF[[x]]-A.PRIME)^2))
    #A.EUC<-mean(A.PRIME /EXP.DIFF[[x]])
    
    #Plot
    PLOT<-data.table(cbind(A.PRIME, EXP.DIFF[[x]]))
    setnames(PLOT, c("PREDICTED","ACTUAL"))
    PLOT$PATIENT.ID<-EXP.DIFF$SAMPLE
    PLOT<-melt(PLOT, id.vars="PATIENT.ID")
    setnames(PLOT, c("PATIENT.ID","RESPONSE","EXP.DIFF"))
    
    my.plot<-ggplot(PLOT, aes(RESPONSE, EXP.DIFF, group=PATIENT.ID, colour=PATIENT.ID)) + geom_point(size=4, shape=21, fill="white") + 
     geom_line() + theme.format + ggtitle(x)
    
    filename=paste("FIGURES/SLE/GENE.FUNC.CNV/","750.TRAIN.",x,".jpeg", sep="")
    jpeg(filename=filename ,width=1200,height=800, quality=100,type="quartz")
    print (my.plot)
    dev.off()
    
    return(A.EUC)
  })
  
  #Build euclidean result table
  main.table<-data.table(Hugo_EXP=genes, EUC=euc.vector)
  
  #Clean up and return
  main.table<-main.table[order(EUC),]
  return(main.table)
}

length(test.patients)
TEST.SLE.TEST.POS<-Function.TEST.SLE.ON.TEST(SOLVE.TABLE$TABLE, EXP.DIFF.TEST, PREPPED.MAF.TEST$CAST) #At position
TEST.SLE.TEST.MUT<-Function.TEST.SLE.ON.TEST(SOLVE.TABLE$TABLE, EXP.DIFF.TEST, PREPPED.MAF.TEST$CAST) #At gene

ggplot(TEST.SLE.TEST.MUT, aes(Hugo_EXP, EUC)) + geom_histogram(stat="identity") + theme.format

BRCA.DIFF.EXP[ID=="CDKN3",] #WHY ARE SOME GENES HARDER TO PREDICT THAN OTHERS?, WHAT SO SPECIAL ABOUT THOSE WERE WE HAVE A PROPORTION OF WIDELY UNMATCHING VALUES?
brca.maf[Hugo_Symbol=="CDKN3",]

#Compare at position and mutation
TEST.SLE.TEST.POS$TYPE<-"POS"
TEST.SLE.TEST.MUT$TYPE<-"MUT"

ggplot(rbind(TEST.SLE.TEST.MUT, TEST.SLE.TEST.POS), aes(Hugo_EXP, EUC, colour=TYPE)) + geom_histogram(stat="identity",position="dodge") + theme.format

####ADDING CNV DATA####
BRCA.CNV<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/021915.BRCA.PROCESSED.CNV.rds")
BRCA.CNV$main.table$CN<- 2*(2^BRCA.CNV$main.table$CNV)
BRCA.CNV$main.table<-BRCA.CNV$main.table[CN>2.46 | CN<1.63,] #Choose as thresholds 0.3 and 0.9 (in terms of CNV)
BRCA.CNV$main.table$CNV.CAT<-ifelse( (BRCA.CNV$main.table$CN>2.46 & BRCA.CNV$main.table$CN<3.73), 1, 
                                       ifelse(BRCA.CNV$main.table$CN>=3.73, 2,
                                              ifelse( (BRCA.CNV$main.table$CN>1.07 & BRCA.CNV$main.table$CN<1.63), -1,
                                                      ifelse(BRCA.CNV$main.table$CN<=1.07, -2,0)))  )
BRCA.CNV$main.table$CNV.CLASS<-paste(BRCA.CNV$main.table$CNV.CAT, BRCA.CNV$main.table$Hugo_Symbol, sep=".")
BRCA.CNV$main.table<-BRCA.CNV$main.table[SAMPLE %in% BRCA.CNV$cancer.patients,] #Filter by cancer samples only

BRCA.CNV$main.table[,CLASS.POP:=length(SAMPLE), by="CNV.CLASS"]
length(unique(BRCA.CNV$main.table[CNV.CAT==2 | CNV.CAT==-2,][CLASS.POP>20,]$CNV.CLASS)) #TOO MANY!! PERHAPS DIVIDE INTO EXTRA CLASS OF HYPER COPY NUMBER i.e. 3, 4

#Find the common groups
BRCA.CNV.GROUPS<-BRCA.CNV$main.table[CLASS.POP>2,]


BRCA.CNV$main.table<-BRCA.CNV$main.table[order(CNV, decreasing=T),]
ggplot(BRCA.CNV$main.table, aes(CN)) + geom_histogram() + coord_trans(y="log1p")

#Load GISTIC
BRCA.GISTIC<-fread("DATABASES/CANCER_DATA/TCGA/CNV/BRCA/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2014041600.0.0/all_thresholded.by_genes.txt", 
                   header=T, sep="\t",stringsAsFactors=F)
setnames(BRCA.GISTIC, c("Hugo_Symbol","Locus_ID","Cytoband", colnames(BRCA.GISTIC)[4:ncol(BRCA.GISTIC)]  ))

BRCA.GISTIC[Hugo_Symbol %in% c("NAALADL2","TP53"),c("Hugo_Symbol", "TCGA-D8-A1JU-01A-11D-A13J-01"), with=F] #0
BRCA.GISTIC[Hugo_Symbol %in% c("LPAL2","TP53"),c("Hugo_Symbol", "TCGA-BH-A0BW-01A-11D-A111-01"), with=F] #-1
BRCA.GISTIC[Hugo_Symbol %in% c("LPAL2","TP53"),c("Hugo_Symbol", "TCGA-S3-A6ZH-01A-22D-A32H-01"), with=F] #-1
BRCA.GISTIC[Hugo_Symbol %in% c("FGFR2","TP53"),c("Hugo_Symbol", "TCGA-D8-A142-01A-11D-A111-01"), with=F] #2
BRCA.GISTIC[Hugo_Symbol %in% c("BCAS3","TP53"),c("Hugo_Symbol", "TCGA-A2-A0CX-01A-21D-A011-01"), with=F] 

BRCA.CNV$main.table[Hugo_Symbol %in% c("BCAS3","TP53") & SAMPLE=="TCGA.A2.A0CX.01A",]

#####Analyze maf####
brca.maf[,MUT.COUNT:=length(TYPE), by=c("SAMPLE")]
brca.maf<-brca.maf[order(MUT.COUNT, decreasing=T),]
brca.maf[SAMPLE=="TCGA.A7.A26I.01B",] #Filer out 01B???
brca.maf[SAMPLE=="TCGA.AC.A2FK.01A",] #Why do you have cancer?
brca.maf[Hugo_Symbol=="PEX5",]

hist(BRCA.CNV$main.table[SAMPLE=="TCGA.AC.A2FK.01A",]$CN)