#####Develop function that describes correlation between MAF and something else in 1000G#####

library(data.table)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(stats)
library(GMD)
heatmap.3<-heatmap.2

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

Function.Prep.MAF<-function(maf.file) {
  
  require(car)
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode",
              "Match_Norm_Seq_Allele1", "Tumor_Seq_Allele2"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #Remove silent mutations
  #maf<-maf[Variant_Classification!="Silent",]
  
  #Equalize REF-ALT pairs
  maf$PRE.REF.ALT<-ifelse(maf$Variant_Type=="SNP", paste(as.vector(maf$Match_Norm_Seq_Allele1), as.vector(maf$Tumor_Seq_Allele2), sep="_"), "-")
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
  maf$CLASS<-ifelse((maf$TYPE=="MISS" | maf$TYPE=="SILENT"), paste(maf$Hugo_Symbol, maf$Start_Position, maf$TYPE, sep="."), paste(maf$Hugo_Symbol, maf$TYPE, sep="."))
  
  #Count how many samples are covered by a type of mutation
  maf[,POP.CLASS:=length(SAMPLE),by="CLASS"]
  
  #Cleand up and Return
  setkey(maf)
  maf<-unique(maf)
  maf<-maf[order(POP.CLASS, decreasing=T),]
  return(maf)
}
############################################################################################################################################
MAF.ANNOVAR<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/ANNOVAR/010715.THOUSAND.SNP.TRIMERS.PLUS.PHAST.45.avinput.exonic_variant_function", drop=c(1,6),
                   sep="\t", header=F, stringsAsFactors=F)
MAF.ANNOVAR<-Function.PROCESS.ANNOVAR(MAF.ANNOVAR)
MAF.ANNOVAR.CLASS<-Function.PREP.MAF.CLASS.ML(MAF.ANNOVAR, CHEN.REP, EXON.NT, c("MAF"), FILTER=T)
MAF.ANNOVAR.CLASS$NT.CUT<-cut(MAF.ANNOVAR.CLASS$NT, quantile(MAF.ANNOVAR.CLASS$NT, seq(0,1,0.20)), include.lowest=T)
MAF.ANNOVAR.CLASS$REP.CLASS<-factor(MAF.ANNOVAR.CLASS$REP.CLASS, levels=c("low", "medium", "high"))
MAF.ANNOVAR.CLASS$PHAST.CLASS<-cut(MAF.ANNOVAR.CLASS$PHAST, quantile(MAF.ANNOVAR.CLASS$PHAST, c(0,0.33,0.66,1)),include.lowest=T, 
                                   labels= as.factor(c("low", "medium", "high")))

ggplot(MAF.ANNOVAR.CLASS, aes(REF.ALT, MAF, colour=REF.ALT)) + geom_boxplot() + geom_jitter(size=0.2) + theme.format +
  facet_wrap(~TYPE) + scale_y_log10()

#Remove more noisy ref_alt (A_G and C_T) and try to correlate to phastcons or reptime
ggplot(MAF.ANNOVAR.CLASS[!(REF.ALT %in% c("A_G","C_T")),], aes(REP.CLASS, MAF, colour=REF.ALT)) + geom_boxplot() +geom_jitter(size=0.2)+theme.format +
  facet_wrap(~TYPE)+ scale_y_log10()

ggplot(MAF.ANNOVAR.CLASS, aes(REF.ALT, MAF, colour=REP.CLASS)) + geom_boxplot() +geom_jitter(size=0.2)+theme.format +
  facet_wrap(~TYPE)+ scale_y_log10() #FOUND: Slight lower maf with respect to rep.time for some ref_alt pairs for non-syn and opposite for syn(for some pairs)!

#Phastcons to NT
ggplot(MAF.ANNOVAR.CLASS[!(REF.ALT %in% c("A_G","C_T")),], aes(NT, MAF, colour=REF.ALT)) + geom_point(size=0.5) +theme.format +
  facet_wrap(REF.ALT~TYPE) + coord_trans(x="log1p")

#Include noisy, so pattern in NT vs MAF 
ggplot(MAF.ANNOVAR.CLASS, aes(NT, MAF, colour=REF.ALT)) + geom_point(size=0.5) +theme.format +
  facet_grid(REF.ALT~TYPE) + coord_trans(x="log1p")

#NOT FOUND no correlation of NT vs MAF
ggplot(MAF.ANNOVAR.CLASS, aes(NT.CUT, MAF, colour=REF.ALT)) + geom_boxplot() + geom_jitter(size=0.5) +theme.format +
  facet_grid(REF.ALT~TYPE) + scale_y_log10()

#Look for correlation between phastcons and maf - FOUND: Negative correlation of MAF with respect to PHASTCons
ggplot(MAF.ANNOVAR.CLASS, aes(PHAST, MAF, colour=REF.ALT)) + geom_point(size=0.3) + theme.format + facet_grid(REF.ALT~TYPE)

ggplot(MAF.ANNOVAR.CLASS, aes(PHAST.CLASS, MAF, colour=REF.ALT)) + geom_boxplot() + geom_jitter(size=0.25) + theme.format + 
  facet_grid(REF.ALT~TYPE) + scale_y_log10()

ggplot(MAF.ANNOVAR.CLASS[TRIMER %in% c("AAG", "GAC"), ], aes(PHAST.CLASS, MAF, colour=REF.ALT)) + geom_boxplot() + geom_jitter(size=0.25) + theme.format + 
  facet_grid(REF.ALT~TYPE) + scale_y_log10()


################################################################################################################################
BRCA.MAF[SAMPLE=="TCGA.A1.A0SH.01A",]###Check on where sample came from (conside mergin mafs...)

##NOTE: Write step to remove pseudogenes!!!
##NOTE: Use MUTSIG File!!!
##NOTE: Wouldn't we have to normalize by CNV?? (For calculation in TCGA mut frequency, that is we might be having the mutation more than once in an individual)
Function.Main<-function(maf, annovar, CANCER.PROB=0.12){
  
  #####Process annovar table##### 
  #Get minimum non-zero maf for 1000G
  min.annovar.maf<-min(annovar$MAF[annovar$MAF!=0])
  
  #Filter out zero mafs
  annovar<-annovar[MAF!=0,]
  
  #First filter silent mutations form annovar file
  annovar<-annovar[TYPE=="synonymous SNV",]
  
  #Classify mutations in annovar table
  annovar$TYPE<-ifelse(annovar$TYPE=="synonymous SNV", "SILENT",
    ifelse(annovar$TYPE=="nonsynonymous SNV", "MISS", "NON.FUNC"))
  
  #Construct mutations classes per gene
  annovar$CLASS<-ifelse((annovar$TYPE=="NON.FUNC" | annovar$TYPE=="SILENT"), paste(annovar$Hugo_Symbol, "NON.FUNC", sep="."),
                        paste(annovar$Hugo_Symbol, annovar$Position, "MISS", sep="."))
  
  #Get total frequencies mutation class
  annovar<-annovar[,list(CLASS.FREQ.NORMAL=sum(MAF)), by="CLASS"]
  
  #####Process cancer maf table#####
  #Filter for silent
  maf<-maf[TYPE=="SILENT",]
  
  #Get total cancer population
  total.pop<-length(unique(maf$SAMPLE))
  
  #Get total frequencies per mutations class
  maf<-maf[,c("CLASS", "POP.CLASS"),with=F]
  setkey(maf)
  maf<-unique(maf)
  maf$CLASS.FREQ.CANCER<-maf$POP.CLASS/total.pop
  maf$POP.CLASS<-NULL
  
  #Get minimum population size to normalize bayesian calculation for cancer sample size
  min.cancer.maf<-min(maf$CLASS.FREQ.CANCER)
  
  #####Combine tables######
  #Merge per site and add minimum annovar maf for non-matches
  main.table<-merge(maf, annovar,by="CLASS", all.x=T)
  
  #Fill those without matches in 1000G with min.maf
  main.table$CLASS.FREQ.NORMAL[is.na(main.table$CLASS.FREQ.NORMAL)]<-min.annovar.maf
  
  #Normalize normal population frequency for cancer sample size minimum threshold ?????- Perhaps normalize at each gene level????
  #main.table$CLASS.FREQ.NORMAL<-main.table$CLASS.FREQ.NORMAL+min.cancer.maf
  
  ######Apply bayes approach######
  main.table$BAYES.PROB<-main.table$CLASS.FREQ.CANCER*CANCER.PROB/(main.table$CLASS.FREQ.CANCER*CANCER.PROB + main.table$CLASS.FREQ.NORMAL*(1-CANCER.PROB))
  
  ######Clean up and return#####
  main.table<-main.table[order(BAYES.PROB, decreasing=T),]
  return(main.table)
}

###BREAST
BRCA.BAYES<-Function.Main(BRCA.MAF, MAF.ANNOVAR, CANCER.PROB=0.12)
hist(BRCA.BAYES$BAYES.PROB)
BRCA.BAYES[BAYES.PROB>0.80,]
BRCA.BAYES[grepl("PIK3CA",CLASS),]

###GBM 
GBM.BAYES<-Function.Main(GBM.MAF, MAF.ANNOVAR, CANCER.PROB=0.10)
GBM.BAYES[BAYES.PROB>0.60,]
GBM.BAYES[grepl("PDGFRA",CLASS),]

###OV
OV.BAYES<-Function.Main(OV.MAF, MAF.ANNOVAR, CANCER.PROB=0.05)
hist(OV.BAYES$BAYES.PROB)
OV.BAYES[BAYES.PROB>0.75,]
OV.BAYES[grepl("RB1",CLASS),]
###NOTE:When testing against CCG, Filter first for genes that we actually have in our set (i.e. GOPC, IDH2) *******

########################################################################################################################################
MAF.ANNOVAR
ggplot(MAF.ANNOVAR[Chrom==1 & TYPE=="nonsynonymous SNV",], aes(Position, MAF, colour=TYPE)) + geom_point(size=0.5) + theme.format +
  facet_wrap(~TYPE) + coord_trans(y="sqrt")

ggplot(MAF.ANNOVAR[Hugo_Symbol=="NEB" & TYPE %in% c("nonsynonymous SNV", "synonymous SNV"),], aes(Position, MAF, colour=REF.ALT)) + geom_point() + theme.format +
  facet_grid(TYPE~REF.ALT) + coord_trans(y="sqrt")

EXON.NT[order(NT),]

Function.ANNOVAR.REF.ALT<-function(annovar, filter=0){
  
  #Filter for synonymous and non-synomous mutations
  annovar<-annovar[TYPE %in% c("nonsynonymous SNV", "synonymous SNV"),]
  
  #Split by types
  annovar.syn<-annovar[TYPE=="synonymous SNV",]
  annovar.nonsyn<-annovar[TYPE=="nonsynonymous SNV",]
  
  #Calculate type frequencies
  annovar.syn[,HUGO.REF.ALT:=sum(MAF), by="Hugo_Symbol"]
  annovar.nonsyn[,HUGO.REF.ALT:=sum(MAF), by="Hugo_Symbol"]
  
  annovar.syn<-annovar.syn[,list(REF.ALT.HUGO=sum(MAF)/unique(HUGO.REF.ALT)), by=c("Hugo_Symbol", "REF.ALT")]
  annovar.syn$REF.ALT<-paste(annovar.syn$REF.ALT, "syn", sep=".")
  annovar.nonsyn<-annovar.nonsyn[,list(REF.ALT.HUGO=sum(MAF)/unique(HUGO.REF.ALT)), by=c("Hugo_Symbol", "REF.ALT")]
  annovar.nonsyn$REF.ALT<-paste(annovar.nonsyn$REF.ALT, "nonsyn", sep=".")
  
  annovar<-rbind(annovar.nonsyn, annovar.syn)
  
  #Cast table
  annovar<-acast(annovar, Hugo_Symbol~REF.ALT, fill=0,value.var="REF.ALT.HUGO")
  
  #Apply filter?
  annovar<-annovar[apply(annovar, 1, function(x) sum(x>0)>filter),]
  
  #Return
  return (annovar)
}

Function.REF.ALT.PLOT<-function(maf.annovar.ref.alt, gene){
  
  #ref.alt names
  ref.alt<-sapply(colnames(maf.annovar.ref.alt)[seq(1,11,2)], function(x) unlist(strsplit(x, "[.]"))[1])
  
  #values
  non.syn<-as.vector(maf.annovar.ref.alt[gene, seq(1,11,2)])
  syn<-as.vector(maf.annovar.ref.alt[gene, seq(2,12,2)])
  
  #Construct matrix
  main.matrix<-matrix(c(non.syn, syn), nrow=2,byrow=T)
  colnames(main.matrix)<-ref.alt
  rownames(main.matrix)<-c("NON.SYN", "SYN")
  
  #Return
  return(main.matrix)
}

MAF.ANNOVAR.REF.ALT<-Function.ANNOVAR.REF.ALT(MAF.ANNOVAR,filter=0)
head(MAF.ANNOVAR.REF.ALT)
heatmap(MAF.ANNOVAR.REF.ALT, scale="none")#REDO
pheatmap(cor(MAF.ANNOVAR.REF.ALT,method="spearman"), scale="none", trace="none", color=brewer.pal(9,"RdYlBu"))

heatmap.2(Function.REF.ALT.PLOT(MAF.ANNOVAR.REF.ALT, "BRCA1"), scale="none", trace="none")

Function.MAF.REF.ALT<-function(cancer.maf, filter=0){
  
  #Filter for synonymous and non-synomous mutations
  maf<-cancer.maf[TYPE %in% c("MISS", "SILENT"),]
  
  #Obtain maf
  maf$MAF<-maf$POP.CLASS/length(unique(maf$SAMPLE))
  
  #Split by types
  maf.syn<-maf[TYPE=="SILENT",]
  maf.nonsyn<-maf[TYPE=="MISS",]
  
  #Calculate type frequencies
  maf.syn[,HUGO.REF.ALT:=sum(MAF), by="Hugo_Symbol"]
  maf.nonsyn[,HUGO.REF.ALT:=sum(MAF), by="Hugo_Symbol"]
  
  maf.syn<-maf.syn[,list(REF.ALT.HUGO=sum(MAF)/unique(HUGO.REF.ALT)), by=c("Hugo_Symbol", "REF.ALT")]
  maf.syn$REF.ALT<-paste(maf.syn$REF.ALT, "syn", sep=".")
  maf.nonsyn<-maf.nonsyn[,list(REF.ALT.HUGO=sum(MAF)/unique(HUGO.REF.ALT)), by=c("Hugo_Symbol", "REF.ALT")]
  maf.nonsyn$REF.ALT<-paste(maf.nonsyn$REF.ALT, "nonsyn", sep=".")
  
  maf<-rbind(maf.nonsyn, maf.syn)
  
  #Cast table
  maf<-acast(maf, Hugo_Symbol~REF.ALT, fill=0,value.var="REF.ALT.HUGO")
  
  #Apply filter?
  maf<-maf[apply(maf, 1, function(x) sum(x>0)>filter),]
  
  #Return
  return (maf)
}

BRCA.MAF.REF.ALT<-Function.MAF.REF.ALT(BRCA.MAF, filter=0)
head(BRCA.MAF.REF.ALT)
dim(BRCA.MAF.REF.ALT)
heatmap.2(BRCA.MAF.REF.ALT, scale="none", trace="none")

#Compare and contrast syn in 1000G vs cancer and also for nonsyn

Function.Compare.REF.ALT<-function(annovar.ref.alt, cancer.ref.alt, filter=0){
  
  ########Work on 1000G first##############
  #annovar.nonsyn<-annovar.ref.alt[,grep(".syn",colnames(annovar.ref.alt), fixed=T)]
  annovar.nonsyn<-annovar.ref.alt[,grep("nonsyn",colnames(annovar.ref.alt), fixed=T)]
  colnames(annovar.nonsyn)<-sapply(colnames(annovar.nonsyn), function(x)  unlist(strsplit(x,"[.]"))[1])
  
  ######Then on cancer ######
  #cancer.nonsyn<-cancer.ref.alt[,grep(".syn",colnames(cancer.ref.alt), fixed=T)]
  cancer.nonsyn<-cancer.ref.alt[,grep("nonsyn",colnames(cancer.ref.alt), fixed=T)]
  colnames(cancer.nonsyn)<-sapply(colnames(cancer.nonsyn), function(x)  unlist(strsplit(x,"[.]"))[1])
  
  #Filter for at least n counts per row
  cancer.nonsyn<-cancer.nonsyn[apply(cancer.nonsyn, 1, function(x) sum(x>0)>filter),]
  
  #####Filter for common genes####
  genes<-intersect(rownames(cancer.nonsyn), rownames(annovar.nonsyn))
  cancer.nonsyn<-cancer.nonsyn[genes,]
  annovar.nonsyn<-annovar.nonsyn[genes,]
  
  ####Build common dendogram based on background(1000G)######
  dendo.r<-hclust(dist(annovar.nonsyn))
  dendo.h<-hclust(dist(t(annovar.nonsyn)))
  
  ######Return as list############
  return(list(cancer=cancer.nonsyn, annovar=annovar.nonsyn, R=dendo.r, C=dendo.h))
}

CANCER.COMPARE.REF.ALT<-Function.Compare.REF.ALT(MAF.ANNOVAR.REF.ALT, BRCA.MAF.REF.ALT, filter=2)
dim(CANCER.COMPARE.REF.ALT$cancer)

layout(rbind(c(8,3,7,6),c(2,1,5,4)), respect = F)
heatmap.3(CANCER.COMPARE.REF.ALT$annovar, scale="none", Rowv=as.dendrogram(CANCER.COMPARE.REF.ALT$R), Colv= as.dendrogram(CANCER.COMPARE.REF.ALT$C),
          main="1000G", trace="none", colsep=c(2,4))
heatmap.3(CANCER.COMPARE.REF.ALT$cancer, scale="none", Rowv=as.dendrogram(CANCER.COMPARE.REF.ALT$R), Colv= as.dendrogram(CANCER.COMPARE.REF.ALT$C),
          main="Breast Cancer", trace="none", colsep=c(2,4))

######Compare using dot product functions between normal and cancer#######
head(CANCER.COMPARE.REF.ALT$cancer)

Function.COS.REF.ALT<-function(compare.ref.alt.cancer, compare.ref.alt.normal){
  
  cos.theta<-function(a,b){
    m<-sum(a*b)/(sqrt(sum(a^2)) * sqrt(sum(b^2)))
    return(m)
  }
  
  ####Obtain comparable gene and ref alt classes####
  genes<-intersect(rownames(compare.ref.alt.cancer), rownames(compare.ref.alt.normal))
  ref.alt<-intersect(colnames(compare.ref.alt.cancer), colnames(compare.ref.alt.normal))
  l.ref.alt<-length(ref.alt)
  
  ####Order matrices and bind####
  cancer<-compare.ref.alt.cancer[genes, ref.alt]
  normal<-compare.ref.alt.normal[genes, ref.alt]
  main<-cbind(cancer,normal)
  
  ####Calculate cos(theta) per gene#####
  cos.thetas<-apply(main, 1, function(x) cos.theta( (x[1:l.ref.alt]+x[(l.ref.alt+1):(2*l.ref.alt)]), x[(l.ref.alt+1):(2*l.ref.alt)]))
  main.table<-data.table(Hugo_Symbol=genes, COS.THETA=cos.thetas)
  
  ####Clean up and return#####
  main.table<-main.table[order(abs(COS.THETA), decreasing=T),]
  return(main.table)
}

COS.THETA.BRCA.ANNOVAR.NONSYN<-Function.COS.REF.ALT(CANCER.COMPARE.REF.ALT$cancer, CANCER.COMPARE.REF.ALT$annovar)
hist(COS.THETA.BRCA.ANNOVAR.NONSYN$COS.THETA)
COS.THETA.BRCA.ANNOVAR.NONSYN[Hugo_Symbol=="TP53",]

CANCER.COMPARE.REF.ALT$cancer[c("TP53", "TTN", "PIK3CA"),]
CANCER.COMPARE.REF.ALT$annovar[c("TP53", "TTN", "PIK3CA"),]

MAF.ANNOVAR[Hugo_Symbol=="AP2M1",]
BRCA.MAF[Hugo_Symbol=="AP2M1",]

############################################################################
######Predict nonsyn ratio using syn ratio######
test<-copy(MAF.ANNOVAR.REF.ALT)
head(test)
dim(test)
#First filter noisy data, that is, those nonsyn that don't have at least 3 non-zero frequencies per gene
test<-test[apply(test, 1, function(x) sum(x[grepl("nonsyn", colnames(test))]>0)>=3), ]

#Then let's make sure we have enough predictive power, that is, keep those syn that have at least 2 non-zero frequencies per gene
test<-test[apply(test, 1, function(x) sum(x[grepl(".syn", colnames(test), fixed=T)]>0)>=3), ]

#Now we are ready to predict, let's try the one with most nonsyn information first
colSums(test)
test<-data.frame(test)

#This appears to be C_T, let's try that first
test.CT<-test[,c("A_C.syn","A_G.syn","A_T.syn","C_A.syn","C_G.syn","C_T.syn","C_T.nonsyn")]

#Predict first using linear model
test.CT.predict<-lm(C_T.nonsyn~., test.CT)
summary(test.CT.predict)

#Then deeplearning using h2o
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '10g', nthreads=-1) 
h2o_test <- as.h2o(localH2O, test.CT, key = 'test.CT') 

#Break into train and test set
set.seed(1234)
ALL_ROWS<-1:nrow(h2o_test)
RAND_FOLDS<-createFolds(ALL_ROWS,5)
TRAIN_ROWS<-unlist(RAND_FOLDS[1:3])
TEST_ROWS<-unlist(RAND_FOLDS[4:5])
TRAIN_CT<-test.CT$C_T.nonsyn[TRAIN_ROWS] #Train with 60%
TEST_CT<-test.CT$C_T.nonsyn[TEST_ROWS] #Test on 40%

#Model without dropout
CT.MODEL<-h2o.deeplearning(x=1:6, y=7, data=h2o_test[TRAIN_ROWS,], classification=F,
                            activation = "MaxoutWithDropout", balance_classes = TRUE, hidden = c(650,650), epochs = 500,
                           input_dropout_ratio = 0.3,hidden_dropout_ratios = c(0.5,0.5))

#Evaluate 
ct_train <- h2o.predict(CT.MODEL, h2o_test[TRAIN_ROWS,])$predict
ct_train <- as.matrix(ct_train)
ct_test <- h2o.predict(CT.MODEL, h2o_test[TEST_ROWS, ])$predict
ct_test <- as.matrix(ct_test)

CT.MODEL.1000.3.100.50.5.E.300<-copy(CT.MODEL) #88.49, -4.55
CT.MODEL.1000.2.200.100.E.500<-copy(CT.MODEL)  #83.31, -3.32
CT.MODEL.1000.2.400.200.E.500<-copy(CT.MODEL)  #85.69, -2.47
CT.MODEL.5000.2.500.100.E.200<-copy(CT.MODEL)  #85.69, -2.47

mean((ct_train-TRAIN_CT)^2)
1-(sum((ct_train-TRAIN_CT)^2) / sum((ct_train-mean(ct_train))^2))
mean((ct_test-TEST_CT)^2)

list.models<-list(c(100,100,100),c(100,100,300),c(100,100,500),
                  c(200,200,100),c(200,200,300),c(200,200,500),
                  c(1000,100,10,100), c(1000,100,10,300), c(1000,100,10,500),
                  c(1000,200,10,100), c(1000,200,10,300), c(1000,200,10,500),
                  c(5000,500,50,100), c(5000,500,50,300), c(5000,500,50,500))

for (n in list.models){
  n.length<-length(n)
  n.h<-n[1:(n.length-1)]
  n.e<-n[n.length]
  
    CT.MODEL<-h2o.deeplearning(x=1:6, y=7, data=h2o_test[TRAIN_ROWS,], classification=F,
                             activation = "Tanh", balance_classes = TRUE, hidden = n.h, epochs = n.e)
  
  par.name<-paste(n.h, collapse=".")
  model.name<-paste("MODEL", par.name, "E", n.e,"h2o", sep=".")
  
  h2o.saveModel(CT.MODEL, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/REF.ALT.PREDICT/",name=model.name)
}

m.folder<-"PIPELINES/METABOLIC.DRIVERS/OBJECTS/1000G.PREDICTION/H2O/REF.ALT.PREDICT/"
MAIN.TABLE<-data.table()
for (m in list.files(m.folder)){
  m.size<-as.numeric(unlist(strsplit(m, "[.]"))[2])
  print(m.size)
  CT.MODEL<-h2o.loadModel(localH2O, paste(m.folder, m, sep=""))
  
  #Evaluate 
  ct_train <- h2o.predict(CT.MODEL, h2o_test[TRAIN_ROWS[1:m.size],])$predict
  ct_train <- as.matrix(ct_train)
  ct_test <- h2o.predict(CT.MODEL, h2o_test[TEST_ROWS, ])$predict
  ct_test <- as.matrix(ct_test)
    
  LS.TRAIN<-sum((ct_train-TRAIN_CT[1:m.size])^2)
  MSE.TRAIN<-mean((ct_train-TRAIN_CT[1:m.size])^2)
  
  LS.TEST<-sum((ct_test-TEST_CT)^2)
  MSE.TEST<-mean((ct_test-TEST_CT)^2)
  
  #Add to table
  MAIN.TABLE<-rbind(MAIN.TABLE, data.table(MODEL=m, SIZE=m.size, LS.TRAIN=LS.TRAIN, MSE.TRAIN=MSE.TRAIN,
                                           LS.TEST=LS.TEST, MSE.TEST=MSE.TEST))
}
  
######NEXT ATTEMPT######
#Can we predict how variable a gene nonsyn is in general based on syn
test<-MAF.ANNOVAR[TYPE %in% c("nonsynonymous SNV", "synonymous SNV"),][,list(HUGO.MAF.TYPE=sum(MAF)),by=c("Hugo_Symbol", "TYPE")]
test<-acast(test, Hugo_Symbol~TYPE, value.var="HUGO.MAF.TYPE", fill=0,drop=T)
test<-data.table(test, keep.rownames=T)
setnames(test, c("Hugo_Symbol", "syn", "nonsyn"))
ggplot(test[nonsyn!=0 & syn!=0,][nonsyn<5 & syn<5,], aes(syn, nonsyn)) + theme.format + geom_point() + coord_trans(x="log1p", y="log1p")
summary(lm(nonsyn~syn, test[nonsyn!=0 & syn!=0,]))

test<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/ANNOVAR/010715.THOUSAND.PHAST.45.SIFT.hg19_ljb26_sift_filtered", header=F, sep="\t", stringsAsFactors=F)
setnames(test, c("Chrom", "Start", "End", "REF", "ALT", "TRIMER", "MT", "EXON", "MAF", "SCORE"))
table(test$EXON)
ggplot(test[EXON==TRUE,], aes(SCORE, MAF, colour=EXON)) + geom_point() + theme.format + facet_wrap(~EXON)
