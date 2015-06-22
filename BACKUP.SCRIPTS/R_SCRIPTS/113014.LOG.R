#TESTING CORRELATION BETWEEN MAF AND PHASTCON SCORES
#120814
# NOT WORKING, Tried phastcons 100, 45 and 45 placental, none show correlation. Tried by filtering read depth and quality of reads in 1000genome dataset, doesn't work either. Tried testing at each gene level, tested a handful including COSMIC breast cancer genes, some show very weak correlation to none, those that have some correlation do not have enough data to do a fit. It is not worth it to continue at this point, will proceed with implementing and improving bidirectional approach to estimate MAF (bakcground variation from 1000G) on sites that we do not have data on by incorporation information on gene length, amount of data we have (sites of data to gene length ratio) and the MEH factor (hypermutability of a patient, this should influece real background mutation rate of individual)

require(data.table)
library(ggplot2)

#Load Exon coordiantes
EXONS<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES", header=T, sep="\t", stringsAsFactors=F)
EXONS[Hugo_Symbol=="TP53",]

#Load 1000G data
THOUSAND<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.SNP.CNV", header=T, sep="\t", stringsAsFactors=F, drop=4:10)
THOUSAND$MAF<-ifelse(THOUSAND$AF>0.5, 1-THOUSAND$AF,THOUSAND$AF)

THOUSAND.Y<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.Y.SNP.CNV",header=T, sep="\t", stringsAsFactors=F, drop=4:10)
THOUSAND.Y$MAF<-ifelse(THOUSAND.Y$AF>0.5, 1-THOUSAND.Y$AF,THOUSAND.Y$AF)

#Load Phastcons (UPDATED 120815 with PHASTCONS.46 version)
PHAST<-fread("DATABASES/PHASTCONS/120314.CHRM.FILTERED.EXONS.45.PARALLEL", header=F, sep="\t", stringsAsFactors=F)
setnames(PHAST, c("Chrom", "Position", "Score"))

#Load Synonymity
CLASS.INFO<-fread("DATABASES/CANCER_DATA/1000GENOME/2013/ANOVAR.RESULTS/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.avinput.exonic_variant_function",
                  header=F, sep="\t", stringsAsFactors=F, drop=c(1, 6:11))
setnames(CLASS.INFO, c("CLASS", "INFO","Chrom","Position" ))
CLASS.INFO<-CLASS.INFO[CLASS %in% c("nonsynonymous SNV", "synonymous SNV") & INFO!="UNKNOWN",]
CLASS.INFO$GENE<-sapply(CLASS.INFO$INFO, function(x) strsplit(x, ":")[[1]][1])
CLASS.INFO$INFO<-NULL

#Correlation of MAF and SCORE across all genes
PHAST.CLASS<-merge(CLASS.INFO, PHAST, by=c("Chrom", "Position"))
PHAST.CLASS<-merge(PHAST.CLASS, THOUSAND[,c(1,2,3,5),with=F],  by=c("Chrom", "Position"))
ggplot(PHAST.CLASS, aes(Score, MAF, colour=CLASS)) + geom_point() + theme.format + facet_wrap(~CLASS) + 
  theme(strip.text.x = element_text(size = 20)) + scale_y_log10()

#Test a few genes
Function.Gene.PhastFit<-function(EXONS.TABLE, THOUSAND.TABLE, PHASTCON.TABLE, CLASS.TABLE , GENE){
  EXONS.GENE<-EXONS.TABLE[Hugo_Symbol==GENE,]
  print (EXONS.GENE)
  THOUSAND.GENE<-THOUSAND.TABLE[Chrom==unique(EXONS.GENE$Chrom) & Position>=unique(EXONS.GENE$START) & Position<=unique(EXONS.GENE$END),]
  print (THOUSAND.GENE)
  PHAST.GENE<-PHASTCON.TABLE[Chrom==unique(EXONS.GENE$Chrom) & Position>=unique(EXONS.GENE$START) & Position<=unique(EXONS.GENE$END),]
  print (PHAST.GENE)
  GENE<-merge(THOUSAND.GENE[,c(1,2,3,5),with=F], PHAST.GENE, by=c("Chrom", "Position"))
  
  GENE$REGION<-ifelse(apply(apply(unique(EXONS.GENE[,c("FEAT_START", "FEAT_END"),with=F]), 1, function(x) x[1]<=GENE$Position & GENE$Position<=x[2]),1,sum)>0, "EXON","INTRON")
  
  #Apply class
  CLASS.GENE<-CLASS.TABLE[GENE==GENE,]
  GENE<-merge(GENE, CLASS.GENE[,1:3,with=F], by=c("Chrom","Position"), all.x=T)
  
  return(GENE)
}

TTN<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST, CLASS.INFO,"TTN")
ggplot(TTN[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="TTN - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(TTN[TYPE=="SNP" & Score<=1.0 & MAF>=0.001,], aes(factor(Score), MAF, colour=REGION)) + geom_boxplot() + theme.format + facet_wrap(~REGION)  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(TTN[TYPE=="SNP",], aes(Score, colour=REGION)) + geom_histogram() + theme.format + facet_wrap(~REGION) +
  ggtitle(label="TTN - PHASTCONS HISTOGRAM (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(TTN[TYPE=="SNP" & CLASS!="NA",], aes(Score, MAF, colour=CLASS)) + geom_point() +theme.format  + facet_wrap(~CLASS) + scale_y_log10() +
  ggtitle(label="TTN - MAF vs PHASTCONS (SYN/NON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

TP53<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST, CLASS.INFO,"TP53")
table(TP53$REGION)
TP53[CLASS!="NA",]
ggplot(TP53[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="TP53 - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(TP53[TYPE=="SNP" & CLASS!="NA",], aes(Score,MAF, colour=CLASS)) + geom_point() +theme.format + facet_wrap(~CLASS) +scale_y_log10() +
  ggtitle(label="TP53 - MAF vs PHASTCONS (SYN/NON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

PIK3CA<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST, CLASS.INFO,"PIK3CA")
table(PIK3CA$REGION)
ggplot(PIK3CA[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="PIK3CA - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(PIK3CA[TYPE=="SNP" & CLASS!="NA",], aes(Score, MAF, colour=CLASS)) + geom_point() +theme.format + facet_wrap(~CLASS) + scale_y_log10() +
  ggtitle(label="PIK3CA - MAF vs PHASTCONS (SYN/NON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

BRCA1<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST, CLASS.INFO, "BRCA1")
table(BRCA1$REGION)
ggplot(BRCA1[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="BRCA1 - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(BRCA1[TYPE=="SNP" & CLASS!="NA",], aes(Score, MAF, colour=CLASS)) + geom_point() +theme.format + facet_wrap(~CLASS) + scale_y_log10() +
  ggtitle(label="BRCA1 - MAF vs PHASTCONS (SYN/NON SYN)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

BRCA2<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST, CLASS.INFO,"BRCA2")
table(BRCA2$REGION)
ggplot(BRCA2[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="BRCA2 - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(BRCA2[TYPE=="SNP" & CLASS!="NA",], aes(Score, MAF, colour=CLASS)) + geom_point() +theme.format + facet_wrap(~CLASS) + scale_y_log10() +
  ggtitle(label="BRCA2 - MAF vs PHASTCONS (SYN/NON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

AKT1<-Function.Gene.PhastFit(EXONS, THOUSAND, PHAST,CLASS.INFO, "AKT1")
table(AKT1$REGION)
ggplot(AKT1[TYPE=="SNP",], aes(Score, MAF, colour=REGION)) + geom_point() +theme.format  + facet_wrap(~REGION)  + scale_y_log10() +
  ggtitle(label="AKT1 - MAF vs PHASTCONS (EXON/INTRON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))
ggplot(AKT1[TYPE=="SNP" & CLASS!="NA",], aes(Score, MAF, colour=CLASS)) + geom_point() +theme.format + facet_wrap(~CLASS) + scale_y_log10() +
  ggtitle(label="AKT11 - MAF vs PHASTCONS (SYN/NON)" )+ 
  theme(strip.text.x = element_text(size = 18), plot.title = element_text(size=24, face="bold"))

#####1000GENOME READ DEPTH
THOUSAND.AVINPUT<-fread("DATABASES/CANCER_DATA/1000GENOME/2013/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.avinput", header=F, sep="\t",
                     stringsAsFactors=F, drop=c(3:6))
setnames(THOUSAND.AVINPUT, c("Chrom", "Position", "QUAL","DP"))
THOUSAND.AVINPUT<-THOUSAND.AVINPUT[DP!=".",] #Remove those whole read depth is undefined (removed 83193 records)
THOUSAND.AVINPUT$DP<-as.numeric(THOUSAND.AVINPUT$DP)

ggplot(THOUSAND.AVINPUT, aes(QUAL)) + geom_histogram()
THOUSAND.AVINPUT[order(QUAL, decreasing=T),]
THOUSAND.AVINPUT[QUAL>2,] # Filter at probabilty that alt call is wrong at 1%

head(THOUSAND.AVINPUT[order(DP, decreasing=T),],1000)
THOUSAND.AVINPUT[DP>25000,]
ggplot(THOUSAND.AVINPUT, aes(DP)) + geom_histogram() + scale_x_log10()

#Filter data by DP and 
FILTERED.PHAST.20000<-merge(PHAST.CLASS, THOUSAND.AVINPUT[QUAL>2 & DP>20000,], by=c("Chrom", "Position"))
FILTERED.PHAST.25000<-merge(PHAST.CLASS, THOUSAND.AVINPUT[QUAL>2 & DP>25000,], by=c("Chrom", "Position"))
FILTERED.PHAST.30000<-merge(PHAST.CLASS, THOUSAND.AVINPUT[QUAL>2 & DP>30000,], by=c("Chrom", "Position"))
ggplot(FILTERED.PHAST.20000, aes(Score)) + geom_histogram()
ggplot(FILTERED.PHAST.25000, aes(factor(Score), MAF, colour=CLASS)) + geom_boxplot() + facet_wrap(~CLASS) + scale_y_log10()
