#METABO ANALYSIS #1
#HISTOGRAM OF GENE CARDINALITY AND SCORE PER PATIENT
#101613 - FROM 150.py, 151.py

TABLE=read.table("NETWORK/R_ANALYSIS/101613_PATIENT_GENE_COVERAGE", header=TRUE, dec=".", check.names=FALSE)
DATA=as.data.frame(TABLE)
CANCERS=as.vector(unique(DATA$CANCER)); CANCERS
head(DATA)

DATA$PERCENT.COVERAGE<-DATA$METABOLIC_MUTATED_GENES/DATA$ALL_MUTATED_GENES

PATIENT.COVERAGE<-function(df, cancer) {
  library(ggplot2)
  
  #Plot coverage first
  TYPE<- df[df$CANCER==cancer,]
  TYPE$PATIENT_COUNT<-1:length(TYPE$PATIENT)
  
  CANCER.COVERAGE<-ggplot(TYPE, aes(x=PATIENT_COUNT)) + geom_line(aes(y=ALL_MUTATED_GENES, colour="ALL_MUTATED_GENES")) +
    geom_line(aes(y=METABOLIC_MUTATED_GENES, colour="METABOLIC_MUTATED_GENES")) + labs(title=cancer)
  
  #Difference between areas as percentage
  PERCENT.DIFF=sum(TYPE$METABOLIC_MUTATED_GENES)/sum(TYPE$ALL_MUTATED_GENES)
  
  #Distribution of non-metabolic mutated genes vs metabolic mutated genes across patients - NORMALIZED TO ALL MUTATED GENES
  NMMG<-data.frame(MUT=(TYPE$NON_METABOLIC_MUTATED_GENES)/TYPE$ALL_MUTATED_GENES) #Non-metabolic mutated genes will be different between all and metabolic
  MMG<-data.frame(MUT=(TYPE$METABOLIC_MUTATED_GENES)/TYPE$ALL_MUTATED_GENES)
  NMMG$OF<-"NMMG"
  MMG$OF<-"MMG"
  COMBINED.2<-rbind(NMMG, MMG)
  
  CANCER.HISTOGRAM<-ggplot(COMBINED.2, aes(x=MUT, fill=OF)) + geom_histogram(position="identity", alpha=0.3,binwidth=0.01) + 
    labs(title=cancer)
  
  #Store and return
  RESULT=list("CANCER.COVERAGE"=CANCER.COVERAGE, "CANCER.HISTOGRAM"=CANCER.HISTOGRAM, "PERCENT.DIFF"=PERCENT.DIFF)
  return(RESULT)  
}

A<-PATIENT.COVERAGE(DATA, "UCEC")
A$CANCER.COVERAGE
A$CANCER.HISTOGRAM
A$PERCENT.DIFF

#BOXPLOT of Distribution of mutations per patient of metabolic vs non-metabolic - THEN PAIRED T-TEST
#Need to reformat into another data frame
DATA$MMG.COVERAGE<-DATA$METABOLIC_MUTATED_GENES/DATA$ALL_MUTATED_GENES
DATA$NMMG.COVERAGE<-DATA$NON_METABOLIC_MUTATED_GENES/DATA$ALL_MUTATED_GENES
head(DATA)
BOXPLOT2A<-data.frame(CANCER=DATA$CANCER, MUTATED.GENES.PROPORTION=DATA$MMG.COVERAGE, TYPE="MMG")
BOXPLOT2B<-data.frame(CANCER=DATA$CANCER, MUTATED.GENES.PROPORTION=DATA$NMMG.COVERAGE, TYPE="NMMG")
BOXPLOT2<-rbind(BOXPLOT2A, BOXPLOT2B)
head(BOXPLOT2)
ggplot(BOXPLOT2, aes(factor(CANCER), MUTATED.GENES.PROPORTION)) + geom_boxplot(aes(fill=factor(TYPE)))

#PAIRED-T-TEST
PAIRED.T.TEST_MUTATIONS<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in CANCERS) {
    T_TEST<-t.test(DATA[DATA$CANCER==cancer,]$MMG.COVERAGE, DATA[DATA$CANCER==cancer,]$NMMG.COVERAGE, paired=TRUE)
    PAIRED.T.TEST_MUTATIONS<-rbind(PAIRED.T.TEST_MUTATIONS, data.frame(CANCER=cancer, P.VALUE=T_TEST$p.value))
}
PAIRED.T.TEST_MUTATIONS

#TODO:  
  #do gene ontology to filter out which are only metabolites(could do by KEGG as well)

#Add column to DATA of patient count proportion, may need to do loop
DATA.ORDERED<-DATA[order(DATA$CANCER, -DATA$ALL_MUTATED_GENES),]
head(DATA.ORDERED)
COUNT.VECTOR<-numeric(0)
for (cancer in as.vector(unique(DATA.ORDERED$CANCER))) {
  DATA.ORDERED.CANCER<-DATA.ORDERED[DATA.ORDERED$CANCER==cancer,]
  COUNT.VECTOR<-c(COUNT.VECTOR, c(1:nrow(DATA.ORDERED.CANCER))/nrow(DATA.ORDERED.CANCER))
}
DATA.ORDERED$PATIENT.COUNT<-COUNT.VECTOR
head(DATA.ORDERED)

#Plot of mutated_genes versus patients - MAKE SURE TO ADD COEFFICIENT OF VARIATION TO EACH PLOT (SD/M)
#Get CV for each cancer
CANCERS.ALL_MUTATED_GENES.CV<-data.frame(CANCERS=c(), CV=c(), MEAN=c())
for (cancer in CANCERS) {
  CANCER.TYPE<-DATA.ORDERED[DATA.ORDERED$CANCER==cancer,]
  CANCER.CV<-data.frame(CANCERS=cancer, CV=sd(CANCER.TYPE$ALL_MUTATED_GENES)/mean(CANCER.TYPE$ALL_MUTATED_GENES), 
                        MEAN=mean(CANCER.TYPE$ALL_MUTATED_GENES))
  CANCERS.ALL_MUTATED_GENES.CV<-rbind(CANCERS.ALL_MUTATED_GENES.CV, CANCER.CV)
}
CANCERS.ALL_MUTATED_GENES.CV
head(DATA.ORDERED)
ggplot(data=DATA.ORDERED, aes(x=PATIENT.COUNT, y=ALL_MUTATED_GENES, col=NUMBER_OF_UNIQUE_METABOLITES)) + geom_point() +
  geom_line(aes(y=118), colour="red") + scale_y_log10() +stat_smooth() + facet_wrap(~CANCER)

#BOXPLOT of average mutations per gene non-metabolic vs metabolic genes #GET T-TEST(MANN.WHITNEY)!!!
#Need to extract into another dataframe
head(DATA)
DATA.BOXPLOT1A<-data.frame(CANCER=DATA$CANCER, AMPG=DATA$AVG_MUTATIONS_PER_GENE_ALL,TYPE="ALL")
DATA.BOXPLOT1B<-data.frame(CANCER=DATA$CANCER, AMPG=DATA$AVG_MUTATIONS_PER_GENE_M,TYPE="M")
DATA.BOXPLOT1C<-data.frame(CANCER=DATA$CANCER, AMPG=DATA$AVG_MUTATIONS_PER_GENE_NM,TYPE="NM")
DATA.BOXPLOT1<-rbind(DATA.BOXPLOT1B, DATA.BOXPLOT1C)
head(DATA.BOXPLOT1)
ggplot(DATA.BOXPLOT1A, aes(factor(CANCER), AMPG)) +geom_boxplot() +scale_y_log10()
ggplot(DATA.BOXPLOT1, aes(factor(CANCER), AMPG)) + geom_boxplot(aes(fill=factor(TYPE))) + scale_y_log10()

MANN.WHITNEY.TEST_AVG_MUTATIONS_PER_GENE<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in CANCERS) {
  #Remove ties first
  WILCOX.DATA<-DATA[DATA$CANCER==cancer,]
  WILCOX.DATA<-WILCOX.DATA[(WILCOX.DATA$AVG_MUTATIONS_PER_GENE_M-WILCOX.DATA$AVG_MUTATIONS_PER_GENE_NM)!=0,]
  #Do test
  M.WHITNEY<-wilcox.test(WILCOX.DATA$AVG_MUTATIONS_PER_GENE_M, WILCOX.DATA$AVG_MUTATIONS_PER_GENE_NM)
  MANN.WHITNEY.TEST_AVG_MUTATIONS_PER_GENE<-rbind(MANN.WHITNEY.TEST_AVG_MUTATIONS_PER_GENE, data.frame(CANCER=cancer, P.VALUE=M.WHITNEY$p.value))
}
MANN.WHITNEY.TEST_AVG_MUTATIONS_PER_GENE

#Plots of percent coverage versus patients
#Need to work on another sorting and loop to get it this time by COVERAGE
DATA.ORDERED.COVERAGE<-DATA[order(DATA$CANCER, -DATA$PERCENT.COVERAGE),]
COUNT.VECTOR.COVERAGE<-numeric(0)
for (cancer in as.vector(unique(DATA.ORDERED.COVERAGE$CANCER))) {
  DATA.ORDERED.COVERAGE.CANCER<-DATA.ORDERED.COVERAGE[DATA.ORDERED.COVERAGE$CANCER==cancer,]
  COUNT.VECTOR.COVERAGE<-c(COUNT.VECTOR.COVERAGE, c(1:nrow(DATA.ORDERED.COVERAGE.CANCER))/nrow(DATA.ORDERED.COVERAGE.CANCER))
}
DATA.ORDERED.COVERAGE$PATIENT.COUNT<-COUNT.VECTOR.COVERAGE
head(DATA.ORDERED.COVERAGE)
ggplot(DATA.ORDERED.COVERAGE, aes(x=PATIENT.COUNT, y=PERCENT.COVERAGE, col=CANCER)) + geom_point() #+ stat_smooth()

#PLOT Metabolic mutated genes vs non-metabolic mutated genes -#GET CORRELATION AND SLOPE!!! - DONE
head(DATA)
ggplot(DATA, aes(x=METABOLIC_MUTATED_GENES, y=NON_METABOLIC_MUTATED_GENES, col=AVG_MUTATIONS_PER_GENE_ALL)) +  
  geom_point() + scale_y_log10() + scale_x_log10() + facet_wrap(~CANCER)

#Correlation betwen metabolic mutated genes and non-metabolic mutated genes per cancer
PEARSON.MUTATED_GENES<-data.frame(CANCER=c(), P.VALUE=c(), RHO=c())
for (cancer in CANCERS) {
  CANCER.TYPE<-DATA[DATA$CANCER==cancer,]
  CORRELATION<-cor.test(CANCER.TYPE$METABOLIC_MUTATED_GENES, CANCER.TYPE$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  CORRELATION.DF<-data.frame(CANCER=cancer,P.VALUE=format.pval(CORRELATION$p.value),RHO=CORRELATION$estimate)
  PEARSON.MUTATED_GENES<-rbind(PEARSON.MUTATED_GENES, CORRELATION.DF)  
}
PEARSON.MUTATED_GENES

#PLOT a bar chart of percentage of common mutated genes across cancers
COMMON.IN.CANCERS<-data.frame(CANCERS=rep(c("LUSC", "READ", "GBM", "KIRC", "UCEC", "OV","BRCA", "COAD"), each=2),
                              FOUND=rep(c("COMMON", "UNCOMMON"), 8),
                              COUNT=c(370, 12275-370, 370, 4082-370, 370, 7028-370, 370, 7362-370, 
                                      370, 16342-370, 370, 3238-370, 370, 11179-370, 370, 4616-370))
COMMON.IN.CANCERS
ggplot(COMMON.IN.CANCERS, aes(CANCERS, COUNT, fill=FOUND)) + geom_bar()

#10/29/13
#DO EXACT BINOMIAL TEST (two-tailed) to compare the proportions of total number of genes in human genome, total metabolic genes, mutated genes per cancer, metabolic mutated genes, get bonferroni correction of p-values
CANCER.GENES.COUNT<-data.frame(CANCER=c("LUSC", "READ", "GBM", "KIRC", "UCEC", "OV","BRCA", "COAD"),
                               ALL=c(12275, 4082, 7028, 7362, 16342, 3238, 11179, 4616),
                               METABOLIC=c(3482,1276,2045,2228, 
                                       4478,967,3193,1521))
CANCER.GENES.COUNT$NON.METABOLIC<-CANCER.GENES.COUNT$ALL-CANCER.GENES.COUNT$METABOLIC
library(reshape)
MELTED.CANCER.GENES.COUNT<-melt(CANCER.GENES.COUNT);MELTED.CANCER.GENES.COUNT
colnames(MELTED.CANCER.GENES.COUNT)<-c("CANCER", "TYPE", "GENES")

MELTED.CANCER.GENES.COUNT                               
ggplot(MELTED.CANCER.GENES.COUNT[MELTED.CANCER.GENES.COUNT$TYPE %in% c("METABOLIC", "NON.METABOLIC"),], aes(CANCER, GENES, fill=TYPE)) + geom_bar()
 
CANCER.GENES.COUNT
CANCER.GENES.COUNT.PVALUE<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in CANCERS) {
  ROW=CANCER.GENES.COUNT[CANCER.GENES.COUNT$CANCER==cancer,]
  ROW.TEST=binom.test(c(ROW$METABOLIC, ROW$NON.METABOLIC), p=5113/21000) #THE ORDER IS IMPORTANT, THE RATIO IS FOR METABOLIC TO ALL
  CANCER.GENES.COUNT.PVALUE<-rbind(CANCER.GENES.COUNT.PVALUE,data.frame(CANCER=cancer, P.VALUE=ROW.TEST$p.value))
}
CANCER.GENES.COUNT.PVALUE$P.CORRECTED<-p.adjust(CANCER.GENES.COUNT.PVALUE$P.VALUE, method="bonferroni")
CANCER.GENES.COUNT.PVALUE

#10/30/13
#PLOT DISTRIBUTION OF GENE LENGTH OF METABOLIC VS OTHER
library(ggplot2)
GENE_LENGTH<-as.data.frame(read.table("NETWORK/R_ANALYSIS/103013_CANCER_MET_NON_AA_LENGTH", header=TRUE, sep="\t"))
head(GENE_LENGTH)
ggplot(GENE_LENGTH,aes(factor(CANCER), LENGTH)) +geom_boxplot(aes(fill=(TYPE))) + scale_y_log10()
#Do Mann whitney
LENGTHS.WHITNEY.TEST<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in as.vector(unique(GENE_LENGTH$CANCER))) {
  MANN.DATA<-GENE_LENGTH[GENE_LENGTH$CANCER==cancer,]
  METABOLIC<-MANN.DATA[MANN.DATA$TYPE=="METABOLIC",]
  NON.METABOLIC<-MANN.DATA[MANN.DATA$TYPE=="NON_METABOLIC",]
  LENGTH.WHITNEY<-wilcox.test(METABOLIC$LENGTH, NON.METABOLIC$LENGTH)
  LENGTHS.WHITNEY.TEST<-rbind(LENGTHS.WHITNEY.TEST, data.frame(CANCER=cancer, P.VALUE=LENGTH.WHITNEY$p.value))
}
LENGTHS.WHITNEY.TEST$P.CORRECTED<-p.adjust(LENGTHS.WHITNEY.TEST$P.VALUE, method="bonferroni")
LENGTHS.WHITNEY.TEST

#BIN GENES IN EACH CANCER BY SIZE, COMPARED METBOLIC BIN VS NON-METABOLIC BIN OF SAME SIZE
head(GENE_LENGTH)
ggplot(GENE_LENGTH, aes(x=LENGTH, fill=TYPE)) +geom_histogram() + facet_wrap(~CANCER) + scale_x_log10()

#Get bins into DECILES, get average across cancers
CANCERS<-as.vector(unique(GENE_LENGTH$CANCER));CANCERS
BINNING<-as.data.frame(sapply(CANCERS, function(x,y) quantile(y[[x]]$LENGTH, prob=seq(0,1,length=11), type=5) , y=split(GENE_LENGTH, GENE_LENGTH$CANCER)))
BINNING
BINNING$PROTEOME<-NULL#Only interested on binning for cancers, not entire affected proteome
BINS<-apply(BINNING, 1, mean) #Getting the averages per decile across canceers
BINS<-as.data.frame(BINS)$BINS; BINS
BINS[1]<-0
BINS
GENE_LENGTH$LENGTH.BINS<-cut(GENE_LENGTH$LENGTH, BINS); head(GENE_LENGTH)
ggplot(GENE_LENGTH,aes(factor(CANCER), LENGTH)) +geom_boxplot(aes(fill=(TYPE))) +facet_grid(LENGTH.BINS~., scales="free")   + scale_y_log10()  #+ coord_flip() #REMEMBER TO FIND SIGNIFICANCE IN "TRUE" SHIFTS WITH SMALLER BINS!!!!

#11/31/13
#ACTUALLY, PROJECT IS RELATED TO COUNTS!!, NOT LENGTHS PROPORTIONS!!
head(GENE_LENGTH)
CANCER_GENE_COUNT_BY_LENGTH<-data.frame(CANCER=c(), TYPE=c(), LENGTH.BINS=c(), COUNT=c())
GENE_LENGTH_SPLIT=split(GENE_LENGTH, GENE_LENGTH$CANCER)

for (cancer in as.vector(unique(GENE_LENGTH$CANCER))) {  
  CANCER_BINS<-split(GENE_LENGTH_SPLIT[[cancer]], GENE_LENGTH_SPLIT[[cancer]]$LENGTH.BINS)

  for (bin in names(CANCER_BINS)) {
    BIN_METABO<-data.frame(CANCER=cancer, TYPE="METABOLIC", 
                           LENGTH.BINS=bin, COUNT=nrow(CANCER_BINS[[bin]][CANCER_BINS[[bin]]$TYPE=="METABOLIC",]))
    BIN_NON_METABO<-data.frame(CANCER=cancer, TYPE="NON_METABOLIC", 
                           LENGTH.BINS=bin, COUNT=nrow(CANCER_BINS[[bin]][CANCER_BINS[[bin]]$TYPE=="NON_METABOLIC",]))
    CANCER_GENE_COUNT_BY_LENGTH<-rbind(CANCER_GENE_COUNT_BY_LENGTH, BIN_METABO)
    CANCER_GENE_COUNT_BY_LENGTH<-rbind(CANCER_GENE_COUNT_BY_LENGTH, BIN_NON_METABO)
  }
}

head(CANCER_GENE_COUNT_BY_LENGTH,15)
ggplot(CANCER_GENE_COUNT_BY_LENGTH, aes(CANCER, COUNT, fill=TYPE)) +geom_bar() +facet_wrap(~LENGTH.BINS)

#Split and join by type (metabolic, non-metbolic)
LENGTH.COUNT.SPLIT<-split(CANCER_GENE_COUNT_BY_LENGTH, CANCER_GENE_COUNT_BY_LENGTH$TYPE)
LENGTH.COUNT.SPLIT$METABOLIC$TYPE<-NULL
colnames(LENGTH.COUNT.SPLIT$METABOLIC)<-c("CANCER", "LENGTH.BINS", "METABOLIC")
head(LENGTH.COUNT.SPLIT$METABOLIC)

LENGTH.COUNT.SPLIT<-cbind(LENGTH.COUNT.SPLIT$METABOLIC, NON_METABOLIC=LENGTH.COUNT.SPLIT$NON_METABOLIC$COUNT) #TOTAL NUMBER OF MUTATED PROTEINS PER BIN 
                                                                                                              #PER CANCER FOR METABOLIC AND NON-METABOLIC

#To obtaine number across all lengths parse BOXPLOT2 (Seen as the ratio seen in BOXPLOT2 plotted)
head(BOXPLOT2)
BOXPLOT2.SPLIT<-split(BOXPLOT2, BOXPLOT2$TYPE)

BOXPLOT2.METABOLIC<-sapply(as.vector(unique(BOXPLOT2$CANCER)), function(x,y)  mean(y[[x]]$MUTATED.GENES.PROPORTION) , #AVERAGES ACROSS PATIENTS
       y=split(BOXPLOT2.SPLIT$MMG, BOXPLOT2.SPLIT$MMG$CANCER)) #add for "PROTEOME" at the end
BOXPLOT2.METABOLIC<-as.data.frame(t(BOXPLOT2.METABOLIC))

BOXPLOT2.NON_METABOLIC<-sapply(as.vector(unique(BOXPLOT2$CANCER)), function(x,y)  mean(y[[x]]$MUTATED.GENES.PROPORTION) , #AVERAGES ACROSS PATIENTS
                           y=split(BOXPLOT2.SPLIT$NMMG, BOXPLOT2.SPLIT$NMMG$CANCER)) 
BOXPLOT2.NON_METABOLIC<-as.data.frame(t(BOXPLOT2.NON_METABOLIC))

BOXPLOT2.NON_METABOLIC
BOXPLOT2.METABOLIC
MUTATED.GENE.RATIO<-BOXPLOT2.NON_METABOLIC/BOXPLOT2.METABOLIC #AVERAGE ACROSS PATIENTS (THE MEAN SEEN IN BOXPLOT ABOVE)
MUTATED.GENE.RATIO$PROTEOME=15143/4964 #Added for "PROTEOME"
MUTATED.GENE.RATIO

#Add ratios to each cancer in LENGTH.COUNT.SPLIT(counts of mutated genes across cancers)
head(LENGTH.COUNT.SPLIT)
LENGTH.COUNT.SPLIT$RATIO<-sapply(as.vector(LENGTH.COUNT.SPLIT$CANCER) , function(x) MUTATED.GENE.RATIO[[x]])

#PLOT ratios versus bins per cancer - AS WE GET CLOSER TO COMPARING LARGER BINS TO EACH OTHER WE CONVERGE AT RATIO SEEN BY PATIENTS
ggplot(data=LENGTH.COUNT.SPLIT, aes(x=LENGTH.BINS, y=NON_METABOLIC/METABOLIC, col=LENGTH.BINS)) + geom_point() + 
  geom_hline(aes(yintercept=RATIO), col="red") + geom_hline(yintercept=15143/4964, col="black", linetype=5)+
  theme(axis.text.x=element_text(angle=45, vjust=0.5))+ facet_wrap(~CANCER)

#11/3/13
#HEATMAP OF PATIENT GENE LENGTH VECTOR FOR METABOLIC AND NON-METABOLIC MUTATED GENES
HEATMAP1.DATA<-DATA
head(HEATMAP1.DATA,3)
HEATMAP1.DATA<-HEATMAP1.DATA[,-6:-9]

#Bin number of genes should be a proportion of ALL_MUTATED_GENES
HEATMAP1.DATA[,6:25]<-HEATMAP1.DATA[,6:25]/HEATMAP1.DATA$ALL_MUTATED_GENES

HEATMAP1.DATA$ALL_MUTATED_GENES<-NULL
HEATMAP1.DATA<-HEATMAP1.DATA[,-ncol(HEATMAP1.DATA)]

#Get P.VALUE per patient for HYPERGEOMETRIC test against expected ratio of 5113/21000
HEATMAP1.PVALUES<-sapply(1:nrow(HEATMAP1.DATA), function(x,y) 
  phyper(y$METABOLIC_MUTATED_GENES[x], 5113, 21000-5113, y$NON_METABOLIC_MUTATED_GENES[x], lower.tail=FALSE),
                         y=HEATMAP1.DATA)
HEATMAP1.DATA$PVALUES<-HEATMAP1.PVALUES

#CONCLUSIONS FOR THE DAY:
#IT APPEARS AT FIRST THAT THE MOST PATIENTS HAVE A LARGE NUMBER OF LONG LENGTH MUTATE GENES IN BOTH CASES FOR METABOLIC AND NON-METABOLIC GENES, THIS WOULD IMPLY THAT THE EFFECT WE SEE IS SOMETHING ELSE AS THERE IS NO PATTERN AS GOING FROM LOW P VALUE TO HIGH P VALUE, THOUGH WE CAN DO SOME THINGS:
#Since longer genes may seem to make no diference, we could get rid of the (1496.05001,35991.0] and see if now we see a pattern
#We could also cluster them first, like rows or columns (for columns to see effect of comparing within similar lengths), and then do       
#the histogram
#CLUSTER BY RATIO, do hclust like in class and pick 2,3 clusters


#11/03/13
#Separate by CANCERS, FDR correct, then order by P.values, convert to matrix to assign p.values as rownames, delete p.value column
head(HEATMAP1.DATA)

###Changed col names to visually appealing bins
colnames(HEATMAP1.DATA)<-c(colnames(HEATMAP1.DATA[1:4]),sapply(strsplit(colnames(HEATMAP1.DATA)[5:24], ","), function(l) paste(l[1], "-->", l[2])), colnames(HEATMAP1.DATA[25]))

HEATMAP1.DATA.SPLIT<-split(HEATMAP1.DATA, HEATMAP1.DATA$CANCER)
summary(HEATMAP1.DATA.SPLIT)

#Bonferroni P.value within cancers, get ratio
for (cancer in as.vector(unique(HEATMAP1.DATA$CANCER))) {
  HEATMAP1.DATA.SPLIT[[cancer]]$CORRECTED.PVALUES<-p.adjust(HEATMAP1.DATA.SPLIT[[cancer]]$PVALUES, method="fdr")
  HEATMAP1.DATA.SPLIT[[cancer]]<-HEATMAP1.DATA.SPLIT[[cancer]][order(HEATMAP1.DATA.SPLIT[[cancer]]$CORRECTED.PVALUES),]
  HEATMAP1.DATA.SPLIT[[cancer]]$RATIO<-
    HEATMAP1.DATA.SPLIT[[cancer]]$NON_METABOLIC_MUTATED_GENES/HEATMAP1.DATA.SPLIT[[cancer]]$METABOLIC
}
head(HEATMAP1.DATA.SPLIT$OV)

#Use FDR corrected values at >0.05 to filter significant patients
HEATMAP1.DATA.SPLIT.LISTS<-list()
  
#Separate by  CORRECTED.P.VALUE<0.05 and CORRECTED.P.VALUE>0.05 per cancer, remove unnecessary lines and set matrix rows as PVALUES
for (cancer in as.vector(unique(HEATMAP1.DATA$CANCER))) {
  
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG<-split(HEATMAP1.DATA.SPLIT[[cancer]], 
                                                              HEATMAP1.DATA.SPLIT[[cancer]]$CORRECTED.PVALUES<0.05)$`TRUE`
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG<-split(HEATMAP1.DATA.SPLIT[[cancer]], 
                                                                   HEATMAP1.DATA.SPLIT[[cancer]]$CORRECTED.PVALUES<0.05)$`FALSE`
  
  #Remove unnecessary columns(CANCER, PATIENT, RATIO) and set as matrix with rows as CORRECTED.PVALUES
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG[,c(5:14,26)])
  rownames(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M)<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M[,ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M)]
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M[,-ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M)]
  
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG[,c(15:24,26)])
  rownames(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM[,ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)]
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM[,-ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)]
  
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG[,c(5:14,26)])
  rownames(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M)<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M[,ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M)]
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M[,-ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M)]  
  
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG[,c(15:24,26)])
  rownames(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM[,ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)]
  HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM<-
    HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM[,-ncol(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)]
}
#11/4/13
#TO DO:
# Try to use smaller number of bins!!!! (i.e 5)
# Find pattern between SIG vs NON_SIG
# Find a way to plot the distribution of counts of bins nm vs m
library(reshape)
library(ggplot2)

#PLOT - First do plot of proportion of significant patients versus non-significant patients per cancer
head(HEATMAP1.DATA.SPLIT.LISTS$COAD$SIG)
SPLIT.SIG<-as.data.frame(as.matrix(sapply(CANCERS, function(x,y) nrow(y[[x]]$SIG), y=HEATMAP1.DATA.SPLIT.LISTS)))
SPLIT.NON.SIG<-as.data.frame(as.matrix(sapply(CANCERS, function(x,y) nrow(y[[x]]$NON_SIG), y=HEATMAP1.DATA.SPLIT.LISTS)))
SPLIT.BOX<-cbind(SPLIT.SIG,SPLIT.NON.SIG)
colnames(SPLIT.BOX)<-c("SIG", "NON_SIG")
SPLIT.BOX$CANCERS<-rownames(SPLIT.BOX)
SPLIT.BOX$TOTAL<-SPLIT.BOX$SIG + SPLIT.BOX$NON_SIG

SPLIT.BOX.MELTED<-melt(SPLIT.BOX, id=c("CANCERS", "TOTAL"))

ggplot(SPLIT.BOX.MELTED, aes(CANCERS, value, fill=variable)) + geom_bar(position="dodge") +ylab("Number of Patients") + 
  geom_text(aes(label=paste(round(100*value/TOTAL,2), "%")), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(title="Significant Number of Patients")

#PLOT
#Get Heatmap per cancer on SIG and NON_SIG
library(RColorBrewer)
hmcol<-brewer.pal(9,"Blues")

#For significant
for (cancer in CANCERS) {
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/110413_", cancer, "_SIG_10.jpeg", sep=""), width=1080, height=1080, quality=100, res=100)
  heatmap(t(cbind(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M, HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)), scale='none', col=hmcol, margins=c(11,11))
  dev.off()
} 

#For Non-significant
for (cancer in CANCERS) {
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/110413_", cancer, "_NON_SIG_10.jpeg", sep=""), width=1080, height=1080, quality=100, res=100)
  heatmap(t(cbind(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M, HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)), scale='none', col=hmcol, margins=c(11,11))
  dev.off()
}

#PLOT
#Binning percent counts in above heatmap

#First for significant
for (cancer in CANCERS) {
  TEMP.PLOT<-ggplot(melt(cbind(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M, HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)),
         aes(factor(X2), value)) + geom_boxplot() +theme(axis.text.x=element_text(angle=90)) +coord_flip() + xlab("Bins")
  ggsave(TEMP.PLOT, file=paste("NETWORK/PICTURES/METABO_ANALYSIS/HEATMAP_BOXPLOTS/110413_", cancer, "_SIG_10_BOXPLOT.jpeg", sep=""), dpi=600)
}

#Then for non-significant
for (cancer in CANCERS) {
  TEMP.PLOT<-ggplot(melt(cbind(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M, HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)),
                    aes(factor(X2), value)) + geom_boxplot() +theme(axis.text.x=element_text(angle=90)) +coord_flip() + xlab("Bins")
  ggsave(TEMP.PLOT, file=paste("NETWORK/PICTURES/METABO_ANALYSIS/HEATMAP_BOXPLOTS/110413_", cancer, "_NON_SIG_10_BOXPLOT.jpeg", sep=""), dpi=600)
}

#PLOT
#Do Heatmap of ratio of same bin NM/M for signicant and non-signicant across cancers
#IMPORTANT - SOME RATIO ARE DIVISIONS BY ZERO, TO AVOID FILTERING THEM OUT, WE ARE GOING TO USE THE APPROXIMATION ln(X+1)~X -> ln(X+1.01)~X

#First for significant, remember to combine into LISTS (one per cancer) for later analysis of ratios by binning
BOXPLOT3.RATIO.HEATMAP<-list()
for (cancer in CANCERS) {
  DIVIDED.HEATMAP<-log(1.01+HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_NM)/log(1.01+HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG_M)
  
  #Get 0/0 Na equal to zero
  DIVIDED.HEATMAP[is.na(DIVIDED.HEATMAP)]<-0
  
  #Get Inf equal to highest ratio (second higest ratio taking Inf into account)
  #SECOND.MAX<-max(DIVIDED.HEATMAP[! DIVIDED.HEATMAP %in% max(DIVIDED.HEATMAP)])
  #DIVIDED.HEATMAP[is.infinite(DIVIDED.HEATMAP)]<-1.5*SECOND.MAX
  
  #INITIALIZE CANCER TO BOXPLOT3.RATIO.HEATMAP DATA FRAME  IN LIST FOR LATER COMPARISSON
  TEMP.FRAME<-melt(DIVIDED.HEATMAP)
  TEMP.FRAME$TYPE<-"SIG"
  BOXPLOT3.RATIO.HEATMAP[[cancer]]<-TEMP.FRAME
  
  #Save plot
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/RATIO/110413_", cancer, "_SIG_BIN_RATIO_10.jpeg", sep=""), width=1080, height=1080, quality=100, res=100)
  heatmap(t(DIVIDED.HEATMAP), scale='none',  margins=c(15,15) ,col=hmcol)
  dev.off()
}

#Then for non-significant
for (cancer in CANCERS) {
  DIVIDED.HEATMAP<-log(1.01+HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_NM)/log(1.01+HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG_M)
  
  #Get 0/0 Na equal to zero
  DIVIDED.HEATMAP[is.na(DIVIDED.HEATMAP)]<-0
  
  #Get Inf equal to highest ratio (second higest ratio taking Inf into account)
  #SECOND.MAX<-max(DIVIDED.HEATMAP[! DIVIDED.HEATMAP %in% max(DIVIDED.HEATMAP)])
  #DIVIDED.HEATMAP[is.infinite(DIVIDED.HEATMAP)]<-1.5*SECOND.MAX

  #ADD TO INITIALIZED BOXPLOT3.RATIO.HEATMAP$CANCER DATA FRAME  IN LIST FOR LATER COMPARISSON
  TEMP.FRAME<-melt(DIVIDED.HEATMAP)
  TEMP.FRAME$TYPE<-"NON_SIG"
  BOXPLOT3.RATIO.HEATMAP[[cancer]]<-rbind(BOXPLOT3.RATIO.HEATMAP[[cancer]], TEMP.FRAME)  
  
  #Save plot
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/RATIO/110413_", cancer, "_NON_SIG_BIN_RATIO_10.jpeg", sep=""), 
       width=1080, height=1080, quality=100, res=100)
  heatmap(t(DIVIDED.HEATMAP), scale='none', col=hmcol, margins=c(15,15))
  dev.off()
}

#PLOT
#DO BOXPLOT COMPARISSON OF RATIOS ACROSS BINS PER CANCER TO FIND SIGNIFICANT CHANGES
#Get significance line in there along with non-significance line

#To get significance line we need to add RATIO from MUTATED.GENE.RATIO described above
for (cancer in CANCERS) {
  BOXPLOT3.RATIO.HEATMAP[[cancer]]$RATIO<-MUTATED.GENE.RATIO[[cancer]]
}

tail(BOXPLOT3.RATIO.HEATMAP$LUSC)
for (cancer in CANCERS) { #"value" are the actual ratios from the distribution
  TEMP.PLOT<-ggplot(BOXPLOT3.RATIO.HEATMAP[[cancer]], aes(factor(X2), value)) + geom_boxplot(aes(fill=factor(TYPE))) + coord_flip() + xlab("Bins") +
    geom_hline(yintercept=15143/4964, col="black", linetype=5) + geom_hline(aes(yintercept=RATIO), col="red") + scale_y_log10()
  ggsave(TEMP.PLOT, file=paste("NETWORK/PICTURES/METABO_ANALYSIS/HEATMAP_BOXPLOTS/RATIO/110513_", cancer, "_BIN_RATIO_COMPARISSON.jpeg", sep=""),
         dpi=800, scale=2)
}
MUTATED.GENE.RATIO
TEMP.PLOT
#STAT TEST
#Do Mann-Whitney to calculate if there is a significance difference between the ratios in each per cancer
#Split by bin in each cancer and get mann whitney sig vs non_sig
head(BOXPLOT3.RATIO.HEATMAP$LUSC)
MANN.WHITNEY.RATIOS.SIG_NONSIG<-list()
for (cancer in CANCERS) {
  MANN.PVAL<-sapply(as.vector(unique(BOXPLOT3.RATIO.HEATMAP[[cancer]]$X2)), function(x,y) 
    wilcox.test(split(y[[x]], y[[x]]$TYPE)$SIG$value,  split(y[[x]], y[[x]]$TYPE)$NON_SIG$value)$p.value, 
             y=split(BOXPLOT3.RATIO.HEATMAP[[cancer]], BOXPLOT3.RATIO.HEATMAP[[cancer]]$X2))
  MANN.PVAL<-as.data.frame(MANN.PVAL)
  MANN.PVAL$CORRECTED.PVALUE<-p.adjust(MANN.PVAL$MANN.PVAL, method="bonferroni")
  MANN.WHITNEY.RATIOS.SIG_NONSIG[[cancer]]<-MANN.PVAL
  write(as.matrix(MANN.PVAL), file=paste("NETWORK/PICTURES/METABO_ANALYSIS/HEATMAP_BOXPLOTS/RATIO/110513_SIG_NONSIG_RATIO_", cancer, sep=""), sep="\t") 
}
View(MANN.WHITNEY.RATIOS.SIG_NONSIG$KIRC[MANN.WHITNEY.RATIOS.SIG_NONSIG$KIRC$CORRECTED.PVALUE<0.05, ])
MANN.WHITNEY.RATIOS.SIG_NONSIG$READ

#HEATMAP of significant ratios (SIG vs NON_SIG) across CANCERS
HEATMAP3.SIGNIFICANT.RATIOS<-data.frame(matrix(nrow=10, ncol=0))
for (cancer in CANCERS) {
  TEMP2<-MANN.WHITNEY.RATIOS.SIG_NONSIG[[cancer]]
  TEMP2[[cancer]]<-as.numeric(TEMP2$CORRECTED.PVALUE<0.05)
  TEMP2$MANN.PVAL<-NULL
  TEMP2$CORRECTED.PVALUE<-NULL
  colnames(TEMP2)<-cancer
  HEATMAP3.SIGNIFICANT.RATIOS<-cbind(HEATMAP3.SIGNIFICANT.RATIOS, TEMP2)
}
HEATMAP3.SIGNIFICANT.RATIOS
heatmap(as.matrix(HEATMAP3.SIGNIFICANT.RATIOS), scale="none", col=hmcol)

####APPLICATION#####
library(mgcv)
#1. Get total NM/M ratio (unibinned) of patients for significant vs non-significant (Uncoupling of plot on 128-129) - COLBIND OF HEATMAP1.DATA.SPLIT
CANCERS.ONLY<-CANCERS[1:8];CANCERS.ONLY
HEATMAP1.DATA.UNSPLIT<-data.frame()
for (cancer in CANCERS.ONLY) {
  HEATMAP1.DATA.UNSPLIT<-rbind(HEATMAP1.DATA.UNSPLIT, HEATMAP1.DATA.SPLIT[[cancer]])
}
tail(HEATMAP1.DATA.UNSPLIT)

#Add a column of TOTAL FDR just to check
HEATMAP1.DATA.UNSPLIT$EXTRAP<-p.adjust(HEATMAP1.DATA.UNSPLIT$PVALUES, method="fdr")
HEATMAP1.DATA.UNSPLIT$EXTRAPTEST<-HEATMAP1.DATA.UNSPLIT$EXTRAP<0.05

#PLOT UNCOUPLED RATIOS above
ggplot(HEATMAP1.DATA.UNSPLIT, aes(x=METABOLIC_MUTATED_GENES, y=NON_METABOLIC_MUTATED_GENES, col=EXTRAPTEST)) +
  geom_point() + scale_y_log10() + scale_x_log10() + facet_wrap(~CANCER)  

#PEARSON CORRELATION of split graph
HDUS<-split(HEATMAP1.DATA.UNSPLIT, HEATMAP1.DATA.UNSPLIT$CANCER)
head(HDUS$BRCA)

PEARSON.HDUS<-data.frame(CANCER=c(), SIG.P.VALUE=c(), SIG.RHO=c(), NON.SIG.P.VALUE=c(), NON.SIG.RHO=c())
for (cancer in CANCERS.ONLY) {
  CANCER.TYPE<-HDUS[[cancer]]
  SIG.CORRELATION<-cor.test(CANCER.TYPE[CANCER.TYPE$EXTRAPTEST==TRUE,]$METABOLIC_MUTATED_GENES, 
                            CANCER.TYPE[CANCER.TYPE$EXTRAPTEST==TRUE,]$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  NON.SIG.CORRELATION<-cor.test(CANCER.TYPE[CANCER.TYPE$EXTRAPTEST==FALSE,]$METABOLIC_MUTATED_GENES, 
                            CANCER.TYPE[CANCER.TYPE$EXTRAPTEST==FALSE,]$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  CORRELATION.DF<-data.frame(CANCER=cancer,SIG.P.VALUE=format.pval(SIG.CORRELATION$p.value),SIG.RHO=SIG.CORRELATION$estimate,
                             NON.SIG.P.VALUE=format.pval(NON.SIG.CORRELATION$p.value), NON.SIG.RHO=NON.SIG.CORRELATION$estimate)
  PEARSON.HDUS<-rbind(PEARSON.HDUS, CORRELATION.DF)  
}
View(PEARSON.HDUS)

#Get slope of split graph - HAD TO DO IT MANUALLY
head(HDUS$BRCA)
coef(lm(NON_METABOLIC_MUTATED_GENES~METABOLIC_MUTATED_GENES, HDUS$KIRC[HDUS$KIRC$EXTRAPTEST==FALSE,]))

#CLUSTERING
head(HEATMAP1.DATA.SPLIT.LISTS$BRCA$SIG)
head(TESTING.MATRIX)
for (cancer in CANCERS.ONLY) {
  #First for SIG
  TESTING.MATRIX<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG[,5:24])
  rownames(TESTING.MATRIX)<-HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$SIG$PATIENT
  
  CLUST<-kmeans(x=TESTING.MATRIX, 2, nstart=100)
  CLUST.MELTED<-melt(CLUST$centers)
  CLUST.MELTED$TYPE<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[1])
  CLUST.MELTED$BIN<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[2])
  
  SIG_PLOT<-ggplot(CLUST.MELTED, aes(x=BIN, y=value, col=TYPE)) +geom_line(aes(group=TYPE)) + facet_wrap(~X1) + 
    theme(axis.text.x=element_text(angle=90))
  ggsave(SIG_PLOT, file=paste("NETWORK/PICTURES/METABO_ANALYSIS/CLUSTERING/110513_", cancer, "_SIG_10_CLUSTER2.jpeg", sep=""), dpi=600)

  #Then for NON_SIG
  TESTING.MATRIX<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG[,5:24])
  rownames(TESTING.MATRIX)<-HEATMAP1.DATA.SPLIT.LISTS[[cancer]]$NON_SIG$PATIENT
  
  CLUST<-kmeans(x=TESTING.MATRIX, 2, nstart=100)
  CLUST.MELTED<-melt(CLUST$centers)
  CLUST.MELTED$TYPE<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[1])
  CLUST.MELTED$BIN<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[2])
  
  NON_SIG_PLOT<-ggplot(CLUST.MELTED, aes(x=BIN, y=value, col=TYPE)) +geom_line(aes(group=TYPE)) + facet_wrap(~X1) + 
    theme(axis.text.x=element_text(angle=90))
  ggsave(NON_SIG_PLOT, file=paste("NETWORK/PICTURES/METABO_ANALYSIS/CLUSTERING/110513_", cancer, "_NON_SIG_10_CLUSTER2.jpeg", sep=""), dpi=600)
}

TESTING.MATRIX<-as.matrix(HEATMAP1.DATA.SPLIT.LISTS$UCEC$NON_SIG[,5:24])
rownames(TESTING.MATRIX)<-HEATMAP1.DATA.SPLIT.LISTS$UCEC$NON_SIG$PATIENT

CLUST<-kmeans(x=TESTING.MATRIX, 2, nstart=100)
CLUST$size
CLUST.MELTED<-melt(CLUST$centers)
head(CLUST.MELTED)
CLUST.MELTED$TYPE<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[1])
CLUST.MELTED$BIN<-sapply(strsplit(as.character(CLUST.MELTED$X2), "_"), function(l) l[2])
ggplot(CLUST.MELTED, aes(x=BIN, y=value, col=TYPE)) +geom_line(aes(group=TYPE)) + facet_wrap(~X1) + theme(axis.text.x=element_text(angle=90))
CLUST.MELTED