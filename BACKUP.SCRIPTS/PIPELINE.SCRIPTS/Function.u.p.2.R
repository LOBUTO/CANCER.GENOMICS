#Function.u.p.2.R
#072514
#Calculates a statistic based on the hypergeometric distribution per gene using the background mutation rate per patient and construct a patient x matrix
#Analogous to Function.u.p.R in terms of using the per patient BMR

#It returns a list of:
# STAT<- data.table of genes with paired.t.stat, wilcoxon.t.stat, affected pateints
# GENES<- List of patients per Hugo_Symbol

#LOGIC:
#   ARE WE DIFFIRENT ENOUGH FROM THE ZERO GROUP IN TERMS OF THIS ONE SIGNIFICANT GENE? IF NOT, THEN EVEN THOUGH THIS GENE IS SIGNIFICANT IN ALL MEMBERS OF THE 
#   ONE GROUP, WE SHOULD NOT TAKE IT INTO CONSIDERATION

############################################################################################################################################################

Function.Outliers.Table.1<-function(table.1) {
  
  #Get mutations counts per patient
  table.1.count<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)),by = "Tumor_Sample_Barcode"]
  
  #Calculate non-outlier range
  summary.stats<-summary(table.1.count$TOTAL.MUTATIONS)
  IR.RANGE<-summary.stats[[5]]-summary.stats[[2]]
  LOW.1.5<-summary.stats[[2]]-1.5*IR.RANGE
  HIGH.1.5<-summary.stats[[5]]+1.5*IR.RANGE
  
  #Remove outliers
  table.1.count<-table.1.count[TOTAL.MUTATIONS>=LOW.1.5,]
  table.1.count<-table.1.count[TOTAL.MUTATIONS<=HIGH.1.5,]
  table.1<-table.1[Tumor_Sample_Barcode %in% table.1.count$Tumor_Sample_Barcode,]
  
  #Return filtered table.1
  return(table.1)
}

Function.paired.t.tests.tailored<-function(vector.a, vector.b, ALT="two.sided", MU=0, var=F, paired=F, type="parametric"){
  #Understand CAREFULLY what you are inputting as MU!!
  #RETURNS LIST!!
  #Default for alternative hypothesis is "two.sided"
  
  if (type=="parametric"){
    #pair t.test
    if(  (length(unique(vector.a))>1) | (length(unique(vector.b))>1)  ) {
      dummy.test<-t.test(vector.a, vector.b, var.equal=var, alternative=ALT, mu=MU, paired=paired)  
    }
    
  } else if (type=="non.parametric") {
    #Wilcoxon-signed rank test
    dummy.test<-wilcox.test(vector.a, vector.b, var.equal=var, alternative=ALT, mu=MU, paired=paired)
  }
  
  #Return results of all tests plus number of patients covered by SIGNIFICANCE vector (vector.a) - LENGTH OF VECTOR USED FOR COUNTING COVERAGE!!!!
  return(list(stat=dummy.test$statistic, p.val=dummy.test$p.value,
              PATIENTS.COVERED=length(vector.a)))
}

Function.u.p.2<-function(table.1, length.table, background.sequenced.genes, sns.table, remove.outliers=F) {
  
  require(data.table)
  require(reshape2)
  require(base)
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  table.1<-TESTING$table.1
  length.table<-table.length
  background.sequenced.genes<-dummy.1$background.genes
  sns.table<-AA.RATIO
  
  #Get total count of patients we are starting with
  TOTAL.PATIENT.COUNT<-length(unique(as.vector(table.1$Tumor_Sample_Barcode)))
  
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }

  #Rename tables
  setnames(length.table, colnames(length.table), c("Hugo_Symbol", "Length"))
  length.table<-as.data.table(length.table) #[Hugo_Symbol, Length]
  
  #Prep SNS.TABLE
  SNS.TABLE<-sns.table[,c("Hugo_Symbol", "SNS_SUM_RATIO", "SEQUENCE_LENGTH"),with=F]
  SNS.TABLE$NORM.LENGTH<-normalize.vector(SNS.TABLE$SEQUENCE_LENGTH)
  SNS.TABLE$SEQUENCE_LENGTH<-NULL
  
  #Filter table.1 for genes that contain length information in length.table [Tumor_Sample_Barcode, Hugo_Symbol, N.MUTATIONS]
  table.1<-table.1[Hugo_Symbol %in% unique(as.vector(length.table$Hugo_Symbol)), ]
  
  #Get total background amino acid length (those that could be measured) - USED FOR TOTAL.BMR
  background.length<-length.table[Hugo_Symbol %in% background.sequenced.genes,]
  background.length<-sum(as.vector(background.length$Length)) #Theoretical length of all translatable amino acids added up
  
  #########DO PRE-FILTERING STEP FOR GENES BY SILENT-TO-NON.SILENT PROBABILITY###########
  #This is Mann Whitney U Test normalize of probability ot non-silent to silent pseudo-normalized by length (length in function)
  PRE.FILTER<-as.data.table(merge(as.data.frame(table.1), as.data.frame(SNS.TABLE), by="Hugo_Symbol"))
  PRE.FILTER<-PRE.FILTER[,Function.paired.t.tests.tailored(Missense, Silent, ALT="greater", var=F,paired=F,type="non.parametric",
                                                           MU=mean(Silent)*(1/mean(SNS_SUM_RATIO))*mean(NORM.LENGTH)-mean(Silent)),
                         by="Hugo_Symbol"]
  
  #Threshold to FDR<0.1(Soft thresholding) - This removes things like TTN and long genes with non-larger than expected mutations by SNS Ratio
  PRE.FILTER$ADJ.P.VAL<-p.adjust(PRE.FILTER$p.val, method="fdr")
  PRE.FILTER<-PRE.FILTER[ADJ.P.VAL<0.1,]
  PRE.FILTER[Hugo_Symbol=="CCDN1",]
  #Apply pre-filtered genes to table.1
  table.1<-table.1[Hugo_Symbol %in% as.vector(PRE.FILTER$Hugo_Symbol),]
  ############################################################################
  
  #Filter table.1 for NON-SILENT mutations only - Need to filter out ZERO MISSENSE MUTATIONS FROM INHERENT Table.1
  table.1<-unique(table.1[,1:3, with=F])
  table.1<-table.1[Missense!=0,]
  setnames(table.1, colnames(table.1), c(colnames(table.1)[1:2], "N.MUTATIONS"))
  
  #Get total count of NON-SILENT mutations per patient (meatabolic and non-metabolic) [Tumor_Sample_Barcode, TOTAL.MUTATIONS]
  patient.all.mutations<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)), by="Tumor_Sample_Barcode"]
  
  #CALCULATE
  #[Tumor_Sample_Barcode, Hugo_Symbol, N.MUTATIONS, TOTAL.MUTATIONS, Length]
  HYPER.TABLE<-as.data.table(merge(as.data.frame(table.1), as.data.frame(patient.all.mutations), by="Tumor_Sample_Barcode")) 
  HYPER.TABLE<-as.data.table(merge(as.data.frame(HYPER.TABLE), as.data.frame(length.table), by="Hugo_Symbol"))
  
  #######.1
  HYPER.TABLE.2<-copy(HYPER.TABLE)
  length(unique(as.vector(HYPER.TABLE.2$Tumor_Sample_Barcode))) #Originally 1006 patients, now 1005
  HYPER.TABLE.2$Hugo_Symbol<-as.character(HYPER.TABLE.2$Hugo_Symbol)
  head(as.data.frame(table(HYPER.TABLE.2$Hugo_Symbol)))
  
  #Drop genes that only appear once since we cannot compute t-tests with just one value - FOR.NOW#########
  HYPER.TABLE.2<-as.data.table(merge(as.data.frame(HYPER.TABLE.2), as.data.frame(table(HYPER.TABLE.2$Hugo_Symbol)) , by.x="Hugo_Symbol",by.y="Var1"))
  HYPER.TABLE.2<-HYPER.TABLE.2[Freq>1,] 
  length(unique(as.vector(HYPER.TABLE.2$Tumor_Sample_Barcode))) #Dropped from 1006 to 1005 patients, could complement with CNV data
  HYPER.TABLE.2$Freq<-NULL
  
  #Continue
  HYPER.TABLE.2$GMR<-HYPER.TABLE.2$N.MUTATIONS/HYPER.TABLE.2$Length
  HYPER.TABLE.2$BMR<-HYPER.TABLE.2$TOTAL.MUTATIONS/background.length
  HYPER.TABLE.2.CALC<-HYPER.TABLE.2[,Function.paired.t.tests.tailored(GMR,BMR,ALT="greater",MU=0,var=F,paired=T,type="non.parametric"), by="Hugo_Symbol"]
  
  #Multiple hypothesis testing
  HYPER.TABLE.2.CALC$P.VAL.ADJ<-p.adjust(HYPER.TABLE.2.CALC$p.val,method="fdr")
  hist(HYPER.TABLE.2.CALC$P.VAL.ADJ)
  HYPER.TABLE.2.CALC[Hugo_Symbol %in% COSMIC.BRCA$Symbol,] #ETV6, NTRK3, PBRM1
  ggplot(HYPER.TABLE.2.CALC, aes(stat, -log(P.VAL.ADJ))) + geom_point()
  
  ########
  
  
  #######.1
  
  #Significance per gene in each patient - [Tumor_Sample_Barcode, Hugo_Symbol, P.VAL]
  HYPER.TABLE.CALC<-HYPER.TABLE[,list(P.VAL=phyper(q=N.MUTATIONS-1, m=TOTAL.MUTATIONS, 
                                          n=background.length-TOTAL.MUTATIONS, k= Length, lower.tail=F)), 
                       by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]
  length(unique(as.vector(HYPER.TABLE.CALC$Hugo_Symbol)))
  
  #Correct for multiple hypothesis [Tumor_Sample_Barcode, Hugo_Symbol, P.VAL, ADJ.P.VAL]
  HYPER.TABLE.CALC$P.VAL.ADJ<-p.adjust(HYPER.TABLE.CALC$P.VAL ,method="fdr")
  
  #Make adjusted P-values binaries - Remove non-significant and convert significant to 1
  HYPER.TABLE.CALC<-HYPER.TABLE.CALC[P.VAL.ADJ<0.05,]
  HYPER.TABLE.CALC$P.VAL.ADJ<-1
  length(unique(as.vector(HYPER.TABLE.CALC$Hugo_Symbol)))
  
  #Get number of patients tested in model
  PATIENTS.IN.MODEL<-length(unique(as.vector(HYPER.TABLE$Tumor_Sample_Barcode)))
  
  #Calculate based on BACKGROUND PROBABILITY OF PICKING AS MANY OR MORE PATIENTS AS THE GENE OF INTEREST
  SIG.GENES<-unique(as.vector(HYPER.TABLE.CALC$Hugo_Symbol))
  SIG.TABLE<-HYPER.TABLE.CALC[,list(P.VAL=sum(as.data.table(table(as.vector(HYPER.TABLE[Hugo_Symbol %in% replicate(1000, sample(SIG.GENES,1)),]$Hugo_Symbol)))$N>=sum(P.VAL.ADJ))/1000),
                   by="Hugo_Symbol"]
  hist(SIG.TABLE$P.VAL)
  SIG.TABLE$P.VAL.ADJ<-p.adjust(SIG.TABLE$P.VAL, method="fdr")
  head(SIG.TABLE[order(P.VAL.ADJ),],20)
  
  #Construct binary genes.x.patients matrix
  genes.x.patients<-as.data.table(dcast(HYPER.TABLE.CALC, Hugo_Symbol~Tumor_Sample_Barcode, value.var="P.VAL.ADJ", fill=0))
  
  ########.2
  head(genes.x.patients[,1:4,with=F])
  dim(genes.x.patients)
  testing.2<-as.data.frame(genes.x.patients)
  rownames(testing.2)<-testing.2$Hugo_Symbol
  testing.2$Hugo_Symbol<-NULL
  testing.2.cor<-cor(t(testing.2), method="pearson")
  head(testing.2.cor[,1:4])
  testing.2.cor<-as.data.table(as.table(testing.2.cor))
  setnames(testing.2.cor, colnames(testing.2.cor), c("Hugo_Symbol.1", "Hugo_Symbol.2", "SPEARMAN.RHO"))
  TTN<-testing.2.cor[Hugo_Symbol.1=="TTN",]
  TTN[order(SPEARMAN.RHO,decreasing=T),]
  #########.2
  
  #########.3
  testing.3<-copy(HYPER.TABLE.CALC)
  testing.3.function<-function(x, table.comp) {
    print (x)
    loop.table<-table.comp[,list(MAIN.COVERAGE=length(intersect(as.vector(Tumor_Sample_Barcode), x))/length(x)), by="Hugo_Symbol"]
    setnames(loop.table, colnames(loop.table), c("Hugo_Compared", "Coverage"))
    print(loop.table)
    return(loop.table)
  }
  testing.3<-testing.3[,testing.3.function(as.vector(Tumor_Sample_Barcode), testing.3), by="Hugo_Symbol"]
  dim(testing.3)
  head(testing.3)
  testing.3.TTN<-testing.3[Hugo_Symbol=="TTN",]
  head(testing.3.TTN[order(Coverage, decreasing=T),],10)
  
  testing.3.PIK3CA<-testing.3[Hugo_Symbol=="PIK3CA",]
  testing.3.PIK3CA[order(Coverage, decreasing=T),]
  #########.3
  
  #########.4
  AA.RATIO<-as.data.table(read.csv("DATABASES/UNIPROT/073114_HUMAN_AA_SILENT_RATIO", header=T, sep="\t"))
  hist(AA.RATIO$SNS_MEDIAN_RATIO)
  AA.RATIO[Hugo_Symbol %in% c("TP53", "TTN", "PIK3CA", "CDH1"),]
  
  #Melt matrix for calculations - [Hugo-Symbol, Tumor_Sample_Barcode, SIGNIFICANCE]
  patients.t.genes<-melt(genes.x.patients, id="Hugo_Symbol") #[Hugo_Symbol, variable, value]
  setnames(patients.t.genes, colnames(patients.t.genes), c("Hugo_Symbol","Tumor_Sample_Barcode","SIGNIFICANCE"))
  patients.t.genes<-as.data.table(patients.t.genes) 
  
  #Test  SIGNIFICANT vs "INSIGNIFICANT" vector (Vector of opposite binary values) - [Hugo_Symbol, paired.t.pvals & paired.t.stats , PATIENTS.COVERED]
  PAIRED.SIGNIFICANCE.TABLE<-patients.t.genes[, Function.paired.t.tests.tailored(SIGNIFICANCE,
                                                                                as.numeric(SIGNIFICANCE!=T), ALT="greater", MU=0), by="Hugo_Symbol"]
  head(PAIRED.SIGNIFICANCE.TABLE[order(PATIENTS.COVERED,decreasing=T),],20)
  
  #Do multiple hypothesis correction
  PAIRED.SIGNIFICANCE.TABLE$PAIRED.P.ADJ<-p.adjust(PAIRED.SIGNIFICANCE.TABLE$paired.t.p.val, method="fdr")
  PAIRED.SIGNIFICANCE.TABLE$WILCOXON.P.ADJ<-p.adjust(PAIRED.SIGNIFICANCE.TABLE$wilcoxon.p.val, method="fdr")
  
  #Create lists per each gene model of patients that have significance mutation of gene
  GENE.T.PATIENT.LIST<-lapply(unique(as.vector(PAIRED.SIGNIFICANCE.TABLE$Hugo_Symbol)) , 
                              function(x) unique(as.vector(patients.t.genes[Hugo_Symbol==x, ]$Tumor_Sample_Barcode)), 
                              USE.NAMES=T)
  
  #Clean and return
  PAIRED.SIGNIFICANCE.TABLE<-PAIRED.SIGNIFICANCE.TABLE[order(WILCOXON.P.ADJ),]
  dummy.return<-list(STATS=PAIRED.SIGNIFICANCE.TABLE, GENE.X.PATIENT.LIST<-GENE.T.PATIENT.LIST,
                     STARTING.PATIENTS=TOTAL.PATIENT.COUNT, PATIENTS.IN.MODEL=PATIENT.IN.MODEL)
  
  return(dummy.return)
  
}


table.length<-as.data.table(read.csv("DATABASES/UNIPROT/042314_GENE_LENGTH",header=T, sep="\t"))
test<-Function.u.p.2(dummy.1$table.1, table.length, dummy.1$background.genes)


PAIRED.SIGNIFICANCE.TABLE[order(PATIENTS.COVERED,decreasing=T),]
as.vector(PAIRED.SIGNIFICANCE.TABLE[order(PATIENTS.COVERED,decreasing=T),]$Hugo_Symbol)[1:3]
dummy.1$table.1
length(unique(as.vector(dummy.1$table.1[Hugo_Symbol %in% as.vector(PAIRED.SIGNIFICANCE.TABLE[order(paired.t.stat,decreasing=T),]$Hugo_Symbol)[1:70], ]$Tumor_Sample_Barcode)))

head(HYPER.TABLE.2.CALC[order(P.VAL.ADJ),],20)
length(unique(as.vector(dummy.1$table.1[Hugo_Symbol %in% as.vector(HYPER.TABLE.2.CALC[order(P.VAL.ADJ),]$Hugo_Symbol)[1:8876], ]$Tumor_Sample_Barcode)))

HYPER.TABLE.2.CALC$RANK<--HYPER.TABLE.2.CALC$P.VAL.ADJ
HYPER.TABLE.2.CALC[order(RANK,decreasing=T),]
HYPER.TABLE.2.CALC[order(stat, decreasing=T),]
ggplot(HYPER.TABLE.2.CALC, aes(-log(P.VAL.ADJ), stat)) + geom_point()

Table.HYPER.2.BRCA.u.rank.enrich.p<-Function.p.rank.enrichment(HYPER.TABLE.2.CALC, c("Hugo_Symbol", "RANK"), as.vector(COSMIC.BRCA$Symbol), normalize=F)
Table.HYPER.2.BRCA.u.rank.enrich.p$NON.CUM.TABLE[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
HYPER.TABLE.2.CALC[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
Function.j.rank.enrich.plot(Table.HYPER.2.BRCA.u.rank.enrich.p)
ggplot(Table.HYPER.2.BRCA.u.rank.enrich.p$CUM.TABLE, aes(x=CUM.RECALL, y=CUM.PRECISION)) + geom_line()

log.1<-Table.BRCA.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.1$METHOD<-"BMR"
log.2<-Table.MUTSIG.BRCA.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.2$METHOD<-"MUTSIG"
log.3<-Table.BRCA.v.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.3$METHOD<-"DIFF.EXP.RATIO"
log.4<-Table.HYPER.2.BRCA.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.4$METHOD<-"HYPER.2"
log.5<-testing.enrichment$CUM.TABLE[,c(1,3:6),with=F]
log.5$METHOD<-"HYPER.2.FILTERED"
log.6<-Table.BRCA.v.rank.enrich.p.2$CUM.TABLE[,c(1,3:6), with=F]
log.6$METHOD<-"v.CNV"
log.7<-Factor.enrichment$CUM.TABLE[,c(1,3:6), with=F]
log.7$METHOD<-"FE"

log<-rbind(log.1, log.2,log.3,log.4,log.5, log.6,log.7)

ggplot(log, aes(x=CUM.RECALL, y=CUM.PRECISION, colour=METHOD)) + geom_line() + theme.format + geom_point(size=4)

test<-as.data.table(merge(as.data.frame(HYPER.TABLE.2.CALC[,c(1,5),with=F]), as.data.frame(Table.v.BRCA.p.2), by="Hugo_Symbol",all=T,))
test$SIG<-test$P.VAL.ADJ<0.05
ggplot(test, aes(v.PROTEIN,normalize.vector(-log(P.VAL.ADJ)), colour=SIG)) + geom_point() + theme.format

#DEVELOP RE-RANKING METHOD THAT COMBINES v(p) and u(p) in order to select a phenotype change threshold "v(p)" that allows for filtering out of u(p) to increase PR curve
test[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
testing<-test[v.PROTEIN>0.63,]
testing.enrichment<-Function.p.rank.enrichment(testing, c("Hugo_Symbol", "RANK"), as.vector(COSMIC.BRCA$Symbol), normalize=T)
testing[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]