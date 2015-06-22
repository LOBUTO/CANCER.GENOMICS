library(data.table)
library(exactRankTests)
library(reshape2)

test<-data.table(labels=seq(10,100,by=10), samples=seq(10,100,by=10))
test


dummy.normal<-sample(1:100, 50 )
lm.log.function<-function(size){
    
  t.vector<-sample(1:300, size)
  #test.vector<-c(dummy.normal, t.vector)
  #model<-c(rep(0,length(dummy.normal)) , rep(1,length(t.vector)))
  #model.lm<-glm(model~test.vector, family=binomial())
  
  #p.val<-summary(model.lm)$coefficients[2,4]
  p.val<-wilcox.exact(dummy.normal, t.vector)$p.value
  return(p.val)
}
test.lm<-test[ , list(p.val=replicate(1000,lm.log.function(samples)))  , by="labels"]
test.lm$adj.p.val<-p.adjust(test.lm$p.val, method="fdr")
test.lm$NEG.LOG<--log(test.lm$adj.p.val)
test.lm$SIG<-test.lm$adj.p.val<0.05

library(ggplot2)
ggplot(test.lm, aes(factor(labels), NEG.LOG, colour=SIG)) + geom_boxplot()

size=50
model
plot(model, test.vector)

#TTN test
Function.RNAseq.Differential.Expression.V2<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  
  require(limma)
  require(data.table)
  
  #Get target combined matrix
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  target.combined.matrix<-normalized.matrices.object$combined.matrices[,c(cancer.samples.in.matrix,normal.samples)]
  
  #Get design matrix
  G.design.matrix<-data.frame(G=c(rep("G1", length(cancer.samples.in.matrix)), rep("G0", length(normal.samples))))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Perform differential expression
  G.fit = lmFit(target.combined.matrix, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Convert to data.table
  if ("ID" %in% colnames(all.G.fit)) {
    all.G.fit<-as.data.table(all.G.fit)
  } else {
    all.G.fit$ID<-rownames(all.G.fit)
    all.G.fit<-as.data.table(all.G.fit)
  }
  
  #Return topTable
  return(all.G.fit)
}

exp.obj<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/102514.CANCER.MATRICES.NORMALIZED.OBJ.rds")
exp.obj$combined.matrices["COL10A1",1090:1100]
median(exp.obj$combined.matrices["COL10A1",1:1092])
median(exp.obj$combined.matrices["COL10A1",1093:1203])

table.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/100514.BRCA.Table.1.rds")
table.1$table.1<-table.1$table.1[PATIENT %in% exp.obj$cancer.patients,][Missense!=0,]

TTN.test<-Function.RNAseq.Differential.Expression.V2(exp.obj, as.vector(table.1$table.1[Hugo_Symbol=="TTN",]$PATIENT))
TTN.test ####Check if we are looking at logFC with respect to normal or tumor (+ or - logFC) - CHECKED!
TTN.test[adj.P.Val<0.05,]
TTN.test[abs(logFC)>=1,]
hist(TTN.test$adj.P.Val)

TP53.test<-Function.RNAseq.Differential.Expression.V2(exp.obj, as.vector(table.1$table.1[Hugo_Symbol=="TP53",]$PATIENT))
TP53.test
TP53.test[adj.P.Val<0.05,]
TP53.test[abs(logFC)>=1,]
hist(TP53.test$adj.P.Val)

#Check if cnv correlates with mutations
gistic<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/092214.BRCA.GISTIC.TH.2.rds")
table.1$table.1

cnv.correlation<-function(patients, gistic.table) {
  require(data.table)
  require(exactRankTests)
  
  internal.function<-function(target.patients,cnv){
    all<-data.table(TARGET=target.patients, CNV=cnv)
    TARGET.CNV<-as.vector(all[TARGET==TRUE,]$CNV)
    NON.TARGET.CNV<-as.vector(all[TARGET==FALSE,]$CNV)
    
    if( (length(TARGET.CNV)!=0) & (length(TARGET.CNV)!=0) ){
      test<-wilcox.exact(TARGET.CNV, NON.TARGET.CNV)
      P.VAL=test$p.value
    } else {
      P.VAL=1
    }

    return(list(P.VAL=P.VAL))    
  }
  
  gistic.table$TARGET.PATIENT<-gistic.table$PATIENT %in% patients
  gistic.table<-gistic.table[,internal.function(TARGET.PATIENT, CNV.TH), by="Hugo_Symbol"]
  
  count<<-count+1
  print (count)
  return(list(HUGO.CNV=as.vector(gistic.table$Hugo_Symbol), P.VAL=as.vector(gistic.table$P.VAL)))
}

length(unique(as.vector(table.1$table.1$Hugo_Symbol)))
count<<-0

cnv.mut.cor<-table.1$table.1[Hugo_Symbol %in% c("TP53", "TTN","BRCA1","COL10A1","PIK3CA", "CDH1", "GATA3","PTEN"),][,cnv.correlation(PATIENT, gistic), by="Hugo_Symbol"]
cnv.mut.cor$P.VAL.ADJ<-p.adjust(cnv.mut.cor$P.VAL, method="fdr")
cnv.mut.cor$SIG<-cnv.mut.cor$P.VAL.ADJ<0.05
ggplot(cnv.mut.cor, aes(x=Hugo_Symbol, y=-log(P.VAL.ADJ), colour=SIG)) + geom_boxplot()
ggplot(cnv.mut.cor, aes(SIG, fill=Hugo_Symbol)) + geom_histogram(position="dodge")

install.packages("cgdsr")
library(cgdsr)
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/public-portal/")

test(mycgds)

# Get list of cancer studies at server
getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[14,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[5,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[1,1]

# Get data slices for a specified list of genes, genetic profile and case list
CNA<-getProfileData(mycgds,c('BRCA1','TTN','TP53','COL10A1'),mygeneticprofile,mycaselist)
CNA<-melt(as.data.table(CNA))
ggplot(CNA, aes(value, fill=variable)) + geom_histogram(position="dodge")

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)

getMutationData(mycgds,"brca_tcga_all", "brca_tcga_mutations", c("BRCA1", "BRCA2"))

####Is there a correlation between mutations and copy number variation for genes

####111114
###Mapping tp53 from 1000Genome processed data
thousand<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.SNP.CNV",header=T,sep="\t",stringsAsFactors=F)
thousand

#Get TP53 only
thousand.TP53<-thousand[Chrom=="17" & 7571719<Position & Position<7590868,]
thousand.TP53$MAF<-ifelse(thousand.TP53$AF>0.5, 1-thousand.TP53$AF, thousand.TP53$AF)
ggplot(thousand.TP53, aes(Position, MAF)) + geom_histogram(stat="identity", position="dodge")
ggplot(thousand.TP53, aes(Position, MAF)) + geom_line() + theme.format

#Map exons
EXONS<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES", header=T, sep="\t", stringsAsFactors=F))
EXONS[Hugo_Symbol=="TP53",]

ggplot(thousand.TP53, aes(Position, MAF)) + geom_point() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) 
ggplot(thousand.TP53, aes(Position, MAF)) + geom_line() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) 

#Map somatic mutations from TCGA
TCGA.maf<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t", stringsAsFactors=F))

ALL.PATIENTS<-length(unique(as.vector(TCGA.maf$Tumor_Sample_Barcode)))

TCGA.maf.TP53<-TCGA.maf[Hugo_Symbol=="TP53",]
TCGA.maf.TP53<-TCGA.maf.TP53[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, TYPE="MUT"), by=c("Start_Position")]

#Combine TCGA and THOUSAND for 1000genome
setnames(TCGA.maf.TP53, c("Position", "MAF","TYPE"))
TP53<-rbind(thousand.TP53[,c("Position", "MAF", "TYPE"), with=F], TCGA.maf.TP53)

ggplot(TP53, aes(Position, MAF, colour=TYPE)) + geom_histogram(stat="identity") + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="lightblue", inherit.aes=F, alpha=0.3) +
  scale_x_continuous(limits = c(7570000, 7580000))

ggplot(TP53, aes(Position, MAF, colour=TYPE)) + geom_histogram(stat="identity") + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="lightblue", inherit.aes=F, alpha=0.3) +
  scale_x_continuous(limits = c(7572500, 7580000))

ggplot(TP53[MAF<0.2,], aes(Position, MAF, colour=TYPE)) + geom_point() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.4) +
  scale_x_continuous(limits = c(7578000, 7580000))

#### TTN 
thousand.TTN<-thousand[Chrom=="2" & 179390715<Position & Position<179672150,]
thousand.TTN$MAF<-ifelse(thousand.TTN$AF>0.5, 1-thousand.TTN$AF, thousand.TTN$AF)
thousand.TTN<-thousand.TTN[TYPE=="SNP",]

TCGA.maf.TTN<-TCGA.maf[Hugo_Symbol=="TTN",]
TCGA.maf.TTN<-TCGA.maf.TTN[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, TYPE="MUT"), by=c("Start_Position")]

#Combine TCGA and THOUSAND for 1000genome
setnames(TCGA.maf.TTN, c("Position", "MAF","TYPE"))
TTN<-rbind(thousand.TTN[,c("Position", "MAF", "TYPE"), with=F], TCGA.maf.TTN)

ggplot(TTN, aes(Position, MAF, colour=TYPE)) + geom_histogram(stat="identity", position="dodge") + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TTN",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) +
  scale_x_continuous(limits = c(179600000, 179625000))
ggplot(TTN, aes(Position, MAF, colour=TYPE)) + geom_point() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TTN",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.4) +
 scale_x_continuous(limits = c(179610000, 179617500))

#### PIK3CA 
thousand.PIK3CA<-thousand[Chrom=="3" & 178866310<Position & Position<178952500,]
thousand.PIK3CA$MAF<-ifelse(thousand.PIK3CA$AF>0.5, 1-thousand.PIK3CA$AF, thousand.PIK3CA$AF)
thousand.PIK3CA<-thousand.PIK3CA[TYPE=="SNP",]

ggplot(thousand.PIK3CA, aes(Position, MAF)) + geom_line() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="PIK3CA",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="lightblue", inherit.aes=F, alpha=0.5) 

TCGA.maf.PIK3CA<-TCGA.maf[Hugo_Symbol=="PIK3CA",]
TCGA.maf.PIK3CA<-TCGA.maf.PIK3CA[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, TYPE="MUT"), by=c("Start_Position")]

#Combine TCGA and THOUSAND for 1000genome
setnames(TCGA.maf.PIK3CA, c("Position", "MAF","TYPE"))
PIK3CA<-rbind(thousand.PIK3CA[,c("Position", "MAF", "TYPE"), with=F], TCGA.maf.PIK3CA)

ggplot(PIK3CA, aes(Position, MAF, colour=TYPE)) + geom_histogram(stat="identity") + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="PIK3CA",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) 

ggplot(PIK3CA, aes(Position, MAF, colour=TYPE)) + geom_histogram(stat="identity") + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="PIK3CA",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) +
  scale_x_continuous(limits = c(178910000, 178955000))

ggplot(PIK3CA, aes(Position, MAF, colour=TYPE)) + geom_point() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="PIK3CA",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) +
  scale_x_continuous(limits = c(178925000, 178930000))

###Develop function to generate artificial mafs in thousand gene matrix
thousand.TP53[,c("Position", "MAF"), with=F]

Function.maf.propagation<-function(thousand.matrix, start, end, ZERO=F) {
  #Propagate minor allele frequencies (maf) across all gene (including introns) borrowing 
  # information of nearby nucleotides to obtain an aritifical maf in nucleotides without information
  #NOTE: 
  #1. For the time being the interpolation is linear, but, may need to look into exponential interpolation
  #2. We are assuming that start and end are the mRNA start and end, that is they should contain mafs before the first and after the last exon to work, 
  #   we may fix this into looking for the closest known nucleotide with maf information to be more accurate about our calculations
  
  require(data.table)
  require(zoo)
  
  #Prep table
  chromosome<-unique(thousand.matrix$Chrom)
  target.table<-thousand.matrix[,c("Position", "MAF"), with=F]
  
  #Obtain minimum maf
  min.maf<-min(as.vector(target.table$MAF))
  
  #Build artifical maf matrix of not known positions
  art.pos<-setdiff(start:end, as.vector(target.table$Position))
  art.matrix<-data.table(Position=art.pos, MAF=NA)
  
  #Integrate to real maf and order for propragation
  target.table<-rbind(target.table, art.matrix)
  target.table<-target.table[order(Position),]
  
  #If interpolation, apply method, otherwise fill with zeroes
  if (ZERO==TRUE){
    target.table$MAF[is.na(target.table$MAF)]<-0
  } else {
    #Propragate maf using cubic spline interpolation in the package zoo - MAXGAP WILL BE FUNCTION OF GENE LENGTH!!! -LONGER GENE, MORE OF A GAP TO FILL OUT, LESS INFO WE HAVE
    target.table$MAF<-na.approx(as.vector(target.table$MAF),na.rm=F,maxgap=250)
    
    #Replace trailing NA obtainig by a linear propagation to unknown maf to minimum maf
    target.table$MAF[is.na(target.table$MAF)]<-0
    
    #Replace negative values produced by interpolation by minimum maf value
    target.table$MAF<-ifelse(target.table$MAF<0, min.maf, target.table$MAF)  
  }
  
  #Clean up and return
  target.table$Chrom<-chromosome
  return(target.table)
}

thousand.TP53.AMAF<-Function.maf.propagation(thousand.TP53, 7571719, 7590868)

ggplot(thousand.TP53.AMAF, aes(Position, MAF)) + geom_point() + theme.format + geom_rect(data=EXONS[Hugo_Symbol=="TP53",], aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.4) 

###Now, develop function to assign probability per mutation through bayesian framework
thousand.TP53.AMAF
TCGA.maf.TP53

Function.maf.bayes<-function(thousand.amaf, tcga.maf, cancer.prob=0.12) {
  #Calculate bayesian probability per site using background artifical maf (amaf) applied to tcga samples
  #The default cancer probability is for breast cancer and was obtained from http://www.breastcancer.org/symptoms/understand_bc/risk/understanding
  #This a per gene calculation, make sure tcga maf chrom and position coordinates match to entered thousand.amaf!!
  
  require(data.table)
  
  #Set up tables
  setnames(tcga.maf, c("Position", "Chrom", "MUT.FREQ", "Sample"))
  target.table<-merge(thousand.amaf, tcga.maf, by=c("Position","Chrom"))
  
  #Calculate the conditional probability for cancer - Probability of having variation given that you have cancer, which is position maf given by MUT.FREQ
  #Calculate the conditional probability for non-cacner - Probability of having variation given that you don't have cancer, which is AMAF given by MAF in thousand 
  
  #Apply bayesian method per position on each patient to obtain posterior probabilities
  target.table$PROB<-(target.table$MUT.FREQ*cancer.prob)/(target.table$MUT.FREQ*cancer.prob  + target.table$MAF*(1-cancer.prob))
  
  #Clean up and return
  target.table<-target.table[order(PROB, decreasing=T),]
  return(target.table)
}

ALL.PATIENTS<-length(unique(as.vector(TCGA.maf$Tumor_Sample_Barcode)))
#For TP53
TCGA.maf.TP53<-TCGA.maf[Hugo_Symbol=="TP53",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
TCGA.maf.TP53<-TCGA.maf.TP53[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.TP53<-Function.maf.bayes(thousand.TP53.AMAF, TCGA.maf.TP53, 0.124)
hist(BAYES.TP53$PROB)

#For TTN
thousand.TTN<-thousand[Chrom=="2" & 179390715<Position & Position<179672150,]
thousand.TTN$MAF<-ifelse(thousand.TTN$AF>0.5, 1-thousand.TTN$AF, thousand.TTN$AF)
thousand.TTN<-thousand.TTN[TYPE=="SNP",]

thousand.TTN.AMAF<-Function.maf.propagation(thousand.TTN, 179390715,179672150,ZERO=FALSE)

TCGA.maf.TTN<-TCGA.maf[Hugo_Symbol=="TTN",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
TCGA.maf.TTN<-TCGA.maf.TTN[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.TTN<-Function.maf.bayes(thousand.TTN.AMAF, TCGA.maf.TTN, 0.124)
hist(BAYES.TTN$PROB)

#For PIK3CA
thousand.PIK3CA.AMAF<-Function.maf.propagation(thousand.PIK3CA, 178866310,178952500)

TCGA.maf.PIK3CA<-TCGA.maf[Hugo_Symbol=="PIK3CA",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
TCGA.maf.PIK3CA<-TCGA.maf.PIK3CA[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.PIK3CA<-Function.maf.bayes(thousand.PIK3CA.AMAF, TCGA.maf.PIK3CA, 0.124)
hist(BAYES.PIK3CA$PROB)

#For MUC4
thousand.MUC4<-thousand[Chrom=="3" & 195473635<Position & Position<195538844,]
thousand.MUC4$MAF<-ifelse(thousand.MUC4$AF>0.5, 1-thousand.MUC4$AF, thousand.MUC4$AF)
thousand.MUC4<-thousand.MUC4[TYPE=="SNP",]

thousand.MUC4.AMAF<-Function.maf.propagation(thousand.MUC4, 195473635,195538844)

TCGA.maf.MUC4<-TCGA.maf[Hugo_Symbol=="MUC4",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
MUC4.PATIENTS<-length(unique(as.vector(TCGA.maf.MUC4$Tumor_Sample_Barcode)))
TCGA.maf.MUC4<-TCGA.maf.MUC4[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.MUC4<-Function.maf.bayes(thousand.MUC4.AMAF, TCGA.maf.MUC4, 0.124)
hist(BAYES.MUC4$PROB)

#For BRCA1
thousand.BRCA1<-thousand[Chrom=="17" & 41196311<Position & Position<41277500,]
thousand.BRCA1$MAF<-ifelse(thousand.BRCA1$AF>0.5, 1-thousand.BRCA1$AF, thousand.BRCA1$AF)
thousand.BRCA1<-thousand.BRCA1[TYPE=="SNP",]

thousand.BRCA1.AMAF<-Function.maf.propagation(thousand.BRCA1, 41196311,41277500)

TCGA.maf.BRCA1<-TCGA.maf[Hugo_Symbol=="BRCA1",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
BRCA1.PATIENTS<-length(unique(as.vector(TCGA.maf.BRCA1$Tumor_Sample_Barcode)))
TCGA.maf.BRCA1<-TCGA.maf.BRCA1[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.BRCA1<-Function.maf.bayes(thousand.BRCA1.AMAF, TCGA.maf.BRCA1, 0.124)
hist(BAYES.BRCA1$PROB)

#For CDH1
thousand.CDH1<-thousand[Chrom=="16" & 68771194<Position & Position<68869444,]
thousand.CDH1$MAF<-ifelse(thousand.CDH1$AF>0.5, 1-thousand.CDH1$AF, thousand.CDH1$AF)
thousand.CDH1<-thousand.CDH1[TYPE=="SNP",]

thousand.CDH1.AMAF<-Function.maf.propagation(thousand.CDH1, 68771194,68869444,ZERO=TRUE)

TCGA.maf.CDH1<-TCGA.maf[Hugo_Symbol=="CDH1",][Variant_Classification %in% c("Missense_Mutation", "Silent"),]
CDH1.PATIENTS<-length(unique(as.vector(TCGA.maf.CDH1$Tumor_Sample_Barcode)))
TCGA.maf.CDH1<-TCGA.maf.CDH1[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, Sample=Tumor_Sample_Barcode), by=c("Start_Position", "Chrom")]
BAYES.CDH1<-Function.maf.bayes(thousand.CDH1.AMAF, TCGA.maf.CDH1, 0.124)
hist(BAYES.CDH1$PROB)

######111414

#Build Likelihood ratio per gene across sites
BAYES.TP53
unique(BAYES.TP53[,1:4,with=F])
Function.Bayes.Posterior.Odds.Hugo<-function(hugo.bayes.table, cancer.prob=0.124){
  #Obtains Posterior Odds per gene 
  
  require(data.table)
  
  #Prep table - We have to obtain the per position probabilities
  hugo.bayes.table<-unique(hugo.bayes.table[,c(1:4,6), with=F])
  
  #Calculate prior odds (Probability of having cancer divided by probability of not having it)
  prior.odds<-cancer.prob/(1-cancer.prob)
  
  #Calculate the likelihood ratio (Ratio of multiplied conditionally independent probabilities of sites given cancer to sites given not cancer)
  #NOTE: Need to split vector length in order to avoid finite zero numbers
  breaks<-ifelse(nrow(hugo.bayes.table)>=200,20,
                 ifelse(nrow(hugo.bayes.table)>=100, 10, 
                 ifelse(nrow(hugo.bayes.table)>=50, 5, 2)))

  
  prob.splits<-split(1:nrow(hugo.bayes.table), 1:breaks)
  
  likelihood.ratio<-lapply(prob.splits, function(x) signif(prod(as.vector(hugo.bayes.table$MUT.FREQ)[x]))/signif(prod(as.vector(hugo.bayes.table$MAF)[x])) )
  likelihood.ratio<-do.call(prod, likelihood.ratio)
  
  #Calculate posterior odds
  Posterior.odds<-prior.odds*likelihood.ratio
  
  #Calculate string integration (noisy - OR)
  link<-lapply(prob.splits, function(x) signif(1-as.vector(hugo.bayes.table$PROB)[x]))
  link<-do.call(prod, link)
  noisy.OR<-signif(1-link)
  
  #Return
  return(data.table(Posterior.Odds=Posterior.odds, noisy.OR=noisy.OR))
}

signif(1-1.16361e-08)
Function.Bayes.Posterior.Odds.Hugo(BAYES.TP53)
Function.Bayes.Posterior.Odds.Hugo(BAYES.TTN)
Function.Bayes.Posterior.Odds.Hugo(BAYES.BRCA1)
Function.Bayes.Posterior.Odds.Hugo(BAYES.CDH1)
Function.Bayes.Posterior.Odds.Hugo(BAYES.PIK3CA)

1-(prod(1-rep(0.2,10)))

plot (1:200, sapply(1:200, function(x)  1-(prod(1-rep(0.01,x))) ))
plot (1:200, sapply(1:200, function(x)  1-(prod(1-rep(0.05,x))) ))
plot (1:200, sapply(1:200, function(x)  1-(prod(1-rep(0.2,x))) ))
plot (1:200, sapply(1:200, function(x)  1-(prod(1-rep(0.5,x))) ))
plot (1:200, sapply(1:200, function(x)  1-(prod(1-rep(0.7,x))) ))

ggplot(melt(unique(BAYES.TTN[,1:4,with=F]), id=c("Position","Chrom")), aes(variable, value)) + geom_boxplot() + scale_y_log10()

wilcox.test(unique(BAYES.TTN[,1:4,with=F])$MUT.FREQ, unique(BAYES.TTN[,1:4,with=F])$MAF, alternative="greater", paired=T)


#Do a patient at a time and update the expection based on the number of patients we accumulate? This will affect the conditional probability until it reaches its true value, the point is that given an

#This is the score for a gene: Continously bayes model site per site based on rank of highest individual to lowest. We count the number of sites we were able to obtain (joint probability of sites) until it passes a low threshold. then divide this number by total number of sites. This is the probability the highest probability that we will obtain out of all (kind of like a random walk, GSEA)