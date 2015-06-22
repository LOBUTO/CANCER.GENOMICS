###120814
#############Develop algorithm to propagate minor allele frequencies based on neightbors######

internal.plotting<-function(GENE, filter=F) {
  EXON.TABLE<-EXONS[Hugo_Symbol==GENE,]
  THOUSAND.TABLE<-THOUSAND[Chrom==unique(EXON.TABLE$Chrom) & unique(EXON.TABLE$START)<Position & Position<unique(EXON.TABLE$END),]
  
  if (filter==T){
    THOUSAND.TABLE<-THOUSAND.TABLE[MAF>=0.000199682,]
  }
  
  #Clean target matrix for position duplicates and aggregate(sum) MAFs contained in duplicates
  target.matrix<-aggregate(MAF~Position,data=THOUSAND.TABLE,FUN=sum)
  target.matrix$TYPE<-"1000G"
  
  #Integrate with TCGA data
  tcga<-TCGA.BRCA.maf[Hugo_Symbol==GENE,]
  tcga<-tcga[,c("Start_Position", "MUT.FREQ"),with=F]
  setnames(tcga, c("Position","MAF"))
  tcga$TYPE<-"TCGA"
  
  #Calculate closest distance
  closest<-sapply(tcga$Position, function(x) min(abs(x-target.matrix$Position)))
  closest<-data.table(Position=tcga$Position, Closest=closest)
  closest<-closest[order(Closest, decreasing=T),]
  print (closest)
  
  target.matrix<-rbind(target.matrix, tcga)
  
  ggplot(target.matrix, aes(Position, MAF, colour=TYPE)) + geom_point() + theme.format + geom_rect(data=EXON.TABLE, aes(xmin=FEAT_START, xmax=FEAT_END, ymin=-Inf, ymax=Inf), fill="yellow", inherit.aes=F, alpha=0.5) 
}

#PREPPING THOUSAND DATA
THOUSAND<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/120814.THOUSAND.SNP.ALL.CNV", header=F, sep="\t", stringsAsFactors=F, drop=6:10)
setnames(THOUSAND, c("Chrom","Position", "TYPE", "REF","ALT", "AF"))
THOUSAND<-THOUSAND[nchar(REF)==1 & nchar(ALT)==1,] #Filter for through SNP (remove insertions and deletions) - Removed 81685538 records ~ 4.21%
THOUSAND<-THOUSAND[AF!=0,] 
THOUSAND<-THOUSAND[TYPE=="SNP",]
THOUSAND$MAF<-ifelse(THOUSAND$AF>0.5, 1-THOUSAND$AF, THOUSAND$AF)

internal.plotting("TTN", filter=T) + scale_x_continuous(limits = c(179500000, 179510000)) + scale_y_log10()

#PREPPING EXON DATA
EXONS<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES", header=T, sep="\t", stringsAsFactors=F))

#PREPPING TCGA DATA
TCGA.BRCA<-unique(fread("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t", stringsAsFactors=F, drop=c(2:4,7:8,10:12,14:15,17:37)))
TCGA.BRCA<-TCGA.BRCA[Variant_Classification %in% c("Missense_Mutation","Silent"),]
ALL.PATIENTS<-length(unique(as.vector(TCGA.BRCA$Tumor_Sample_Barcode)))
TCGA.BRCA.maf<-TCGA.BRCA[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, TYPE="MUT", Sample=Tumor_Sample_Barcode, Chrom=Chrom), 
                         by=c("Start_Position","Hugo_Symbol")]

ggplot(TCGA.BRCA.maf, aes(MUT.FREQ)) + geom_histogram() + scale_x_log10() + scale_y_log10()

###Propagation problem
Function.maf.propagation<-function(thousand.matrix, exons.table, gene, ZERO=F) {
  #Propagate minor allele frequencies (maf) across all gene (including introns) borrowing 
  # information of nearby nucleotides to obtain an aritifical maf in nucleotides without information
  #NOTE: 
  #1. For the time being the interpolation is linear, but, may need to look into exponential interpolation
  #2. We are assuming that start and end are the mRNA start and end, that is they should contain mafs before the first and after the last exon to work, 
  #   we may fix this into looking for the closest known nucleotide with maf information to be more accurate about our calculations
  
  require(data.table)
  require(zoo)
  
  #Prep table
  EXON.TABLE<-exons.table[Hugo_Symbol==gene,]
  chromosome<-unique(EXON.TABLE$Chrom)
  start<-unique(EXON.TABLE$START)
  end<-unique(EXON.TABLE$END)
  target.table<-thousand.matrix[Chrom==chromosome & Position>=start & Position<=end ,c("Position", "MAF"), with=F]
  
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
    target.table$MAF<-na.approx(as.vector(target.table$MAF),na.rm=F,maxgap=245)
    
    #Replace trailing NA obtainig by a linear propagation to unknown maf to minimum maf
    target.table$MAF[is.na(target.table$MAF)]<-0
    
    #Replace negative values produced by interpolation by minimum maf value
    target.table$MAF<-ifelse(target.table$MAF<0, min.maf, target.table$MAF)  
  }
  
  #Clean up and return
  target.table$Chrom<-chromosome
  return(target.table)
}

Function.maf.bayes<-function(thousand.amaf, tcga.maf, gene,cancer.prob=0.12) {
  #Calculate bayesian probability per site using background artifical maf (amaf) applied to tcga samples
  #The default cancer probability is for breast cancer and was obtained from http://www.breastcancer.org/symptoms/understand_bc/risk/understanding
  #This a per gene calculation, make sure tcga maf chrom and position coordinates match to entered thousand.amaf!!
  
  require(data.table)
  
  #Set up tables
  tcga.maf<-tcga.maf[Hugo_Symbol==gene, c("Start_Position", "Chrom", "MUT.FREQ", "Sample", "Hugo_Symbol"), with=F]
  setnames(tcga.maf, c("Position", "Chrom", "MUT.FREQ", "Sample", "Hugo_Symbol"))
  
  target.table<-merge(thousand.amaf, tcga.maf, by=c("Position","Chrom"))
  
  #Calculate the conditional probability for cancer - Probability of having variation given that you have cancer, which is position maf given by MUT.FREQ
  #Calculate the conditional probability for non-cacner - Probability of having variation given that you don't have cancer, which is AMAF given by MAF in thousand 
  
  #Apply bayesian method per position on each patient to obtain posterior probabilities
  target.table$PROB<-(target.table$MUT.FREQ*cancer.prob)/(target.table$MUT.FREQ*cancer.prob  + target.table$MAF*(1-cancer.prob))
  
  #Clean up and return
  target.table<-target.table[order(PROB, decreasing=T),]
  return(target.table)
}

TP53.AMAF<-Function.maf.propagation(THOUSAND, EXONS, "TP53")
TP53.BAYES<-Function.maf.bayes(TP53.AMAF, TCGA.BRCA.maf, "TP53")
hist(TP53.BAYES$PROB)

TTN.BAYES<-Function.maf.bayes(Function.maf.propagation(THOUSAND, EXONS, "TTN"), TCGA.BRCA.maf, "TTN")
PIK3CA.BAYES<-Function.maf.bayes(Function.maf.propagation(THOUSAND, EXONS, "PIK3CA"), TCGA.BRCA.maf, "PIK3CA")
BRCA1.BAYES<-Function.maf.bayes(Function.maf.propagation(THOUSAND, EXONS, "BRCA1"), TCGA.BRCA.maf, "BRCA1")
hist(PIK3CA.BAYES$PROB)


