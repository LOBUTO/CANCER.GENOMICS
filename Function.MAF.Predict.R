#Function.MAF.Predict.R
#121714
#Develop and algorithm to predict MAF potentially based on:
#   Codon context
#   Phastcon (?)
#   Gene variation score (sum of variation for all gene (exons or exons+introns))
#   Make sure to filter out small MAF (close to noise)

THOUSAND
normalize.vector<-function(x){
  y=(x-min(x))/(max(x)-min(x))
  return(y)
} 

Function.Gene.MAF.Variation<-function(thousand.table, exons.table, gene){
  #Function calculates the sum of MAFs in whole gene or exon divided by nt length
  
  require(data.table)
  require(IRanges)
  
  #Calculate total exon nt count
  gene.exon<-exons.table[Hugo_Symbol==gene,]
  gene.ranges<-IRanges(gene.exon$FEAT_START, gene.exon$FEAT_END)
  gene.ranges<-reduce(gene.ranges)
  gene.exon.nt<-sum(end(gene.ranges) - start(gene.ranges))
  
  #Calculate total gene nt count
  gene.total.nt<-unique(gene.exon$END)-unique(gene.exon$START)
  
  #Get processed thousand table
  start<-unique(gene.exon$START)
  end<-unique(gene.exon$END)
  chrom<-unique(gene.exon$Chrom)
  target.matrix<-thousand.table[Chrom==chrom & Position>start & Position<end,]
  target.matrix<-aggregate(MAF~Position,data=target.matrix,FUN=sum)
  
  #Filter low count MAF (1/5008)
  print (as.data.table(target.matrix)[order(MAF),])
  target.matrix<-target.matrix[target.matrix$MAF>=0.000199682,]
  print (as.data.table(target.matrix)[order(MAF),])
  
  #Calculate maf coverage sum for exons
  mafs.exon.sum<-apply(gene.exon[,c("FEAT_START", "FEAT_END"), with=F], 1, 
                      function(x) as.vector(target.matrix$Position) %in% x[1]:x[2])
  
  mafs.exon.sum<-apply(mafs.exon.sum, 1, sum)
  target.matrix$exon.maf<-mafs.exon.sum
  target.matrix<-target.matrix[target.matrix$exon.maf==T,]
  mafs.exon.sum<-sum(target.matrix$MAF)
  
  exon.maf.ratio<-mafs.exon.sum/gene.exon.nt
  
  #Calculate maf coverage sum for whole gene 
  gene.maf.ratio<-sum(target.matrix$MAF)/gene.total.nt

  #Return as list
  return(list(exon.maf.ratio=exon.maf.ratio, gene.maf.ratio=gene.maf.ratio))  
}

#Look for an average maf ratio sample from the population
cov.sum.maf<-data.frame(a=c(), b=c(), c=c())
for (gene in c("TP53", "TTN", "PIK3CA", "BRCA1", "CDH1", "ZZZ3")){
  print (gene)
  gene.list<-Function.Gene.MAF.Variation(THOUSAND, EXONS, gene)
  print (c(gene, gene.list$exon.maf.ratio, gene.list$gene.maf.ratio))
  cov.sum.maf<-rbind(cov.sum.maf, c(gene, gene.list$exon.maf.ratio, gene.list$gene.maf.ratio))
  #cov.sum.maf<-c(cov.sum.maf,list(c(gene, gene.list$exon.maf.ratio, gene.list$gene.maf.ratio)))
}
cov.sum.maf<-as.data.table(cov.sum.maf)
setnames(cov.sum.maf, c("Hugo_Symbol", "EXON.MAF.RATIO", "GENE.MAF.RATIO"))

ggplot(unique(cov.sum.maf), aes(Hugo_Symbol, as.numeric(EXON.MAF.RATIO))) + geom_histogram(stat="identity") + theme.format + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_text(size=10)) + coord_flip() 
ggplot(unique(cov.sum.maf), aes(as.numeric(GENE.MAF.RATIO), as.numeric(EXON.MAF.RATIO))) + geom_point() + theme.format
median(cov.sum.maf$EXON.MAF.RATIO)

######120814######
#Create a function to optimize the function that classifies most mutations found in TP53,BRCA1, PIK3CA and CDH1 as cancer and TTN,ZZZ3, A1CF and MUC4.
#   This function aims to increase the difference between cancer-genes mutations to their MAFs and reducing into between TTN,ZZZ3 and MUC.
#   This done by optimizing the predicted TP53 MAFs
Function.MAF.Mapping<-function(thousand.table, exons.table, gene, filter=F){
  #Map MAFs to gene at the exon and gene level
  
  require(data.table)
  require(IRanges)
  
  #####Map MAFs to EXONS ranges first######
  gene.exon<-exons.table[Hugo_Symbol==gene,]
  exon.ranges<-IRanges(gene.exon$FEAT_START, gene.exon$FEAT_END)
  exon.ranges<-reduce(exon.ranges)
  exon.ranges<-data.table(FEAT_START=start(exon.ranges), FEAT_END=end(exon.ranges))  
  
  #Get processed thousand table
  start<-unique(gene.exon$START)
  end<-unique(gene.exon$END)
  chrom<-unique(gene.exon$Chrom)
  target.matrix<-thousand.table[Chrom==chrom & Position>=start & Position<=end,]
  target.matrix<-as.data.table(aggregate(MAF~Position,data=target.matrix,FUN=sum))
  
  #Filter low count MAF (1/5008) if desired
  if (filter==T){
    target.matrix<-target.matrix[MAF>=0.000199682,]  #<-- MAFs at gene level
  }
  gene.mafs<-copy(target.matrix)
  
  #Obtain MAFs in exon ranges
  mafs.exon.sum<-apply(exon.ranges[,c("FEAT_START", "FEAT_END"), with=F], 1, 
                       function(x) as.vector(target.matrix$Position) %in% x[1]:x[2])
  
  mafs.exon.sum<-apply(mafs.exon.sum, 1, sum)
  target.matrix$exon.maf<-mafs.exon.sum
  target.matrix<-target.matrix[exon.maf==T,] #<--MAFs at exon level
  
  #Return
  return(list(gene.mafs=gene.mafs, exon.mafs=target.matrix))    
}

TP53.MAP<-Function.MAF.Mapping(THOUSAND,EXONS,"TP53")
TTN.MAP<-Function.MAF.Mapping(THOUSAND,EXONS,"TTN",filter=T)

#Look at MAF trimers
MAF.TRIMERS<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/121914.THOUSAND.SNP.TRIMERS", header=T, sep="\t", stringsAsFactors=F)
MAF.TRIMERS<-MAF.TRIMERS[AF!=0,] 
MAF.TRIMERS$MAF<-ifelse(MAF.TRIMERS$AF>0.5, 1-MAF.TRIMERS$AF, MAF.TRIMERS$AF)
MAF.TRIMERS<-MAF.TRIMERS[MAF>=0.000199682,] #Filter smallest counts
ggplot(MAF.TRIMERS, aes(TRIMER, MAF)) + geom_boxplot() + theme.format + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggplot(MAF.TRIMERS[TRIMER=="ACG",], aes(Chrom, MAF)) + geom_boxplot() + theme.format + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

#Now look at exon primer only
MAF.TRIMERS.EXON<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/122014.THOUSAND.SNP.TRIMERS.EXON.FILT", header=T, sep="\t", stringsAsFactors=F)
MAF.TRIMERS.EXON<-MAF.TRIMERS.EXON[AF!=0,] 
MAF.TRIMERS.EXON$MAF<-ifelse(MAF.TRIMERS.EXON$AF>0.5, 1-MAF.TRIMERS.EXON$AF, MAF.TRIMERS.EXON$AF)
MAF.TRIMERS.EXON<-MAF.TRIMERS.EXON[MAF>=0.000199682,] #Filter smallest counts
ggplot(MAF.TRIMERS.EXON, aes(TRIMER, MAF)) + geom_boxplot() + theme.format + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggplot(MAF.TRIMERS.EXON[Chrom %in% c("1","5","X","Y"),], aes(TRIMER, MAF)) + geom_boxplot() + theme.format + scale_y_log10() +
   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(Chrom~.) 

#####Load processed trimers with temp and phastcon scores
MAF.TRIMERS.EXON.PHAST<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST", sep="\t", stringsAsFactors=F)
setnames(MAF.TRIMERS.EXON.PHAST, c("Chrom", "Position", "REF","ALT","TRIMER","MT","EXON","AF", "PHAST.100"))
MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[AF!=0,] 
MAF.TRIMERS.EXON.PHAST$MAF<-ifelse(MAF.TRIMERS.EXON.PHAST$AF>0.5, 1-MAF.TRIMERS.EXON.PHAST$AF, MAF.TRIMERS.EXON.PHAST$AF)
MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[!is.na(MT),]#clear of potential NXN trimers MT of NA!!!
MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[MAF>=0.000199682,] #Filter smallest counts

######122914######
#Testing linear models
MAF.TRIMERS.EXON.PHAST$Chrom<-as.factor(MAF.TRIMERS.EXON.PHAST$Chrom)
MAF.TRIMERS.EXON.PHAST$TRIMER<-as.factor(MAF.TRIMERS.EXON.PHAST$TRIMER)
MAF.TRIMERS.EXON.PHAST$REF<-as.factor(MAF.TRIMERS.EXON.PHAST$REF)
MAF.TRIMERS.EXON.PHAST$ALT<-as.factor(MAF.TRIMERS.EXON.PHAST$ALT)
MAF.TRIMERS.EXON.PHAST<-MAF.TRIMERS.EXON.PHAST[EXON==TRUE,]
mod<-lm(formula=MAF~exp(PHAST.100)+REF+ALT+MT+TRIMER+Chrom, data=MAF.TRIMERS.EXON.PHAST)
summary(mod)

test.copy<-copy(MAF.TRIMERS.EXON.PHAST)
test.copy$predict.maf<-predict(mod, test.copy)


#Test the model
set.seed(42)

training.sample<-sample(nrow(MAF.TRIMERS.EXON.PHAST), 750000) #70%
testing.sample<-setdiff(1:nrow(MAF.TRIMERS.EXON.PHAST), training.sample) #10%
training.sample<-MAF.TRIMERS.EXON.PHAST[training.sample,]
testing.sample<-MAF.TRIMERS.EXON.PHAST[testing.sample,]

mod.1<-lm(formula=MAF~exp(PHAST.100)+REF+ALT+MT+TRIMER+Chrom, data=training.sample)
summary(mod.1)
testing.sample$predict.maf<-predict(mod.1, testing.sample)
testing.sample$predict.score<-ifelse(testing.sample$MAF>testing.sample$predict.maf, 
                                     testing.sample$predict.maf/testing.sample$MAF, testing.sample$MAF/testing.sample$predict.maf)
ggplot(testing.sample, aes(predict.score)) + geom_histogram() + theme.format 

#Model based on lm won't work, try decission trees
library(rpart)
m.rpart<-rpart(MAF~., data=training.sample[,c("PHAST.100", "REF","ALT","TRIMER","MT","Chrom", "MAF"),with=F])
summary(m.rpart)
testing.sample$predict.maf<-predict(m.rpart, testing.sample)

#Model based on model trees
library(XLConnect)
xlcFreeMemory()
options( java.parameters = "-Xmx16384m" )
library(RWeka)
m.m5p<-M5P(MAF~., data=training.sample[,c("PHAST.100", "REF","ALT","TRIMER","MT","Chrom", "MAF"),with=F])

library(tree)
m.tree<-tree(MAF~exp(PHAST.100)+REF+ALT+MT+Chrom, data=training.sample[,c("PHAST.100", "REF","ALT","TRIMER","MT","Chrom", "MAF"),with=F])
summary(m.tree)

#010615
########Test with syn/non_syn information########-1000G
MAF.ANNOVAR<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/ANNOVAR/122114.THOUSAND.SNP.TRIMERS.PLUS.PHAST.avinput.exonic_variant_function", drop=c(1,6),
                   sep="\t", header=F, stringsAsFactors=F)
MAF.ANNOVAR$V3<-sapply(MAF.ANNOVAR$V3, function(x) strsplit(x,":")[[1]][1])
setnames(MAF.ANNOVAR,  c("TYPE", "Hugo_Symbol","Chrom", "Position", "REF","ALT","TRIMER","MT","EXON","AF", "PHAST.100"))
MAF.ANNOVAR<-MAF.ANNOVAR[TYPE %in% c("nonsynonymous SNV", "synonymous SNV"),] #Syn and non-syn only
MAF.ANNOVAR$MAF<-ifelse(MAF.ANNOVAR$AF>0.5, 1-MAF.ANNOVAR$AF, MAF.ANNOVAR$AF)
MAF.ANNOVAR<-MAF.ANNOVAR[!is.na(MT),]#clear of potential NXN trimers MT of NA!!!
MAF.ANNOVAR<-MAF.ANNOVAR[EXON==TRUE,]
ggplot(MAF.ANNOVAR, aes(PHAST.100, MAF, colour=TYPE)) + geom_point() + facet_wrap(~TYPE) + scale_y_log10()

annovar.mod<-lm(MAF~exp(PHAST.100)+TYPE+REF+ALT, MAF.ANNOVAR)
summary(annovar.mod)

MAF.ANNOVAR.FILTERED<-MAF.ANNOVAR[MAF>=0.000199682,]
MAF.ANNOVAR.FILTERED$PHAST.CUT<-cut(MAF.ANNOVAR.FILTERED$PHAST.100, seq(0,1,0.1),include.lowest=T)

ggplot(MAF.ANNOVAR.FILTERED, aes(PHAST.CUT, MAF, colour=TYPE)) + geom_boxplot() + facet_wrap(~TYPE) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme.format

#010715
#######Try with Phastcons 45####### NOTE: CORRECTED REF and ALT nts!!!!
THOUSAND.PHAST.45<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/ANNOVAR/010715.THOUSAND.SNP.TRIMERS.PLUS.PHAST.45.avinput.exonic_variant_function", 
                         header=F, sep="\t", stringsAsFactors=F, drop=c(1,6:8))
THOUSAND.PHAST.45$Hugo_Symbol<-sapply(THOUSAND.PHAST.45$V3, function(x) strsplit(x,":")[[1]][1])
THOUSAND.PHAST.45$REF<-sapply(THOUSAND.PHAST.45$V3, function(x) {
  t<-unlist(strsplit(unlist(strsplit(x,":c."))[2], ":p."))[1]
  t<-unlist(strsplit(t,""))
  return(t[1])
  })
THOUSAND.PHAST.45$ALT<-sapply(THOUSAND.PHAST.45$V3, function(x) {
  t<-unlist(strsplit(unlist(strsplit(x,c(":c.")))[2], ":p."))[1]
  t<-unlist(strsplit(t,""))
  return(t[length(t)])
})
THOUSAND.PHAST.45$V3<-NULL
setnames(THOUSAND.PHAST.45, c("TYPE", "Chrom","Position","TRIMER","MT", "EXON","AF","PHAST.45","Hugo_Symbol","REF","ALT"))
THOUSAND.PHAST.45<-THOUSAND.PHAST.45[TYPE %in% c("nonsynonymous SNV", "synonymous SNV"),]
THOUSAND.PHAST.45<-THOUSAND.PHAST.45[EXON==TRUE,] #Excluded for some analysis
THOUSAND.PHAST.45<-THOUSAND.PHAST.45[!is.na(MT),]
THOUSAND.PHAST.45$MAF<-ifelse(THOUSAND.PHAST.45$AF>0.5, 1-THOUSAND.PHAST.45$AF,THOUSAND.PHAST.45$AF)

THOUSAND.PHAST.45.FILTERED<-THOUSAND.PHAST.45[MAF>=0.000199682,] 
THOUSAND.PHAST.45.FILTERED$PHAST.CUT<-cut(THOUSAND.PHAST.45.FILTERED$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
ggplot(THOUSAND.PHAST.45.FILTERED, aes(PHAST.CUT, MAF, colour=TYPE)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme.format + facet_wrap(~TYPE)
ggplot(THOUSAND.PHAST.45.FILTERED, aes(PHAST.CUT, MAF, fill=ALT)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme.format + facet_wrap(~TYPE)
summary(lm(MAF~exp(PHAST.45)+TYPE+REF+ALT+Chrom,data=THOUSAND.PHAST.45.FILTERED))

table(THOUSAND.PHAST.45.FILTERED$PHAST.CUT)
THOUSAND.PHAST.45.FILTERED[PHAST.45==1,]
#Try including normal expression levels of all genes (need to find info)
length(unique(as.vector(THOUSAND.PHAST.45$Hugo_Symbol)))

######Inclusde replication time
CHEN.REP<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/010815.CHEN.REP.TIMES", header=T, sep="\t",stringsAsFactors=F, drop=2)
ggplot(CHEN.REP, aes(REP.CLASS, REP.TIME)) + geom_boxplot() + theme.format
table(CHEN.REP$REP.CLASS)

THOUSAND.PHAST.45<-merge(THOUSAND.PHAST.45, CHEN.REP, by="Hugo_Symbol")
THOUSAND.PHAST.45$REP.CLASS<-factor(THOUSAND.PHAST.45$REP.CLASS, levels=c("high","medium","low"))
summary(lm(MAF~exp(PHAST.45)+REP.TIME+TYPE+REF+ALT+Chrom, data=THOUSAND.PHAST.45 ))

ggplot(THOUSAND.PHAST.45[TYPE=="nonsynonymous SNV" & EXON==TRUE,], aes(REP.CLASS, MAF, colour=TYPE)) + geom_boxplot() + 
  facet_grid(REF~ALT) + theme.format + scale_y_log10()

THOUSAND.PHAST.45.FILTERED<-THOUSAND.PHAST.45#[MAF>=0.000199682,]
THOUSAND.PHAST.45.FILTERED$PHAST.CUT<-cut(THOUSAND.PHAST.45.FILTERED$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)

ggplot(THOUSAND.PHAST.45.FILTERED, aes(PHAST.CUT, MAF, colour=REP.CLASS)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme.format + facet_grid(~TYPE)
ggplot(THOUSAND.PHAST.45.FILTERED, aes(PHAST.CUT, MAF, colour=REP.CLASS)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme.format + facet_grid(TYPE~REF)
ggplot(THOUSAND.PHAST.45.FILTERED, aes(REP.CLASS, MAF, colour=TYPE)) + geom_boxplot() + theme.format + facet_wrap(~TYPE) +
  scale_y_log10()

summary(lm(MAF~exp(PHAST.45)+REP.CLASS+TYPE+REF:ALT+Chrom, data=THOUSAND.PHAST.45.FILTERED )) 

#MODIFY BAYES
#Given that |replication time(class) - Increase number of replication classes

#########011015########
#Learn how noise is represented in genes
THOUSAND.PHAST.45
ggplot(THOUSAND.PHAST.45[Hugo_Symbol=="TP53",], aes(Position, MAF, colour=REF)) + geom_point() +
  theme.format + scale_y_log10() + facet_grid(REF~TYPE)
ggplot(THOUSAND.PHAST.45[Hugo_Symbol=="BRCA1",], aes(Position, MAF, colour=TYPE)) + geom_point() +
  theme.format + scale_y_log10() + facet_grid(REF~ALT)
ggplot(THOUSAND.PHAST.45[Hugo_Symbol=="TTN",], aes(Position, MAF, colour=TYPE)) + geom_point() +
  theme.format + scale_y_log10() + facet_grid(REF~ALT)

ggplot(THOUSAND.PHAST.45.FILTERED[EXON==TRUE,], aes(PHAST.CUT, MAF,colour=TYPE)) + geom_boxplot() + theme.format +
  facet_grid(REF~ALT) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()

THOUSAND.PHAST.45.FILTERED[EXON==TRUE,][,list(low.count=sum(MAF<0.000199682),
                                              high.count=sum(MAF>0.000199682)), by=PHAST.CUT][order(PHAST.CUT),]
table(THOUSAND.PHAST.45.FILTERED[EXON==TRUE,]$PHAST.CUT)

#Formal function based on Pr(Cancer | REF->ALT n REP.TIME n PHAST.45)
Function.THOUSAND.Prob<-function(THOUSAND.PHAST.CHEN.TABLE){
  
  require(data.table)
  
  #Filter table
  main.table<-THOUSAND.PHAST.CHEN.TABLE[]
  
}