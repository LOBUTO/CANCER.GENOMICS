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

#Try including normal expression levels of all genes (need to find info)
length(unique(as.vector(THOUSAND.PHAST.45$Hugo_Symbol)))

######Inclusde replication time
CHEN.REP<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/010815.CHEN.REP.TIMES.MID.4.EXP", header=T, sep="\t",stringsAsFactors=F, drop=2)
setkey(CHEN.REP)
ggplot(CHEN.REP, aes(REP.CLASS, REP.TIME)) + geom_boxplot() + theme.format
table(CHEN.REP$REP.CLASS)

THOUSAND.PHAST.45.MID.1<-merge(THOUSAND.PHAST.45, CHEN.REP, by="Hugo_Symbol")
setkey(THOUSAND.PHAST.45.MID.1)
THOUSAND.PHAST.45.MID.1<-unique(THOUSAND.PHAST.45.MID.1)

THOUSAND.PHAST.45.FILTERED<-THOUSAND.PHAST.45[MAF>=0.000199682,]
THOUSAND.PHAST.45.FILTERED$PHAST.CUT<-cut(THOUSAND.PHAST.45.FILTERED$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
THOUSAND.PHAST.45.MID.1.FILTERED<-merge(THOUSAND.PHAST.45.FILTERED, CHEN.REP, by="Hugo_Symbol")
setkey(THOUSAND.PHAST.45.MID.1.FILTERED)
THOUSAND.PHAST.45.MID.1.FILTERED<-unique(THOUSAND.PHAST.45.MID.1.FILTERED)

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

#Formal function based on Pr(Cancer | REF->ALT n REP.TIME n PHAST.45)
Function.THOUSAND.Prob<-function(THOUSAND.PHAST.CHEN.TABLE, CUTS=F, PHAST=T){
  #Keep in mind that pre-processing needs to be done (mt, merging with chen, setnames...etc)
  #FIX!!!!
  require(data.table)
  
  #Filter table
  main.table<-THOUSAND.PHAST.CHEN.TABLE[EXON==TRUE,]
  
  #Prep classes
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  
  if (CUTS==T){
    main.table$PHAST.CLASS<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)  
    main.table$REP.CLASS<-cut(main.table$REP.TIME, seq(0,1,0.05),include.lowest=T)
  }
  
  #Calculate probabilities - APPLY CHANGES TO TCGA FUNCTION!!!!!
  n.experiments<-sum(main.table$MAF)
  main.table[,REF.ALT.PROB:=sum(MAF)/n.experiments,by=REF.ALT]
  main.table[,TYPE.PROB:=sum(MAF)/n.experiments,by=TYPE]
  
  if (CUTS==T){
    main.table[,REP.TIME.PROB:=sum(MAF)/n.experiments,by=REP.CLASS]
    main.table[,PHAST.45.PROB:=sum(MAF)/n.experiments,by=PHAST.CLASS]
  } else{
    main.table[,REP.TIME.PROB:=sum(MAF)/n.experiments,by=REP.TIME]
    main.table[,PHAST.45.PROB:=sum(MAF)/n.experiments,by=PHAST.45]
  }
  
  #Set key to all columns so we don't remove non-duplicates with unique later on
  setkey(main.table)
  
  if (PHAST==T){
    main.table$THOUSAND.PROB<-main.table$REF.ALT.PROB*main.table$REP.TIME.PROB*main.table$TYPE.PROB*main.table$PHAST.45.PROB
    
    #Clean up and return
    if (CUTS==T){
      prop.table<-unique(main.table[,c("REF.ALT","REP.CLASS", "TYPE", "PHAST.CLASS", "THOUSAND.PROB"),with=F])
    } else{
      prop.table<-unique(main.table[,c("REF.ALT","REP.TIME", "TYPE", "PHAST.45", "THOUSAND.PROB"),with=F]) 
    }
  } else {
    main.table$THOUSAND.PROB<-main.table$REF.ALT.PROB*main.table$REP.TIME.PROB*main.table$TYPE.PROB
    
    #Clean up and return
    if (CUTS==T){
      prop.table<-unique(main.table[,c("REF.ALT","REP.CLASS", "TYPE", "THOUSAND.PROB"),with=F])
    } else{
      prop.table<-unique(main.table[,c("REF.ALT","REP.TIME", "TYPE", "THOUSAND.PROB"),with=F]) 
    }
  }
  return(prop.table)
}

Function.TCGA.Prob<-function(TCGA.PHAST.CHEN.TABLE, CUTS=F, PHAST=T){
  #Keep in mind that pre-processing needs to be done (mt, mergint with chen, setnames....etc)
  
  require(data.table)
  
  #Filter table - change labels for compatibility with thousand table
  main.table<-TCGA.PHAST.CHEN.TABLE[TYPE %in% c("Missense_Mutation", "Silent"),]
  main.table$TYPE<-ifelse(main.table$TYPE=="Missense_Mutation", "nonsynonymous SNV", "synonymous SNV")
  
  #Prep class
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  if (CUTS==T){
    main.table$PHAST.CLASS<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)  
    main.table$REP.CLASS<-cut(main.table$REP.TIME, seq(0,1,0.05),include.lowest=T)
  }
  
  #Calculate probabilities and add to each line
  n.experiments<-sum(main.table$MUT.FREQ)
  main.table[,REF.ALT.PROB:=sum(MUT.FREQ)/n.experiments,by=REF.ALT]
  main.table[,TYPE.PROB:=sum(MUT.FREQ)/n.experiments,by=TYPE]
  
  if (CUTS==T){
    main.table[,REP.TIME.PROB:=sum(MUT.FREQ)/n.experiments,by=REP.CLASS]
    main.table[,PHAST.45.PROB:=sum(MUT.FREQ)/n.experiments,by=PHAST.CLASS]
  } else{
    main.table[,REP.TIME.PROB:=sum(MUT.FREQ)/n.experiments,by=REP.TIME]
    main.table[,PHAST.45.PROB:=sum(MUT.FREQ)/n.experiments,by=PHAST.45]
  }
  
  #Use phastcons scores?
  if (PHAST==T){
    main.table$TCGA.PROB<-main.table$REF.ALT.PROB*main.table$REP.TIME.PROB*main.table$TYPE.PROB*main.table$PHAST.45.PROB  
  } else {
    main.table$PHAST.45.PROB<-NULL
    main.table$TCGA.PROB<-main.table$REF.ALT.PROB*main.table$REP.TIME.PROB*main.table$TYPE.PROB  
  }
  
  #Return
  return (main.table)
}

Function.Main.Bayes<-function(thousand.prop.table, tcga.main.table, cancer.prob=0.12, threshold=F, all=T, CUTS=F, PHAST=T){
  #Merges cancer and non-cancer probabilities to calculate non-naive bayesian probability per site
  
  library(data.table)
  
  #Merge probabilites (TCGA.PROB, THOUSAND.PROB)
  if (all==T){
    threshold.prob<-min(thousand.prop.table$THOUSAND.PROB)
    
    #Analyse by cuts?
    if (CUTS==T){
      if(PHAST==T){
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.CLASS","TYPE","PHAST.CLASS"), all.x=TRUE)  
      } else{
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.CLASS","TYPE"), all.x=TRUE)  
      }
      
    } else{
      if (PHAST==T){
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.TIME","TYPE","PHAST.45"), all.x=TRUE)    
      } else{
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.TIME","TYPE"), all.x=TRUE)    
      }
    }
    
    #Apply lowest threshold?
    if (threshold==T){
      main.table[is.na(main.table),]<-threshold.prob
    }  
  } else {
    
    #Analyze by cuts?
    if (CUTS==T){
      if (PHAST==T){
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.CLASS","TYPE","PHAST.CLASS"))  
      } else{
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.CLASS","TYPE"))  
      } 
      
    } else{
      if (PHAST==T){
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.TIME","TYPE","PHAST.45"))  
      } else{
        main.table<-merge(tcga.main.table, thousand.prop.table, by=c("REF.ALT","REP.TIME","TYPE"))  
      }
      
    }
  }
  
  #Calculate bayes prob per site based on features
  main.table$BAYES.PROB<-(main.table$TCGA.PROB*cancer.prob)/(main.table$TCGA.PROB*cancer.prob + main.table$THOUSAND.PROB*(1-cancer.prob))
  
  #Clean up and return
  main.table<-main.table[order(BAYES.PROB, decreasing=T),]
  return(main.table)
}

Function.Bayes.Test<-function(main.table, cosmic.genes){
  
  require (data.table)
  
  #Classify
  main.table$cancer<-main.table$Hugo_Symbol %in% cosmic.genes
  
  #Analyze...
  
  #Return
  return (main.table) 
}

########Test bayes########
BRCA.COSMIC.GENES<-c("AKT1", "BAP1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "EP300", "ERBB2",
                     "FOXA1", "MAP2K4","PALB2","PBRM1", "PIK3CA", "RB1","TP53")

#Load TCGA file
TCGA.MUT.45<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/011015.BRCA.MAF.TRIMERS.45.csv", header=F, sep="\t",
                   stringsAsFactors=F, drop=c(8:9,11))
setnames(TCGA.MUT.45, c("Chrom","Position","Hugo_Symbol", "MUT.FREQ", "TYPE", "REF","ALT","PHAST.45"))
TCGA.MUT.45.MID.4<-merge(TCGA.MUT.45, CHEN.REP, by="Hugo_Symbol")
setkey(TCGA.MUT.45.MID.1)
TCGA.MUT.45.MID.1<-unique(TCGA.MUT.45.MID.1)

#Apply functions
thousand.prop.table<-Function.THOUSAND.Prob(THOUSAND.PHAST.45.MID.4, CUTS=F, PHAST=T)
tcga.main.table<-Function.TCGA.Prob(TCGA.MUT.45.MID.4, CUTS=F, PHAST=T)
tcga.bayes<-Function.Main.Bayes(thousand.prop.table,tcga.main.table,0.12,threshold=F,all=F,CUTS=F,PHAST=T)
tcga.bayes.test<-Function.Bayes.Test(tcga.bayes,BRCA.COSMIC.GENES)

ggplot(tcga.bayes.test, aes(cancer, BAYES.PROB, colour=cancer)) + geom_boxplot() + geom_jitter(size=1) + theme.format + facet_wrap(~TYPE)

#DON'T USE UNIQUE()!!! IT IGNORES POSITION, UNLESS U DO SETKEY() BEFOREHAND!!!
#TRY ALTERNATIVE TO POSITIVE STRAND!!!! SAME ARE NOT ANNOTATED WELL 
#COULD ALSO DO PROBABILITIES PER GENE!!! (by=c(""..,"Hugo_Symbol"))
#COULD ALSO USE EXPRESSION LEVELS AS PREDICTORS????
bayes.box.test<-data.table()
for (prob in c(0.95, 0.9, 0.85, 0.80, 0.75, 0.7, 0.6, 0.5,0.4,0.3,0.2)){
  n.cancer<-length(unique(as.vector(tcga.bayes.test[BAYES.PROB>=prob & cancer==TRUE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  n.non.cancer<-length(unique(as.vector(tcga.bayes.test[BAYES.PROB>=prob & cancer==FALSE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  bayes.box.test<-rbind(bayes.box.test, data.table(PROB=prob, CANCER=n.cancer, NON.CANCER=n.non.cancer))
}
ggplot(melt(bayes.box.test,id.vars="PROB"), aes(PROB, value, colour=variable)) + geom_bar(stat="identity", position="dodge") + theme.format +
  geom_text(aes(label=value), position=position_dodge(width=0.04), vjust=-0.25) 
 
tcga.bayes[Hugo_Symbol=="TP53",]

#########011215########
######Modify functions  so that they use joint probabilities rather than conditional independence per gene

#Functions
Function.THOUSAND.Prob.Joint<-function(THOUSAND.PHAST.CHEN.TABLE, exp.table,hic.table, biogrid, noise=1, rounded=3){
  
  require(data.table)
  require(car)
  
  #Filter table
  main.table<-THOUSAND.PHAST.CHEN.TABLE[EXON==TRUE,]
  main.table$REP.TIME<-round(main.table$REP.TIME, rounded)
  
  #Prep classes - Recode ref.alt to reflect binomial chance on either strand
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  
  #Introduce expression info
  main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  main.table$logFC<-abs(main.table$logFC)
  main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce PPI info
  main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  
  #Introduce chromatin open state info
  main.table<-merge(main.table, hic.table[,c("Hugo_Symbol","HIC.CUT"),with=F],by="Hugo_Symbol")
  
  #Modified joint
  #main.table$PHAST.CUT<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
  main.table[,BACK.PROB:=sum(MAF)*noise, by=c("TYPE", "REF.ALT","REP.CLASS","HIC.CUT","EXP.CUTS", "DEGREE.CUT")]
  main.table[,THOUSAND.PROB:=sum(MAF)/BACK.PROB, by=c("Hugo_Symbol", "TYPE","REF.ALT","REP.CLASS","HIC.CUT","EXP.CUTS","DEGREE.CUT")]
  
  #Get joint probabilites assuming conditional independence
#   n.maf<-sum(main.table$MAF)/noise
#   
#   main.table[,HUGO.PROB:=sum(MAF)/n.maf, by="Hugo_Symbol"]
#   main.table[,TYPE.PROB:=sum(MAF)/n.maf, by="TYPE"]
#   main.table[,REF.ALT.PROB:=sum(MAF)/n.maf, by="REF.ALT"]
#   main.table[,REP.PROB:=sum(MAF)/n.maf, by="REP.TIME"]
#   main.table[,PHAST.PROB:=sum(MAF)/n.maf, by="PHAST.45"]
#   main.table[,EXP.PROB:=sum(MAF)/n.maf, by="EXP.CUTS"]
#   main.table[,REP.CLASS.PROB:=sum(MAF)/n.maf, by="REP.CLASS"]
#   
#   #Calculate probability 
#   main.table<-main.table[,list(THOUSAND.PROB=TYPE.PROB*REF.ALT.PROB*REP.CLASS.PROB*EXP.PROB*HUGO.PROB*PHAST.PROB),by=c("Hugo_Symbol", "REF.ALT","TYPE","REP.CLASS", "PHAST.45", "EXP.CUTS")]
  
  #Clean up and return
  main.table<-main.table[,c("Hugo_Symbol","REF.ALT","TYPE","REP.CLASS","HIC.CUT","EXP.CUTS", "DEGREE.CUT", "THOUSAND.PROB"),with=F]
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.TCGA.Prob.Joint<-function(TCGA.PHAST.CHEN.TABLE, exp.table, hic.table,biogrid, rounded=3){
  require(data.table)
  
  #Filter table - change labels for compatibility with thousand table
  main.table<-TCGA.PHAST.CHEN.TABLE[TYPE %in% c("Missense_Mutation", "Silent"),]
  main.table$TYPE<-ifelse(main.table$TYPE=="Missense_Mutation", "nonsynonymous SNV", "synonymous SNV")
  
  #Prep class
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  main.table$REP.TIME<-round(main.table$REP.TIME, rounded)
  
  #Introduce expression info
  main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  main.table$logFC<-abs(main.table$logFC)
  main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce chromatin open state info
  main.table<-merge(main.table, hic.table[,c("Hugo_Symbol","HIC.CUT"),with=F],by="Hugo_Symbol")
  
  #Introduce PPI info
  main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  
  #Modified joint
  #main.table$PHAST.CUT<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
  main.table[,BACK.PROB:=sum(MUT.FREQ), by=c("TYPE", "REF.ALT","REP.CLASS","HIC.CUT","EXP.CUTS", "DEGREE.CUT")]
  main.table[,TCGA.PROB:=sum(MUT.FREQ)/BACK.PROB, by=c("Hugo_Symbol", "TYPE","REF.ALT","REP.CLASS","HIC.CUT","EXP.CUTS", "DEGREE.CUT")]
  
#   #Get joint probabilites assuming conditional independence
#   n.maf<-sum(main.table$MUT.FREQ)
#   main.table[,HUGO.PROB:=sum(MUT.FREQ)/n.maf, by="Hugo_Symbol"]
#   main.table[,TYPE.PROB:=sum(MUT.FREQ)/n.maf, by="TYPE"]
#   main.table[,REF.ALT.PROB:=sum(MUT.FREQ)/n.maf, by="REF.ALT"]
#   main.table[,REP.PROB:=sum(MUT.FREQ)/n.maf, by="REP.TIME"]
#   main.table[,PHAST.PROB:=sum(MUT.FREQ)/n.maf, by="PHAST.45"]
#   main.table[,EXP.PROB:=sum(MUT.FREQ)/n.maf, by="EXP.CUTS"]
#   main.table[,REP.CLASS.PROB:=sum(MUT.FREQ)/n.maf, by="REP.CLASS"]
#   
#   #Calculate probability
#   main.table[,TCGA.PROB:=TYPE.PROB*REF.ALT.PROB*REP.CLASS.PROB*EXP.PROB*HUGO.PROB*PHAST.PROB ,by=c("Hugo_Symbol", "REF.ALT","TYPE","REP.CLASS","PHAST.45", "EXP.CUTS")]
  
  #Clean up and return
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.Main.Bayes.Joint<-function(thousand.prop.joint, tcga.prob.joint, cancer.prob=0.12, threshold=F, all=T){
  #Merges cancer and non-cancer probabilities to calculate non-naive bayesian probability per site
  
  library(data.table)
  
  #Merge probabilites (TCGA.PROB, THOUSAND.PROB)
  main.table<-merge(tcga.prob.joint, thousand.prop.joint, by=c("REF.ALT","TYPE","Hugo_Symbol","REP.CLASS","HIC.CUT","EXP.CUTS"))
  
  #Calculate bayes prob per site based on features
  main.table$BAYES.PROB<-(main.table$TCGA.PROB*cancer.prob)/(main.table$TCGA.PROB*cancer.prob + main.table$THOUSAND.PROB*(1-cancer.prob))
  
  #Clean up and return
  main.table<-main.table[order(BAYES.PROB, decreasing=T),]
  return(main.table)
}

#Apply
thousand.test<-Function.THOUSAND.Prob.Joint(THOUSAND.PHAST.45.MID.4[MAF!=0,], brca.exp,hic,biogrid.degree,noise=1,rounded=2)
tcga.test<-Function.TCGA.Prob.Joint(TCGA.MUT.45.MID.4, brca.exp,hic,biogrid.degree,rounded=2)
test.bayes<-Function.Main.Bayes.Joint(thousand.test, tcga.test,cancer.prob=0.08)
test.bayes.plot<-Function.Bayes.Test(test.bayes, BRCA.COSMIC.GENES)

ggplot(test.bayes.plot, aes(cancer, BAYES.PROB, colour=cancer)) + geom_boxplot() + geom_jitter(size=0.5) + theme.format + facet_wrap(~TYPE)

ggplot(test.bayes.plot, aes(REP.TIME, BAYES.PROB, colour=cancer)) + geom_point()  + theme.format + facet_grid(cancer~TYPE) + 
  geom_text(data=subset(test.bayes.plot, cancer==TRUE), aes(REP.TIME,BAYES.PROB,label=Hugo_Symbol), size=4)
ggplot(test.bayes.plot, aes(PHAST.45, BAYES.PROB, colour=cancer, label=Hugo_Symbol)) + geom_point()  + theme.format + facet_grid(cancer~TYPE) + 
  geom_text(data=subset(test.bayes.plot, cancer==TRUE), aes(REP.TIME,BAYES.PROB,label=Hugo_Symbol), size=3)

bayes.box.test<-data.table()
for (prob in c(0.99,0.95, 0.9, 0.85, 0.80, 0.75, 0.7, 0.6, 0.5,0.4,0.3,0.2,0)){
  n.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & cancer==TRUE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  n.non.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & cancer==FALSE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  bayes.box.test<-rbind(bayes.box.test, data.table(PROB=prob, CANCER=n.cancer, NON.CANCER=n.non.cancer))
  print (c(prob, n.cancer, n.non.cancer))
}
ggplot(melt(bayes.box.test,id.vars="PROB"), aes(PROB, value, colour=variable)) + geom_bar(stat="identity", position="dodge") + theme.format +
  geom_text(aes(label=value), position=position_dodge(width=0.04), vjust=-0.25) 

tcga.test[Hugo_Symbol=="A1CF",]
thousand.test[Hugo_Symbol=="A1CF",]
test.bayes[Hugo_Symbol=="BRCA2",]
unique(test.bayes.plot[cancer==TRUE & BAYES.PROB>=0.85 & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)
test.bayes.plot

ggplot(THOUSAND.PHAST.45.MID.2[MAF!=0,], aes(REP.CLASS, MAF, colour=TYPE)) + geom_boxplot() + geom_jitter(size=2) + theme.format + facet_wrap(~TYPE) + scale_y_log10()
ggplot(TCGA.MUT.45.MID.2, aes(REP.CLASS, MUT.FREQ, colour=TYPE)) + geom_boxplot() + geom_jitter(size=2) + theme.format + facet_wrap(~TYPE) + scale_y_log10()

###Try with expression values
brca.exp<-readRDS("LOGS/111114.BRCA.DIFF.EXP.rds")
setnames(brca.exp, c("logFC","AveExpr","t","P.Value","adj.P.Val","B","Hugo_Symbol"))
brca.exp[Hugo_Symbol %in% BRCA.COSMIC.GENES,]
hist(brca.exp$logFC)

###Try with chromatin open close states
hic<-fread("http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/gene.covariates.txt",header=T, sep="\t",
             stringsAsFactors=F)
setnames(hic, c("Hugo_Symbol", "EXPR","REPLICATION","hic"))
hic<-hic[!is.na(hic),]
hic$HIC.CUT<-cut(hic$hic, quantile(hic$hic, c(0,0.33,0.67,1)), labels=c("low","medium","high"))

##############011715##############
#Modified functions
#Functions
Function.THOUSAND.Prob.Joint.2<-function(THOUSAND.PHAST.CHEN.TABLE, exp.table,hic.table, biogrid, noise=1, rounded=3){
  
  require(data.table)
  require(car)
  
  #Filter table
  main.table<-THOUSAND.PHAST.CHEN.TABLE[EXON==TRUE,]
  main.table$REP.TIME<-round(main.table$REP.TIME, rounded)
  
  #Prep classes - Recode ref.alt to reflect binomial chance on either strand
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  
  #Introduce expression info
  main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  main.table$logFC<-abs(main.table$logFC)
  main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce PPI info
  #main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  
  #Introduce chromatin open state info
  #main.table<-merge(main.table, hic.table[,c("Hugo_Symbol","HIC.CUT"),with=F],by="Hugo_Symbol")
  
  #Modified joint
  n.maf<-sum(main.table$MAF)
  main.table<-main.table[,list(THOUSAND.PROB=sum(MAF)/n.maf), by=c("REF.ALT","TYPE","Chrom","REP.TIME","EXP.CUTS")]
  
  #Clean up and return
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.TCGA.Prob.Joint.2<-function(TCGA.PHAST.CHEN.TABLE, exp.table, hic.table,biogrid, rounded=3){
  require(data.table)
  
  #Filter table - change labels for compatibility with thousand table
  main.table<-TCGA.PHAST.CHEN.TABLE[TYPE %in% c("Missense_Mutation", "Silent"),]
  main.table$TYPE<-ifelse(main.table$TYPE=="Missense_Mutation", "nonsynonymous SNV", "synonymous SNV")
  
  #Prep class
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  main.table$REP.TIME<-round(main.table$REP.TIME, rounded)
  
  #Introduce expression info
  main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  main.table$logFC<-abs(main.table$logFC)
  main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce chromatin open state info
  #main.table<-merge(main.table, hic.table[,c("Hugo_Symbol","HIC.CUT"),with=F],by="Hugo_Symbol")
  
  #Introduce PPI info
  #main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  
  #Modified joint
  #n.maf<-sum(main.table$MUT.FREQ)
  #main.table[,TCGA.PROB:=sum(MUT.FREQ)/n.maf, by=c("TYPE","REF.ALT","Chrom","REP.TIME","EXP.CUTS")]
  
  #Clean up and return
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.Main.Bayes.Joint.2<-function(thousand.prop.joint, tcga.prob.joint, cancer.prob=0.12, threshold=F, all=T){
  #Merges cancer and non-cancer probabilities to calculate non-naive bayesian probability per site
  
  library(data.table)
  
  #Merge probabilites (TCGA.PROB, THOUSAND.PROB)
  main.table<-merge(tcga.prob.joint, thousand.prop.joint, by=c("TYPE","REF.ALT","Chrom","REP.TIME","EXP.CUTS"))
  
  #Calculate bayes prob per site based on features
  main.table$BAYES.PROB<-(main.table$MUT.FREQ*cancer.prob)/(main.table$MUT.FREQ*cancer.prob + main.table$THOUSAND.PROB*(1-cancer.prob))
  
  #Clean up and return
  main.table<-main.table[order(BAYES.PROB, decreasing=T),]
  return(main.table)
}

#Apply
thousand.test<-Function.THOUSAND.Prob.Joint.2(THOUSAND.PHAST.45.MID.2[MAF!=0,], brca.exp,hic,biogrid.degree,noise=1,rounded=2)
tcga.test<-Function.TCGA.Prob.Joint.2(TCGA.MUT.45.MID.2, brca.exp,hic,biogrid.degree,rounded=2)
test.bayes<-Function.Main.Bayes.Joint.2(thousand.test, tcga.test,cancer.prob=0.12)
test.bayes.plot<-Function.Bayes.Test(test.bayes, BRCA.COSMIC.GENES)

ggplot(test.bayes.plot, aes(cancer, BAYES.PROB, colour=cancer)) + geom_boxplot() + geom_jitter(size=0.5) + theme.format + facet_wrap(~TYPE)

bayes.box.test<-data.table()
for (prob in c(0.95, 0.9, 0.85, 0.80, 0.75, 0.7, 0.6, 0.5,0.4,0.3,0.2,0)){
  n.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & cancer==TRUE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  n.non.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & cancer==FALSE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
  bayes.box.test<-rbind(bayes.box.test, data.table(PROB=prob, CANCER=n.cancer, NON.CANCER=n.non.cancer))
  print (c(prob, n.cancer, n.non.cancer))
}
ggplot(melt(bayes.box.test,id.vars="PROB"), aes(PROB, value, colour=variable)) + geom_bar(stat="identity", position="dodge") + theme.format +
  geom_text(aes(label=value), position=position_dodge(width=0.04), vjust=-0.25) 

test.bayes[Hugo_Symbol=="TTN",]

TCGA.MUT.45.MID.2
Function.BAYES.SELF<-function(main.table, cancer.prob=0.12, rounded=3, cancer.genes){
  require(data.table) 
  require(car)
  
  #Prep class
  main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
  main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')
  main.table$REP.TIME<-round(main.table$REP.TIME, rounded)
  
  #Calculate background probabilities based on synonymous mutations
  n.sil<-sum(main.table[TYPE=="Silent",]$MUT.FREQ)
  back.table<-main.table[TYPE=="Silent",][,list(BACK.PROB=sum(MUT.FREQ)/n.sil), by="REF.ALT"]
  
  #Integrate to calculate bayesian probabilities
  mis.table<-main.table[TYPE=="Missense_Mutation",]
  mis.table<-merge(mis.table, back.table, by="REF.ALT")
  
  #Calculate bayesian prob
  mis.table$BAYES.PROB<-mis.table$MUT.FREQ*cancer.prob/((mis.table$MUT.FREQ*cancer.prob)+mis.table$BACK.PROB*(1-cancer.prob))
  
  #Clean up and Return
  mis.table$cancer<-mis.table$Hugo_Symbol %in% cancer.genes
  mis.table<-mis.table[order(BAYES.PROB, decreasing=T),]
  return(mis.table)
  
}

test.self<-Function.BAYES.SELF(TCGA.MUT.45.MID.2, 0.12, cancer.genes=BRCA.COSMIC.GENES)

ggplot(test.self, aes(cancer, BAYES.PROB, colour=cancer)) + geom_boxplot() + geom_jitter(size=2) + theme.format