######COPY NUMBER VARATION FUNCTION
#080614

#Process CNV Thresholded file to Table 1 type
BRCA.CNV.THRESHOLDED<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/CNV/BRCA/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2014041600.0.0/all_thresholded.by_genes.txt",header=T, sep="\t",stringsAsFactors=F))
dim(BRCA.CNV.THRESHOLDED)
BRCA.CNV.THRESHOLDED[,1:7,with=F]
BRCA.CNV.TABLE<-BRCA.CNV.THRESHOLDED[,c(1,4:ncol(BRCA.CNV.THRESHOLDED)),with=F]
BRCA.CNV.TABLE<-melt(BRCA.CNV.TABLE, id.vars="Gene.Symbol",variable.name="Sample")
BRCA.CNV.TABLE<-BRCA.CNV.TABLE[abs(value)==2,]
BRCA.CNV.TABLE$PATIENT<-substr(BRCA.CNV.TABLE$Sample,1,16)
BRCA.CNV.TABLE$Sample<-NULL
BRCA.CNV.TABLE$value<-NULL

length(unique(BRCA.CNV.TABLE$Hugo_Symbol))

CNV.1<-BRCA.CNV.TABLE[,list(PATIENTS=length(Sample)), by="Gene.Symbol"]
CNV.2<-BRCA.CNV.TABLE[,list(GENES=length(Gene.Symbol)), by="Sample"]
ggplot(CNV.1, aes(PATIENTS)) + geom_histogram() + theme.format + ylab("Number of CNV Genes in Patient")
ggplot(CNV.2, aes(GENES)) + geom_histogram() + theme.format + ylab("Number of Patients with Number of CNVs")
length(unique(CNV.2$Sample))
CNV.1[order(PATIENTS, decreasing=T),]
CNV.1[Gene.Symbol=="ERLIN2",]

###
Table.v.BRCA.p.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.Table.v.p.rds")
Table.v.BRCA.p.2$v.PROTEIN.NORM<-Table.v.BRCA.p.2$v.PROTEIN * Table.v.BRCA.p.2$SAMPLE.COVERAGE
Table.v.BRCA.p.2[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

Table.BRCA.v.rank.enrich.p.2<-Function.p.rank.enrichment(Table.v.BRCA.p.2, c("Hugo_Symbol", "v.PROTEIN"), as.vector(COSMIC.BRCA$Symbol), normalize=T)
Function.j.rank.enrich.plot(Table.BRCA.v.rank.enrich.p.2)
Table.v.BRCA.p.2[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
Table.BRCA.v.rank.enrich.p.2$NON.CUM.TABLE[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
ggplot(Table.BRCA.v.rank.enrich.p.2$CUM.TABLE, aes(CUM.RECALL,CUM.PRECISION)) + geom_line()

length(unique(as.vector(dummy.1.cnv[Hugo_Symbol %in% as.vector(Table.v.BRCA.p.2[order(v.PROTEIN, decreasing=T),]$Hugo_Symbol)[1:100], ]$Tumor_Sample_Barcode)))

###combining
HYPER.TABLE.2.CALC
Table.v.BRCA.p.2
test<-as.data.table(merge(as.data.frame(HYPER.TABLE.2.CALC[,c(1,5),with=F]), as.data.frame(Table.v.BRCA.p.2), by="Hugo_Symbol",all=T,))
test<-test[P.VAL.ADJ<0.05 | is.na(P.VAL.ADJ),]

test$P.VAL.ADJ[is.na(test$P.VAL.ADJ)]<-1
test$v.PROTEIN[is.na(test$v.PROTEIN)]<-0

test$P.VAL.ADJ.NORM<-normalize.vector(-log(test$P.VAL.ADJ))
test<-test[P.VAL.ADJ<0.05,]
test<-test[P.VAL.ADJ<0.05 | v.PROTEIN>0.5,]
test[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
test$FACTOR<-test$v.PROTEIN*(1-0.8) + test$P.VAL.ADJ.NORM*0.8
Factor.enrichment<-Function.p.rank.enrichment(test, c("Hugo_Symbol", "v.PROTEIN"), as.vector(COSMIC.BRCA$Symbol), normalize=F)
Function.j.rank.enrich.plot(Factor.enrichment)
ggplot(Factor.enrichment$CUM.TABLE, aes(CUM.RECALL, CUM.PRECISION)) + geom_line()

dummy.1.cnv<-dummy.1$table.1[,1:2,with=F]
dummy.1.cnv$Tumor_Sample_Barcode<-sapply(as.character(dummy.1.cnv$Tumor_Sample_Barcode), function(x) paste0(strsplit(x,"-")[[1]][1:4],collapse=".")  )
dummy.1.cnv<-unique(rbind(dummy.1.cnv, BRCA.CNV.TABLE))

length(unique(as.vector(dummy.1.cnv[Hugo_Symbol %in% as.vector(test[order(v.PROTEIN, decreasing=T),]$Hugo_Symbol)[1:100], ]$Tumor_Sample_Barcode)))

######Function for PR for Patient Coverage#######

Function.PR.PATIENT.COVERAGE<-function(table.p, interest.columns, CANCER, c.table.1 ){
  
  require(data.table)
  require(parallel)
  
  #Process scoring table
  table.p<-table.p[,interest.columns, with=F] 
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  
  #CALCULATIONS
  background.genes<-unique(as.vector(table.p$Hugo_Symbol))
  background.patients<-unique(as.vector(c.table.1$Tumor_Sample_Barcode)) #Patients that we can choose from
  target.patients<-unique(as.vector(c.table.1[TYPE==CANCER & Hugo_Symbol %in% background.genes,]$Tumor_Sample_Barcode)) #Targets we could have chosen (White in urn) - ALL POSITIVES
  
  #Breaks for tables 0-1
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.0001)
  PR.BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  c.table.1<-copy(c.table.1)
  clusterExport(cl, varlist=c("table.p", "target.patients", "c.table.1"),envir=environment())
  
  #PR
  SCORE.RECALL<-parSapply(cl, PR.BREAKS, function (x)
    length(unique(as.vector(c.table.1[c.table.1$TYPE==CANCER & c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) / length(target.patients))
  
  SCORE.PRECISION<-parSapply(cl, PR.BREAKS, function(x)
    length(unique(as.vector(c.table.1[c.table.1$TYPE==CANCER & c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) / 
                            length(unique(as.vector(c.table.1[c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) )
  
  #Stop parallelization
  stopCluster(cl)
  
  #Clean up and Return
  dummy.return<-as.data.table(do.call(cbind, list(PR.BREAKS, SCORE.RECALL, SCORE.PRECISION)))
  setnames(dummy.return, colnames(dummy.return), c("PR.BREAKS", "RECALL", "PRECISION"))
  return(dummy.return)
  
}

Function.composite.table.1<-function(tables.1, cnv.tables, cancers) {
  
  require(base)
  require(data.table)
  
  dummy.list<-list()
  
  for (type in 1:length(cancers)){
    cancer.cnv<-readRDS(cnv.tables[type])
    setnames(cancer.cnv, colnames(cancer.cnv), c("Tumor_Sample_Barcode", "Hugo_Symbol"))
    
    cancer.table.1<-readRDS(tables.1[type])
    cancer.table.1<-cancer.table.1$table.1[,1:2,with=F]
    cancer.table.1$Tumor_Sample_Barcode<-sapply(as.character(cancer.table.1$Tumor_Sample_Barcode), function(x) paste0(strsplit(x,"-")[[1]][1:4],collapse=".")  )
    
    cancer.table.1<-unique(rbind(cancer.table.1, cancer.cnv))
    cancer.table.1$TYPE<-cancers[type]
    
    dummy.list[[type]]<-cancer.table.1
  }
  
  dummy.return<-as.data.table(do.call(rbind, dummy.list))
  return(dummy.return)
}

composite.tables.1<-Function.composite.table.1(c("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds",
                                                  "PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/072214_Table1_COAD.rds",
                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/072314_Table1_UCEC.rds",
                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/072214_Table1_GBM.rds"),
                                               c("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.GISTIC.TH.2.rds",
                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/080714.COAD.GISTIC.TH.2.rds",
                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/080714.GBM.GISTIC.TH.2.rds",
                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/080714.UCEC.GISTIC.TH.2.rds"),
                                               c("BRCA", "COAD","UCEC","GBM"))

test.FACTOR.PATIENT.ENRICHMENT<-Function.PR.PATIENT.COVERAGE(test, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
test.FACTOR.PATIENT.ENRICHMENT
ggplot(test.FACTOR.PATIENT.ENRICHMENT, aes(RECALL, PRECISION)) + geom_line() + geom_point()+theme.format

MUTSIG.PATIENT.ENRICHMENT<-Function.PR.PATIENT.COVERAGE(MUTSIG.BRCA, c("genegg", "NEG.LOG.q"), "BRCA", composite.tables.1)
ggplot(MUTSIG.PATIENT.ENRICHMENT,  aes(RECALL, PRECISION)) + geom_line() + geom_point()+theme.format

Table.u.BRCA.p$STATS$RANK<--log()
p.FACTOR.PATIENT.ENRICHMENT<-Function.PR.PATIENT.COVERAGE(Table.u.BRCA.p$STATS, c("Hugo_Symbol", "paired.t.stat"),"BRCA", composite.tables.1)

HYPER.2.CALC.PATIENT.ENRICHMENT<-Function.PR.PATIENT.COVERAGE(HYPER.TABLE.2.CALC, c("Hugo_Symbol", "RANK"),"BRCA", composite.tables.1)

log.1.P<-copy(test.FACTOR.PATIENT.ENRICHMENT)
log.1.P$METHOD<-"v.CNV"
log.2.P<-copy(MUTSIG.PATIENT.ENRICHMENT)
log.2.P$METHOD<-"MUTSIG"
log.3.P<-copy(p.FACTOR.PATIENT.ENRICHMENT)
log.3.P$METHOD<-"all.paired"
log.4.P<-copy(HYPER.2.CALC.PATIENT.ENRICHMENT)
log.4.P$METHOD<-"HYPER.2"

log.P<-rbind(log.1.P, log.2.P, log.3.P, log.4.P)

ggplot(log.P, aes(x=RECALL, y=PRECISION, colour=METHOD)) + geom_line() + theme.format + geom_point(size=4)

########081014#######
BRCA.v.p.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080914.BRCA.Table.v.p.rds")
BRCA.v.p.2$v.table
BRCA.v.p.2$neg.list[1:2]
BRCA.v.p.2$neg.list$A1CF[1:10]

write.table(BRCA.v.p.2$v.table, file="PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p",quote=F,sep="\t",col.names=T,row.names=F)
lapply(BRCA.v.p.2$neg.list, cat,"\n",file="PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p.neg", append=T)
write.table(names(BRCA.v.p.2$neg.list), file="PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p.neg.names", sep="\t", col.names=F, quote=F,row.names=F)

a<-BRCA.v.p.2$v.table[,c(2,5,6),with=F]
a$NEG<-a$SUM.LogFC-a$POS.SUM.LogFC
a<-melt(a, id.vars="v.PROTEIN")
b<-BRCA.v.p.2$v.table[,c(2,3), with=F]
b$v.PROTEIN.NEG<-b$v.PROTEIN-b$v.PROTEIN.POS
b$DIFF<-b$v.PROTEIN.POS-b$v.PROTEIN.NEG
b<-melt(b, id.vars="v.PROTEIN")
b$variable<-as.character(b$variable)

ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN, POS.SUM.LogFC)) + geom_point() + theme.format
ggplot(a, aes(v.PROTEIN, abs(value), colour=variable)) + geom_point() + theme.format +theme(legend.position="bottom")
ggplot(b, aes(v.PROTEIN, value, colour=variable)) + geom_point() + theme.format
cor.test(x=subset(b,variable=="DIFF")$v.PROTEIN, y=subset(b,variable=="DIFF")$value, method="pearson")

ggplot(BRCA.v.p.2$v.table, aes(y=POS.SUM.LogFC/v.PROTEIN.POS, x=SUM.LogFC/v.PROTEIN)) + geom_point() + theme.format 


test<-BRCA.v.p.2$neg.list
test.comb<-t(combn(names(test), 2, simplify=T))
library(RBGL)
library(gRbase)
test.comb.1<-t(combnPrim(names(test),2,simplify=T))
test.comb.1<-as.data.table(test.comb.1)
test.comb.1<-unique(test.comb.1)

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("test.comb", "test"),envir=environment())

test.ints<-system.time(parApply(cl,test.comb[1:100,], 1, function(x) length(intersect(test[[x[1]]],test[[x[2]]]))))
stopCluster(cl)

##########081414##############

#Analyze pairs of 20MILLLION BRCA v.diff samples
BRCA.NEG.20<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081414.BRCA.Table.v.p.neg.20.pairs",header=T, sep="\t"))
ggplot(BRCA.NEG.20, aes(NEG.INTERSECTION)) + geom_histogram() + theme.format

BRCA.v.p.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080914.BRCA.Table.v.p.rds")
BRCA.v.p.2$v.table

test<-as.data.table(merge(as.data.frame(BRCA.NEG.20), as.data.frame(BRCA.v.p.2$v.table), by.x="GENE.1", by.y="Hugo_Symbol"))
test<-as.data.table(merge(as.data.frame(test[,c(1,2,3,4),with=F]), as.data.frame(BRCA.v.p.2$v.table[,1:2,with=F]), by.x="GENE.2", by.y="Hugo_Symbol"))
test$v.diff<-abs(test$v.PROTEIN.x- test$v.PROTEIN.y)
test$DIFF.GROUP<-as.vector(cut(test$v.diff, 5))

ggplot(test, aes(v.diff)) + geom_histogram() + theme.format
ggplot(test, aes(x=DIFF.GROUP, y=NEG.INTERSECTION)) + geom_boxplot() + theme.format +
  theme(axis.text.x=element_text(angle=90))

#Randomly pick 100K pairs x 100, rbind and plot
library(RBGL)
library(gRbase)
BRCA.NEG.LIST<-BRCA.v.p.2$neg.list
test.comb<-t(combnPrim(names(BRCA.NEG.LIST),2,simplify=T))
test.comb<-as.data.table(test.comb)

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("test.comb", "BRCA.NEG.LIST"),envir=environment())

test.ints<-data.frame(V1=character(), V2=character, NEG.INTERSECTION=numeric())

for (rep in 1:100) {
  random.rows<-sample(1:nrow(test.comb),100000)
  PAIRS<-test.comb[random.rows,]
  PAIRS.INTERSECTION<-parApply(cl,test.comb[random.rows,], 1, function(x) length(intersect(BRCA.NEG.LIST[[x[1]]],BRCA.NEG.LIST[[x[2]]])))
  PAIRS$NEG.INTERSECTION<-PAIRS.INTERSECTION
  test.ints<-rbind(test.ints, PAIRS)
}
stopCluster(cl)

ggplot(test.ints, aes(NEG.INTERSECTION)) + geom_histogram() + theme.format
test.ints.plus<-as.data.table(merge(as.data.frame(test.ints), as.data.frame(BRCA.v.p.2$v.table[,1:2,with=F]), by.x="V1", by.y="Hugo_Symbol"))
test.ints.plus<-as.data.table(merge(as.data.frame(test.ints.plus), as.data.frame(BRCA.v.p.2$v.table[,1:2,with=F]), by.x="V2", by.y="Hugo_Symbol"))
test.ints.plus$v.diff<-abs(test.ints.plus$v.PROTEIN.x - test.ints.plus$v.PROTEIN.y)
ggplot(test.ints.plus, aes(v.diff)) + geom_histogram() + theme.format

test.ints.plus$DIFF.GROUP<-as.vector(cut(test.ints.plus$v.diff, 10))
ggplot(test.ints.plus, aes(DIFF.GROUP, NEG.INTERSECTION)) + geom_boxplot() + theme.format +
  theme(axis.text.x=element_text(angle=90))
ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN)) + geom_histogram() + theme.format

#Bin by location of both v.PROTEIN pair
test.ints.plus$v.PROTEIN.LOC<-ifelse(test.ints.plus$v.PROTEIN.x<0.6 & test.ints.plus$v.PROTEIN.y<0.6, "(0.5,0.6]",
                                     ifelse(test.ints.plus$v.PROTEIN.x>0.6 & test.ints.plus$v.PROTEIN.y>0.6 & 
                                              test.ints.plus$v.PROTEIN.x<0.65 & test.ints.plus$v.PROTEIN.y<0.65,"(0.6,0.65]",
                                            ifelse(test.ints.plus$v.PROTEIN.x>0.65 & test.ints.plus$v.PROTEIN.y>0.65 &
                                                     test.ints.plus$v.PROTEIN.x<0.70 & test.ints.plus$v.PROTEIN.y<0.70, "(0.65,0.70]",
                                                   ifelse(test.ints.plus$v.PROTEIN.x>0.70 & test.ints.plus$v.PROTEIN.y>0.70 &
                                                            test.ints.plus$v.PROTEIN.x<0.75 & test.ints.plus$v.PROTEIN.y<0.75, "(0.70,0.75]",
                                                          ifelse(test.ints.plus$v.PROTEIN.x>0.75 & test.ints.plus$v.PROTEIN.y>0.75, "(0.70,0.85)",
                                                                 "None")))))

test.ints.plus$bin.x<-as.vector(cut(test.ints.plus$v.PROTEIN.x,4))
test.ints.plus$bin.y<-as.vector(cut(test.ints.plus$v.PROTEIN.y,4))
intersect(as.vector(unique(test.ints.plus$bin.x)) , as.vector(unique(test.ints.plus$bin.y)))
test.ints.plus[bin.x==bin.y,]
test.ints.plus[v.PROTEIN.LOC!="None",]

ggplot(test.ints.plus[v.PROTEIN.LOC!="None",], aes(DIFF.GROUP, NEG.INTERSECTION)) + geom_boxplot() + theme.format +
  theme(axis.text.x=element_text(angle=90)) + facet_wrap(~v.PROTEIN.LOC)

ggplot(test.ints.plus[bin.x==bin.y,], aes(DIFF.GROUP, NEG.INTERSECTION)) + geom_boxplot() + theme.format + geom_jitter(size=0.5)+
  theme(axis.text.x=element_text(angle=90)) + facet_wrap(~bin.x) + theme(strip.text.x = element_text(size = 12))

ggplot(test.ints.plus[bin.x==bin.y,], aes(v.diff, NEG.INTERSECTION)) + geom_point(size=0.3) + theme.format +
  theme(axis.text.x=element_text(angle=90)) + facet_wrap(~bin.x)

#Check this accumulation as we increase v(p) - Calculated in python
BRCA.NEG.ACCUM<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081514.BRCA.Table.v.p.neg.cums", header=T, sep="\t", stringsAsFactors=F))

ggplot(BRCA.NEG.ACCUM, aes(v.PROTEIN, CUM_NEG))+ geom_point() +theme.format

BRCA.NEG.ACCUM.UNSORTED<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081514.BRCA.Table.v.p.neg.cums.unsorted", header=T, sep="\t", stringsAsFactors=F))
ggplot(BRCA.NEG.ACCUM.UNSORTED, aes(v.PROTEIN, CUM_NEG)) + geom_point() + theme.format

#######081514#######
BRCA.v.p.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080914.BRCA.Table.v.p.rds")

BRCA.NEG.CUM.INTERSECTION<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081514.BRCA.Table.v.p.neg.cums.intersection", header=T, sep="\t", stringsAsFactors=F))

ggplot(BRCA.NEG.CUM.INTERSECTION, aes(v.PROTEIN, CUM_NEG_INTERSECTION)) + geom_point() + theme.format

#######081614#######
BRCA.Table.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds")
BRCA.Table.1$table.1[Hugo_Symbol=="TP53",]

BRCA.CANCER.MATRICES<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES")
BRCA.CANCER.MATRICES$tumor[1:3,1:3]

BRCA.TP53.PATIENTS<-as.vector(unique(BRCA.Table.1$table.1[Hugo_Symbol=="TP53",]$Tumor_Sample_Barcode))
BRCA.TP53.PATIENTS<-unlist(lapply(BRCA.TP53.PATIENTS, function(x) paste0(strsplit(x,"-")[[1]][1:4],collapse="." ) ))
BRCA.TP53.PATIENTS<-BRCA.TP53.PATIENTS[BRCA.TP53.PATIENTS %in% colnames(BRCA.CANCER.MATRICES$tumor)]
dim(BRCA.CANCER.MATRICES$tumor)
dim(BRCA.CANCER.MATRICES$tumor[,BRCA.TP53.PATIENTS])

BRCA.TP53.CANCER.MATRIX<-BRCA.CANCER.MATRICES$tumor[,BRCA.TP53.PATIENTS]
head(BRCA.TP53.CANCER.MATRIX[,sample(1:305,5)])


nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "BRCA.CANCER.MATRICES", "BRCA.TP53.CANCER.MATRIX","BRCA.TP53.PATIENTS"),envir=environment())

BRCA.TP53.DIFF.RAN.300<-parLapply(cl,1:100, function(x)  
  Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.TP53.CANCER.MATRIX[,sample(BRCA.TP53.PATIENTS, 300)]) )
stopCluster(cl)

BRCA.TP53.DIFF.RAN.20.POS<-sapply(BRCA.TP53.DIFF.RAN.20, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC>0,]))
hist(BRCA.TP53.DIFF.RAN.20.POS)
BRCA.TP53.DIFF.RAN.20.NEG<-sapply(BRCA.TP53.DIFF.RAN.20, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC<0,]))
hist(BRCA.TP53.DIFF.RAN.20.NEG)

BRCA.TP53.DIFF.RAN.50.POS<-sapply(BRCA.TP53.DIFF.RAN.50, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC>0,]))
hist(BRCA.TP53.DIFF.RAN.50.POS)
BRCA.TP53.DIFF.RAN.50.NEG<-sapply(BRCA.TP53.DIFF.RAN.50, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC<0,]))
hist(BRCA.TP53.DIFF.RAN.50.NEG)

BRCA.TP53.DIFF.RAN.100.POS<-sapply(BRCA.TP53.DIFF.RAN.100, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC>0,]))
hist(BRCA.TP53.DIFF.RAN.100.POS)
BRCA.TP53.DIFF.RAN.100.NEG<-sapply(BRCA.TP53.DIFF.RAN.100, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC<0,]))
hist(BRCA.TP53.DIFF.RAN.100.NEG)

BRCA.TP53.DIFF.RAN.200.POS<-sapply(BRCA.TP53.DIFF.RAN.200, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC>0,]))
hist(BRCA.TP53.DIFF.RAN.200.POS)
BRCA.TP53.DIFF.RAN.200.NEG<-sapply(BRCA.TP53.DIFF.RAN.200, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC<0,]))
hist(BRCA.TP53.DIFF.RAN.200.NEG)

BRCA.TP53.DIFF.RAN.300.POS<-sapply(BRCA.TP53.DIFF.RAN.300, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC>0,]))
hist(BRCA.TP53.DIFF.RAN.300.POS)
BRCA.TP53.DIFF.RAN.300.NEG<-sapply(BRCA.TP53.DIFF.RAN.300, function(x) nrow(x[x$adj.P.Val<0.05 &  x$logFC<0,]))
hist(BRCA.TP53.DIFF.RAN.300.NEG)

BRCA.TP53.DIFF.RAN.POS<-as.data.table(do.call(cbind, lapply(c("20","50","100","200","300"), function(x) get(paste0("BRCA.TP53.DIFF.RAN.",x,".POS")))))
setnames(BRCA.TP53.DIFF.RAN.POS, c("20","50","100","200","300"))
BRCA.TP53.DIFF.RAN.POS<-melt(BRCA.TP53.DIFF.RAN.POS)
BRCA.TP53.DIFF.RAN.NEG<-as.data.table(do.call(cbind, lapply(c("20","50","100","200","300"), function(x) get(paste0("BRCA.TP53.DIFF.RAN.",x,".NEG")))))
setnames(BRCA.TP53.DIFF.RAN.NEG, c("20","50","100","200","300"))
BRCA.TP53.DIFF.RAN.NEG<-melt(BRCA.TP53.DIFF.RAN.NEG)

ggplot(BRCA.TP53.DIFF.RAN.POS, aes(value, fill=variable)) + geom_histogram(binwidth = 10) + theme.format + 
  geom_vline(xintercept=nrow(BRCA.TP53.DIFF.POS), colour="red")
ggplot(BRCA.TP53.DIFF.RAN.NEG, aes(value, fill=variable)) + geom_histogram(binwidth = 10) + theme.format +
  geom_vline(xintercept=nrow(BRCA.TP53.DIFF.NEG), colour="red")

BRCA.TP53.DIFF<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.TP53.CANCER.MATRIX)
BRCA.TP53.DIFF.POS<-BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05 & BRCA.TP53.DIFF$logFC>0,]
BRCA.TP53.DIFF.NEG<-BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05 & BRCA.TP53.DIFF$logFC<0,]
nrow(BRCA.TP53.DIFF.POS)
nrow(BRCA.TP53.DIFF.NEG)

####Make function for gene discovery as above per gene at 10%, 30% and 90% of number of samples covered per gene
#Apply this function to 10 random COSMIC genes and 10 random non-COSMIC genes and analyze behavior of POS and particularly NEG closeness of random sampled distributions to true value

ggplot(BRCA.v.p.2$v.table, aes(SAMPLE.COVERAGE)) + geom_histogram(binwidth=.005) + theme.format
ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN, SAMPLE.COVERAGE)) + geom_point() + theme.format
ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN.POS, SAMPLE.COVERAGE)) + geom_point() + theme.format
ggplot(BRCA.v.p.2$v.table, aes(SUM.LogFC, SAMPLE.COVERAGE)) + geom_point() + theme.format
ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN-v.PROTEIN.POS, SAMPLE.COVERAGE)) + geom_point() + theme.format
ggplot(BRCA.v.p.2$v.table, aes(v.PROTEIN, (abs(SUM.LogFC-POS.SUM.LogFC) + POS.SUM.LogFC))) + geom_point() + theme.format

BRCA.v.test<-BRCA.v.p.2$v.table
BRCA.v.test$NORM.v<-BRCA.v.test$v.PROTEIN / BRCA.v.test$SAMPLE.COVERAGE
BRCA.v.test[order(NORM.v, decreasing=T),]
BRCA.v.test[order(v.PROTEIN, decreasing=T),]
BRCA.v.test[Hugo_Symbol=="CCND1",]

BRCA.v.test.enrich<-Function.p.rank.enrichment(BRCA.v.test, c("Hugo_Symbol", "NORM.v"), target.genes=as.vector(COSMIC.BRCA$Symbol),normalize=F)
ggplot(BRCA.v.test.enrich$CUM.TABLE, aes(x=CUM.RECALL, y=CUM.PRECISION)) + geom_point() + geom_line() + theme.format
BRCA.v.test

hist(sapply(BRCA.v.p.2$neg.list, function(x) length(  BRCA.v.p.2$neg.list[x] )))

BRCA.v.p.filtered.2000<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.2000",header=F, sep="\t",drop=3)
BRCA.v.p.filtered.2000[V1 %in% as.vector(COSMIC.BRCA$Symbol),]
BRCA.v.p.filtered.2000[order(V2, decreasing=T),]
as.vector(COSMIC.BRCA$Symbol)
hist(BRCA.v.p.filtered.2000$V2)

BRCA.v.p.filtered.2000.enrich<-Function.p.rank.enrichment(BRCA.v.p.filtered.2000, c("V1","V2"), target.genes=as.vector(COSMIC.BRCA$Symbol),normalize=F)
ggplot(BRCA.v.p.filtered.2000.enrich$CUM.TABLE, aes(x=CUM.RECALL, y=CUM.PRECISION)) + geom_point() + geom_line() + theme.format

BRCA.v.p.filtered.2000.enrich.patients<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.2000, c("V1","V2"),"BRCA", composite.tables.1)

ggplot(BRCA.v.p.filtered.2000.enrich.patients, aes(RECALL, PRECISION)) + geom_line() + geom_point()+theme.format

BRCA.u.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/072014.BRCA.Table.u.p.rds")
BRCA.u.p$STATS
BRCA.u.p$STATS["Hugo_Symbol"=="CCND1",]


#####
#Load RDS files 
dummy.table.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds")
dummy.table.1<-dummy.table.1$table.1[,1:2, with=F] #[Tumor_Sample_Barcode, Hugo_Symbol]

dummy.table.1$PATIENT<-sapply(as.character(dummy.table.1$Tumor_Sample_Barcode), function(x) paste(strsplit(x, "-")[[1]][1:4] , collapse="."))
dummy.table.1$Tumor_Sample_Barcode<-NULL

#Integrate CNV data and keep only those that have expression information
cnv.table<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.GISTIC.TH.2.rds")
dummy.table.1<-unique(as.data.table(rbind(dummy.table.1, cnv.table)))
dummy.table.1.split<-split(dummy.table.1, dummy.table.1$Hugo_Symbol)

dummy.table.1.split$A1CF
writeLines(sapply(names(dummy.table.1.split),function(x) paste(x,paste(dummy.table.1.split[[x]]$PATIENT,collapse=" "))),sep="\n", "PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081914.BRCA.Table.1.GENE.PATIENT" )

#####v(p) with patient pre-treatment for COSMIC prediction
BRCA.v.p.filtered.patients.5<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.5", header=T, sep="\t",
                                                   stringsAsFactors=F))
BRCA.v.p.filtered.patients.10<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.10", header=T, sep="\t",
                                                     stringsAsFactors=F))
BRCA.v.p.filtered.patients.15<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.15", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.20<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.20", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.25<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.25", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.30<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.30", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.35<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.35", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.40<-as.data.table(read.csv("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.40", header=T, sep="\t",
                                                      stringsAsFactors=F))
BRCA.v.p.filtered.patients.40[Hugo_Symbol %in% as.vector(COSMIC.BRCA$Symbol),]
head(BRCA.v.p.filtered.patients.40,10)

BRCA.v.p.filtered.patients.enrich.40<-Function.p.rank.enrichment(BRCA.v.p.filtered.patients.40, c("Hugo_Symbol","v.PROTEIN"), target.genes=as.vector(COSMIC.BRCA$Symbol),normalize=F)
ggplot(BRCA.v.p.filtered.patients.enrich.40$CUM.TABLE, aes(x=CUM.RECALL, y=CUM.PRECISION)) + geom_point() + geom_line() + theme.format

log.5<-BRCA.v.p.filtered.patients.enrich.5$CUM.TABLE[,c(5,6), with=F]
log.5$TYPE<-"5%"
log.10<-BRCA.v.p.filtered.patients.enrich.10$CUM.TABLE[,c(5,6), with=F]
log.10$TYPE<-"10%"
log.15<-BRCA.v.p.filtered.patients.enrich.15$CUM.TABLE[,c(5,6), with=F]
log.15$TYPE<-"15%"
log.20<-BRCA.v.p.filtered.patients.enrich.20$CUM.TABLE[,c(5,6), with=F]
log.20$TYPE<-"20%"
log.25<-BRCA.v.p.filtered.patients.enrich.25$CUM.TABLE[,c(5,6), with=F]
log.25$TYPE<-"25%"
log.30<-BRCA.v.p.filtered.patients.enrich.30$CUM.TABLE[,c(5,6), with=F]
log.30$TYPE<-"30%"
log.35<-BRCA.v.p.filtered.patients.enrich.35$CUM.TABLE[,c(5,6), with=F]
log.35$TYPE<-"35%"
log.40<-BRCA.v.p.filtered.patients.enrich.40$CUM.TABLE[,c(5,6), with=F]
log.40$TYPE<-"40%"
log.cosmic<-rbind(log.5, log.10,log.15, log.20, log.25, log.30, log.35, log.40)
ggplot(log.cosmic, aes(CUM.RECALL, CUM.PRECISION, colour=TYPE)) + geom_point() + geom_line() + theme.format

#####v(p) with patient pre-treatment for BREAST CANCER PREDICTION
BRCA.v.p.filtered.patients.enrich.patients.5<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.5, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.10<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.10, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.15<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.15, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.20<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.20, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.25<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.25, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.30<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.30, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)
BRCA.v.p.filtered.patients.enrich.patients.40<-Function.PR.PATIENT.COVERAGE(BRCA.v.p.filtered.patients.40, c("Hugo_Symbol", "v.PROTEIN"), "BRCA", composite.tables.1)

ggplot(BRCA.v.p.filtered.patients.enrich.patients.40, aes(RECALL, PRECISION)) + geom_point() + geom_line() + theme.format

log.5.patient<-BRCA.v.p.filtered.patients.enrich.patients.5[,2:3, with=F]
log.5.patient$TYPE<-"5%"
log.10.patient<-BRCA.v.p.filtered.patients.enrich.patients.10[,2:3, with=F]
log.10.patient$TYPE<-"10%"
log.15.patient<-BRCA.v.p.filtered.patients.enrich.patients.15[,2:3, with=F]
log.15.patient$TYPE<-"15%"
log.20.patient<-BRCA.v.p.filtered.patients.enrich.patients.20[,2:3, with=F]
log.20.patient$TYPE<-"20%"
log.25.patient<-BRCA.v.p.filtered.patients.enrich.patients.25[,2:3, with=F]
log.25.patient$TYPE<-"25%"
log.30.patient<-BRCA.v.p.filtered.patients.enrich.patients.30[,2:3, with=F]
log.30.patient$TYPE<-"30%"
log.40.patient<-BRCA.v.p.filtered.patients.enrich.patients.40[,2:3, with=F]
log.40.patient$TYPE<-"40%"
log.patient.pred<-rbind(log.5.patient, log.10.patient, log.15.patient, log.20.patient, log.25.patient, log.30.patient, log.40.patient)
ggplot(log.patient.pred, aes(RECALL, PRECISION, colour=TYPE)) + geom_point() + geom_line() + theme.format

######AFTER GM#####
BRCA.v.p.2$v.table

BRCA.CANCER.MATRICES<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES")
BRCA.DIFF.EXP<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor)
BRCA.DIFF.EXP.GENES<-rownames(BRCA.DIFF.EXP[BRCA.DIFF.EXP$adj.P.Val<0.05,])

head(BRCA.DIFF.EXP)
head(BRCA.TP53.DIFF)
length(rownames(BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05,]))
length(rownames(BRCA.DIFF.EXP[BRCA.DIFF.EXP$adj.P.Val<0.05,]))
length(intersect(rownames(BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05,]),rownames(BRCA.DIFF.EXP[BRCA.DIFF.EXP$adj.P.Val<0.05,])))

#Get patient coverage per gene
BRCA.v.test
BRCA.Table.1.PLUS<-copy(BRCA.Table.1$table.1)
BRCA.Table.1.PLUS$Tumor_Sample_Barcode<-as.character(BRCA.Table.1.PLUS$Tumor_Sample_Barcode)
BRCA.Table.1.PLUS$PATIENT<-sapply(BRCA.Table.1.PLUS$Tumor_Sample_Barcode, function(x) paste0(strsplit(x, "-")[[1]][1:4], collapse="."))
BRCA.Table.1.PLUS$Tumor_Sample_Barcode<-NULL

BRCA.CNV.Table.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.GISTIC.TH.2.rds")
BRCA.Table.1.PLUS<-unique(rbind(BRCA.Table.1.PLUS[,c(1,3), with=F], BRCA.CNV.Table.1))
BRCA.Table.1.PATIENT.COVERAGE<-BRCA.Table.1.PLUS[,list(PATIENT.COVERAGE=length(PATIENT)), by="Hugo_Symbol"]
as.data.table(table(BRCA.Table.1.PATIENT.COVERAGE$PATIENT.COVERAGE)) #A POTENTIAL 221 X100 PERMUTATIONS 

#Use number of patients from TP53 to create random distribution sampling from all cancer patients
BRCA.TP53.DIFF<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, 
                                                        BRCA.CANCER.MATRICES$tumor[, intersect(colnames(BRCA.CANCER.MATRICES$tumor),as.vector(BRCA.Table.1.PLUS[Hugo_Symbol=="TP53",]$PATIENT))])

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "BRCA.CANCER.MATRICES"),envir=environment())
BRCA.DIFF.EXP.SAMPLING.TP53<-parLapply(cl, 1:100, function(x)
  Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[,sample(colnames(BRCA.CANCER.MATRICES$tumor),319)] ))
stopCluster(cl)

BRCA.DIFF.EXP.SAMPLING.TP53.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.TP53, function(x) nrow(x[x$adj.P.Val<0.05,])/nrow(x))
hist(BRCA.DIFF.EXP.SAMPLING.TP53.DIST)
length(rownames(BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05,]))/nrow(BRCA.TP53.DIFF)

BRCA.DIFF.EXP.SAMPLING.TP53.TRUE.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.TP53, function(x) 
  length(intersect(rownames(x[x$adj.P.Val<0.05,]), BRCA.DIFF.EXP.GENES))/nrow(x[x$adj.P.Val<0.05,])   )
hist(BRCA.DIFF.EXP.SAMPLING.TP53.TRUE.DIST)
length(intersect(rownames(BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05,]),BRCA.DIFF.EXP.GENES)) /
  length(rownames(BRCA.TP53.DIFF[BRCA.TP53.DIFF$adj.P.Val<0.05,]))

#Use number of patients from FGF19 to create random distribution sampling from all cancer patients
BRCA.FGF19.DIFF<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[, intersect(colnames(BRCA.CANCER.MATRICES$tumor),as.vector(BRCA.Table.1.PLUS[Hugo_Symbol=="FGF19",]$PATIENT))])

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "BRCA.CANCER.MATRICES"),envir=environment())

BRCA.DIFF.EXP.SAMPLING.FGF19<-parLapply(cl, 1:100, function(x)
  Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[,sample(colnames(BRCA.CANCER.MATRICES$tumor),159)] ))
stopCluster(cl)

BRCA.DIFF.EXP.SAMPLING.FGF19.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.FGF19, function(x) nrow(x[x$adj.P.Val<0.05,])/nrow(x))
hist(BRCA.DIFF.EXP.SAMPLING.FGF19.DIST)
length(rownames(BRCA.FGF19.DIFF[BRCA.FGF19.DIFF$adj.P.Val<0.05,]))/nrow(BRCA.FGF19.DIFF)

BRCA.DIFF.EXP.SAMPLING.FGF19.TRUE.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.FGF19, function(x) 
  length(intersect(rownames(x[x$adj.P.Val<0.05,]), BRCA.DIFF.EXP.GENES))/nrow(x[x$adj.P.Val<0.05,])   )
hist(BRCA.DIFF.EXP.SAMPLING.FGF19.TRUE.DIST)
length(intersect(rownames(BRCA.FGF19.DIFF[BRCA.FGF19.DIFF$adj.P.Val<0.05,]),BRCA.DIFF.EXP.GENES)) /
  length(rownames(BRCA.FGF19.DIFF[BRCA.FGF19.DIFF$adj.P.Val<0.05,]))

#Use number of patients from TITIN to create random distribution sampling from all cancer patients
BRCA.TTN.DIFF<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[, intersect(colnames(BRCA.CANCER.MATRICES$tumor),as.vector(BRCA.Table.1.PLUS[Hugo_Symbol=="TTN",]$PATIENT))])

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "BRCA.CANCER.MATRICES"),envir=environment())

BRCA.DIFF.EXP.SAMPLING.TTN<-parLapply(cl, 1:100, function(x)
  Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[,sample(colnames(BRCA.CANCER.MATRICES$tumor),180)] ))
stopCluster(cl)

BRCA.DIFF.EXP.SAMPLING.TTN.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.TTN, function(x) nrow(x[x$adj.P.Val<0.05,])/nrow(x))
hist(BRCA.DIFF.EXP.SAMPLING.TTN.DIST)
length(rownames(BRCA.TTN.DIFF[BRCA.TTN.DIFF$adj.P.Val<0.05,]))/nrow(BRCA.TTN.DIFF)

BRCA.DIFF.EXP.SAMPLING.TTN.TRUE.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.TTN, function(x) 
  length(intersect(rownames(x[x$adj.P.Val<0.05,]), BRCA.DIFF.EXP.GENES))/nrow(x[x$adj.P.Val<0.05,])   )
hist(BRCA.DIFF.EXP.SAMPLING.TTN.TRUE.DIST)
length(intersect(rownames(BRCA.TTN.DIFF[BRCA.TTN.DIFF$adj.P.Val<0.05,]),BRCA.DIFF.EXP.GENES)) /
  length(rownames(BRCA.TTN.DIFF[BRCA.TTN.DIFF$adj.P.Val<0.05,]))

#Use number of patients from A1CF to create random distribution sampling from all cancer patients
BRCA.A1CF.DIFF<-Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[, intersect(colnames(BRCA.CANCER.MATRICES$tumor),as.vector(BRCA.Table.1.PLUS[Hugo_Symbol=="A1CF",]$PATIENT))])

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression", "BRCA.CANCER.MATRICES"),envir=environment())

BRCA.DIFF.EXP.SAMPLING.A1CF<-parLapply(cl, 1:100, function(x)
  Function.RNAseq.Differential.Expression(BRCA.CANCER.MATRICES$normal, BRCA.CANCER.MATRICES$tumor[,sample(colnames(BRCA.CANCER.MATRICES$tumor),15)] ))
stopCluster(cl)

BRCA.DIFF.EXP.SAMPLING.A1CF.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.A1CF, function(x) nrow(x[x$adj.P.Val<0.05,])/nrow(x))
hist(BRCA.DIFF.EXP.SAMPLING.A1CF.DIST)
length(rownames(BRCA.A1CF.DIFF[BRCA.A1CF.DIFF$adj.P.Val<0.05,]))/nrow(BRCA.A1CF.DIFF)

BRCA.DIFF.EXP.SAMPLING.A1CF.TRUE.DIST<-sapply(BRCA.DIFF.EXP.SAMPLING.A1CF, function(x) 
  length(intersect(rownames(x[x$adj.P.Val<0.05,]), BRCA.DIFF.EXP.GENES))/nrow(x[x$adj.P.Val<0.05,])   )
hist(BRCA.DIFF.EXP.SAMPLING.A1CF.TRUE.DIST)
length(intersect(rownames(BRCA.A1CF.DIFF[BRCA.A1CF.DIFF$adj.P.Val<0.05,]),BRCA.DIFF.EXP.GENES)) /
  length(rownames(BRCA.TTN.DIFF[BRCA.A1CF.DIFF$adj.P.Val<0.05,]))

