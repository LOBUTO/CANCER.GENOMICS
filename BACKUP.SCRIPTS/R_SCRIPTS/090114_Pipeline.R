history.database<-read.table("~/.rstudio-desktop/history_database", sep=":",fill=T,stringsAsFactors=F)
head(history.database)
history.database$V1

as.numeric("1401863172962")
as.numeric(as.POSIXct("2014-08-01 10:00:00 CET"))*1000


v.p.object.test<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")
length(names(v.p.object.test))

COSMIC<-as.data.table(read.csv("DATABASES/CANCER_DATA/COSMIC/cancer_gene_census.csv", header=T,sep=","))
COSMIC.BRCA<-COSMIC[apply(COSMIC, 1, function(x) any(grepl("breast",x,ignore.case=T))),]

#Split v.p.object.test into 10 files
saveRDS(v.p.object.test[1:2340], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.10.rds")
saveRDS(v.p.object.test[2341:4679], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.20.rds")
saveRDS(v.p.object.test[4680:7017], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.30.rds")
saveRDS(v.p.object.test[7018:9356], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.40.rds")
saveRDS(v.p.object.test[9357:11695], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.50.rds")
saveRDS(v.p.object.test[11696:14034], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.60.rds")
saveRDS(v.p.object.test[14035:16373], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.70.rds")
saveRDS(v.p.object.test[16374:18711], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.80.rds")
saveRDS(v.p.object.test[18712:21050], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.90.rds")
saveRDS(v.p.object.test[21051:23389], "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082714.BRCA.v3.tables.list.100.rds")

v.p.object.test[1:5]
v.p.quantiles<-as.data.table(quantile(1:23389, probs=seq(0,1, by=0.1)), keep.rownames=T)
v.p.quantiles$V2<-round(v.p.quantiles$V2)

######AUGUST RECOVERY#######
BRCA.NORM.MATRICES.OBJ<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082614.CANCER.MATRICES.NORMALIZED.OBJ.rds") #As batch.test
BRCA.diff.exp.table<-Function.RNAseq.Differential.Expression.V2(BRCA.NORM.MATRICES.OBJ, BRCA.NORM.MATRICES.OBJ$cancer.patients)
saveRDS(BRCA.diff.exp.table, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/090114.BRCA.DIFF.EXP.TABLE.rds")

#Get TH counts per genes v(p)
v.p.table.mean.abs.logFC<-as.data.table(sapply(v.p.object.test, function(x) mean(abs(x$logFC))  ,USE.NAMES=T),keep.rownames=T)

hist(v.p.table.mean.abs.logFC$V2)
setnames(v.p.table.mean.abs.logFC, c("Hugo_Symbol", "mean.abs.logFC"))

BRCA.Table.1.PATIENT.COVERAGE
v.p.table.mean.abs.logFC<-as.data.table(merge(as.data.frame(v.p.table.mean.abs.logFC), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE), by="Hugo_Symbol"))

ggplot(v.p.table.mean.abs.logFC, aes(mean.abs.logFC, PATIENT.COVERAGE)) + geom_point() + theme.format

v.p.table.logFC.TH.count<-as.data.table(t(sapply(v.p.object.test, function(x) {
  y1.5<-sum(abs(x$logFC)>1.5)
  y2<-sum(abs(x$logFC)>2)
  y2.5<-sum(abs(x$logFC)>2.5)
  y3<-sum(abs(x$logFC)>3)
  y3.5<-sum(abs(x$logFC)>3.5)
  y4<-sum(abs(x$logFC)>4)
  return(c(y1.5,y2, y2.5, y3, y3.5, y4))
} , USE.NAMES=T)), keep.rownames=T)

setnames(v.p.table.logFC.TH.count, c("Hugo_Symbol", "1.5", "2.0","2.5","3.0", "3.5","4.0"))
saveRDS(v.p.table.logFC.TH.count, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082814.BRCA.v.p.logFC.abs.TH.count")

v.p.table.logFC.TH.count<-melt(v.p.table.logFC.TH.count, id="Hugo_Symbol")
BRCA.Table.1.PATIENT.COVERAGE<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PATIENT.COVERAGE.rds")
v.p.table.logFC.TH.count<-as.data.table(merge(as.data.frame(v.p.table.logFC.TH.count), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE),by="Hugo_Symbol"))

ggplot(v.p.table.logFC.TH.count, aes(PATIENT.COVERAGE, value, colour=variable)) + geom_point() + theme.format
ggplot(v.p.table.logFC.TH.count, aes(as.factor(PATIENT.COVERAGE), value, colour=variable)) + geom_boxplot() + theme.format+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Getting Null to apply to TH count found in genes v(p)
v.p.object.NULL<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v.p.object.NULL.rds")

v.p.table.null.logFC.TH.count<-lapply(names(v.p.object.NULL), function(x) {
  z=t(sapply(v.p.object.NULL[[x]], function(y) {
    y.all=abs(y$logFC)
    y.1.5<-sum(y.all>1.5)
    y.2<-sum(y.all>2)
    y.2.5<-sum(y.all>2.5)
    y.3<-sum(y.all>3)
    y.3.5<-sum(y.all>3.5)
    y.4<-sum(y.all>4)
    return(c(y.1.5, y.2, y.2.5, y.3, y.3.5, y.4))
  }))
  z<-as.data.table(z)
  setnames(z, c("1.5", "2.0", "2.5", "3.0", "3.5", "4.0"))
  z$PATIENT.COVERAGE<-x
  z<-melt(z, id="PATIENT.COVERAGE")
  return(z)
})

v.p.table.logFC.TH.count<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/LOGS/082814.BRCA.v.p.logFC.abs.TH.count")
names(v.p.table.null.logFC.TH.count)<-names(v.p.object.NULL)
v.p.table.null.logFC.TH.count[["19"]]
v.p.table.null.logFC.TH.count[["227"]]

v.p.table.null.logFC.TH.count.table<-do.call(rbind, v.p.table.null.logFC.TH.count)

ggplot(v.p.table.null.logFC.TH.count.table, aes(PATIENT.COVERAGE, value, colour=variable)) + geom_boxplot() + theme.format +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

v.p.table.null.logFC.TH.count.table$PATIENT.COVERAGE<-as.numeric(v.p.table.null.logFC.TH.count.table$PATIENT.COVERAGE)

ggplot(v.p.table.null.logFC.TH.count.table, aes(as.factor(PATIENT.COVERAGE), value, colour=variable)) + geom_boxplot() + theme.format +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Look at variance distribution of NULL
ggplot(v.p.table.null.logFC.TH.count.table[, list(variance=var(value)), by=c("PATIENT.COVERAGE", "variable")], aes(PATIENT.COVERAGE, variance, colour=variable)) +
  geom_point()+ theme.format

ggplot(v.p.table.null.logFC.TH.count.table[, list(variance=var(value)), by=c("PATIENT.COVERAGE", "variable")], aes(PATIENT.COVERAGE, variance)) +
  geom_point()+ theme.format + facet_wrap(~variable,ncol=1,scales="free_y") + theme(strip.text.x=element_text(size=20))

#Now compare null distribution per sample coverage against obtained count thresholded
v.p.table.logFC.TH.count
v.p.table.null.logFC.TH.count.table

ggplot(v.p.table.logFC.TH.count, aes(factor(PATIENT.COVERAGE), value)) + geom_point(shape=8) + geom_line(colour="red")+ theme.format +
  geom_boxplot(data=v.p.table.null.logFC.TH.count.table, aes(factor(PATIENT.COVERAGE), value)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + facet_wrap(~variable,scales="free_y") + theme(strip.text.x=element_text(size=20)) +
  ylab("Number of Differentially Expressed Genes Beyond Threshold")

#Now obtain actual p-values per gene based on applied null distribution
setnames(v.p.table.null.logFC.TH.count.table, c("NULL.PATIENT.COVERAGE", "null.variable", "null.value"))
v.p.table.logFC.TH.count
v.p.table.null.logFC.TH.count.table

v.p.table.logFC.TH.count.p.val.upper<-v.p.table.logFC.TH.count[Hugo_Symbol %in% COSMIC.BRCA$Symbol,][,list(p.val= mean(as.vector(v.p.table.null.logFC.TH.count.table[as.numeric(null.variable)==as.numeric(variable) &
                                                                                                                           as.numeric(NULL.PATIENT.COVERAGE)==PATIENT.COVERAGE,]$null.value)>=value)),
                                                         by=c("Hugo_Symbol","variable")]
v.p.table.logFC.TH.count.p.val.upper$p.val.adj<-p.adjust(v.p.table.logFC.TH.count.p.val.upper$p.val, method="fdr")
v.p.table.logFC.TH.count.p.val.upper[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

v.p.table.logFC.TH.count.p.val.lower<-v.p.table.logFC.TH.count[Hugo_Symbol %in% COSMIC.BRCA$Symbol,][,list(p.val= mean(as.vector(v.p.table.null.logFC.TH.count.table[as.numeric(null.variable)==as.numeric(variable) &
                                                                                                                                 as.numeric(NULL.PATIENT.COVERAGE)==PATIENT.COVERAGE,]$null.value)<=value)),
                                                               by=c("Hugo_Symbol","variable")]
v.p.table.logFC.TH.count.p.val.lower$p.val.adj<-p.adjust(v.p.table.logFC.TH.count.p.val.lower$p.val, method="fdr")
v.p.table.logFC.TH.count.p.val.lower[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]

ggplot(v.p.table.logFC.TH.count.p.val.upper[Hugo_Symbol %in% COSMIC.BRCA$Symbol,], aes(Hugo_Symbol, p.val.adj, colour=variable)) + geom_bar(stat="identity", position="dodge") +
  theme.format +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0.025, color="red")

ggplot(v.p.table.logFC.TH.count.p.val.lower[Hugo_Symbol %in% COSMIC.BRCA$Symbol,], aes(Hugo_Symbol, p.val.adj, colour=variable)) + geom_bar(stat="identity", position="dodge") +
  theme.format +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0.025, color="red")

#Look at logFC expressions at different thresholds of house-keeping genes
HK.genes<-as.data.table(read.csv("DATABASES/HK_genes.EISENBERG.txt", header=F, sep="\t", stringsAsFactors=F))
setnames(HK.genes, c("Hugo_Symbol", "Refseq"))
HK.genes

HK.genes$Hugo_Symbol<-sapply(HK.genes$Hugo_Symbol, function(x) strsplit(x," ")[[1]][1])
BRCA.diff.exp.table$House.keeping<-as.character(BRCA.diff.exp.table$ID) %in% as.character(as.vector(HK.genes$Hugo_Symbol))

ggplot(BRCA.diff.exp.table, aes(logFC)) + geom_histogram(, binwidth=0.01) + theme.format + facet_wrap(~House.keeping,ncol=1)
var.test(as.vector(BRCA.diff.exp.table[House.keeping==FALSE,]$logFC), as.vector(BRCA.diff.exp.table[House.keeping==TRUE,]$logFC))

#Analyze variability (variance) in expression for all v(p) by looking at mean variance for all genes across all patients affected by gene v(p)
Table.1.PLUS<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PLUS.rds")
v.p.corr.matrix.obj<-Function.v.p.corr.matrix.obj(BRCA.NORM.MATRICES.OBJ,
                                                  Table.1.PLUS[Hugo_Symbol %in% c("TTN",
                                                                                  sample(unique(as.vector(Table.1.PLUS$Hugo_Symbol)),1000),
                                                                                  as.vector(COSMIC.BRCA$Symbol)),], as.vector(HK.genes$Hugo_Symbol) )
as.data.table(v.p.corr.matrix.obj, keep.rownames=T)[rn %in% COSMIC.BRCA$Symbol,]
as.data.table(v.p.corr.matrix.obj, keep.rownames=T) 
ggplot(as.data.table(merge(as.data.frame(as.data.table(v.p.corr.matrix.obj, keep.rownames=T)), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE), by.x="rn", by.y="Hugo_Symbol")),
       aes(v.p.corr.matrix.obj, PATIENT.COVERAGE)) + geom_point() + theme.format

#Get intersection recall and coverage per v(p) gene at different thresholds to whole cancer threshold sets
v.p.object.test<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")

v.p.prediction.table<-Function.v.p.prediction(v.p.object.test[1:100], BRCA.diff.exp.table) #DID ALL IN CLUSTER

#######090214#######
#Get new v(p) table based on fully normalized matrices for breast cancer
v.p.object.BRCA<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")
v.p.table.BRCA.V2<-Function.v.p.table.V2(v.p.object.BRCA) #DONE IN CLUSTER

#Instead of v(p) + mean var (v.p.corr.matrix.obj), thresholded + mean var??

####TO DO
#1. Get intersection recall and coverage per v(p) gene at different thresholds to whole cancer threshold sets - NEED NULL!!! MAYBE?? LOOK AT COSMIC PATTERN
BRCA.PRED.INT.TH.v.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/090114.BRCA.PRED.INT.TH.rds")
BRCA.PRED.INT.TH.v.p[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
BRCA.PRED.INT.TH.v.p$a1.5<-BRCA.PRED.INT.TH.v.p$p1.5 + BRCA.PRED.INT.TH.v.p$c1.5
BRCA.PRED.INT.TH.v.p$a2.0<-BRCA.PRED.INT.TH.v.p$p2.0 + BRCA.PRED.INT.TH.v.p$c2.0
BRCA.PRED.INT.TH.v.p$a2.5<-BRCA.PRED.INT.TH.v.p$p2.5 + BRCA.PRED.INT.TH.v.p$c2.5
BRCA.PRED.INT.TH.v.p$a3.0<-BRCA.PRED.INT.TH.v.p$p3.0 + BRCA.PRED.INT.TH.v.p$c3.0
BRCA.PRED.INT.TH.v.p$a3.5<-BRCA.PRED.INT.TH.v.p$p3.5 + BRCA.PRED.INT.TH.v.p$c3.5
BRCA.PRED.INT.TH.v.p$a4.0<-BRCA.PRED.INT.TH.v.p$p4.0 + BRCA.PRED.INT.TH.v.p$c4.0

BRCA.PRED.INT.TH.v.p.melt<-melt(BRCA.PRED.INT.TH.v.p[,c(1,14:19),with=F], id="Hugo_Symbol")
BRCA.PRED.INT.TH.v.p.melt<-as.data.table(merge(as.data.frame(BRCA.PRED.INT.TH.v.p.melt), as.data.frame(BRCA.Table.1.PATIENT.COVERAGE), by="Hugo_Symbol"))

ggplot(BRCA.PRED.INT.TH.v.p.melt, aes(PATIENT.COVERAGE, value, colour=variable)) + geom_point() + theme.format #UNFORTUNATELY IT DEPENDS ON SAMPLE SIZE - NOT USEFUL!

BRCA.PRED.INT.TH.v.p.melt[Hugo_Symbol %in% COSMIC.BRCA$Symbol,] #DO NULL FOR ONLY COSMIC GENES, MAKE FUNCTION

#Get NULL
BRCA.diff.exp.table<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/090114.BRCA.DIFF.EXP.TABLE.rds")
BRCA.PRED.INT.TH.v.p.NULL<-Function.v.p.prediction.NULL(v.p.object.NULL, BRCA.diff.exp.table)

#2. Given the mean variances per gene, how likely is it that it will be lower than null (since lower mean variance implies that patients chosen by gene are more similar
#   to each other), NEED TO BUILD NULL FOR THIS - SAMPLE SIZE INDEPENDENT 

