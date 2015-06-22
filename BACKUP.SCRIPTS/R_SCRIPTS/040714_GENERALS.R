#Figures for generals proposal

#Non-silent mutations only
library(data.table)
library(grid)
BRCA.table<-read.csv("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t")
head(BRCA.table)
dim(BRCA.table)

BRCA.table<-BRCA.table[BRCA.table$Variant_Classification!="Silent",]
BRCA.table$Line_Number<-NULL
BRCA.table<-unique(BRCA.table)
BRCA.table<-BRCA.table[,c(1,16)]
BRCA.table<-as.data.table(BRCA.table)
length(unique(as.vector(BRCA.table$Hugo_Symbol)))

#Make table of mutated non-synonymous gene per patient
BRCA.gene.per.patient<-unique(BRCA.table)[,list(N.GENES=length(Hugo_Symbol)), by="Tumor_Sample_Barcode"] #Unique to account for multiple mutations per gene
ggplot(BRCA.gene.per.patient, aes(x=Tumor_Sample_Barcode, y=N.GENES)) + geom_bar(stat="identity", position="dodge") +
  scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

#Load KEGG enzymes and select patients that only have enzymes mutated in this dataset 
KEGG.FILE<-read.csv("DATABASES/KEGG/032614_ENZYME_PRODUCT", header=T, sep="\t")
KEGG.enzymes<-unique(KEGG.FILE$Enzyme)
length(KEGG.enzymes) #2702 enzymes, not filtered out with HMDB yet
KEGG.samples<-unique(BRCA.table)
KEGG.samples<-KEGG.samples[Hugo_Symbol %in% KEGG.enzymes,]
KEGG.samples<-KEGG.samples[,list(N.GENES.KEGG=length(Hugo_Symbol)), by="Tumor_Sample_Barcode"]
length(unique(KEGG.samples$Tumor_Sample_Barcode))

#Plot all genes per sample vs KEGG enzymes per sample
BRCA.KEGG.per.patient<-merge(BRCA.gene.per.patient, KEGG.samples, by="Tumor_Sample_Barcode")
BRCA.KEGG.per.patient<-BRCA.KEGG.per.patient[order(N.GENES, decreasing=T),]
BRCA.KEGG.per.patient.melted<-as.data.table(melt(BRCA.KEGG.per.patient))
BRCA.KEGG.per.patient.melted<-BRCA.KEGG.per.patient.melted[order(value, decreasing=T),]
BRCA.KEGG.per.patient.melted$Tumor_Sample_Barcode<-factor(BRCA.KEGG.per.patient.melted$Tumor_Sample_Barcode, levels=BRCA.KEGG.per.patient.melted$Tumor_Sample_Barcode,
                                                          ordered=T)
main.plot.2<-ggplot(BRCA.KEGG.per.patient.melted, aes(x=Tumor_Sample_Barcode, y=value, fill=variable)) + geom_histogram(stat="identity",binwidth=2,position="identity") + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title=element_blank(),
        legend.position="top", 
        axis.text.y=element_text(size=rel(2.0)), plot.title=element_text(size=rel(2)), legend.text = element_text(size = 22),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25)) +
  labs(title="KEGG Breast Cancer Coverage", size=9) + ylab("Number of Mutated Genes") + xlab("Patient Samples") +
  opts(axis.text.x = theme_blank()) +
  scale_fill_discrete(breaks=c("N.GENES.KEGG", "N.GENES"), labels=c("Mutated KEGG enzymes", "All mutated genes")) 

COVERAGE<-copy(BRCA.KEGG.per.patient)
COVERAGE$Coverage<-COVERAGE$N.GENES.KEGG/COVERAGE$N.GENES

sub.plot.2<-ggplot(COVERAGE, aes(x="Patients",y=Coverage)) + geom_boxplot() +
  theme(axis.title.x=element_blank(), axis.text.y=element_text(size=rel(2.0)),
        axis.title.x=element_text(size=20), axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=20))

vp <- viewport(width = 0.1, height = 0.4, x = 0.93, y = 0.68)

print (main.plot.2)
print (sub.plot.2, vp=vp)
  dev.off()

#Do differential expression in BRCA of tumor samples against normal
read.expression.files = function(path="DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/Level_3", files) {
  #There are "null" in table that function reads as NA instead to deal with easier
  tables = lapply(files, function(f) read.csv(paste0(path, "/", f), header=FALSE, sep="\t", skip=2, na.strings=c("null"),
                                                              col.names=c("gene", f)))
  
  m = tables[[1]]
  for (i in 2:length(tables)) {
    m = merge(m, tables[[i]], all=FALSE) #So that only common genes in tables are merged
  }
  rownames(m) = m$gene
  m$gene = NULL
  return (as.matrix(m))
}

BRCA.sample.map<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/FILE_SAMPLE_MAP.txt", header=T, sep="\t"))
BRCA.sample.map$type<-substr(BRCA.sample.map$barcode.s., 14,15)
BRCA.sample.map<-BRCA.sample.map[type %in% c("01","11")] #Remove metastatic for now
BRCA.exp.type<-split(BRCA.sample.map, BRCA.sample.map$type)
head(BRCA.sample.map)
head(BRCA.exp.type$'11')

breast.tumor<-read.expression.files(files=BRCA.exp.type$`01`$filename)
breast.normal<-read.expression.files(files=BRCA.exp.type$`11`$filename)
dim(breast.tumor)
breast.tumor[c(1,2),1:5]
head(colnames(breast.tumor))
head(BRCA.exp.type$`01`$filename)

colnames(breast.tumor) = paste0(colnames(breast.tumor), ".tumor")
colnames(breast.normal) = paste0(colnames(breast.normal), ".normal")

combined = cbind(breast.tumor, breast.normal[rownames(breast.tumor), ])

library(limma)
length(colnames(breast.tumor)); length(colnames(breast.normal))
d = data.frame(cancer=c(rep("tumor", 526), rep("normal", 61)))
library(GSEAL)
library(org.Hs.eg.db)

combined.filtered = combined[rownames(combined) %in% names(as.list(org.Hs.egSYMBOL2EG)), ] #Keep only those that have an entrez gene identifier
rownames(combined.filtered) = as.character(org.Hs.egSYMBOL2EG)[rownames(combined.filtered)] #change names in actual matrix

combined.filtered<-combined.filtered[complete.cases(combined.filtered),] #Remove NAs

boxplot(colSums(combined.filtered) ~d$cancer)
norm = normalizeBetweenArrays(combined.filtered, method="quantile") 
boxplot(colSums(norm) ~d$cancer)
dim(norm)
norm[1:5,1:5]
norm["366",524:528]
x.366<-mean(norm["1300",1:526])-mean(norm["1300",527:587]) #be wary of the sign in the logFC results

fit = lmFit(norm, model.matrix(~ cancer, d)) #fitting data 
dim(model.matrix(~cancer,d))
eb = eBayes(fit)
hist(eb$p.value[, 2])
volcanoplot(eb, highlight=20)

y = eb$coefficients[, 2]
head(eb$coefficients)
hist(y)
topTable(eb,coef=2,sort.by="p")

all.brca.fit<-topTable(eb, coef=2, n=Inf)
head(all.brca.fit)
nrow(all.brca.fit)
dim(all.brca.fit[all.brca.fit$adj.P.Val<0.01,])
tail(all.brca.fit,10)
head(all.brca.fit[all.brca.fit$adj.P.Val<0.05,],10)
all.brca.fit["1",]

mar<-all.brca.fit[all.brca.fit$adj.P.Val<0.05,]
mean(exp(mar$B)/(1+exp(mar$B))  )
mean(exp(all.brca.fit$B)/(1+exp(all.brca.fit$B))  )
nrow(all.brca.fit[all.brca.fit$adj.P.Val<0.01,])/ (nrow(all.brca.fit[all.brca.fit$adj.P.Val>=0.01,]) + nrow(all.brca.fit[all.brca.fit$adj.P.Val<0.01,]))
nrow(all.brca.fit[all.brca.fit$adj.P.Val<0.01,])/nrow(all.brca.fit)

#Distance to mutation plot - Comes from 013114_EXP_MUT_CORR.R
head(ALL_PATIENTS_DIST_MUT_MELTED)
ggplot(ALL_PATIENTS_DIST_MUT_MELTED[ALL_PATIENTS_DIST_MUT_MELTED$variable %in% c("GENES", "ALL_MUTATIONS") ,], aes(x=DISTANCE, y=value)) + geom_point(size=2.5) +
  facet_wrap(~variable) + scale_y_log10() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 30) , legend.text = element_text(size = 22),
        axis.text.y=element_text(size=rel(2.0)), axis.text.x=element_text(size=rel(2.0)),
        axis.title.x=element_text(size=25), axis.title.y=element_text(size=25)) +
  ylab("Mutation Counts")



