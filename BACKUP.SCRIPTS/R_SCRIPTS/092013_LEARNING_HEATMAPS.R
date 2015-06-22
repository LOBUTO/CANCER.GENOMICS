#LEARNING HEATMAPS
#092013 - MADE UP, DATE UNKNOWN
library("ALL")

data("ALL")
ALL
ALL$mol.biol

esset<-ALL[, ALL$mol.biol %in% c("BCR/ABL", "ALL1/AF4")] #filter 128 for criteria
head(exprs(esset[1:100,])) #matrix of gene expression (12625 genes) vs patient number (<128 patients)
exprs(esset[1:100,]) #<- class=MATRIX
heatmap(exprs(esset[1:100,])) #draw heatmap of MATRIX
        
#Look for signicance
library(limma)
esset$mol.biol
f<-factor(as.character(esset$mol.biol)) #Convert to factors
design<-model.matrix(~f) #Creates matrix containing column of 0/1 for presence or absence of one of the two factors in f
fit<-eBayes(lmFit(esset,design)) #Do statistics using whole filtered data set, constructs linear model fit
topTable(fit, coef=2)# Extract table from top ranked genes from linear model fit

#Select those who have p-value less than 0.05
selected<- p.adjust(fit$p.value[,2]) <0.05
essetSEL<-esset[selected,]

heatmap(exprs(essetSEL))
        
#LOAD CLUSTER FILE
FILE=read.table("NETWORK/TSF_FILES/CLUSTER_CANCER_AVG_LOGP",header=TRUE)
DATA=FILE[,-1]; DATA
CLUSTERS=FILE[,1]
CLUSTERS[1]="10A"
CLUSTERS[2]="10B"
CLUSTERS[3]="10C"
CLUSTERS[4]="11A"
CLUSTERS[4]="11B"
library(RColorBrewer)
hmcol<-brewer.pal(9,"Blues")
heatmap(as.matrix(DATA[50:62,]), col=hmcol)
FILE[,1]