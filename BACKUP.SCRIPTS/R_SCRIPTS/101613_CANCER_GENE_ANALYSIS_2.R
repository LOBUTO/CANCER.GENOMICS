#CANCER_GENES_ANALYSIS 2
#08/28/13
#FREQUENCY BARPLOTS OF CANCER MUTATIONS
#DON'T WORRY ABOUT ZEROES, THEY WERE REMOVED DURING PYTHON CODING AS ZERO FREQUENCY OF GENE MEANS THAT THE GENE DOES NOT BELONG TO CANCER CATEGORY

CANCER.FREQ=read.table("NETWORK/R_ANALYSIS/082813_CANCER_GENE_FREQUENCY_UNIPARTITE", stringsAsFactors=TRUE,
                       fill=TRUE,sep="\t")
CANCER.FREQ.T=t(CANCER.FREQ)
CANCER.FREQ.T=data.frame(CANCER.FREQ.T)
head(CANCER.FREQ.T)

colnames(CANCER.FREQ.T)<-c("COAD", "LUSC","READ","GBM","UCEC","KIRC","BRCA","OV")
CANCER.FREQ.T=CANCER.FREQ.T[-1,]

COAD=as.vector(CANCER.FREQ.T$COAD)
COAD=as.numeric(COAD[!is.na(COAD)])

LUSC=as.vector(CANCER.FREQ.T$LUSC)
LUSC=as.numeric(LUSC[!is.na(LUSC)])

READ=as.vector(CANCER.FREQ.T$READ)
READ=as.numeric(READ[!is.na(READ)])

GBM=as.vector(CANCER.FREQ.T$GBM)
GBM=as.numeric(GBM[!is.na(GBM)])

UCEC=as.vector(CANCER.FREQ.T$UCEC)
UCEC=as.numeric(UCEC[!is.na(UCEC)])

KIRC=as.vector(CANCER.FREQ.T$KIRC)
KIRC=as.numeric(KIRC[!is.na(KIRC)])

BRCA=as.vector(CANCER.FREQ.T$BRCA)
BRCA=as.numeric(BRCA[!is.na(BRCA)])

OV=as.vector(CANCER.FREQ.T$OV)
OV=as.numeric(OV[!is.na(OV)])

COLORS<-c("black", "blue", "green", "red", "purple", "yellow", "brown", "grey")
boxplot(list(LUSC,READ,OV,GBM,UCEC,KIRC,BRCA,COAD), cex=0.7,
        names=c("LUSC","READ","OV","GBM","UCEC","KIRC","BRCA","COAD"), cex.axis=0.6,
        log="y", col=COLORS,
        main="Cancer Frequency of Mutation",xlab="Cancer Type", ylab="Raw Mutation Frequency")

#GET CLUSTERS 41, 43B and 637 (43A had only one value)
C.41=c(23, 15, 16, 10, 34, 8, 17, 13, 24, 15, 17, 24, 22, 5, 24, 35, 12, 12, 3, 8)
C.43=c(27, 48, 7, 96, 12, 9, 20, 8, 2, 12, 28, 26, 6, 65, 1, 3, 5, 12, 27)
C.637=c(12, 20, 1, 7, 4, 12, 5, 13, 30, 34, 18, 49, 20, 20, 26, 24, 32, 10, 52, 107, 23, 12, 11, 7, 11, 68, 16, 8, 10, 7, 13, 49, 24, 30, 10, 16, 6, 26, 6, 8, 7, 46, 90, 41, 45, 17, 302, 30, 17, 20, 38, 18, 56, 27, 13, 16, 36, 17, 18, 3, 15, 92, 36, 11, 6, 6, 6, 15, 16, 28, 11, 14, 37, 8, 20, 17, 9, 109, 4, 23, 3, 9, 55, 50, 12, 33, 24, 18, 2, 6, 19, 6, 54, 24, 8, 13, 23, 12, 37, 10, 39, 68, 10, 8, 4, 14, 24, 36, 17, 34, 39, 11, 17, 3, 4, 4, 6, 10, 9, 3, 6, 19, 12, 38, 37, 11, 24, 12, 6, 8, 2, 5, 29, 44, 18, 5, 10, 26, 9, 58, 33, 58, 39, 12, 99, 9, 120, 45, 461, 29, 78, 22, 86, 168, 30, 55, 35, 52, 20, 46, 13, 1, 2, 55, 99, 10, 44, 20, 13, 14, 24, 312, 1, 90, 9, 12, 58, 27, 8, 25, 68, 48, 48, 16, 26, 22, 4, 10, 15, 37, 30, 37, 58, 9, 1, 24, 12, 57, 27, 29, 7, 74, 19, 9, 8, 60, 7, 16, 4, 4, 9, 14, 5, 23, 16, 2, 11, 24, 81, 31, 48, 15, 36, 14, 30, 18, 24, 66, 2, 2, 43, 9, 15, 28, 21, 10, 18, 44, 15, 7, 10, 19, 28, 23, 31, 14, 22, 22, 46, 1, 12, 45, 28, 110, 23, 36, 11, 15, 45, 4, 2, 9, 32, 11, 255, 2, 14, 19, 10, 48, 58, 16, 2, 80, 88, 8, 18, 62, 14, 14, 41, 7, 30, 3, 13, 18, 4, 4, 9, 5, 36, 5, 26, 23, 13, 11)

#Compare them against UCEC and BRCA (HIGHEST AND LOWEST FREQUENCY SCORES BY CLUSTERING ANALYSIS)
COLORS2=c("Brown", "White", "White", "White", "Purple")
boxplot(list(BRCA, C.41, C.43, C.637, UCEC), cex=0.7,
        names=c("BRCA", "41", "43", "637", "UCEC"), cex.axis=0.6,
        log="y",col=COLORS2,
        main="Cancer Frequency of Mutation vs 41,43,637", xlab="Cancer Types and Clusters",
        ylab="Raw Mutation Frequency")