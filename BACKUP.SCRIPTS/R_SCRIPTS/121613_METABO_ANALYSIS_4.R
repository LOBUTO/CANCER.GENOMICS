#METABO ANALYSIS #4
#12/16/13
#1. Look for optimal TC filtering to get better enrichemnt of metabolites across patients
#2. If cannot find optimal TC without loosing much genomic information then combine metabolites into metabolic modules with TC=1.0
#CONCLUSION = Will use a TC=1.0 so that we don't loose any information

#####
#Plot metabolic genes and non-metabolic genes per cancer
CANCER_FILTER80<-data.frame(CANCER=c("ACC", "CESC", "LUAD", "LUSC", "READ", "KIRP", "GBM", "LGG", "BRCA", "STAD", "UCEC",
                                       "PRAD", "PAAD", "KIRC", "THCA", "HNSC", "KICH", "BLCA", "SKCM"), 
                              CANCER_GENES=c(3516, 3800, 15441, 12193, 8034, 3499, 7004, 8237, 10907, 14988, 11186,
                                             5967, 7824, 8344, 3175, 11933, 3010, 3484, 15461),
                              METABOLIC_MUTATED_GENES= c(523, 508, 2220, 1758, 1140, 509, 1018, 1236, 1587, 2195, 1639,
                                                         856, 1154, 1268, 458, 1722, 425, 499, 2203),
                            METABOLITES=c(346, 380, 682, 618, 515, 372, 507, 528, 580, 670, 602, 
                                          443, 532, 558, 339, 599, 350, 360, 676)
                            )
CANCER_FILTER80
library(ggplot2)
MELTED_CANCER_GENE_COUNT<-melt(CANCER_FILTER80[,1:3])
MELTED_CANCER_GENE_COUNT
ggplot(MELTED_CANCER_GENE_COUNT, aes(CANCER, value, fill=variable)) +geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0)))

#Do hypergeometric per metabolite in table and fdr
METABOLITE_HYPER_TABLE<-read.table("NETWORK/R_ANALYSIS/121813_METABOLITES_FOR_HYPER", sep="\t", header=T, quote="")
head(METABOLITE_HYPER_TABLE)
METABOLITE_HYPER_TABLE$P.VALUE<-phyper(METABOLITE_HYPER_TABLE$TOTAL_INTERACTING_CANCER_P, METABOLITE_HYPER_TABLE$TOTAL_INTERACTING_P,
                                       METABOLITE_HYPER_TABLE$TOTAL_NOT_INTERACTING_P, METABOLITE_HYPER_TABLE$TOTAL_NOT_INTERACTING_CANCER_P,
                                       lower.tail=FALSE)
METABOLITE_HYPER_TABLE$P.VALUE.CORRECTED<-p.adjust(METABOLITE_HYPER_TABLE$P.VALUE, method="fdr")

#Keep those that have an FDR<0.05
length(unique(METABOLITE_HYPER_TABLE$METABOLITE))
METABOLITE_HYPER_SIG<-METABOLITE_HYPER_TABLE[METABOLITE_HYPER_TABLE$P.VALUE.CORRECTED<0.05, ]
length(unique(METABOLITE_HYPER_SIG$METABOLITE))

head(METABOLITE_HYPER_SIG)
METABOLITE_HYPER_SIG<-as.data.table(METABOLITE_HYPER_SIG)
METABOLITE_SIG_PER_CANCER<-METABOLITE_HYPER_SIG[, list(SIG_METABOLITES=length(unique(METABOLITE))), by="CANCER"]

MERGED_PER_CANCER_MET<-merge(as.data.frame(METABOLITE_SIG_PER_CANCER), CANCER_FILTER80, by=c("CANCER", "CANCER"))

#Plot number of metabolites per cancer (total and sig)
MELTED_MET_PER_CANCER<-melt(MERGED_PER_CANCER_MET[,c(1,2,5)])
MELTED_MET_PER_CANCER
ggplot(MELTED_MET_PER_CANCER, aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0)))

#Output table of those metabolites that are significant per cancer - Only need metabolites name and cancer type
write.table(METABOLITE_HYPER_SIG, file="NETWORK/FROM_R_ANALYSIS/121813_SIG_MET_PER_CANCER", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


PATIENT_SIG<-data.frame(CANCER=c("ACC", "CESC", "LUAD", "LUSC", "READ", "KIRP", "GBM", "LGG", "BRCA", "STAD", "UCEC",
                                 "PRAD", "PAAD", "KIRC", "THCA", "HNSC", "KICH", "BLCA", "SKCM"), 
                        ALL_PATIENTS=c(90, 39, 519, 176, 80, 112, 291, 220, 774, 243, 194,
                                       250, 57, 293, 403, 306, 66, 28, 345),
                        SIG_PATIENTS=c(60, 37, 506, 176, 78, 99, 277, 173, 682, 239, 193,
                                       189, 42, 283, 162, 296, 59, 28, 337)
                        )
PATIENT_SIG
MELTED_PATIENT_SIG<-melt(PATIENT_SIG)
ggplot(MELTED_PATIENT_SIG, aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0)))

####FOR 100
CANCER_FILTER100<-data.frame(CANCER=c("ACC", "CESC", "LUAD", "LUSC", "READ", "KIRP", "GBM", "LGG", "BRCA", "STAD", "UCEC",
                                     "PRAD", "PAAD", "KIRC", "THCA", "HNSC", "KICH", "BLCA", "SKCM"), 
                            CANCER_GENES=c(3516, 3800, 15441, 12193, 8034, 3499, 7004, 8237, 10907, 14988, 11186,
                                           5967, 7824, 8344, 3175, 11933, 3010, 3484, 15461),
                            METABOLIC_MUTATED_GENES= c(773, 833, 3346, 2691, 1773, 828, 1558, 1876, 2436, 3307, 2516,
                                                       1339, 1760, 1930, 716, 2665, 663, 809, 3332),
                            METABOLITES=c(846, 914, 1481, 1357, 1198, 915, 1146, 1219, 1320, 1467, 1337, 
                                          1045, 1208, 1261, 893, 1336, 861, 861, 1470)
)

MELTED_CANCER_GENE_COUNT100<-melt(CANCER_FILTER100[,1:3])
MELTED_CANCER_GENE_COUNT100
ggplot(MELTED_CANCER_GENE_COUNT100, aes(CANCER, value, fill=variable)) +geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0)))+
  labs(title="FILTER_100")

#Do hypergeometric per metabolite in table and fdr
METABOLITE_HYPER_TABLE100<-read.table("NETWORK/R_ANALYSIS/121813_METABOLITES_FOR_HYPER_100", sep="\t", header=T, quote="")
head(METABOLITE_HYPER_TABLE)
METABOLITE_HYPER_TABLE100$P.VALUE<-phyper(METABOLITE_HYPER_TABLE100$TOTAL_INTERACTING_CANCER_P, METABOLITE_HYPER_TABLE100$TOTAL_INTERACTING_P,
                                       METABOLITE_HYPER_TABLE100$TOTAL_NOT_INTERACTING_P, METABOLITE_HYPER_TABLE100$TOTAL_NOT_INTERACTING_CANCER_P,
                                       lower.tail=FALSE)
METABOLITE_HYPER_TABLE100$P.VALUE.CORRECTED<-p.adjust(METABOLITE_HYPER_TABLE100$P.VALUE, method="fdr")

#Keep those that have an FDR<0.05
length(unique(METABOLITE_HYPER_TABLE100$METABOLITE))
METABOLITE_HYPER_SIG100<-METABOLITE_HYPER_TABLE100[METABOLITE_HYPER_TABLE100$P.VALUE.CORRECTED<0.05, ]
length(unique(METABOLITE_HYPER_SIG100$METABOLITE))

head(METABOLITE_HYPER_SIG100)
METABOLITE_HYPER_SIG100<-as.data.table(METABOLITE_HYPER_SIG100)
METABOLITE_SIG_PER_CANCER100<-METABOLITE_HYPER_SIG100[, list(SIG_METABOLITES=length(unique(METABOLITE))), by="CANCER"]

MERGED_PER_CANCER_MET100<-merge(as.data.frame(METABOLITE_SIG_PER_CANCER100), CANCER_FILTER100, by=c("CANCER", "CANCER"))

#Plot number of metabolites per cancer (total and sig)
MELTED_MET_PER_CANCER100<-melt(MERGED_PER_CANCER_MET100[,c(1,2,5)])
MELTED_MET_PER_CANCER
ggplot(MELTED_MET_PER_CANCER100, aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0))) +
  labs(title="FILTER_100")

#Output table of those metabolites that are significant per cancer - Only need metabolites name and cancer type
write.table(METABOLITE_HYPER_SIG100, file="NETWORK/FROM_R_ANALYSIS/121813_SIG_MET_PER_CANCER_100", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

PATIENT_SIG100<-data.frame(CANCER=c("ACC", "CESC", "LUAD", "LUSC", "READ", "KIRP", "GBM", "LGG", "BRCA", "STAD", "UCEC",
                                 "PRAD", "PAAD", "KIRC", "THCA", "HNSC", "KICH", "BLCA", "SKCM"), 
                        ALL_PATIENTS=c(90, 39, 519, 176, 80, 112, 291, 220, 774, 243, 194,
                                       250, 57, 293, 403, 306, 66, 28, 345),
                        SIG_PATIENTS=c(78, 39, 514, 176, 80, 109, 288, 215, 745, 242, 193,
                                       237, 51, 292, 367, 305, 65, 28, 344)
)
PATIENT_SIG100
MELTED_PATIENT_SIG100<-melt(PATIENT_SIG100)
ggplot(MELTED_PATIENT_SIG100, aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") + theme(axis.text.x=element_text(size=rel(2.0))) +
  labs(title="FILTER_100")