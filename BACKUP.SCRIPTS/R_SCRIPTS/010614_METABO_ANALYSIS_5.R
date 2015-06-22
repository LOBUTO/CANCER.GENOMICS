#METABO ANALYSIS #5
#1/6/14
#1. Use table produced in python with mutation count per metabolite to new hypergeowmetric approach
#From 122.py

#Load table
MET_MUT<-read.table("NETWORK/R_ANALYSIS/121813_METABOLITE_TABLE_MUTATIONS", sep="\t", header=T, quote="")
MET_MUT<-as.data.frame(MET_MUT)
head(MET_MUT)

#FILTER OUT THOSE METABOLITES THAT HAVE NO RELATED CANCER MUTATIONS
dim(MET_MUT)
MET_MUT<-MET_MUT[MET_MUT$CANCER_MET_MUTATIONS!=0,]

#Hypergeometric calculation
MET_MUT$P_VALUE<-phyper(MET_MUT$CANCER_MET_MUTATIONS, MET_MUT$TOTAL_MET_MUTATIONS, 
                        MET_MUT$TOTAL_NON_MET_MUTATIONS, MET_MUT$CANCER_NON_MET_MUTATIONS, lower.tail=FALSE)
head(MET_MUT)

#Correct for multiple hypothesis testing
MET_MUT$P_VALUE_CORRECTED<-p.adjust(MET_MUT$P_VALUE, method="fdr")
MET_MUT$SIGNIFICANT<-MET_MUT$P_VALUE_CORRECTED<0.05

#Extract those metabolites that are significant in each cancer type
MET_MUT_SIGNIFICANT<-MET_MUT[MET_MUT$SIGNIFICANT==TRUE, ]
head(MET_MUT_SIGNIFICANT)
dim(MET_MUT_SIGNIFICANT)

#Count how many metabolites per cancer remained
MET_MUT_SIGNIFICANT<-as.data.table(MET_MUT_SIGNIFICANT)
head(MET_MUT_SIGNIFICANT)

SIG_CANCER_MET_COUNT<- MET_MUT_SIGNIFICANT[, list(N_SIG_MET_CANCER=length(METABOLITE)) , by="CANCER"]
SIG_CANCER_MET_COUNT<-SIG_CANCER_MET_COUNT[order(CANCER),]
SIG_CANCER_MET_COUNT$N_TOTAL_MET_CANCER<-c(805, 824, 1222,862 , 1326, 1067,
                                           1232, 810, 1162, 865, 1135, 1336,
                                           1245, 1124, 970, 1115, 1324, 
                                           1344, 854,1222)
SIG_CANCER_MET_COUNT

#PLOT counts of total vs significant metabolites per cancer - SAVE GRAPHS !!!!!!!!
ggplot(melt(SIG_CANCER_MET_COUNT), aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") + 
  theme(axis.text.x=element_text(size=rel(2.0)), axis.text.y=element_text(size=rel(2.0))) + ylab("Number of Metabolites")

ggplot(SIG_CANCER_MET_COUNT, aes(x=N_SIG_MET_CANCER, y=N_TOTAL_MET_CANCER, label=CANCER)) + geom_point() + scale_y_log10() + scale_x_log10() +
  geom_text(colour="black", alpha=0.9, size=8) + theme(axis.text.x=element_text(size=rel(2.0)), axis.text.y=element_text(size=rel(2.0)))

#PEARSON correlation of graph above
cor.test(SIG_CANCER_MET_COUNT$N_SIG_MET_CANCER, SIG_CANCER_MET_COUNT$N_TOTAL_MET_CANCER, method=c("pearson"))

#Write out to file significant metabolites
write.table(MET_MUT_SIGNIFICANT, file="NETWORK/FROM_R_ANALYSIS/010614_SIGNIFICANT_METABOLITES", sep="\t", quote=F, row.names=FALSE)

#Continue analysis in PYTHON.....
