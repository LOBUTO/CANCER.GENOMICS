#METABO ANALYSIS #3
#11/19/13 - FROM 154.py
#Analysis of enrichment of metabolites across patients per cancer

library(data.table)
library(ggplot2)

#Phyper table per patient in cancers
ENRICH.DATA<-read.table("NETWORK/R_ANALYSIS/111913_PATIENT_ENRICHMENT_PREFDR", header=TRUE, sep="\t",quote="")
ENRICH.DATA<-as.data.table(ENRICH.DATA)
head(ENRICH.DATA,15)

#Patients per cancer that show enrichment
ENRICH.DATA[, list(N_PATIENTS=length(unique(PATIENT))), by="CANCER"]

#FDR pvalues
ENRICH.DATA$FDR_CORRECTED<-p.adjust(ENRICH.DATA$P_VALUE, method="fdr")
head(ENRICH.DATA)
nrow(ENRICH.DATA)
ENRICH.DATA[ENRICH.DATA$P_VALUE==0.0475886953026,]

#Filter for FDR<0.05
ENRICH.DATA.FILTERED<-ENRICH.DATA[ENRICH.DATA$FDR_CORRECTED<0.05,]
nrow(ENRICH.DATA.FILTERED)
head(ENRICH.DATA.FILTERED)
ENRICH.DATA.FILTERED[,list(N_PATIENTS=length(unique(PATIENT))), by="CANCER"] #Number of patients that are significant

#COUNT TOTAL NUMBER OF UNIQUE METABOLITES PER PATIENT VS SIGNIFICANT NUMBER OF METABOLITES PER PATIENT
ENRICH.COUNT.DATA<-read.table("NETWORK/R_ANALYSIS/111913_PATIENT_ENRICHMENT_COUNT", header=TRUE, sep="\t")
ENRICH.COUNT.DATA<-ENRICH.COUNT.DATA[order(ENRICH.COUNT.DATA$CANCER, -ENRICH.COUNT.DATA$TOTAL_UNIQUE_METABOLITES),]
head(ENRICH.COUNT.DATA)

#Get average percentage of metabolites lost per patient
ENRICH.COUNT.DATA.TABLE<-as.data.table(ENRICH.COUNT.DATA)
ENRICHED<-ENRICH.COUNT.DATA.TABLE[,list(PER.REMAIN=sum(SIGNIFICANT_METABOLITES)/sum(TOTAL_UNIQUE_METABOLITES)), by="CANCER"]

ENRICHED<-dt[, list(mean=mean(price)), by="cut"]

#SAVE PLOTS for coverage
for (cancer in CANCERS.ONLY) {
  TYPE1<-ENRICH.COUNT.DATA[ENRICH.COUNT.DATA$CANCER==cancer,]
  TYPE1$COUNT<-1:length(TYPE1$PATIENT)
  
  TYPE1.PLOT<-ggplot(TYPE1, aes(x=COUNT)) +geom_line(aes(y=TOTAL_UNIQUE_METABOLITES, colour="TOTAL_UNIQUE_METABOLITES")) +
    geom_line(aes(y=SIGNIFICANT_METABOLITES, colour="SIGNIFICANT_METABOLITES")) + labs(title=cancer) + 
    theme(legend.position="bottom", legend.text=element_text(size=12)) + xlab(label="PATIENTS")
  ggsave(TYPE1.PLOT, 
         file=paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/111913_", cancer, "MET_SIG_COVERAGE.jpeg", sep=""),
         dpi=600, scale=2)
}

#HEATMAP DATA PATIENT X METABOLITE
head(ENRICH.DATA.FILTERED)
HEATMAP.ENRICH.DATA<-ENRICH.DATA.FILTERED[,c(1:4),with=FALSE]
head(HEATMAP.ENRICH.DATA)

#Build heatmap including all metabolites in all patients per cancer
HEATMAP.ENRICH.DATA.LISTS=list()
for (cancer in CANCERS.ONLY) {
  
  #Separate by cancer
  TYPE2<-HEATMAP.ENRICH.DATA[CANCER==cancer,]
  
  #Get all metabolites per cancer
  CANCER_METABOLITES<-as.vector(unique(TYPE2$METABOLITE))
  
  #Get all patients
  CANCER_PATIENTS<-as.vector(unique(TYPE2$PATIENT))
  
  #Start patient data frame with empty rows per metabolite column
  HEATMAP.PATIENT.DF=data.frame(matrix(nrow=length(CANCER_PATIENTS), ncol=length(CANCER_METABOLITES)))
  rownames(HEATMAP.PATIENT.DF)<-CANCER_PATIENTS
  colnames(HEATMAP.PATIENT.DF)<-CANCER_METABOLITES
  
  #Initialize with zeroes
  HEATMAP.PATIENT.DF[is.na(HEATMAP.PATIENT.DF)]<-0
  
  #Fill with values
  for (patient in CANCER_PATIENTS) {
    
    #Get patient metabolites
    PATIENT_METABOLITES<-as.vector(TYPE2[TYPE2$PATIENT==patient,]$METABOLITE)
    
    #Fill with metabolite counts #HEATMAP NEED NAME OF PATIENTS!!!
    for (met in PATIENT_METABOLITES) {
      HEATMAP.PATIENT.DF[patient,met]<-TYPE2[PATIENT==patient & METABOLITE==met,]$METABOLITE_COUNT
    }      
  }
  HEATMAP.ENRICH.DATA.LISTS[[cancer]]<-HEATMAP.PATIENT.DF
}

HEATMAP.ENRICH.DATA.LISTS$LUSC[c(1:5), c(1:10)]

#PLOT HEATMAP - PICK TOP 20%
library(RColorBrewer)
hmcol<-brewer.pal(9,"Blues")

for (cancer in CANCERS.ONLY) {

  TOP20CUTOFF<-as.vector(quantile(rowMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]]), c(0.90)))
  
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/111913_", cancer, "_SIG_METABOLITES.jpeg", sep=""),
       width=1080, height=1080, quality=100, res=100)
  heatmap(as.matrix(HEATMAP.ENRICH.DATA.LISTS[[cancer]][colMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]])>TOP20CUTOFF]), col=hmcol, 
          labRow=NA, margins=c(11,11))
  dev.off()
}

## LEARNING PCA
#Example 1
DAT=read.table("http://astrostatistics.psu.edu/su09/lecturenotes/marks.dat", head=T)
head(DAT)
plot(DAT$Phys, DAT$Stat)
PC=princomp(DAT,)
PC$loadings #Components order in terms of which gives maximum spread first (Comp1), then Comp2 gives maximum spread perpendicular to it and so on
            #These are unit vectors
PC #Spread along each direction (magnitudes)
plot(PC)
names(PC)
PC$scores #Projections of data points in principal component frame

#Example 2
QUASAR=read.table("http://astrostatistics.psu.edu/su09/lecturenotes/SDSS_quasar.dat", head=T)
head(QUASAR)
dim(QUASAR)
QUASAR=na.omit(QUASAR)
PC2=princomp(QUASAR[,-1], scores=T)
PC2
plot(PC2, type="lines") #As seen most of the bulk (variation) is seen from 2 components

#Look for the place which has this 2 planes
M<-PC2$loadings #These are 2 PC that account for bulk
t(M) %*% M #This should be perpendicular, check, should produce Identity Matrix

head(PC2$scores)
#Project data onto plane
plot(PC2$scores[,1], PC2$scores[,2], pch=".")

#PCA with function prcomp()
head(USArrests)
PCA1=prcomp(USArrests, scale.=TRUE,)
PCA1$rotation #Loadings (directions of unit vectors of PC)
head(PCA1$x) #PCs (aka scores) #Projection of data along principal components
plot(PCA1)
#PLOT so it makes sense
library(ggplot2)
PCA.SCORES<-as.data.frame(PCA1$x)
head(PCA.SCORES)
ggplot(PCA.SCORES, aes(x=PC1, y=PC2, label=rownames(PCA.SCORES))) +
  geom_point() + geom_text(colour="tomato", alpha=0.8, size=4) #WOULD JUST TELL ME HOW PATIENTS CLUSTER, NOT UNDER WHAT CONDITIONS

correlations = as.data.frame(cor(USArrests, PCA1$x))
correlations
cor(USArrests$Murder, PCA1$x[,2]) #Correlation between PCA and original vector, if it was causal for shown variance, then it will show correlation
library(FactoMineR)

#PLOT correlations of actual variables
library(corrplot)
corrplot(as.matrix(correlations, order="hclust"))

##TRY PCA and CORR with BRCA
TOP20CUTOFF.BRCA<-as.vector(quantile(rowMeans(HEATMAP.ENRICH.DATA.LISTS$BRCA), c(0.90)))
BRCA.TOP20<-as.matrix(HEATMAP.ENRICH.DATA.LISTS$BRCA[colMeans(HEATMAP.ENRICH.DATA.LISTS$BRCA)>TOP20CUTOFF.BRCA])

PCA.BRCA<-prcomp(BRCA.TOP20)#, scale.=TRUE)
plot(PCA.BRCA) #Use TOP 3 PCAs
PCA.BRCA.SCORES<-as.data.frame(PCA.BRCA$x)

#PLOT patients by first 2 PCAs
head(PCA.BRCA.SCORES)
ggplot(PCA.BRCA.SCORES, aes(x=PC1, y=PC2, label=rownames(PCA.BRCA.SCORES)))+
  geom_point() + geom_text(colour="tomato", alpha=0.8, size=4)

#PLOT variables (metabolites) as correlations matrices by first 3 PCAs
BRCA.CORRELATIONS<-as.data.frame(cor(BRCA.TOP20, PCA.BRCA$x))
dim(BRCA.CORRELATIONS)
BRCA.CORRELATIONS[1:2,1:2]
corrplot(as.matrix(BRCA.CORRELATIONS[120:177,1:6], order="hclust"),tl.cex=0.5,cl.pos="b") 

#PLOT AS ABOVE (CIRCLE PLOT)
# function to create a circle
circle <- function(center=c(0,0), npoints=100)
{
  r = 1
  tt = seq(0, 2*pi, length=npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0,0), npoints = 100)

# data frame with arrows coordinates
BRCA.CORRELATIONS.FILTERED<-BRCA.CORRELATIONS[,1:2]
BRCA.CORRELATIONS.FILTERED<-BRCA.CORRELATIONS.FILTERED[sqrt(BRCA.CORRELATIONS.FILTERED$PC1^2 + BRCA.CORRELATIONS.FILTERED$PC2^2)>0.3, ]
dim(BRCA.CORRELATIONS.FILTERED)

arrows = data.frame(x1=rep(0,nrow(BRCA.CORRELATIONS.FILTERED)), y1=rep(0,nrow(BRCA.CORRELATIONS.FILTERED)),
                    x2=BRCA.CORRELATIONS.FILTERED$PC1, y2=BRCA.CORRELATIONS.FILTERED$PC2)

# geom_path will do open circles
ggplot() +
  geom_path(data=corcir, aes(x=x, y=y), colour="gray65") +
  geom_segment(data=arrows, aes(x=x1, y=y1, xend=x2, yend=y2), colour="gray65") +
  geom_text(data=BRCA.CORRELATIONS.FILTERED, aes(x=PC1, y=PC2, label=rownames(BRCA.CORRELATIONS.FILTERED))) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  xlim(-1.1,1.1) + ylim(-1.1,1.1) +
  labs(x="pc1 aixs", y="pc2 axis") #+opts(title="Circle of correlations")

##APPLY PCA TO ALL CANCERS - APPLIED TO TOP 15%
# function to create a circle -NEED FOR CIR PLOT
circle <- function(center=c(0,0), npoints=100)
{
  r = 1
  tt = seq(0, 2*pi, length=npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

for (cancer in CANCERS.ONLY) {
  
  #Get top 15%
  TOP15CUTOFF.CANCER<-as.vector(quantile(rowMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]]), c(0.90)))
  CANCER.TOP15<-as.matrix(HEATMAP.ENRICH.DATA.LISTS[[cancer]][colMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]])>TOP15CUTOFF.CANCER]) 
  
  #DO PCA
  PCA.CANCER<-prcomp(CANCER.TOP15)
  PCA.CANCER.SCALED<-prcomp(CANCER.TOP15, scale.=TRUE) #FOR COMPARISSON IN PCA.PLOT!
  
  #PLOT PCA
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/PCA/PCAS/111913_",cancer,"_PCAS.jpeg",sep=""),
       width=1080, height=1080, quality=100, res=100)
  plot(PCA.CANCER)
  dev.off()
  
  #PLOT Patient clustering by PC1 and PC2 using SCORES ($x, score of patient, mapping of patient into PCA frames)
  PCA.CANCER.SCORES<-as.data.frame(PCA.CANCER$x)
  PCA.PLOT<-ggplot(PCA.CANCER.SCORES, aes(x=PC1, y=PC2, label=rownames(PCA.CANCER.SCORES))) +
    geom_point() + geom_text(colour="tomato", alpha=0.8, size=4)
  ggsave(PCA.PLOT,
         file=paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/PCA/PATIENT_BY_2DPCA/111913_",
                    cancer,"_PATIENT_PCA.jpeg", sep=""),
         dpi=600, scale=2) 
  
  #SCALED VERSION OF ABOVE
  PCA.CANCER.SCORES.SCALED<-as.data.frame(PCA.CANCER.SCALED$x)
  PCA.PLOT.SCALED<-ggplot(PCA.CANCER.SCORES.SCALED, aes(x=PC1, y=PC2, label=rownames(PCA.CANCER.SCORES.SCALED))) +
    geom_point() + geom_text(colour="tomato", alpha=0.8, size=4)
  ggsave(PCA.PLOT.SCALED,
         file=paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/PCA/PATIENT_BY_2DPCA/111913_",
                    cancer,"_PATIENT_PCA_SCALED.jpeg", sep=""),
         dpi=600, scale=2)    
  
  #PLOT variable(metabolites) as correlations matrices by !!FIRST 6 PCAs!! (CORRELATION OF METABO COLUMN WITH PCA VECTOR) - SAMPLE OF 1:40 METABO
  CANCER.CORRELATIONS<-as.data.frame(cor(CANCER.TOP15, PCA.CANCER$x))
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/PCA/METABO_PCA_CORR/111913_", cancer,"_SAMPLE40.jpeg", sep=""),
       width=1080,height=1080, quality=100, res=100)
  corrplot(as.matrix(CANCER.CORRELATIONS[1:40,1:6], order="hclust"), tl.cex=0.5, cl.pos="b")
  dev.off()
  
  #PLOT ABOVE as CIRCLE PLOT - Use PC1 and PC2 only at CUTOFF FOR PCA VECTOR AT >0.3 MAGNITUDE
  corcir = circle(c(0,0), npoints = 100)
  
  #data frame with arrow coordinates
  CANCER.CORRELATIONS.FILTERED<-CANCER.CORRELATIONS[,1:2] #USING PC1 and PC2 ONLY
  CANCER.CORRELATIONS.FILTERED<-CANCER.CORRELATIONS.FILTERED[sqrt(CANCER.CORRELATIONS.FILTERED$PC1^2 + CANCER.CORRELATIONS.FILTERED$PC2^2)>0.3, ]
  
  arrows = data.frame(x1=rep(0,nrow(CANCER.CORRELATIONS.FILTERED)), y1=rep(0,nrow(CANCER.CORRELATIONS.FILTERED)),
                      x2=CANCER.CORRELATIONS.FILTERED$PC1, y2=CANCER.CORRELATIONS.FILTERED$PC2)
  
  #geom_path will do open circles
  CIRCLE.COOR<-ggplot() +
    geom_path(data=corcir, aes(x=x, y=y), colour="gray65") +
    geom_segment(data=arrows, aes(x=x1, y=y1, xend=x2, yend=y2), colour="gray65") +
    geom_text(data=CANCER.CORRELATIONS.FILTERED, aes(x=PC1, y=PC2, label=rownames(CANCER.CORRELATIONS.FILTERED))) +
    geom_hline(yintercept=0, colour="gray65") +
    geom_vline(xintercept=0, colour="gray65") +
    xlim(-1.1,1.1) + ylim(-1.1,1.1) +
    labs(x="pc1 aixs", y="pc2 axis", title=cancer)
  
  ggsave(CIRCLE.COOR,
         file=paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/PCA/METABO_PCA_CORR_CIR/111913_",
                    cancer,"_CIR_2DPCA.jpeg", sep=""),
         dpi=600, scale=2)
}

##K-MEAN CLUSTERING ON ABOVE PCAS (UNSCALED)

#Try for BRCA.TOP20%
TOP20CUTOFF.BRCA<-as.vector(quantile(rowMeans(HEATMAP.ENRICH.DATA.LISTS$BRCA), c(0.90)))
BRCA.TOP20<-as.matrix(HEATMAP.ENRICH.DATA.LISTS$BRCA[colMeans(HEATMAP.ENRICH.DATA.LISTS$BRCA)>TOP20CUTOFF.BRCA])
head(BRCA.TOP20)

PCA.BRCA<-prcomp(BRCA.TOP20)
plot(PCA.BRCA) #Use TOP 5 PCAs
PCA.BRCA.SCORES<-as.data.frame(PCA.BRCA$x) #get coordinate of patients for all PCAs
dim(PCA.BRCA.SCORES)
PCA.BRCA.SCORES[1:3,1:3]

#IF K-MEANS CLUSTERING IS APPLIED TO PCA.BRCA.SCORES AT THIS POINT USING 5-6 PCAs I WOULD PROBABLY GET PATIENT CLUSTERS SIMILAR TO SPACES DIFFERENCES SHOWN BY HEATMAP 

#FOR NOW FOCUS ON METABOLITES CLUSTERING - For this, I need correlations so I have in PCA IN VARIABLE TERMS
BRCA.CORRELATIONS<-as.data.frame(cor(BRCA.TOP20, PCA.BRCA$x))
dim(BRCA.CORRELATIONS)

BRCA.CORRELATIONS.FILTERED<-BRCA.CORRELATIONS[,1:6]
head(BRCA.CORRELATIONS.FILTERED)
#Filtering by magnitude
BRCA.CORRELATIONS.FILTERED<-BRCA.CORRELATIONS.FILTERED[sqrt(BRCA.CORRELATIONS.FILTERED$PC1^2 + BRCA.CORRELATIONS.FILTERED$PC2^2 +
                                                            BRCA.CORRELATIONS.FILTERED$PC3^2 + BRCA.CORRELATIONS.FILTERED$PC4^2 +
                                                            BRCA.CORRELATIONS.FILTERED$PC5^2 + BRCA.CORRELATIONS.FILTERED$PC6^2  )>0.3, ]

#K-MEAN clustering using 6 PCAs
KM.BRCA<-kmeans(BRCA.CORRELATIONS.FILTERED, 6, nstart=100)
KM.BRCA$size

KM.BRCA$cluster
KM.BRCA$centers
plot(BRCA.CORRELATIONS.FILTERED, col=KM.BRCA$cluster)

BRCA.KM.CLUSTERS<-as.data.frame(KM.BRCA$cluster)
BRCA.KM.CLUSTERS$METABOLITE<-rownames(BRCA.KM.CLUSTERS)
colnames(BRCA.KM.CLUSTERS)<-c("CLUSTER", "METABOLITES")
View(BRCA.KM.CLUSTERS[order(BRCA.KM.CLUSTERS$CLUSTER),])

dim(BRCA.TOP20)

TEST.HM<-heatmap.2(BRCA.TOP20[, colSums()] , col=hmcol)
HC<-as.hclust(TEST.HM$roDendrogram)


##
head(BRCA.TOP20)
d=dist(t(BRCA.TOP20)) #Need a distance matrix
h=hclust(d) #clusters distance matrix
plot(h) #We can see the dendogram

D.C<-as.data.frame(cutree(h,50)) #Cuts dendogram into pre-determined clusters
D.C$MET<-rownames(D.C)
colnames(D.C)<-c("CLUSTER", "MET")
head(D.C)

D.C13<-D.C[D.C$CLUSTER==13,]
D.C13
write.table(D.C13[,2], file="NETWORK/FROM_R_ANALYSIS/112013_BRCA_CLUSTER13",sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  