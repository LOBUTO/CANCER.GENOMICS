library(data.table)
library(Hmisc)
library(pheatmap)
library(gplots)

PATIENT.COVERAGE<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PATIENT.COVERAGE.rds")
Table.1.PLUS<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.Table.1.PLUS.rds")

BRCA.NORM.MATRICES<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082614.CANCER.MATRICES.NORMALIZED.OBJ.rds")
head(BRCA.NORM.MATRICES$combined.matrices[,1:3])

##
earlier.test<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES")
##

TP53.PATIENTS<-as.vector(Table.1.PLUS[Hugo_Symbol=="TP53",]$PATIENT)
TTN.PATIENTS<-as.vector(Table.1.PLUS[Hugo_Symbol=="TTN",]$PATIENT)
CCND1.PATIENTS<-as.vector(Table.1.PLUS[Hugo_Symbol=="CCND1"]$PATIENT)
PIK3CA.PATIENTS<-as.vector(Table.1.PLUS[Hugo_Symbol=="PIK3CA"]$PATIENT)
TP53.TTN.PATIENTS<-intersect(TP53.PATIENTS, TTN.PATIENTS)
ONLY.TP53.PATIENTS<-TP53.PATIENTS[!(TP53.PATIENTS %in% TP53.TTN.PATIENTS)]
ONLY.TTN.PATIENTS<-TTN.PATIENTS[!(TTN.PATIENTS %in% TP53.TTN.PATIENTS)]

##
MATRIX.TP53<-BRCA.NORM.MATRICES$combined.matrices[,intersect(colnames(BRCA.NORM.MATRICES$combined.matrices), ONLY.TP53.PATIENTS)]
MATRIX.TTN<-BRCA.NORM.MATRICES$combined.matrices[,intersect(colnames(BRCA.NORM.MATRICES$combined.matrices), ONLY.TTN.PATIENTS)]
MATRIX.CCND1<-BRCA.NORM.MATRICES$combined.matrices[,intersect(colnames(BRCA.NORM.MATRICES$combined.matrices), CCND1.PATIENTS)]
MATRIX.BOTH<-BRCA.NORM.MATRICES$combined.matrices[,intersect(colnames(BRCA.NORM.MATRICES$combined.matrices), TP53.TTN.PATIENTS)]
NOPE.PATIENTS<-colnames(BRCA.NORM.MATRICES$combined.matrices)[!(colnames(BRCA.NORM.MATRICES$combined.matrices) %in% c(TP53.PATIENTS, TTN.PATIENTS, CCND1.PATIENTS))]
MATRIX.NOPE<-BRCA.NORM.MATRICES$combined.matrices[,NOPE.PATIENTS]

##
MATRIX.TP53<-earlier.test$tumor[,intersect(colnames(earlier.test$tumor), ONLY.TP53.PATIENTS)]
MATRIX.TTN<-earlier.test$tumor[,intersect(colnames(earlier.test$tumor), ONLY.TTN.PATIENTS)]
MATRIX.CCND1<-earlier.test$tumor[,intersect(colnames(earlier.test$tumor), CCND1.PATIENTS)]
MATRIX.PIK3CA<-earlier.test$tumor[,intersect(colnames(earlier.test$tumor), PIK3CA.PATIENTS)]
MATRIX.BOTH<-earlier.test$tumor[,intersect(colnames(earlier.test$tumor), TP53.TTN.PATIENTS)]
NOPE.PATIENTS<-colnames(earlier.test$tumor)[!(colnames(earlier.test$tumor) %in% c(TP53.PATIENTS, TTN.PATIENTS,CCND1.PATIENTS)) ]
MATRIX.NOPE<-earlier.test$tumor[,NOPE.PATIENTS]
##

colnames(MATRIX.TTN)<- as.vector(sapply(1:length(colnames(MATRIX.TTN)), function(x) paste0("TTN.",x)))
colnames(MATRIX.TP53)<- as.vector(sapply(colnames(MATRIX.TP53), function(x) paste0("TP53.",x)))
colnames(MATRIX.CCND1)<-as.vector(sapply(1:length(colnames(MATRIX.CCND1)), function(x) paste0("CCND1.",x)))
colnames(MATRIX.PIK3CA)<-as.vector(sapply(1:length(colnames(MATRIX.PIK3CA)), function(x) paste0("PIK3CA.",x)))
colnames(MATRIX.BOTH)<- as.vector(sapply(1:length(colnames(MATRIX.BOTH)), function(x) paste0("BOTH!!.",x)))
colnames(MATRIX.NOPE)<- as.vector(sapply(1:length(colnames(MATRIX.NOPE)), function(x) paste0("--RANDOM--.",x)))

WHOLE.MATRIX<-cbind(MATRIX.TP53[,sample(colnames(MATRIX.TP53),40)], MATRIX.TTN[,sample(colnames(MATRIX.TTN),40)], MATRIX.BOTH[,sample(colnames(MATRIX.BOTH),10)])
WHOLE.MATRIX<-cbind(MATRIX.TP53[,sample(colnames(MATRIX.TP53),50)], MATRIX.TTN[,sample(colnames(MATRIX.TTN),50)], MATRIX.NOPE[,sample(colnames(MATRIX.NOPE),20)])
WHOLE.MATRIX<-cbind(MATRIX.TP53[,sample(colnames(MATRIX.TP53),50)],  MATRIX.NOPE[,sample(colnames(MATRIX.NOPE),50)])
WHOLE.MATRIX<-cbind(MATRIX.TTN[,sample(colnames(MATRIX.TTN),50)],  MATRIX.NOPE[,sample(colnames(MATRIX.NOPE),50)])
WHOLE.MATRIX<-cbind(MATRIX.CCND1[,sample(colnames(MATRIX.CCND1),50)],  MATRIX.NOPE[,sample(colnames(MATRIX.NOPE),50)] )
WHOLE.MATRIX<-cbind(MATRIX.PIK3CA[,sample(colnames(MATRIX.PIK3CA),100)],  MATRIX.NOPE[,sample(colnames(MATRIX.NOPE),100)] )
head(WHOLE.MATRIX[,1:3])
WHOLE.MATRIX.CORR<-rcorr(as.matrix(WHOLE.MATRIX), type="spearman")
head(WHOLE.MATRIX.CORR$r[,1:3])

heatmap(WHOLE.MATRIX.CORR$r, scale="none")
pheatmap(WHOLE.MATRIX.CORR$r, scale="none")
heatmap.2(WHOLE.MATRIX.CORR$r, scale="none",trace = "none")

#K-means clustering
TP53.RANDOM<-kmeans(t(WHOLE.MATRIX),3,nstart=100)
as.data.table(TP53.RANDOM$cluster,keep.rownames=T)
