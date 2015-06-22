#METABO ANALYSIS #6
#010714
#1 Check how many patients are lost per cancer
#2. Use significant PMV cancer tables produced in python to display by heatmap
#3. Do PCA analysis
#From 123.py

#Load tables to list
TABLE_FILES<-list.files("NETWORK/R_ANALYSIS", "PMV_TABLE")

PMV_CANCER_TABLES<-list()

for (file in TABLE_FILES) {
    CANCER_NAME<-strsplit(file,"_")[[1]][3]
    
    PMV_CANCER_TABLES[[CANCER_NAME]]<-read.table(paste("NETWORK/R_ANALYSIS/", file,sep="") , sep="\t", quote="",check.names=F,header=T)
    
    rownames(PMV_CANCER_TABLES[[CANCER_NAME]])<-PMV_CANCER_TABLES[[CANCER_NAME]][,1]
    
    PMV_CANCER_TABLES[[CANCER_NAME]]<-as.data.frame(PMV_CANCER_TABLES[[CANCER_NAME]][,-1])    
}

#Patient count per cancer
PMV_PATIENT_COUNT<-data.frame(CANCER=c(), SIG_PATIENTS=c())

for (cancer in names(PMV_CANCER_TABLES)) {
  PMV_PATIENT_COUNT<-rbind(PMV_PATIENT_COUNT, data.frame(CANCER=cancer, SIG_PATIENTS=nrow(PMV_CANCER_TABLES[[cancer]])) )  
}
PMV_PATIENT_COUNT
PMV_PATIENT_COUNT$ALL_PATIENTS<-c(90, 28, 774, 39, 268, 291,
                                  306, 66, 293, 112, 220, 519,
                                  176, 57, 250 , 80, 345, 243 , 403,
                                  194)

#PLOT Patient coverage
ggplot(melt(PMV_PATIENT_COUNT), aes(CANCER, value, fill=variable)) + geom_bar(position="dodge") +
  theme(axis.text.x=element_text(size=rel(2.0)), axis.text.y=element_text(size=rel(2.0))) + ylab("Number of Patients")

ggplot(PMV_PATIENT_COUNT, aes(x=SIG_PATIENTS, y=ALL_PATIENTS, label=CANCER)) + geom_point() + scale_y_continuous() + scale_x_continuous() +
  geom_text(colour="black", alpha=0.9, size=8) + theme(axis.text.x=element_text(size=rel(2.0)), axis.text.y=element_text(size=rel(2.0)))

#Pearson correlation of significant patients vs all patients per cancer
cor.test(PMV_PATIENT_COUNT$SIG_PATIENTS, PMV_PATIENT_COUNT$ALL_PATIENTS, method=c("pearson"))

#HEATMAP cancers
library(RColorBrewer)
hmcol<-brewer.pal(9,"Blues")

for (cancer in names(PMV_CANCER_TABLES)) {
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/HEATMAP/010713_",cancer,"_SIG_METABOLITES.jpeg", sep=""),
       width=1080, height=1080, quality=100, res=100)
  heatmap(as.matrix(PMV_CANCER_TABLES[[cancer]]), col=hmcol, margins=c(16,3),labRow=NA ,ylab="PATIENTS")  
  dev.off()
}


for (cancer in CANCERS.ONLY) {
  
  TOP20CUTOFF<-as.vector(quantile(rowMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]]), c(0.90)))
  
  jpeg(paste("NETWORK/PICTURES/METABO_ANALYSIS/METABO_VECTORS/METABOLITE_SIGNIFICANCE/111913_", cancer, "_SIG_METABOLITES.jpeg", sep=""),
       width=1080, height=1080, quality=100, res=100)
  heatmap(as.matrix(HEATMAP.ENRICH.DATA.LISTS[[cancer]][colMeans(HEATMAP.ENRICH.DATA.LISTS[[cancer]])>TOP20CUTOFF]), col=hmcol, 
          labRow=NA, margins=c(11,11))
  dev.off()
}

citation('Biobase')
