#031714
#Check mutations across families

#Load BRCA.mi
BRCA.mirna.1<-read.table("DATABASES/CANCER_DATA/miRNA/031714_BRCA.maf.mirna", sep="\t", header=T,stringsAsFactors=F)
head(BRCA.mirna.1)
length(unique(BRCA.mirna.1$Tumor_Sample_Barcode)) #135 tumor samples affected

#Load family information
mirna.fam<-read.table("DATABASES/CANCER_DATA/miRNA/family/miR_Family_Info.txt", header=T, sep="\t",stringsAsFactors=F)
mirna.fam<-mirna.fam[grepl("hsa",mirna.fam$MiRBase.ID),]
head(mirna.fam)
dim(mirna.fam)

#Load miRs mature info
mirna.coor<-read.table("DATABASES/CANCER_DATA/miRNA/hsa._coordinatesgff3.gff3", header=F, sep="\t", skip=13,stringsAsFactors=F)
head(mirna.coor)
mirna.mature<-mirna.coor[mirna.coor$V3=="miRNA",]
mirna.mature$MiRBase.Accession<-sapply(mirna.mature$V9, function(x)  strsplit(strsplit(x,";")[[1]][1],"=")[[1]][2] )
mirna.mature$Accession<-sapply(mirna.mature$V9, function(x)  strsplit(strsplit(x,";")[[1]][4],"=")[[1]][2] )
mirna.mature$mature<-sapply(mirna.mature$V9, function(x)  strsplit(strsplit(x,";")[[1]][3],"=")[[1]][2] )
mirna.mature<-mirna.mature[,c(1,4,5,7,10,11,12)]
dim(mirna.mature)
mirna.mature<-merge(mirna.mature, mirna.fam[,c(1,7)], by="MiRBase.Accession")

#Merge
head(BRCA.mirna.1)
head(mirna.mature)
dim(BRCA.mirna.1)
BRCA.mirna.1<-merge(BRCA.mirna.1, mirna.mature[,c(6,7,8)], by="Accession")

#Plot family frequencies
library(ggplot2)
ggplot(BRCA.mirna.1, aes(x=miR.family)) + geom_histogram(binwidth=1.0,position="identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_text(size=rel(1.5)), plot.title = element_text(size = rel(2))) + 
  annotate("text",label=paste0("Samples=",length(BRCA.mirna.1$miR.family)),x=10,y=10, size=7)+
  labs(title="BRCA", size=7)
  

#031814
#Function to map families to mirna
#Get MI to MIMA file
write.table(file="DATABASES/CANCER_DATA/miRNA/031814_MI_TO_MIMA", mirna.mature[,c(6, 1)], quote=F, row.names=F, sep="\t")

plot.mirna.family<-function (maf.mirna.mut.file, mima.file, fam.file, cancer) {
  #Clean fam file for hsa
  fam<-read.table(fam.file, header=T, sep="\t",stringsAsFactors=F)
  fam<-fam[grepl("hsa",fam$MiRBase.ID),]
  
  #Load maf - of the format 031814_BRCA.maf.mature.mutations
  maf<-read.table(maf.mirna.mut.file, header=T, sep="\t", stringsAsFactors=F)
  
  #Load mima file
  mima<-read.table(mima.file, header=T, sep="\t", stringsAsFactors=F)
  
  #Obtain mima identifiers for maf file
  maf<-merge(maf, mima, by="Accession")
  
  #Assign family based on mima identifier - Keep in mind that not all may have an mirna family
  maf<-merge(maf, fam[,c(7,1)], by="MiRBase.Accession")
  
  #Plot
  ggplot(maf, aes(x=miR.family)) + geom_histogram(binwidth=1.0) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2))) +
    annotate("text",label=paste0("Samples=",length(maf$miR.family)),x=10,y=10, size=7) +
    labs(title=cancer, size=7)  
}

plot.mirna.family("DATABASES/CANCER_DATA/miRNA/031614_UCECmaf.mature.mutations", 
                  "DATABASES/CANCER_DATA/miRNA/031814_MI_TO_MIMA", 
                  "DATABASES/CANCER_DATA/miRNA/family/miR_Family_Info.txt", "UCEC")

#Plot mature vs mirna #Values obtained from running python mirna_maf_mut
mature.vs.mirna.cov<-read.table(col.names=c("Cancer", "Type","Coverage"), text="BLCA 	mature 	4
                                BLCA 	mirna 	9
                                BRCA 	mature 	58
                                BRCA 	mirna 	181
                                CESC 	mature 	4
                                CESC 	mirna 	15
                                COAD 	mature 	7
                                COAD 	mirna 	34
                                HNSC 	mature 	44
                                HNSC 	mirna 	112
                                KICH 	mature 	0
                                KICH 	mirna 	14
                                KIRC 	mature 	10
                                KIRC 	mirna 	35
                                LGG 	mature 	0
                                LGG 	mirna 	5
                                LUAD 	mature 	34
                                LUAD 	mirna 	111
                                LUSC 	mature 	4
                                LUSC 	mirna 	19
                                PAAD 	mature 	8
                                PAAD 	mirna 	29
                                PRAD 	mature 	0
                                PRAD 	mirna 	2
                                READ 	mature 	0
                                READ 	mirna 	3
                                SKCM 	mature 	32
                                SKCM 	mirna 	98
                                THCA 	mature 	0
                                THCA 	mirna 	4
                                UCEC 	mature 	127
                                UCEC 	mirna 	296")
ggplot(mature.vs.mirna.cov, aes(x=Cancer, y=Coverage, fill=Type))+geom_bar(stat="identity", position="dodge") +
  theme(axis.text.y=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.text.x=element_text(size=rel(2.0)), 
        strip.text.x = element_text(size = 30) , legend.text = element_text(size = 22)) +
  labs(title="mature vs mirna", size=7) 
library(reshape)
mature.vs.mirna.cov.casted<-cast(mature.vs.mirna.cov, Cancer~Type )
ggplot(mature.vs.mirna.cov.casted, aes(x=mature, y=mirna, label=Cancer)) + geom_point(size=4) + stat_smooth(method="lm",se=F) +
  theme(axis.text.y=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.text.x=element_text(size=rel(2.0)), 
        strip.text.x = element_text(size = 30) , legend.text = element_text(size = 22)) +
  labs(title="mature vs mirna", size=7) +geom_text(aes(label=Cancer),hjust=0, vjust=0)
lm(mature.vs.mirna.cov.casted$mirna~mature.vs.mirna.cov.casted$mature)