#Map mutations to pri- and mature mirna

#Load BRCA.maf
BRCA.maf<-read.table("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/031714/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t", stringsAsFactors=F)
BRCA.maf<-BRCA.maf[,c(1,5,6,7,11,13,16,17)]
BRCA.maf<-unique(BRCA.maf)
head(BRCA.maf)

#Load pri-miRNA
pri.mirna<-read.table("DATABASES/CANCER_DATA/miRNA/031614_primirna_coor_seq", header=T, sep="\t",stringsAsFactors=F)
colnames(pri.mirna)[2]<-"Hugo_Symbol"
head(pri.mirna)

#Filter BRCA.maf for miRNA 
BRCA.mirna<-BRCA.maf[BRCA.maf$Hugo_Symbol %in% pri.mirna$Hugo_Symbol,]
head(BRCA.mirna)
length(unique(BRCA.mirna$Hugo_Symbol)) #152 miRNA recovered

####Verify -### Done in 031614
length(unique(BRCA.maf$Hugo_Symbol))
TEST<-BRCA.maf[grepl("MIR", BRCA.maf$Hugo_Symbol), ]
TEST<-TEST[!grepl("HG", TEST$Hugo_Symbol),]
head(TEST)
dim(TEST)
setdiff(unique(TEST$Hugo_Symbol), unique(BRCA.mirna$Hugo_Symbol))
#####

#Incorporate strand information to BRCA.mirna
BRCA.mirna<-merge(BRCA.mirna, pri.mirna[,c(2,4,5,6,7,8)], by="Hugo_Symbol")
head(BRCA.mirna)
write.table(file="DATABASES/CANCER_DATA/miRNA/031714_BRCA.maf.mirna", BRCA.mirna, quote=F, row.names=F, sep="\t")

#031814
#Add Hugo Symbol info to mature
mature.coor.seq<-read.table("DATABASES/CANCER_DATA/miRNA/031614_mature_coor_seq", header=T, sep="\t", stringsAsFactors=F)
mature.coor.seq<-merge(mature.coor.seq, pri.mirna[,c(7,2)], by="Accession")
mature.coor.seq<-mature.coor.seq[,c(6,8,2,3,4,5,1,7)]
head(mature.coor.seq)
write.table(file="DATABASES/CANCER_DATA/miRNA/031814_mature_coor_seq.hugo", mature.coor.seq, quote=F, row.names=F, sep="\t")

#Function map pri and miRNA to .maf
map.maf.mirna<-function (maf.file, mirna.file, outfile,type) {
  #Prep maf file
  maf<-read.csv(maf.file, header=T, sep="\t")
  maf<-maf[,c(1,5,6,7,11,13,16,17)]
  maf<-unique(maf)
  
  #Get number of samples
  all.samples<-length(unique(maf$Tumor_Sample_Barcode))
  
  #Load mirna.file (column order of 031614_primirna_coor_seq, make sure it has the hugo_symbol column)
  mirna<-read.table(mirna.file, header=T, sep="\t",stringsAsFactors=F)
  colnames(mirna)[2]<-"Hugo_Symbol"
  
  #Filter maf file for mirna
  maf<-maf[maf$Hugo_Symbol %in% mirna$Hugo_Symbol,]
  
  #Only continue if there are actually mirna in maf file
  mirna.samples<-length(unique(maf$Tumor_Sample_Barcode))
  
  #Output for info
  cat (type, all.samples, mirna.samples,sep="\t")
  cat ("\n")
  
  if (mirna.samples>0) {
    #Incorporate mirna info to maf
    maf<-merge(maf, mirna[,c(2,4,5,6,7,8,1)], by="Hugo_Symbol")
    colnames(maf)[14]<-"Name"
    
    #Write to file
    write.table(file=outfile, maf, quote=F, row.names=F, sep="\t")
  }
}

map.maf.mirna("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/031714/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", "DATABASES/CANCER_DATA/miRNA/031814_mature_coor_seq.hugo", "DATABASES/CANCER_DATA/miRNA/031814_BRCA.maf.mature")

#For pri-mirna
FOLDERS<-list.files(path="DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS")
for (cancer in FOLDERS) {
  file.name<-list.files(path=paste0("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/",cancer),pattern=".maf" )
  map.maf.mirna(paste0("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/",cancer,"/",file.name),
                "DATABASES/CANCER_DATA/miRNA/031614_primirna_coor_seq",
                paste0("DATABASES/CANCER_DATA/miRNA/031614_",cancer,"maf.mature"),
                cancer)
}

mirna.sample.table<-read.table(col.names=c("Cancer","All.samples", "miRNA.samples"),text="ACC  90	0
AML	197	6
BLCA	28	8
BRCA	992	135
CESC	39	11
COAD	219	29
GBM	291	0
HNSC	306	92
KICH	66	13
KIRC	491	43
KIRP	115	1
LGG	220	4
LUAD	538	80
LUSC	178	35
OV	91	0
PAAD	91	13
PRAD	251	2
READ	81	3
SKCM	346	77
STAD	289	0
THCA	405	3
UCEC	248	72")
library(reshape2)
mirna.sample.table.melted<-melt(mirna.sample.table)
ggplot(mirna.sample.table.melted, aes(x=Cancer, y=value, fill=variable)) + geom_bar(stat="identity", position="dodge") + 
  theme(axis.text.y=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.text.x=element_text(size=rel(2.0)), 
        strip.text.x = element_text(size = 30) , legend.text = element_text(size = 22)) +
  labs(title="All samples vs mirna", size=7)  


cor.test(mirna.sample.table$All.samples, mirna.sample.table$miRNA.samples,method="pearson")

#For mature miRNA
for (cancer in FOLDERS) {
  file.name<-list.files(path=paste0("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/",cancer),pattern=".maf" )
  map.maf.mirna(paste0("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/",cancer,"/",file.name),
                "DATABASES/CANCER_DATA/miRNA/031814_mature_coor_seq.hugo",
                paste0("DATABASES/CANCER_DATA/miRNA/031614_",cancer,"maf.mature"),
                cancer)
}
