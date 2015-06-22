#TCGA.PRE.PROCESS.R
#123114
#Pre-process TCGA.maf file for TRIMER assignment

library(data.table)

args<-commandArgs(trailingOnly=T)

tcga.maf<-args[1]
output.file<-args[2]

TCGA.BRCA<-unique(fread(tcga.maf, header=T, sep="\t", stringsAsFactors=F, drop=c(2:4,7:8,10,12,14:15,17:37))) 
#As in "DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"

setkey(TCGA.BRCA)
TCGA.BRCA<-unique(TCGA.BRCA[Variant_Classification %in% c("Missense_Mutation","Silent"),])

ALL.PATIENTS<-length(unique(as.vector(TCGA.BRCA$Tumor_Sample_Barcode)))

TCGA.BRCA.maf<-TCGA.BRCA[,list(MUT.FREQ=length(Tumor_Seq_Allele2)/ALL.PATIENTS, 
                              Chrom=Chrom, TYPE=Variant_Classification, REF=Reference_Allele, 
                               ALT=Tumor_Seq_Allele2, Tumor_Sample_Barcode=Tumor_Sample_Barcode), 
                         by=c("Start_Position","Hugo_Symbol")]
colnames(TCGA.BRCA.maf)[1]<-"Position"

write.table(file=output.file, TCGA.BRCA.maf, sep="\t", quote=F, row.names=F, col.names=T)
cat ("done pre-processing maf file")


