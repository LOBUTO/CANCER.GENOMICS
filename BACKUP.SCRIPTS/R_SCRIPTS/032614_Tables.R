########032614######

library(data.table)
BRCA.table<-read.csv("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header=T, sep="\t")
head(BRCA.table)
dim(BRCA.table)

BRCA.table<-BRCA.table[BRCA.table$Variant_Classification!="Silent",]
BRCA.table$Line_Number<-NULL
BRCA.table<-unique(BRCA.table)
BRCA.table<-BRCA.table[,c(1,16)]
BRCA.table<-as.data.table(BRCA.table)
length(unique(BRCA.table$Tumor_Sample_Barcode))

#Making TABLE 1 for Metabolic network
#Sample - Gene - Number of mutations
#Non-silent mutation only
BRCA.table.1<-copy(BRCA.table)
BRCA.table.1$fill<-1
BRCA.table.1<-BRCA.table.1[,list(N.MUTATIONS=length(fill)), by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]
BRCA.table.1[order(N.MUTATIONS, decreasing=T),]
write.table(file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/032614_BRCA.TCGA_TABLE1", BRCA.table.1, sep="\t", quote=F, row.names=F)

########042114#########

#Load Table 2 - Processed in PC
#Metabolite - Kegg ID - Binding protein (This is the HMDB data that has been cross-referenced against KEGG, meaning that all metabolites have a KEGG ID)
BRCA.table.2<-as.data.table(read.csv("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/032614_BRCA.TCGA_TABLE2", header=T, sep="\t"))
BRCA.table.2

dim(BRCA.table.1); length(unique(BRCA.table.1$Tumor_Sample_Barcode)) #992 samples
dim(BRCA.table.1[Hugo_Symbol %in% BRCA.table.2$GENE,]); length(unique(BRCA.table.1[Hugo_Symbol %in% BRCA.table.2$GENE,]$Tumor_Sample_Barcode)) #976 samples

#Load Table 3 
#Enzyme - Kegg ID - Product
BRCA.table.3<-as.data.table(read.csv("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/032914_ENZYME_PRODUCT", header=T, sep="\t"))
length(unique(BRCA.table.3$KEGG_ID))

#Process Tables
#Have to limit all metabolites and genes in Tables 1 and 2 by the KEGG metabolites found in Table 3
#Leave metabolites in KEGG Table with KEGG_ID only
BRCA.table.3.processed<-copy(BRCA.table.3)
BRCA.table.3.processed<-unique(BRCA.table.3.processed[,c(1,2),with=F]) #DIFFERENT
length(unique(BRCA.table.3.processed$Enzyme)) #2702 enzymes

BRCA.table.2.processed<-copy(BRCA.table.2)
BRCA.table.2.processed<-BRCA.table.2.processed[KEGG_ID %in% BRCA.table.3.processed$KEGG_ID,] #DIFFERENT, need table 3 to do table 2
length(unique(BRCA.table.2.processed$KEGG_ID)) #758 metabolites for network
length(unique(BRCA.table.2.processed$GENE))

BRCA.table.1.processed<-unique(BRCA.table.1[,c(1,2),with=F]) #992 samples
BRCA.table.1.processed<-BRCA.table.1.processed[Hugo_Symbol %in% BRCA.table.2.processed$GENE,] #Samples - Mutated genes found in KEGG and HMDB tables (2 and 3)
length(unique(BRCA.table.1.processed$Tumor_Sample_Barcode)) #974 samples
length(unique(BRCA.table.1.processed$Hugo_Symbol))

length(unique(BRCA.table.1.processed[Hugo_Symbol %in% BRCA.table.3.processed$Enzyme,]$Tumor_Sample_Barcode)) #954 samples that have at least 1 KEGG enzyme

#Save processed Tables
save(list=c("BRCA.table.3.processed", "BRCA.table.1.processed", "BRCA.table.2.processed"),file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042214_MTables.RData")

######042214#########

#Plot Tables
Table.3.plot<-BRCA.table.3.processed[,list(Enzymes=length(unique(Enzyme))), by ="KEGG_ID"]
Table.3.plot<-Table.3.plot[order(Enzymes, decreasing=T),]
Table.3.plot$KEGG_ID<-factor(x=Table.3.plot$KEGG_ID,levels=Table.3.plot$KEGG_ID,ordered=T)
ggplot(Table.3.plot, aes(x=KEGG_ID, y=Enzymes)) + geom_histogram(stat="identity",binwidth=5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="KEGG Dataset") + xlab("Metabolites") + ylab("Enzymes") +
  opts(axis.text.x = theme_blank())

Table.2.plot<-BRCA.table.2.processed[,list(Proteins=length(unique(GENE))), by="METABOLITE"]
Table.2.plot<-Table.2.plot[order(Proteins, decreasing=T),]
Table.2.plot$METABOLITE<-factor(x=Table.2.plot$METABOLITE,levels=Table.2.plot$METABOLITE, ordered=T )
ggplot(Table.2.plot, aes(x=METABOLITE, y=Proteins)) + geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

Table.1.plot<-BRCA.table.1.processed[,list(Genes=length(unique(Hugo_Symbol))), by="Tumor_Sample_Barcode"]
Table.1.plot<-Table.1.plot[order(Genes, decreasing=T),]
Table.1.plot$Tumor_Sample_Barcode<-factor(x=Table.1.plot$Tumor_Sample_Barcode,levels=Table.1.plot$Tumor_Sample_Barcode,ordered=T)
ggplot(Table.1.plot, aes(x=Tumor_Sample_Barcode, y=Genes)) + geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="TCGA-BRCA Dataset") + xlab("Patient Samples") + ylab("Number of Mutated Genes")
