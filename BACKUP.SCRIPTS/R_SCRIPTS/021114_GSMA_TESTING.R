#FUNCTIONS
GENE.SET.TYPES<-c("GSEA/CANONICAL_PATHWAYS.v4.0.symbols.csv", "GSEA/c3.all.v4.0.symbols.csv" )
MUTATION.DT<-copy(CORE.MUTATIONS[,1:3, with=F])
GENE.SET.BACKGROUND<-"GSEA/msigdb.v4.0.symbols.csv"

GSEA.HYPER<-function (GENE.SET.BACKGROUND, GENE.SET.TYPES, MUTATION.DT, UNIVERSE.CHOICE) {
  library (plyr)
  library (data.table)
  library (qvalue)
  
  #Inner function to compute Gene sets Hypergeometric per each patient using GSEA TABLE of type chosen
  PATIENT.SET.HYPER<-function(PATIENT.DT, GSEA.TYPE.FILE, BACKGROUND.UNIVERSE) {
    
    #Background universe - The set of genes that will be used to determine how many patient genes are in the actual comparisson
    MUTATION.COUNT<-nrow(PATIENT.DT) 
    MUTATED.GENES<-unique(as.vector(PATIENT.DT[["GENE"]]))
    
    PATIENT.GSEA<-as.data.table(t(as.matrix(apply(GSEA.TYPE.FILE, 1, 
                                                  function(AG) c(P.VALUE=phyper(q= length(intersect(as.vector(AG)[as.vector(AG)!=""], MUTATED.GENES))-1,
                                                                                m= length(as.vector(AG)[as.vector(AG)!=""]),
                                                                                n= 45956-length(as.vector(AG)[as.vector(AG)!=""]),
                                                                                k= length(intersect(MUTATED.GENES, BACKGROUND.UNIVERSE )),
                                                                                lower.tail=F),
                                                                 GENE.OVERLAP=length(intersect(as.vector(AG)[as.vector(AG)!=""], MUTATED.GENES)),
                                                                 MUTATION.OVERLAP=nrow(PATIENT.DT[GENE %in% intersect(as.vector(AG)[as.vector(AG)!=""],
                                                                                                                      MUTATED.GENES),]))
                                                  
    )
    )), keep.rownames=T)
    
    #FDR correction
    PATIENT.GSEA$FDR<-p.adjust(PATIENT.GSEA$P.VALUE, method="fdr")
    
    #Q.VALUE correction
    PATIENT.GSEA$Q.VALUE.S<-qvalue(PATIENT.GSEA$P.VALUE, pi0.method="smoother")$qvalues
    PATIENT.GSEA$Q.VALUE.B<-qvalue(PATIENT.GSEA$P.VALUE, pi0.method="bootstrap")$qvalues
    
    #Return
    return (PATIENT.GSEA)
  }
  
  #Process and merge csv files of GSEA types - COMBINED GSEA (Table of GS as row.names and Genes as rows, one gene set per row)
  GSEA.FILES<- lapply(GENE.SET.TYPES, function(g) as.matrix(read.csv(g, header=F, sep=",", fill=T, row.names=1)[,-1]))
  COMBINED.GSEA<-GSEA.FILES[[1]]
  if (length(GSEA.FILES)>1) {
    for (i in 2:length(GSEA.FILES)) {
      COMBINED.GSEA<- as.matrix(rbind.fill(as.data.table(COMBINED.GSEA, keep.rownames=T), as.data.table(GSEA.FILES[[i]], keep.rownames=T)))
      COMBINED.GSEA[is.na(COMBINED.GSEA)]=""
      rownames(COMBINED.GSEA)<-COMBINED.GSEA[,1]
      COMBINED.GSEA<-COMBINED.GSEA[,-1]
    }
  }
  
  #Filter GSEA.types for sets that only contain genes in the patient samples (INTERNAL FILTERING!!!)
  PATIENT.TOTAL.GENES<-unique(as.vector(MUTATION.DT$GENE))
  COMBINED.GSEA<-COMBINED.GSEA[apply(COMBINED.GSEA, 1, function(FG) length(intersect(FG, PATIENT.TOTAL.GENES))!=0),] #FIX!!!???? Try rm() above
  
  #Get GSEA universes (Genes in universes)
  INTERNAL.UNIVERSE<-unique(c(COMBINED.GSEA))
  INTERNAL.UNIVERSE<-INTERNAL.UNIVERSE[INTERNAL.UNIVERSE!=""] #Genes present only in GSEA types provided
  COMPLETE.UNIVERSE<-unique(c(as.matrix(read.csv(GENE.SET.BACKGROUND, header=F, sep=",", fill=T, row.names=1)[,-1])))
  COMPLETE.UNIVERSE<-COMPLETE.UNIVERSE[COMPLETE.UNIVERSE!=""]
  
  #Format mutation file and filter out SILENT mutations - KEEP IN MIND WE ARE REMOVING THEM!
  colnames(MUTATION.DT)<-c("PATIENT", "GENE", "MUT.CLASS")
  MUTATION.DT<-MUTATION.DT[toupper(MUT.CLASS)!="SILENT",]
  
  #Perform test per patient
  PATIENT.LIST<-split(MUTATION.DT, MUTATION.DT$PATIENT)
  PATIENT.GSEA.LISTS<-sapply(names(PATIENT.LIST), function(HG) PATIENT.SET.HYPER(PATIENT.LIST[[HG]], COMBINED.GSEA, COMPLETE.UNIVERSE),
                             simplify=F,USE.NAMES=T)
  
  #Return
  return(PATIENT.GSEA.LISTS)
}

#####################################################################################################################################
BRCA<-read.csv("CANCER_DATA/TCGA/FROM_R_ANALYSIS/012514_BRCA_TOP_CORR_PATIENTS_1_7", sep="\t", header=F, col.names=c("filename")) #47 CORE PATIENTS files
BRCA_IDS<-read.csv("CANCER_DATA/TCGA/EXP_GENE/BRCA/9078bd08-7662-45b9-ab3b-25ddf04f2897/FILE_SAMPLE_MAP.txt",header=T, sep="\t", 
                   col.names=c("filename", "Tumor_Sample_Barcode")) #File identifiers

BRCA.CORE<-merge(BRCA,BRCA_IDS, by="filename")
BRCA.CORE$filename<-NULL
BRCA.CORE$PATIENT<-substr(BRCA.CORE$Tumor_Sample_Barcode, 1,16)

BRCA.MAF<-read.csv("CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/012714_level2.maf", header=T, sep="\t")
BRCA.MAF<-unique(BRCA.MAF[,c(1,9,10,11,13,16)])
BRCA.MAF$PATIENT<-substr(BRCA.MAF$Tumor_Sample_Barcode,1,16)
CORE.MUTATIONS<-unique(merge(BRCA.CORE, BRCA.MAF, by="PATIENT")[,c(1,3,4,5,6,7)]) #46 CORE PATIENTS WITH MUTATIONS INFO!!!

#PLOT_1 - Distribution of mutated genes across core patients
ggplot(as.data.table(CORE.MUTATIONS)[,list(N_GENES=length(unique(Hugo_Symbol))), by=PATIENT], 
       aes(PATIENT, N_GENES)) + geom_histogram() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

#PLOT_2 - Distribution of mutations across core patients
ggplot(as.data.table(CORE.MUTATIONS)[,list(N_MUTATIONS=length(Tumor_Seq_Allele2)), by=PATIENT], 
       aes(PATIENT, N_MUTATIONS)) + geom_histogram() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

CORE.MUTATIONS<-as.data.table(CORE.MUTATIONS)

#No common mutated genes across CORE patients
CORE.GENES<-unique(as.vector(split(CORE.MUTATIONS, CORE.MUTATIONS$PATIENT)[[1]]$Hugo_Symbol))
for (i in 2:length(split(CORE.MUTATIONS, CORE.MUTATIONS$PATIENT))) {
  CORE.GENES<-intersect (CORE.GENES,  unique(as.vector(split(CORE.MUTATIONS, CORE.MUTATIONS$PATIENT)[[i]]$Hugo_Symbol)) )
}
CORE.GENES

#Find common Gene sets across them (FDR Corrected) - Use GSEA.HYPER function with CANONICAL gene sets
GSEA.FILE.TYPES<-c("GSEA/CANONICAL_PATHWAYS.v4.0.symbols.csv", "GSEA/c3.all.v4.0.symbols.csv")
GSEA.FILE.BACKGROUND<-"GSEA/msigdb.v4.0.symbols.csv"
CORE.GSEA.C2P.C3.UPDATED<-GSEA.HYPER(GSEA.FILE.BACKGROUND, GSEA.FILE.TYPES, CORE.MUTATIONS[,1:3, with=F], 3)


CORE.GSEA.C2P.C3.C7[["TCGA-BH-A0DT-01A"]][order(FDR),]
CORE.GSEA.C2P.C7[["TCGA-BH-A0DT-01A"]][order(FDR),]
CORE.GSEA.C2P.C7[[1]][order(FDR),]
CORE.GSEA.C7[["TCGA-BH-A0DT-01A"]][order(FDR),]
CORE.GSEA.C7[[1]][order(FDR),]

CORE.RESULTS<-copy(PATIENT.GSEA.LISTS)
CORE.RESULTS.FREQ<-as.data.table(as.matrix(sapply(CORE.RESULTS, function(CORE) nrow(CORE[Q.VALUE.B<0.05]) , USE.NAMES=T)), keep.rownames=T) 
colnames(CORE.RESULTS.FREQ)<-c("PATIENT", "SETS")
ggplot(CORE.RESULTS.FREQ, aes(x=PATIENT, y=SETS)) + geom_bar(position="identity") + 
  scale_y_log10()+theme(axis.text.x = element_text(angle = 90, hjust = 1))

#CORE.GSEA.CANONICAL
################################################ NOTES ###############################################################################################
GENE_SET<-as.matrix(read.csv("GSEA/CANONICAL_PATHWAYS.v4.0.symbols.csv", header=F, sep=",", fill=T, stringsAsFactors=F, row.names=1)[,-1])
GENE.UNIVERSE<-unique(c(GENE_SET))

TEST.VECTOR<-as.vector(GENE_SET["BIOCARTA_HCMV_PATHWAY",])
TEST.VECTOR<-TEST.VECTOR[TEST.VECTOR!=""]
TEST.HUGO<-as.vector(unique(CORE.MUTATIONS[CORE.MUTATIONS$PATIENT=="TCGA-BH-A0DT-01A",]$Hugo_Symbol))
phyper(q=length(intersect(TEST.HUGO,TEST.VECTOR))-1,m=17   ,n=45956-17 , 
       k=length(intersect(TEST.HUGO,ALL.GENE.UNIVERSE))-1, lower.tail=F) #Probability of drawing more or as extreme as we see (closer to 17 from a very small number of gene we overlap), that is why is we do -1 (P[X>x] vs P[X<=x])

#Find minimum number of sets needed to cover all patients
A<-stack(lapply(CORE.GSEA.C2P.C3, function(SS) as.vector(SS[P.VALUE<0.05,]$rn)))
head(A)
dim(A)
B<-data.frame(table(A$values))
dim(B[B$Freq!=0,])
head(B[order(B$Freq, decreasing=T),])
ggplot(A, aes(x=values)) + geom_histogram() +theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(B, aes(x=Var1, y=Freq)) + geom_bar()

#To PC
save(CORE.GSEA.C2P.C3, file="GSEA/TO_PC/CORE.GSEA.C2P.C3.RData")