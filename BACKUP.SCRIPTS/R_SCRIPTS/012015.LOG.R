####Analsyis of enrichment scores for cancer metabolites for grant####
#012015
library(data.table)
library(pheatmap)
library(gplots)

###FUNCTIONS###
Function.Process.MAF<-function(maf.file){
  #Function to clean maf files straight from TCGA
  
  maf.record<-fread(maf.file, sep="\t", header=T, stringsAsFactors=F, drop=c(2:4,7:8, 10,12,14:15,17:37))
  maf.record<-maf.record[Variant_Classification %in% c("Missense_Mutation","Silent"),]
  maf.record<-maf.record[nchar(Tumor_Seq_Allele2)==1 & nchar(Reference_Allele)==1,]
  setkey(maf.record)
  maf.record<-unique(maf.record)
  
  #Return
  return(maf.record)
}

###DOCUMENTS###
#Cancer files
AML="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/AML/Somatic_Mutations/WUSM__IlluminaGA_DNASeq/Level_2/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
BRCA="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/100514/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
COAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/COAD/100914/Somatic_Mutations/BCM__IlluminaGA_DNASeq/Level_2/hgsc.bcm.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
GBM="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/GBM/100814/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
LUAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/LUAD/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
LUSC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/LUSC/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
HNSC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/HNSC/8998ada8-b0cb-4cf2-bdd0-7321f8904f13/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
KIRC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/KIRC/062014/Somatic_Mutations/BI__IlluminaGA_DNASeq_automated/Level_2/broad.mit.edu__IlluminaGA_automated_DNA_sequencing_level2.maf"
OV="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/OV/101014/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
PRAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/PRAD/970b92f0-04d7-4021-980c-9a0c32ecadf3/Somatic_Mutations/BI__IlluminaGA_DNASeq_curated/Level_2/broad.mit.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
READ="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/READ/101014/Somatic_Mutations/BCM__SOLiD_DNASeq/Level_2/hgsc.bcm.edu__ABI_SOLiD_DNA_Sequencing_level2.maf"
SKCM="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/SKCM/608af7c2-400c-4102-bb99-eaced152f04f/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"
STAD="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/STAD/broad.mit.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
UCEC="/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/UCEC/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"

###ANALYSIS###
#Process metabolites 
met.enriched<-data.table()
cancers=c("AML","GBM","LUSC","SKCM","HNSC","OV","STAD","BRCA","KIRC","PRAD","UCEC","COAD","LUAD","READ")
for (cancer in cancers){
  cancer.in<-fread(paste0("PIPELINES/METABOLIC.DRIVERS/TABLES/",cancer,"/012015.",cancer,".MET.HUGO.ENRICHMENT"), header=T, sep="\t",stringsAsFactors=F)
  cancer.in$CANCER<-cancer
  
  #Filter for significant scores
  cancer.in<-cancer.in[P.VAL.ADJ<0.05,]
  
  if (nrow(cancer.in)>0){
    met.enriched<-rbind(met.enriched, cancer.in)
  }
}

#Find common metabolites
met.enriched[,MET.COUNT:=length(CANCER), by=METABOLITE]
met.enriched[MET.COUNT>1,][order(METABOLITE),]

#Find common genes using binding info in table.2
table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
met.enriched<-merge(met.enriched, table.2, by=c("METABOLITE","KEGG_ID"))

#Filter for those genes actually found mutated in the cancer
met.enriched.filtered<-list()
for (cancer in unique(met.enriched$CANCER)){
  cancer.maf<-Function.Process.MAF(get(cancer))
  cancer.maf<-cancer.maf[Variant_Classification=="Missense_Mutation",]
  cancer.genes<-unique(as.vector(cancer.maf$Hugo_Symbol))
  met.enriched.filtered[[cancer]]<- met.enriched[CANCER==cancer & GENE %in% cancer.genes,]
}
met.enriched.filtered<-do.call(rbind, met.enriched.filtered)

#Count genes across cancer
met.enriched.filtered[,GENE.COUNT:=length(unique(CANCER)), by=GENE]

unique(met.enriched.filtered[GENE.COUNT>=3,][order(GENE),][,c("CANCER","GENE","GENE.COUNT"),with=F])
met.enriched.filtered[CANCER=="LUAD",]

#Gene count per cancer
met.enriched.filtered[,list(N.GENES=length(unique(GENE))), by=CANCER]

saveRDS(met.enriched.filtered,file="PIPELINES/METABOLIC.DRIVERS/OBJECTS/012115.MET.ENRICHED.FILTERED.rds")

#What do patients with metabolite associated mutations have in common? - Look at breast cancer  mutations - NOTHING APPARENT IN TERMS OF MUTATIONS
maf.record
met.enriched.filtered[CANCER=="BRCA",]
maf.record.brca.met<-maf.record[PATIENT %in% unique(maf.record[Hugo_Symbol %in% met.enriched.filtered[CANCER=="BRCA",]$GENE, ]$PATIENT),]
maf.record.brca.met$count<-1
maf.record.brca.met<-as.data.frame(dcast(maf.record.brca.met, PATIENT~Hugo_Symbol,sum,fill=0))
rownames(maf.record.brca.met)<-maf.record.brca.met$PATIENT
maf.record.brca.met$PATIENT<-NULL
maf.record.brca.met<-as.matrix(maf.record.brca.met)
maf.record.brca.met[1:3,1:3]

pheatmap(maf.record.brca.met,scale="none")
heatmap.2(as.matrix(maf.record.brca.met), scale="none",trace="none")

#What do patients with metabolite associated mutations have in common? - Look at breast cancer expressiong - NOTHING APPARENT IN TERMS OF MUTATIONS OR EXPRESSION
BRCA.EXP.OBJ<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/102514.CANCER.MATRICES.NORMALIZED.OBJ.rds")
BRCA.EXP.OBJ$combined.matrices[, unique(maf.record[Hugo_Symbol %in% met.enriched.filtered[CANCER=="BRCA",]$PATIENT)]
BRCA.MET.PATIENTS<-unique(maf.record[Hugo_Symbol %in% met.enriched.filtered[CANCER=="BRCA",]$GENE,]$PATIENT)

BRCA.MET.PATIENTS.ANNO<-data.frame(row.names=colnames(BRCA.EXP.OBJ$combined.matrices), 
                                   BRCA.MET=ifelse(colnames(BRCA.EXP.OBJ$combined.matrices) %in% BRCA.MET.PATIENTS, "METABOLIC", "NON.METABOLIC" ))
BRCA.MET.PATIENTS.ANNO<-BRCA.MET.PATIENTS.ANNO[order(BRCA.MET.PATIENTS.ANNO$BRCA.MET),,drop=F]
pheatmap(BRCA.EXP.OBJ$combined.matrices[intersect(rownames(BRCA.EXP.OBJ$combined.matrices), met.enriched.filtered[CANCER=="BRCA",]$GENE ),
                                        rownames(BRCA.MET.PATIENTS.ANNO)], cluster_cols=F,
                                       scale="none",trace="none", annotation=BRCA.MET.PATIENTS.ANNO)

unique(met.enriched.filtered[CANCER=="UCEC",]$METABOLITE)

######012715#######
#Plotting enrichment found per cancer type
library(grid)
mr.enrich<-fread("~/Desktop/FOLDER/GRANTS/JANUARY.15/Expressio.results.csv", header=T, stringsAsFactors=F, sep=",")
mr.enrich$SIGNIFICANT<-ifelse(mr.enrich$N.DIFF.PVAL<=0.05, "SIGNIFICANT", "NOT SIGNIFICANT")
mr.enrich$CANCER<-as.factor(mr.enrich$CANCER)
mr.enrich<-mr.enrich[order(CANCER, METABOLITE),]
mr.enrich$LIMITS<-1:nrow(mr.enrich)
mr.enrich[,LIMIT.MAX:=max(LIMITS), by=CANCER]
mr.enrich[,LIMIT.MIN:=min(LIMITS), by=CANCER]
mr.enrich$SIGNIFICANCE<--log(mr.enrich$N.DIFF.PVAL)
mr.enrich$SIGNIFICANCE<-ifelse(mr.enrich$SIGNIFICANCE==Inf, 7, mr.enrich$SIGNIFICANCE)

ggplot(mr.enrich, aes(LIMITS, SIGNIFICANCE)) + geom_bar(stat="identity", position="dodge", width=1)  + 
    coord_polar(theta="x") + geom_hline(y=-log(0.05), colour="red") +
  theme(panel.background = element_rect(fill = "white"))  +
  geom_rect(aes(fill=CANCER, xmax=LIMIT.MAX+0.5, xmin=LIMIT.MIN-0.5, ymax=10, ymin=11), colour="white") +
  geom_text(aes(y=10.5, label=METABOLITE), size=4, angle=-80 - 360/length(unique(mr.enrich$LIMITS)) * seq_along(mr.enrich$LIMITS) ) + 
  theme(plot.margin=unit(c(0.1,0.1,0.1,0), "lines")) +
  scale_x_continuous(minor_breaks=unique(c(mr.enrich$LIMIT.MAX+0.5, unique(mr.enrich$LIMIT.MIN)[-1]-0.5))) +
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
  xlab("Driver Metabolites Beyond Red Threshold of Significance for Each Cancer Type (P-value <0.05)") + 
  ylab("Significance Measured by -log(P-value)")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
