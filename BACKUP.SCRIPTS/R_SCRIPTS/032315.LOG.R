###########FUNCTIONS############
library(data.table)
library(reshape2)
library(ggplot2)
library(gplots)
library(pheatmap)
library(made4)
library(gtools)
library(parallel)
library(network)
library(XML)
library(combinat)

Function.Prep.MAF<-function(maf.file) {
  
  require(car)
  
  #Load cancer data
  maf<-fread(maf.file, header=T, sep="\t",stringsAsFactors=F)
  maf<-maf[,c("Hugo_Symbol","Chrom","Start_Position","Variant_Classification","Variant_Type", "Tumor_Sample_Barcode",
              "Reference_Allele", "Tumor_Seq_Allele2"),with=F]
  
  #Filter for "Unknown" gene
  maf<-maf[Hugo_Symbol!="Unknown",]
  
  #Unique
  setkey(maf)
  maf<-unique(maf)
  
  #Remove silent mutations
  #maf<-maf[Variant_Classification!="Silent",]
  
  #Equalize REF-ALT pairs
  maf$PRE.REF.ALT<-ifelse(maf$Variant_Type=="SNP", paste(as.vector(maf$Reference_Allele), as.vector(maf$Tumor_Seq_Allele2), sep="_"), "-")
  maf$REF.ALT<-ifelse(maf$PRE.REF.ALT!="-",recode(maf$PRE.REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" '),"-")
  maf$PRE.REF.ALT<-NULL
  
  #Separate by type of mutation
  non.func.class<-c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Nonsense_Mutation","Splice_Site", "In_Frame_Del","In_Frame_Ins")
  non.func.type<-c("DEL","INS")
  maf$TYPE<-ifelse(maf$Variant_Classification=="Silent", "SILENT",
                   ifelse(maf$Variant_Classification %in% non.func.class, "NON.FUNC", 
                          ifelse(maf$Variant_Type %in% non.func.type, "NON.FUNC", "MISS")))
  
  #Convert IDs to expression sample IDs
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste0(unlist(strsplit(x, "-"))[1:4],collapse="." ))
  maf$Tumor_Sample_Barcode<-NULL
  
  #Get one sample representing all vials
  maf$VIAL<-substr(maf$SAMPLE, 14,16) #Get vial
  maf$SAMPLE<-substr(maf$SAMPLE,1,12)
  maf[,REP.VIAL:= names(sort(table(VIAL),decreasing=T))[1] , by=SAMPLE] #Get highest counting vial ID
  maf$SAMPLE<-paste( maf$SAMPLE, maf$REP.VIAL, sep="." ) #Rename sample with highest counting vial ID
  maf$VIAL<-NULL
  maf$REP.VIAL<-NULL
  
  #Classify mutations by classes
  maf$CLASS<-ifelse((maf$TYPE=="MISS" | maf$TYPE=="SILENT"), 
                    paste(maf$Hugo_Symbol, maf$Start_Position, maf$REF.ALT, maf$TYPE, sep="."), 
                    paste(maf$Hugo_Symbol, maf$TYPE, sep="."))
  
  #Count how many samples are covered by a type of mutation
  maf[,POP.CLASS:=length(SAMPLE),by="CLASS"]
  
  #Cleand up and Return
  setkey(maf)
  maf<-unique(maf)
  maf<-maf[order(POP.CLASS, decreasing=T),]
  return(maf)
}

normalize.vector<-function(x){
  y=(x-min(x))/(max(x)- min(x))
  return(y)
}

#########LOAD FILES############
BRCA.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf")

GBM.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/GBM/031315/Somatic_Mutations/BI__IlluminaGA_DNASeq/Level_2/broad.mit.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf")

AML.MAF<-Function.Prep.MAF("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/AML/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf")

BRCA.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/022615.BRCA.CANCER.MATRICES.NORMALIZED.OBJ.NB.rds")

PATHWAY<-fread("DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND", header=T, sep="\t", stringsAsFactors=F, colClasses=rep("character",3))
PATHWAY<-PATHWAY[PATHWAY!="01100",]

##########ANALYSIS###########
KRAS.SAMPLES<-BRCA.MAF[Hugo_Symbol %in% c("BRAF") ,]$SAMPLE
#KRAS.SAMPLES<-BRCA.MAF[Hugo_Symbol %in% c("KRAS") & CLASS=="KRAS.25398284.C_A.MISS" ,]$SAMPLE
#KRAS.SAMPLES<-BRCA.MAF[Hugo_Symbol=="KRAS",]$SAMPLE
NON.KRAS.SAMPLES<-setdiff(BRCA.EXP$cancer.patients, KRAS.SAMPLES)
TESTED.GENES<-c("SLC5A1","SLC2A2", "SLC2A1","SLC2A5") #SGLT1, ,GLUT2, GLUT1, GLUT5
TESTED.GENES %in% rownames(BRCA.EXP$combined.matrices)

pheatmap(BRCA.EXP$combined.matrices[TESTED.GENES, c(KRAS.SAMPLES, NON.KRAS.SAMPLES)], scale="none", cluster_cols=F,
         annotation=data.frame(MUTATION=c(rep("KRAS", length(KRAS.SAMPLES)), rep("NON.KRAS", length(NON.KRAS.SAMPLES))),
                               row.names=c(KRAS.SAMPLES, NON.KRAS.SAMPLES)))

KRAS.TABLE<-data.table(t(BRCA.EXP$combined.matrices[TESTED.GENES,]),keep.rownames=T)
KRAS.TABLE$KRAS<-KRAS.TABLE$rn %in% KRAS.SAMPLES
KRAS.TABLE<-melt(KRAS.TABLE, id.vars=c("rn", "KRAS"))

ggplot(KRAS.TABLE, aes(KRAS, value, colour=KRAS)) + geom_boxplot() + geom_jitter(size=0.7) + theme.format + facet_wrap(~variable,scales="free") +
  ylab("Expression")

KRAS.PVAL=data.table()
for (gene in TESTED.GENES){
  P.VAL=wilcox.test(KRAS.TABLE[variable==gene & KRAS==T,]$value, KRAS.TABLE[variable==gene & KRAS==F,]$value, alternative="greater", paired=F)$p.value
  KRAS.PVAL=rbind(KRAS.PVAL, data.table(gene=gene, p.val=P.VAL))
}
KRAS.PVAL

BRCA.MAF[Hugo_Symbol=="IDH1",]
BRCA.MAF[Hugo_Symbol=="PHGDH",]
AML.MAF[Hugo_Symbol=="IDH1",]
GBM.MAF[Hugo_Symbol=="IDH1",]

BRCA.EXP.MI.NORMALS<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/030415.BRCA.EXP.MI.NORMALS", header=T, sep="\t", stringsAsFactors=F)
BRCA.EXP.MI.NORMALS<-BRCA.EXP.MI.NORMALS[order(MI, decreasing=T),]
hist(BRCA.EXP.MI.NORMALS$MI)
BRCA.EXP.MI.NORMALS[Hugo.1=="BRAF" & Hugo.2=="KRAS",]

BRCA.EXP.COR.NORMAL<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$normal.patients]), method="spearman")
BRCA.EXP.COR.CANCER<-cor(t(BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients]), method="spearman")
genes<-c("TP53","PIK3CA","TTN",
         "CDH1","CDH2","CDH3","CDH4","CDH5",
         "GATA1","GATA2","GATA3","GATA4","GATA5","GATA6",
         "SLC5A1", "SLC2A2","TAS1R2","TAS1R3",
         "BRCA1", "BRCA2", "BRIP1", "PALB2", "FANCD2",
         "BRAF","KRAS")
heatplot(BRCA.EXP.COR.NORMAL[genes,genes],scale="none",dualScale=FALSE)
heatplot(BRCA.EXP.COR.CANCER[genes,genes],scale="none",dualScale=FALSE)

#Test home-made network ()
glu.genes<-c("SLC5A1", "SLC2A2","TAS1R2","TAS1R3")
BRCA.MAF[Hugo_Symbol %in% glu.genes,]
BRCA.MAF[Hugo_Symbol=="KRAS",]

BRCA.EXP.COR.NORMAL["SLC5A1",1:3]
SLC5A1<-data.table(Hugo.1="SLC5A1", Hugo.2=colnames(BRCA.EXP.COR.NORMAL) , COR=BRCA.EXP.COR.NORMAL["SLC5A1",])
SLC5A1<-SLC5A1[abs(COR)>0.35,]
SLC5A1<-SLC5A1[order(abs(COR), decreasing=T),]
hist(SLC5A1$COR)

TP53<-data.table(Hugo.1="TP53", Hugo.2=colnames(BRCA.EXP.COR.NORMAL) , COR=BRCA.EXP.COR.NORMAL["TP53",])
TP53<-TP53[abs(COR)>0.35,]
TP53<-TP53[order(abs(COR), decreasing=T),]

#To make enzymatic network, load enzyme information
ENZYME<-fread("DATABASES/KEGG/061214_ENZYME_PRODUCT",header=T, sep="\t", stringsAsFactors=F)
ENZYME<-ENZYME[Enzyme %in% rownames(BRCA.EXP.COR.NORMAL),]
heatplot(BRCA.EXP.COR.NORMAL[unique(ENZYME$Enzyme), unique(ENZYME$Enzyme)], scale="none", dualScale=FALSE)
heatplot(BRCA.EXP.COR.CANCER[unique(ENZYME$Enzyme), unique(ENZYME$Enzyme)], scale="none", dualScale=FALSE)

sum(unique(ENZYME$Enzyme) %in% rownames(BRCA.EXP.COR.CANCER))
ENZYME[!(Enzyme %in% rownames(BRCA.EXP.COR.NORMAL)),]
"RTC1" %in% rownames(BRCA.EXP.COR.NORMAL)

Function.COR.EXP.MUT<-function(cancer.exp, mut.class, maf, enzymes){
  #Construct expression correlation matrix for mut class using cancer expression data
  
  #Filter out silent mutations
  maf<-maf[TYPE!="SILENT",]
  
  #Filter for mut class
  maf<-maf[CLASS==mut.class,]
  
  #Filter exp matrix for target samples
  cancer.exp<-cancer.exp[,intersect(colnames(cancer.exp), unique(maf$SAMPLE))]
  
  #Filter exp matrix for target enzymes
  cancer.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  
  #Construct correlation matrix
  cancer.cor<-cor(t(cancer.exp), method="spearman")
  
  #Return
  return(cancer.cor)
}
BRCA.BRCA1.NON.FUNC<-Function.COR.EXP.MUT(BRCA.EXP$combined.matrices, "BRCA1.NON.FUNC", BRCA.MAF, unique(ENZYME$Enzyme))
BRCA.BRCA1.NON.FUNC[1:3,1:3]
heatplot(BRCA.BRCA1.NON.FUNC, scale="none",dualScale=F)

Function.COR.EXTRACT<-function(main.cor, back.cor){
  #Function to substract background correlation signal to cancer one
  
  #Order of genes in main.cor
  genes.order<-rownames(main.cor)
  
  #Substract
  signal.cor<-main.cor-back.cor[genes.order, genes.order]
  
  #Cluster in terms of original back.cor
  #dendo.r<-hclust(dist(back.cor))
  #dendo.h<-hclust(dist(t(back.cor)))
  
  #Return
  return(signal.cor)
}

BRCA1.NON.FUNC.SUBSTRACTED<-Function.COR.EXTRACT(BRCA.BRCA1.NON.FUNC, BRCA.EXP.COR.NORMAL)
heatplot(BRCA1.NON.FUNC.SUBSTRACTED, scale="none", dualScale=F)

hist(BRCA1.NON.FUNC.SUBSTRACTED)

Function.COR.ENZYME.DIFF<-function(main.cor, back.cor){
  #Function to calculate the average change in enzyme in network, that is change in functionality with respect to functional module
  
  #Order of genes in main.cor
  genes.order<-rownames(main.cor)
  
  #Substract
  signal.cor<-main.cor-back.cor[genes.order, genes.order]
  
  #Get sum of absolute value of change of enzyme with respect to all (LOG-CONVERTED)
  sum.delta.abs<-log(apply(signal.cor, 1, function(x) sum(abs(x))))
  
  #Organize to table
  main.table<-matrix(sum.delta.abs)
  rownames(main.table)<-rownames(signal.cor)
  #main.table<-data.table(Hugo_Symbol=rownames(signal.cor), SUM.DELTA.ABS=sum.delta.abs)
  
  #Clean up and return
  #main.table<-main.table[order(SUM.DELTA.ABS, decreasing=T),]
  return(main.table)
}

BRCA1.NON.FUNC.ENZYME.CHANGE<-Function.COR.ENZYME.DIFF(BRCA.BRCA1.NON.FUNC, BRCA.EXP.COR.NORMAL)
head(BRCA1.NON.FUNC.ENZYME.CHANGE)

Function.COR.ENZYME.BACKGROUND<-function(maf, mut.class, cancer.exp, back.cor, enzymes){
  #Function to build background enzyme correlation distribution based on random sampling of n samples corresponding to mutation class
  
  #Filter cancer exp matrix for common patient
  cancer.exp<-cancer.exp[,intersect(colnames(cancer.exp), unique(maf$SAMPLE))]
  
  #Filter out silent mutations
  maf<-maf[TYPE!="SILENT",]
  
  #Filter for mut class
  maf<-maf[CLASS==mut.class,]
  
  #Get number of samples for mut class from intersect with exp data
  n.samples<-length(intersect(unique(maf$SAMPLE), colnames(cancer.exp)))
  
  #Filter exp matrix for target enzymes
  cancer.exp<-cancer.exp[intersect(enzymes,rownames(cancer.exp)),]
  
  #Filter background correlation for enzyme genes found in cancer matrix
  back.cor<-back.cor[rownames(cancer.exp), rownames(cancer.exp)]
  
  #Random sampling from expression matrix (100x)
  print (dim(cancer.exp))
  ran.samples<-replicate(1000,sample(colnames(cancer.exp), n.samples), simplify=F)
  
  #Obtain gene order for back correlation
  gene.order<-rownames(back.cor)
  
  #Ran correlation simulations
  main.list<-lapply(ran.samples, function(x) {
    
    #Construct correlation matrix for random samples
    main.cor<-cor(t(cancer.exp[,x]), method="spearman")
    
    #Substract
    signal.cor<-main.cor[gene.order, gene.order]-back.cor
    
    #Get sum of absolute value of change of enzyme with respect to all - LOG-TRANSFORMED
    sum.delta.abs<-log(apply(signal.cor, 1, function(x) sum(abs(x))))
    
    #Return
    return(sum.delta.abs)
  })
  
  #Assign to genes
  main.table<-do.call(cbind, main.list)
  rownames(main.table)<-gene.order
  
  #Clean up of NAs
  main.table<-t(main.table)
  main.table<-main.table[complete.cases(main.table),]
  main.table<-t(main.table)
  
  #Return
  return(main.table)
}

BRCA1.NON.FUNC.RANDOM.COR.SUBS<-Function.COR.ENZYME.BACKGROUND(BRCA.MAF, "BRCA1.NON.FUNC", BRCA.EXP$combined.matrices[, BRCA.EXP$cancer.patients],
                                                               BRCA.EXP.COR.NORMAL, unique(ENZYME$Enzyme))
BRCA1.NON.FUNC.RANDOM.COR.SUBS[1:3,1:3]
dim(BRCA1.NON.FUNC.RANDOM.COR.SUBS)

Function.PVALUE<-function(change.matrix, ran.matrix){
  #Calculate p value
  
  #Interest genes
  genes<-rownames(change.matrix)
  
  #Calculate empirical p-value
  p.vals<-sapply(genes, function(x) mean(as.vector(change.matrix[x,])<=ran.matrix[x,])  )
  
  #Build results
  main.table<-data.table(Hugo_Symbol=genes, P.VAL=p.vals)
  
  #Correct for multiple hypothesis testing
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Clean up and Return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}
BRCA1.NON.FUNC.PVAL<-Function.PVALUE(BRCA1.NON.FUNC.ENZYME.CHANGE, BRCA1.NON.FUNC.RANDOM.COR.SUBS )

######DO IT AT ENZYME LEVEL#######
setnames(ENZYME, c("Hugo_Symbol", "KEGG_ID", "Product"))
setkey(ENZYME)
ENZYME<-unique(ENZYME)

ENZYME[KEGG_ID=="C00024",]

Function.KEGG.PATH.TO.ENZYME<-function(kegg.hsa.url, kegg.path.hsa.url, pathway.url){
  
  ####Construct kegg hugo symbol to kegg path id descriptor####
  #Load kegg gene list
  genes<-fread(kegg.hsa.url, header=F, sep="\t", stringsAsFactors=F)
  
  #Obtain hugo from kegg id
  genes$Hugo_Symbol<-sapply(genes$V2, function(x) unlist(strsplit(unlist(strsplit(x, ","))[1], ";"))[1])
  genes$V2<-NULL
  genes<-genes[!(grepl("uncharacterized",Hugo_Symbol)),]
  
  #Load kegg gene to path list
  genestopath<-fread(kegg.path.hsa.url, header=F, sep="\t", stringsAsFactors=F)
  
  #Merge to obtain genes in annotated paths
  main.table<-merge(genes, genestopath, by="V1")
  main.table$V1<-NULL
  setnames(main.table, c("Hugo_Symbol", "PATHWAY"))
  
  ####Attach pathway description to kegg path id####
  #Load pathway description file
  keggpathdescr<-fread(pathway.url, header=F, sep="\t", stringsAsFactors=F)
  setnames(keggpathdescr, c("PATHWAY", "DESCRIPTION"))
  
  #Merge
  main.table<-merge(main.table, keggpathdescr, by="PATHWAY")
  
  ####Clean up and return####
  #Clen path id
  main.table$PATHWAY<-sapply(main.table$PATHWAY, function(x) unlist(strsplit(x, ":", fixed=T))[2])
  
  #Filter out irrelevant pathways (i.e. metabolic pathways 01100)
  main.table<-main.table[PATHWAY!="hsa01100",]
  
  #Remove duplicates
  setkey(main.table)
  main.table<-unique(main.table)
  
  #Return
  return(main.table)
}

PATHWAY<-Function.KEGG.PATH.TO.ENZYME("http://rest.kegg.jp/list/hsa", "http://rest.kegg.jp/link/pathway/hsa", "http://rest.kegg.jp/list/pathway/hsa")
length(unique(PATHWAY$Hugo_Symbol))

BRCA.DIFF.EXP<-Function.RNAseq.Differential.Expression.V2(BRCA.EXP, BRCA.EXP$cancer.patients)

Function.test.1<-function(cancer.cor, normal.cor, enzyme.table, pathway.table, diff.exp.table){
  
  internal.function<-function(genes, pathway){
    
    ######First look at level of level of dyscorrelation of metabolite with respect to pathway#####
    #Find difference in correlation
    non.genes<-setdiff(unique(pathway.list[[pathway]]$Hugo_Symbol), genes)
    cancer<-cancer.cor[genes, non.genes, drop=F]
    normal<-normal.cor[genes, non.genes, drop=F]
    diff.cor<-abs(cancer-normal) #degree of difference in correlation between cancer and normal
    
    #Look for median difference
    diff.cor<-apply(diff.cor, 1, function(x) median(x))
    diff.cor<-median(diff.cor)
    
    #####Then level of diff exp for metabolite#####
    diff.exp<-diff.exp.table[Hugo_Symbol %in% genes,]
    diff.exp<-median(abs(diff.exp$LOG.FC))
    
    #####Return#####
    return(list(COR.SCORE=diff.cor, EXP.SCORE=diff.exp))
  }
  
  ######Obtain annotations#######
  kegg.anno<-enzyme.table[,list(PRODUCT=Product[1]), by="KEGG_ID"]
  path.anno<-unique(pathway.table[,c("PATHWAY","DESCRIPTION"), with=F])
  
  ######Process enzyme table#####
  enzyme<-enzyme.table[,c("Hugo_Symbol", "KEGG_ID"), with=F]
  
  #Filter out small molecules and co-factors (NEEDS UPDATE)
  enzyme<-enzyme[!(KEGG_ID %in% c("C00001", "C00003","C00004", "C00080", "C00005","C00006","C00008","C00011","C00014",
                                  "C00015","C00035","C00013","C00009")),]
  
  setkey(enzyme)
  enzyme<-unique(enzyme)
  
  #Filter out on pathway
  pathway.table<-pathway.table[PATHWAY!="hsa05146",]
  
  #Combine enzyme and pathway information#
  enzyme<-merge(enzyme, pathway.table, by="Hugo_Symbol", allow.cartesian=T)
  
  #####Co-filter enzyme, cor and diff tables for genes#####
  diff.exp.table<-diff.exp.table[,c("logFC","adj.P.Val","ID"), with=F]
  setnames(diff.exp.table, c("LOG.FC","ADJ.PVAL", "Hugo_Symbol"))
  all.genes<-intersect(intersect(unique(enzyme$Hugo_Symbol), rownames(cancer.cor)), diff.exp.table$Hugo_Symbol)
  
  diff.exp.table<-diff.exp.table[Hugo_Symbol %in% all.genes,]
  enzyme<-enzyme[Hugo_Symbol %in% all.genes,]
  cancer.cor<-cancer.cor[all.genes, all.genes]
  normal.cor<-normal.cor[all.genes, all.genes]
  
  #####Get pathway gene list#####
  #Filter out pathways that don't have at least one kegg
  enzyme[,POP:=length(unique(KEGG_ID)), by="PATHWAY"]
  enzyme<-enzyme[POP>=2,]
  enzyme$POP<-NULL
  pathway.list<-split(enzyme, enzyme$PATHWAY, drop=T)
  
  ####Obtain dysregulation and dyscorrelation score for cancer vs normal per metabolite in pathway
  main.table<-enzyme[,internal.function(Hugo_Symbol, PATHWAY),by=c("PATHWAY", "KEGG_ID")]
  
  #Obtain total score
  main.table$EXP.SCORE<-normalize.vector(main.table$EXP.SCORE)
  main.table$TOTAL.SCORE<-main.table$COR.SCORE + main.table$EXP.SCORE
  
  ####Clean up and return
  main.table<-merge(main.table, kegg.anno, by="KEGG_ID")
  main.table<-merge(main.table, path.anno, by="PATHWAY")
  main.table<-main.table[order(TOTAL.SCORE, decreasing=T),]
  return(main.table)
}

test.1<-Function.test.1(BRCA.EXP.COR.CANCER, BRCA.EXP.COR.NORMAL, ENZYME, PATHWAY, BRCA.DIFF.EXP)
test.1[1:100,]
hist(test.1$TOTAL.SCORE)
write.table(test.1, "LOGS/040815.BRCA.MET.SCORES", quote=F,sep="\t",row.names=F, col.names=T)

BRCA.MAF[grepl("B6.A0RD",SAMPLE),]
length(unique(ENZYME$KEGG_ID))

######Process metabolomic profiles from TANG et al. paper######
TANG.PROF<-fread("DATABASES/METABOLOMICS/TANG.2014/CLEANED.PROFILES.csv", header=T, sep=",", stringsAsFactors=F)
colnames(TANG.PROF)
data.frame(TANG.PROF, row.names=as.vector(TANG.PROF$METABOLITE))[1:3,1:3]

Function.Process.TANG<-function(tang.cleaned, kegg.cpd.url, hmdb, extra.ann){
  require(stringr)
  #Function to convert tang cleaned file into processed matrix
  
  ######Process tang file into matrix########
  #Load file
  main.table<-fread(tang.cleaned, header=T, sep=",", stringsAsFactors=F,)
  
  #Process into TCGA ID format
  samples<-colnames(main.table)[2:ncol(main.table)]
  ids<-c()
  toll<-1
  for (record in samples){
    if (grepl("-", record)==T){
      ids<-c(ids, paste("TCGA.", str_trim(gsub("-",".", record)), sep=""))
    } else {
      ids<-c(ids, paste("Normal", toll, sep="."))
      toll<-toll+1
    }
  }
  setnames(main.table, c("METABOLITE", ids))
  
  ######Feature filter tang metabolite names######
  main.table$METABOLITE<-sapply(main.table$METABOLITE, function(x) unlist(strsplit(unlist(strsplit(x, "*", fixed=T))[1], " \\("))[1])
  all.metabolites<-as.vector(toupper(unique(main.table$METABOLITE)))
  
  ######Obtain kegg identifiers ############
  print ("annotating kegg_ids from kegg")
  #Load kegg compounds
  cpds<-fread(kegg.cpd.url, header=F, sep="\t", stringsAsFactors=F)
  setnames(cpds, c("CPD", "PRODUCT"))
  
  #Process product names
  cpds<-cpds[,list(METABOLITE=unlist(strsplit(PRODUCT, "; "))), by="CPD"]
  
  ######Merge to get KEGG identifiers and convert to matrix######
  main.table$METABOLITE<-toupper(main.table$METABOLITE)
  cpds$METABOLITE<-toupper(cpds$METABOLITE)
  main.table<-merge(main.table, cpds, by="METABOLITE",all.x=T)
  
  ######Separate for hmdb annotation######
  main.table.hmdb<-main.table[is.na(CPD),]
  main.table.kegg<-main.table[!(is.na(CPD)),]
  
  #Apply HMDB annotation - add extra potential synonyms
  print ("annotating kegg_ids from hmdb")
  hmdb$METABOLITE<-toupper(hmdb$METABOLITE) 
  hmdb.extra<-hmdb[grepl("-", METABOLITE),]
  hmdb.extra$METABOLITE<-sapply(hmdb.extra$METABOLITE, function(x) paste(unlist(strsplit(x, "-")), collapse=" "))
  hmdb<-unique(rbind(hmdb, hmdb.extra))
  
  main.table.hmdb$CPD<-NULL
  main.table.hmdb<-merge(main.table.hmdb, hmdb, by="METABOLITE", all.x=T)
  
  ######Separate for independent annotation#####
  main.table.extra<-main.table.hmdb[is.na(KEGG_ID),]
  main.table.hmdb<-main.table.hmdb[!(is.na(KEGG_ID)),]
  
  #Apply extra annotation
  setnames(extra.ann, c("METABOLITE", "KEGG_ID"))
  extra.ann$METABOLITE<-toupper(extra.ann$METABOLITE)
  extra.ann<-extra.ann[KEGG_ID!="NONE",]
  
  main.table.extra$KEGG_ID<-NULL
  main.table.extra<-merge(main.table.extra, extra.ann, by="METABOLITE")

  ######List those unidentifiable by a KEGG ID metabolites#######
  not.found<-setdiff(all.metabolites, c(main.table.hmdb$METABOLITE, main.table.kegg$METABOLITE, main.table.extra$METABOLITE))
  main.table.hmdb$METABOLITE<-NULL
  main.table.extra$METABOLITE<-NULL
  
  ######Clean up into matrices######
  print ("Merging matrices")
  print (main.table.kegg[,1:5,with=F])
  print (main.table.hmdb[,1:5, with=F])
  #Extract KEGG IDs from kegg table
  main.table.kegg$KEGG_ID<-sapply(main.table.kegg$CPD, function(x) unlist(strsplit(x, ":"))[2])
  main.table.kegg$CPD<-NULL
  main.table.kegg$METABOLITE<-NULL
  
  #Combine kegg and hmdb tables
  main.table<-rbind(main.table.kegg, main.table.hmdb[,colnames(main.table.kegg), with=F])
  
  #Melt and cast to mean expression over kegg duplciates
  main.table<-melt(main.table, id.vars="KEGG_ID")
  main.table<-dcast.data.table(main.table, KEGG_ID~variable,fun.aggregate=mean)
  print (main.table)
  
  #Convert to matrix with rownames as metabolites
  print (main.table[KEGG_ID %in% c("C00003", "C04677", "C16513"),c("KEGG_ID", colnames(main.table)[1:12]), with=F])
  print (length(main.table$KEGG_ID))
  print (length(unique(main.table$KEGG_ID)))
  
  main.table<-data.frame(main.table, row.names=as.vector(main.table$KEGG_ID))
  main.table$KEGG_ID<-NULL
  main.table<-as.matrix(main.table)
  
  #######Return######
  return(list(MAIN.TABLE=main.table, NOT.FOUND=not.found))
}

TANG.PROF<-Function.Process.TANG("DATABASES/METABOLOMICS/TANG.2014/CLEANED.PROFILES.csv", "http://rest.kegg.jp/list/cpd/", HMDB, NOT.FOUND.ANN)
dim(TANG.PROF$MAIN.TABLE)
heatplot(log(TANG.PROF$MAIN.TABLE[apply(TANG.PROF$MAIN.TABLE, 1, sd)!=0,]),scale="none",dualScale=FALSE, margins=c(7,5))

length(colnames(TANG.PROF$MAIN.TABLE))
setdiff(colnames(TANG.PROF$MAIN.TABLE), substr(BRCA.EXP$cancer.patients, 1, 12))
setdiff(colnames(TANG.PROF$MAIN.TABLE), substr( unique(BRCA.MAF$SAMPLE) , 1, 12))

#Manually annotate not found
write.table(data.table(TANG.PROF$NOT.FOUND), "LOGS/NOT.FOUND", quote=F, sep="\t",row.names=F, col.names=F)
NOT.FOUND.ANN<-fread("DATABASES/METABOLOMICS/TANG.2014/NOT.FOUND.ANNOTATED", header=F, sep="\t", stringsAsFactors=F)

#Load HMDB xml
Function.HMDB<-function(hmdb.folder){
  require(XML)
  
  #Obtain xml files and filter out large XML if possible
  xmls<-list.files(hmdb.folder)
  xmls<-xmls[grepl("HMDB", xmls, ignore.case=F)]
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "xmls", "xmlToList") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-parLapply(cl, xmls, function(x) {
    
    #Load file
    xml.file<-paste(hmdb.folder, x, sep="/")
    xml.table<-xmlToList(xml.file)
    
    #Only store if there is an actual kegg id
    if (is.null(xml.table$kegg_id)==F){
      xml.table<-data.table(METABOLITE=c(xml.table$name, as.vector(sapply(xml.table$synonyms, function(y) y))), KEGG_ID=xml.table$kegg_id)
      
      #Return
      return(xml.table)
    }
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Converge list
  main.table<-do.call(rbind, main.list)
  
  #Clean up and Return
  main.table<-unique(main.table)
  main.table<-main.table[METABOLITE!="\n  ",]
  return(main.table)
}

HMDB<-Function.HMDB("DATABASES/HMDB_SOURCE/HMDB.3.6/hmdb_metabolites")
BRCA.MAF[Hugo_Symbol=="PEPD",]

#####Analyze annotated metabolites from Tang et al.######
dim(TANG.PROF$MAIN.TABLE)
Function.TANG.DIFF.MET<-function(tang.table){
  
  #Separate constant data and assign p.val=1 later
  cons.met<-rownames(tang.table[apply(tang.table, 1, function(x) sd(x))==0,])
  tang.table<-tang.table[apply(tang.table, 1, function(x) sd(x))!=0,]
  
  #Obtain patients
  normals<-paste("Normal", seq(1,5,1), sep=".")
  cancer<-setdiff(colnames(tang.table), normals)
  
  #Obtain diff metabolites
  p.vals<-apply(tang.table, 1, function(x) t.test(x[normals], x[cancer], paired=F, var.equal=F)$p.value)
  
  #Construct table
  main.table<-data.table(KEGG_ID=rownames(tang.table), P.VAL=p.vals)
  
  #Add constant data
  main.table<-rbind(main.table, data.table(KEGG_ID=cons.met, P.VAL=1))
  
  #Correct for multiple hypothesis testing
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #Clean up and return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}
TANG.DIFF.MET<-Function.TANG.DIFF.MET(TANG.PROF$MAIN.TABLE)
ggplot(TANG.DIFF.MET, aes(P.VAL.ADJ)) + geom_histogram() + theme.format + ylab("Metabolite count") + xlab("Corrected p-values") +
  geom_vline(xintercept=0.05, colour="red", linetype="dashed")

####Can we predict metabolite levels from enzyme levels?#####

#Combinatorial approach
Function.Comb.Met.Pred<-function(tang.proc.table, exp.obj, enzyme.table){
  
  #####Obtain average expression to normal for cancer patients and re-construct matrix######
  gene.exp.avg<-t(apply(exp.obj$combined.matrices, 1, function(x) x[exp.obj$cancer.patients] - mean(x[exp.obj$normal.patients])))
  colnames(gene.exp.avg)<-sapply(colnames(gene.exp.avg), function(x) substr(x,1,12))
  print (gene.exp.avg[1:3,1:3])
  
  #####Obtain average metabolic expression to normal for cancer patients####
  met.normal<-paste("Normal.", seq(1,5,1), sep="")
  met.cancer<-setdiff(colnames(tang.proc.table), met.normal)
  met.exp.avg<-t(apply(tang.proc.table, 1, function(x) x[met.cancer] - mean(x[met.normal])))
  met.exp.avg<-data.table(melt(met.exp.avg))
  setnames(met.exp.avg, c("KEGG_ID", "SAMPLE", "MET.EXP"))
  print (met.exp.avg)
  
  #####Convert expression table into met table by obtaining the median expression of enzymes to metabolite#####
  #Reduce for common patients and enzymes to metabolite expression table
  common.patients<-intersect(colnames(gene.exp.avg), met.cancer)
  print(common.patients)
  common.enzymes<-unique(enzyme.table[KEGG_ID %in% unique(met.exp.avg$KEGG_ID), ]$Hugo_Symbol)
  print (length(common.enzymes))
  gene.exp.avg<-gene.exp.avg[intersect(common.enzymes, rownames(gene.exp.avg)) , common.patients]
  print (gene.exp.avg[1:3,1:3])
  
  #Convert to data table and merge with enzyme table to obtain common KEGG IDs
  gene.met<-data.table(gene.exp.avg, keep.rownames=T)
  setnames(gene.met, c("Hugo_Symbol", colnames(gene.met)[2:ncol(gene.met)]))
  gene.met<-merge(gene.met, unique(enzyme.table[,1:2,with=F]), by="Hugo_Symbol")
  gene.met$Hugo_Symbol<-NULL
  print(gene.met[1:3,1:3, with=F])
  
  #Melt to obtain median expression of enzymes to predict metabolite expression
  gene.met<-melt(gene.met, id.vars="KEGG_ID")
  gene.met<-gene.met[,list(KEGG.PRED.EXP=median(value)), by=c("KEGG_ID", "variable")]
  setnames(gene.met, c("KEGG_ID", "SAMPLE", "MET.PRED.EXP"))
  print (gene.met)
  
  #####Obtain correlation table of metabolite and average enzyme fold expression##### 
  main.table<-merge(met.exp.avg, gene.met, by=c("KEGG_ID", "SAMPLE"))
  main.table<-main.table[,list(RHO=cor.test(MET.EXP, MET.PRED.EXP, method="spearman")$estimate,
                               P.VAL=cor.test(MET.EXP, MET.PRED.EXP, method="spearman")$p.value) ,by="KEGG_ID"]
  
  #Correct for multiple hypothesis testing
  main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
  
  #####Clean up and return#####
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}

COMB.MET.PRED<-Function.Comb.Met.Pred(TANG.PROF$MAIN.TABLE, BRCA.EXP, ENZYME)
ggplot(COMB.MET.PRED, aes(RHO)) + geom_histogram(binwidth=0.05) + theme.format + 
  ggtitle("Correlation Distribution of Enzyme prediction to Actual Metabolic Levels")
ggplot(COMB.MET.PRED, aes(P.VAL.ADJ)) + geom_histogram(binwidth=0.05) + theme.format + 
  geom_vline(xintercept=0.05, colour="red", linetype="dashed") +
  ggtitle("P-value Distribution of Correlation of Metabolic Prediction")

setdiff(rownames(TANG.PROF$MAIN.TABLE), unique(ENZYME$KEGG_ID))[1:10]
ENZYME[Hugo_Symbol=="DBH",]