#PROPOSAL.R
library(data.table)
library(reshape2)
library(ggplot2)

kegg.edges<-Function.kegg.filtered("DATABASES/KEGG/071415.ENZYME.SUB.PROD.MAIN.FILT", "DATABASES/RECON/042215.PROCESSED.METABOLITES", weight.filter = 60, n.edge.filter = 60)

enz.prod<-fread("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_ENZYME_PRODUCT", header = T)

gene.length<-fread("DATABASES/UNIPROT/042314_GENE_LENGTH", header=T, sep="\t", stringsAsFactors = F)

tcga.mut<-Function.tcga.mut.prep("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", mut.filter=5) #at least 5 people have all mutations
length(unique(tcga.mut[Hugo_Symbol %in% unique(kegg.edges$Hugo_Symbol), ]$SAMPLE)) #833@1, 726@5

Function.D.3.3<-function(tcga.mut, kegg.edges, gene.length, enz.prod=F){
  
  #Get associated metabolites per gene 
  if (enz.prod==T){
    kegg<-unique(kegg.edges[,c("Enzyme", "KEGG_ID"), with=F])
    setnames(kegg, c("Hugo_Symbol", "MET"))
  } else {
    kegg<-rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT))
    kegg<-unique(kegg) 
  }
  
  #Obtain mutation rate per metabolite
  main.table<-merge(merge(kegg, tcga.mut, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")
  filtered.genes<-unique(main.table[,c("Hugo_Symbol","MUTATION"), with=F])
  main.table<-main.table[, list(RATE=length(Hugo_Symbol)/(length(SAMPLE) * sum(LENGTH)), GENE.COUNT=length(unique(Hugo_Symbol))), by=c("MET", "MUTATION")]
  
  #Clean up and Return
  main.table<-main.table[!is.na(RATE),]
  return(list(MET.TABLE=main.table, GENE.TABLE=filtered.genes))
}

main.D.3.3<-Function.D.3.3(tcga.mut, enz.prod, gene.length, enz.prod = T)
main.D.3.3$GENE.TABLE

ggplot(main.D.3.3$MET.TABLE, aes(GENE.COUNT, RATE)) + geom_point() + facet_wrap(~MUTATION) + scale_x_log10() #Have to draw from same background
ggplot(main.D.3.3$MET.TABLE, aes(RATE)) + geom_histogram() + facet_wrap(~MUTATION) 
ggplot(main.D.3.3$MET.TABLE, aes(MUTATION,RATE)) + geom_boxplot()

Function.D.3.3.NULL<-function(main.table, gene.length, tcga.mut){
  
  #Filter tcga to correct pools
  tcga.pool<-merge(tcga.mut, gene.length, by="Hugo_Symbol")
  tcga.pool<-unique(tcga.pool[,c("SAMPLE", "Hugo_Symbol", "MUTATION"), with=F])
  #tcga.pool<-tcga.pool[Hugo_Symbol %in% unique(kegg.edges$Hugo_Symbol),]
  
  #Obtain metabolite gene count based on D.3.3 results
  types<-unique(main.table$MUTATION)
  
  score.table<-data.table()
  for (t in types){
    gene.counts<-unique(main.table[MUTATION==t,]$GENE.COUNT)
    tcga.genes<-unique(tcga.pool[MUTATION==t,]$Hugo_Symbol)
    
    print (t)
    print (gene.counts)
    for (gc in gene.counts){
      print (gc)
      random.pool<-replicate(100, sample(tcga.genes, gc))
       
      #Calculate score per random pool
      SCORES=sapply(random.pool, function(x) gc/(length(tcga.pool[MUTATION==t,][Hugo_Symbol %in% x,]$SAMPLE)*sum(gene.length[Hugo_Symbol %in% x,]$LENGTH)))
      
      #Store results
      score.table<-rbind(score.table, data.table(MUTATION=t, GENE.COUNT=gc, SCORES=SCORES))
    }
  }
  
  #Return
  return(score.table)
}

D.3.3.NULL<-Function.D.3.3.NULL(main.D.3.3$MET.TABLE[MUTATION!="SILENT",], gene.length, tcga.mut)

Function.D.3.3.PVAL<-function(main.table, null.table){
  
  setnames(null.table, c("MUT", "GC", "SCORES"))
  
  main.pvals<-main.table[,list(PVAL=mean(null.table[MUT==MUTATION & GC==GENE.COUNT,]$SCORES >= RATE)), by=c("MET", "MUTATION")]
  main.pvals$PVAL.ADJ<-p.adjust(main.pvals$PVAL, method = "fdr")
  
  #Clean up and Return
  return(main.pvals)
}

main.D.3.3.PVAL<-Function.D.3.3.PVAL(main.D.3.3$MET.TABLE[MUTATION!="SILENT",], D.3.3.NULL)

main.D.3.3.PVAL[order(PVAL.ADJ),]
ggplot(main.D.3.3.PVAL, aes(PVAL.ADJ)) + geom_histogram() + facet_wrap(~MUTATION)

##Find significant metabolite in CORE paper##
library(gdata)
nci.60.core<-as.matrix(read.table("DATABASES/CANCER_DATA/NCI.60/NCI60.CORE.S1.csv", header=T))
dim(nci.60.core)
nci.60.core[1:3, 1:3]
data.matrix(nci.60.core)[1:3, 1:3]

nci.doubling<-fread("DATABASES/CANCER_DATA/NCI.60/doubling.time.csv", header=T)
nci.doubling$Cell.line<-sapply(nci.doubling$Cell.line, function(x) gsub("-", ".", x))

nci.doubling[!(Cell.line %in% colnames(nci.60.core)),]
colnames(nci.60.core)
