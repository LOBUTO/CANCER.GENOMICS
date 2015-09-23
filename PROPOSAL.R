#PROPOSAL.R
library(data.table)
library(reshape2)
library(ggplot2)
library(gplots)

kegg.edges<-Function.kegg.filtered("DATABASES/KEGG/071415.ENZYME.SUB.PROD.MAIN.FILT", "DATABASES/RECON/042215.PROCESSED.METABOLITES", weight.filter = 60, n.edge.filter = 60)

table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")

gene.length<-fread("DATABASES/UNIPROT/042314_GENE_LENGTH", header=T, sep="\t", stringsAsFactors = F)

tcga.mut<-Function.tcga.mut.prep("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", mut.filter=1, combine = F) #at least 5 people have all mutations
length(unique(tcga.mut[Hugo_Symbol %in% unique(kegg.edges$Hugo_Symbol), ]$SAMPLE)) #833@1, 726@5

#
main.D.3.3<-Function.D.3.3(tcga.mut, enz.prod , gene.length, enz.prod = T)
main.D.3.3$MET.TABLE[order(PVAL.ADJ),]

ggplot(main.D.3.3$MET.TABLE, aes(GENE.COUNT, RATE)) + geom_point() + facet_wrap(~MUTATION) + scale_x_log10() #Have to draw from same background
ggplot(main.D.3.3$MET.TABLE, aes(PVAL.ADJ)) + geom_histogram() + facet_wrap(~MUTATION) 
ggplot(main.D.3.3$MET.TABLE, aes(MUTATION,RATE)) + geom_boxplot()

#Calculate NULL
Function.D.3.3.NULL<-function(main.table, gene.length, tcga.mut){

  #Filter tcga to correct pools
  tcga.pool<-merge(tcga.mut, gene.length, by="Hugo_Symbol")
  setkey(tcga.pool)
  tcga.pool<-unique(tcga.pool[,c("SAMPLE", "Hugo_Symbol", "MUTATION"), with=F])
  #tcga.pool<-tcga.pool[Hugo_Symbol %in% unique(kegg.edges$Hugo_Symbol),]
  
  #Mutation types
  types<-unique(main.table$MUTATION)
  
  #Mutation count per gene for tcga
  tcga.mut<-tcga.mut[,list(N.MUT=length(SAMPLE)), by=c("Hugo_Symbol", "MUTATION")]
  all.filt.samples<-length(unique(merge(merge(kegg.edges, tcga.mut, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")[["SAMPLE"]]))
  
  #Obtain null distribution
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
      SCORES=sapply(random.pool, function(x) (sum(tcga.mut[MUTATION==t,][Hugo_Symbol %in% x,]$N.MUT)/(sum(gene.length[Hugo_Symbol %in% x,]$LENGTH)))*(1/all.filt.samples))
      
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

main.D.3.3.PVAL[order(PVAL),]
ggplot(main.D.3.3.PVAL, aes(PVAL.ADJ)) + geom_histogram() + facet_wrap(~MUTATION)

##Find significant metabolite in CORE paper##
library(gdata)
nci.60.core<-as.matrix(read.csv("DATABASES/CANCER_DATA/NCI.60/NCI60.CORE.S1.csv", header=F))
dim(nci.60.core)

colnames(nci.60.core)<-nci.60.core[1,]
met.names<-nci.60.core[,1]
nci.60.core<-nci.60.core[2:nrow(nci.60.core),2:ncol(nci.60.core)]
nci.60.core<-apply(nci.60.core, 2, as.numeric)
rownames(nci.60.core)<-met.names[2:length(met.names)]
nci.60.core[1:3, 1:3]

nci.doubling<-fread("DATABASES/CANCER_DATA/NCI.60/doubling.time.csv", header=T)

nci.60.cor<-Function.nci60.cor(nci.60.core, nci.doubling)
nci.60.cor[PVAL.ADJ<0.01,]

nci.60.cor.cancer<-Function.nci60.cor.cancers(nci.60.core, nci.doubling)
nci.60.cor.cancer[order(PVAL.ADJ),]
nci.60.cor.cancer[PVAL.ADJ<0.05,]
nci.60.cor.cancer[CANCER=="Colon",][order(PVAL.ADJ),]
nci.60.cor.cancer[PVAL.ADJ<0.1,]
nci.60.cor.cancer[CANCER=="Breast",][order(PVAL.ADJ),]

#Use method across stages to find significantly mutated metabolites
TCGA.BRCA.CLINICAL<-Function.tcga.brca.clinical("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
tcga.mut

D.3.STAGES<-data.table()
for (stage in unique(TCGA.BRCA.CLINICAL$STAGE)){
  print(stage)
  
  #Obtain stage samples
  stage.samples<-unique(TCGA.BRCA.CLINICAL[STAGE==stage,]$SAMPLE)
  stage.samples<-intersect(stage.samples, tcga.mut$SAMPLE)
  
  stage.table<-Function.D.3.3(tcga.mut[SAMPLE %in% stage.samples,], enz.prod , gene.length, enz.prod = T)
  stage.table$MET.TABLE$STAGE<-stage
  D.3.STAGES<-rbind(D.3.STAGES, stage.table$MET.TABLE)
}
ggplot(D.3.STAGES[MUTATION!="SILENT",][STAGE!="Stage.IV",], aes(RATE, colour=STAGE)) + geom_density(position = "dodge") + facet_wrap(~MUTATION) + scale_y_sqrt()

D.3.STAGES$STAGE<-ifelse(D.3.STAGES$STAGE=="Stage.I", 1, 
                         ifelse(D.3.STAGES$STAGE=="Stage.II", 2,
                                ifelse(D.3.STAGES$STAGE=="Stage.III", 3,4)))

D.3.STAGES<-D.3.STAGES[STAGE!=4,]
D.3.STAGES<-data.table(acast(D.3.STAGES, MET + MUTATION~STAGE, fill = 0, value.var = "RATE" ), keep.rownames = T)
D.3.STAGES$MET<-sapply(D.3.STAGES$rn, function(x) strsplit(x, "_")[[1]][1] )
D.3.STAGES$MUTATION<-sapply(D.3.STAGES$rn, function(x) strsplit(x, "_")[[1]][2] )
D.3.STAGES$rn<-NULL

D.3.STAGES<-D.3.STAGES[apply(D.3.STAGES[,1:3,with=F], 1, function(x)  mean(c(x[1]>0,x[2]>0,x[3]>0)))>0.66,]
D.3.STAGES<-D.3.STAGES[apply(D.3.STAGES[,1:3,with=F], 1, function(x) sd(x))!=0,]

D.3.STAGES<-melt(D.3.STAGES, id.vars = c("MET", "MUTATION"))
setnames(D.3.STAGES, c("MET", "MUTATION", "STAGE", "RATE"))
D.3.STAGES$STAGE<-as.numeric(D.3.STAGES$STAGE)

ggplot(D.3.STAGES, aes(factor(STAGE), RATE)) + geom_boxplot() + facet_wrap(~MUTATION)

D.3.STAGES<-D.3.STAGES[,list(COR=cor.test(STAGE, RATE, method="pearson")$estimate,
                             PVAL=cor.test(STAGE, RATE, method="pearson")$p.value), by=c("MET", "MUTATION") ]
D.3.STAGES[,PVAL.ADJ:=p.adjust(PVAL, method="fdr"), by="MUTATION"]
D.3.STAGES[order(PVAL.ADJ),]
D.3.STAGES[order(PVAL),]

hist(D.3.STAGES$PVAL.ADJ)
D.3.STAGES[PVAL.ADJ<0.05,]

D.3.STAGES[MET=="C00164",]

#Use combinatorial method
brca.exp<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041715.BRCA.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
brca.icgc.mut<-Function.icgc.mut.prep("DATABASES/CANCER_DATA/ICGC/BRCA/ssm.tsv")
table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
kegg.path<-fread("DATABASES/KEGG/060415_PATHWAY_TO_COMPOUND", header=T, sep="\t")
saveRDS(brca.icgc.mut, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/092115.BRCA.ICGC.MUT.rds")
saveRDS(kegg.edges, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/092115.KEGG.EDGES.rds")

icgc.brca.enrich.1<-Function.icgc.enrich.1(brca.icgc.mut, brca.exp, table.2)
icgc.brca.enrich.1[PVAL.ADJ<0.05,]
icgc.brca.enrich.1$PVAL.ADJ.2<-p.adjust(icgc.brca.enrich.1$PVAL, method="fdr")
icgc.brca.enrich.1$PVAL.ADJ.3<-p.adjust(icgc.brca.enrich.1$PVAL, method="bonferroni")
icgc.brca.enrich.1[PVAL.ADJ.2<0.05,]
icgc.brca.enrich.1[PVAL.ADJ.3<0.05,]
saveRDS(icgc.brca.enrich.1, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/092215.icgc.brca.enrich.1.table.2.rds")

length(unique(icgc.brca.enrich.1[PVAL.ADJ<0.05,][KEGG.ID %in% table.2[,list(N=length(unique(GENE))), by="KEGG_ID"][N<40,]$KEGG_ID,]$KEGG.ID))

icgc.brca.enrich.1[PVAL.ADJ<0.05,]

icgc.hub<-Function.icgc.enrich.coverage(icgc.brca.enrich.1, brca.icgc.mut, table.2, pval.th = 0.05, kegg.th = 4, pval.col="PVAL.ADJ", gene.th=165, met.hub.th = 0.8, table.2=T)

plot(icgc.hub$MET.HUBS.GRAPH)
heatmap.2(icgc.hub$KEGG.COV.MET, scale="none", trace = "none")
icgc.hub$KEGG.COV.SAMPLES
icgc.hub$HUB.COV.SAMPLES
icgc.hub$HUB.COV.GENES$C00390.C00399.C00721.C01134

icgc.brca.path<-Function.kegg.path.enrich(kegg.path, unique(icgc.brca.enrich.1[PVAL.ADJ.2<0.05,]$KEGG.ID), table.2, table.2 = T, th = 40)
icgc.brca.path[PVAL.ADJ<0.05,]

table.2.dist<-Function.met.node.distance(table.2, table.2 = T, th = 40, degree.th = 60)

icgc.brca.assign.distance<-Function.met.assign.distance(icgc.brca.enrich.1, table.2.dist$KEGG.DISTANCE, pval.th = 0.05, pval.column = "PVAL.ADJ")

ggplot(icgc.brca.assign.distance[EXP.MET.COUNT>=5,], aes(KEGG.ID, DIST)) + geom_boxplot()

icgc.brca.enrich.distance<-Function.enrich.distance(icgc.brca.assign.distance, table.2.dist$KEGG.DISTANCE, exp.met.count.th = 5)
hist(icgc.brca.enrich.distance$PVAL.ADJ)
icgc.brca.enrich.distance

icgc.brca.assign.distance[KEGG.ID=="C00236",]
x<-table.2.dist$KEGG.GRAPH
V(x)$color<-ifelse(V(x)$name=="C05498", "red",
                   ifelse(V(x)$name %in% icgc.brca.assign.distance[KEGG.ID=="C05498",]$EXP.METS, "green", "yellow"))
plot(x, layout=layout.fruchterman.reingold.grid, vertex.size=2, vertex.label.cex=0.3)

y<-unique(icgc.brca.assign.distance[KEGG.ID=="C05498",]$EXP.METS)
z<-matrix(nrow=length(y), ncol=length(y), dimnames = list(y, y))
for (i in y){
  for (j in y){
    z[i,j]<-length(intersect(unique(table.2[KEGG_ID==i,]$GENE),unique(table.2[KEGG_ID==j,]$GENE)))/
      length(union(unique(table.2[KEGG_ID==i,]$GENE),unique(table.2[KEGG_ID==j,]$GENE)))
  }
}
