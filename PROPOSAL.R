#PROPOSAL.R
library(data.table)
library(reshape2)
library(ggplot2)
library(gplots)
library(pheatmap)
library(made4)
library(igraph)

kegg.edges<-Function.kegg.filtered("DATABASES/KEGG/071415.ENZYME.SUB.PROD.MAIN.FILT", "DATABASES/RECON/042215.PROCESSED.METABOLITES", weight.filter = 60, n.edge.filter = 200)

table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
write.table(table.2, "PIPELINES/METABOLIC.DRIVERS/TABLES/092415.TABLE.2", quote = F, sep = "\t", row.names = F, col.names = F)

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

table.2.filt<-Function.kegg.table.filter(table.2, table.2 = T,gene.th = 40, degree.th = 100)

icgc.hub<-Function.icgc.enrich.coverage(icgc.brca.enrich.1, brca.icgc.mut, table.2.filt, pval.th=0.05, kegg.th=5, pval.col="PVAL.ADJ", gene.th=60, met.hub.th=0.8, table.2=T)

plot(icgc.hub$MET.HUBS.GRAPH)
heatmap.2(icgc.hub$KEGG.COV.MET, scale="none", trace = "none")
icgc.hub$KEGG.COV.SAMPLES
icgc.hub$HUB.COV.SAMPLES
heatmap.2(icgc.hub$HUB.COV.GENES$C02140.C05497.C05498, scale="none", trace="none", margins=c(8,8))
icgc.hub$ALL.COVERAGE

icgc.brca.path<-Function.kegg.path.enrich(kegg.path, unique(icgc.brca.enrich.1[PVAL.ADJ.2<0.05,]$KEGG.ID), table.2.filt, table.2 = T, gene.th = 40)
icgc.brca.path[PVAL.ADJ<0.05,]

table.2.dist<-Function.met.node.distance(table.2, table.2 = T, gene.th = 40, degree.th = 100)

icgc.brca.assign.distance<-Function.met.assign.distance(icgc.brca.enrich.1, table.2.dist$KEGG.DISTANCE, pval.th = 0.05, pval.column = "PVAL.ADJ")

ggplot(icgc.brca.assign.distance[EXP.MET.COUNT>=5,], aes(KEGG.ID, DIST)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

icgc.brca.enrich.distance<-Function.enrich.distance(icgc.brca.assign.distance, table.2.dist$KEGG.DISTANCE, exp.met.count.th = 5)
ggplot(icgc.brca.enrich.distance, aes(PVAL.ADJ)) + geom_histogram() + theme.format + geom_vline(xintercept=c(0.05), linetype="dashed", colour="red")
icgc.brca.enrich.distance


#
icgc.brca.assign.distance[KEGG.ID=="C02075",]
x<-table.2.dist$KEGG.GRAPH
V(x)$color<-ifelse(V(x)$name=="C05284", "red",
                   ifelse(V(x)$name %in% icgc.brca.assign.distance[KEGG.ID=="C05284",]$EXP.METS, "green", "yellow"))
par(mai=c(0,0,1,0))
plot(x, layout=layout.fruchterman.reingold, vertex.size=2, vertex.label.cex=0.3, edge.width=0.3)
dev.off()

y<-unique(icgc.brca.assign.distance[KEGG.ID=="C05284",]$EXP.METS)
z<-matrix(nrow=length(y), ncol=length(y), dimnames = list(y, y))
for (i in y){
  for (j in y){
    z[i,j]<-length(intersect(unique(table.2[KEGG_ID==i,]$GENE),unique(table.2[KEGG_ID==j,]$GENE)))/
      length(union(unique(table.2[KEGG_ID==i,]$GENE),unique(table.2[KEGG_ID==j,]$GENE)))
  }
}

heatmap.2(z, scale = "none", trace="none", margins=c(8,8))

#Use kegg.edges results and apply pipeline to it
icgc.brca.enrich.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/092115.ICGC.KEGG.IMPACT.rds")
kegg.edges.filt<-Function.kegg.table.filter(kegg.edges, table.2 = F,gene.th = 40, degree.th = 60)
icgc.brca.enrich.2$PVAL.ADJ.2<-p.adjust(icgc.brca.enrich.2$PVAL, method="fdr")
icgc.brca.enrich.2$PVAL.ADJ.3<-p.adjust(icgc.brca.enrich.2$PVAL, method="bonferroni")
icgc.brca.enrich.2[PVAL.ADJ<0.05,]
icgc.brca.enrich.2[PVAL.ADJ.2<0.05,]
icgc.brca.enrich.2[PVAL.ADJ.3<0.05,]

icgc.brca.enrich.combn<-icgc.brca.enrich.2[KEGG.ID %in% unique(table.2$KEGG_ID) & EXP.METS %in% unique(table.2$KEGG_ID),]
icgc.brca.enrich.combn[,PVAL:=p.adjust(PVAL, method="fdr"), by="KEGG.ID"]
icgc.brca.enrich.combn$PVAL.ADJ.2<-p.adjust(icgc.brca.enrich.combn$PVAL, method="fdr")
icgc.brca.enrich.combn$PVAL.ADJ.3<-p.adjust(icgc.brca.enrich.combn$PVAL, method="bonferroni")
icgc.brca.enrich.combn[PVAL.ADJ<0.05,]


icgc.hub.2<-Function.icgc.enrich.coverage(icgc.brca.enrich.combn, brca.icgc.mut, kegg.edges.filt, pval.th=0.05, kegg.th=5, pval.col="PVAL.ADJ",gene.th=40, met.hub.th=0.8, table.2=F)

plot(icgc.hub.2$MET.HUBS.GRAPH)
heatmap.2(icgc.hub.2$KEGG.COV.MET, scale="none", trace = "none")
icgc.hub.2$KEGG.COV.SAMPLES
icgc.hub.2$HUB.COV.SAMPLES
heatmap.2(icgc.hub.2$HUB.COV.GENES$C02934.C00655.C00077.C00065.C00181.C00090.C02266.C00674.C00570.C12448.C16834.C00364.C03150.C01272.C00184.C00475.C00307.C00311, scale="none", trace="none", margins=c(8,8))
icgc.hub.2$ALL.COVERAGE

icgc.brca.path.2<-Function.kegg.path.enrich(kegg.path, unique(icgc.brca.enrich.combn[PVAL.ADJ<0.05,]$KEGG.ID), kegg.edges.filt, table.2 = F, gene.th = 40)
icgc.brca.path.2[PVAL.ADJ<0.05,]

kegg.dist<-Function.met.node.distance(kegg.edges, table.2 = F, gene.th = 40, degree.th = 60)

icgc.brca.assign.distance.combn<-Function.met.assign.distance(icgc.brca.enrich.combn, kegg.dist$KEGG.DISTANCE, pval.th = 0.05, pval.column = "PVAL.ADJ")

ggplot(icgc.brca.assign.distance.combn[EXP.MET.COUNT>=5,], aes(KEGG.ID, DIST)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

icgc.brca.enrich.distance.combn<-Function.enrich.distance(icgc.brca.assign.distance.combn, kegg.dist$KEGG.DISTANCE, exp.met.count.th = 5)
ggplot(icgc.brca.enrich.distance.combn, aes(PVAL.ADJ)) + geom_histogram() + theme.format + geom_vline(xintercept=c(0.05), linetype="dashed", colour="red")

#NEED TO REPROCESS KEGG.EDGES
write.table(kegg.edges, "PIPELINES/METABOLIC.DRIVERS/TABLES/092315.KEGG.EDGES", quote = F, sep = "\t", row.names = F, col.names = F)
kegg.edges.refilt<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/092315.KEGG.EDGES.FILT", header=T) #done in python
kegg.edges.refilt<-Function.kegg.table.filter(kegg.edges.refilt, table.2 = T,gene.th = 40, degree.th = 60) #40, 60

length(unique(kegg.edges.refilt$KEGG_ID)) #897, 754
length(unique(kegg.edges.refilt$GENE)) #1313, 1175

icgc.brca.enrich.filtered<-icgc.brca.enrich.2[KEGG.ID %in% kegg.edges.refilt$KEGG_ID & EXP.METS %in% kegg.edges.refilt$KEGG_ID,][,c("KEGG.ID", "EXP.METS", "PVAL"),with=F]
icgc.brca.enrich.filtered[,PVAL.ADJ:=p.adjust(PVAL, method="fdr"), by="KEGG.ID"]
icgc.brca.enrich.filtered$PVAL.ADJ.2<-p.adjust(icgc.brca.enrich.filtered$PVAL, "fdr")

icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]
icgc.brca.enrich.filtered[PVAL.ADJ.2<0.05,]

icgc.hub.3<-Function.icgc.enrich.coverage(icgc.brca.enrich.filtered, brca.icgc.mut, kegg.edges.refilt, pval.th=0.05, kegg.th=5, pval.col="PVAL.ADJ",gene.th=40, met.hub.th=0.8, table.2=T)

table.2<-rbind(table.2, data.table(METABOLITE=c("S-Adenosylmethionine", "L-Cysteate", "Phenylacetic acid", "Benzenecarboxylic acid", "Butanoic acid", "5'-Deoxy-5-fluorocytidine"),
                                   KEGG_ID=c("C00019", "C00506", "C07086", "C00180", "C00246","C16635"),
                                   GENE=c(NA, NA, NA, NA, NA, NA)))
V(icgc.hub.3$MET.HUBS.GRAPH)$name<-sapply(V(icgc.hub.3$MET.HUBS.GRAPH)$name, function(x)  as.character(table.2[KEGG_ID==x,]$METABOLITE[1]))
plot(icgc.hub.3$MET.HUBS.GRAPH, layout=layout.fruchterman.reingold, vertex.size=6, vertex.label.cex=0.8, edge.width=0.8, edge.arrow.size=0.15)
table.2

heatmap.2(icgc.hub.3$KEGG.COV.MET, scale="none", trace = "none")
icgc.hub.3$KEGG.COV.SAMPLES
icgc.hub.3$HUB.COV.SAMPLES
heatmap.2(icgc.hub.3$HUB.COV.GENES$C06429.C00606.C00506.C01120.C00510.C05580.C00412.C00245.C00319.C00836.C07086, scale="none", trace="none", margins=c(8,8))
icgc.hub.3$ALL.COVERAGE

icgc.brca.path.3<-Function.kegg.path.enrich(kegg.path, unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID), kegg.edges.refilt, table.2 = T, gene.th = 40)
icgc.brca.path.3[PVAL.ADJ<0.05,]

kegg.refilt.dist<-Function.met.node.distance(kegg.edges.refilt, table.2 = T, gene.th = 40, degree.th = 60)

icgc.brca.assign.distance.refilt<-Function.met.assign.distance(icgc.brca.enrich.filtered, kegg.refilt.dist$KEGG.DISTANCE, pval.th = 0.05, pval.column = "PVAL.ADJ")

ggplot(icgc.brca.assign.distance.refilt[EXP.MET.COUNT>=5,], aes(KEGG.ID, DIST)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

icgc.brca.enrich.distance.refilt<-Function.enrich.distance(icgc.brca.assign.distance.refilt, kegg.refilt.dist$KEGG.DISTANCE, exp.met.count.th = 5)
ggplot(icgc.brca.enrich.distance.refilt, aes(PVAL.ADJ)) + geom_histogram() + theme.format + geom_vline(xintercept=c(0.05), linetype="dashed", colour="red")

icgc.brca.enrich.distance.refilt

#
icgc.brca.assign.distance.refilt[KEGG.ID=="C00095",]
x<-kegg.refilt.dist$KEGG.GRAPH
V(x)$color<-ifelse(V(x)$name=="C00095", "red",
                   ifelse(V(x)$name %in% icgc.brca.assign.distance.refilt[KEGG.ID=="C00095",]$EXP.METS, "green", "yellow"))
par(mai=c(0,0,1,0))
plot(x, layout=layout.fruchterman.reingold, vertex.size=2, vertex.label.cex=0.3, edge.width=0.3)
dev.off()

y<-unique(icgc.brca.assign.distance.refilt[KEGG.ID=="C00095",]$EXP.METS)
z<-matrix(nrow=length(y), ncol=length(y), dimnames = list(y, y))
for (i in y){
  for (j in y){
    z[i,j]<-length(intersect(unique(kegg.edges.refilt[KEGG_ID==i,]$GENE),unique(kegg.edges.refilt[KEGG_ID==j,]$GENE)))/
      length(union(unique(kegg.edges.refilt[KEGG_ID==i,]$GENE),unique(kegg.edges.refilt[KEGG_ID==j,]$GENE)))
  }
}

heatmap.2(z, scale = "none", trace="none", margins=c(8,8))

#Validate against Terunuma and Tang differentially expressed metabolites
length(unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID))

#TERU
teru.diff.met<-sapply(rownames(teru.cancer.matrix$MATRIX), function(x) {
  PVAL=wilcox.test(teru.cancer.matrix$MATRIX[x, ], teru.normal.matrix$MATRIX[x,])$p.value
})
teru.diff.met<-data.table(MET=rownames(teru.cancer.matrix$MATRIX), PVAL=teru.diff.met)
teru.diff.met<-teru.diff.met[!is.na(PVAL),]
teru.diff.met$PVAL.ADJ<-p.adjust(teru.diff.met$PVAL, method="fdr")
teru.diff.met[PVAL.ADJ<0.05,]

hist(teru.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID),]$PVAL.ADJ)
binom.test( 23 ,24,p = 196/241, alternative = "greater") #KEGG.ID
binom.test( 8,9,p = 196/241, alternative = "greater") #EXP.METS
phyper(23-1,196, 241-196, 24, lower.tail = F) #KEGG.ID
phyper(8-1, 196, 241-196, 9, lower.tail = F) #EXP.METS

length(teru.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ)
sum(teru.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ<0.05)
binom.test(20,22,p = 196/241, alternative = "greater") #EXP.METS

#TANG
colnames(tang.matrix)<-colnames(data.frame(tang.matrix))
tang.diff.met<-apply(tang.matrix, 1,  function(x) {
  normal=c("NORMAL", "NORMAL.1", "NORMAL.2", "NORMAL.3", "NORMAL.4")
  cancer=setdiff(colnames(tang.matrix), normal)
  PVAL=wilcox.test(x[normal], x[cancer])$p.value
})
tang.diff.met<-data.table(MET=rownames(tang.matrix), PVAL=tang.diff.met)
tang.diff.met<-tang.diff.met[!is.na(PVAL),]
tang.diff.met$PVAL.ADJ<-p.adjust(tang.diff.met$PVAL, method="fdr")
tang.diff.met[PVAL.ADJ<0.05,]

hist(tang.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID),]$PVAL.ADJ)
binom.test(16, 22,143/207, alternative = "greater") #KEGG.ID
binom.test(8, 9,143/207, alternative = "greater") #EXP.METS

length(tang.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID),]$PVAL.ADJ)
sum(tang.diff.met[MET %in% unique(icgc.brca.enrich.filtered[PVAL.ADJ<0.05,]$KEGG.ID),]$PVAL.ADJ<0.05)
binom.test(25, 34,143/207, alternative = "greater") #EXP.METS

#Previous proposal validation with tang and terunuma
teru.diff.met[MET %in% c("C01194", "C05489", "C00518", "C07490", "C05299", "C11136", "C03033", "C14869", "C05302", "C01194", "C05981", "C00334", "C00404"),]
tang.diff.met[MET %in% c("C01194", "C05489", "C00518", "C07490", "C05299", "C11136", "C03033", "C14869", "C05302", "C01194", "C05981", "C00334", "C00404"),]

teru.diff.met[MET %in% c("C02297", "C02714"),]
tang.diff.met[MET %in% c("C02297", "C02714"),]

x<-fread("DATABASES/KEGG/063015.ENZYME.SUB.PROD", header=T)
length(union(x$SUBSTRATE, x$PRODUCT))
length(unique(x$Hugo_Symbol))

######APPLY METHOD TO TANG DATASET########
tang.brca.enrich.1<-Function.icgc.enrich.1(brca.icgc.mut[SAMPLE %in% colnames(tang.matrix),], brca.exp, kegg.edges.refilt, n.samples = 2)
hist(tang.brca.enrich.1$PVAL.ADJ)
tang.brca.enrich.1[PVAL.ADJ<0.1,]

tang.diff.met[MET %in% unique(tang.brca.enrich.1[PVAL.ADJ<0.1,]$EXP.METS),]

icgc.brca.enrich.filtered[PVAL.ADJ<0.05,][,list(N.MET=length(unique(EXP.METS))), by="KEGG.ID"][order(N.MET),]
x<-unique(kegg.edges.refilt[GENE %in% brca.icgc.mut[MUTATION=="MISSENSE",]$Hugo_Symbol,][KEGG_ID %in% icgc.brca.enrich.filtered[PVAL.ADJ<0.05,][,list(N.MET=length(unique(EXP.METS))), by="KEGG.ID"][order(N.MET),]$KEGG.ID,]$GENE)
ccg.genes<-c("AKT1", "ARID1A", "ARID1B", "BAP1", "BRCA1", "BRCA2", "BRIP1", "CASP8", "CCND1", "CDH1", "CDKN1B", "CHEK2", "EP300","ERBB2", "ESR1", "ETV6",
  "FOXA1", "GATA3", "MAP2K4", "MAP3K1", "MAP3K13", "NCOR1", "NOTCH1", "NTRK3", "PALB2", "PBRM1", "RB1", "SMARCD1","TBX3", "TP53", "PIK3CA") 
sum(ccg.genes %in% table.2.refilt$Hugo_Symbol)

######NEED TO REFILT TABLE.2 SO WE CAN HAVE GENES ASSOCIATED WITH CANCER PROGRESSION#######
#Genes in latest original version of kegg.edges do belong to CCG list for breast cancer

table.2.refilt<-fread("PIPELINES/METABOLIC.DRIVERS/TABLES/092415.TABLE.2.FILT", header=T) #Done in python
length(unique(table.2.refilt$KEGG_ID))
length(unique(table.2.refilt$GENE))
table.2.refilt[KEGG_ID=="C00005",]
table.2.refilt[KEGG_ID=="C00037",]
table.2.refilt[GENE=="PIK3CA",]


table.2.refilt<-Function.kegg.table.filter(table.2.refilt, table.2 = T,gene.th = 195, degree.th = 285) 
icgc.brca.enrich.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/092215.icgc.brca.enrich.1.table.2.rds")


icgc.brca.enrich.1.table.2.refilt<-icgc.brca.enrich.1[KEGG.ID %in% table.2.refilt$KEGG_ID & EXP.METS %in% table.2.refilt$KEGG_ID,][,c("KEGG.ID", "EXP.METS", "PVAL"),with=F]
icgc.brca.enrich.1.table.2.refilt[,PVAL.ADJ:=p.adjust(PVAL, method="fdr"), "KEGG.ID"]
icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]


icgc.hub.table.2.refilt<-Function.icgc.enrich.coverage(icgc.brca.enrich.1.table.2.refilt, 
                                                       brca.icgc.mut, table.2.refilt, 
                                                       pval.th=0.05, kegg.th=3, pval.col="PVAL.ADJ", 
                                                       gene.th=500, met.hub.th=0.7, table.2=T)

heatmap.2(icgc.hub.table.2.refilt$KEGG.COV.MET, scale="none", trace = "none")
icgc.hub.table.2.refilt$KEGG.COV.SAMPLES
plot(icgc.hub.table.2.refilt$MET.HUBS.GRAPH)
icgc.hub.table.2.refilt$HUB.COV.SAMPLES
heatmap.2(icgc.hub.table.2.refilt$HUB.COV.GENES$C01762.C00387.C00559.C00239, scale="none", trace="none", margins=c(8,8))
icgc.hub.table.2.refilt$ALL.COVERAGE




hist(teru.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ)
length(teru.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ)
sum(teru.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ<0.05)
binom.test(29, 31, 196/241)

hist(tang.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ)
length(tang.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ)
sum(tang.diff.met[MET %in% unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS),]$PVAL.ADJ<0.05)
binom.test(23,30, 143/207)

method.ccg.1<-Function.cumulative.ccg(icgc.brca.enrich.1.table.2.refilt, table.2.refilt, pval.th = 0.05, table.2 = T)
ggplot(method.ccg.1$EXP.MET.CCG, aes(LEVEL, CCG.RATIO)) + geom_point() + geom_line()
method.ccg.1$EXP.MET.CCG
icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,][,list(N=length(unique(KEGG.ID))), by="EXP.METS"][order(N, decreasing = T),]

length(unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$KEGG.ID))
length(unique(icgc.brca.enrich.1.table.2.refilt[PVAL.ADJ<0.05,]$EXP.METS))

grid.search.1<-Function.grid.search.1()
grid.search.1[TERU.PVAL<0.05 & TANG.PVAL<0.05,]
grid.search.1[TERU.PVAL<0.1,]
grid.search.1[order(TERU.PVAL, TANG.PVAL)]
grid.search.1[order(TANG.PVAL, TERU.PVAL)]

hist(grid.search.1$TERU.PVAL)
grid.search.1[TERU.PVAL<0.1,][C01194==T & C05981==T,][,list(GENE.TH=GENE.TH[1], DEGREE.TH=DEGREE.TH[1], SIG.KEGG.ID=SIG.KEGG.ID[1], SIG.EXP.METS=SIG.EXP.METS[1]), by="TERU.PVAL"]
