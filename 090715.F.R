Function.Feature.Met.Table<-function(kegg.edges, met.cancer.matrix, met.normal.matrix, teru.gene.matrix, kegg.path){
  #Build feature table for metabolites based on direct gene expression data for Terunuma et al. dataset
  
  #Get common cancer and normal patients for gene and met data
  normal.samples<-intersect(colnames(met.normal.matrix), colnames(teru.gene.matrix))
  cancer.samples<-intersect(colnames(met.cancer.matrix), colnames(teru.gene.matrix))
  
  #Normalize expression matrices to z-scores
  teru.gene.matrix<-scale(teru.gene.matrix, scale = T, center = T)
  
  #Filter metabolite matrix by those that have gene associated to it in both substrate and product and that we have gene expression information
  kegg.edges<-kegg.edges[Hugo_Symbol %in% rownames(teru.gene.matrix),]
  keggs<-rownames(met.cancer.matrix)
  keggs.filt<-sapply(keggs, function(x) ((x %in% kegg.edges$SUBSTRATE) + (x %in% kegg.edges$PRODUCT))==2)
  keggs<-keggs[keggs.filt]
  
  #Filter gene and met matrices for common samples
  met.cancer.matrix<-met.cancer.matrix[, cancer.samples]
  met.normal.matrix<-met.normal.matrix[, normal.samples]
  teru.gene.matrix<-teru.gene.matrix[, c(normal.samples, cancer.samples)]
  
  #Build feature table for all filtered mets
  all.mets<-unique(keggs)
  main.table<-data.table()
  for (met in all.mets){
    print (met)
    
    #Obtain met normal median
    normal.median<-median(met.normal.matrix[met,])
    
    #Obtain met log2 to median normal
    met.table<-data.table(t(log2( cbind(met.cancer.matrix[met,, drop=F], met.normal.matrix[met,,drop=F]) / normal.median)), keep.rownames = T)
    setnames(met.table, c("SAMPLE", "MET.LEVELS"))
    met.table$NAME<-met
    
    #Find correlation to non-related genes and choose TOP 10 genes (5 Negative and 5 Positvely correlated)
    met.genes<-union(kegg.edges[SUBSTRATE==met,]$Hugo_Symbol, kegg.edges[PRODUCT==met,]$Hugo_Symbol)
    non.met.genes<-setdiff(rownames(teru.gene.matrix), met.genes)
    
    cor.table<-data.table(GENE=non.met.genes, COR=sapply(non.met.genes, function(x) cor(met.table$MET.LEVELS, teru.gene.matrix[x, met.table$SAMPLE])))
    top.pos.genes<-cor.table[order(COR, decreasing = T),]$GENE[1:5]
    top.neg.genes<-cor.table[order(COR, decreasing = F),]$GENE[1:5]
    
    #Construct table
    sub.genes<-unique(kegg.edges[SUBSTRATE==met,]$Hugo_Symbol)
    prod.genes<-unique(kegg.edges[PRODUCT==met,]$Hugo_Symbol)
    
    sub.table<-data.table(t(t(colMeans(teru.gene.matrix[sub.genes,,drop=F]))), keep.rownames = T)
    setnames(sub.table, c("SAMPLE", "MEAN.SUB"))
    prod.table<-data.table(t(t(colMeans(teru.gene.matrix[prod.genes,, drop=F]))), keep.rownames = T)
    setnames(prod.table, c("SAMPLE", "MEAN.PROD"))
    
    top.pos.table<-data.table(t(teru.gene.matrix[top.pos.genes,]), keep.rownames = T)
    setnames(top.pos.table, c("SAMPLE", letters[1:5]))
    top.neg.table<-data.table(t(teru.gene.matrix[top.neg.genes,]), keep.rownames = T)
    setnames(top.neg.table, c("SAMPLE", letters[6:10]))
    
    top.table<-merge(top.pos.table, top.neg.table, by="SAMPLE")
    gene.table<-merge(merge(sub.table, prod.table, by="SAMPLE"), top.table, by="SAMPLE")
    main.met.table<-merge(met.table, gene.table, by="SAMPLE")
    
    #Store
    main.table<-rbind(main.table, main.met.table)
  }
  
  #Add path info
  main.table$GLYCO.TCA.LIPID<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Glycol", DESCRIPTION),]$COMPOUND, 1,
                                               ifelse(main.table$NAME %in% kegg.path[grepl("fatty acid", DESCRIPTION, ignore.case = T),]$COMPOUND, 1,
                                                      ifelse(main.table$NAME %in% kegg.path[grepl("TCA", DESCRIPTION),]$COMPOUND, 1, 0))))
  main.table$CANCER<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("cancer", DESCRIPTION, ignore.case = T),]$COMPOUND, 1, 
                                      ifelse(main.table$NAME %in% kegg.path[grepl("Glioma", DESCRIPTION),]$COMPOUND, 1, 
                                             ifelse(main.table$NAME %in% kegg.path[grepl("Melanoma", DESCRIPTION),]$COMPOUND, 1, 0))))
  main.table$STEROID<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Steroid", DESCRIPTION, ignore.case = T),]$COMPOUND, 1, 0))
  main.table$ABC<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("ABC transporters", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ANTIBIO<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Biosynthesis of antibiotics", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$AA<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Biosynthesis of amino acids", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$PROTEIN.ABS<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Protein digestion and absorption", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$AMINOACYL.TRNA<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Aminoacyl-tRNA biosynthesis", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$TWO.OXO<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("2-Oxocarboxylic acid metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$PURINE<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Purine metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$GST<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Glycine, serine and threonine", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$CARBON<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Carbon metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ARG.PRO<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Arginine and proline metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ALA.ASP.GLU<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Alanine, aspartate and glutamate metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$CYS.MET<-as.factor(ifelse(main.table$NAME %in% kegg.path[grepl("Cysteine and methionine metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  
  #Return
  return(main.table)
}

normalize.vector<-function(x){
  y=(x-min(x))/(max(x)- min(x))
  return(y)
}

Function.Feat.Table.Std<-function(met.feat, normalize=F){
  #Takes output from Function.Feature.Met.Table() and standarizes metabolite
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  met.feat$MET.LEVELS<-as.factor(ifelse(abs(met.feat$MET.LEVELS)>1, 1, 0))
  
  #Normalize?
  if (normalize==T){
    norm.col<-c("MEAN.SUB", "MEAN.PROD", colnames(met.feat)[colnames(met.feat) %in% letters]  )
    
    for (n in norm.col){
      met.feat[[n]]<-normalize.vector(met.feat[[n]])
    }
  }
  
  #Return
  return(met.feat)
}

Function.nci60.cor<-function(nci.60.core, nci.doubling){
  
  common.lines<-intersect(colnames(nci.60.core), nci.doubling$Cell.line)
  common.lines<-common.lines[order(common.lines,decreasing = T)]
  common.index<-which(colnames(nci.60.core) %in% common.lines)
  
  nci.60.pvalue<-apply(nci.60.core, 1, function(x) 
    cor.test(x[common.index], sapply(colnames(nci.60.core[,common.index]), function(y) nci.doubling[Cell.line==y,]$Doubling.time) , method = "pearson")$p.value)
  nci.60.rho<-apply(nci.60.core, 1, function(x) 
    cor.test(x[common.index], sapply(colnames(nci.60.core[,common.index]), function(y) nci.doubling[Cell.line==y,]$Doubling.time) , method = "pearson")$estimate)
  nci.60.cor<-data.table(MET=names(nci.60.pvalue), COR=nci.60.rho, PVAL=nci.60.pvalue)
  nci.60.cor$PVAL.ADJ<-p.adjust(nci.60.cor$PVAL, method = "fdr")
  
  nci.60.cor<-nci.60.cor[order(PVAL.ADJ),]  
  return(nci.60.cor)
}

Function.nci60.cor.cancers<-function(nci.60.core, nci.doubling){
  
  nci.cancers<-unique(nci.doubling[Cell.line!="MDA-N" & Cancer!="Prostate",]$Cancer)
  
  main.table<-data.table()
  for (cancer in nci.cancers){
    
    common.lines<-intersect(colnames(nci.60.core), nci.doubling[Cancer==cancer,]$Cell.line)
    common.lines<-common.lines[order(common.lines,decreasing = T)]
    common.index<-which(colnames(nci.60.core) %in% common.lines)
    doubling.vector<-sapply(colnames(nci.60.core[,common.index]), function(y) nci.doubling[Cell.line==y,]$Doubling.time)
    print (c(cancer, common.lines))
    
    nci.60.pvalue<-apply(nci.60.core, 1, function(x) 
      cor.test(x[common.index], doubling.vector, method = "pearson")$p.value)
    nci.60.rho<-apply(nci.60.core, 1, function(x) 
      cor.test(x[common.index], doubling.vector, method = "pearson")$estimate)
    nci.60.main<-data.table(MET=names(nci.60.pvalue), COR=nci.60.rho, PVAL=nci.60.pvalue)
    nci.60.main$PVAL.ADJ<-p.adjust(nci.60.main$PVAL, method = "fdr")
    
    nci.60.main<-nci.60.main[order(PVAL.ADJ),]  
    nci.60.main$CANCER<-cancer
    main.table<-rbind(main.table, nci.60.main)
  }
  
  return(main.table)
}

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
  
  #Obtain posible mutation length per metabolite
  met.length<-merge(kegg, gene.length, by="Hugo_Symbol")
  setkey(met.length)
  met.length<-unique(met.length)
  met.length<-met.length[,list(MET.LENGTH=sum(LENGTH)), by="MET"]
  
  #Obtain total mutation length for tcga sample
  tcga.length<-merge(tcga.mut, gene.length, by="Hugo_Symbol")
  tcga.length<-sum(unique(tcga.length[,c("Hugo_Symbol", "LENGTH"),with=F]$LENGTH))
  
  #Calculate individuals mutation count and total mutation length
  main.table<-merge(merge(kegg, tcga.mut, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")
  main.table[,SAMPLE.N.MUT:=length(Hugo_Symbol), by=c("SAMPLE")]
  
#   sample.length<-main.table[,c("Hugo_Symbol", "SAMPLE", "LENGTH"), with=F]
#   setkey(sample.length)
#   sample.length<-unique(sample.length)
#   sample.length<-sample.length[,list(SAMPLE.LENGTH=sum(LENGTH)), by="SAMPLE"]
#  main.table<-merge(main.table, sample.length, by="SAMPLE")
  main.table<-merge(main.table, met.length, by="MET")
  print (main.table)
  
  main.table<-main.table[,list(PVAL=binom.test(length(Hugo_Symbol), sum(SAMPLE.N.MUT), unique(MET.LENGTH)/tcga.length, alternative = "greater" )$p.value), 
                         by=c("MET", "MUTATION")]
  
  #Obtain mutation rate per metabolite
  #   all.filt.samples<-length(unique(merge(merge(kegg, tcga.mut, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")[["SAMPLE"]]))
  #   
  #   tcga.mut<-tcga.mut[,list(N.MUT=length(SAMPLE)), by=c("Hugo_Symbol", "MUTATION")]
  #   main.table<-merge(merge(kegg, tcga.mut, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")
  #   filtered.genes<-unique(main.table[,c("Hugo_Symbol","MUTATION"), with=F])
  #   main.table<-main.table[, list(RATE=   (sum(N.MUT)/sum(LENGTH)) * (1/all.filt.samples), GENE.COUNT=length(unique(Hugo_Symbol))), by=c("MET", "MUTATION")]
  #   
  #   #Clean up and Return
  #   main.table<-main.table[!is.na(RATE),]
  
  #Correct each metabolic case for FDR
  #main.table[,PVAL.ADJ:=p.adjust(PVAL, method="fdr"), by=c("MET", "MUTATION")]
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
  
  #Obtain score
  #main.table<-main.table[,list(RATE=mean(PVAL.ADJ<0.05)), by=c("MET", "MUTATION")]
  
  #return(list(MET.TABLE=main.table, GENE.TABLE=filtered.genes))
  return(list(MET.TABLE=main.table))
}

Function.Met.TCGA.1000G<-function(tcga.mut, kegg.edges, thousand.mut, gene.length, enz.prod=F){
  
  #Get associated metabolites per gene 
  if (enz.prod==T){
    kegg<-unique(kegg.edges[,c("Enzyme", "KEGG_ID"), with=F])
    setnames(kegg, c("Hugo_Symbol", "MET"))
  } else {
    kegg<-rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT))
    kegg<-unique(kegg) 
  }
  
  #Obtain mutation classification for thousand genomes
  t<-thousand.mut[,c("Hugo_Symbol", "TOTAL.ALLELES", "ALT.ALLELES", "Variant_Class"),with=F][,list(N.RATE=sum(ALT.ALLELES)/sum(TOTAL.ALLELES)), by=c("Hugo_Symbol", "Variant_Class")]
  setnames(t, c("Hugo_Symbol", "MUTATION", "THOUSAND.RATE"))
  t$MUTATION<-ifelse(t$MUTATION=="synonymous", "SILENT",
                     ifelse(t$MUTATION=="nonsynonymous", "GOF", "LOF"))
  
  #Calculate mutation rates in tcga per gene
  n.tcga<-length(unique(tcga.mut[MUTATION=="GOF",]$SAMPLE))
  m<-tcga.mut[MUTATION=="GOF",][,list(TCGA.RATE=length(SAMPLE)/n.tcga), by=c("Hugo_Symbol", "MUTATION")]
  
  #Combine metabolic and length info info
  m<-merge(merge(m, kegg, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")
  t<-merge(merge(t, kegg, by="Hugo_Symbol"), gene.length, by="Hugo_Symbol")
  m$TCGA.RATE<-m$TCGA.RATE/m$LENGTH
  t$THOUSAND.RATE<-t$THOUSAND.RATE/t$LENGTH
  
  #Obtain all metabolites that have at least 3 gene associations in both tcga and thousand panels
  t.met.info<-t[,list(N=length(Hugo_Symbol)), by=c("MUTATION", "MET")]
  m.met.info<-m[,list(N=length(Hugo_Symbol)), by=c("MUTATION", "MET")]
  t.met.info<-t.met.info[N>=3,]
  m.met.info<-m.met.info[N>=3,]
  met.info<-merge(m.met.info, t.met.info, by=c("MUTATION", "MET"))
  met.info$STATUS<-TRUE
  met.info<-met.info[,c("MUTATION", "MET", "STATUS"),with=F]
  
  #Filter both tables and perform statistic
  t<-merge(t, met.info, by=c("MUTATION", "MET"))
  m<-merge(m, met.info, by=c("MUTATION", "MET"))
  setnames(t, c("MUTATION.T", "MET.T", "Hugo_Symbol", "THOUSAND.RATE", "LENGTH", "STATUS"))
  
  main.table<-m[,list(PVAL=t.test(TCGA.RATE, t[MUTATION.T==MUTATION & MET.T==MET, ]$THOUSAND.RATE, alternative = "greater", paired = F)$p.value), by=c("MUTATION", "MET")]
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
  main.table[order(PVAL.ADJ),]
  main.table[order(PVAL)]
}

Function.icgc.mut.prep<-function(ssm.file){
  #Functions to process icgc somatic mutations for missense and synonymous variants
  #NOTE: One isoform per gene
  
  require(biomaRt)
  require(data.table)
  require(reshape2)
  
  #Clean up icgc table
  brca.icgc.mut<-fread(ssm.file, header=T)
  brca.icgc.mut<-brca.icgc.mut[consequence_type %in% c("missense_variant", "synonymous_variant"),]
  brca.icgc.mut<-brca.icgc.mut[,c("submitted_sample_id", "aa_mutation", "gene_affected", "consequence_type"), with=F]
  brca.icgc.mut$SAMPLE<-sapply(brca.icgc.mut$submitted_sample_id, function(x) paste(strsplit(x, "-")[[1]][1:4], collapse = "."))
  
  #Replace ensembl with hugo identifiers
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  test.genes<-unique(brca.icgc.mut$gene_affected)
  G_list <- data.table(getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=test.genes, mart= mart))
  setnames(G_list, c("ensembl", "Hugo_Symbol"))
  G_list<-G_list[Hugo_Symbol!="",]
  
  #Clean up
  setnames(brca.icgc.mut, c("old", "AA.SITE", "ensembl","MUTATION", "SAMPLE"))
  main.table<-merge(brca.icgc.mut, G_list, by="ensembl")
  main.table<-main.table[,c("SAMPLE", "Hugo_Symbol", "AA.SITE", "MUTATION"), with=F]
  main.table$MUTATION<-ifelse(main.table$MUTATION=="missense_variant", "MISSENSE", "SILENT")
  
  #Watch out for multiple replicates of transcripts within the same patient for the sample gene (Keep only one isoform per gene)
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.icgc.enrich.1<-function(icgc.mut, brca.exp, table.2){
  #Finds enriched driver metabolites based on mutation and expression analysis
  
  #Filter expression and mutation info against each other for common samples
  common.samples<-intersect(unique(icgc.mut$SAMPLE), unique(colnames(brca.exp$tumor)))
  icgc.mut<-icgc.mut[SAMPLE %in% common.samples,]
  brca.exp$tumor<-brca.exp$tumor[,common.samples]
  
  #Filter expression matrix by non-zero expression
  brca.exp$tumor<-brca.exp$tumor[rowSums(brca.exp$tumor)!=0,]
  
  #Assign sample mutations to kegg identifiers based on associated genes
  icgc.mut<-icgc.mut[MUTATION=="MISSENSE",]
  table.2<-unique(table.2[,c("KEGG_ID", "GENE"),with=F])
  setnames(table.2, c("KEGG.ID", "Hugo_Symbol"))
  mut.table<-merge(icgc.mut, table.2, by="Hugo_Symbol")
  
  #Filter mutation table by those metabolites that have at least 10 samples assocaited to them
  mut.table[,N.SAMPLES:=length(unique(SAMPLE)), by="KEGG.ID"]
  mut.table<-mut.table[N.SAMPLES>=10,]
  
  #Filter expression metabolites by those that have at least 4 genes associated to them
  exp.table<-table.2[Hugo_Symbol %in% rownames(brca.exp$tumor),]
  exp.table[,N.GENES:=length(unique(Hugo_Symbol)), by="KEGG.ID"]
  exp.table<-exp.table[N.GENES>=4,]
  exp.mets<-unique(exp.table$KEGG.ID)
  
  #Iteration through filtered mutation mets and obtain significantly differentiated values for expression mets
  internal.function<-function(samples, KEGG.ID){
    print (KEGG.ID)
    
    not.samples<-setdiff(colnames(brca.exp$tumor), samples)
    met.pvals<-sapply(exp.mets, function(x) {
      met.hugos<-unique(exp.table[KEGG.ID==x, ]$Hugo_Symbol)
      PVAL=wilcox.test(rowMeans(brca.exp$tumor[met.hugos, samples]), rowMeans(brca.exp$tumor[met.hugos, not.samples]), paired = T)$p.value
      return(PVAL)
    })
    main.mets<-data.table(EXP.METS=exp.mets, PVAL=met.pvals)
    main.mets$PVAL.ADJ<-p.adjust(main.mets$PVAL, method="fdr")
    
    return(list(EXP.METS=main.mets$EXP.METS, PVAL=main.mets$PVAL, PVAL.ADJ=main.mets$PVAL.ADJ))
  }
  
  print (length(unique(mut.table$KEGG.ID)))
  main.table<-mut.table[,internal.function(unique(SAMPLE), unique(KEGG.ID)),by="KEGG.ID"]
  
  #Return
  return(main.table)
}

Function.kegg.path.enrich<-function(kegg.path, sig.mets, original.kegg, table.2=F, th=40){
  #Calculates enrichment for kegg pathway. Filters out high degree met from original table
  
  #Filter sig.mets
  if (table.2==T){
    kegg.filt<-unique(original.kegg[, c("KEGG_ID", "GENE"), with=F])
    kegg.filt<-kegg.filt[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<th,]$KEGG_ID
  } else {
    kegg.filt<-unique(rbind(data.table(KEGG_ID=original.kegg$SUBSTRATE, GENE=original.kegg$Hugo_Symbol),
                            data.table(KEGG_ID=original.kegg$PRODUCT, GENE=original.kegg$Hugo_Symbol)))
    kegg.filt<-kegg.filt[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<th,]$KEGG_ID
  }
  sig.mets<-sig.mets[sig.mets %in% kegg.filt]
  
  #Remove "- Homo sapiens" from path description
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - Homo"))[1])
  
  #Remove "Metabolic pathways" from kegg.path
  kegg.path<-kegg.path[!grepl("metabolic", DESCRIPTION, ignore.case = T),]
  
  #Remove non-direct pathways of other disease
  kegg.path<-kegg.path[!DESCRIPTION %in% c("Asthma","Amoebiasis"),]
  
  #Filter kegg.path for those pathways present only in significant mets
  sig.mets<-sig.mets[sig.mets %in% unique(kegg.path$COMPOUND)]
  filt.paths<-unique(kegg.path[COMPOUND %in% sig.mets, ]$DESCRIPTION)
  kegg.path<-kegg.path[DESCRIPTION %in% filt.paths,]
  
  #Calculate enrichment
  all.path.mets<-unique(kegg.path$COMPOUND)
  main.table<-kegg.path[,list(PVAL=phyper(q = length(intersect(sig.mets, unique(COMPOUND)))-1,
                                          m = length(unique(COMPOUND)),
                                          n = length(all.path.mets) - length(unique(COMPOUND)),
                                          k = length(sig.mets), lower.tail=F)), by="DESCRIPTION"]
  
  #Correct for fdr
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method = "fdr")
  
  #Clean up and return
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.met.node.distance<-function(original.kegg, table.2=F, th=40, degree.th=40){
  #Gene treshold should be same as for path enrichment
  
  #Filter table metabolites by high associations
  if (table.2==T){
    kegg.filt<-unique(original.kegg[, c("KEGG_ID", "GENE"), with=F])
    kegg.filt<-kegg.filt[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<th,]$KEGG_ID
    
    original.kegg<-original.kegg[KEGG_ID %in% kegg.filt,]
    kegg.1<-unique(original.kegg[,c("KEGG_ID", "GENE"),with=F])
    kegg.2<-unique(original.kegg[,c("KEGG_ID", "GENE"),with=F])
    setnames(kegg.1, c("SUBSTRATE", "GENE"))
    setnames(kegg.2, c("PRODUCT", "GENE"))
    original.kegg<-merge(kegg.1, kegg.2, by="GENE", allow.cartesian=T)
    
  } else {
    kegg.filt<-unique(rbind(data.table(KEGG_ID=original.kegg$SUBSTRATE, GENE=original.kegg$Hugo_Symbol),
                            data.table(KEGG_ID=original.kegg$PRODUCT, GENE=original.kegg$Hugo_Symbol)))
    kegg.filt<-kegg.filt[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<th,]$KEGG_ID
    original.kegg<-original.kett[SUBSTRATE %in% kegg.filt & PRODUCT %in% kegg.filt,]
  }
  
  #Build graph and obtain distance - Threshold by degree
  original.kegg<-original.kegg[SUBSTRATE!=PRODUCT,]
  kegg.graph<-graph.data.frame(unique(original.kegg[,c("SUBSTRATE", "PRODUCT"), with=F]), directed = F)
  kegg.degree<-data.table(as.matrix(degree(kegg.graph, mode = "all",normalized = F)), keep.rownames = T)
  kegg.degree.filt<-unique(kegg.degree[V1<=degree.th,]$rn)
  kegg.graph<-graph.data.frame(unique(original.kegg[SUBSTRATE %in% kegg.degree.filt & PRODUCT %in% kegg.degree.filt,][,c("SUBSTRATE", "PRODUCT"), with=F]), directed = F)
  kegg.distance<-shortest.paths(kegg.graph)
  
  #Return distance
  return(list(KEGG.DISTANCE=kegg.distance, KEGG.GRAPH=kegg.graph))
}

Function.met.assign.distance<-function(icgc.enrich.met, dist.matrix, pval.th=0.05, pval.column="PVAL.ADJ"){
  #MAKE SURE DISTANCE MATRIX WAS BUILT USING THE SAME METABOLITE DATA TABLE USED TO CONSTRUCT ENRICHMENT TABLE!!
  
  #Threshold enriched table
  icgc.met.enrich<-icgc.enrich.met[icgc.enrich.met[[pval.column]]<pval.th,]
  
  #Filter by presence in distance matrix (Presumably correct since we filtered both before)
  icgc.met.enrich<-icgc.met.enrich[KEGG.ID %in% rownames(dist.matrix) & EXP.METS %in% rownames(dist.matrix)]
  
  #Find distances
  icgc.met.enrich$DIST<-apply(icgc.met.enrich, 1, function(x) dist.matrix[x[1], x[2]])
  
  #Clean up
  icgc.met.enrich<-icgc.met.enrich[DIST!=Inf,]
  
  #Add EXP.METS count for each KEGG.ID MUT causal
  icgc.met.enrich[,EXP.MET.COUNT:=length(unique(EXP.METS)),by="KEGG.ID"] 
  
  #Clean up and return
  icgc.met.enrich$KEGG.ID<-as.character(icgc.met.enrich$KEGG.ID)
  icgc.met.enrich$EXP.METS<-as.character(icgc.met.enrich$EXP.METS)
  return(icgc.met.enrich)
}

Function.enrich.distance<-function(assigned.distance, dist.matrix, exp.met.count.th=5){
  #MAKE SURE DISTANCE MATRIX WAS BUILT USING THE SAME METABOLITE DATA TABLE USED TO CONSTRUCT ENRICHMENT TABLE!!
  #Performs test to find if enriched metabolites (EXP.METs) are closer to KEGG.ID than rest of metaboiltes in network
  
  #Filter for metabolites that have and exp.met count great than threshold 
  assigned.distance<-assigned.distance[EXP.MET.COUNT>=exp.met.count.th,]
  
  #Calculate statisitcs
  main.keggs<-unique(assigned.distance$KEGG.ID)
  main.pvals<-sapply(main.keggs, function(x) {
    self.keggs<-assigned.distance[KEGG.ID==x,]$EXP.METS
    other.keggs<-setdiff(colnames(dist.matrix), c(self.keggs, x))
    self.dist<-assigned.distance[EXP.METS %in% self.keggs,]$DIST
    other.dist<-dist.matrix[x, other.keggs]
    other.dist<-other.dist[other.dist!=Inf]
    
    pval<-t.test(self.dist, other.dist, alternative="less")$p.value
    return(pval)
  })
  
  #Clean up and return
  main.table<-data.table(KEGG.ID=main.keggs, PVAL=main.pvals)
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.icgc.enrich.coverage<-function(enrich.table, icgc.mut, original.kegg, pval.th=0.05, kegg.th=5, pval.col="PVAL.ADJ", gene.th=40, met.hub.th=0.6, table.2=T){
  #Obtain coverage per significantly influencial metabolite and metabolic hubs depending on influenced metabolites jaccard threshold
  #"original.kegg" has to be the same dataset used to build enrich.table
  
  #PRODUCES:
  #   Cast matrix of kegg.id jaccard for all influenced exp.mets
  #   Table of covered samples per kegg.id
  #   Graph of metabolic hubs
  #   Table of covered samples per metabolic kegg.id hub (based on met.hub.th)
  #   Cast matrix of kegg.id jaccard for all associated genes (per hub)
  #       - This serves to establish the point that kegg.ids that are in metabolic hubs (influencing same metabolites) do no do so because of same mutated genes
  
  #Threshold by p-value
  enrich.table<-enrich.table[enrich.table[[pval.col]]<0.05,]
  
  #Filter by associated genes by original table
  if (table.2==T){
    original.kegg<-unique(original.kegg[, c("KEGG_ID", "GENE"), with=F])
    kegg.filt<-original.kegg[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<gene.th,]$KEGG_ID
    
  } else {
    original.kegg<-unique(rbind(data.table(KEGG_ID=original.kegg$SUBSTRATE, GENE=original.kegg$Hugo_Symbol),
                                data.table(KEGG_ID=original.kegg$PRODUCT, GENE=original.kegg$Hugo_Symbol)))
    kegg.filt<-original.kegg[,list(N=length(GENE)),by="KEGG_ID"]
    kegg.filt<-kegg.filt[N<gene.th,]$KEGG_ID
  }
  enrich.table<-enrich.table[KEGG.ID %in% kegg.filt & EXP.METS %in% kegg.filt,]
  
  #Filter for at least exp.mets per kegg.id
  enrich.table[,N:=length(unique(EXP.METS)), by="KEGG.ID"]
  enrich.table<-enrich.table[N>=kegg.th,]
  
  #Obtain table of covered samples per kegg id
  all.mets<-unique(enrich.table$KEGG.ID)
  mut.cover<-original.kegg[KEGG_ID %in% all.mets,][,c("KEGG_ID", "GENE"),with=F]
  setkey(mut.cover)
  mut.cover<-unique(mut.cover)
  setnames(mut.cover, c("KEGG.ID", "Hugo_Symbol"))
  
  icgc.mut<-icgc.mut[MUTATION=="MISSENSE",]
  mut.cover<-merge(mut.cover, icgc.mut, by="Hugo_Symbol")
  mut.coverage<-mut.cover[,list(N.SAMPLES=length(unique(SAMPLE))), by="KEGG.ID"] ### 2
  
  #Obtain kegg affected coverage but causal kegg
  kegg.cov<-matrix(ncol=length(all.mets), nrow=length(all.mets), dimnames = list(all.mets, all.mets)) ### 1
  for (i in all.mets){
    for (j in all.mets){
      kegg.cov[i,j]<-length(intersect(enrich.table[KEGG.ID==i,]$EXP.METS, enrich.table[KEGG.ID==j,]$EXP.METS))/
        length(union(enrich.table[KEGG.ID==i,]$EXP.METS, enrich.table[KEGG.ID==j,]$EXP.METS))
    }
  }
  
  #Obtain metabolic hubs
  met.hubs<-data.table(melt(kegg.cov))[Var1!=Var2,][value>=met.hub.th,]
  met.hubs<-graph.data.frame(met.hubs) ### 3
  
  hub.list<-lapply(decompose.graph(met.hubs), function(x) V(x)$name)
  hub.coverage<-data.table(MET.HUB=sapply(hub.list, function(x) paste(x, collapse = ".")),
                           HUB.COVERAGE=sapply(hub.list, function(x)  length(unique(mut.cover[KEGG.ID %in% x,]$SAMPLE)) )) ### 4
  
  #Obtain gene jaccard coverage per hub
  hub.gene.cover<-lapply(hub.list, function(y) {
    
    hub.met.gene.jaccard<-matrix(nrow=length(y), ncol=length(y), dimnames = list(y,y))
    for (i in y){
      for (j in y){
        hub.met.gene.jaccard[i,j]<-length(intersect(original.kegg[KEGG_ID==i,]$GENE, original.kegg[KEGG_ID==j,]$GENE))/
          length(union(original.kegg[KEGG_ID==i,]$GENE, original.kegg[KEGG_ID==j,]$GENE))
      }
    }
    return(hub.met.gene.jaccard)
  })
  names(hub.gene.cover)<-sapply(hub.list, function(x) paste(x, collapse = ".")) ### 5
  
  #Return
  return(list(KEGG.COV.MET=kegg.cov,  KEGG.COV.SAMPLES=mut.coverage, MET.HUBS.GRAPH=met.hubs, HUB.COV.SAMPLES=hub.coverage, HUB.COV.GENES=hub.gene.cover))
}