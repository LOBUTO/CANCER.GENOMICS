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

