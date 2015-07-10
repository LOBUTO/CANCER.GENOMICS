Function.teru.diff.gene.exp<-function(teru.obj){
  
  #Clean matrix
  gene.exp<-teru.obj$MATRIX[complete.cases(teru.obj$MATRIX),]
  
  #Find differentially expressed genes
  cancer<-unique(teru.obj$CLASS[CLASS=="Tumor",]$SAMPLE)
  normal<-unique(teru.obj$CLASS[CLASS=="Normal",]$SAMPLE)
  
  print ("Calculating p-values")
  pvals<-apply(gene.exp, 1, function(x) wilcox.test(x[cancer], x[normal], paired=F)$p.value)
  lgfc<-apply(gene.exp, 1, function(x) log2(median(x[cancer])/median(x[normal])))
  
  #Correct for fdr
  main.table<-data.table(Hugo_Symbol=rownames(gene.exp), PVAL=pvals, LG.FC=lgfc)
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
  
  #Clean up and return
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.recon.edges<-function(recon.table) {
  #Produces substrate-hugo-product table for recon dataset
  
  #Remove Acyl carrier protein (acp) - keep in mind when doing enrichment!!!
  recon.product<-recon.table$PRODUCT[!grepl("acp", NAME, ignore.case = T), ]
  recon.substrate<-recon.table$SUBSTRATE[!grepl("acp", NAME, ignore.case = T), ]
  
  #Pre-clean up
  recon.product<-recon.product[!grepl("_", Hugo_Symbol),]
  recon.product<-recon.product[KEGG_ID!="NONE",]
  recon.product<-recon.product[,c("REACTION.ID", "Hugo_Symbol", "KEGG_ID"), with=F]
  setnames(recon.product, c("REACTION.ID", "Hugo_Symbol", "PRODUCT"))
  
  recon.substrate<-recon.substrate[!grepl("_", Hugo_Symbol),]
  recon.substrate<-recon.substrate[KEGG_ID!="NONE",]
  recon.substrate<-recon.substrate[,c("REACTION.ID", "Hugo_Symbol", "KEGG_ID"), with=F]
  setnames(recon.substrate, c("REACTION.ID", "Hugo_Symbol", "SUBSTRATE"))
  
  #Manually insert reaction information (for R-2HG)
  recon.product<-rbind(recon.product, data.table(REACTION.ID=c(rep("ONCO.REACT.1",2), rep("ONCO.REACT.2", 2), rep("ONCO.REACT.3",1)),
                                                 Hugo_Symbol=c("ADHFE1", "ADHFE1", 
                                                               "PHGDH",  "PHGDH",
                                                               "D2HGDH"),
                                                 PRODUCT=c("C00164", "C01087",
                                                           "C03232",  "C00026",
                                                           "C00026")))
  
  recon.substrate<-rbind(recon.substrate, data.table(REACTION.ID=c(rep("ONCO.REACT.1",2), rep("ONCO.REACT.2", 2), rep("ONCO.REACT.3",1)),
                                                     Hugo_Symbol=c("ADHFE1", "ADHFE1", 
                                                                   "PHGDH",  "PHGDH",
                                                                   "D2HGDH"),
                                                     SUBSTRATE=c("C03197", "C00026",
                                                               "C01087", "C00197",
                                                               "C01087")))
  setkey(recon.product)
  setkey(recon.substrate)
  recon.substrate<-unique(recon.substrate)
  recon.product<-unique(recon.product)
  
  #Obtain edge table
  main.table<-merge(recon.substrate, recon.product, by=c("REACTION.ID", "Hugo_Symbol"), allow.cartesian = T)
  main.table$REACTION.ID<-NULL
  
  #Clean up and return
  main.table<-main.table[SUBSTRATE!=PRODUCT,]
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Funcion.directed.network<-function(recon.table) {
  
  #Remove Acyl carrier protein (acp) - keep in mind when doing enrichment!!!
  recon.product<-recon.table$PRODUCT[!grepl("acp", NAME, ignore.case = T), ]
  recon.substrate<-recon.table$SUBSTRATE[!grepl("acp", NAME, ignore.case = T), ]
  
  #Pre-clean up
  recon.product<-recon.product[!grepl("_", Hugo_Symbol),]
  recon.product<-recon.product[KEGG_ID!="NONE",]
  recon.product<-recon.product[,c("REACTION.ID", "Hugo_Symbol", "KEGG_ID"), with=F]
  setnames(recon.product, c("REACTION.ID", "Hugo_Symbol", "PRODUCT"))
  
  recon.substrate<-recon.substrate[!grepl("_", Hugo_Symbol),]
  recon.substrate<-recon.substrate[KEGG_ID!="NONE",]
  recon.substrate<-recon.substrate[,c("REACTION.ID", "Hugo_Symbol", "KEGG_ID"), with=F]
  setnames(recon.substrate, c("REACTION.ID", "Hugo_Symbol", "SUBSTRATE"))
  
  #Manually insert reaction information (for R-2HG)
  recon.product<-rbind(recon.product, data.table(REACTION.ID=c(rep("ONCO.REACT.1",2), rep("ONCO.REACT.2", 2), rep("ONCO.REACT.3",1)),
                                                 Hugo_Symbol=c("ADHFE1", "ADHFE1", 
                                                               "PHGDH",  "PHGDH",
                                                               "D2HGDH"),
                                                 PRODUCT=c("C00164", "C01087",
                                                           "C03232",  "C00026",
                                                           "C00026")))
  
  recon.substrate<-rbind(recon.substrate, data.table(REACTION.ID=c(rep("ONCO.REACT.1",2), rep("ONCO.REACT.2", 2), rep("ONCO.REACT.3",1)),
                                                     Hugo_Symbol=c("ADHFE1", "ADHFE1", 
                                                                   "PHGDH", "PHGDH",
                                                                   "D2HGDH"),
                                                     SUBSTRATE=c("C03197", "C00026",
                                                                 "C01087", "C00197",
                                                                 "C01087")))
  setkey(recon.product)
  setkey(recon.substrate)
  recon.substrate<-unique(recon.substrate)
  recon.product<-unique(recon.product)
  
  #Obtain edge table
  main.table<-merge(recon.substrate, recon.product, by=c("REACTION.ID", "Hugo_Symbol"), allow.cartesian = T)
  main.table<-main.table[SUBSTRATE!=PRODUCT,]
  main.table$REACTION.ID<-NULL
  setkey(main.table)
  main.table<-unique(main.table)
  
  #Obtain Jaccard per two metabolites
  main.copy<-copy(main.table)
  setnames(main.copy, c("gene", "sub", "prod"))
  main.table<-main.table[,list(JACCARD= length(unique(Hugo_Symbol))/length(unique(c(main.copy[sub==SUBSTRATE,]$gene, main.copy[prod==PRODUCT,]$gene)))), 
                         by=c("SUBSTRATE", "PRODUCT")]
  
  #Clean up and return
  setkey(main.table)
  main.table<-unique(main.table)
  main.table<-main.table[order(JACCARD, decreasing = T),]
  return(main.table)
}

Function.kegg.directed<-function(kegg.file, recon.met.file, weight.filter=1000) {
  
  #Load file
  kegg.table<-fread(kegg.file, header=T, sep="\t", stringsAsFactors = F)
  setkey(kegg.table)
  kegg.table<-unique(kegg.table)
  
  #Manually insert reaction information (for R-2HG)
  kegg.table<-rbind(kegg.table, data.table(Hugo_Symbol=c("ADHFE1", "ADHFE1","ADHFE1", "ADHFE1",
                                                         "PHGDH", "PHGDH", "PHGDH", "PHGDH",
                                                         "D2HGDH"),
                                           SUBSTRATE=c("C03197","C03197", "C00026", "C00026",
                                                       "C01087","C00197","C01087", "C00197",
                                                       "C01087"),
                                           PRODUCT=c("C00164", "C01087", "C00164", "C01087",
                                                     "C03232","C03232","C00026","C00026",
                                                     "C00026")))
  
  #Pre-clean up
  kegg.table<-kegg.table[SUBSTRATE!=PRODUCT,]
  setkey(kegg.table)
  kegg.table<-unique(kegg.table)
  
  #Form edge network
  kegg.copy<-copy(kegg.table)
  setnames(kegg.copy, c("gene", "sub", "prod"))
  edge.net<-kegg.table[,list(JACCARD=length(unique(Hugo_Symbol)) / length(union(kegg.copy[sub==SUBSTRATE,]$gene, kegg.copy[prod==PRODUCT,]$gene))), 
                       by=c("SUBSTRATE", "PRODUCT")]
  
  #Filter by weight using recon.table
  recon.met<-fread(recon.met.file, header=T, sep="\t", stringsAsFactors = F)
  recon.met<-recon.met[KEGG_ID!="NONE",]
  recon.met<-recon.met[WEIGHT!="NONE",]
  recon.met$WEIGHT<-as.numeric(recon.met$WEIGHT)
  recon.met<-recon.met[WEIGHT<weight.filter,]
  recon.met<-unique(recon.met$KEGG_ID)
  edge.net<-edge.net[!(SUBSTRATE %in% recon.met),][!(PRODUCT %in% recon.met),]
  edge.net<-edge.net[SUBSTRATE!="C00087",][PRODUCT!="C00087",] #Manual substrate
  
  #Clean up and return
  edge.net<-edge.net[order(JACCARD,decreasing = T),]
  return(edge.net)
}

Function.kegg.edges<-function(kegg.file, recon.met.file, weight.filter=1000){
  
  #Load file
  kegg.table<-fread(kegg.file, header=T, sep="\t", stringsAsFactors = F)
  setkey(kegg.table)
  kegg.table<-unique(kegg.table)
  
  #Manually insert reaction information (for R-2HG)
  kegg.table<-rbind(kegg.table, data.table(Hugo_Symbol=c("ADHFE1", "ADHFE1","ADHFE1", "ADHFE1",
                                                         "PHGDH", "PHGDH", "PHGDH", "PHGDH",
                                                         "D2HGDH"),
                                           SUBSTRATE=c("C03197","C03197", "C00026", "C00026",
                                                       "C01087","C00197","C01087", "C00197",
                                                       "C01087"),
                                           PRODUCT=c("C00164", "C01087", "C00164", "C01087",
                                                     "C03232","C03232","C00026","C00026",
                                                     "C00026")))
  
  #Pre-clean up
  kegg.table<-kegg.table[SUBSTRATE!=PRODUCT,]
  setkey(kegg.table)
  kegg.table<-unique(kegg.table)
  
  #Filter by weight using recon.table
  recon.met<-fread(recon.met.file, header=T, sep="\t", stringsAsFactors = F)
  recon.met<-recon.met[KEGG_ID!="NONE",]
  recon.met<-recon.met[WEIGHT!="NONE",]
  recon.met$WEIGHT<-as.numeric(recon.met$WEIGHT)
  recon.met<-recon.met[WEIGHT<weight.filter,]
  recon.met<-unique(recon.met$KEGG_ID)
  kegg.table<-kegg.table[!(SUBSTRATE %in% recon.met),][!(PRODUCT %in% recon.met),]
  kegg.table<-kegg.table[SUBSTRATE!="C00087",][PRODUCT!="C00087",] #Manual substrate
  
  #Clean up and return
  return(kegg.table)
}

Function.kegg.path.hm<-function(kegg.path, target.mets){
  
  #Remove "- Homo sapiens" from path description
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - Homo"))[1])
  
  #Remove "Metabolic pathways" from kegg.path
  kegg.path<-kegg.path[!grepl("metabolic", DESCRIPTION, ignore.case = T),]
  
  #Filter kegg path for mets of interest
  kegg.path<-kegg.path[COMPOUND %in% target.mets,]
  
  #Cast into matrix ready for plotting
  kegg.path$PATHWAY<-NULL
  kegg.path$N<-1
  setkey(kegg.path)
  kegg.path<-unique(kegg.path)
  kegg.count<-kegg.path[,list(COUNT=length(unique(COMPOUND))), by="DESCRIPTION"] #Filter low counts
  kegg.path<-kegg.path[DESCRIPTION %in% kegg.count[COUNT>=2,]$DESCRIPTION,]
  main.cast<-acast(kegg.path, COMPOUND~DESCRIPTION, fill = 0, value.var = "N")
  
  #Return
  return(main.cast)
}

Function.kegg.path.enrich<-function(kegg.path, prr.pval.results){
  
  #Remove "- Homo sapiens" from path description
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - Homo"))[1])
  
  #Remove "Metabolic pathways" from kegg.path
  kegg.path<-kegg.path[!grepl("metabolic", DESCRIPTION, ignore.case = T),]
  
  #Filter kegg path for those metabolites that we actually tested
  kegg.path<-kegg.path[COMPOUND %in% unique(prr.pval.results$SUBSTRATE),]
  
  #Filter kegg.path for those pathways present only in significant mets
  sig.mets<-prr.pval.results[PVAL.ADJ<0.05,]$SUBSTRATE
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

Function.network.graph<-function(network.table, d.filter=1000){
  
  require(igraph)
  
  #Pre-clean network from transposable mets
  network.table<-network.table[SUBSTRATE!=PRODUCT,]
  
  #Filter table by degree
  filt.table<-rbind(data.table(ID.1=network.table$SUBSTRATE, ID.2=network.table$PRODUCT),
                    data.table(ID.1=network.table$PRODUCT, ID.2=network.table$SUBSTRATE))
  filt.table<-unique(filt.table)
  filt.table<-filt.table[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
  filt.table<-unique(filt.table[DEGREE<d.filter,]$ID.1)
  network.table<-network.table[SUBSTRATE %in% filt.table,][PRODUCT %in% filt.table,]
  
  #Create directed graph
  graph.table<-graph.data.frame(network.table, directed = T)
  
  #Return
  return(graph.table)
}

Function.met.scoring<-function(network.table, met.targets, d=0.85, d.filter=1000) {
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-clean network.table from transposable mets
  network.table<-unique(network.table[,c("SUBSTRATE", "PRODUCT"), with=F])
  network.table<-network.table[SUBSTRATE!=PRODUCT,]
  
  #Filter network table by degree (filter out currency metabolites)
  degree.filter<-unique(rbind(data.table(ID.1=network.table$SUBSTRATE, ID.2=network.table$PRODUCT),
                              data.table(ID.1=network.table$PRODUCT, ID.2=network.table$SUBSTRATE)))
  degree.filter<-degree.filter[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
  degree.filter<-degree.filter[DEGREE<d.filter,]
  print (degree.filter[ID.1 %in% met.targets,])
  print (degree.filter[DEGREE %in% degree.filter[ID.1 %in% met.targets,]$DEGREE,])
  print (table(degree.filter[DEGREE %in% degree.filter[ID.1 %in% met.targets,]$DEGREE,]$DEGREE))
  print (table(degree.filter$DEGREE))
  
  network.table<-network.table[SUBSTRATE %in% unique(degree.filter$ID.1),][PRODUCT %in% unique(degree.filter$ID.1),]
  print (dim(network.table))
  print (met.targets %in% network.table$SUBSTRATE)
  print (met.targets %in% network.table$PRODUCT)
  
  #Create transition matrix
  t.m<-network.table[,list(DEGREE=length(unique(PRODUCT))), by="SUBSTRATE"]
  t.m<-merge(t.m, network.table[,c("SUBSTRATE", "PRODUCT"), with=F], by = "SUBSTRATE")
  t.m$trans<-1/t.m$DEGREE #may modify depending on weighted jaccard edges
  degree.info<-unique(data.table(SUBSTRATE=t.m$SUBSTRATE, DEGREE=t.m$DEGREE))
  degree.info<-rbind(degree.info, data.table(SUBSTRATE=unique(t.m$PRODUCT)[!(unique(t.m$PRODUCT) %in% degree.info$SUBSTRATE)], DEGREE=0))
  t.m$DEGREE<-NULL
  t.m<-rbind(t.m, data.table(SUBSTRATE=unique(t.m$SUBSTRATE), PRODUCT=unique(t.m$SUBSTRATE), trans=0))
  
  all.mets<-unique(c(t.m$SUBSTRATE, t.m$PRODUCT))
  adj.t.m<-matrix(ncol=length(all.mets), nrow=length(all.mets), dimnames = list(all.mets, all.mets))
  adj.t.m[is.na(adj.t.m)]<-0
  apply(as.matrix(t.m), 1, function(x) {
    adj.t.m[x[2], x[1]]<<-as.numeric(x[3])
  })
  
  #Obtain personalized vector for path
  #personalized.vector<-unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  personalized.vector<-met.targets[met.targets %in% rownames(adj.t.m)]
  personalized.vector<-ifelse(rownames(adj.t.m) %in% personalized.vector, 1/length(personalized.vector), 0)
  
  #####Iterate to minimize mean change per row and converge to stable network####
  initial.rank<-rep(1/nrow(adj.t.m), nrow(adj.t.m))
  
  #Obtain first iteration page rank
  A.1<-(adj.t.m*d)%*%initial.rank + (1-d)*personalized.vector #Implementing personalized vector
  
  #Iterate to converge on optimal page rank
  delta<-Inf
  while(delta>0.0000000001){
    delta<-mean(abs(A.1 - ((adj.t.m*d)%*%A.1 + (1-d)*personalized.vector)))
    A.1<-((adj.t.m*d)%*%A.1 + (1-d)*personalized.vector)
  }
  
  #####Obtained the non-personalized version####
  personalized.vector<-rep(1/nrow(adj.t.m), nrow(adj.t.m)) #all equal for all
  
  #Obtain first iteration page rank
  A.2<-(adj.t.m*d)%*%initial.rank + (1-d)*personalized.vector #Implementing personalized vector
  
  #Iterate to converge on optimal page rank
  delta<-Inf
  while(delta>0.0000000001){
    delta<-mean(abs(A.2 - ((adj.t.m*d)%*%A.2 + (1-d)*personalized.vector)))
    A.2<-((adj.t.m*d)%*%A.2 + (1-d)*personalized.vector)
  }
  ##############################################
  
  #Normalize for unpersonalized rank (degree normalization)
  A.1<-as.vector(A.1)
  A.2<-as.vector(A.2)
  PAGE.RANK<-A.1
  
  #Clean up and return
  main.table<-data.table(SUBSTRATE=rownames(adj.t.m), PAGE.RANK, A.1, A.2)
  main.table$PAGE.RANK<-normalize.vector(main.table$PAGE.RANK)
  main.table$PAGE.RANK<-main.table$PAGE.RANK / sum(main.table$PAGE.RANK) #Normalizing page rank to 1
  main.table<-merge(main.table, degree.info, by="SUBSTRATE")
  main.table$TARGET.PATH<-main.table$SUBSTRATE %in% met.targets
  #main.table$TARGET.PATH<-main.table$SUBSTRATE %in% unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  main.table<-main.table[order(PAGE.RANK, decreasing = T),]
  return(main.table)
}

Function.met.scoring.pvalue<-function(prr, network.table, met.targets, d=0.9, d.filter=100){
  #NOTE: d, d.filter, met.targers have to specific to prr
  
  require(parallel)
  
  #Clean network
  network.table<-network.table[SUBSTRATE!=PRODUCT,]
  
  #Obtain ground ranks from the prr
  prr<-prr[order(PAGE.RANK, decreasing=T),]
  prr.ranks<-data.table(SUBSTRATE=prr$SUBSTRATE, RANK=1:nrow(prr))
  
  #Filter network by degree
  degree.filter<-unique(rbind(data.table(ID.1=network.table$SUBSTRATE, ID.2=network.table$PRODUCT),
                              data.table(ID.1=network.table$PRODUCT, ID.2=network.table$SUBSTRATE)))
  degree.filter<-degree.filter[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
  degree.filter<-degree.filter[DEGREE<d.filter,]
  network.table<-network.table[SUBSTRATE %in% unique(degree.filter$ID.1),][PRODUCT %in% unique(degree.filter$ID.1),]
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "Function.met.scoring", "d", "d.filter") ,envir=environment())
  
  #Obtain shuffled candidate mets of the same size as target that meet degree sum condition
  sub.mets<-setdiff(unique(network.table$SUBSTRATE), met.targets)
  degree.sum<-sum(unique(degree.filter[ID.1 %in% met.targets,]$DEGREE))
  shuffled<-replicate(25000, sample(sub.mets, length(met.targets)), simplify = F)
  shuffled.deg.sum<-parSapply(cl, shuffled, function(x) sum(degree.filter[ID.1 %in% x,]$DEGREE))
  shuffled<-shuffled[shuffled.deg.sum==degree.sum]
  shuffled<-lapply(shuffled, function(x) sort(x)) #Fix for duplicates
  print (length(shuffled))
  print (nrow(unique(data.table(do.call(rbind, shuffled)))))
  
  #Calculate p-values
  c<-1
  main.list<-parLapply(cl, shuffled, function(x) {
    
    print (c)
    temp.prr<-Function.met.scoring(network.table, x, d, d.filter)
    
    #Obtain ranks for shuffle
    shuffled.ranks<-data.table(SHUFFLED.SUBSTRATE=temp.prr$SUBSTRATE, SHUFFLED.RANK=1:nrow(temp.prr))
    c<<-c+1
    return(shuffled.ranks)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Obtain p-values (By looking if we can obtain lower than or equal ranks in shuffled networks), correct for fdr
  SHUFFLED.RANKS<-do.call(rbind, main.list)
  main.table<-prr.ranks[,list(PVAL= mean(SHUFFLED.RANKS[SHUFFLED.SUBSTRATE==SUBSTRATE,]$SHUFFLED.RANK<=RANK)), by="SUBSTRATE"]
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method = "fdr")
  
  #Clean up and return
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.BEST.PRR.PARAM<-function(met.diff, network.table, target.mets){
  
  #Test parameters
  PRR.DIR<-data.table()
  for (d.t in c(50, 70, 80, 100,110,120, 150, 200, 300, 900)){
    for (c in seq(0.1,0.95,0.05)){
      for (j in c(T,F)){
        cat(d.t, c, j)
        temp.table<-Function.PRR.PDISCOVER(met.diff, network.table, d=c, jaccard=j, d.filter=d.t)
        temp.table<-Function.KEGG.ROC(temp.table, target.mets)
        temp.table$JACCARD<-j
        temp.table$DEGREE.TH<-d.t
        temp.table$d<-c
        PRR.DIR<-rbind(PRR.DIR, temp.table)
      }
    }
  }
  
  #Clean up and return
  PRR.DIR<-PRR.DIR[order(AUC, decreasing = T),]
  return(PRR.DIR)
}

Function.PRR.REGRESSION<-function(diff.exp, edge.table, met.diff, network.table, target.mets, d=0.70, jaccard = F, d.filter = 100){
  
  #Clean network table
  network.table<-unique(network.table[,list(JACCARD=mean(JACCARD)), by=c("SUBSTRATE", "PRODUCT")])
  network.table<-network.table[SUBSTRATE!=PRODUCT,]
  
  #Obtain a score per edge based on genes that are significantly dysregulated (for those we have information on)
  common.genes<-intersect(unique(diff.exp$Hugo_Symbol), unique(edge.table$Hugo_Symbol))
  diff.exp<-diff.exp[Hugo_Symbol %in% common.genes,] 
  edge.table<-edge.table[Hugo_Symbol %in% common.genes,]
  
  #Score is average of absolute value log fold change for significantly dysregulated genes
  diff.exp$LG.FC<-ifelse(diff.exp$PVAL.ADJ<0.05, diff.exp$LG.FC, 0)
  setnames(diff.exp, c("gene", "PVAL", "LG.FC", "PVAL.ADJ"))
  edge.score<-edge.table[, list(SCORE=mean(  abs(diff.exp[gene %in% Hugo_Symbol,]$LG.FC) )), by=c("SUBSTRATE", "PRODUCT")]
  edge.score<-edge.score[order(SCORE, decreasing = T),]
  ordered.scores<-unique(edge.score$SCORE)
  print (ordered.scores)
  
  #Calculate initial conditions of PRR Reverse
  initial.cond<-Function.PRR.PDISCOVER(met.diff, network.table, d=d, jaccard = jaccard, d.filter = d.filter)
  current.ranks<-which(initial.cond$KEGG.ID %in% target.mets)
  
  #Delete network edges(network.table) progressively in order of decreasing gene expression dysregulation score (edge.score)
  network.matrix<-as.matrix(network.table)
  print (dim(network.matrix))
  print (current.ranks)
  
  delta<-Inf #delta will be measure at every point until there is no more gain in ranks
  del.paths<-list()
  pro.paths<-list()
  for (s in ordered.scores){
    
    network.matrix<-as.matrix(network.matrix)
    
    #Edges to be modified
    s.edges<-edge.score[SCORE==s,]
    #print (s.edges)
    
    #Convert to matrix
    s.edges<-as.matrix(s.edges[, c("SUBSTRATE", "PRODUCT"), with=F]) #same order as in network table
    
    #Create dummy matrix to test effect of deletion
    dummy.matrix<-copy(network.matrix)
    
    #Delete edges from network of present score
    for (n in 1:nrow(s.edges)){
      del.edge<-as.vector(s.edges[n,])
      del.rows<-apply(dummy.matrix, 1, function(x) x[1]==del.edge[1] & x[2]==del.edge[2] )
      dummy.matrix<-dummy.matrix[del.rows!=1,]
    }
    
    #Apply R-PRR algorithm to new network after deletion of edges belonging to "order score"
    dummy.matrix<-as.data.table(dummy.matrix)
    dummy.matrix$JACCARD<-as.numeric(dummy.matrix$JACCARD)
    PRR.P.DISCOVER.DEL<-Function.PRR.PDISCOVER(met.diff, dummy.matrix, d=d, jaccard = jaccard, d.filter = d.filter)
    new.ranks<-which(PRR.P.DISCOVER.DEL$KEGG.ID %in% target.mets)
    print (new.ranks)
    
    #Obtain new delta
    delta<-mean(current.ranks) - mean(new.ranks)
    print (delta)
    
    #Evaluate conditions
    if (delta< -3 | sum(target.mets %in% unique(dummy.matrix$SUBSTRATE))<3){
      break
      
    } else if (delta>0){
      #We have improved due to deletion of edges: update matrix, delete corrupting edges and upadte to new rank
      network.matrix<-dummy.matrix
      del.paths<-c(del.paths, list(as.data.table(s.edges))) 
      current.ranks<-new.ranks
      
    } else if (delta<0 & delta>-3){
      #We have regressed in our ranks: should not update matrix, but store separatedly those edges that are benefitial to us (due to loss of rank upon their deletion)
      pro.paths<-c(pro.paths, list(as.data.table(s.edges)))
      
    } else if (delta==0){
      #Nothing has changed, not update of matrix or ranks
      network.matrix<-network.matrix
    }
    print (dim(network.matrix))
    print (current.ranks)
  }
  
  #Clean up and return resultant matrix, deleted paths, and beneficial paths
  PRO.PATHS<-do.call(rbind, pro.paths)
  DEL.PATHS<-do.call(rbind, del.paths)
  network.table<-as.data.table(network.matrix)
  return(list(NETWORK=network.table, DEL.PATHS=DEL.PATHS, PRO.PATHS=PRO.PATHS ))
}

Function.teru.met.matrix<-function(teru.clean, teru.kegg){
  
  #Clean up teru met matrix
  teru<-read.csv(teru.clean, header=F, stringsAsFactors=F)
  setnames(teru, as.vector(as.matrix(teru[2,])) )
  teru<-teru[4:nrow(teru),]
  rownames(teru)<-teru$STATUS
  teru$STATUS<-NULL
  colnames(teru)<-ifelse(colnames(teru)!="NORMAL", "CANCER", "NORMAL")
  teru<-teru[!(grepl("X - ", rownames(teru))),]
  
  teru<-data.matrix(teru)
  
  #Introduce kegg identifiers file
  teru.kegg<-fread(teru.kegg, header=F, sep=",", na.strings = "")
  teru.kegg<-teru.kegg[!is.na(V2),]
  teru.kegg$V3<-NULL
  setnames(teru.kegg, c("MET", "KEGG.ID"))
  teru.kegg<-teru.kegg[!grepl(",", KEGG.ID),]
  
  #Introduce kegg identifiers
  teru<-teru[(rownames(teru) %in% unique(teru.kegg$MET))==1,]
  rownames(teru)<-sapply(rownames(teru), function(x) unique(teru.kegg[MET==x,]$KEGG.ID))
  
  #Return
  return(teru)
}

Function.kegg.path.plot<-function(kegg.path, keywords=c("fatty acid biosynthesis", "glycolisis", "citrate"), met.matrix, teru.kegg){
  
  require(reshape2)
  
  #Prep path table
  kegg.path<-kegg.path[grepl("fatty acid biosynthesis", DESCRIPTION, ignore.case = T) | grepl("glycolysis", DESCRIPTION, ignore.case = T) | grepl("citrate", DESCRIPTION, ignore.case = T),]
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) strsplit(x, " - ")[[1]][1])
  kegg.path$DESCRIPTION<-ifelse(kegg.path$DESCRIPTION=="Citrate cycle (TCA cycle)", "TCA Cycle", kegg.path$DESCRIPTION) #minor fix
  kegg.path$PATHWAY<-NULL
  
  #Introduce available kegg identifiers
  met.matrix<-met.matrix[(rownames(met.matrix) %in% unique(kegg.path$COMPOUND))==1, ]
  
  teru.kegg<-fread(teru.kegg, header=F, sep=",", na.strings = "")
  teru.kegg<-teru.kegg[!is.na(V2),]
  teru.kegg$V3<-NULL
  setnames(teru.kegg, c("MET", "KEGG.ID"))
  teru.kegg<-teru.kegg[!grepl(",", KEGG.ID),]
  
  rownames(met.matrix)<-sapply(rownames(met.matrix), function(x) unique(teru.kegg[KEGG.ID==x,]$MET))
  setnames(kegg.path, c("KEGG.ID", "DESCRIPTION"))
  kegg.path<-merge(kegg.path, teru.kegg, by="KEGG.ID")
  
  #Return annotated heatmap
  colnames(met.matrix)<-paste(colnames(met.matrix), 1:ncol(met.matrix), sep=".")
  anno.col<-data.frame(type=sapply(colnames(met.matrix), function(x) strsplit(x,"[.]")[[1]][1] ),
                       row.names = colnames(met.matrix))
  
  anno.row<-data.frame(path=sapply(rownames(met.matrix), function(x) paste(unique(kegg.path[MET==x,]$DESCRIPTION), collapse=".")),
                   row.names=rownames(met.matrix))
  
  met.matrix<-met.matrix[order(anno.row$path),]
  pheatmap(log2(met.matrix), scale="none", size=14, cluster_cols = F, cluster_rows=F ,annotation_row = anno.row, annotation_col = anno.col, fontsize = 12,
           color=colorRampPalette(c("white", "green","green4","violet","purple"))(100))
}

