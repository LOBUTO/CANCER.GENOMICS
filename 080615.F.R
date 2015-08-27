Function.Tang<-function(tang.clean, teru.kegg, kegg.extra){
  #Constructs clean Tang matrix
  
  #Load tang matrix
  tang<-read.csv(tang.clean, header=F, stringsAsFactors=F)
  setnames(tang, as.vector(as.matrix(tang[1,])) )
  tang<-tang[2:nrow(tang),]
  rownames(tang)<-tang$METABOLITE
  tang$METABOLITE<-NULL
  colnames(tang)<-sapply(colnames(tang), function(x) {
    if (x=="Normal Breast"){
      cn<-"NORMAL"
    } else {
      cn<-paste0(c("TCGA",strsplit(x, "-")[[1]], "01A"),collapse = ".")
    }
    return(cn)
  })
  
  tang<-data.matrix(tang)
  rownames(tang)<-toupper(rownames(tang))
  
  #Load kegg identifiers
  tang.kegg<-fread(teru.kegg, header=F, sep=",", na.strings = "")
  tang.kegg<-tang.kegg[!is.na(V2),]
  tang.kegg$V3<-NULL
  setnames(tang.kegg, c("MET", "KEGG.ID"))
  tang.kegg$MET<-toupper(tang.kegg$MET)
  
  #Load extra kegg identifiers
  extra.kegg<-fread(kegg.extra, head=F, sep="\t", stringsAsFactors = F)
  setnames(extra.kegg, c("MET", "KEGG.ID"))
  extra.kegg<-extra.kegg[KEGG.ID!="NONE",]
  
  tang.kegg<-unique(rbind(tang.kegg, extra.kegg))
  
  #Use first annottated identifier
  tang.kegg<-tang.kegg[,list(KEGG.ID=KEGG.ID[1]), by="MET"]
  
  #Apply identifiers
  mets<-rownames(tang)[rownames(tang) %in% tang.kegg$MET]
  tang<-tang[mets,]
  rownames(tang)<-sapply(rownames(tang), function(x) unique(tang.kegg[MET==x,]$KEGG.ID))
  
  return(tang)
}

Function.teru.met.matrix<-function(teru.met.file, teru.id.to.gsm, teru.gene.obj, teru.kegg, kegg.extra, CANCER=T, teru.id.to.er){
  #Obtain terunuma correlation matrix depending on cancer/non-cancer patients
  
  #Prep met expression table
  x<-fread(teru.met.file, header=T, skip=2, stringsAsFactors = F)
  x<-data.frame(x, row.names = "ID")
  x<-data.matrix(x)
  y<-fread(teru.id.to.gsm, header=F, stringsAsFactors = F)
  colnames(x)<-sapply(colnames(x), function(x) y[V1==strsplit(x, "X")[[1]][2],]$V2)
  met.matrix<-x[,sapply(colnames(x), function(x) x!="character(0)")]
  
  #Filter by cancer/non-cancer
  teru.gene.obj<-readRDS(teru.gene.obj)
  if (CANCER==F){
    met.matrix<-met.matrix[,teru.gene.obj$CLASS[CLASS=="Normal",]$SAMPLE]  
  } else{
    met.matrix<-met.matrix[,teru.gene.obj$CLASS[CLASS=="Tumor",]$SAMPLE]
  }
  rownames(met.matrix)<-toupper(rownames(met.matrix))
  
  #Load kegg identifiers
  tang.kegg<-fread(teru.kegg, header=F, sep=",", na.strings = "")
  tang.kegg<-tang.kegg[!is.na(V2),]
  tang.kegg$V3<-NULL
  setnames(tang.kegg, c("MET", "KEGG.ID"))
  tang.kegg$MET<-toupper(tang.kegg$MET)
  
  #Load extra kegg identifiers
  extra.kegg<-fread(kegg.extra, head=F, sep="\t", stringsAsFactors = F)
  setnames(extra.kegg, c("MET", "KEGG.ID"))
  extra.kegg<-extra.kegg[KEGG.ID!="NONE",]
  
  tang.kegg<-unique(rbind(tang.kegg, extra.kegg))
  
  #Use first annottated identifier
  tang.kegg<-tang.kegg[,list(KEGG.ID=KEGG.ID[1]), by="MET"]
  
  #Apply identifiers
  mets<-rownames(met.matrix)[rownames(met.matrix) %in% tang.kegg$MET]
  met.matrix<-met.matrix[mets,]
  rownames(met.matrix)<-sapply(rownames(met.matrix), function(x) unique(tang.kegg[MET==x,]$KEGG.ID))
  
  #Obtain ER.STATUS if cancer
  if (CANCER==T){
    id.to.gsm<-fread(teru.id.to.gsm, header=F, stringsAsFactors = F)
    id.to.er<-fread(teru.id.to.er, header=F, stringsAsFactors = F, select=c(1,4))
    id.er<-merge(id.to.gsm, id.to.er, by="V1")
    id.er$V1<-NULL
    setnames(id.er, c("SAMPLE", "ER.STATUS"))
    id.er<-id.er[grepl("TUMOR", ER.STATUS),]
    id.er$ER.STATUS<-ifelse(id.er$ER.STATUS=="NEG TUMOR", "NEG", "POS")
    
    #Return
    return(list(MATRIX=met.matrix, ER.STATUS=id.er))
    
  } else {
    #Return
    return(list(MATRIX=met.matrix))
  }
}

Function.met.gene.pre.classify.teru<-function(met.matrix, gene.matrix, met.target=c(), met.th=2){
  #Build Terunuma gene/met matrix for classification
  
  #Find common patients
  common.samples<-intersect(colnames(met.matrix), colnames(gene.matrix))
  
  #Filter for target met and apply threshold
  met.matrix<-met.matrix[met.target,,drop=F]
  met.matrix<-ifelse(met.matrix<2, 0, 1)
  rownames(met.matrix)<-met.target
  print (met.matrix)
  
  #Combine matrices
  all.names<-c(rownames(gene.matrix), met.target)
  main.matrix<-rbind(gene.matrix[,common.samples], met.matrix[,common.samples])
  rownames(main.matrix)<-all.names
  
  #Clean up and return
  main.matrix<-t(main.matrix)
  return(main.matrix)
}

Function.tcga.mut.prep<-function(tcga.file, mut.filter=0){
  #Preps tcga mutation file for processing, simple mapping
  
  #Load file
  x<-fread(tcga.file, header=T, sep="\t",stringsAsFactors = F, select = c(1, 9,10, 16))
  
  #Fix sample identifier to 4 letter code
  x$Tumor_Sample_Barcode<-sapply(x$Tumor_Sample_Barcode, function(y) paste(strsplit(y, "-")[[1]][1:4], collapse = ".") )
  setkey(x)
  x<-unique(x)
  
  #Classify mutations into gain-of-function (GOF) or loss-of-function (LOF)
  x$MUTATION<-ifelse(x$Variant_Classification=="Silent", "SILENT",
                     ifelse(x$Variant_Classification %in% c("Missense_Mutation","RNA"), "GOF", "LOF"))
  x$MUTATION<-ifelse(x$Variant_Type %in% c("INS", "DEL"), "LOF", x$MUTATION)
  
  #Clean up names
  x$Variant_Classification<-NULL
  x$Variant_Type<-NULL
  setnames(x, c("Hugo_Symbol", "SAMPLE", "MUTATION"))
  
  #Filter samples by mutation number for GOF and LOF
  x.mut<-x[MUTATION %in% c("GOF", "LOF"),]
  x.sil<-x[MUTATION=="SILENT",]
  
  x.mut[,MUT.COUNT:=length(unique(SAMPLE)), by="Hugo_Symbol"] #At least number of samples to have mutated gene
  x.mut<-x.mut[MUT.COUNT>=mut.filter,]
  x.mut$MUT.COUNT<-NULL
  
  x<-rbind(x.sil, x.mut)
  
  #Clean up
  x<-x[!grepl("-Sep", Hugo_Symbol),][!grepl("-Mar", Hugo_Symbol),]
  
  #Return
  return(x)
}

Function.kegg.filtered<-function(kegg.file, recon.met.file, weight.filter=1000, n.edge.filter=25){
  #Builds the kegg directed network and filters for reactions between a substrate and product 
  #   that posses less than a threshold number of participating genes
  
  #Load file
  kegg.table<-fread(kegg.file, header=T, sep="\t", stringsAsFactors = F)
  setkey(kegg.table)
  kegg.table<-unique(kegg.table)
  
  #Manually insert reaction information (for R-2HG)
  kegg.table<-rbind(kegg.table, data.table(Hugo_Symbol=c("ADHFE1", "ADHFE1","ADHFE1", "ADHFE1",
                                                         "PHGDH", "PHGDH", "PHGDH", "PHGDH",
                                                         "D2HGDH",
                                                         "EGLN1", "EGLN2", "EGLN3"),
                                           SUBSTRATE=c("C03197","C03197", "C00026", "C00026",
                                                       "C01087","C00197","C01087", "C00197",
                                                       "C01087", 
                                                       "C00026", "C00026", "C00026"),
                                           PRODUCT=c("C00164", "C01087", "C00164", "C01087",
                                                     "C03232","C03232","C00026","C00026",
                                                     "C00026",
                                                     "C00042", "C00042", "C00042")))
  
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
  
  #Filter by number of genes present in edge
  kegg.table[,EDGE.N:=length(unique(Hugo_Symbol)), by=c("SUBSTRATE", "PRODUCT")]
  kegg.table<-kegg.table[EDGE.N<n.edge.filter,]
  kegg.table$EDGE.N<-NULL
  
  #Clean up and return
  return(kegg.table)
}

Function.tcga.brca.clinical<-function(tcga.file){
  
  #Load file and extract staging info
  clinical<-fread(tcga.file, sep="\t", header=F, stringsAsFactors = F, skip=3, select = c(2, 41, 44))
  setnames(clinical, c("SAMPLE", "STAGE", "ER.STATUS"))
  
  #Clean up
  clinical<-clinical[!(STAGE %in% c("[Not Available]", "Stage Tis", "Stage X")),]
  clinical<-clinical[!(ER.STATUS %in% c("[Not Evaluated]", "Indeterminate")),]
  
  #Combined stage per level
  clinical$STAGE<-ifelse(clinical$STAGE %in% c("Stage I", "Stage IA", "Stage IB"), "Stage.I", 
                         ifelse(clinical$STAGE %in% c("Stage II", "Stage IIA", "Stage IIB"), "Stage.II",
                                ifelse(clinical$STAGE %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), "Stage.III", 
                                       ifelse(clinical$STAGE %in% c("Stage IV", "Stage IVA"), "Stage.IV", clinical$STAGE))))
  
  #Adapt SAMPLE names to TCGA data
  clinical$SAMPLE<-sapply(clinical$SAMPLE, function(x) paste(strsplit(x, "-")[[1]], collapse="."))
  clinical<-clinical[,list(SAMPLE=c(paste(SAMPLE, "01A", sep="."), paste(SAMPLE,"01B", sep="."))), by=c("STAGE", "ER.STATUS")]
  
  #Return
  return(clinical)
}

Function.met.gold.std<-function(teru.normal.obj, teru.cancer.obj, tang.matrix, tcga.clinical, p.val.th=0.1, cutoff.th=10, lfc.th=2){
  #Finds differentially expressed metabolites between cancer and normal in both terunuma and tang datasets with respect to ER status
  #cutoff.th is the top number of metabolites to return per ER class ranked based on the average absolute log fold between datasets
  
  #First find differentially expressed metabolites in Tang with respect to er
  tang.er.pos<-intersect(colnames(tang.matrix), unique(tcga.clinical[ER.STATUS=="Positive",]$SAMPLE))
  tang.er.neg<-intersect(colnames(tang.matrix), unique(tcga.clinical[ER.STATUS=="Negative",]$SAMPLE))
  colnames(tang.matrix)<-colnames(data.frame(tang.matrix))
  tang.normal<-colnames(tang.matrix)[grepl("NORMAL", colnames(tang.matrix))]
  
  tang.pos<-data.table(MET=rownames(tang.matrix), PVAL=apply(tang.matrix, 1, function(x)  wilcox.test(x[tang.normal], x[tang.er.pos])$p.value),
                       FOLD=apply(tang.matrix, 1, function(x) log2(median(x[tang.er.pos])/median(x[tang.normal]))) )
  tang.pos$PVAL.ADJ<-p.adjust(tang.pos$PVAL, method="fdr")
  tang.pos<-unique(tang.pos[PVAL.ADJ<p.val.th,])
  tang.neg<-data.table(MET=rownames(tang.matrix), PVAL=apply(tang.matrix, 1, function(x)  wilcox.test(x[tang.normal], x[tang.er.neg])$p.value),
                       FOLD=apply(tang.matrix, 1, function(x) log2(median(x[tang.er.neg])/median(x[tang.normal]))) )
  tang.neg$PVAL.ADJ<-p.adjust(tang.neg$PVAL, method="fdr")
  tang.neg<-unique(tang.neg[PVAL.ADJ<p.val.th,])
  
  #Then for Terunuma
  teru.er.pos<-intersect(colnames(teru.cancer.obj$MATRIX), teru.cancer.obj$ER.STATUS[ER.STATUS=="POS",]$SAMPLE)
  teru.er.neg<-intersect(colnames(teru.cancer.obj$MATRIX), teru.cancer.obj$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE)
  teru.normal<-colnames(teru.normal.obj$MATRIX)
  
  teru.mets<-intersect(rownames(teru.normal.obj$MATRIX), rownames(teru.cancer.obj$MATRIX))
  teru.matrix<-cbind(teru.normal.obj$MATRIX[teru.mets,], teru.cancer.obj$MATRIX[teru.mets,])
  
  teru.pos<-data.table(MET=teru.mets, PVAL=apply(teru.matrix,1, function(x) wilcox.test(x[teru.normal], x[teru.er.pos])$p.value),
                       FOLD=apply(teru.matrix, 1, function(x) log2(median(x[teru.er.pos])/median(x[teru.normal]))) )
  teru.pos$PVAL.ADJ<-p.adjust(teru.pos$PVAL, method="fdr")
  teru.pos<-unique(teru.pos[PVAL.ADJ<p.val.th,])
  teru.neg<-data.table(MET=teru.mets, PVAL=apply(teru.matrix,1, function(x) wilcox.test(x[teru.normal], x[teru.er.neg])$p.value),
                       FOLD=apply(teru.matrix, 1, function(x) log2(median(x[teru.er.neg])/median(x[teru.normal]))) )
  teru.neg$PVAL.ADJ<-p.adjust(teru.neg$PVAL, method="fdr")
  teru.neg<-unique(teru.neg[PVAL.ADJ<p.val.th,])
  
  #Combine matrices and find intersection of differentially expressed metabolites from both sets
  #Top lfc across both datasets per ER class
  MET.POS<-merge(teru.pos, tang.pos, by="MET")
  MET.POS$AVE.LFC<-(MET.POS$FOLD.x + MET.POS$FOLD.y)/2
  MET.NEG<-merge(teru.neg, tang.neg, by="MET")
  MET.NEG$AVE.LFC<-(MET.NEG$FOLD.x + MET.NEG$FOLD.y)/2
  
  if (cutoff.th==Inf){
    MET.POS<-MET.POS[order(abs(AVE.LFC), decreasing = T),]$MET
    MET.NEG<-MET.NEG[order(abs(AVE.LFC), decreasing = T),]$MET  
  } else if (cutoff.th=="LFC") {
    MET.POS<-MET.POS[order(abs(AVE.LFC), decreasing = T),][abs(AVE.LFC)>=lfc.th,]$MET
    MET.NEG<-MET.NEG[order(abs(AVE.LFC), decreasing = T),][abs(AVE.LFC)>=lfc.th,]$MET
  } else {
    MET.POS<-MET.POS[order(abs(AVE.LFC), decreasing = T),][1:cutoff.th,]$MET
    MET.NEG<-MET.NEG[order(abs(AVE.LFC), decreasing = T),][1:cutoff.th,]$MET  
  }
  
  #Clean up and return
  MAIN.TABLE=merge(data.table(INDEX=MET.POS, MET.POS=MET.POS), data.table(INDEX=MET.NEG, MET.NEG=MET.NEG), by="INDEX", all=T)
  return(MAIN.TABLE)
}

Function.met.score.1<-function(tcga.mut, kegg.edges, aa.length, tcga.clinical){
  #Initial score will not discriminate between GOF/LOF  mutations
  #Scores based on ER-/+ subpopulations
  
  internal.scores<-function(mut){
    
    #Only consider silent mutations, for the time being, treat GOF or LOF the same
    mut<-mut[MUTATION!="SILENT",]
    
    #Normalize genes by amino acid length
    main.table<-merge(mut, aa.length, by="Hugo_Symbol")
    main.table[,MUT.SAMPLE.COUNT:=length(Hugo_Symbol), by=c("SAMPLE", "Hugo_Symbol")]
    main.table$MUT.NORM<-main.table$MUT.SAMPLE.COUNT/main.table$LENGTH
    
    #Normalize genes per patient mutation weight [Hugo_Symbol, SAMPLE, MUT.NORM]
    main.table[,SAMPLE.WEIGHT:=sum(MUT.NORM), by="SAMPLE"]
    main.table<-main.table[,list(MUT.NORM=unique(MUT.NORM)/unique(SAMPLE.WEIGHT)), by=c("SAMPLE", "Hugo_Symbol")]
    
    #Simplify kegg table [Hugo_Symbol, MET]
    kegg<-unique(rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                       data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT)))
    
    #SCORING
    #NOTE: We will do median for now of gene scores across patients
    main.table<-main.table[,list(MED.SCORE=median(MUT.NORM)), by="Hugo_Symbol"] #[Hugo_Symbol, MED.SCORE]
    
    #NOTE: Median across genes per metabolite
    main.table<-rbind(main.table, data.table(Hugo_Symbol= setdiff(unique(kegg$Hugo_Symbol), main.table$Hugo_Symbol), MED.SCORE=0))
    main.table<-merge(kegg, main.table, by="Hugo_Symbol")
    main.table<-main.table[,list(MET.SCORE=median(MED.SCORE)), by="MET"]
    
    #Clean up and return
    main.table<-main.table[order(MET.SCORE, decreasing = T),]
    return(main.table)
  }
  
  #Do for both er+/-
  mut.pos<-tcga.mut[SAMPLE %in% unique(tcga.clinical[ER.STATUS=="Positive"]$SAMPLE),]
  mut.neg<-tcga.mut[SAMPLE %in% unique(tcga.clinical[ER.STATUS=="Negative"]$SAMPLE),]
  
  mut.pos<-internal.scores(mut.pos)
  mut.neg<-internal.scores(mut.neg)
  
  #Return both
  return(list(MUT.POS=mut.pos, MUT.NEG=mut.neg))
}

Function.met.pr<-function(score.table, target.mets){
  #Score table should have columns "MET" and "MET.SCORE"
  
  require(ROCR)
  require(pROC)
  
  #Prep data
  score.table<-score.table[order(MET.SCORE,decreasing = T),]
  pred.scores<-score.table$MET.SCORE
  pred.mets<-score.table$MET
  
  #Filter for those mets we can actually predict
  target.mets<-target.mets[target.mets %in% pred.mets]
  
  #Calculate tpr and fpr 
  target.pred<-prediction(pred.scores, pred.mets %in% target.mets)
  target.perf<-performance(target.pred, "tpr", "fpr")
  
  #Convert to ready to ggplot format
  pred.table<-data.table(FPR=target.perf@x.values[[1]], TPR=target.perf@y.values[[1]])
  pred.auc<-auc(pred.mets %in% target.mets, pred.scores)
  
  #Return
  return(list(PRED.TABLE=pred.table, AUC=pred.auc)) 
}

Function.met.score.2<-function(exp.obj, kegg.edges, tcga.clinical, pval.th=0.1, fold.th=1){
  #Calculate metabolite score based on expression data
  #Initial score will not discriminate between GOF/LOF  mutations
  #Scores based on ER-/+ subpopulations
  #NOTE: exp.obj is object containing both normal and cancer matrices
  
  internal.scores<-function(exp.samples){
    
    #Filter for expression samples in cancer matrix
    cancer.matrix<-exp.obj$tumor[,intersect(colnames(exp.obj$tumor), exp.samples)]
    
    #Combine normal and cancer matrices by genes of interest
    common.genes<-intersect(rownames(cancer.matrix), unique(kegg.edges$Hugo_Symbol))
    main.matrix<-cbind(exp.obj$normal[common.genes,], cancer.matrix[common.genes,])
    
    #Calculate differential expresion for hugos of interest - WITHOUT Normalization
    normal.samples<-colnames(exp.obj$normal)
    cancer.samples<-colnames(cancer.matrix)
    gene.pvals<-apply(main.matrix, 1, function(x) wilcox.test(x[cancer.samples], x[normal.samples])$p.value)
    gene.folds<-apply(main.matrix, 1, function(x) log2(mean(x[cancer.samples])/mean(x[normal.samples])))
    
    #Clean up gene table [Hugo_Symbol, PVAL.ADJ, LFC]
    main.table<-data.table(Hugo_Symbol=rownames(main.matrix), PVAL=gene.pvals, LFC=gene.folds)
    main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
    
    #Obtain met score based on significance cut off - PVAL and ABSOLUTE value of log fold change
    main.table$PASS<-main.table$PVAL.ADJ<pval.th & abs(main.table$LFC)>=fold.th
    main.table$PASS<-as.numeric(main.table$PASS)
    
    kegg<-unique(rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                       data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT)))
    
    main.table<-merge(main.table, kegg, by="Hugo_Symbol")
    main.table<-main.table[,list(MET.SCORE=mean(PASS)), by="MET"]
    
    #Clean up and return
    main.table<-main.table[order(MET.SCORE, decreasing = T),]
    return(main.table)
  }
  
  #Pre-filter expression matrices for non-variable genes (sd==0)
  exp.obj$tumor<-data.matrix(exp.obj$tumor)
  exp.obj$tumor<-exp.obj$tumor[apply(exp.obj$tumor, 1, sd)!=0,]
  exp.obj$normal<-data.matrix(exp.obj$normal)
  exp.obj$normal<-exp.obj$normal[apply(exp.obj$normal, 1, sd)!=0,]
  common.genes<-intersect(rownames(exp.obj$tumor), rownames(exp.obj$normal))
  exp.obj$tumor<-exp.obj$tumor[common.genes,]
  exp.obj$normal<-exp.obj$normal[common.genes,]
  
  #Do for both er+/-
  pos.samples<-unique(tcga.clinical[ER.STATUS=="Positive",]$SAMPLE)
  neg.samples<-unique(tcga.clinical[ER.STATUS=="Negative",]$SAMPLE)
  
  exp.pos<-internal.scores(pos.samples)
  exp.neg<-internal.scores(neg.samples)
  
  #Return both
  return(list(EXP.POS=exp.pos, EXP.NEG=exp.neg))
}

Function.met.score.3<-function(tcga.mut, exp.obj, tcga.clinical, kegg.edges, pval.th=0.1, fold.th=1, met.sample.th=5){
  #Calculate metabolite score based on BOTH mutation and expression data. CONTRIBUTION of mutation to phenotype.
  #   NOTE: Contribution is determined based on number of differentially expression genes for samples that have metabolite mutations vs rest of samples in cancer
  #Initial score will not discriminate between GOF/LOF  mutations
  #Scores based on ER-/+ subpopulations
  #NOTE: exp.obj is object containing both normal and cancer matrices
  
  require(parallel)
  require(data.table)
  require(reshape2)
  
  #Set up parallelization
  print ("Importing values for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("tcga.mut", "exp.obj", "tcga.clinical", "pval.th", "fold.th", "met.sample.th", "kegg.edges",
                              "as.data.table", "data.table"),envir=environment())
  
  internal.scores<-function(exp.samples){
    
    #Consolidate mutation data
    kegg<-unique(rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                       data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT)))
    
    ###Obtain score for each metabolite based on phenotype contribution to exp.samples cohort###
    #Given that there are actually metabolites with associated mutation
    merged.kegg<-merge(kegg, tcga.mut, by="Hugo_Symbol")
    merged.kegg<-merged.kegg[SAMPLE %in% exp.samples,] #Based on SAMPLE COHORT
    
    #Filter for those metabolites that have samples in more than met.sample.th
    merged.kegg[,COUNT.TH:=length(unique(SAMPLE)), by="MET"]
    merged.kegg<-merged.kegg[COUNT.TH>=met.sample.th,]
    all.mets<-unique(merged.kegg$MET)
    
    count<-0
    met.scores<-sapply(all.mets, function(x) {
      
      #Separate those samples affected by metabolite associated mutations and the rest
      met.samples<-unique(merged.kegg[MET==x,]$SAMPLE)
      rest.samples<-setdiff(exp.samples, met.samples)
      
      #Calculate number of differentially expressed genes between the 2 groups (Phenotype contribution)
      gene.pvals<-parApply(cl, exp.obj$tumor, 1, function(y) wilcox.test(y[met.samples], y[rest.samples])$p.value)
      gene.lfc<-parApply(cl, exp.obj$tumor, 1, function(y) log2(mean(y[met.samples])/mean(y[rest.samples])))
      met.table<-data.table(PVAL=gene.pvals, LFC=gene.lfc)
      met.table$PVAL.ADJ<-p.adjust(met.table$PVAL, method="fdr")
      
      #Return score
      met.table<-met.table[!(is.na(LFC)),][!(is.na(PVAL.ADJ)),] #NOTE: Temporary clean up to account for inconsistencies!
      met.score<-mean(met.table$PVAL.ADJ<pval.th & abs(met.table$LFC)>=fold.th)
      
      count<<-count+1
      print(count/length(all.mets))
      return(met.score)
    })
    
    #Clean up and return
    main.table<-data.table(MET=all.mets, MET.SCORE=met.scores)
    main.table<-main.table[order(MET.SCORE, decreasing = T),]
    return(main.table)
  }
  
  #Only consider silent mutations, for the time being, treat GOF or LOF the same
  tcga.mut<-tcga.mut[MUTATION!="SILENT",]
  
  #Pre-filter expression matrices for non-variable genes (sd==0)
  exp.obj$tumor<-data.matrix(exp.obj$tumor)
  exp.obj$tumor<-exp.obj$tumor[apply(exp.obj$tumor, 1, sd)!=0,]
  
  #Do for both er+/- For those that have both expression and mutation data
  common.samples<-intersect(unique(tcga.mut$SAMPLE), colnames(exp.obj$tumor))
  pos.samples<-intersect(unique(tcga.clinical[ER.STATUS=="Positive",]$SAMPLE), common.samples)
  neg.samples<-intersect(unique(tcga.clinical[ER.STATUS=="Negative",]$SAMPLE), common.samples)
  
  exp.pos<-internal.scores(pos.samples)
  exp.neg<-internal.scores(neg.samples)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return both
  return(list(EXP.POS=exp.pos, EXP.NEG=exp.neg)) 
}

Function.gene.diff.exp<-function(exp.obj, target.samples, type="tcga"){
  #Obtain log fold change per gene of cancer vs normal
  require(limma)
  
  #Prep main matrix
  if (type=="tcga"){
    normal.samples<-colnames(exp.obj$normal)
    cancer.samples<-intersect(colnames(exp.obj$tumor), target.samples)  
    common.genes<-intersect(rownames(exp.obj$normal), rownames(exp.obj$tumor))
    main.matrix<-cbind(exp.obj$normal[common.genes, ], exp.obj$tumor[common.genes, cancer.samples])
    
  } else if (type=="teru"){
    normal.samples<-exp.obj$CLASS[CLASS=="Normal",]$SAMPLE
    cancer.samples<-intersect(exp.obj$CLASS[CLASS=="Tumor",]$SAMPLE, target.samples)
    main.matrix<-normalizeBetweenArrays(exp.obj$MATRIX, method = "quantile")
    main.matrix<-2^main.matrix
  }
  
  #Calculate log fold expression
  lfc<-apply(main.matrix, 1, function(x) log2(median(x[cancer.samples])/median(x[normal.samples])))
  
  #Clean up and return
  main.table<-data.table(Hugo_Symbol=rownames(main.matrix), LFC=lfc)
  main.table<-main.table[!is.na(LFC),][LFC!=Inf,][LFC!=-Inf,]
  main.table<-main.table[order(abs(LFC),decreasing = T),]
  return(main.table)
}

Function.met.diff.exp<-function(met.matrix, target.samples, type="tang", normal.matrix){
  #Calculates wether metabolite is differentially expressed or not 
  #NOTE: type is either "teru" or "tang". If type is "teru", then normal matrix needs to be provided
  
  #Prep data sets
  if (type=="tang"){
    colnames(met.matrix)<-colnames(data.frame(met.matrix))
    normal.samples<-colnames(met.matrix)[grepl("NORMAL", colnames(met.matrix))]
    cancer.samples<-intersect(setdiff(colnames(met.matrix), normal.samples) , target.samples)
    main.matrix<-met.matrix
    
  } else if (type=="teru"){
    normal.samples<-colnames(normal.matrix)
    cancer.samples<-intersect(colnames(met.matrix), target.samples)
    common.genes<-intersect(rownames(normal.matrix), rownames(met.matrix))
    main.matrix<-cbind(met.matrix[common.genes, cancer.samples], normal.matrix[common.genes,])
  }
  
  #Calculate log fold change
  pvals<-apply(main.matrix, 1, function(x) wilcox.test(x[cancer.samples], x[normal.samples])$p.value)
  lfc<-apply(main.matrix, 1, function(x) log2(median(x[cancer.samples])/median(x[normal.samples])))
  pvals.adj<-p.adjust(pvals, method="fdr")
  
  #Clean up and return
  main.table<-data.table(MET=rownames(main.matrix), PVAL.ADJ=pvals.adj, LFC=lfc)
  main.table<-main.table[!is.na(PVAL.ADJ),][PVAL.ADJ!=Inf,][PVAL.ADJ!=-Inf,]
  main.table<-main.table[!is.na(LFC),][LFC!=Inf,][LFC!=-Inf,]
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.met.gene.cor.diff<-function(met.pval.table, exp.obj, target.samples, kegg.edges, type="tcga"){
  #Find differences in gene expression correlation between normal and cancer samples for metabolites
  #NOTE: type is either "teru" or "tang"
  
  #Convert to single kegg.table
  kegg<-unique(data.table(MET=c(kegg.edges$SUBSTRATE, kegg.edges$PRODUCT), Hugo_Symbol=c(kegg.edges$Hugo_Symbol, kegg.edges$Hugo_Symbol)))
  
  #Prep main matrix
  if (type=="tcga"){
    normal.samples<-colnames(exp.obj$normal)
    cancer.samples<-intersect(colnames(exp.obj$tumor), target.samples)  
    common.genes<-intersect(rownames(exp.obj$normal), rownames(exp.obj$tumor))
    main.matrix<-cbind(exp.obj$normal[common.genes, ], exp.obj$tumor[common.genes, cancer.samples])
    
  } else if (type=="teru"){
    normal.samples<-exp.obj$CLASS[CLASS=="Normal",]$SAMPLE
    cancer.samples<-intersect(exp.obj$CLASS[CLASS=="Tumor",]$SAMPLE, target.samples)
    main.matrix<-exp.obj$MATRIX
    #main.matrix<-normalizeBetweenArrays(exp.obj$MATRIX, method = "quantile")
  }
  
  #Filter kegg table for mets of interest and filter out those mets that have less than 2 associated genes
  kegg<-kegg[MET %in% unique(met.pval.table$MET),]
  common.hugos<-intersect(unique(kegg$Hugo_Symbol), rownames(main.matrix))
  kegg<-kegg[Hugo_Symbol %in% common.hugos,]
  kegg[,H.COUNT:=length(unique(Hugo_Symbol)), by="MET"]
  kegg<-kegg[H.COUNT>1,]
  kegg$H.COUND<-NULL
  
  #Filter expression matrix for genes of interest
  kegg.hugos<-unique(kegg$Hugo_Symbol)
  main.matrix<-main.matrix[kegg.hugos,]
  
  #Obtain cancer and normal correlations
  cancer.cor<-cor(t(main.matrix[,cancer.samples]), method = "spearman")
  normal.cor<-cor(t(main.matrix[,normal.samples]), method = "spearman")
  
  #Extract mean correlation difference
  diag(cancer.cor)<-0
  diag(normal.cor)<-0
  target.mets<-unique(kegg$MET)
  met.gene.cor<-kegg[,list(MEAN.COR.DIFF=  mean(cancer.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(cancer.cor[Hugo_Symbol, Hugo_Symbol])]) - 
                             mean(normal.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(normal.cor[Hugo_Symbol, Hugo_Symbol])])),
                     by="MET"]
  
  #Add zero difference for mets that had only one gene associated, since we could not calculate correlations for it
  met.gene.cor<-met.gene.cor[!is.na(MEAN.COR.DIFF),]
  main.table<-rbind(met.gene.cor, data.table(MET=setdiff(met.pval.table$MET, met.gene.cor$MET), MEAN.COR.DIFF=0))
  
  #ABSOLUTE MEAN.COR.DIFF
  main.table$ABS.MEAN.COR.DIFF<-abs(main.table$MEAN.COR.DIFF)
  main.table$MEAN.COR.DIFF<-NULL
  
  #Clean up and return
  main.table<-main.table[order(abs(ABS.MEAN.COR.DIFF), decreasing = T),]
  return(main.table)
}

Function.prep.kegg.pred.table<-function(kegg.edges, gene.diff.exp, met.diff.exp, kegg.path, met.gene.mcd, pval.th=0.05, lfc.th=1){
  #Constructs table for classification of differetially expressed metabolites based on thresholds
  
  require(igraph)
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Obtain metabolite degree (based on both in and out degree)
  kegg.graph<-graph.data.frame(unique(kegg.edges[,2:3, with=F]), directed = T)
  kegg.degree<-data.table(as.matrix(degree(kegg.graph, mode = "all",normalized = T)), keep.rownames = T)
  setnames(kegg.degree, c("MET", "NORM.DEGREE"))
  
  #Obtain metabolite betweeness centrality
  kegg.between<-data.table(as.matrix(betweenness(kegg.graph, directed = T, normalized = T)), keep.rownames = T)
  setnames(kegg.between, c("MET", "NORM.BC"))
  
  #Classify metabolites into differentially expressed or not (based on thresholds)
  met.sig<-met.diff.exp[PVAL.ADJ<pval.th & abs(LFC)>lfc.th,]$MET
  met.diff.exp$DIFF<-as.factor(ifelse(met.diff.exp$MET %in% met.sig, 1,0))
  
  #Obtain mean absolute lfc change of enzymes for each metabolite
  kegg.hugos<-unique(data.table(MET=c(kegg.edges$SUBSTRATE, kegg.edges$PRODUCT), Hugo_Symbol=c(kegg.edges$Hugo_Symbol, kegg.edges$Hugo_Symbol)))
  kegg.hugos<-merge(kegg.hugos, gene.diff.exp, by="Hugo_Symbol")
  kegg.hugos<-kegg.hugos[,list(HUGO.MED.LFC=mean(abs(LFC))), by="MET"]
  
  #Combine all data
  main.table<-merge(met.diff.exp[,c("MET", "DIFF"), with=F], kegg.degree, by="MET") #Add degree info
  main.table<-merge(main.table, kegg.hugos, by="MET") # Add Hugo LFC
  main.table<-merge(main.table, kegg.between, by="MET") # Add BC
  main.table<-merge(main.table, met.gene.mcd, by="MET") #Add met gene mean correlation difference 
  
  #Add path info - NOTE: For the time being it will be [# of paths, presence in Fats, Glyco, TCA and CANCER]
  kegg.path.count<-kegg.path[,list(PATH.COUNT=length(unique(DESCRIPTION))), by="COMPOUND"]
  kegg.path.count<-rbind(kegg.path.count, 
                         data.table(COMPOUND=setdiff(unique(main.table$MET), kegg.path.count$COMPOUND), PATH.COUNT=0)
                         )
  setnames(kegg.path.count, c("MET", "PATH.COUNT"))
  
  main.table$GLYCO.TCA.LIPID<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Glycol", DESCRIPTION),]$COMPOUND, 1,
                                     ifelse(main.table$MET %in% kegg.path[grepl("fatty acid", DESCRIPTION, ignore.case = T),]$COMPOUND, 1,
                                            ifelse(main.table$MET %in% kegg.path[grepl("TCA", DESCRIPTION),]$COMPOUND, 1, 0))))
  main.table$CANCER<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("cancer", DESCRIPTION, ignore.case = T),]$COMPOUND, 1, 
                                      ifelse(main.table$MET %in% kegg.path[grepl("Glioma", DESCRIPTION),]$COMPOUND, 1, 
                                             ifelse(main.table$MET %in% kegg.path[grepl("Melanoma", DESCRIPTION),]$COMPOUND, 1, 0))))
  main.table$STEROID<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Steroid", DESCRIPTION, ignore.case = T),]$COMPOUND, 1, 0))
  main.table$ABC<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("ABC transporters", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ANTIBIO<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Biosynthesis of antibiotics", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$AA<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Biosynthesis of amino acids", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$PROTEIN.ABS<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Protein digestion and absorption", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$AMINOACYL.TRNA<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Aminoacyl-tRNA biosynthesis", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$TWO.OXO<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("2-Oxocarboxylic acid metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$PURINE<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Purine metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$GST<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Glycine, serine and threonine", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$CARBON<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Carbon metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ARG.PRO<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Arginine and proline metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$ALA.ASP.GLU<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Alanine, aspartate and glutamate metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  main.table$CYS.MET<-as.factor(ifelse(main.table$MET %in% kegg.path[grepl("Cysteine and methionine metabolism", DESCRIPTION),]$COMPOUND, 1, 0))
  
  main.table<-merge(main.table, kegg.path.count, by="MET")
  
  #Normalize needed variables between 0 and 1
  main.table$HUGO.MED.LFC<-normalize.vector(main.table$HUGO.MED.LFC)
  main.table$ABS.MEAN.COR.DIFF<-normalize.vector(main.table$ABS.MEAN.COR.DIFF)
  main.table$PATH.COUNT<-normalize.vector(main.table$PATH.COUNT)
  
  #Return
  return(main.table)
}

Function.class.master.sampler<-function(target.samples, sample.size=10, n= 10, type="tcga", pval.th=0.05, lfc.th=1){
  #Resamples from population to obtain increased independent master.tables for classification
  #n is number of replications
  
  #Obtain random samples - Pseudobootstrapping
  random.list<-replicate(n, sample(target.samples, sample.size), simplify = F)
  
  #Apply functions
  master.table<-data.table()
  for (i in random.list){
    
    ####Depending on type (tcga/teru)###
    if (type=="tcga"){
      
      #Obtain gene log fold changes
      lfc.table<-Function.gene.diff.exp(brca.exp, i, type=type)  
      
      #Obtain met diff expression
      met.pval<-Function.met.diff.exp(tang.matrix, i, "tang")
      
      #Obtain met's gene mean correlation difference
      met.mcd<-Function.met.gene.cor.diff(met.pval, brca.exp, i, kegg.edges, type = "tcga")
      
    } else if (type=="teru"){
      
      #Obtain gene log fold changes
      lfc.table<-Function.gene.diff.exp(teru.gene.exp, i , type=type)
      
      #Obtain met diff expression
      met.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, i , "teru", teru.normal.matrix$MATRIX)
      
      #Obtain met's gene mean correlation difference
      met.mcd<-Function.met.gene.cor.diff(met.pval, teru.gene.exp, i, kegg.edges, type = "teru")
      
    }
    ###################################
    
    #Obtain class table and store
    class.table<-Function.prep.kegg.pred.table(kegg.edges, lfc.table, met.pval, kegg.path, met.mcd, pval.th = pval.th, lfc.th = lfc.th)
    master.table<-rbind(master.table, class.table)
  }
  
  #Return master table
  return (master.table)
}

Function.deep.met.learning<-function(met.class.table){
  
  library(h2o)
  
  #Initialize h2o node
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 
  
  #Prep data
  
}

Function.tcga.pred.met<-function(tcga.exp, h2o.model, target.samples, filter.pval.met, kegg.edges, kegg.path, localH2O){
  #Obtains prediction of significantly differentiated metabolites in cancer dataset
  #NOTE: 
  #   tcga.exp is object with "tumor" and "normal" expression matrices 
  #   filter.pval.met is table of type "tang.plus.pval.met" to filter for those metabolites that we initially had information on 
  #   h2o.model expected to be in the localh2o()
  #   localH2O is the h2o.init() object
  
  #Prep tcga prediction table
  print ("prepping data")
  table.lfc.gene<-Function.gene.diff.exp(tcga.exp, target.samples , type = "tcga")
  table.mcd<-Function.met.gene.cor.diff(filter.pval.met, tcga.exp, target.samples, kegg.edges, type = "tcga")
  class.table<-Function.prep.kegg.pred.table(kegg.edges, table.lfc.gene, filter.pval.met, kegg.path, table.mcd, pval.th = 0.05, lfc.th = 1)
  
  #Apply model in h2o mode
  print ("building predictions")
  class.table.h2o<-as.h2o(localH2O, class.table, "class.table")
  tcga.pred<-data.table(MET=class.table$MET, PRED=as.data.frame(h2o.predict(h2o.model, class.table.h2o))[["p1"]])
  h2o.rm(localH2O, c("class.table"))
  
  #Return
  return(tcga.pred)
}

Function.kegg.path.enrich<-function(kegg.path, sig.mets){
  
  #Remove "- Homo sapiens" from path description
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - Homo"))[1])
  
  #Remove "Metabolic pathways" from kegg.path
  kegg.path<-kegg.path[!grepl("metabolic", DESCRIPTION, ignore.case = T),]
  
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

Function.process.icgc.exp.raw<-function(icgc.matrix, icgc_info, target.cancer.samples){
  #Finds differential expression between groups of: "Primary Tumor" and "Normal"
  
  require(edgeR)
  require(qvalue)
  
  #Read and clean info file
  icgc_info<-fread(icgc_info, header=T)
  icgc_info<-icgc_info[specimen_type %in% c("Normal - tissue adjacent to primary", "Primary tumour - solid tissue"),]
  icgc_info$SAMPLE<-sapply(icgc_info$submitted_specimen_id, function(x) paste0(strsplit(x, "-")[[1]], collapse = "." ) )
  icgc_info<-icgc_info[,c("SAMPLE", "specimen_type"), with=F]
  
  #Filter by cancer samples
  icgc.normals<-icgc_info[specimen_type=="Normal - tissue adjacent to primary",]
  icgc.tumors<-icgc_info[specimen_type=="Primary tumour - solid tissue",]
  icgc.tumors<-icgc.tumors[SAMPLE %in% target.cancer.samples,]
  icgc_info<-rbind(icgc.normals, icgc.tumors)
  
  #Read and clean expression file
  cast.matrix<-readRDS(icgc.matrix)
  common.samples<-intersect(colnames(cast.matrix), unique(icgc_info$SAMPLE))
  cast.matrix<-cast.matrix[,common.samples]
  icgc_info<-icgc_info[SAMPLE %in% common.samples,]
  
  #Prep cpm, filter out genes that don't have at least 1 cpm in 1/4 of all samples 
  matrix_cpm<-cpm(cast.matrix) #For sample ordering in DGEList later
  min_count<-ncol(matrix_cpm)/4
  keep<-rowSums(matrix_cpm>1)>=min_count
  matrix_filt<-matrix_cpm[keep, icgc_info$SAMPLE]
  
  #Obtain differentially expressed genes
  icgc_DGE <- DGEList(matrix_filt, group = icgc_info$specimen_type)
  icgc_DGE <- calcNormFactors(icgc_DGE)
  icgc_DGE <- estimateCommonDisp(icgc_DGE)
  icgc_DGE <- estimateTagwiseDisp(icgc_DGE)
  icgc_ET <- exactTest(icgc_DGE)
  
  #Correct for multiple hypothesis testing
  main.table<-data.table(icgc_ET$table, keep.rownames = T)
  setnames(main.table, c("Hugo_Symbol", "LFC", "LOG.CPM", "PVAL"))
  main.table$PVAL.ADJ<-qvalue(main.table$PVAL)$qvalues
  
  #Clean up and return
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

Function.teru.diff.limma<-function(teru.obj, target.samples, norm=T){
  #Finds differentially expressed genes according to ebayesfit in the limma package after quantile normalization
  
  require(data.table)
  require(limma)
  require(edgeR)
  
  #Filter for target samples
  tumor.samples<-intersect(target.samples, teru.obj$CLASS[CLASS=="Tumor",]$SAMPLE)
  normal.samples<-teru.obj$CLASS[CLASS=="Normal",]$SAMPLE
  exp.matrix<-teru.obj$MATRIX[,c(normal.samples, tumor.samples)]
  
  #Normalize?
  if (norm==T){
    exp.matrix<-normalizeBetweenArrays(exp.matrix, method = "quantile")
  }
  
  #Obtain design matrix
  design<-cbind(WT=1, MTvsWT=colnames(exp.matrix) %in% tumor.samples )
  fit<-lmFit(exp.matrix, design)
  fit<-eBayes(fit)
  
  #Clean up and Return
  top.fit<-topTable(fit, coef = 2, number = 25000)
  top.fit<-data.table(top.fit, keep.rownames = T)
  top.fit<-top.fit[,c(1,2,6),with=F]
  setnames(top.fit, c("Hugo_Symbol", "LFC", "PVAL.ADJ"))
  top.fit<-top.fit[order(abs(LFC), decreasing = T),]
  return(top.fit)
}