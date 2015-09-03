library(h2o)
library(data.table)
library(reshape2)
library(caret)

Function.met.split<-function(limit.met, seed=23, pval.met, pval.th = 0.05, lfc.th = 1, folds=10){
  require(data.table)
  require(caret)
  
  set.seed(seed)
  
  #Filter table by limit.met
  pval.met<-pval.met[MET %in% limit.met,]
  met.sig<-pval.met[PVAL.ADJ<pval.th & abs(LFC)>lfc.th,]$MET
  pval.met$DIFF<-as.factor(ifelse(pval.met$MET %in% met.sig, 1,0))
  
  #Split table
  met.splits<-createFolds(pval.met$DIFF, folds)
  met.splits<-lapply(met.splits, function(x) pval.met$MET[x])
  
  #Return
  return(met.splits)
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
  kegg$H.COUNT<-NULL
  
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
  #   met.gene.cor<-kegg[,list(MEAN.COR.DIFF=  mean(cancer.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(cancer.cor[Hugo_Symbol, Hugo_Symbol])]) - 
  #                              mean(normal.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(normal.cor[Hugo_Symbol, Hugo_Symbol])])),
  #                      by="MET"]
  met.gene.cor<-kegg[,list(MEAN.CANCER.COR=mean(cancer.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(cancer.cor[Hugo_Symbol, Hugo_Symbol])]),
                           MEAN.NORMAL.COR=mean(normal.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(normal.cor[Hugo_Symbol, Hugo_Symbol])]),
                           RATIO.CANCER.COR.POS=mean(cancer.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(cancer.cor[Hugo_Symbol, Hugo_Symbol])]>0),
                           RATIO.CANCER.COR.NEG=mean(cancer.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(cancer.cor[Hugo_Symbol, Hugo_Symbol])]<0),
                           RATIO.NORMAL.COR.POS=mean(normal.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(normal.cor[Hugo_Symbol, Hugo_Symbol])]>0),
                           RATIO.NORMAL.COR.NEG=mean(normal.cor[Hugo_Symbol, Hugo_Symbol][upper.tri(normal.cor[Hugo_Symbol, Hugo_Symbol])]<0)),
                     by="MET"]
  
  #Add zero difference for mets that had only one gene associated, since we could not calculate correlations for it
  #   met.gene.cor<-met.gene.cor[!is.na(MEAN.COR.DIFF),]
  #   main.table<-rbind(met.gene.cor, data.table(MET=setdiff(met.pval.table$MET, met.gene.cor$MET), MEAN.COR.DIFF=0))  
  met.added<-met.gene.cor[is.na(MEAN.CANCER.COR),]
  if (nrow(met.added)>0){
    met.added<-data.table(MET=met.added$MET, MEAN.CANCER.COR=0, MEAN.NORMAL.COR=0, 
                          RATIO.CANCER.COR.POS=0, RATIO.CANCER.COR.NEG=0, 
                          RATIO.NORMAL.COR.POS=0, RATIO.NORMAL.COR.NEG=0)
    met.gene.cor<-met.gene.cor[!is.na(MEAN.CANCER.COR),]
    met.gene.cor<-rbind(met.gene.cor, met.added)  
  }
  
  #main.table<-rbind(met.gene.cor, data.table(MET=setdiff(met.pval.table$MET, met.gene.cor$MET), 
  
  #ABSOLUTE MEAN.COR.DIFF
  #   main.table$ABS.MEAN.COR.DIFF<-abs(main.table$MEAN.COR.DIFF)
  #   main.table$MEAN.COR.DIFF<-NULL
  
  #Clean up and return
  #   main.table<-main.table[order(abs(ABS.MEAN.COR.DIFF), decreasing = T),]
  return(met.gene.cor)
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
  main.table<-merge(main.table, met.gene.mcd, by="MET") #Add met gene mean correlation information
  
  #Add path info (if necessary) - NOTE: For the time being it will be [# of paths, presence in Fats, Glyco, TCA and CANCER]
  kegg.path.count<-kegg.path[,list(PATH.COUNT=length(unique(DESCRIPTION))), by="COMPOUND"]
  if (length(setdiff(unique(main.table$MET), kegg.path.count$COMPOUND))>0){
    kegg.path.count<-rbind(kegg.path.count, 
                           data.table(COMPOUND=setdiff(unique(main.table$MET), kegg.path.count$COMPOUND), PATH.COUNT=0)
    )
  }
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
  #main.table$ABS.MEAN.COR.DIFF<-normalize.vector(main.table$ABS.MEAN.COR.DIFF)
  main.table$PATH.COUNT<-normalize.vector(main.table$PATH.COUNT)
  
  #Return
  return(main.table)
}


#Arguments
args<-commandArgs(trailingOnly=T)
teru.plus.mcd<-readRDS(args[1])
teru.plus.pval.met<-readRDS(args[2])
teru.normal.matrix<-readRDS(args[3])
teru.cancer.matrix<-readRDS(args[4])
teru.gene.exp<-readRDS(args[5])
kegg.edges<-readRDS(args[6])
kegg.path<-readRDS(args[7])
h2o.folder<-args[8]

output.file<-args[9]


plus.met.splits<-Function.met.split(teru.plus.mcd$MET, 23, teru.plus.pval.met, folds = 10)
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 

metrics.bootstrap.2<-data.table()
for (i in 1:length(plus.met.splits)){
  
  #Split mets
  testing.mets<-plus.met.splits[[i]]
  training.mets<-setdiff(unlist(plus.met.splits), testing.mets)
  
  #Obtain bootstrapped samples
  er.plus.samples<-teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE
  bootstrap.train<-replicate(10, sample(er.plus.samples, 10), simplify = F)
  
  #Create master table per bootsrapped set (for bootstrapped samples and training mets)
  master.table<-data.table()
  
  for (j in bootstrap.train){
    
    train.limma<-Function.teru.diff.limma(teru.gene.exp, j, norm = F)
    train.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  j , "teru", teru.normal.matrix$MATRIX)
    train.mcd<-Function.met.gene.cor.diff(train.pval, teru.gene.exp, j, kegg.edges, type = "teru")
    train.table<-Function.prep.kegg.pred.table(kegg.edges, train.limma, train.pval, kegg.path, train.mcd[MET %in% training.mets, ], 
                                               pval.th = 0.05, lfc.th = 1)
    
    #Store
    master.table<-rbind(master.table, train.table)
  }
  
  #Create train table (with all samples and testing mets, no bootstrap)
  test.limma<-Function.teru.diff.limma(teru.gene.exp, er.plus.samples, norm = F)
  test.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  er.plus.samples, "teru", teru.normal.matrix$MATRIX)
  test.mcd<-Function.met.gene.cor.diff(test.pval, teru.gene.exp, er.plus.samples, kegg.edges, type = "teru")
  test.table<-Function.prep.kegg.pred.table(kegg.edges, test.limma, test.pval, kegg.path, test.mcd[MET %in% testing.mets,], 
                                            pval.th = 0.05, lfc.th = 1)
  
  #Build model
  master.table<-master.table[sample(1:nrow(master.table)),]
  key.z<-sample(letters,1)
  train.h2o<-as.h2o(localH2O, data.frame(master.table), destination_frame = key.z)
  test.h2o<-as.h2o(localH2O, data.frame(test.table), destination_frame = "test")
  FEATURES<-setdiff(colnames(train.h2o), c("MET", "DIFF", "HUGO.MED.LFC"))
  
  #Model
  #   deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame =  train.h2o, use_all_factor_levels = T,
  #                                input_dropout_ratio = 0.0, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
  #                               activation = "RectifierWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)
  
  hidden<-c(20,30,50,80,100,200)
  input.dropout<-c(0,0.1,0.2)
  hidden.dropout<-c(0.1,0.3,0.5)
  for (id in input.dropout){
    for (h in hidden){
      for (hd in hidden.dropout){
        
        deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame = train.h2o, activation = "RectifierWithtDropout", 
                                     epochs = 500, hidden = rep(h, 3),
                                     input_dropout_ratio = id, hidden_dropout_ratios = rep(hd,3))
        #gbm.model<-h2o.gbm(x=FEATURES, y="DIFF", training_frame =  train.h2o, learn_rate = l, ntrees = n, max_depth = m)
        
        #DEEP model predictions 
        gbm.conf.matrix<-h2o.confusionMatrix(deep.model, train.h2o)
        gbm.train.acc<-gbm.conf.matrix$Error[3]
        gbm.train.adj.acc<-var(gbm.conf.matrix$Error[1:2]) + gbm.train.acc
        
        #Predict
        gbm.pred.conf.matrix<-h2o.confusionMatrix(deep.model, test.h2o)
        gbm.test.acc<-gbm.pred.conf.matrix$Error[3]
        gbm.test.adj.acc<-var(gbm.pred.conf.matrix$Error[1:2]) + gbm.test.acc
        
        #Store
        metrics.bootstrap.2<-rbind(metrics.bootstrap.2, data.table(FOLD=i, HIDDEN=h, INPUT.DROPOUT=id, HIDDEN.DROPOUT=hd,
                                                                   GBM.TRAIN.ACC=gbm.train.acc, GBM.TRAIN.ADJ.ACC=gbm.train.adj.acc,
                                                                   GBM.TEST.ACC=gbm.test.acc, GBM.TEST.ADJ.ACC=gbm.test.adj.acc))
        #Save models
        model.id<-paste0("090315.TERUNUMA.PLUS","_",i,"_",h, "_", id, "_",hd,"_", round(gbm.train.acc,3), "_", round(gbm.test.acc,3))
        h2o.saveModel(gbm.model, 
                      paste0("//", h2o.folder, "/"), 
                      model.id)
        
        print (metrics.bootstrap.2)
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z, "test"))) 
      }
    }
  }
}


saveRDS(object = metrics.bootstrap.2 , file = output.file)
print ("Done writing output")
