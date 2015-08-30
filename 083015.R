#BOOTSTRAP.2
#Boostrapping samples for 9/10 of metabolites across all individuals and testing on the rest of metabolites

#The limit will be at teru.minus.mcd
Function.met.split<-function(limit.met, seed, pval.met, pval.th = 0.05, lfc.th = 1, folds=10){
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

minus.met.splits<-Function.met.split(teru.minus.mcd$MET, 23, teru.minus.pval.met, folds = 10)

metrics.bootstrap.2<-data.table()


for (i in 1:length(minus.met.splits)){
  
  #Split mets
  testing.mets<-minus.met.splits[[i]]
  training.mets<-setdiff(unlist(minus.met.splits), testing.mets)

  #Obtain bootstrapped samples
  er.minus.samples<-teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE
  bootstrap.train<-replicate(10, sample(er.minus.samples, 10), simplify = F)
  
  #Create master table per bootsrapped set (for bootstrapped samples and training mets)
  master.table<-data.table()
  
  for (j in bootstrap.train){
    
    train.limma<-Function.teru.diff.limma(teru.gene.exp, j, norm = F)
    train.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  j , "teru", teru.normal.matrix$MATRIX)
    train.mcd<-Function.met.gene.cor.diff(train.pval, teru.gene.exp, j, kegg.edges, type = "teru")
    train.table<-Function.prep.kegg.pred.table(kegg.edges, train.limma, train.pval, kegg.path, train.mcd[MET %in% training.mets, ], pval.th = 0.05, lfc.th = 1)
    
    #Store
    master.table<-rbind(master.table, train.table)
  }
  
  #Create train table (with all samples and testing mets, no bootstrap)
  test.limma<-Function.teru.diff.limma(teru.gene.exp, er.minus.samples, norm = F)
  test.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  er.minus.samples, "teru", teru.normal.matrix$MATRIX)
  test.mcd<-Function.met.gene.cor.diff(test.pval, teru.gene.exp, er.minus.samples, kegg.edges, type = "teru")
  test.table<-Function.prep.kegg.pred.table(kegg.edges, test.limma, test.pval, kegg.path, test.mcd[MET %in% testing.mets,], pval.th = 0.05, lfc.th = 1)
  
  #Build model
  master.table<-master.table[sample(1:nrow(master)),]
  key.z<-sample(letters,1)
  train.h2o<-as.h2o(localH2O, data.frame(master.table), destination_frame = key.z)
  test.h2o<-as.h2o(localH2O, data.frame(test.table), destination_frame = "test")
  FEATURES<-setdiff(colnames(train.h2o), c("MET", "DIFF"))
  
  #Model
  deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame =  train.h2o, use_all_factor_levels = T,
                               input_dropout_ratio = 0.0, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
                              activation = "RectifierWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)
  
  gbm.model<-h2o.gbm(x=FEATURES, y="DIFF", training_frame =  train.h2o, learn_rate = 0.05, ntrees = 1000)
  
  
  #DEEP model predictions
  #deep.f1<-deep.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[1] 
  deep.conf.matrix<-h2o.confusionMatrix(deep.model, train.h2o )
  deep.train.acc<-deep.conf.matrix$Error[3]
  deep.train.adj.acc<-var(c(deep.conf.matrix$Error[1],deep.conf.matrix$Error[2])) + deep.train.acc
  
  #Predict
  deep.pred.conf.matrix<-h2o.confusionMatrix(deep.model, test.h2o)
  deep.test.acc<-deep.pred.conf.matrix$Error[3]
  deep.test.adj.acc<-var(c(deep.pred.conf.matrix$Error[1], deep.pred.conf.matrix$Error[2])) + deep.test.acc
  
  #GBM model predictions 
  gbm.conf.matrix<-h2o.confusionMatrix(gbm.model, train.h2o )
  gbm.train.acc<-gbm.conf.matrix$Error[3]
  gbm.train.adj.acc<-var(c(gbm.conf.matrix$Error[1], gbm.conf.matrix$Error[2])) + gbm.train.acc
  
  #Predict
  gbm.pred.conf.matrix<-h2o.confusionMatrix(gbm.model, test.h2o)
  gbm.test.acc<-gbm.pred.conf.matrix$Error[3]
  gbm.test.adj.acc<-var(c(gbm.pred.conf.matrix$Error[1], gbm.pred.conf.matrix$Error[2])) + gbm.test.acc
  
  #Store
  metrics.bootstrap.2<-rbind(metrics.bootstrap.2, data.table(FOLD=i,
                                                             TRAIN.ACC=deep.train.acc, TRAIN.ADJ.ACC=deep.train.adj.acc,
                                                             TEST.ACC=deep.test.acc, TEST.ADJ.ACC=deep.test.adj.acc,
                                                             GBM.TRAIN.ACC=gbm.train.acc, GBM.TRAIN.ADJ.ACC=gbm.train.adj.acc,
                                                             GBM.TEST.ACC=gbm.test.acc, GBM.TEST.ADJ.ACC=gbm.test.adj.acc))
  print (metrics.bootstrap.2)
  h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z))) 
  
}
