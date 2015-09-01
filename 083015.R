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
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 

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
  master.table<-master.table[sample(1:nrow(master.table)),]
  key.z<-sample(letters,1)
  train.h2o<-as.h2o(localH2O, data.frame(master.table), destination_frame = key.z)
  test.h2o<-as.h2o(localH2O, data.frame(test.table), destination_frame = "test")
  FEATURES<-setdiff(colnames(train.h2o), c("MET", "DIFF"))
  
  #Model
#   deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame =  train.h2o, use_all_factor_levels = T,
#                                input_dropout_ratio = 0.0, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
#                               activation = "RectifierWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)
  
  learn.rate<-c(0.005,0.01,0.05,0.1,0.2)
  max.depth<-c(2,4,6,8,10,20)
  n.trees<-c(5,10,15,20,50,100,200)
  for (l in learn.rate){
    for (m in max.depth){
      for (n in n.trees){
        
        gbm.model<-h2o.gbm(x=FEATURES, y="DIFF", training_frame =  train.h2o, learn_rate = l, ntrees = n, max_depth = m)
        
        #GBM model predictions 
        gbm.conf.matrix<-h2o.confusionMatrix(gbm.model, train.h2o)
        gbm.train.acc<-gbm.conf.matrix$Error[3]
        gbm.train.adj.acc<-var(gbm.conf.matrix$Error[1:2]) + gbm.train.acc
        
        #Predict
        gbm.pred.conf.matrix<-h2o.confusionMatrix(gbm.model, test.h2o)
        gbm.test.acc<-gbm.pred.conf.matrix$Error[3]
        gbm.test.adj.acc<-var(gbm.pred.conf.matrix$Error[1:2]) + gbm.test.acc
        
        #Store
        metrics.bootstrap.2<-rbind(metrics.bootstrap.2, data.table(FOLD=i, LEARN.RATE=l, MAX.DEPTH=m, N.TREES=n,
                                                                   GBM.TRAIN.ACC=gbm.train.acc, GBM.TRAIN.ADJ.ACC=gbm.train.adj.acc,
                                                                   GBM.TEST.ACC=gbm.test.acc, GBM.TEST.ADJ.ACC=gbm.test.adj.acc))
        #Save models
        model.id<-paste0("083015.TERUNUMA.MINUS","_",i,"_",l, "_", m, "_",n,"_", round(gbm.train.acc,3), "_", round(gbm.test.acc,3))
        h2o.saveModel(gbm.model, 
                      paste0("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS.2/GBM/BOOTSTRAP/"), 
                      model.id)
        
        print (metrics.bootstrap.2)
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z, "test"))) 
      }
    }
  }
}

metrics.bootstrap.2
ggplot(melt(metrics.bootstrap.2[,c(2:4,5,7),with=F], measure.vars = c("GBM.TRAIN.ACC", "GBM.TEST.ACC")), aes(factor(LEARN.RATE), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(MAX.DEPTH~N.TREES)

metrics.bootstrap.2[order(GBM.TEST.ACC),]

best.minus.model<-h2o.loadModel("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS.2/GBM/BOOTSTRAP/083015.TERUNUMA.MINUS_8_0.1_4_200_0.034_0.062",localH2O)
best.minus.model@parameters$training_frame<-"holder.4"
h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.549989)
