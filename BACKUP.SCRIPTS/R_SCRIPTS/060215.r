##MAIN TEST - PREDICT 2HG EXP LEVELS USING COMBINED TANG AND TERU MATRICES
TERU.2HG<-Function.classify.teru.target.met("DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", c("2-hydroxyglutarate"), 
                                            "DATABASES/METABOLOMICS/TERUNUMA.2014/ID.TO.GSM.csv")
TERU.2HG.EXP<-Function.match.teru.expression(TERU.2HG, "DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds", 
                                             "DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", 
                                             "all", normal=F)
TANG.2HG<-Function.classify.tang.target.met("DATABASES/METABOLOMICS/TANG.2014/clean.met.er.csv", c("2-hydroxyglutarate") )
TANG.2HG.EXP<-Function.match.tang.expression(TANG.2HG, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041815.BRCA.MATRICES.AGILENT.rds",
                                             "DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt",
                                             "all", normal=F, exp.version="ma")
TERU.2HG.EXP$ER.STATUS<-factor(ifelse(TERU.2HG.EXP$ER.STATUS=="POS", 1, 2))
TANG.2HG.EXP$ER.STATUS<-factor(ifelse(TANG.2HG.EXP$ER.STATUS=="POS", 1, 2))

#Need to find best features for combined datasets
candidate.features<-setdiff(intersect(colnames(TANG.2HG.EXP), colnames(TERU.2HG.EXP)), c("SAMPLE", "SUBTYPE"))
MAIN.2HG<-rbind(TERU.2HG.EXP[,candidate.features,with=F], TANG.2HG.EXP[,candidate.features,with=F])

#Look for most correlated features
MAIN.2HG.COR<-apply(data.matrix(MAIN.2HG[,!c("ER.STATUS","METABOLITE"), with=F]), 2, function(x) {
  pval<-cor.test(x, MAIN.2HG$METABOLITE, type="spearman")$p.value
  return(pval)
})
MAIN.2HG.COR<-data.table(Hugo_Symbol=colnames(MAIN.2HG[,!c("ER.STATUS","METABOLITE"), with=F]), PVAL=MAIN.2HG.COR)
MAIN.2HG.COR$PVAL.ADJ<-p.adjust(MAIN.2HG.COR$PVAL, method="fdr") #1094 genes are significantly correlated to 2HG
MAIN.2HG.COR<-MAIN.2HG.COR[order(PVAL.ADJ),] #Ordered in matrix

#Send to cluster to look for best features
saveRDS(MAIN.2HG, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/060215.MAIN.2HG.EXP.MET.rds")
BEST.FEATURES<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/O60215.BEST.FEATURES.MAIN.2HG.FILTERED.rds")

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("as.data.table","data.table", "BEST.FEATURES", "MAIN.2HG") ,envir=environment())
BEST.AIC<-parSapply(cl, BEST.FEATURES, function(x) extractAIC(lm(METABOLITE~., data.frame(MAIN.2HG)[,c("METABOLITE", x)]))[2])
BEST.ADJ.SQR<-parSapply(cl, BEST.FEATURES, function(y) summary(lm(METABOLITE~., data.frame(MAIN.2HG)[,c("METABOLITE", y)]))$adj.r.squared)
stopCluster(cl)
BEST.FEATURES.PLOT<-data.table(AIC=BEST.AIC, ADJ.R.SQR=BEST.ADJ.SQR)

ggplot(BEST.FEATURES.PLOT, aes(AIC, ADJ.R.SQR)) + geom_point() + theme.format
order(BEST.FEATURES.PLOT$AIC)[1:5]

#Now do regularization with lm based on best models
set.seed(10)
main.folds<-createFolds(1:nrow(MAIN.2HG), k=5)

models.lm.2hg<-data.table()
BEST.FEATURES.PLOT[489,]
for (aic in 1:length(BEST.FEATURES) ){
  
  AIC=as.vector(BEST.FEATURES.PLOT[aic,]$AIC)
  
  for (i in main.folds){
    
    it.folds<-i
    other.folds<-setdiff(unlist(main.folds), it.folds)
    feat.limit<-length(other.folds)-2
    
    training.matrix<-data.matrix(data.frame(MAIN.2HG)[other.folds, BEST.FEATURES[[aic]][1:feat.limit]])
    training.target<-data.matrix(MAIN.2HG[other.folds,"METABOLITE", with=F])
    testing.matrix<-data.matrix(data.frame(MAIN.2HG)[it.folds, BEST.FEATURES[[aic]][1:feat.limit]])
    testing.target<-data.matrix(MAIN.2HG[it.folds,"METABOLITE",with=F])
    tr.met.max<-max(as.vector(training.target))
    tr.met.min<-min(as.vector(training.target))
    test.met.max<-max(as.vector(testing.target))
    test.met.min<-min(as.vector(testing.target))
    
    #Non-regularized
    nr.model<-lm(METABOLITE~., data.frame(cbind(training.matrix,training.target)))
    tr.predict<-predict.lm(nr.model, data.frame(training.matrix))
    test.predict<-predict.lm(nr.model, data.frame(testing.matrix))
    tr.nrmsd<-sqrt(mean( (as.vector(training.target)-tr.predict)^2)) / (tr.met.max - tr.met.min)
    test.nrmsd<-sqrt(mean( (as.vector(testing.target)-test.predict)^2)) / (test.met.max - test.met.min)
    
    #Regularized with Lasso, Ridge and Elastic net
    penalty<-c(1, 0, 0.5)
    names(penalty)<-c("Lasso", "Ridge", "ElNET")
    for (p in names(penalty)){
      
      #Build penalty model and apply to train and test
      reg.fit<-cv.glmnet(x=training.matrix, y=training.target, nfolds=5, alpha=penalty[[p]]) #to find best lambda
      reg.tr.pred<-predict(reg.fit, newx=training.matrix, s="lambda.min")
      reg.test.pred<-predict(reg.fit, newx=testing.matrix, s="lambda.min")
      reg.tr.nrmsd<-sqrt(mean( (as.vector(training.target)-reg.tr.pred)^2)) / (tr.met.max - tr.met.min)
      reg.test.nrmsd<-sqrt(mean( (as.vector(testing.target)-reg.test.pred)^2)) / (test.met.max - test.met.min)
      assign(paste(p, "reg.tr.nrmsd", sep="."), reg.tr.nrmsd)
      assign(paste(p, "reg.test.nrmsd", sep="."), reg.test.nrmsd)
    }
    
    #Return MSE
    models.lm.2hg<-rbind(models.lm.2hg, data.table(AIC=AIC, tr.nrmsd, test.nrmsd, 
                                                   Lasso.reg.tr.nrmsd, Lasso.reg.test.nrmsd, 
                                                   Ridge.reg.tr.nrmsd, Ridge.reg.test.nrmsd, 
                                                   ElNET.reg.tr.nrmsd, ElNET.reg.test.nrmsd))  
  }
}
saveRDS(models.lm.2hg, "060215.MAIN.2HG.BEST.FEATURES.AIC.MODELS.rds")
models.lm.2hg.melt<-melt(models.lm.2hg, id.vars = "AIC")
models.lm.2hg.melt$SET<-rep(rep(c("TRAIN", "TEST"), 4), each=5000)
models.lm.2hg.melt$REGULARIZATION<-rep(c("None", "Lasso", "Ridge", "Elastic.Net"), each = 10000)
ggplot(models.lm.2hg.melt[,list(MEAN.NRMSE=mean(value)), by=c("AIC", "SET", "REGULARIZATION")], 
       aes(AIC, MEAN.NRMSE, colour=REGULARIZATION)) + geom_line() + facet_wrap(~SET) + theme.format
ggplot(models.lm.2hg.melt[,list(MEAN.NRMSE=mean(value)), by=c("AIC", "SET", "REGULARIZATION")][REGULARIZATION!="None",], 
       aes(AIC, MEAN.NRMSE, colour=REGULARIZATION)) + geom_line() + facet_wrap(~SET,ncol = 1, scales = "free") + theme.format

#Test best model!
models.lm.2hg.melt[,list(MEAN.NRMSE=mean(value)), by=c("AIC", "SET", "REGULARIZATION")][SET=="TEST",][order(MEAN.NRMSE),]$AIC[1] # @AIC=-541.1382
FEAT.N<-which(BEST.FEATURES.PLOT$AIC==models.lm.2hg.melt[,list(MEAN.NRMSE=mean(value)), by=c("AIC", "SET", "REGULARIZATION")][SET=="TEST",][order(MEAN.NRMSE),]$AIC[1])
BEST.FEATURES[[FEAT.N]]
x<-copy(MAIN.2HG)

set.seed(10)
main.folds<-createFolds(1:nrow(MAIN.2HG), k=5)
best.models.lm.2hg<-data.table()

for (i in main.folds){
  
  it.folds<-i
  other.folds<-setdiff(unlist(main.folds), it.folds)
  feat.limit<-length(other.folds)-2
  
  training.matrix<-data.matrix(data.frame(MAIN.2HG)[other.folds, BEST.FEATURES[[FEAT.N]][1:feat.limit]])
  training.target<-data.matrix(MAIN.2HG[other.folds,"METABOLITE", with=F])
  testing.matrix<-data.matrix(data.frame(MAIN.2HG)[it.folds, BEST.FEATURES[[FEAT.N]][1:feat.limit]])
  testing.target<-data.matrix(MAIN.2HG[it.folds,"METABOLITE",with=F])
  tr.met.max<-max(as.vector(training.target))
  tr.met.min<-min(as.vector(training.target))
  test.met.max<-max(as.vector(testing.target))
  test.met.min<-min(as.vector(testing.target))
  
  #Regularized with Lasso  
  #Build penalty model and apply to train and test
  reg.fit<-cv.glmnet(x=training.matrix, y=training.target, nfolds=5, alpha=1) #LASSO penalty
  reg.tr.pred<-predict(reg.fit, newx=training.matrix, s="lambda.min")
  reg.test.pred<-predict(reg.fit, newx=testing.matrix, s="lambda.min")
  Lasso.reg.tr.nrmsd<-sqrt(mean( (as.vector(training.target)-reg.tr.pred)^2)) / (tr.met.max - tr.met.min)
  Lasso.reg.test.nrmsd<-sqrt(mean( (as.vector(testing.target)-reg.test.pred)^2)) / (test.met.max - test.met.min)
  
  #Return MSE
  best.models.lm.2hg<-rbind(best.models.lm.2hg, data.table(Lasso.reg.tr.nrmsd, Lasso.reg.test.nrmsd))  
}
best.models.lm.2hg

x$PREDICT<-predict(reg.fit, newx=data.matrix(data.frame(MAIN.2HG)[,BEST.FEATURES[[FEAT.N]][1:feat.limit]]), s="lambda.min") #Chosen from best model in iteartion at best.models.lm.2hg
ggplot(x, aes(METABOLITE, METABOLITE)) + geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="predict")) + theme.format


#Test with deep learning (build and predict without er.status in model)
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
key.z<-sample(letters,1)
MAIN.2HG<-MAIN.2HG[sample(nrow(MAIN.2HG)),]
h2o_2HG<-as.h2o(localH2O, data.frame(MAIN.2HG), key=key.z) 

DEEP.2HG<-data.table()
METHOD="Tanh"
HIDDEN<-c(80,40)
INPUT.DR<-c(0)
HIDDEN.DR<-rep(0.5, length(HIDDEN))
FEATURES<-80

for (id in INPUT.DR){
  for (n in 1:10) {
    #Model with/out dropout
    print (c("building model", id, n))
    
    MODEL.2HG<-h2o.deeplearning(x=BEST.FEATURES[[FEAT.N]], y="METABOLITE", data=h2o_2HG, classification = F, nfolds = 5,
                                activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500)
                               #input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )
    
    #     MODEL.2HG@xval
    #     MODEL.2HG@model
    #     MODEL.2HG@model$train_class_error
    #     MODEL.2HG@model$valid_class_error
    
    CUR.PRED<-data.table(TR.SQR.ERROR=MODEL.2HG@model$train_sqr_error, VALID.SQR.ERROR=MODEL.2HG@model$valid_sqr_error, FEATURES=FEATURES,
                         ITER=n, METHOD=METHOD, INPUT.DR=id, HIDDEN=paste(HIDDEN,collapse="."), HIDDEN.DR=paste(HIDDEN.DR,collapse="."))
    
    #Assign predictors
    DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
    
    #Clean H2o memory
    h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.z))
  }  
}
DEEP.2HG

ggplot(DEEP.2HG, aes(HIDDEN, sqrt(TR.SQR.ERROR)/(max(MAIN.2HG$METABOLITE)-min(MAIN.2HG$METABOLITE)), colour=factor(INPUT.DR))) + geom_boxplot() + facet_wrap(~METHOD)
ggplot(DEEP.2HG, aes(HIDDEN, sqrt(VALID.SQR.ERROR)/(max(MAIN.2HG$METABOLITE)-min(MAIN.2HG$METABOLITE)), colour=factor(INPUT.DR))) + geom_boxplot() + facet_wrap(~METHOD)

MODEL.2HG<-h2o.deeplearning(x=BEST.FEATURES[[FEAT.N]], y="METABOLITE", data=h2o_2HG, classification = F, nfolds = 5,
                            activation = "TanhWithDropout", balance_classes = TRUE, hidden = c(80,40), epochs = 500,
                            input_dropout_ratio = 0.2 , hidden_dropout_ratios =c(0.5,0.5) )

x<-copy(MAIN.2HG)
x$PREDICT<-as.numeric(as.matrix(h2o.predict(MODEL.2HG, h2o_2HG)))
ggplot(x, aes(METABOLITE, METABOLITE)) + geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="predict"))

hist(replicate(100,wilcox.test(runif(100), runif(100), paired=T)$p.value))
hist(replicate(100,wilcox.test(runif(1000), runif(1000), paired=T)$p.value))
hist(replicate(100,wilcox.test(runif(10), runif(10), paired=T)$p.value))

