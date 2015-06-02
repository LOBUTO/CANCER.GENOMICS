library(parallel)
library(data.table)
library(caret)
library(glmnet)

Function.fit.models<-function(MAIN.2HG, BEST.FEATURES){
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "BEST.FEATURES", "MAIN.2HG") ,envir=environment())
  
  #Prep best features
  print ("Calculating best features")
  BEST.AIC<-parSapply(cl, BEST.FEATURES, function(x) extractAIC(lm(METABOLITE~., data.frame(MAIN.2HG)[,c("METABOLITE", x)]))[2])
  BEST.ADJ.SQR<-parSapply(cl, BEST.FEATURES, function(y) summary(lm(METABOLITE~., data.frame(MAIN.2HG)[,c("METABOLITE", y)]))$adj.r.squared)
  BEST.FEATURES.PLOT<-data.table(AIC=BEST.AIC, ADJ.R.SQR=BEST.ADJ.SQR)
  print ("done calculating best features")
  
  #Calculate penalized and unpenalized models
  print ("Calculating best models")
  set.seed(10)
  main.folds<-createFolds(1:nrow(MAIN.2HG), k=5)
  models.lm.2hg<-data.table()
  
  models.lm.2hg<-parLapply(cl, 1:length(BEST.FEATURES), function(aic) {
    
    AIC=as.vector(BEST.FEATURES.PLOT[aic,]$AIC)
    model.2hg<-data.table()
    
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
      model.2hg<-rbind(model.2hg, data.table(AIC=AIC, tr.nrmsd, test.nrmsd, 
                            Lasso.reg.tr.nrmsd, Lasso.reg.test.nrmsd, 
                            Ridge.reg.tr.nrmsd, Ridge.reg.test.nrmsd, 
                            ElNET.reg.tr.nrmsd, ElNET.reg.test.nrmsd))
    }
    #Return
    return(model.2hg)
  }
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Construct all AIC table
  main.table<-do.call(rbind, models.lm.2hg)
  
  #Return
  return(main.table)
}

#Arguments
args<-commandArgs(trailingOnly=T)
main.obj<-args[1]
BEST.FEATURES<-args[2]
output.file<-args[3]

MAIN.2HG<-readRDS(main.obj)
MAIN.FIT.MODELS<-Function.fit.models(MAIN.2HG, BEST.FEATURES)

saveRDS(object = MAIN.FIT.MODELS, file = output.file)
print ("Done writing output")
