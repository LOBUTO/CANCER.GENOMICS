#DEEP.ITER.2

library(h2o)
library(data.table)
library(reshape2)
library(caret)

#Arguments
args<-commandArgs(trailingOnly=T)
teru.met.feat.std<-readRDS(args[1])
h2o.folder<-args[2]

output.file<-args[3]

#Initialize h2o
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 

teru.met.h2o<-as.h2o(localH2O, data.frame(teru.met.feat.std), "teru.met")
FEATURES<-setdiff(colnames(teru.met.h2o), c("SAMPLE", "NAME", "MET.LEVELS"))

n_folds<-3
rand_folds<-createFolds(as.factor(as.matrix(teru.met.h2o$MET.LEVELS)), k=n_folds)


#Iterate through sets and build models
metrics.bootstrap<-data.table()
for (fold in 1:n_folds){
  
  test_rows<-as.numeric(unlist(rand_folds[fold]))  
  train_rows<-setdiff(as.numeric(unlist(rand_folds)), test_rows)
  
  iterations<-1:10
  hidden<-c(50,100,200)
  input.dropout<-c(0,0.1,0.2)
  hidden.dropout<-c(0.2, 0.5)
  
  for (iter in iterations){
    for (h in hidden){
      for (id in input.dropout){
        for (hd in hidden.dropout){
          
          deep.model<-h2o.deeplearning(x = FEATURES, y="MET.LEVELS", training_frame = teru.met.h2o[train_rows,], 
                                       activation = "TanhWithDropout", hidden = rep(h, 3), epochs = 500, 
                                       input_dropout_ratio = id, hidden_dropout_ratios = rep(hd,3))
          
          #DEEP model predictions 
          conf.matrix<-h2o.confusionMatrix(deep.model, teru.met.h2o[train_rows,])
          train.acc<-conf.matrix$Error[3]
          
          #Predict
          conf.matrix<-h2o.confusionMatrix(deep.model, teru.met.h2o[test_rows,])
          test.acc<-conf.matrix$Error[3]
          
          #Store
          metrics.bootstrap<-rbind(metrics.bootstrap, data.table(FOLD=fold, ITERATION=iter, HIDDEN=h, INPUT.DROPOUT=id, HIDDEN.DROPOUT=hd,
                                                                 TRAIN.ACC=train.acc, TEST.ACC=test.acc))
          #Save models
          model.id<-paste0("090715.TERUNUMA.PLUS","_",fold,"_",iter, "_",h, "_", id, "_",hd,"_", round(train.acc,3), "_", round(test.acc,3))
          h2o.saveModel(deep.model, 
                        paste0("//", h2o.folder, "/"), 
                        model.id)
          
          print (metrics.bootstrap)
          h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c("teru.met"))) 
        }
      }
    }  
  }
}

#Store output
saveRDS(object = metrics.bootstrap , file = output.file)
print ("Done writing output")