
#Arguments
args<-commandArgs(trailingOnly=T)
met.obj<-args[1]
output.file<-args[2]

library(h2o)
library(data.table)

TERU.2HG.EXP<-readRDS(met.obj)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 

#Convert data into h2o
key.v<-sample(letters,1)
TERU.2HG.EXP<-TERU.2HG.EXP[sample(nrow(TERU.2HG.EXP)),]
h2o_2HG<-as.h2o(localH2O, TERU.2HG.EXP, key=key.v) 

DEEP.2HG<-data.table()

FEATURES<-c(10,20,30,40,50,100,200,400)
METHODS=c("Tanh", "TanhWithDropout")
INPUT.DR<-c(0.1,0.2)

for (feat in FEATURES){
  HIDDEN<-c(feat, round(feat/2), round(feat/10))
  HIDDEN.DR<-rep(0.5, length(HIDDEN))
  for (method in METHODS){
    for (id in INPUT.DR){
      for (n in 1:10) {
        #Model with/out dropout
        print (c("building model", method, feat, id, n))
        
        if (method=="TanhWithDropout"){
          MODEL.2HG<-h2o.deeplearning(x=2:feat, y=1, data=h2o_2HG, classification = F, nfolds = 5,
                                      activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 300,
                                      input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )  
        } else {
          MODEL.2HG<-h2o.deeplearning(x=2:feat, y=1, data=h2o_2HG, classification = F, nfolds = 5,
                                      activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 300)
        }
        
        
        CUR.PRED<-data.table(TR.SQR.ERROR=MODEL.2HG@model$train_sqr_error, VALID.SQR.ERROR=MODEL.2HG@model$valid_sqr_error,
                             ITER=n, METHOD=method, INPUT.DR=id, HIDDEN=paste(HIDDEN,collapse="."), HIDDEN.DR=paste(HIDDEN.DR,collapse="."))
        
        #Assign predictors
        DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
        
        #Clean H2o memory
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))
      }  
    }
  }
}

saveRDS(object = DEEP.2HG , file = output.file)