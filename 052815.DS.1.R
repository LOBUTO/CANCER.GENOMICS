
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
METHOD="TanhWithDropout"
HIDDEN<-c(1600,800,400)
INPUT.DR<-seq(0.1,0.3,0.1)
HIDDEN.DR<-rep(0.5, length(HIDDEN))

for (id in INPUT.DR){
  for (n in 1:10) {
    #Model with/out dropout
    print (c("building model", id, n))
    
    MODEL.2HG<-h2o.deeplearning(x=2:ncol(TERU.2HG.EXP), y=1, data=h2o_2HG, classification = F, nfolds = 5,
                                activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 300,
                                input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )
    
    #     MODEL.2HG@xval
    #     MODEL.2HG@model
    #     MODEL.2HG@model$train_class_error
    #     MODEL.2HG@model$valid_class_error
    
    CUR.PRED<-data.table(TR.SQR.ERROR=MODEL.2HG@model$train_sqr_error, VALID.SQR.ERROR=MODEL.2HG@model$valid_sqr_error,
                          ITER=n, METHOD=METHOD, INPUT.DR=id, HIDDEN=paste(HIDDEN,collapse="."), HIDDEN.DR=paste(HIDDEN.DR,collapse="."))
    
    #Assign predictors
    DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
    
    #Clean H2o memory
    h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))
  }  
}

saveRDS(object = DEEP.2HG , file = output.file)