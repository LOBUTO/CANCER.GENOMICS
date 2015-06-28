library(parallel)
library(data.table)
library(h2o)

#Arguments
args<-commandArgs(trailingOnly=T)
digit.train<-fread(args[1], header=T, sep=",", stringsAsFactors = F)
output.file<-args[2]
print ("done loading files")

######Execute######
digit.train<-digit.train[sample(1:nrow(digit.train)),] #randomize
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
key.z<-sample(letters,1)
train.h2o<-as.h2o(localH2O, data.frame(digit.train), key=key.z)

#MODEL
DEEP.2HG<-data.table()

#No dropout model
print ("modeling without dropout")
MODEL.2HG<-h2o.deeplearning(x=2:785, y=1, data=train.h2o, classification = T, nfolds = 10,
                            activation = "Tanh", balance_classes = TRUE, hidden = c(800,800,800), epochs = 500)

DEEP.2HG<-rbind(DEEP.2HG, data.table(METHOD="Tanh", INPUT.DR=0.5, HIDDEN="800.800.800", HIDDEN.DR=0,
                                     TRAIN.CLASS.ERROR=MODEL.2HG@model$train_class_error,
                                     VALID.CLASS.ERROR=MODEL.2HG@model$valid_class_error))

#Dropout model
print ("modeling with dropout")
METHOD="TanhWithDropout"
HIDDEN<-c(800, 800, 800)
INPUT.DR<-c(0, 0.1, 0.2)
HIDDEN.DR<-rep(0.5, length(HIDDEN))

for (id in INPUT.DR){
  
  #Model with/out dropout
  print (c("building model", id))
  
  MODEL.2HG<-h2o.deeplearning(x=2:785, y=1, data=train.h2o, classification = T, nfolds = 10,
                              activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                              input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR)
  
  CUR.PRED<-data.table(METHOD=METHOD, INPUT.DR=id,  HIDDEN=paste(HIDDEN, collapse="."), HIDDEN.DR=paste(HIDDEN.DR,collapse="."),
                       TRAIN.CLASS.ERROR=MODEL.2HG@model$train_class_error,
                       VALID.CLASS.ERROR=MODEL.2HG@model$valid_class_error))
  
  #Assign predictors
  DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
  
  #Clean H2o memory
  h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.z))

}
###################

#Write to output
saveRDS(object = DEEP.2HG, file = output.file)
print ("Done writing output")