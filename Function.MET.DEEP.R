#Function.MET.DEEP.R
#050915
#Calculate prediciton errors of metabolite levels using h2o's deep learning module

library(data.table)
library(reshape2)
library(h2o)
library(caret)

#Functions
Function.Main<-function(met.obj, method, hidden, hidden.dr, input.dr){
  
  #Load met table
  MET<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  h2o_MET<-as.h2o(localH2O, MET, key=key.v) 
  
  #Introduce parameters for training
  DEEP.MET<-data.table()
  METHOD=method
  HIDDEN<-as.numeric(unlist(strsplit(hidden,"[.]")))
  INPUT.DR<-as.numeric(input.dr)
  HIDDEN.DR<-1/as.numeric(unlist(strsplit(hidden.dr,"[.]")))
  FEATURES<-"ALL.RECON.EXP"
  
  for (n in 1:10) {
    
    #Model with/out dropout
    print (c("building model:",n))
    
    if (METHOD=="TanhWithDropout"){
      MODEL.MET<-h2o.deeplearning(x=5:ncol(MET), y=2, data=h2o_MET, classification = F, nfolds = 5,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                  input_dropout_ratio = INPUT.DR , hidden_dropout_ratios = HIDDEN.DR)  
    } else if (METHOD=="Tanh"){
      MODEL.MET<-h2o.deeplearning(x=5:ncol(MET), y=2, data=h2o_MET, classification = F, nfolds = 5,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500)
    }
    
    
    
    TRAIN.ERROR<-MODEL.MET@model$train_sqr_error
    TEST.ERROR<-MODEL.MET@model$valid_sqr_error
    
    CUR.PRED<-data.table(ITER=n, TRAIN.ACC=1-sqrt(TRAIN.ERROR), TEST.ACC=1-sqrt(TEST.ERROR), METHOD=METHOD, FEATURES=FEATURES,
                         HIDDEN=paste(HIDDEN,collapse = "."), INPUT.DR=INPUT.DR, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
    
    #Assign predictors
    DEEP.MET<-rbind(DEEP.MET,CUR.PRED)
  }
  
  #Return
  return(DEEP.MET)
}

Function.Main.Bin<-function(met.obj, method, hidden, hidden.dr, input.dr){
  
  
  #Load met table
  MET<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  h2o_MET<-as.h2o(localH2O, MET, key=key.v) 
  
  #Introduce parameters for training
  DEEP.MET<-data.table()
  METHOD=method
  HIDDEN<-as.numeric(unlist(strsplit(hidden,"[.]")))
  INPUT.DR<-as.numeric(input.dr)
  HIDDEN.DR<-1/as.numeric(unlist(strsplit(hidden.dr,"[.]")))
  FEATURES<-"ALL.RECON.EXP.BIN"
  
  #Break into equal representative folds of 1/0 dysregulated cohorts
  x<-ifelse(MET$MET>2, 1, 0)
  RAND_FOLDS.1<-createFolds(which(x==1),5)
  RAND_FOLDS.0<-createFolds(which(x==0),5)
  TEST_ROWS<-c( which(x==1)[RAND_FOLDS.1[[1]]], which(x==0)[RAND_FOLDS.0[[1]]])
  TRAIN_ROWS<-setdiff(1:nrow(MET), TEST_ROWS)
  
  for (n in 1:10) {
    
    #Model with/out dropout
    print (c("building model",n))
    
    if (METHOD=="TanhWithDropout"){
      MODEL.MET<-h2o.deeplearning(x=5:200, y=2, data=h2o_MET[TRAIN_ROWS,], classification = F, nfolds = 4,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                  input_dropout_ratio = INPUT.DR , hidden_dropout_ratios = HIDDEN.DR)  
    } else {
      MODEL.MET<-h2o.deeplearning(x=5:200, y=2, data=h2o_MET[TRAIN_ROWS,], classification = F, nfolds = 4,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500)  
    }
    
    TRAIN.ERROR<-MODEL.MET@model$train_sqr_error
    TEST.ERROR<-MODEL.MET@model$valid_sqr_error
    
    TRAIN.ACC<-1-sqrt(TRAIN.ERROR)/(max(MET$MET)- min(MET$MET))
    TEST.ACC<-1-sqrt(TEST.ERROR)/(max(MET$MET)- min(MET$MET))
    
    TEST.SAMPLE<-MET[TEST_ROWS, c("SAMPLE", "MET"), with=F]
    TEST.SAMPLE$BIN.MET<-ifelse(TEST.SAMPLE$MET>2, 1, 0)
    TEST.SAMPLE$PREDICT<-as.numeric(as.matrix(h2o.predict(MODEL.MET, h2o_MET[TEST_ROWS,])$predict))
    TEST.SAMPLE$BIN.PRED<-ifelse(TEST.SAMPLE$PREDICT>2, 1, 0)
    VALID.ACC<-mean(TEST.SAMPLE$BIN.MET==TEST.SAMPLE$BIN.PRED)
    
    CUR.PRED<-data.table(TRAIN.ACC=TRAIN.ACC, TEST.ACC=TEST.ACC, VALID.ACC=VALID.ACC, METHOD=METHOD, FEATURES=FEATURES,
                         HIDDEN=paste(HIDDEN,collapse="."), INPUT.DR=INPUT.DR, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
    
    #Assign predictors
    DEEP.MET<-rbind(DEEP.MET, CUR.PRED)
  }
  
  #Return
  return(DEEP.MET)
}

Function.Main.Class<-function(met.obj, method, hidden, hidden.dr, input.dr){
  
  #Load met table
  MET<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  h2o_MET<-as.h2o(localH2O, MET, key=key.v) 
  
  #Introduce parameters for training
  DEEP.MET<-data.table()
  METHOD=method
  HIDDEN<-as.numeric(unlist(strsplit(hidden,"[.]")))
  INPUT.DR<-as.numeric(input.dr)
  HIDDEN.DR<-1/as.numeric(unlist(strsplit(hidden.dr,"[.]")))
  FEATURES<-c(7,10,20,50,100,200,400,800,ncol(MET))
  
  #Distribute into folds evenly for cross-validation
  RAND_FOLDS.1<-createFolds(which(MET$MET==1),5)
  RAND_FOLDS.0<-createFolds(which(MET$MET==0),5)
  
  for (f in FEATURES){
    
    for (r in 1:5){
      
      TEST_ROWS<-c( which(MET$MET==1)[RAND_FOLDS.1[[r]]], which(MET$MET==0)[RAND_FOLDS.0[[r]]])
      TRAIN_ROWS<-setdiff(1:nrow(MET), TEST_ROWS)
      
      for (n in 1:10) {
        
        print (c("building model", f, r, n))
        
        #Model with/out dropout
        if (METHOD=="TanhWithDropout"){
          MODEL.MET<-h2o.deeplearning(x=5:f, y=2, data=h2o_MET[TRAIN_ROWS,], classification = T,
                                      activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                      input_dropout_ratio = INPUT.DR , hidden_dropout_ratios = HIDDEN.DR)  
        } else {
          MODEL.MET<-h2o.deeplearning(x=5:f, y=2, data=h2o_MET[TRAIN_ROWS,], classification = T, 
                                      activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500)  
        }
        
        PREDICT.MET.TRAIN<-h2o.predict(MODEL.MET, h2o_MET[TRAIN_ROWS,])$predict
        PREDICT.MET.TRAIN<-as.factor(as.matrix(PREDICT.MET.TRAIN))
        PREDICT.MET.TEST<-h2o.predict(MODEL.MET, h2o_MET[TEST_ROWS,])$predict
        PREDICT.MET.TEST<-as.factor(as.matrix(PREDICT.MET.TEST))
        
        TRAIN.ACC<-round(confusionMatrix(PREDICT.MET.TRAIN, MET$MET[TRAIN_ROWS])$overall[1],4)
        TEST.ACC<-round(confusionMatrix(PREDICT.MET.TEST, MET$MET[TEST_ROWS])$overall[1],4)  
        
        CUR.PRED<-data.table(TRAIN.ACC=TRAIN.ACC, TEST.ACC=TEST.ACC, ITER=n,  FOLD=r, FEATURES=f-5+1,  METHOD=METHOD,
                             HIDDEN=paste(HIDDEN,collapse="."), INPUT.DR=INPUT.DR, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
        
        #Assign predictors
        DEEP.MET<-rbind(DEEP.MET,CUR.PRED)
        
        #Clean H2o memory
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))
      }
    }
  }
  
  #Return
  return(DEEP.MET)
}

Function.Main.Class.2<-function(met.obj, method, hidden, hidden.dr, input.dr){
  
  #Load met table
  MET.TABLE<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  MET.TABLE<-MET.TABLE[sample(nrow(MET.TABLE)),] #Too balanace classes!!!!
  h2o_MET<-as.h2o(localH2O, MET.TABLE, key=key.v) 
  
  #Introduce parameters for training
  DEEP.MET<-data.table()
  METHOD=method
  HIDDEN<-as.numeric(unlist(strsplit(hidden,"[.]")))
  INPUT.DR<-as.numeric(input.dr)
  HIDDEN.DR<-1/as.numeric(unlist(strsplit(hidden.dr,"[.]")))
  FEATURES<-c(7,10,20,50,100,200,400,800,ncol(MET.TABLE))
  
  for (f in FEATURES){
    for (n in 1:10) {
      
      #Model with/out dropout
      print (c("building model", f, n))
      
      #Model with/out dropout
      if (METHOD=="TanhWithDropout"){
        MODEL.MET<-h2o.deeplearning(x=5:f, y=2, data=h2o_MET, classification = T, nfolds=5,
                                    activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                    input_dropout_ratio = INPUT.DR , hidden_dropout_ratios = HIDDEN.DR)  
      } else if (METHOD=="Tanh"){
        MODEL.MET<-h2o.deeplearning(x=5:f, y=2, data=h2o_MET, classification = T, nfolds=5, 
                                    activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500)  
      }
      
      CUR.PRED<-data.table(TRAIN.ACC=1-MODEL.MET@model$train_class_error, TEST.ACC=1-MODEL.MET@model$valid_class_error, ITER=n, METHOD=METHOD, FEATURES=f,
                           HIDDEN=paste(HIDDEN,collapse="."), INPUT.DR=INPUT.DR, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
      
      #Assign predictors
      DEEP.MET<-rbind(DEEP.MET, CUR.PRED)
      
      #Clean H2o memory
      h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))
    } 
  }
  
  #Return
  return(DEEP.MET)
}

Function.MET.RANDOM.DEEP<-function(met.obj) {
  
  ##########SET UP###########
  #Load met table
  MET.TABLE<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  MET.TABLE<-MET.TABLE[sample(nrow(MET.TABLE)),] #Too balanace classes!!!!
  h2o_MET<-as.h2o(localH2O, MET.TABLE, key=key.v) 
  ###########################
  
  ########EXECUTE###########
  FEAT.COLUMNS<-colnames(MET.TABLE)[5:ncol(MET.TABLE)]
  TARGET.COLUMN<-"MET"
  
  #Randomize feature columns 
  seed<-sample(1:10000,1)
  set.seed(seed)
  ran.feat<-sample(FEAT.COLUMNS)
  
  #Set up input frame
  DEEP.2HG<-data.table()
  
  #Go through all random samples every 5 samples
  STD.SEQ<-c(seq(1,50,5), seq(51,100,10), seq(101,1000,100), seq(1001,20000,1000))
  FIVE.FEAT<-intersect(STD.SEQ, 1:length(FEAT.COLUMNS))
  if (!(length(FEAT.COLUMNS) %in% FIVE.FEAT)){
    FIVE.FEAT<-c(FIVE.FEAT, length(FEAT.COLUMNS))
  }
  
  for (f in FIVE.FEAT[2:length(FIVE.FEAT)]){
    
    hidden.1<-f
    hidden.2<-round(f/2)
    HIDDEN<-c(hidden.1, hidden.2)
    id=0.3
    HIDDEN.DR<-c(0.5,0.5)
    METHOD="TanhWithDropout"
    
    for (n in 1:5){
      print (c("current:", f, n ))
      
      MODEL.2HG<-h2o.deeplearning(x=1:f, y=TARGET.COLUMN, data=h2o_MET[,c(ran.feat, TARGET.COLUMN)], classification = T, nfolds = 5,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                  input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )
      
      CUR.PRED<-data.table(TRAIN.ACC=1-MODEL.2HG@model$train_class_error, TEST.ACC=1-MODEL.2HG@model$valid_class_error, ITER=n, METHOD=METHOD, FEATURES=f,
                           HIDDEN=paste(HIDDEN,collapse="."), INPUT.DR=id, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
      
      #Assign predictors
      DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
      
      #Clean H2o memory
      h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))  
    } 
  }
  
  #Clean up and return
  DEEP.2HG$SEED<-seed
  return(DEEP.2HG)
}

Function.MET.COR.DEEP<-function(met.obj, cor.obj){
  
  #Load cor.obj
  met.cor<-readRDS(cor.obj)
  met.cor<-met.cor[order(COR,decreasing = T),]
 
  ##########SET UP###########
  #Load met table
  MET.TABLE<-readRDS(met.obj)
  
  #Open h2o connection
  localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '32g', nthreads=-1) 
  
  #Introduce met obj as h2o object
  key.v<-sample(letters,1)
  MET.TABLE<-MET.TABLE[sample(nrow(MET.TABLE)),] #Too balanace classes!!!!
  h2o_MET<-as.h2o(localH2O, MET.TABLE, key=key.v) 
  ###########################
  
  ########EXECUTE###########
  TARGET.COLUMN<-"MET"
  
  #Get gene walk order and filter for those genes that are actually in feature table
  gene.order<-unique(met.cor$Hugo_Symbol)
  gene.order<-gene.order[gene.order %in% colnames(MET.TABLE)]
  
  #Go through all correlation genes every n genes
  STD.SEQ<-c(seq(1,50,5), seq(51,100,10), seq(101,1000,100), seq(1001,20000,1000))
  FIVE.FEAT<-intersect(STD.SEQ, 1:length(gene.order))
  if (!(length(gene.order) %in% FIVE.FEAT)){
    FIVE.FEAT<-c(FIVE.FEAT, length(gene.order))
  }
  
  #Set up input frame
  DEEP.2HG<-data.table()
  
  #Iterate through ordered genes in correlation
  for (f in FIVE.FEAT[2:length(FIVE.FEAT)]){
    
    hidden.1<-f
    hidden.2<-round(f/2)
    HIDDEN<-c(hidden.1, hidden.2)
    id=0.3
    HIDDEN.DR<-c(0.5,0.5)
    METHOD="TanhWithDropout"
    
    for (n in 1:5){
      print (c("current:", f, n ))
      
      MODEL.2HG<-h2o.deeplearning(x=1:f, y=TARGET.COLUMN, data=h2o_MET[,c(gene.order, TARGET.COLUMN)], classification = T, nfolds = 5,
                                  activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 500,
                                  input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )
      
      CUR.PRED<-data.table(TRAIN.ACC=1-MODEL.2HG@model$train_class_error, TEST.ACC=1-MODEL.2HG@model$valid_class_error, ITER=n, METHOD=METHOD, FEATURES=f,
                           HIDDEN=paste(HIDDEN,collapse="."), INPUT.DR=id, HIDDEN.DR=paste(HIDDEN.DR,collapse="."))  
      
      #Assign predictors
      DEEP.2HG<-rbind(DEEP.2HG, CUR.PRED)
      
      #Clean H2o memory
      h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))  
    } 
  }
  return(DEEP.2HG)
}

#Arguments
args<-commandArgs(trailingOnly=T)
met.obj<-args[1]
output.file<-args[2]
hidden<-args[3] #hidden layers as character string "80.40.10"
hidden.dr<-args[4] #hidden.dr layers as denomitor character strings "2.2.2" for "0.5 0.5 0.5" (1/2)
input.dr<-args[5] #input dropout rate for input layer
method<-args[6] #choose between "Tanh" or ""TanhWithDropout"
cor.obj<-args[7]
print("opened files")

print (met.obj)
#Execute
#main.obj<-Function.MET.RANDOM.DEEP(met.obj, method=method, hidden=hidden, hidden.dr=hidden.dr, input.dr=input.dr)
#main.obj<-Function.MET.RANDOM.DEEP(met.obj)
main.obj<-Function.MET.COR.DEEP(met.obj, cor.obj)

#Write out
saveRDS(main.obj, output.file)
print ("done writing to output")