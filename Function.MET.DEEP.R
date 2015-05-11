#Function.MET.DEEP.R
#050915
#Calculate prediciton errors of metabolite levels using h2o's deep learning module

library(data.table)
library(reshape2)
library(h2o)

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
  
  for (n in c(1)) {
    
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

#Arguments
args<-commandArgs(trailingOnly=T)
met.obj<-args[1]
output.file<-args[2]
hidden<-args[3] #hidden layers as character string "80.40.10"
hidden.dr<-args[4] #hidden.dr layers as denomitor character strings "2.2.2" for "0.5 0.5 0.5" (1/2)
input.dr<-args[5] #input dropout rate for input layer
method=args[6] #choose between "Tanh" or ""TanhWithDropout"
print("opened files")

#Execute
main.obj<-Function.Main(met.obj, method=method, hidden=hidden, hidden.dr=hidden.dr, input.dr=input.dr)

#Write out
saveRDS(main.obj, output.file)
print ("done writing to output")