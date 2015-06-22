library(data.table)
library(ggplot2)
library(reshape2)
library(parallel)
library(h2o)
library(MASS)
library(glmnet)
library(caret)

Function.brca.maf.simple<-function(maf.file){
  
  maf<-fread(maf.file, header=T, sep="\t", stringsAsFactors = F)
  maf<-maf[,c("Hugo_Symbol", "Start_Position","End_Position", "Variant_Classification", "Variant_Type" ,
              "Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"), with=F]
  
  #Classify mutations
  maf$MUT.TYPE<-ifelse(maf$Variant_Type %in% c("INS", "DEL"), "LOF", 
                       ifelse(maf$Variant_Classification %in% c("Missense_Mutation", "RNA") , "MISSENSE", 
                       ifelse(maf$Variant_Classification=="Silent", "SILENT", "LOF")))
  maf$Variant_Classification<-NULL
  maf$Variant_Type<-NULL
  
  #Remove patients with metastatic info
  maf<-maf[!(grepl("-06A", Tumor_Sample_Barcode)),]
  
  #Constrict sample name
  maf$SAMPLE<-sapply(maf$Tumor_Sample_Barcode, function(x) paste(unlist(strsplit(x, "-"))[1:3], collapse = ".") )
  
  #Clean up and return
  setkey(maf)
  maf<-unique(maf)
  return(maf)
}

brca.maf<-Function.brca.maf.simple("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf")

TANG.2HG<-data.table(t(TANG.2HG), keep.rownames = T)
setnames(TANG.2HG, c("SAMPLE", "R2HG"))
TANG.2HG$SAMPLE<-paste("TCGA",TANG.2HG$SAMPLE, sep = ".")
TANG.2HG<-merge(brca.maf, TANG.2HG, by="SAMPLE")

TANG.2HG[Hugo_Symbol %in% c("PHGDH","ADHFE1", "D2HGDH"),]
TANG.2HG$HIGH<-TANG.2HG$R2HG>2.5

#TEST: Correlate non-synonymous mutated genes with 2HG levels
tang.pvals<-sapply(unique(TANG.2HG$Hugo_Symbol), function(x) {
  print (x)
  
  hugo.samples<-unique(TANG.2HG[Hugo_Symbol==x & MUT.TYPE!="SILENT",]$SAMPLE)
  
  if (length(hugo.samples)>0){
    other.samples<-setdiff(unique(TANG.2HG$SAMPLE), hugo.samples)
    
    r2hg.hugo<-unique(TANG.2HG[SAMPLE %in% hugo.samples,][,c("SAMPLE", "R2HG"),with=F])$R2HG
    r2hg.other<-unique(TANG.2HG[SAMPLE %in% other.samples,][,c("SAMPLE", "R2HG"),with=F])$R2HG
    
    w<-wilcox.test(r2hg.hugo, r2hg.other)$p.value
    return(w)  
  } else{
    w<-1
    return(w)
  }
})
tang.pvals<-data.table(Hugo_Symbol=unique(TANG.2HG$Hugo_Symbol), PVAL=tang.pvals)
tang.pvals<-tang.pvals[PVAL!=1,]
tang.pvals$PVAL.ADJ<-p.adjust(tang.pvals$PVAL, method="fdr")
tang.pvals<-tang.pvals[order(PVAL.ADJ),]
#CONCLUSION: No association of mutated genes with tang 2-hg levels

#TEST: Correlate expresion levels in terunuma dataset with 2HG leves using DEEP LEARNING
TERU.2HG<-Function.classify.teru.target.met("DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", c("2-hydroxyglutarate"), 
                                            "DATABASES/METABOLOMICS/TERUNUMA.2014/ID.TO.GSM.csv")
TERU.2HG.EXP<-Function.match.teru.expression(TERU.2HG, "DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds", 
                                             "DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", 
                                             "all", normal=F)
ggplot(TERU.2HG.EXP, aes(ER.STATUS, METABOLITE))+geom_boxplot()

#Pick "smart features"
TERU.2HG.EXP[1:3,1:9, with=F]
TERU.2HG.COR<-apply(data.matrix(TERU.2HG.EXP[,5:ncol(TERU.2HG.EXP), with=F]), 2, function(x) {
  pval<-cor.test(x, TERU.2HG.EXP$METABOLITE, type="spearman")$p.value
  return(pval)
  })
TERU.2HG.COR<-data.table(Hugo_Symbol=colnames(TERU.2HG.EXP[,5:ncol(TERU.2HG.EXP), with=F]), PVAL=TERU.2HG.COR)
TERU.2HG.COR$PVAL.ADJ<-p.adjust(TERU.2HG.COR$PVAL, method="fdr") #1618 genes are significantly correlated to 2HG
TERU.2HG.COR<-TERU.2HG.COR[order(PVAL.ADJ),] #Ordered in matrix
TERU.2HG.COR<-TERU.2HG.COR[Hugo_Symbol %in% colnames(TANG.2HG.EXP),]
TERU.2HG.EXP<-TERU.2HG.EXP[,c("METABOLITE", "ER.STATUS", TERU.2HG.COR[PVAL.ADJ<0.05,]$Hugo_Symbol), with=F]

saveRDS(object = TERU.2HG.EXP, file = "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/053015.TERU.2HG.EXP.MA.TANG.FILT.rds")
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
key.z<-sample(letters,1)
TERU.2HG.EXP<-TERU.2HG.EXP[sample(nrow(TERU.2HG.EXP)),]
h2o_2HG<-as.h2o(localH2O, TERU.2HG.EXP, key=key.z) 

DEEP.2HG<-data.table()
METHOD="Tanh"
HIDDEN<-c(50,25,5)
INPUT.DR<-0.1
HIDDEN.DR<-rep(0.5, length(HIDDEN))
FEATURES<-50

for (id in INPUT.DR){
  for (n in 1:10) {
    #Model with/out dropout
    print (c("building model", id, n))
    
    MODEL.2HG<-h2o.deeplearning(x=2:ncol(TERU.2HG.EXP), y=1, data=h2o_2HG, classification = F, nfolds = 5,
                                activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 300)
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
    h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$Key,key.v))
  }  
}

ggplot(DEEP.2HG, aes(factor(FEATURES), TR.SQR.ERROR, colour=factor(INPUT.DR) )) + geom_boxplot() + geom_jitter() + facet_wrap(~METHOD)
ggplot(DEEP.2HG, aes(factor(FEATURES), VALID.SQR.ERROR, colour=factor(INPUT.DR) )) + geom_boxplot() + geom_jitter() + facet_wrap(~METHOD)

CLUSTER.2HG.RESULTS<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/052815.TESTIGN.2HG.2.rds")
CLUSTER.2HG.RESULTS$FEATURES<-rep(c(10,20,30,40,50,100,200,400),each = 40)
ggplot(CLUSTER.2HG.RESULTS, aes(factor(FEATURES), TR.SQR.ERROR, colour=factor(INPUT.DR) )) + geom_boxplot() + geom_jitter() + facet_wrap(~METHOD)
ggplot(CLUSTER.2HG.RESULTS, aes(factor(FEATURES), VALID.SQR.ERROR, colour=factor(INPUT.DR) )) + geom_boxplot() + geom_jitter() + facet_wrap(~METHOD)

#TEST: Use the best model and test it on the Tanh dataset (30)
#Picture of residuals between model prediciton and actual values
TANG.2HG<-Function.classify.tang.target.met("DATABASES/METABOLOMICS/TANG.2014/clean.met.er.csv", c("2-hydroxyglutarate") )
TANG.2HG.EXP<-Function.match.tang.expression(TANG.2HG, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041815.BRCA.MATRICES.AGILENT.rds",
                                             "DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt",
                                             "all", normal=F, exp.version="ma")
TANG.2HG.EXP$ER.STATUS<-factor(ifelse(TANG.2HG.EXP$ER.STATUS=="POS", 1, 2))

TANG_2HG<-as.h2o(localH2O, TANG.2HG.EXP, key=key.v) 
MODEL.2HG<-h2o.deeplearning(x=2:100, y=1, data=h2o_2HG, classification = F, nfolds = 5,
                            activation = "TanhWithDropout", balance_classes = TRUE, hidden = c(100,50,10), epochs = 500,
                            input_dropout_ratio = 0.5 , hidden_dropout_ratios =c(0.5,0.5,0.5) )

TANG.2HG.EXP$PREDICT<-as.numeric(as.matrix(h2o.predict(MODEL.2HG, TANG_2HG)$predict))
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="green"))
#CONCLUSSION: DEEP LEARNING may not suitable for small datasets such as this one as it is not robust to outliers

###TEST: Use regression models to predict 2HG levels with L1 regularizatin to account for outliers
#TEST 1a - Select best features
TERU.2HG.EXP[1:3,1:9, with=F]
TERU.2HG.LM<-lm(METABOLITE~., TERU.2HG.EXP[,c(1,2:59),with=F])
extractAIC(TERU.2HG.LM)[2]
summary(TERU.2HG.LM)
plot(TERU.2HG.LM)

teru.dt<-data.table()
for (f in 3:59){
  teru.lm<-lm(METABOLITE~., TERU.2HG.EXP[,c(1,2:f), with=F])
  teru.aic<-extractAIC(teru.lm)[2]
  teru.adj.r<-summary(teru.lm)$adj.r.squared
  teru.dt<-rbind(teru.dt, data.table(feat=f, aic=teru.aic, adj.r=teru.adj.r))
}

ggplot(teru.dt, aes(adj.r, aic)) + geom_point() + theme.format + ggtitle("AIC vs ADJ.R.SQUARE Comparisson")

Function.predict.features.lm<-function(target.matrix, feat, target.feat){
  
  feat.max<-nrow(target.matrix)-1
  main.list<-list()
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "target.matrix", "feat", "target.feat") ,envir=environment())
  print ("Done exporting values")
  
  main.list<-parLapply(cl, 1:1000, function(x){
    print (x)
    
    #Randomize feature samples after each iteration
    features<-sample(feat)
    
    #Reset AIC
    current.AIC<-Inf
    current.feat<-c("ER.STATUS")
    count=1
    f.adj.sq<-1
    
    #Continue step-wise AIC search till we reach maximum allowed features
    while(is.na(f.adj.sq)==FALSE | count<length(features)){
      testing.feat<-features[count]
      
      f.formula<-as.formula(paste(target.feat ," ~ ", paste(c(current.feat, testing.feat), collapse= "+")))
      print (f.formula)
      f.lm<-lm(f.formula, target.matrix)
      f.aic<-extractAIC(f.lm)[2]
      f.adj.sq<-summary(f.lm)$adj.r.squared
      
      #Keep feature if lower than current AIC and update AIC
      if ((f.aic<current.AIC) & (is.na(f.adj.sq)==FALSE)){
        current.feat<-c(current.feat, testing.feat)
        current.AIC<-f.aic
      }
      print (c(current.feat, current.AIC))
      #Update count
      count<-count+1
    }
    return(current.feat)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  return(main.list)
}

TERU.2HG.SELECTED.FEATURES<-Function.predict.features.lm(data.frame(TERU.2HG.EXP), setdiff(colnames(data.frame(TERU.2HG.EXP)),"METABOLITE") ,"METABOLITE")

hist(sapply(TERU.2HG.SELECTED.FEATURES, function(x) length(x)))
sapply(TERU.2HG.SELECTED.FEATURES, function(x) extractAIC(lm(METABOLITE~.,TERU.2HG.EXP[,c("METABOLITE",x),with=F]))[2])
sort(table(unlist(TERU.2HG.SELECTED.FEATURES)))

x<-lm(METABOLITE~., TERU.2HG.EXP[,c("METABOLITE",TERU.2HG.SELECTED.FEATURES[[11]] ),with=F])
plot(lm(METABOLITE~., TERU.2HG.EXP[,c("METABOLITE",TERU.2HG.SELECTED.FEATURES[[1]]),with=F]))

y<-copy(TERU.2HG.EXP)
y$PREDICT<-predict.lm(x, data.frame(TERU.2HG.EXP))
ggplot(y, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT), colour="red") + theme.format

TANG.2HG.EXP$PREDICT<-predict.lm(x, data.frame(TANG.2HG.EXP))
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="green"))

CLUSTER.2HG.SELECTED.FEATURES<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/O53015.BEST.FEATURES.TERU.FILTERED.rds") #run 1000 interations online
CLUSTER.2HG.SELECTED.FEATURES<-data.table(t(table(unlist(CLUSTER.2HG.SELECTED.FEATURES))), keep.rownames = T)
CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:59] #Use top 58 to build model

TERU.FILTERED.MODEL<-lm(METABOLITE~., TERU.2HG.EXP[,c("METABOLITE",CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:58]),with=F])
summary(TERU.FILTERED.MODEL)

y$PREDICT<-predict.lm(TERU.FILTERED.MODEL, data.frame(TERU.2HG.EXP))
ggplot(y, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT), colour="red") + theme.format

TANG.2HG.EXP$PREDICT<-predict.lm(TERU.FILTERED.MODEL, data.frame(TANG.2HG.EXP))
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="green"))
#TEST 1a: Conclussions: We found the 58 best features that we will use to train the model

#TEST 1b - Use glmnet to find best lambda and use crossvalidation on TERU test and training data, and validate on Tang
dim(TERU.2HG.EXP)
teru.folds<-createFolds(1:nrow(TERU.2HG.EXP), k=5)
TERU.2HG.EXP$ER.STATUS<-factor(ifelse(TERU.2HG.EXP$ER.STATUS=="POS", 1, 2))
cv.fit<-cv.glmnet(x= data.matrix(TERU.2HG.EXP[unlist(teru.folds[1:4]),CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:58], with=F]), 
                 y= data.matrix(TERU.2HG.EXP[unlist(teru.folds[1:4]),"METABOLITE",with=F]), nfolds=5)

plot(cv.fit)
coef(cv.fit, s="lambda.min")
teru.train<-TERU.2HG.EXP[unlist(teru.folds[5]),]
teru.train$PREDICT<-predict(cv.fit, newx=data.matrix(teru.train[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:58],with=F]), s="lambda.1se")
ggplot(teru.train, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT), colour="red") + theme.format
predict.cv.glmnet()

TANG.2HG.EXP$PREDICT<-predict(cv.fit, newx=data.matrix(TANG.2HG.EXP[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:58],with=F]), s="lambda.1se")
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="green"))

set.seed(10)
teru.folds<-createFolds(1:nrow(TERU.2HG.EXP), k=5)

models.lm.2hg<-data.table()
for (i in teru.folds){
  
  it.folds<-i
  other.folds<-setdiff(unlist(teru.folds), it.folds)
  feat.limit<-length(other.folds)-2
  
  training.matrix<-data.matrix(TERU.2HG.EXP[other.folds,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit], with=F])
  training.target<-data.matrix(TERU.2HG.EXP[other.folds,"METABOLITE",with=F])
  testing.matrix<-data.matrix(TERU.2HG.EXP[it.folds,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit], with=F])
  testing.target<-data.matrix(TERU.2HG.EXP[it.folds,"METABOLITE",with=F])
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
  valid.predict<-predict.lm(nr.model, data.frame(data.matrix(TANG.2HG.EXP[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit],with=F])))
  valid.nrmsd<-sqrt(mean( (TANG.2HG.EXP$METABOLITE-valid.predict)^2)) / (max(TANG.2HG.EXP$METABOLITE) - min(TANG.2HG.EXP$METABOLITE))
  
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
    
    #Apply model to validation data.set (tang)
    reg.valid.pred<-predict(reg.fit, newx=data.matrix(TANG.2HG.EXP[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit],with=F]), s="lambda.min")
    reg.valid.nrmsd<- sqrt(mean((TANG.2HG.EXP$METABOLITE - reg.valid.pred)^2)) / (max(TANG.2HG.EXP$METABOLITE) - min(TANG.2HG.EXP$METABOLITE))
    assign(paste(p, "reg.valid.nrmsd", sep="."), reg.valid.nrmsd)
  }
  
  #Return MSE
  models.lm.2hg<-rbind(models.lm.2hg, data.table(tr.nrmsd, test.nrmsd, valid.nrmsd,
                                                 Lasso.reg.tr.nrmsd, Lasso.reg.test.nrmsd, Lasso.reg.valid.nrmsd,
                                                 Ridge.reg.tr.nrmsd, Ridge.reg.test.nrmsd, Ridge.reg.valid.nrmsd,
                                                 ElNET.reg.tr.nrmsd, ElNET.reg.test.nrmsd, ElNET.reg.valid.nrmsd))  
}
models.lm.2hg

#check on tang's er.status so the factor is the same as teru's er.status
TANG.2HG.EXP$PREDICT<-predict(reg.fit, newx=data.matrix(TANG.2HG.EXP[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit],with=F]), s="lambda.min")
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="Prediction")) + theme.format

TERU.2HG.EXP$PREDICT<-predict(reg.fit, newx=data.matrix(TERU.2HG.EXP[,CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:feat.limit],with=F]), s="lambda.min")
ggplot(TERU.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="Prediction"))

model.lm.2hg<-melt(models.lm.2hg)
model.lm.2hg$SET<-rep(rep(c("Train", "Test", "Tang"),each = 5 ), 4)
model.lm.2hg$REGULARIZATION<-rep(c("None", "Lasso", "Ridge", "Elastic Net"), each=15)
ggplot(model.lm.2hg, aes(REGULARIZATION, value, colour=SET )) + geom_boxplot() + theme.format
ggplot(model.lm.2hg[REGULARIZATION!="None",], aes(REGULARIZATION, value, colour=SET )) + geom_boxplot() + theme.format

#TEST 1b: Conclusions: It is useful for training and testing on Teru, but even with lm regularization, it is difficult to predict on Tang

#TEST: Try building deep learning model again on Teru and test on Tang, but this time use best features found
BEST.FEATURES<-CLUSTER.2HG.SELECTED.FEATURES[order(N,decreasing = T),]$V2[1:58]

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 

#Convert data into h2o
key.z<-sample(letters,1)
TERU.2HG.EXP<-TERU.2HG.EXP[sample(nrow(TERU.2HG.EXP)),]
h2o_2HG<-as.h2o(localH2O, TERU.2HG.EXP, key=key.z) 

DEEP.2HG<-data.table()
METHOD="TanhWithDropout"
HIDDEN<-c(60,30,10)
INPUT.DR<-c(0,0.1,0.2)
HIDDEN.DR<-rep(0.5, length(HIDDEN))
#FEATURES<-50

for (id in INPUT.DR){
  for (n in 1:10) {
    #Model with/out dropout
    print (c("building model", id, n))
    
    MODEL.2HG<-h2o.deeplearning(x=BEST.FEATURES, y="METABOLITE", data=h2o_2HG, classification = F, nfolds = 5,
                                activation = METHOD, balance_classes = TRUE, hidden = HIDDEN, epochs = 300,
                                input_dropout_ratio = id , hidden_dropout_ratios =HIDDEN.DR )
    
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
sqrt(DEEP.2HG$VALID.SQR.ERROR)/ (max(TERU.2HG.EXP$METABOLITE) - min(TERU.2HG.EXP$METABOLITE))

TANG_2HG<-as.h2o(localH2O, TANG.2HG.EXP, key=key.v) 
MODEL.2HG<-h2o.glm(x=setdiff(colnames(TERU.2HG.EXP),"METABOLITE") , y="METABOLITE", data=h2o_2HG, 
                   nfolds = 5,family = "gaussian",lambda_search = T, nlambda = 10, alpha = seq(0,1,0.25) )
MODEL.2HG.REG<-h2o.glm(x=setdiff(colnames(TERU.2HG.EXP),"METABOLITE"), y="METABOLITE", data=h2o_2HG, 
                   nfolds = 5,family = "gaussian", lambda =  0.2754207, alpha=1 )

TANG.2HG.EXP$PREDICT<-as.numeric(as.matrix(h2o.predict(MODEL.2HG.REG, TANG_2HG)$predict))
ggplot(TANG.2HG.EXP, aes(METABOLITE, METABOLITE))+ geom_point() + geom_point(aes(METABOLITE, PREDICT, colour="prediction"))

#MAIN CONCLUSION: PREDICTING TANG 2HG EXPRESSION ON TERU WILL BE VERY VERY DIFFICULT