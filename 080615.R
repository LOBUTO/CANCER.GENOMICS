library(data.table)
library(reshape2)
library(ggplot2)
library(igraph)
library(pheatmap)
library(gplots)
library(dplyr)

teru.gene.exp<-readRDS("DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds")

tang.matrix<-Function.Tang("DATABASES/METABOLOMICS/TANG.2014/CLEANED.PROFILES.csv", 
                           "DATABASES/METABOLOMICS/TERUNUMA.2014/MET.KEGG.HMDB.csv",
                           "DATABASES/METABOLOMICS/TANG.2014/NOT.FOUND.ANNOTATED")

teru.normal.matrix<-Function.teru.met.matrix("DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", "DATABASES/METABOLOMICS/TERUNUMA.2014/ID.TO.GSM.csv",
                                             "DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds", "DATABASES/METABOLOMICS/TERUNUMA.2014/MET.KEGG.HMDB.csv",
                                             "DATABASES/METABOLOMICS/TANG.2014/NOT.FOUND.ANNOTATED", CANCER=F)

teru.cancer.matrix<-Function.teru.met.matrix("DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", "DATABASES/METABOLOMICS/TERUNUMA.2014/ID.TO.GSM.csv",
                                             "DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds", "DATABASES/METABOLOMICS/TERUNUMA.2014/MET.KEGG.HMDB.csv",
                                             "DATABASES/METABOLOMICS/TANG.2014/NOT.FOUND.ANNOTATED", CANCER=T, "DATABASES/METABOLOMICS/TERUNUMA.2014/TERU.ID.STATUS.csv")

dim(teru.normal.matrix$MATRIX)
dim(teru.cancer.matrix$MATRIX)

teru.matrix<-cbind(teru.cancer.matrix[intersect(rownames(teru.cancer.matrix), rownames(teru.normal.matrix)),],teru.normal.matrix[intersect(rownames(teru.cancer.matrix),rownames(teru.normal.matrix)),] )
teru.class.C00042<-Function.met.gene.pre.classify.teru(teru.matrix, teru.gene.exp$MATRIX, "C00042", 2)

#
TCGA.BRCA.CLINICAL<-Function.tcga.brca.clinical("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
saveRDS(TCGA.BRCA.CLINICAL, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/081415.BRCA.CLINICAL.rds")

#
met.gold.std<-Function.met.gold.std(teru.normal.matrix, teru.cancer.matrix, tang.matrix, TCGA.BRCA.CLINICAL, 0.05, cutoff.th ="LFC", lfc.th = 2)

tcga.mut<-Function.tcga.mut.prep("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/030415/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", mut.filter=1) #at least 1 people have all mutations
tcga.mut

kegg.edges<-Function.kegg.filtered("DATABASES/KEGG/063015.ENZYME.SUB.PROD", "DATABASES/RECON/042215.PROCESSED.METABOLITES", weight.filter = 60, n.edge.filter = 100)
saveRDS(kegg.edges, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/081415.KEGG.EDGES.rds")

gene.length<-fread("DATABASES/UNIPROT/042314_GENE_LENGTH", header=T, sep="\t", stringsAsFactors = F)

#
met.score.1<-Function.met.score.1(tcga.mut, kegg.edges, gene.length, TCGA.BRCA.CLINICAL)

met.score.1$MUT.POS[MET %in% c("C00122", "C00042", "C01087"),]
ggplot(met.score.1$MUT.POS, aes(MET.SCORE)) + geom_histogram() + theme.format

met.score.1.neg.pr<-Function.met.pr(met.score.1$MUT.NEG[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.NEG[!is.na(met.gold.std$MET.NEG)])
met.score.1.pos.pr<-Function.met.pr(met.score.1$MUT.POS[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.POS[!is.na(met.gold.std$MET.POS)])

ggplot(met.score.1.neg.pr, aes(FPR,TPR)) + geom_point() + geom_line() + theme.format
ggplot(met.score.1.pos.pr, aes(FPR,TPR)) + geom_point() + geom_line() + theme.format

#
brca.exp<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041715.BRCA.RNASEQ.MATRICES.V2.RSEM.UQ.rds")

met.score.2<-Function.met.score.2(brca.exp, kegg.edges, TCGA.BRCA.CLINICAL, 0.1, 2)

met.score.2.neg.pr<-Function.met.pr(met.score.2$EXP.NEG[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.NEG[!is.na(met.gold.std$MET.NEG)])
met.score.2.pos.pr<-Function.met.pr(met.score.2$EXP.POS[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.POS[!is.na(met.gold.std$MET.POS)])

ggplot(met.score.2.neg.pr$PRED.TABLE, aes(FPR,TPR)) + geom_point() + geom_line()
ggplot(met.score.2.pos.pr$PRED.TABLE, aes(FPR,TPR)) + geom_point() + geom_line()

#
met.score.3<-Function.met.score.3(tcga.mut, brca.exp, TCGA.BRCA.CLINICAL, pval.th = 0.1, fold.th = 1,met.sample.th = 5) #DONE IN CLUSTER
met.score.3<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/081415.BRCA.MET.SCORE.3.TH.5.rds")

met.score.3.neg.pr<-Function.met.pr(met.score.3$EXP.NEG[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.NEG[!is.na(met.gold.std$MET.NEG)])
met.score.3.pos.pr<-Function.met.pr(met.score.3$EXP.POS[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.POS[!is.na(met.gold.std$MET.POS)])

ggplot(met.score.3.neg.pr$PRED.TABLE, aes(FPR,TPR)) + geom_point() + geom_line()
ggplot(met.score.3.pos.pr$PRED.TABLE, aes(FPR,TPR)) + geom_point() + geom_line()
met.score.3.neg.pr$AUC
met.score.3.pos.pr$AUC

#
met.score.1
met.score.2
met.score.3

alpha<- 0
met.score.1.3.neg<-merge(met.score.2$EXP.NEG, met.score.3$EXP.NEG, by="MET")
met.score.1.3.neg$MET.SCORE<-alpha*met.score.1.3.neg$MET.SCORE.x + (1-alpha)*met.score.1.3.neg$MET.SCORE.y

met.score.1.3.neg.pr<-Function.met.pr(met.score.1.3.neg[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                    met.gold.std$MET.NEG[!is.na(met.gold.std$MET.NEG)])
ggplot(met.score.1.3.neg.pr$PRED.TABLE, aes(FPR,TPR)) + geom_point() + geom_line()
met.score.1.3.neg.pr$AUC

f <- function(alpha)  {
  met.score.1.3.neg$MET.SCORE<-alpha*met.score.1.3.neg$MET.SCORE.x + (1-alpha)*met.score.1.3.neg$MET.SCORE.y
  
  met.score.1.3.neg.pr<-Function.met.pr(met.score.1.3.neg[MET %in% union(rownames(tang.matrix), rownames(teru.cancer.matrix$MATRIX)),], 
                                        met.gold.std$MET.NEG[!is.na(met.gold.std$MET.NEG)])
  return (met.score.1.3.neg.pr$AUC)
}
plot(seq(-1,1,0.1), sapply(seq(-1,1,0.1), function(y) f(y)))

#
tang.plus.lfc.gene<-Function.gene.diff.exp(brca.exp, intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE) , type = "tcga")
tang.minus.lfc.gene<-Function.gene.diff.exp(brca.exp, intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE) , type = "tcga")

teru.gene.exp<-readRDS("DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds")
teru.plus.lfc.gene<-Function.gene.diff.exp(teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, "teru"  )
teru.minus.lfc.gene<-Function.gene.diff.exp(teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE, "teru"  )

tang.plus.pval.met<-Function.met.diff.exp(tang.matrix, TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, "tang")
tang.minus.pval.met<-Function.met.diff.exp(tang.matrix, TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, "tang")
tang.pval.met<-Function.met.diff.exp(tang.matrix, TCGA.BRCA.CLINICAL$SAMPLE, "tang")
teru.plus.pval.met<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE , "teru", teru.normal.matrix$MATRIX)
teru.minus.pval.met<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE , "teru", teru.normal.matrix$MATRIX)

kegg.path<-fread("DATABASES/KEGG/060415_PATHWAY_TO_COMPOUND", header=T, sep="\t")

tang.plus.mcd<-Function.met.gene.cor.diff(tang.plus.pval.met, brca.exp, TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, kegg.edges, type = "tcga") #FIX SAMPLES
tang.minus.mcd<-Function.met.gene.cor.diff(tang.minus.pval.met, brca.exp, TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, kegg.edges, type = "tcga") #FIX SAMPLES
teru.plus.mcd<-Function.met.gene.cor.diff(teru.plus.pval.met, teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, kegg.edges, type = "teru")
teru.minus.mcd<-Function.met.gene.cor.diff(teru.minus.pval.met, teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE, kegg.edges, type = "teru")

tang.plus.class.table<-Function.prep.kegg.pred.table(kegg.edges, tang.plus.lfc.gene, tang.plus.pval.met, kegg.path, tang.plus.mcd, pval.th = 0.05, lfc.th = 1)
tang.minus.class.table<-Function.prep.kegg.pred.table(kegg.edges, tang.minus.lfc.gene, tang.minus.pval.met, kegg.path, tang.minus.mcd, pval.th = 0.05, lfc.th = 1)

ggplot(tang.plus.class.table, aes(DIFF, NORM.DEGREE)) + geom_boxplot() + scale_y_log10()
ggplot(tang.plus.class.table, aes(DIFF, HUGO.MED.LFC)) + geom_boxplot() + scale_y_log10()
ggplot(tang.plus.class.table, aes(DIFF, BC)) + geom_boxplot() + scale_y_sqrt()
ggplot(tang.plus.class.table, aes(DIFF, PATH.COUNT)) + geom_boxplot() + scale_y_sqrt()
ggplot(tang.plus.class.table, aes(DIFF, ABS.MEAN.COR.DIFF)) + geom_boxplot() + scale_y_sqrt()

teru.plus.class.table<-Function.prep.kegg.pred.table(kegg.edges, teru.plus.lfc.gene, teru.plus.pval.met, kegg.path, teru.plus.mcd, pval.th = 0.05, lfc.th = 1)
teru.minus.class.table<-Function.prep.kegg.pred.table(kegg.edges, teru.minus.lfc.gene, teru.minus.pval.met, kegg.path, teru.minus.mcd, pval.th = 0.05, lfc.th = 1)

ggplot(teru.plus.class.table, aes(DIFF, NORM.DEGREE)) + geom_boxplot() + scale_y_log10()
ggplot(teru.plus.class.table, aes(DIFF, HUGO.MED.LFC)) + geom_boxplot() + scale_y_log10()
ggplot(teru.plus.class.table, aes(DIFF, BC)) + geom_boxplot() + scale_y_sqrt()

tang.master.table<-Function.class.master.sampler(intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL$SAMPLE), sample.size = 10, n = 10, type = "tcga",pval.th = 0.05, lfc.th = 1)
tang.master.table.minus<-Function.class.master.sampler(intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE), 
                                                       sample.size = 5, n = 10, type = "tcga",pval.th = 0.05, lfc.th = 1)
tang.master.table.plus<-Function.class.master.sampler(intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE), 
                                                       sample.size = 5, n = 10, type = "tcga",pval.th = 0.05, lfc.th = 1)

teru.master.table<-Function.class.master.sampler(teru.cancer.matrix$ER.STATUS$SAMPLE, sample.size = 10, n = 10, type = "teru",pval.th = 0.05, lfc.th = 1)
#
library(h2o)
library(caret)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 
tang.plus.class.table

key.z<-sample(letters,1)
MAIN.MET<-do.call(rbind, list(tang.plus.class.table))
MAIN.MET<-MAIN.MET[,2:ncol(MAIN.MET), with=F][sample(nrow(MAIN.MET)),]
h2o_MET<-as.h2o(localH2O, data.frame(MAIN.MET), destination_frame = key.z) 
FEATURES<-setdiff(colnames(h2o_MET), c("MET", "DIFF"))

h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, key.z))

n_folds<-5
rand_folds<-createFolds(as.factor(as.matrix(h2o_MET$DIFF)), k=n_folds)
train_rows<-as.numeric(unlist(rand_folds[1:4]))
test_rows<-as.numeric(unlist(rand_folds[5]))

AUC<-0.5
while(AUC<0.75){
  MODEL.MET<-h2o.deeplearning(x=FEATURES, y="DIFF",  h2o_MET[train_rows[1:length(train_rows)],], use_all_factor_levels = T,
                              input_dropout_ratio = 0.1, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
                              activation = "RectifierWithDropout", balance_classes = F, hidden=c(200,200,200), epochs=500)
  
  #x<-h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],])
  #AUC<-x@metrics$AUC
  AUC<-auc(as.data.frame(h2o_MET[test_rows[1:length(test_rows)],])[["DIFF"]], as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["p1"]])
  print (AUC)
}

plot(h2o.performance(MODEL.MET, h2o_MET[train_rows[1:length(train_rows)],]), type="roc")
plot(h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]), type="roc")
confusionMatrix(as.data.frame(h2o_MET[test_rows[1:length(test_rows)],])[["DIFF"]], as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["predict"]])
table(as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["p1"]]>0.544591)

table(as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["predict"]])


h2o_TEST<-as.h2o(localH2O, as.data.frame(tang.minus.class.table[,2:ncol(tang.minus.class.table),with=F]), destination_frame = "zz")

h2o_MET[,2]<-as.factor(h2o_MET[,2])
h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],])
x<-h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],])
h2o:::h2o.find_threshold_by_max_metric(x, "f2")

x@metrics$max_criteria_and_metric_scores
x@metrics$thresholds_and_metric_scores

h2o.confusionMatrix(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],])
table(as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["predict"]])

table(as.data.frame(h2o.predict(MODEL.MET, h2o_TEST[1:nrow(h2o_TEST),]))[["predict"]])

#
target.pred<-prediction(as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["p1"]], as.data.frame(h2o_MET[test_rows[1:length(test_rows)],])[["DIFF"]] )
target.perf<-performance(target.pred, "tpr", "fpr")
pred.table<-data.table(FPR=target.perf@x.values[[1]], TPR=target.perf@y.values[[1]])
ggplot(pred.table, aes(FPR, TPR)) + geom_point() + geom_line()

pred.auc<-auc(as.data.frame(h2o_MET[test_rows[1:length(test_rows)],])[["DIFF"]], as.data.frame(h2o.predict(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]))[["p1"]])
plot(h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],]), type="roc")
#

h2o.saveModel(MODEL.MET, "file:///Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/", "tang.plus.class.best.3")

h2o_TEST<-as.h2o(localH2O, data.frame(tang.minus.class.table))
h2o.performance(MODEL.MET, h2o_TEST)
plot(h2o.performance(MODEL.MET, h2o_TEST, type="roc"))

master.sampler.perf<-data.table()

data.sets<-c("tang.master.table", "tang.master.table.minus", "tang.master.table.plus")
for (dataset in data.sets){
  
  #Prep dataset
  MAIN.MET<-get(dataset)
  MAIN.MET<-MAIN.MET[sample(nrow(MAIN.MET)),]
  h2o_MET<-as.h2o(localH2O, data.frame(MAIN.MET), destination_frame = key.z) 
  FEATURES<-setdiff(colnames(h2o_MET), c("MET", "DIFF"))
  
  n_folds<-10
  rand_folds<-createFolds(as.factor(as.matrix(h2o_MET$DIFF)), k=n_folds)
  train_rows<-as.numeric(unlist(rand_folds[1:9]))
  test_rows<-as.numeric(unlist(rand_folds[10]))
  
  INPUT.DROPOUT<-c(0,0.1,0.2)
  HIDDEN.DROPOUT.RATIOS<-list(rep(0.1,3), rep(0.3,3), rep(0.5,3))
  ITERATIONS=1:10
  
  #Learn
  for (i in ITERATIONS){
    for (id in INPUT.DROPOUT){
      for (hdp in HIDDEN.DROPOUT.RATIOS){
        print (c(i, id, hdp))
        
        MODEL.TRAIN<-h2o.deeplearning(x=FEATURES, y="DIFF",
                                      input_dropout_ratio = id, hidden_dropout_ratios =hdp,
                                      h2o_MET[train_rows[1:length(train_rows)],], use_all_factor_levels = T,
                                      activation = "RectifierWithDropout", balance_classes = F, hidden=c(200,200,200), epochs=500)  
        
        #Get performance
        train.perf<-h2o.performance(MODEL.TRAIN, h2o_MET[train_rows[1:length(train_rows)],])@metrics$AUC
        test.perf<-h2o.performance(MODEL.TRAIN, h2o_MET[test_rows[1:length(test_rows)],])@metrics$AUC
        master.sampler.perf<-rbind(master.sampler.perf, 
                                   data.table(DATA=dataset, ITERATION=i, TRAIN.AUC=train.perf, TEST.AUC=test.perf, HIDDEN.DROPOUT=paste0(hdp, collapse ="."),
                                              INPUT.DROPOUT=id))
        
        #Save model
        model.id<-paste0(dataset,"_",i,"_",id,"_",paste0(hdp, collapse ="."), "_", round(test.perf,3))
        h2o.saveModel(MODEL.TRAIN, 
                      "file:///Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/", 
                      model.id)
        
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, key.z))   
      }
    }
  }
}

saveRDS(master.sampler.perf, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/081815.TANG.MET.PREDICT.RESULTS.rds")
melt(master.sampler.perf, measure.vars = c("TRAIN.AUC", "TEST.AUC"))
ggplot(melt(master.sampler.perf, measure.vars = c("TRAIN.AUC", "TEST.AUC")), aes(HIDDEN.DROPOUT, value, colour=variable)) + geom_boxplot() + theme.format +
  facet_grid(DATA~INPUT.DROPOUT) + theme(strip.text.x = element_text(size = 14), strip.text.y = element_text(size = 14))

       
MODEL.TRAIN<-h2o.deeplearning(x=FEATURES, y="DIFF", 
                              h2o_MET[train_rows[1:length(train_rows)],], use_all_factor_levels = T,
                              activation = "Rectifier", balance_classes = F, hidden=c(100,100,100), epochs=500)  

MODEL.TRAIN<-h2o.deeplearning(x=FEATURES, y="DIFF", model_id = "td", 
                              input_dropout_ratio = 0, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
                              h2o_MET[train_rows[1:length(train_rows)],], use_all_factor_levels = T,
                              activation = "RectifierWithDropout", balance_classes = F, hidden=c(200,200,200), epochs=500)  

h2o.performance(MODEL.TRAIN, h2o_MET[test_rows[1:length(test_rows)],])

plot(h2o.performance(MODEL.TRAIN, h2o_MET[train_rows[1:length(train_rows)],]), type="roc")
plot(h2o.performance(MODEL.TRAIN, h2o_MET[test_rows[1:length(test_rows)],]), type="roc")

#
BEST.TANG.MASTER.PLUS.MODEL<-
  h2o.loadModel("file:///Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/tang.master.table.plus_9_0_0.1.0.1.0.1_0.918",
                                           localH2O)
BEST.TANG.MASTER.PLUS.MODEL@parameters$training_frame<-"holder.1"
plot(h2o.performance(BEST.TANG.MASTER.PLUS.MODEL, h2o_TEST), type="roc")
h2o.performance(BEST.TANG.MASTER.PLUS.MODEL, h2o_TEST)

BEST.TANG.MASTER.MINUS.MODEL<-
  h2o.loadModel("file:///Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/tang.master.table.minus_3_0_0.1.0.1.0.1_0.943",
                                            localH2O)
BEST.TANG.MASTER.MINUS.MODEL@parameters$training_frame<-"holder.2"

BEST.TANG.MASTER.MODEL<-
  h2o.loadModel("file:///Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/tang.master.table_10_0_0.1.0.1.0.1_0.971",
                localH2O)
BEST.TANG.MASTER.MODEL@parameters$training_frame<-"holder.3"

h2o_TEST<-as.h2o(localH2O, data.frame(teru.minus.class.table))
h2o.performance(BEST.TANG.MASTER.MINUS.MODEL, h2o_TEST)

#Test best model on TCGA non-tang samples and get predicted metabolites
tcga.breast.plus.pred<-Function.tcga.pred.met(brca.exp, BEST.TANG.MASTER.PLUS.MODEL, setdiff(TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, colnames(tang.matrix)), tang.pval.met,
                                       kegg.edges, kegg.path, localH2O)
tcga.breast.minus.pred<-Function.tcga.pred.met(brca.exp, BEST.TANG.MASTER.MINUS.MODEL, setdiff(TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, colnames(tang.matrix)), tang.pval.met,
                                        kegg.edges, kegg.path, localH2O)

tcga.breast.plus.pred.enrich<-Function.kegg.path.enrich(kegg.path, tcga.breast.plus.pred[PRED==1,]$MET)
tcga.breast.minus.pred.enrich<-Function.kegg.path.enrich(kegg.path, tcga.breast.minus.pred[PRED==1,]$MET)
tcga.breast.plus.pred.enrich[PVAL.ADJ<0.01,]
tcga.breast.minus.pred.enrich[PVAL.ADJ<0.01,]

GBM.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/081915.GBM.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
COAD.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/081915.COAD.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
READ.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/READ/081915.READ.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
UCEC.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/081915.UCEC.GA.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
LUAD.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/LUAD/081915.LUAD.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
LUSC.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/LUSC/081915.LUSC.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
OV.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/OV/081915.OV.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
KIRC.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/KIRC/081915.KIRC.RNASEQ.MATRICES.V2.RSEM.UQ.rds")
KIRP.EXP<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/KIRP/081915.KIRP.RNASEQ.MATRICES.V2.RSEM.UQ.rds")

tcga.gbm.pred<-Function.tcga.pred.met(GBM.EXP, BEST.TANG.MASTER.MODEL, colnames(GBM.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.coad.pred<-Function.tcga.pred.met(COAD.EXP, BEST.TANG.MASTER.MODEL, colnames(COAD.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.read.pred<-Function.tcga.pred.met(READ.EXP, BEST.TANG.MASTER.MODEL, colnames(READ.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.ucec.pred<-Function.tcga.pred.met(UCEC.EXP, BEST.TANG.MASTER.MODEL, colnames(UCEC.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.luad.pred<-Function.tcga.pred.met(LUAD.EXP, BEST.TANG.MASTER.MODEL, colnames(LUAD.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.lusc.pred<-Function.tcga.pred.met(LUSC.EXP, BEST.TANG.MASTER.MODEL, colnames(LUSC.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.ov.pred<-Function.tcga.pred.met(list(tumor=OV.EXP$tumor, normal=brca.exp$normal), BEST.TANG.MASTER.MODEL, colnames(OV.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.kirc.pred<-Function.tcga.pred.met(KIRC.EXP, BEST.TANG.MASTER.MODEL, colnames(KIRC.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)
tcga.kirp.pred<-Function.tcga.pred.met(KIRP.EXP, BEST.TANG.MASTER.MODEL, colnames(KIRP.EXP$tumor), tang.pval.met, kegg.edges, kegg.path, localH2O)

Function.met.pred.combine<-function(pred.list){
  #Takes list of outputs from Function.tcga.pred.met() and merges them into a single matrix
  
  #Obtain tables in list format
  main.list<-lapply(pred.list, function(x) {
    
    cancer.table<-get(x)
    setnames(cancer.table, c("MET", toupper(strsplit(x, "tcga.")[[1]][2])))
    cancer.table<-unique(cancer.table)
    return(cancer.table)
  })
  
  #Combine them
  main.pred.matrix<-Reduce(function(x, y) merge(x, y, by="MET"), main.list)
  
  #Clean up and Return
  main.pred.matrix<-data.frame(main.pred.matrix, row.names = "MET")
  main.pred.matrix<-data.matrix(main.pred.matrix)
  return(main.pred.matrix)
}

tcga.pred<-Function.met.pred.combine(c("tcga.breast.plus.pred", "tcga.breast.minus.pred", "tcga.gbm.pred", "tcga.coad.pred",
                                       "tcga.read.pred", "tcga.ucec.pred", "tcga.luad.pred", "tcga.lusc.pred", "tcga.ov.pred",
                                       "tcga.kirc.pred", "tcga.kirp.pred"))
saveRDS(tcga.pred, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/082015.TCGA.PRED.MET.1.rds")
pheatmap(tcga.pred, scale = "none")
pheatmap(cor(tcga.pred, method = "spearman"), scale="none", color=colorRampPalette(c("white","green","green4","violet","purple"))(100))

tcga.pred.sample<-melt(data.table(tcga.pred[c("C00042", "C00122", "C01087"),], keep.rownames = T), id.vars = "rn")
setnames(tcga.pred.sample, c("MET", "CANCER", "PRED"))
ggplot(tcga.pred.sample, aes(MET, PRED, fill=CANCER)) + geom_histogram(stat = "identity", position = "dodge") + scale_fill_brewer(palette="Set2") +
  theme.format

#Renormalize expression matrices
brca.diff.exp.all<-Function.process.icgc.exp.raw("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv")
saveRDS(brca.diff.exp.all, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.BRCA.DIFF.EXP.ALL.rds")

teru.gene.exp<-readRDS("DATABASES/METABOLOMICS/TERUNUMA.2014/NORMALIZED.AFFY.EXP.MATRIX.rds")
hist(teru.gene.exp$MATRIX)

teru.diff.limma<-Function.teru.diff.limma(teru.gene.exp, unique(teru.gene.exp$CLASS$SAMPLE), norm = T)

ggplot(merge(teru.diff.limma, brca.diff.exp.all, by="Hugo_Symbol"), aes(LOG.FC.x, LOG.FC.y)) + geom_point()
cor.test(merge(teru.diff.limma, brca.diff.exp.all, by="Hugo_Symbol")[["LOG.FC.x"]], merge(teru.diff.limma, brca.diff.exp.all, by="Hugo_Symbol")[["LOG.FC.y"]], method = "spearman")


brca.diff.exp.erplus<-Function.process.icgc.exp.raw("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv",
                                                    TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE)
teru.diff.limma.erplus<-Function.teru.diff.limma(teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, norm = F)

ggplot(merge(brca.diff.exp.erplus, teru.diff.limma.erplus, by="Hugo_Symbol"), aes(LOG.FC.x, LOG.FC.y)) + geom_point()
cor.test(merge(teru.diff.limma.erplus, brca.diff.exp.erplus, by="Hugo_Symbol")[["LOG.FC.x"]], merge(teru.diff.limma.erplus, brca.diff.exp.erplus, by="Hugo_Symbol")[["LOG.FC.y"]], method = "spearman")


brca.diff.exp.erminus<-Function.process.icgc.exp.raw("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv",
                                                    TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE)
teru.diff.limma.erminus<-Function.teru.diff.limma(teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE, norm = F)

ggplot(merge(brca.diff.exp.erminus, teru.diff.limma.erminus, by="Hugo_Symbol"), aes(LOG.FC.x, LOG.FC.y)) + geom_point()
cor.test(merge(teru.diff.limma.erminus, brca.diff.exp.erminus, by="Hugo_Symbol")[["LOG.FC.x"]], merge(teru.diff.limma.erminus, brca.diff.exp.erminus, by="Hugo_Symbol")[["LOG.FC.y"]], method = "spearman")

#BUILD MODEL BASED ON TERUNUMA DATASET
teru.diff.limma.erplus

teru.plus.pval.met<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE , "teru", teru.normal.matrix$MATRIX)
teru.minus.pval.met<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE , "teru", teru.normal.matrix$MATRIX)

teru.plus.mcd<-Function.met.gene.cor.diff(teru.plus.pval.met, teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, kegg.edges, type = "teru")
teru.minus.mcd<-Function.met.gene.cor.diff(teru.minus.pval.met, teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE, kegg.edges, type = "teru")

teru.plus.class.table<-Function.prep.kegg.pred.table(kegg.edges, teru.diff.limma.erplus, teru.plus.pval.met, kegg.path, teru.plus.mcd, pval.th = 0.05, lfc.th = 1)
teru.minus.class.table<-Function.prep.kegg.pred.table(kegg.edges, teru.diff.limma.erminus, teru.minus.pval.met, kegg.path, teru.minus.mcd, pval.th = 0.05, lfc.th = 1)
table(teru.minus.class.table$DIFF)

library(h2o)
library(caret)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 
saveRDS(teru.minus.class.table, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/082815.TERU.MINUS.CLASS.TABLE.rds")
saveRDS(teru.plus.class.table, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/082815.TERU.PLUS.CLASS.TABLE.rds")

#Do Non-bootstrap first
#Prep dataset
key.z<-sample(letters,1)
MAIN.MET<-teru.minus.class.table[sample(nrow(teru.minus.class.table)),]
h2o_MET<-as.h2o(localH2O, data.frame(MAIN.MET), destination_frame = key.z) 
FEATURES<-setdiff(colnames(h2o_MET), c("MET", "DIFF"))

set.seed(43)
n_folds<-5
rand_folds<-createFolds(as.factor(as.matrix(h2o_MET$DIFF)), k=n_folds)
train_rows<-as.numeric(unlist(rand_folds[1:4]))
test_rows<-as.numeric(unlist(rand_folds[5]))

train.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame =  h2o_MET[train_rows[1:length(train_rows)],],
                              input_dropout_ratio = 0.2, hidden_dropout_ratios =c(0.5,0.5,0.5) ,
                              use_all_factor_levels = T,
                              activation = "RectifierWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)   
train.gbm<-h2o.gbm(x=FEATURES, y="DIFF", training_frame =  h2o_MET[train_rows[1:length(train_rows)],], learn_rate = 0.01, ntrees = 2000)

h2o.performance(train.model, h2o_MET[test_rows[1:length(test_rows)],])
h2o.performance(train.gbm, h2o_MET[test_rows[1:length(test_rows)],])

metrics.non.bootstrap<-data.table()
set.seed(43)
n_folds<-5
rand_folds<-createFolds(as.factor(as.matrix(h2o_MET$DIFF)), k=n_folds)
input.dropout<-c("0.0","0.1","0.2")


for (ip in input.dropout){
  
  for (i in 1:length(rand_folds)){
    
    for (iter in 1:10){
      
      test_rows<-as.numeric(unlist(rand_folds[i]))
      train_rows<-setdiff(as.numeric(unlist(rand_folds)), test_rows)
      
      #Obtain deep learning model and parameters
      deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", 
                                   h2o_MET[train_rows[1:length(train_rows)],], use_all_factor_levels = T,
                                   activation = "Rectifier", balance_classes = F, hidden=c(100,100,100), epochs=500)   
                                    #input_dropout_ratio = as.numeric(ip), hidden_dropout_ratios =c(0.5,0.5,0.5) ,
      
      deep.f1<-deep.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[1] 
      
      deep.conf.matrix<-h2o.confusionMatrix(deep.model, h2o_MET[train_rows[1:length(train_rows)],] , thresholds=deep.f1)
      deep.test.acc<-deep.conf.matrix$Error[3]
      deep.test.adj.acc<-var(c(deep.conf.matrix$Error[1],deep.conf.matrix$Error[2])) + deep.test.acc
      
      #Predict
      deep.pred.conf.matrix<-h2o.confusionMatrix(deep.model, h2o_MET[test_rows[1:length(test_rows)],], thresholds=deep.f1)
      deep.train.acc<-deep.pred.conf.matrix$Error[3]
      deep.train.adj.acc<-var(c(deep.pred.conf.matrix$Error[1], deep.pred.conf.matrix$Error[2])) + deep.train.acc
      
      #Store
      metrics.non.bootstrap<-rbind(metrics.non.bootstrap, data.table(FOLD=i, ITERATION=iter, 
                                                                     TEST.ACC=deep.test.acc, TEST.ADJ.ACC=deep.test.adj.acc,
                                                                     TRAIN.ACC=deep.train.acc, TRAIN.ADJ.ACC=deep.train.adj.acc,
                                                                     INPUT.DROPOUT=0, METHOD="RECTIFIER", REGULARIZED=FALSE,
                                                                     SEED=43))
      #Save model
      model.id<-paste0("082915.TERUNUMA.MINUS","_",43,"_", i, "_", iter, "_", round(deep.test.adj.acc,3), "_", round(deep.train.adj.acc,3))
      h2o.saveModel(deep.model, 
                    paste0("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS/RECTIFIER/NO.BOOTSTRAP/NO.REG/"), 
                    model.id)
      
      h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z)))  
    }
  }
}

ggplot(melt(metrics.non.bootstrap[,c(3,5,7:10),with=F], measure.vars =  c("TEST.ACC", "TRAIN.ACC")), aes(factor(INPUT.DROPOUT), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(METHOD~REGULARIZED)
ggplot(melt(metrics.non.bootstrap[,c(4,6,7:10),with=F], measure.vars =  c("TEST.ADJ.ACC", "TRAIN.ADJ.ACC")), aes(factor(INPUT.DROPOUT), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(METHOD~REGULARIZED)

best.minus.model<-h2o.loadModel("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS/TANH/NO.BOOTSTRAP/INPUT.DROPOUT.0.1/082915.TERUNUMA.MINUS_43_3_9_0.16_0.145/", localH2O)
best.minus.model@parameters$training_frame<-"holder.1"
plot(h2o.performance(best.minus.model, tang.minus.h2o), type="roc")
h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.191111)

#Do BOOTSTRAP.1 now

minus.samples<-teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE
minus.split<-split(minus.samples, 1:5)

metrics.bootstrap<-data.table()
metrics.bootstrap.gbm<-data.table()

for (i in 1:length(minus.split)){
  
  #Prep master sampler table for train.samples
  test.samples<-minus.split[[i]]
  train.samples<-setdiff(unlist(minus.samples), test.samples)
  
  #Create master table per bootsrapped set
  bootstrap.train<-replicate(10, sample(train.samples, 10), simplify = F)
  master.table<-data.table()
  
  for (j in bootstrap.train){
    
    train.limma<-Function.teru.diff.limma(teru.gene.exp, j, norm = F)
    train.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  j , "teru", teru.normal.matrix$MATRIX)
    train.mcd<-Function.met.gene.cor.diff(train.pval, teru.gene.exp, j, kegg.edges, type = "teru")
    train.table<-Function.prep.kegg.pred.table(kegg.edges, train.limma, train.pval, kegg.path, train.mcd, pval.th = 0.05, lfc.th = 1)
    
    #Store
    master.table<-rbind(master.table, train.table)
  }
  
  #Create testing table
  test.limma<-Function.teru.diff.limma(teru.gene.exp, test.samples, norm = F)
  test.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  test.samples , "teru", teru.normal.matrix$MATRIX)
  test.mcd<-Function.met.gene.cor.diff(test.pval, teru.gene.exp, test.samples, kegg.edges, type = "teru")
  test.table<-Function.prep.kegg.pred.table(kegg.edges, test.limma, test.pval, kegg.path, test.mcd, pval.th = 0.05, lfc.th = 1)
  
  #Build model out of 10 iterations
  master.table<-master.table[sample(1:nrow(master.table)),]
  key.z<-sample(letters,1)
  train.h2o<-as.h2o(localH2O, data.frame(master.table), destination_frame = key.z)
  test.h2o<-as.h2o(localH2O, data.frame(test.table), destination_frame = "test")
  FEATURES<-setdiff(colnames(train.h2o), c("MET", "DIFF"))
  
  learn.rate<-c(0.01,0.05,0.1,0.2)
  max.depth<-c(2,4,6,8,10,20)
  n.trees<-c(5,10,15,20,50,100)
  for (l in learn.rate){
    for (m in max.depth){
      for (n in n.trees){
        #Model
        #       deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF",  train.h2o, use_all_factor_levels = T,
        #                                    input_dropout_ratio = as.numeric(ip), hidden_dropout_ratios =c(0.5,0.5,0.5) ,
        #                                    activation = "TanhWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)
        
        gbm.model<-h2o.gbm(x=FEATURES, y="DIFF",  train.h2o, learn_rate = l, max_depth =m , ntrees = n)
        gbm.train.acc<-h2o.confusionMatrix(gbm.model, train.h2o)$Error[3]
        gbm.train.adj.acc<-gbm.train.acc + var(h2o.confusionMatrix(gbm.model, train.h2o)$Error[1:2])
        
        gbm.test.acc<-h2o.confusionMatrix(gbm.model, test.h2o)$Error[3]
        gbm.test.adj.acc<-gbm.test.acc + var(h2o.confusionMatrix(gbm.model, test.h2o)$Error[1:2])
        
        #       #Obtain accuracy
        #       deep.f1<-deep.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[1] 
        #       
        #       deep.conf.matrix<-h2o.confusionMatrix(deep.model, test.h2o , thresholds=deep.f1)
        #       deep.test.acc<-deep.conf.matrix$Error[3]
        #       deep.test.adj.acc<-var(c(deep.conf.matrix$Error[1],deep.conf.matrix$Error[2])) + deep.test.acc
        #       
        #       #Predict
        #       deep.pred.conf.matrix<-h2o.confusionMatrix(deep.model, train.h2o, thresholds=deep.f1)
        #       deep.train.acc<-deep.pred.conf.matrix$Error[3]
        #       deep.train.adj.acc<-var(c(deep.pred.conf.matrix$Error[1], deep.pred.conf.matrix$Error[2])) + deep.train.acc
        
        #Save models
        #i.split<-i
        #iteration<-n
        model.id<-paste0("082815.TERUNUMA.MINUS","_",l, "_", m, "_",n,"_", round(gbm.train.acc,3), "_", round(gbm.test.acc,3))
        h2o.saveModel(gbm.model, 
                      paste0("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS/GBM/BOOTSTRAP/"), 
                      model.id)
        
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z, "test")))
        
        #Store for reference
        metrics.bootstrap.gbm<-rbind(metrics.bootstrap.gbm, data.table(FOLD=i, LEARN.RATE=l, MAX.DEPTH=m, N.TREES=n,
                                                               TRAIN.ACC=gbm.train.acc, TRAIN.ADJ.ACC=gbm.train.adj.acc,
                                                               TEST.ACC=gbm.test.acc, TEST.ADJ.ACC=gbm.test.adj.acc,
                                                                METHOD="GBM", REGULARIZED=TRUE))
      }
    }
  }
}

ggplot(melt(metrics.bootstrap[,c(3,5,7:9), with=F], measure.vars = c("TRAIN.ACC", "TEST.ACC")), aes(factor(INPUT.DROPOUT), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(METHOD~REGULARIZED) + theme.format +  theme(strip.text.x = element_text(size = 14))
ggsave("~/Documents/FOLDER/LAB/NOTES.FIGURES/083015.TERUNUMA.ERMINUS.MODEL.DL.tiff", width = 12, height = 16, dpi = 600)

ggplot(melt(metrics.bootstrap.gbm[,c(2:4, 5,7),with=F], measure.vars = c("TRAIN.ACC", "TEST.ACC")), aes(factor(LEARN.RATE), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(MAX.DEPTH~N.TREES)
metrics.bootstrap.gbm[order(TEST.ACC),]

best.minus.model<-h2o.loadModel("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS/GBM/BOOTSTRAP/082815.TERUNUMA.MINUS_0.1_8_20_0.068_0.124",localH2O)
best.minus.model@parameters$training_frame<-"holder.4"
h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.555767)


#####Test best model####
icgc.obj<-Function.process.icgc.matrix.to.obj("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv")
tang.diff.exp.minus<-Function.process.icgc.exp.raw("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", 
                                                   "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv",
                                                   intersect(TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, colnames(tang.matrix)))
tang.minus.pval.met<-Function.met.diff.exp(tang.matrix, intersect(colnames(tang.matrix),TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE), "tang")
tang.minus.mcd<-Function.met.gene.cor.diff(tang.minus.pval.met, icgc.obj, intersect(colnames(tang.matrix),TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE), 
                                           kegg.edges, type = "tcga") 
tang.minus.class.table<-Function.prep.kegg.pred.table(kegg.edges, tang.diff.exp.minus, tang.minus.pval.met, kegg.path, tang.minus.mcd, pval.th = 0.05, lfc.th = 1)
tang.minus.h2o<-as.h2o(tang.minus.class.table, localH2O, destination_frame = "tang.minus")

best.minus.model<-h2o.loadModel("file:///Users/jzamalloa/Dropbox/Lab/PIPELINES/METABOLIC.DRIVERS/OBJECTS/H2O/NEW.MODELS/RECTIFIER/BOOTSTRAP.10/INPUT.DROPOUT.0.1/082515.TERUNUMA.MINUS.X_3_1_0.965_0.833",localH2O)
best.minus.model@parameters$training_frame<-"holder.4"
best.minus.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[1]

h2o.performance(best.minus.model, tang.minus.h2o)
h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.470639)
plot(h2o.performance(best.minus.model, tang.minus.h2o), type="roc")


c<-h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.470639)
var(c(c$Error[1],c$Error[2])) + c$Error[3]

list.files("~/Dropbox/Lab/PIPELINES/METABOLIC.DRIVERS/OBJECTS/H2O/NEW.MODELS/RECTIFIER/BOOTSTRAP.10/INPUT.DROPOUT.0.1/")

minus.model.data<-data.table()
for (h2o.model in list.files("~/Dropbox/Lab/PIPELINES/METABOLIC.DRIVERS/OBJECTS/H2O/NEW.MODELS/RECTIFIER/BOOTSTRAP.10/NO.REG/")){
  
  #Obtain model stats
  model.file<-paste0("file:///Users/jzamalloa/Dropbox/Lab/PIPELINES/METABOLIC.DRIVERS/OBJECTS/H2O/NEW.MODELS/RECTIFIER/BOOTSTRAP.10/NO.REG/",h2o.model )
  best.minus.model<-h2o.loadModel(model.file,localH2O)
  best.minus.model@parameters$training_frame<-"holder.4"
  model.f1<-best.minus.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[1] 
  model.f2<-best.minus.model@model$training_metrics@metrics$max_criteria_and_metric_scores$threshold[2] 
  
  #Predict
  c.f1<-h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=model.f1)
  c.f2<-h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=model.f2)
  c.f1.score<-var(c(c.f1$Error[1],c.f1$Error[2])) + c.f1$Error[3]
  c.f2.score<-var(c(c.f2$Error[1],c.f2$Error[2])) + c.f2$Error[3]
  
  #Store
  minus.model.data<-rbind(minus.model.data, data.table(METHOD="Rectifier", REGULARIZED=FALSE, INPUT.DROPOUT=0.0, MODEL.ID=h2o.model, F1.SCORE=c.f1.score, F2.SCORE=c.f2.score))
}
minus.model.data[order(F1.SCORE),]

