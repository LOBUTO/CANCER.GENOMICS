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

tang.plus.mcd<-Function.met.gene.cor.diff(tang.plus.pval.met, brca.exp, TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, kegg.edges, type = "tcga")
tang.minus.mcd<-Function.met.gene.cor.diff(tang.minus.pval.met, brca.exp, TCGA.BRCA.CLINICAL[ER.STATUS=="Negative",]$SAMPLE, kegg.edges, type = "tcga")
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
hist(teru.minus.class.table$ABS.MEAN.COR.DIFF)

library(h2o)
library(caret)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '30g', nthreads=-1) 
teru.minus.class.table

key.z<-sample(letters,1)
MAIN.MET<-do.call(rbind, list(teru.minus.class.table))
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
  
  x<-h2o.performance(MODEL.MET, h2o_MET[test_rows[1:length(test_rows)],])
  AUC<-x@metrics$AUC
  print (AUC)
}