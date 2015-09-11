#BOOTSTRAP.2
#Boostrapping samples for 9/10 of metabolites across all individuals and testing on the rest of metabolites
library(h2o)
library(caret)
library(randomForest)

#The limit will be at teru.minus.mcd
Function.met.split<-function(limit.met, seed, pval.met, pval.th = 0.05, lfc.th = 1, folds=10){
  require(data.table)
  require(caret)
  
  set.seed(seed)
  
  #Filter table by limit.met
  pval.met<-pval.met[MET %in% limit.met,]
  met.sig<-pval.met[PVAL.ADJ<pval.th & abs(LFC)>lfc.th,]$MET
  pval.met$DIFF<-as.factor(ifelse(pval.met$MET %in% met.sig, 1,0))
  
  #Split table
  met.splits<-createFolds(pval.met$DIFF, folds)
  met.splits<-lapply(met.splits, function(x) pval.met$MET[x])
  
  #Return
  return(met.splits)
}

minus.met.splits<-Function.met.split(teru.minus.mcd$MET, 23, teru.minus.pval.met, folds = 10)
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 

metrics.bootstrap.2<-data.table()
for (i in 1:length(minus.met.splits)){
  
  #Split mets
  testing.mets<-minus.met.splits[[i]]
  training.mets<-setdiff(unlist(minus.met.splits), testing.mets)

  #Obtain bootstrapped samples
  er.minus.samples<-teru.cancer.matrix$ER.STATUS[ER.STATUS=="NEG",]$SAMPLE
  bootstrap.train<-replicate(10, sample(er.minus.samples, 10), simplify = F)
  
  #Create master table per bootsrapped set (for bootstrapped samples and training mets)
  master.table<-data.table()
  
  for (j in bootstrap.train){
    
    train.limma<-Function.teru.diff.limma(teru.gene.exp, j, norm = F)
    train.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  j , "teru", teru.normal.matrix$MATRIX)
    train.mcd<-Function.met.gene.cor.diff(train.pval, teru.gene.exp, j, kegg.edges, type = "teru")
    train.table<-Function.prep.kegg.pred.table(kegg.edges, train.limma, train.pval, kegg.path, train.mcd[MET %in% training.mets, ], pval.th = 0.05, lfc.th = 1)
    
    #Store
    master.table<-rbind(master.table, train.table)
  }
  
  #Create train table (with all samples and testing mets, no bootstrap)
  test.limma<-Function.teru.diff.limma(teru.gene.exp, er.minus.samples, norm = F)
  test.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  er.minus.samples, "teru", teru.normal.matrix$MATRIX)
  test.mcd<-Function.met.gene.cor.diff(test.pval, teru.gene.exp, er.minus.samples, kegg.edges, type = "teru")
  test.table<-Function.prep.kegg.pred.table(kegg.edges, test.limma, test.pval, kegg.path, test.mcd[MET %in% testing.mets,], pval.th = 0.05, lfc.th = 1)
  
  #Build model
  master.table<-master.table[sample(1:nrow(master.table)),]
  key.z<-sample(letters,1)
  train.h2o<-as.h2o(localH2O, data.frame(master.table), destination_frame = key.z)
  test.h2o<-as.h2o(localH2O, data.frame(test.table), destination_frame = "test")
  FEATURES<-setdiff(colnames(train.h2o), c("MET", "DIFF"))
  
  #Model
#   deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame =  train.h2o, use_all_factor_levels = T,
#                                input_dropout_ratio = 0.0, hidden_dropout_ratios =c(0.2,0.2,0.2) ,
#                               activation = "RectifierWithDropout", balance_classes = F, hidden=c(100,100,100), epochs=500)
  
  learn.rate<-c(0.005,0.01,0.05,0.1,0.2)
  max.depth<-c(2,4,6,8,10,20)
  n.trees<-c(5,10,15,20,50,100,200)
  for (l in learn.rate){
    for (m in max.depth){
      for (n in n.trees){
        
        gbm.model<-h2o.gbm(x=FEATURES, y="DIFF", training_frame =  train.h2o, learn_rate = l, ntrees = n, max_depth = m)
        
        #GBM model predictions 
        gbm.conf.matrix<-h2o.confusionMatrix(gbm.model, train.h2o)
        gbm.train.acc<-gbm.conf.matrix$Error[3]
        gbm.train.adj.acc<-var(gbm.conf.matrix$Error[1:2]) + gbm.train.acc
        
        #Predict
        gbm.pred.conf.matrix<-h2o.confusionMatrix(gbm.model, test.h2o)
        gbm.test.acc<-gbm.pred.conf.matrix$Error[3]
        gbm.test.adj.acc<-var(gbm.pred.conf.matrix$Error[1:2]) + gbm.test.acc
        
        #Store
        metrics.bootstrap.2<-rbind(metrics.bootstrap.2, data.table(FOLD=i, LEARN.RATE=l, MAX.DEPTH=m, N.TREES=n,
                                                                   GBM.TRAIN.ACC=gbm.train.acc, GBM.TRAIN.ADJ.ACC=gbm.train.adj.acc,
                                                                   GBM.TEST.ACC=gbm.test.acc, GBM.TEST.ADJ.ACC=gbm.test.adj.acc))
        #Save models
        model.id<-paste0("083015.TERUNUMA.MINUS","_",i,"_",l, "_", m, "_",n,"_", round(gbm.train.acc,3), "_", round(gbm.test.acc,3))
        h2o.saveModel(gbm.model, 
                      paste0("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS.2/GBM/BOOTSTRAP/"), 
                      model.id)
        
        print (metrics.bootstrap.2)
        h2o.rm(localH2O, setdiff(h2o.ls(localH2O)$key, c(key.z, "test"))) 
      }
    }
  }
}

metrics.bootstrap.2
ggplot(melt(metrics.bootstrap.2[,c(2:4,5,7),with=F], measure.vars = c("GBM.TRAIN.ACC", "GBM.TEST.ACC")), aes(factor(LEARN.RATE), 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(MAX.DEPTH~N.TREES)

metrics.bootstrap.2[order(GBM.TEST.ACC),]

best.minus.model<-h2o.loadModel("//home/lobuto/Documents/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/MET.PREDICTION/H2O/CORRECTED.MODELS.2/GBM/BOOTSTRAP/083015.TERUNUMA.MINUS_8_0.1_4_200_0.034_0.062",localH2O)
best.minus.model@parameters$training_frame<-"holder.4"
h2o.confusionMatrix(best.minus.model, tang.minus.h2o, thresholds=0.549989)

##########REMOVE DIFFERENTIAL FEATURES AND TRAIN AGAIN###########
teru.diff.limma.erplus<-Function.teru.diff.limma(teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, norm = F)

teru.plus.pval.met<-Function.met.diff.exp(teru.cancer.matrix$MATRIX, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE , "teru", teru.normal.matrix$MATRIX)
teru.plus.mcd<-Function.met.gene.cor.diff(teru.plus.pval.met, teru.gene.exp, teru.cancer.matrix$ER.STATUS[ER.STATUS=="POS",]$SAMPLE, kegg.edges, type = "teru")
teru.plus.class.table<-Function.prep.kegg.pred.table(kegg.edges, teru.diff.limma.erplus, teru.plus.pval.met, kegg.path, teru.plus.mcd, pval.th = 0.05, lfc.th = 1)

teru.plus.class.table

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 

rand.rows<-createFolds(teru.plus.class.table$DIFF, 5)
plus.h2o<-as.h2o(localH2O, data.frame(teru.plus.class.table[sample(nrow(teru.plus.class.table)),]), "teru.plus")
FEATURES<-setdiff(colnames(teru.plus.class.table), c("DIFF", "MET", "HUGO.MED.LFC", "ABC", "CYS.MET", "ALA.ASP.GLU", "AA",
                                                     "GST", "GLYCO.TCA.LIPID", "PURINE",
                                                     "TWO.OXO", "ANTIBIO", "ARG.PRO", "CANCER", 
                                                     "CARBON", "PROTEIN.ABS", "AMINOACYL.TRNA", "STEROID"))
FEATURES<-setdiff(colnames(teru.plus.class.table), c("DIFF", "MET", "HUGO.MED.LFC"))

n_folds<-5
rand_folds<-createFolds(as.factor(as.matrix(plus.h2o$DIFF)), k=n_folds)
train_rows<-as.numeric(unlist(rand_folds[1:4]))
test_rows<-as.numeric(unlist(rand_folds[5]))

deep.model<-h2o.deeplearning(x=FEATURES, y="DIFF", training_frame = plus.h2o[train_rows,],
                             activation = "Tanh", hidden = c(100,100,100), epochs = 100)
h2o.performance(deep.model, plus.h2o[test_rows,])

h2o.gbm(x=FEATURES, y="DIFF", training_frame = plus.h2o[train_rows,],
        ntrees = 16, max_depth = 40, learn_rate = 0.005)

h2o.rm(setdiff(h2o.ls(localH2O)$key, "teru.plus"), localH2O)
h2o.ls(localH2O)

#Test best plus model
teru.plus.gbm.models<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/DEEP.ITER.1/090215.TERU.PLUS.GBM.MODEL.rds")
ggplot(melt(teru.plus.gbm.models[,c(2:4, 5,7), with=F], measure.vars = c("GBM.TRAIN.ACC", "GBM.TEST.ACC") ), aes(N.TREES, 1-value, colour=variable)) +
  geom_boxplot() + facet_grid(LEARN.RATE~MAX.DEPTH)
teru.plus.gbm.models[order(GBM.TEST.ADJ.ACC),]
teru.plus.best.model<-h2o.loadModel("file:///Users/jzamalloa/FOLDER/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/DEEP.ITER.1/090215.TERUNUMA.PLUS_5_0.05_6_100_0.085_0.133", localH2O)
teru.plus.best.model@parameters$training_frame<-"holder.1"

icgc.obj<-Function.process.icgc.matrix.to.obj("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv")
tang.diff.exp.plus<-Function.process.icgc.exp.raw("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082315.ICGC.RAW.COUNTS.MATRIX.rds", 
                                                   "DATABASES/CANCER_DATA/ICGC/BRCA/specimen.BRCA-US.tsv",
                                                   intersect(TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE, colnames(tang.matrix)))
tang.plus.pval.met<-Function.met.diff.exp(tang.matrix, intersect(colnames(tang.matrix),TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE), "tang")
tang.plus.mcd<-Function.met.gene.cor.diff(tang.plus.pval.met, icgc.obj, 
                                          intersect(colnames(tang.matrix),TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE), 
                                          kegg.edges, type = "tcga") 
tang.plus.class.table<-Function.prep.kegg.pred.table(kegg.edges, tang.diff.exp.plus, tang.plus.pval.met, kegg.path, tang.plus.mcd, pval.th = 0.05, lfc.th = 1)

h2o.performance(teru.plus.best.model, as.h2o(data.frame(tang.plus.class.table), localH2O))

teru.plus.deep.models<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/DEEP.ITER.1/090315.TERU.PLUS.DEEP.MODEL.rds")
ggplot(melt(teru.plus.deep.models[,c(2:4, 5,7), with=F], measure.vars = c("GBM.TRAIN.ACC", "GBM.TEST.ACC") ), 
       aes(factor(HIDDEN), 1-value, colour=variable)) +geom_boxplot() + facet_grid(INPUT.DROPOUT~HIDDEN.DROPOUT)
teru.plus.deep.models[order(GBM.TEST.ACC),]

teru.plus.best.model<-h2o.loadModel("file:///Users/jzamalloa/FOLDER/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/DEEP.ITER.1/090315.TERUNUMA.PLUS_5_80_0_0.3_0.15_0.133", localH2O)
teru.plus.best.model@parameters$training_frame<-"holder.1"
h2o.performance(teru.plus.best.model, as.h2o(data.frame(tang.plus.class.table), localH2O))

#Exploratory analysis - Do metabolites that are closer in metabolic network exhibit higher correlation (distance vs cor)
dim (teru.cancer.matrix)
kegg.edges

kegg.rel<-unique(rbind(data.table(MET.1=kegg.edges$SUBSTRATE, MET.2=kegg.edges$PRODUCT),
                       data.table(MET.1=kegg.edges$PRODUCT, MET.2=kegg.edges$SUBSTRATE)))
kegg.rel<-kegg.rel[MET.1 %in% rownames(teru.cancer.matrix$MATRIX) & MET.2 %in% rownames(teru.cancer.matrix$MATRIX),]

kegg.test<-kegg.edges[,list(MET=unique(c(SUBSTRATE, PRODUCT))), by="Hugo_Symbol"]
kegg.test<-kegg.test[MET %in% rownames(teru.cancer.matrix$MATRIX),]

kegg.layers<-kegg.test[,list(HUGO.1=length(unique(Hugo_Symbol)), 
                             HUGO.2= length(union(kegg.test[MET==kegg.rel[MET.1==MET,]$MET.2,]$Hugo_Symbol, Hugo_Symbol))  ), by="MET"]
ggplot(kegg.layers, aes(HUGO.1, HUGO.2)) + geom_point()
ggplot(kegg.layers, aes(HUGO.2/HUGO.1)) + geom_histogram() + scale_x_log10()

kegg.layers[order(HUGO.2/HUGO.1, decreasing = F),][1:60,]

##
teru.cancer.matrix$MATRIX[1:3,1:3]
teru.cancer.matrix$ER.STATUS

teru.cancer.samples<-intersect(colnames(teru.cancer.matrix$MATRIX), teru.cancer.matrix$ER.STATUS[]$SAMPLE)
teru.stage.cor<-
  data.table(MET=rownames(teru.cancer.matrix$MATRIX),
             COR=apply(teru.cancer.matrix$MATRIX, 1, function(x) {
               x.cor<-cor.test(log2(x[teru.cancer.samples]),teru.cancer.matrix$ER.STATUS[SAMPLE == teru.cancer.samples,]$STAGE, method = "spearman")$estimate
               return(x.cor)
             }),
             PVAL=apply(teru.cancer.matrix$MATRIX, 1, function(x) {
               x.cor<-cor.test(log2(x[teru.cancer.samples]),teru.cancer.matrix$ER.STATUS[SAMPLE == teru.cancer.samples,]$STAGE, method = "spearman")$p.value
               return(x.cor)
             }))

teru.stage.cor<-teru.stage.cor[!is.na(COR),]
teru.stage.cor$PVAL.ADJ<-p.adjust(teru.stage.cor$PVAL, method = "fdr")

hist(teru.stage.cor$PVAL.ADJ)
teru.stage.cor[PVAL.ADJ<0.2,]
teru.stage.cor[order(abs(COR), decreasing = T),][21:30,]
hist(teru.stage.cor$COR)

teru.cancer.matrix$ER.STATUS

C08281<-as.data.table(t(teru.cancer.matrix$MATRIX["C00042", teru.cancer.samples, drop=F]), keep.rownames = T)
setnames(C08281, c("SAMPLE", "MET"))
ggplot(merge(teru.cancer.matrix$ER.STATUS,C08281, by="SAMPLE"), aes(factor(STAGE), log2(MET))) + geom_boxplot()

tang.matrix[1:3,1:3]
TCGA.BRCA.CLINICAL$STAGE.2<-ifelse(TCGA.BRCA.CLINICAL$STAGE=="Stage.I", 1, 
                                   ifelse(TCGA.BRCA.CLINICAL$STAGE=="Stage.II", 2,
                                          ifelse(TCGA.BRCA.CLINICAL$STAGE=="Stage.III", 3, 4)))
tang.cancer.samples<-intersect(colnames(tang.matrix), TCGA.BRCA.CLINICAL[ER.STATUS=="Positive",]$SAMPLE)

tang.stage.cor<-
  data.table(MET=rownames(tang.matrix),
             COR=apply(tang.matrix, 1, function(x) {
               x.cor<-cor.test(log2(x[tang.cancer.samples]), TCGA.BRCA.CLINICAL[SAMPLE==tang.cancer.samples,]$STAGE.2, method = "spearman")$estimate
               return(x.cor)
             }),
             PVAL=apply(tang.matrix, 1, function(x) {
               x.cor<-cor.test(log2(x[tang.cancer.samples]), TCGA.BRCA.CLINICAL[SAMPLE==tang.cancer.samples,]$STAGE.2, method = "spearman")$p.value
               return(x.cor)
             }))

tang.stage.cor<-tang.stage.cor[!is.na(COR),]
tang.stage.cor$PVAL.ADJ<-p.adjust(tang.stage.cor$PVAL, method="fdr")

hist(tang.stage.cor$PVAL.ADJ)
tang.stage.cor[order(abs(COR), decreasing = T),][1:20,]

ggplot(merge(teru.stage.cor[,1:2,with=F], tang.stage.cor[,1:2, with=F], by="MET"), aes(abs(COR.x- COR.y))) + geom_histogram()
merge(teru.stage.cor[,1:2,with=F], tang.stage.cor[,1:2, with=F], by="MET")[order(abs(COR.x - COR.y)),][1:30,]

####
teru.cancer.matrix$MATRIX[1:3,1:3]
teru.gene.exp$MATRIX[1:3,1:3]

Function.construct.class.matrix<-function(met.matrix, gene.matrix,target.samples, met, type="teru"){
  
  #Obtain samples of interest
  common.samples<-intersect(intersect(colnames(met.matrix), colnames(gene.matrix)), target.samples)
  
  #Process met data
  met.data<-data.table(t(met.matrix[met, common.samples,drop=F]), keep.rownames = T)
  setnames(met.data, c("SAMPLE", "METABOLITE"))
  
  #Process gene data
  gene.data<-data.table(t(gene.matrix[,common.samples, drop=F]), keep.rownames = T)
  gene.col<-colnames(gene.data)
  setnames(gene.data, c("SAMPLE", gene.col[2:length(gene.col)]))
  
  #Combine data
  main.table<-merge(met.data, gene.data, by="SAMPLE")
  
  #Return
  return(main.table)
}

teru.feat.table<-Function.construct.class.matrix(teru.cancer.matrix$MATRIX, teru.gene.exp$MATRIX, colnames(teru.cancer.matrix$MATRIX), "C00327")
kegg.edges[PRODUCT=="C00327",]$Hugo_Symbol
teru.feat.table[,c("METABOLITE", intersect(unique(kegg.edges[PRODUCT=="C01921",]$Hugo_Symbol), colnames(teru.feat.table))),with=F]

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 

teru.feat.table.h2o<-as.h2o(teru.feat.table[,c("METABOLITE", intersect(union(kegg.edges[PRODUCT=="C00327",]$Hugo_Symbol, kegg.edges[SUBSTRATE=="C00327",]$Hugo_Symbol), colnames(teru.feat.table))),with=F], 
                            localH2O, destination_frame = "z")

test.deep<-h2o.deeplearning(setdiff(colnames(teru.feat.table.h2o), "METABOLITE"), "METABOLITE", teru.feat.table.h2o, 
                            activation = "Tanh",hidden = c(50,50), epochs = 400)

hist(teru.feat.table$METABOLITE)

hist(teru.cancer.matrix$MATRIX["C00989",])
hist(teru.normal.matrix$MATRIX["C00989",])

common.mets<-intersect(rownames(teru.cancer.matrix$MATRIX), rownames(teru.normal.matrix$MATRIX))
filt.mets<-sapply(common.mets, function(z) sd(teru.cancer.matrix$MATRIX[z,])==0 & sd(teru.normal.matrix$MATRIX[z, ])==0)
common.mets<-common.mets[!filt.mets]
x<-data.table(MET=common.mets, PVAL=sapply(common.mets, function(y) {
  t.met<-t.test(teru.cancer.matrix$MATRIX[y, ], teru.normal.matrix$MATRIX[y,], paired = F, var.equal = F)$p.value
  return(t.met)
}), FC=sapply(common.mets, function(y){
  fc.met<-log2(median(teru.cancer.matrix$MATRIX[y, ])/ median(teru.normal.matrix$MATRIX[y,]))
  return(fc.met)
}))

x$PVAL.ADJ<-p.adjust(x$PVAL, method = "bonferroni")
x[PVAL.ADJ<0.05,][order(FC),]

hist(teru.gene.exp$MATRIX)
hist(scale(teru.gene.exp$MATRIX, center = T, scale = T ))
hist(rowMeans(scale(teru.gene.exp$MATRIX, center = T, scale = T)))

hist(icgc.obj$tumor)
hist(scale(icgc.obj$tumor, center = T, scale=T))

dim(teru.gene.exp$MATRIX)
dim(icgc.obj$tumor)
dim(brca.exp$tumor)
length(intersect(rownames(teru.gene.exp$MATRIX), rownames(icgc.obj$tumor)))
length(intersect(rownames(teru.gene.exp$MATRIX), rownames(brca.exp$tumor)))


