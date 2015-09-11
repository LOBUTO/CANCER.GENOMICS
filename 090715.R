library(data.table)
library(reshape2)
library(ggplot2)
library(h2o)
library(parallel)
library(dplyr)
library(caret)
kegg.path<-fread("DATABASES/KEGG/060415_PATHWAY_TO_COMPOUND", header=T, sep="\t")

common.mets<-intersect(rownames(teru.cancer.matrix$MATRIX), rownames(teru.normal.matrix$MATRIX))

y<-log2(teru.cancer.matrix$MATRIX[common.mets,]/ apply(teru.normal.matrix$MATRIX[common.mets,],1, median))
y.samples<-data.table(MET=colnames(y), PERC=apply(y,2,  function(q) mean(q>1) ))
y.samples[order(PERC),]
hist(y.samples$PERC)

y.met<-data.table(MET=rownames(y), PERC=apply(y,1,  function(q) mean(q>1) ))
y.met[order(PERC,decreasing = T),][1:10,]
hist(y.met$PERC)

#Testing reformatted features
teru.plus.limma<-Function.teru.diff.limma(teru.gene.exp, teru.cancer.matrix$ER.STATUS[]$SAMPLE , norm = F)
teru.plus.pval<-Function.met.diff.exp(teru.cancer.matrix$MATRIX,  teru.cancer.matrix$ER.STATUS[]$SAMPLE , "teru", teru.normal.matrix$MATRIX)
teru.plus.mcd<-Function.met.gene.cor.diff(teru.plus.pval, teru.gene.exp, teru.cancer.matrix$ER.STATUS[]$SAMPLE, kegg.edges, type = "teru")
teru.plus.table<-Function.prep.kegg.pred.table(kegg.edges, teru.plus.limma, teru.plus.pval, kegg.path, teru.plus.mcd, pval.th = 0.05, lfc.th = 1)
table(teru.plus.table$DIFF)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 
teru.plus.h2o<-as.h2o(teru.plus.table, localH2O, "teru.plus.h2o")
FEATURES<-setdiff(colnames(teru.plus.table), c("DIFF", "MET", "HUGO.MED.LFC", "STEROID", "GST", "CYS.MET", "ALA.ASP.GLU"))

n_folds<-3
rand_folds<-createFolds(as.factor(as.matrix(teru.plus.h2o$DIFF)), k=n_folds)
train_rows<-as.numeric(unlist(rand_folds[1:2]))
test_rows<-as.numeric(unlist(rand_folds[3]))

deep.model<-h2o.deeplearning(x = FEATURES, y="DIFF", training_frame = teru.plus.h2o[train_rows,],
                             activation = "TanhWithDropout", hidden = c(20,20), epochs = 600,
                             input_dropout_ratio = 0, hidden_dropout_ratios = c(0.5,0.5))
deep.model
h2o.performance(deep.model, teru.plus.h2o[test_rows,])
plot(h2o.performance(deep.model, teru.plus.h2o[test_rows,]), type = "roc")

h2o.rm(setdiff(h2o.ls(localH2O)$key, c("teru.plus.h2o") ))

#####
nci60.met<-fread("DATABASES/NCI.60/WEB_DATA_METABOLON_ALL.TXT", header=T)
nci60.mut<-fread("DATABASES/NCI.60/WEB_DATA_ALL_MT.TXT", header=T)

table(nci60.mut$UNITS)
table(nci60.mut[grepl("mutat",UNITS),]$UNITS)
unique(nci60.mut[UNITS=="Presence of mutation (2=homo mut; 1=het mut; 0=wt)",]$GENE)
unique(nci60.mut[UNITS=="Presence of mutation (2=homo mut; 1=het mut; 0=wt)",]$pname)
table(nci60.mut$pname)

table(nci60.met$cellname)

ccle.drug<-fread("DATABASES/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", header=T)
length(unique(ccle.drug[[2]]))

intersect(unique(ccle.drug[[2]]), unique(nci60.met$cellname))
setdiff(unique(nci60.met$cellname), unique(ccle.drug[[2]]) )

ccle.drug[[2]][grepl("CAKI", ccle.drug[[2]])]

gi50<-data.table(read.csv("DATABASES/NCI.60/GI50_SEPT2014", head=T, strip.white=T))
intersect(unique(gi50$CELL), unique(nci60.met$cellname))
gi50
ggplot(gi50, aes(LCONC, NLOGGI50)) + geom_point()

Function.nci60.met.cast<-function(nci60.met, LOG2=T){
  require(reshape2)
  
  main.cast<-acast(nci60.met, TITLE~cellname, fun.aggregate = mean, value.var = "Value")
  
  if (LOG2==T){
    main.cast<-log2(main.cast)
  }
  
  #Return
  return(main.cast)
}
nci60.met.cast<-Function.nci60.met.cast(nci60.met)

Function.nci60.gi50.met<-function(met.cast, gi50, filter=0){
  
  #Convert to table for mergin
  met.table<-data.table(t(met.cast), keep.rownames = T)
  setnames(met.table, c("CELL", colnames(met.table)[2:ncol(met.table)]))
  
  #Extract important infor from gi50
  gi50<-gi50[,c("NSC", "LCONC", "CELL", "NLOGGI50", "PANEL"),with=F]
  
  #Merge
  main.table<-merge(gi50, met.table, by="CELL")
  
  #Do we filter by cell count?
  cell.count<-main.table[,list(N.CELL=length(unique(CELL))), by="NSC"]
  main.table<-main.table[NSC %in% cell.count[N.CELL>=filter,]$NSC,]
  
  #Return
  return(main.table)
  
}

gi50[CELL=="NCI-H23",]$NLOGGI50
ggplot(gi50[CELL=="NCI-H23",], aes(NLOGGI50)) + geom_histogram()

ggplot(gi50[NSC==17,], aes(NLOGGI50)) + geom_histogram()

nci60.gi50.met<-Function.nci60.gi50.met(nci60.met.cast, gi50, filter = 50)
dim(nci60.gi50.met)
nci60.gi50.met[1:3,1:10,with=F]

library(h2o)

localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, max_mem_size= '8g', nthreads=-1) 

nci60.1<-as.h2o(nci60.gi50.met[NSC %in% c(1),], localH2O, "test.1")

FEATURES<-setdiff(colnames(nci60.1), c("NLOGGI50"))

n_folds<-5
rand_folds<-createFolds(as.factor(as.matrix(nci60.1$NLOGGI50)), k=n_folds)
train_rows<-as.numeric(unlist(rand_folds[1:4]))
test_rows<-as.numeric(unlist(rand_folds[5]))

deep.model.1<-h2o.deeplearning(FEATURES, "NLOGGI50", nci60.1[train_rows,],
                               activation = "Rectifier", hidden = c(200,200), epochs = 500)
                               #input_dropout_ratio = 0, hidden_dropout_ratios = c(0.5,0.5))

plot(as.numeric(as.matrix(nci60.1[train_rows,]$NLOGGI50)), as.matrix(as.numeric(h2o.predict(deep.model.1, nci60.1[train_rows,]))))
h2o.performance(deep.model.1, nci60.1[test_rows,])

plot(as.numeric(as.matrix(nci60.1[test_rows,]$NLOGGI50)), as.matrix(as.numeric(h2o.predict(deep.model.1, nci60.1[test_rows,]))))

h2o.rm(setdiff(h2o.ls(localH2O)$key, "test.1"))
