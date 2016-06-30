library(data.table)
library(reshape2)

Function.scale.data.table<-function(d.t, col.protect=1:3, criteria=c(), scaling="SCALE"){
  #Scaling can be regular z-scoring ("SCALE") or 0-1 range ("ZERO")

  y.colnames <- colnames(d.t[,-col.protect,with=F])

  if (length(criteria)>0){
    all.criteria <- unique(d.t[[criteria]])

    all.y <- lapply(all.criteria, function(x)  {
      z <- d.t[d.t[[criteria]]==x,]
      z <- data.frame(z[,-col.protect,with=F])
      col.except <- apply(z, 2, sd)==0
      not.col.except <- apply(z, 2, sd)!=0

      if (scaling=="SCALE"){
        z <- cbind(z[,col.except], scale(z[,not.col.except]))
      } else if (scaling=="ZERO"){
        z <- cbind(z[,col.except], apply(z[,not.col.except], 2 ,function(zz) Function.range.0.1(zz)) )
      }

      return (z)
    })
    y <- do.call(rbind, all.y)

  } else {
    y<-data.frame(d.t[,-col.protect,with=F])

    if (scaling=="SCALE"){
      y<-scale(y)
    } else if (scaling=="ZERO"){
      y <- apply(y, 2, function(zz) Function.range.0.1(zz))
    }

  }

  x.colnames<-colnames(d.t[,col.protect,with=F])
  x<-cbind(d.t[,col.protect,with=F], y)
  setnames(x, c(x.colnames, y.colnames))

  return(x)
}

OBJ_FOLDER <- "~/Documents/FOLDER/OBJECTS/" #For Lab
TCGA.TRAINING <- "~/Documents/FOLDER/TABLES/TCGA.TRAINING/"
OBJ_FOLDER <- "/tigress/zamalloa/OBJECTS/" #For tigress
TCGA.TRAINING <- "/tigress/zamalloa/TABLES/TCGA.TRAINING/"

#tcga.drug.feat <- readRDS("~/Documents/FOLDER/OBJECTS/050916.TCGA.DRUG.FEAT.rds")
#tcga.drug.feat.filtered <- readRDS("~/Documents/FOLDER/OBJECTS/061616.TCGA.DRUG.FEAT.FILTERED.2.rds")
#nci.cgp.feat.class <- readRDS("~/Documents/FOLDER/OBJECTS/061716.NCI.CGP.FEAT.CLASS.rds")
all.tcga.clinical <- readRDS(paste0(OBJ_FOLDER,"050916.ALL.TCGA.CLINICAL.rds"))
cancer.samples <- readRDS(paste0(OBJ_FOLDER,"052016.CANCER.SAMPLES.rds"))
cgp.cor.AUC <- readRDS(paste0(OBJ_FOLDER,"052516.CGP.COR.AUC.rds"))

tcga.drug.feat <- readRDS(paste0(OBJ_FOLDER,"050916.TCGA.DRUG.FEAT.UNSCALED_FEAT.rds"))
nci.cgp.feat.class <- readRDS(paste0(OBJ_FOLDER,"061816.NCI.CGP.FEAT.CLASS.UNSCALED_FEAT.rds"))

target.samples <- unique(cancer.samples[CANCER=="Positive",]$SAMPLE)

#Establish tables
target.table <- tcga.drug.feat$FEAT.TABLE[SAMPLE %in% target.samples,] #or tcga.drug.feat$FEAT.TABLE
training.table <- nci.cgp.feat.class
setnames(target.table, colnames(cgp.cor.AUC))

#Is scaling necessary?
main.scaling=F
if (main.scaling==T){

  pre.table <- rbind(training.table, target.table)
  pre.table <- Function.scale.data.table(pre.table, col.protect=1:3)

  training.table <- pre.table[cell_name %in% unique(training.table$cell_name),]
  target.table <- pre.table[cell_name %in% unique(target.table$cell_name),]
}

#training.table <- Function.scale.data.table(training.table, col.protect=1:3)

for (i in c(0)){
  print (i)

  #ITERATE THROUGH RANDOM SAMPLING
  for (iter in 1){

    test.table <- target.table
    test.table$NORM_AUC <- scale(test.table$NORM_AUC)
    print (dim(test.table))

    temp.table <- training.table
    print (dim(temp.table))
    print (0.9)
    train.rows <- sample(1:nrow(temp.table), nrow(temp.table)*0.9)
    valid.rows <- setdiff(1:nrow(temp.table), train.rows)
    train.table <- temp.table[train.rows,]
    valid.table <- temp.table[valid.rows,]

    test.table <- test.table[order(cell_name),] #Need to order to get survival curves later!!!

    target.drug <- "ALL"
    file.name <- paste0("CGP.AUC_DRUG.TH_", i ,"_",target.drug, "_TARGET",".",iter)
    file.name <- paste0(TCGA.TRAINING, file.name)

    write.table(valid.table, paste0(file.name, ".VALID", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )
    write.table(test.table, paste0(file.name, ".TEST", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )

    #Write indices
    # batch.size <- 100
    # train.index <- rep(batch.size, nrow(train.table)%/%batch.size)
    #
    # l.batch <- nrow(train.table)%%batch.size
    # if (l.batch>80){
    #   train.index <- c(train.index, l.batch)
    # }
    #train.index <- as.vector(table(train.table$Compound))
    write.table(train.table, paste0(file.name, ".TRAIN", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )
    #write.table(train.index, paste0(file.name, ".TRAIN.INDEX.1"), row.names = F, quote = F, col.names = F)

    #Write clinical table in valid table order
    write.table(all.tcga.clinical[SAMPLE %in% test.table$cell_name,][order(SAMPLE),],
                paste0(TCGA.TRAINING,"master.clinical.txt"),
                sep="\t", quote=F, row.names = F, col.names=T )

  }
}

print ("Done writing tables")
