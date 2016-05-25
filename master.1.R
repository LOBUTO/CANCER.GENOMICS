library(data.table)
library(reshape2)

tcga.drug.feat <- readRDS("~/Documents/FOLDER/OBJECTS/050916.TCGA.DRUG.FEAT.rds")
nci.cgp.feat.class <- readRDS("~/Documents/FOLDER/OBJECTS/052516.NCI.CGP.FEAT.CLASS.rds")
all.tcga.clinical <- readRDS("~/Documents/FOLDER/OBJECTS/050916.ALL.TCGA.CLINICAL.rds")
cancer.samples <- readRDS("~/Documents/FOLDER/OBJECTS/052016.CANCER.SAMPLES.rds")

target.samples <- unique(cancer.samples[CANCER=="Positive",]$SAMPLE)

for (i in c(0)){
  print (i)

  #ITERATE THROUGH RANDOM SAMPLING
  for (iter in 1){

    test.table <- tcga.drug.feat$FEAT.TABLE[SAMPLE %in% target.samples,]
    setnames(test.table, colnames(cgp.cor.AUC))
    temp.table <- nci.cgp.feat.class
    train.rows <- sample(1:nrow(temp.table), nrow(temp.table)*0.8)
    valid.rows <- setdiff(1:nrow(temp.table), train.rows)
    train.table <- temp.table[train.rows,]
    valid.table <- temp.table[train.rows,]

    test.table <- test.table[order(cell_name),] #Need to order to get survival curves later!!!

    target.drug <- "ALL"
    file.name <- paste0("CGP.AUC_DRUG.TH_", i ,"_",target.drug, "_TARGET",".",iter)
    file.name <- paste0("~/Documents/FOLDER/TABLES/TCGA.TRAINING/", file.name)

    write.table(test.table, paste0(file.name, ".VALID", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )
    write.table(test.table, paste0(file.name, ".TEST", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )

    #Write indices
    batch.size <- 100
    train.index <- rep(batch.size, nrow(train.table)%/%batch.size)

    l.batch <- nrow(train.table)%%batch.size
    if (l.batch>80){
      train.index <- c(train.index, l.batch)
    }
    #train.index <- as.vector(table(train.table$Compound))
    write.table(train.table, paste0(file.name, ".TRAIN", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )
    write.table(train.index, paste0(file.name, ".TRAIN.INDEX.1"), row.names = F, quote = F, col.names = F)

    #Write clinical table in valid table order
    write.table(all.tcga.clinical[SAMPLE %in% test.table$cell_name,][order(SAMPLE),],
                "~/Documents/FOLDER/TABLES/TCGA.TRAINING/master.clinical.txt",
                sep="\t", quote=F, row.names = F, col.names=T )
  }
}
