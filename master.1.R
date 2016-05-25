library(data.table)
library(reshape2)

for (i in c(0.6)){
  print (i)

  selected.drugs <- unique(nci.gi50$NSC)

  #ITERATE THROUGH RANDOM SAMPLING
  for (iter in 1){

    test.table <- tcga.drug.feat$FEAT.TABLE[DRUG %in% target.drug,][SAMPLE %in% target.samples,]
    setnames(test.table, colnames(cgp.cor.AUC))
    train.table <- nci.cgp.feat[Compound %in% selected.drugs, ]

    test.table <- test.table[order(cell_name),] #Need to order to get survival curves later!!!
    #test.table$NORM_AUC <- -test.table$NORM_AUC

    #train.table <- train.table[sample(nrow(train.table)),]
    #train.table$NORM_AUC <- scale(train.table$NORM_AUC)
    train.table <- train.table[order(Compound),]

    target.drug <- "ALL"
    file.name <- paste0("CGP.AUC_DRUG.TH_", i ,"_",target.drug, "_TARGET",".",iter)
    file.name <- paste0("PIPELINES/METABOLIC.DRIVERS/TABLES/DRUG.PRED/", "SMART.BATCH.COR", "/", file.name)

    write.table(test.table, paste0(file.name, ".TEST", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )

    #Write indices
    #     batch.size <- 100
    #     train.index <- rep(batch.size, nrow(train.table)%/%batch.size)
    #
    #     l.batch <- nrow(train.table)%%batch.size
    #     if (l.batch>40){
    #       train.index <- c(train.index, l.batch)
    #     }
    train.index <- as.vector(table(train.table$Compound))
    write.table(train.table, paste0(file.name, ".TRAIN", ".1"), sep = "\t", quote = F, row.names = F, col.names = T )
    write.table(train.index, paste0(file.name, ".TRAIN.INDEX.1"), row.names = F, quote = F, col.names = F)

    #Write clinical table in valid table order
    write.table(all.tcga.clinical[SAMPLE %in% test.table$cell_name,][order(SAMPLE),],
                "PIPELINES/METABOLIC.DRIVERS/TABLES/DRUG_THEANO_PRED/MLP/REGRESSION/TESTING/TRYING.FIGURES/master.clinical.txt",
                sep="\t", quote=F, row.names = F, col.names=T )
  }
}
