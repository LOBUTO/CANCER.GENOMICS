#nci60_feat_build_pca.R
#Function to reduce dimensionality of nci60 data through patience_increase

library(data.table)
library(reshape2)
library(dplyr)

Function.fast.merging <- function(data.list, nparts=10, BY=""){

  while(length(data.list) != 1){ #Loop until everything is merged
       #define sections according to nparts
       if(length(data.list) > nparts){
            starts = seq(1, length(data.list), nparts)
            ends = seq(nparts, length(data.list), nparts)
            #starts and ends are of equal size
            #if length(data.list) divides nparts.
       if(length(ends) < length(starts)) ends = c(ends, length(data.list))
                #making sure things are even
            sections = matrix(c(starts, ends), ncol=2, byrow=FALSE)
            sections = apply(sections, 1, list)
       }else{
            sections = list(c(1, length(data.list)))
       }

 #We have the standard way inside lapply
 #lapply merges a length(sections) amount of nparts sections

 data.list = lapply(sections, function(x, data.list){
  if(is.list(x)) x = x[[1]]
  #the standard way starts ->
  part = data.list[[x[1]]]
  for(i in x[1]:x[2]){
   part = inner_join(part, data.list[[i]], by=BY, allow.cartersian=T)
  }
  #<- standard way ends
  return(part)
 }, data.list = data.list)
 #rinse and repeat until all is done
 }
 return(data.list[[1]]) #returning the merged data frame
}


Function.nci60.pca.feat <- function(nci60.hmdb.file, hmdb, nci.gi50, nci60.exp, all.tcga.exp, classify = F,
  split=0.9, OBJ_FOLDER){

  #Build metabolite feature table
  x <- fread(nci60.hmdb.file)
  x$DRUG <- as.character(x$DRUG)

  mets <- intersect(unique(hmdb$METABOLITE), unique(x$METABOLITE))

  x <- x[METABOLITE %in% mets,]
  x[,SD:=sd(TC), by="DRUG"]
  x <- x[SD!=0,]

  x[,SD:=sd(TC), by="METABOLITE"]
  x <- x[SD!=0,]

  x <- data.table(acast(x,DRUG~METABOLITE, value.var = "TC"), keep.rownames=T)
  setnames(x, c("NSC", colnames(x)[2:ncol(x)]))
  print ("processed casting of DRUG~METABOLITE matrix")

  #Build expression table
  common.genes <- intersect(rownames(nci60.exp), rownames(all.tcga.exp))
  nci.exp <- nci60.exp[common.genes, ]
  nci.exp <- data.table(t(nci.exp), keep.rownames = T)
  setnames(nci.exp, c("CELL", colnames(nci.exp)[2:ncol(nci.exp)]))

  #Is this for a classification problem
  if (classify==T){
    main.act <- nci.gi50[,c("CELL", "NSC", "SCALE.ACT"),with=F]
    main.act$CLASS <- ifelse(main.act$SCALE.ACT>=0.8, 1,
                             ifelse(main.act$SCALE.ACT<=-0.8, 0, 2))
    main.act <- main.act[CLASS!=2,]
    main.act$SCALE.ACT <- main.act$CLASS
  } else{
    main.act <- nci.gi50[,c("CELL", "NSC", "SCALE.ACT"),with=F]
  }
  print ("chose classfication vs regression")

  #Merge all with expression
  common.nsc <- intersect(unique(main.act$NSC), unique(x$NSC))
  common.cell <- intersect(unique(main.act$CELL), unique(nci.exp$CELL))
  main.act <- main.act[NSC %in% common.nsc, ][CELL %in% common.cell,]
  x <- x[NSC %in% common.nsc,]
  nci.exp <- nci.exp[CELL %in% common.cell,]
  print (dim(x))
  print (dim(nci.exp))

  main.table <- inner_join(main.act[,c("NSC", "CELL", "SCALE.ACT"),with=F], x, by="NSC", allow.cartesian = T)
  #main.table <- Function.fast.merging(list(main.act[,c("NSC", "CELL", "SCALE.ACT"),with=F],
  #                                         x), nparts=10, BY="NSC")
  print (dim(main.table))
  print (common.cell)
  main.list <- lapply(common.cell, function(cell){
    print (cell)
    temp.table <- merge(main.table[CELL==cell,], nci.exp[CELL==cell,], by="CELL",  allow.cartesian=T)
    print (object.size(temp.table)/1024^3)
    #saveRDS(temp.table, paste0(OBJ_FOLDER, "nci60_all_feat_table.", cell, ".rds"))
    return(temp.table)
  })
  main.table <- do.call(rbind, main.list)
  rm(main.list)
  print (dim(main.table))
  #main.table <- inner_join(main.table, nci.exp, by="CELL", allow.cartersian = T)
  #main.table <- Function.fast.merging(list(main.table, nci.exp), nparts=10, BY="CELL")
  saveRDS(main.table, paste0(OBJ_FOLDER, "nci60_all_feat_table.rds"))
  print ("done merging all tables")

  #Split into training and validation
  all.rows <- 1:nrow(main.table)
  train.rows <- sample(all.rows, length(all.rows)*split)
  valid.rows <- setdiff(all.rows, train.rows)
  train.table <- main.table[train.rows,]
  valid.table <- main.table[valid.rows,]
  print ("done splitting tables")

  #Do PCA on train.table
  main.pca <- train.table[,-c("CELL", "NSC", "SCALE.ACT"),with=F]
  sd.cols <- colnames(main.pca)[apply(main.pca, 2, sd)!=0]  #USE same columns on valid.table!!!
  main.pca <- main.pca[,sd.cols,with=F]
  main.prim <- prcomp(main.pca, scale. = T)

  #Obtain feature matrix with PCAs that have features which explain 98% of the variance
  pr_var <- main.prim$sdev^2
  pr_var <- pr_var/sum(pr_var)
  plot(pr_var[1:300], type="b")

  needed_pcas <- Function.pca.sd.exp(pr_var, filter=0.98)
  filtered_pca <- main.prim$x[,1:length(needed_pcas)]

  filtered_pca <- cbind(train.table[,c("CELL", "NSC", "SCALE.ACT"),with=F],
                        data.table(filtered_pca))

  print ("done with PCA")
  print (dim(filtered_pca))

  #To obtain a new set using rotation then do:
  # scale(as.matrix(data)) %*% rotation
  #Obtain validation table with PCA features
  valid_pca <- valid.table[ , sd.cols ,with=F]
  valid_pca <- scale(as.matrix(valid_pca)) %*% main.prim$rotation
  valid_pca <- cbind(valid.table[,c("CELL", "NSC", "SCALE.ACT"), with=F],
                     data.table(valid_pca))
  print ("done with validation table")
  print (dim(valid_pca))

  return(list(TRAIN.DATA = filtered_pca,
              VALID.DATA = valid_pca,
              LOADINGS = main.prim$rotation,
              needed_cols = c("CELL", "NSC", "SCALE.ACT", sd.cols)))
}

Function.pca.sd.exp <- function(pr, filter=0.95){

  sd.cum <- c()
  i <- 1
  while(sum(sd.cum)<filter){
    sd.cum <- c(sd.cum, pr[i])
    i <- i + 1
  }

  return(sd.cum)
}

Function.pca.tcga.data <- function(tcga.data, loadings, needed_cols){

  #Filter data
  tcga.data <- tcga.data[,c("SAMPLE","DRUG", "LIVED", needed_cols[4:length(needed_cols)]),with=F]
  tcga.matrix <- tcga.data[,4:ncol(tcga.data),with=F]

  tcga.pca <- scale(as.matrix(tcga.matrix)) %*% loadings
  tcga.pca <- cbind(tcga.data[,c("SAMPLE", "DRUG", "LIVED"),with=F],
                    data.table(tcga.pca))

  #Return
  return(tcga.pca)
}


#################################READ FILES#####################################

OBJ_FOLDER <- "/home/zamalloa/Documents/FOLDER/OBJECTS/" #Lab
TABLES_FOLDER <- "/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/" #Lab
OBJ_FOLDER <- "/tigress/zamalloa/OBJECTS/" #tigress
TABLES_FOLDER <- "/tigress/zamalloa/TABLES/TCGA.TRAINING/"#tigress

nci60.hmdb.file <- paste0(TABLES_FOLDER, "NCI60.TC.HMDB.FP4")
hmdb <- readRDS(paste0(OBJ_FOLDER, "061916.DRUGS.MET.PROFILE.rds"))
nci.gi50 <- readRDS(paste0(OBJ_FOLDER, "061916.NCI.GI50.rds"))
nci60.exp <- readRDS(paste0(OBJ_FOLDER, "061916.NCI60.EXP.rds"))
all.tcga.exp <- readRDS(paste0(OBJ_FOLDER, "062316.ALL.TCGA.EXP.rds"))
tcga.drug.feat.all <- readRDS(paste0(OBJ_FOLDER, "062416.TCGA.DRUG.FEAT.ALL.rds"))

#################################EXECUTE########################################
feat.obj <- Function.nci60.pca.feat(nci60.hmdb.file, hmdb, nci.gi50, nci60.exp, all.tcga.exp, classify = T, split=0.9, OBJ_FOLDER=OBJ_FOLDER)

#################################WRITE TABLES###################################
file.name <- paste0("CGP.AUC_DRUG.TH_", 0, "_", "ALL", "_TARGET", ".", 1)
file.name <- paste0(TABLES_FOLDER, file.name)

write.table(feat.obj$TRAIN.DATA, paste0(file.name, ".TRAIN", ".1"), sep="\t", quote=F, row.names=F, col.names=T)
write.table(feat.obj$VALID.DATA, paste0(file.name, ".VALID", ".1"), sep="\t", quote=F, row.names=F, col.names=T)
saveRDS(feat.obj$LOADINGS, paste0(file.name, ".LOADINGS", ".1"))
saveRDS(feat.obj$needed_cols, paste0(file.name, ".needed_cols"))

#PCA tcga data and write it
tcga_pca <- Function.pca.tcga.data(tcga.drug.feat.all, feat.obj$LOADINGS, feat.obj$needed_cols)

write.table(tcga_pca, paste0(file.name, ".TEST", ".1", sep="\t", quote=F, row.names=F, col.names=T))
