#cgp_new_prep_mf.R

library(data.table)
library(reshape2)

# FUNCTIONS
Function_top_cell_morgan_bits_features_extracted_mf <- function(feat_table, exp_table, morgan_table, max_cells=10, max_bits=10, scaled=T,
                                                                class=F, genes=F, pca=F, lower_th=0.4, higher_th=0.4){

  # Extract most variable cell features
  if (pca==F){
    print("pca==F")

    if (genes==F){
      print("genes==F")

      cell_feat <- cor(exp_table, method = "spearman")
      cell_var  <- data.table(cells = colnames(cell_feat),
                              VAR   = apply(cell_feat, 2, var))
      top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

      cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)

    } else {
      print("genes==T")

      exp_table <- t(exp_table)
      gene_var  <- data.table(genes = colnames(exp_table),
                              VAR   = apply(exp_table, 2, var))
      top_genes <- gene_var[order(-VAR),]$genes[1:max_cells]

      cell_feat <- data.table(exp_table[,top_genes], keep.rownames = T)
    }

  } else{
    print("pca==T")

    if (genes==F){
      print("genes==F")

      cell_feat   <- cor(exp_table, method = "spearman")
      cell_feat   <- prcomp(cell_feat, center = T, scale. = T)
      cell_feat   <- cell_feat$x[,1:max_cells]
      cell_feat   <- data.table(cell_feat, keep.rownames = T)

    } else {
      print("genes==T")

      ctrp_pca    <- prcomp(t(exp_table), center=T, scale. = T) #NOTE: Scaling previously shown to be essential for out of set accuracy
      cell_feat   <- ctrp_pca$x[,1:max_cells]
      cell_feat   <- data.table(cell_feat, keep.rownames = T)
      print("dim cell feat:")
      print(dim(cell_feat))
    }
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features
  if (max_bits > 0){
    if (pca == F){

      top_bits  <- morgan_table[,list(VAR = var(value)), by="bit_pos"][order(-VAR)]$bit_pos[1:max_bits]
      drug_feat <- morgan_table[bit_pos %in% top_bits, ]
      drug_feat <- acast(drug_feat, Compound~bit_pos, value.var="value")
      drug_feat <- data.table(drug_feat, keep.rownames = T)

    } else {
      drug_feat <- acast(morgan_table, Compound~bit_pos, value.var="value")
      drug_pca  <- prcomp(drug_feat, center=F, scale. = F) #Bits, being binary don't need to be scaled (Assumption)
      drug_feat <- data.table(drug_pca$x[, 1:max_bits], keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  } else {
    drug_feat <- data.table()
  }

  # Construct feature table based on classification/regression
  if(class==F){

    if(scaled==T){
      feat_table <- feat_table[,list(cell_name=cell_name, target=scale(target)), by="Compound"]

    } else {
      feat_table <- feat_table
    }

  } else{
    feat_table <- feat_table[,list(cell_name=cell_name, target=scale(target)), by="Compound"] # Need to scale for binary classification
    feat_table[,med:=median(target), by="Compound"]
    feat_table[,sd:= sd(target), by="Compound"]

    feat_table$target <- ifelse(feat_table$target < (feat_table$med - lower_th*abs(feat_table$sd)), 1,
                                ifelse(feat_table$target >  (feat_table$med + higher_th*abs(feat_table$sd)), 0, 2)
                                )
    feat_table <- feat_table[target!=2,]
    feat_table <- feat[,-c("med", "sd"),with=F]
    print("Class average - 1:Sensitive, 0:Resistant")
    print(mean(feat_table$target))
  }

  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  # Merge to obtain total number of possible combinations
  if (max_bits > 0){
    feat_table <- merge(feat_table, drug_feat[, 1:2, with=F], by="Compound", allow.cartesian=TRUE)
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name", allow.cartesian=TRUE)

  setkey(feat_table) #NOTE, may need to remove when "all_rebalance"
  feat_table <- unique(feat_table) #NOTE, may need to remove when "all_rebalance"
  print(paste0("actual total number of samples ", nrow(feat_table)))

  # Extract indices for combinations
  if (max_bits > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else {
    drug_index <- c()
  }
  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$NORM_pIC50

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))
}

Function_top_cell_morgan_counts_features_extracted_mf <- function(feat_table, exp_table, morgan_table, max_cells=10, max_counts=10, scaled=T,
                                                                  class=F, genes=F, pca=F, lower_th=0.4, higher_th=0.4){

  # Extract most variable cell features
  if (pca==F){
    print("pca==F")

    if (genes==F){
      print("genes==F")

      cell_feat <- cor(exp_table, method = "spearman")
      cell_var  <- data.table(cells = colnames(cell_feat),
                              VAR   = apply(cell_feat, 2, var))
      top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

      cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)

    } else {
      print("genes==T")

      exp_table <- t(cgp_exp)
      gene_var  <- data.table(genes = colnames(cgp_exp),
                              VAR   = apply(cgp_exp, 2, var))
      top_genes <- gene_var[order(-VAR),]$genes[1:max_cells]

      cell_feat <- data.table(exp_table[,top_genes], keep.rownames = T)
    }

  } else {
    print("pca==T")

    if (genes==F){
      print("genes==F")

      cell_feat   <- cor(exp_table, method = "spearman")
      cell_feat   <- prcomp(cell_feat, center = T, scale. = T)
      cell_feat   <- cell_feat$x[,1:max_cells]
      cell_feat   <- data.table(cell_feat, keep.rownames = T)

    } else {
      print("genes==T")

      ctrp_pca    <- prcomp(t(exp_table), center=T, scale. = T) #NOTE: Scaling previously shown to be essential for out of set accuracy
      cell_feat   <- ctrp_pca$x[,1:max_cells]
      cell_feat   <- data.table(cell_feat, keep.rownames = T)
      print("dim cell feat:")
      print(dim(cell_feat))
    }
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features (if needed)
  if (max_counts > 0 ){
    if (pca==F) {

      top_counts <- apply(acast(morgan_table, Compound~Substructure, value.var="Counts", fill=0), 2, var)
      top_counts <- names(sort(top_counts, decreasing = T))[1:max_counts]

      drug_feat  <- morgan_table[Substructure %in% top_counts, ]
      drug_feat  <- acast(drug_feat, Compound~Substructure, value.var="Counts", fill=0)
      drug_feat  <- data.table(drug_feat, keep.rownames = T)

    } else {

      drug_feat  <- acast(morgan_table, Compound~Substructure, value.var="Counts", fill=0)
      drug_pca   <- prcomp(drug_feat, center=F, scale. = F) #NOTE: Not scaling needed for discrete variables (Assumed)
      drug_feat  <- data.table(drug_pca$x[, 1:max_counts], keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  } else {
    drug_feat  <- data.table()
  }

  # Construct feature table based on classification/regression
  if(class==F){

    if(scaled==T){
      feat_table <- feat_table[,list(cell_name=cell_name, target=scale(target)), by="Compound"]

    } else {
      feat_table <- feat_table
    }

  } else{

    feat_table <- feat_table[,list(cell_name=cell_name, target=scale(target)), by="Compound"] # Need to scale for binary classification
    feat_table[,med:=median(target), by="Compound"]
    feat_table[,sd:= sd(target), by="Compound"]

    feat_table$target <- ifelse(feat_table$target < (feat_table$med - lower_th*abs(feat_table$sd)), 1,
                                ifelse(feat_table$target >  (feat_table$med + higher_th*abs(feat_table$sd)), 0, 2)
                                )
    feat_table <- feat_table[target!=2,]
    feat_table <- feat[,-c("med", "sd"),with=F]

    print("Class average - 1:Sensitive, 0:Resistant")
    print(mean(feat_table$target))
  }

  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  # Merge to obtain total number of possible combinations
  if (max_counts > 0){
    feat_table <- merge(feat_table, drug_feat[, 1:2, with=F], by="Compound")
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name")

  setkey(feat_table) #NOTE: May need to remove when doing rebalancing
  feat_table <- unique(feat_table) #NOTE: May need to remove when doing rebalancing

  # Extract indices for combinations
  if (max_counts > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else{
    drug_index <- c()
  }

  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$NORM_pIC50

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))

}

Function.range.0.1 <- function(x){
  scale.1 <-  (x-min(x))/(max(x)-min(x))
  return(scale.1)
}

Function_load_morgan_bits <- function(morgan=T){
  if (morgan==T){
    morgan_nci_bits   <- fread(paste0(in_morgan, "NCI_MORGAN_BITS_r_16_b_2048.txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==16 & bits==2048,]
    morgan_nci_bits$position <- paste0("mcf_", morgan_nci_bits$position)
  } else{
    morgan_nci_bits   <- c()
  }
  return(morgan_nci_bits)
}

Function_load_morgan_counts <- function(morgan=T){
  if (morgan==T){
    morgan_nci_counts <- fread(paste0(in_morgan, "NCI_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "numeric", "numeric"))[radius==12,]
  } else{
    morgan_nci_counts <- c()
  }
}

Function_gene_target <- function(gene_target, cgp_exp){

  gee_folder <- "/tigress/zamalloa/OBJECTS/GEELEHER/"

  if (gene_target=="ccle"){
    print("gene_target: ccle")
    selected_genes <- intersect(rownames(cgp_exp), rownames(ccle_exp))

  } else if (gene_target=="geeleher_cisplatin"){
    print("gene_target: geeleher_cisplatin")
    cis            <- readRDS(paste0(gee_folder, "030217_GEE_CISPLATIN.rds"))
    selected_genes <- intersect(rownames(cgp_exp), rownames(cis$exp_table))

  } else if (gene_target=="geeleher_docetaxel"){
    print("gene_target: geeleher_docetaxel")
    doc            <- readRDS(paste0(gee_folder, "030217_GEE_DOCETAXEL.rds"))
    selected_genes <- intersect(rownames(cgp_exp), rownames(doc$exp_table))

  } else if (gene_target=="geeleher_bortezomib_a"){
    print("gene_target: geeleher_bortezomib_a")
    bor            <- readRDS(paste0(gee_folder, "030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
    selected_genes <- intersect(rownames(cgp_exp), rownames(bor$exp_table_a))

  } else if (gene_target=="geeleher_bortezomib_b"){
    print("gene_target: geeleher_bortezomib_b")
    bor            <- readRDS(paste0(gee_folder, "030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
    selected_genes <- intersect(rownames(cgp_exp), rownames(bor$exp_table_b))

  } else if (gene_target=="geeleher_erlotinib"){
    print("gene_target: geeleher_erlotinib")
    erl            <- readRDS(paste0(gee_folder, "030417_GEE_ERLOTINIB.rds"))
    selected_genes <- intersect(rownames(cgp_exp), rownames(erl$exp_table))

  } else {
    print("gene_target: None")
    selected_genes <- rownames(cgp_exp)
  }

  return(cgp_exp[selected_genes,])
}

Function_write_tables <- function(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder, file_name, scale_tables=F){

  train_drug       <- unique(names(feat_table$drug_index[train_rows]))
  valid_drug       <- unique(names(feat_table$drug_index[valid_rows]))
  test_drug        <- unique(names(feat_table$drug_index[test_rows]))

  train_cell       <- unique(names(feat_table$cell_index[train_rows]))
  valid_cell       <- unique(names(feat_table$cell_index[valid_rows]))
  test_cell        <- unique(names(feat_table$cell_index[test_rows]))

  if (max_drugs > 0){
    train_drug_table <- feat_table$drug_feat[Compound %in% train_drug,]
    valid_drug_table <- feat_table$drug_feat[Compound %in% valid_drug,]
    test_drug_table  <- feat_table$drug_feat[Compound %in% test_drug,]
  } else{
    train_drug_table <- data.table()
    valid_drug_table <- data.table()
    test_drug_table  <- data.table()
  }

  train_cell_table <- feat_table$cell_feat[cell_name %in% train_cell,]
  valid_cell_table <- feat_table$cell_feat[cell_name %in% valid_cell,]
  test_cell_table  <- feat_table$cell_feat[cell_name %in% test_cell,]

  # Obtain new index
  train_feat_table <- feat_table$feat_table[train_rows,]
  valid_feat_table <- feat_table$feat_table[valid_rows,]
  test_feat_table  <- feat_table$feat_table[test_rows,]

  if (max_drugs > 0){
    train_drug_index <- unlist(sapply(train_feat_table$Compound, function(x) which(x==train_drug_table$Compound)))
    valid_drug_index <- unlist(sapply(valid_feat_table$Compound, function(x) which(x==valid_drug_table$Compound)))
    test_drug_index  <- unlist(sapply(test_feat_table$Compound,  function(x) which(x==test_drug_table$Compound)))
  } else {
    train_drug_index <- 0
    valid_drug_index <- 0
    test_drug_index  <- 0
  }

  train_cell_index <- unlist(sapply(train_feat_table$cell_name, function(x) which(x==train_cell_table$cell_name)))
  valid_cell_index <- unlist(sapply(valid_feat_table$cell_name, function(x) which(x==valid_cell_table$cell_name)))
  test_cell_index  <- unlist(sapply(test_feat_table$cell_name,  function(x) which(x==test_cell_table$cell_name)))

  train_target     <- train_feat_table$NORM_pIC50
  valid_target     <- valid_feat_table$NORM_pIC50
  test_target      <- test_feat_table$NORM_pIC50

  train_index      <- data.table(drug = train_drug_index - 1,
                                 cell = train_cell_index - 1,
                                 NORM_pIC50 = train_target)

  valid_index      <- data.table(drug = valid_drug_index - 1,
                                 cell = valid_cell_index - 1,
                                 NORM_pIC50 = valid_target)

  test_index       <- data.table(drug = test_drug_index - 1,
                                 cell = test_cell_index - 1,
                                 NORM_pIC50 = test_target)

  # Scale tables with respect to training set (for Multiplicative_fusion, only scale drug features if they are continous variables)
  if (met_type=="drug_cor"){
    not_scaling <- 1:3
    scaling     <- (max(not_scaling) + 1):ncol(feat_table)

  } else if (met_type=="morgan_bits"){
    not_scaling <- 1
    scaling     <- 2:ncol(train_cell_table)

  } else if (met_type=="morgan_counts"){
    not_scaling <- 1
    scaling     <- 2:ncol(train_cell_table)
  }

  if (scale_tables==T){
    print("scaling training cell tables")
    cell_scaling     <- Function_scaling_tables(train_cell_table, valid_cell_table, test_cell_table,
                                                scaling, not_scaling)

    train_cell_table <- cell_scaling[["train_table"]]
    valid_cell_table <- cell_scaling[["valid_table"]]
    test_cell_table  <- cell_scaling[["test_table"]]

    if (pca==T & max_drugs>0){
      print("scaling training drug tables")
      drug_not_scaling <- 1
      drug_scaling     <- 2:ncol(train_drug_table)

      drug_scaling     <- Function_scaling_tables(train_drug_table, valid_drug_table, test_drug_table,
                                                  drug_scaling, drug_not_scaling)

      train_drug_table <- drug_scaling[["train_table"]]
      valid_drug_table <- drug_scaling[["valid_table"]]
      test_drug_table  <- drug_scaling[["test_table"]]
    }
  }

  # WRITE TABLES
  write.table(train_index, paste0(out_folder, file_name, "_train_index"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(valid_index, paste0(out_folder, file_name, "_valid_index"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(test_index,  paste0(out_folder, file_name, "_test_index"),
              quote=F, sep="\t", row.names=F, col.names=T)

  write.table(train_drug_table, paste0(out_folder, file_name, "_train_drug"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(valid_drug_table, paste0(out_folder, file_name, "_valid_drug"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(test_drug_table, paste0(out_folder, file_name,  "_test_drug"),
              quote=F, sep="\t", row.names=F, col.names=T)

  write.table(train_cell_table, paste0(out_folder, file_name, "_train_cell"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(valid_cell_table, paste0(out_folder, file_name, "_valid_cell"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(test_cell_table, paste0(out_folder, file_name,  "_test_cell"),
              quote=F, sep="\t", row.names=F, col.names=T)

  print("Done writing sel tables")

}

Function_write_train_tables <- function(feat_table, train_rows, max_drugs, pca, out_folder, file_name, scale_tables=F){

  train_drug       <- unique(names(feat_table$drug_index[train_rows]))
  train_cell       <- unique(names(feat_table$cell_index[train_rows]))

  if (max_drugs > 0){
    train_drug_table <- feat_table$drug_feat[Compound %in% train_drug,]
  } else{
    train_drug_table <- data.table()
  }

  train_cell_table <- feat_table$cell_feat[cell_name %in% train_cell,]

  # Obtain new index
  train_feat_table <- feat_table$feat_table[train_rows,]

  if (max_drugs > 0){
    train_drug_index <- unlist(sapply(train_feat_table$Compound, function(x) which(x==train_drug_table$Compound)))
  } else {
    train_drug_index <- 0
  }

  train_cell_index <- unlist(sapply(train_feat_table$cell_name, function(x) which(x==train_cell_table$cell_name)))
  train_target     <- train_feat_table$NORM_pIC50

  train_index      <- data.table(drug = train_drug_index - 1,
                                 cell = train_cell_index - 1,
                                 NORM_pIC50 = train_target)

  # Do we need to scale them?
  if (met_type=="drug_cor"){
    not_scaling <- 1:3
    scaling     <- (max(not_scaling) + 1):ncol(feat_table)

  } else if (met_type=="morgan_bits"){
    not_scaling <- 1
    scaling     <- 2:ncol(train_cell_table)

  } else if (met_type=="morgan_counts"){
    not_scaling <- 1
    scaling     <- 2:ncol(train_cell_table)
  }

  if (scale_tables==T){
    print("scaling training cell tables")
    scale_train      <- scale(train_cell_table[, scaling, with=F])
    train_cell_table <- cbind(train_cell_table[, not_scaling, with=F], scale_train)


    if (pca==T & max_drugs>0){
      print("scaling training drug tables")
      drug_not_scaling <- 1
      drug_scaling     <- 2:ncol(train_drug_table)

      scale_train      <- scale(train_drug_table[,drug_scaling,with=F])
      train_drug_table <- cbind(train_drug_table[,drug_not_scaling,with=F], scale_train)
    }
  }


  # WRITE TABLES
  write.table(train_index, paste0(out_folder, file_name, "_train_index"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(train_drug_table, paste0(out_folder, file_name, "_train_drug"),
              quote=F, sep="\t", row.names=F, col.names=T)
  write.table(train_cell_table, paste0(out_folder, file_name, "_train_cell"),
              quote=F, sep="\t", row.names=F, col.names=T)

  print("Done writing sel train only tables")
}

Function_scaling_tables <- function(train_table, valid_table, test_table, scaling, not_scaling){

  scale_train <- scale(train_table[, scaling, with=F])
  train_table <- cbind(train_table[, not_scaling, with=F], scale_train)

  train_mean  <- attributes(scale_train)$`scaled:center`
  train_sd    <- attributes(scale_train)$`scaled:scale`

  valid_scale <- sweep(valid_table[, scaling, with=F], 2, train_mean, "-")
  valid_scale <- sweep(valid_scale, 2, train_sd, "/")
  valid_table <- cbind(valid_table[, not_scaling, with=F],
                       valid_scale)

  test_scale  <- sweep(test_table[, scaling, with=F], 2, train_mean, "-")
  test_scale  <- sweep(test_scale, 2, train_sd, "/")
  test_table  <- cbind(test_table[, not_scaling, with=F],
                       test_scale)

  return(list(train_table = train_table,
              valid_table = valid_table,
              test_table  = test_table))

}

Function_rebalance <- function(cgp_table, th){
  # Rebalance compound count in term of number of samples and activity

  # First rebalance by activity
  print("balancing by activity")
  cgp_table     <- cgp_table[,-c("Site"),with=F] #To reduce overhead
  all_compounds <- unique(cgp_table$Compound)
  reb_1         <- data.table()

  for (drug in all_compounds){
    temp_table  <- cgp_table[Compound==drug,]

    above_table <- temp_table[((1-AUC)*10)>=th,]
    below_table <- temp_table[((1-AUC)*10)<th,]

    if (nrow(above_table)>0){

      if (nrow(below_table) > nrow(above_table)){
        needed_rows <- nrow(below_table) - nrow(above_table)
        set.seed(1234)
        above_table <- rbind(above_table,
                             above_table[sample(1:nrow(above_table), needed_rows, replace=T),])
      } else{
        needed_rows <- nrow(above_table) - nrow(below_table)
        set.seed(1234)
        below_table <- rbind(below_table,
                             below_table[sample(1:nrow(below_table), needed_rows, replace=T),])
      }
      reb_1       <- rbind(reb_1,
                           rbind(below_table, above_table))
    } else{
      reb_1       <- rbind(reb_1, below_table)
    }
  }

  # Then rebalance by count
  print("balancing by count")
  drug_count    <- reb_1[,list(N=length(cell_name)), by="Compound"]
  max_count     <- max(drug_count$N)

  reb_2         <- data.table()

  for (drug in all_compounds){
    temp_table  <- reb_1[Compound==drug,]

    if (nrow(temp_table) < max_count){
      needed_count    <- max_count - nrow(temp_table)
      set.seed(1234)
      temp_table      <- rbind(temp_table,
                               temp_table[sample(1:nrow(temp_table), needed_count, replace=T),])
    }

    reb_2       <- rbind(reb_2, temp_table)
  }

  # Return
  print(paste0("Approximate total number of samples = ", nrow(reb_2)))
  return(reb_2)
}

# LOAD INPUT
args        <- commandArgs(trailingOnly = TRUE)

max_cells   <- as.numeric(args[1])
max_drugs   <- as.numeric(args[2])
file_name   <- args[3]
met_type    <- args[4]
class_mlp   <- as.logical(args[5])
samples     <- args[6]
genes       <- as.logical(args[7])
batch_norm  <- args[8]
pca         <- as.logical(args[9])
rebalance   <- as.logical(args[10])
radii_set   <- as.numeric(args[11])
bit_set     <- as.numeric(args[12])
gene_target <- args[13]
fold        <- args[14]
lower_th    <- as.numeric(args[15])
higher_th   <- as.numeric(args[16])

in_folder    <- "/tigress/zamalloa/CGP_FILES/" #For tigress
in_morgan    <- "/tigress/zamalloa/MORGAN_FILES/"
in_objects   <- "/tigress/zamalloa/OBJECTS/"
tcga_objects <- "/tigress/zamalloa/TCGA_FILES/"
out_folder   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/"

# LOAD DATA
DRUGS.MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
ctrp_new          <- readRDS(paste0(in_objects, "ctrp_tables.rds"))[["ctrp_2.0"]][,c("cpd_name","ccl_name","area_under_curve"),with=F] #More data points in 2.0
setnames(ctrp_new, c("Compound", "cell_name", "target"))

if (batch_norm=="ctrp_cgp") {
  print("ctrp_cgp")
  ctrp_exp        <- readRDS(paste0(in_objects, "041817_cgp_ctrp_exp_norm.rds"))[["EXP.2"]]

} else{
  print("None")
  ctrp_exp        <- readRDS(paste0(in_object, "ctrp_exp_2.1.rds")) #But expression information found in version 2.1, total of ~370K
}

ctrp_bits         <- fread(paste0(in_morgan, "CTRP_MORGAN_BITS_r_",radii_set,"_b_",bit_set,".txt"))
setnames(ctrp_bits, c("radius", "bits", "Compound", "bit_pos", "value"))
ctrp_bits$bit_pos <- paste0("mcf_", ctrp_bits$bit_pos)

# morgan_counts     <- fread(paste0(in_morgan, "CGP_MORGAN_COUNTS.txt"),
#                        colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,] #12 normal setting

#################################################### EXECUTE ####################################################
# Apply
if (met_type=="drug_cor"){
  feat_table <- Function_top_cell_drug_features_extracted_mf(ctrp_new, ctrp_exp, DRUGS.MET.PROFILE,
                                                          max_cells = max_cells, max_drugs = max_drugs,
                                                          class = class_mlp, scaled = T)

  drug_sim   <- data.table(melt(cor(acast(DRUGS.MET.PROFILE, METABOLITE~DRUG, value.var = "TC"))))
  drug_sim   <- drug_sim[Var1!=Var2,]

} else if (met_type=="morgan_bits"){
  feat_table <- Function_top_cell_morgan_bits_features_extracted_mf(ctrp_new, ctrp_exp, morgan_bits,
                                                          max_cells = max_cells, max_bits = max_drugs,
                                                          class = class_mlp, scaled = F, #MODIFIED
                                                          genes = genes, pca = pca, lower_th, higher_th)

  drug_sim   <- data.table(melt(cor(acast(morgan_bits, bit_pos~Compound, value.var = "value"))))
  drug_sim   <- drug_sim[Var1!=Var2,]

} else if (met_type=="morgan_counts"){
  feat_table <- Function_top_cell_morgan_counts_features_extracted_mf(ctrp_new, ctrp_exp, morgan_counts,
                                                          max_cells = max_cells, max_counts = max_drugs,
                                                          class = class_mlp, scaled = F, #MODIFIED
                                                          genes = genes, pca = pca, lower_th, higher_th)

  drug_sim   <- data.table(melt(cor(acast(morgan_counts, Substructure~Compound, value.var = "Counts", fill = 0))))
  drug_sim   <- drug_sim[Var1!=Var2,]

}

# Split tables depending on all or drug samples
#UPDATE: Take folding into account
if ( (samples == "all") | (grepl("all_rebalance_", samples)==T)){
  #NOTE: If "all_rebalance_", this type of rebalance has been already applied prior to building feat_table and indexes

  temp_rows          <- 1:length(feat_table$target)

  if (fold=="fold_early"){
    # NOTE:Train with all data (For the time being, use early stopping of 25%)
    set.seed(1234)
    train_rows       <- sample(temp_rows, length(temp_rows)*0.75)
    valid_rows       <- setdiff(temp_rows, train_rows) #For early stopping
    test_rows        <- valid_rows  #Same as valid

    Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                          file_name, scale_tables=F) #NOTE:CHANGED SCALING TEMPORARILY!!!

  } else if(fold=="fold_none"){
    # NOTE:Train with all data, literally, no split done for early stopping
    train_rows       <- 1:length(feat_table$target)

    Function_write_train_tables(feat_table, train_rows, max_drugs, pca, out_folder,
                                file_name, scale_tables=F) #NOTE:CHANGED SCALING TEMPORARILY!!!

  } else if (fold=="fold_all"){
    for (split_fold in 1:5){
      print(fold)
      print(split_fold)
      splits           <- split(temp_rows, 1:5)
      test_rows        <- splits[[split_fold]]

      testing_rows     <- setdiff(temp_rows, test_rows)
      set.seed(1234)
      train_rows       <- sample(testing_rows, length(testing_rows)*0.75)
      valid_rows       <- setdiff(testing_rows, train_rows)

      print(file_name)
      print(gsub("fold_all", split_fold, file_name))
      Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                            gsub("fold_all", split_fold, file_name), scale_tables=F)
    }
  }

} else if ( (grepl("act_rebalance_", samples)==T) & (grepl("zero_", samples)==F) ){

  th_rebalance <- as.numeric(gsub("act_rebalance_", "", samples))

  if (fold=="fold_early"){
    # NOTE:Train with all data (For the time being, use early stopping of 25%)

    extra_rows   <- which(feat_table$target >= th_rebalance)
    non_rows     <- which(feat_table$target < th_rebalance)

    set.seed(1234)
    extra_rows   <- sample(extra_rows, length(non_rows), replace=T)

    all_rows     <- c(non_rows, extra_rows)

    set.seed(1234)
    train_rows   <- sample(all_rows, length(all_rows)*0.75)
    valid_rows   <- setdiff(all_rows, train_rows)
    test_rows    <- valid_rows  #Same as valid

    Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                          file_name, scale_tables=F) #NOTE:CHANGED SCALING TEMPORARILY!!!

  } else if(fold=="fold_none"){
    # NOTE:Train with all data, literally, no split done for early stopping
    extra_rows   <- which(feat_table$target >= th_rebalance)
    non_rows     <- which(feat_table$target < th_rebalance)

    set.seed(1234)
    extra_rows   <- sample(extra_rows, length(non_rows), replace=T)

    train_rows       <- c(non_rows, extra_rows)

    Function_write_train_tables(feat_table, train_rows, max_drugs, pca, out_folder,
                                file_name, scale_tables=F) #NOTE:CHANGED SCALING TEMPORARILY!!!

  }

} else if (grepl("act_rebalancetop_", samples)==T){

  print("act_rebalancetop_")
  th_rebalance <- as.numeric(gsub("act_rebalancetop_", "", samples))

  set.seed(1234)
  test_rows    <- sample(1:length(feat_table$target), length(feat_table$target)*0.2)

  drugs        <- unique(feat_table$feat_table$Compound)

  extra_rows   <- lapply(drugs, function(x) {
    temp_th    <- feat_table$feat_table[Compound==x,][order(-NORM_pIC50),]

    non_cells  <- feat_table$feat_table[test_rows,][Compound==x,]$cell_name

    temp_th    <- temp_th[!cell_name %in% non_cells,]$NORM_pIC50[1:th_rebalance]
    temp_th    <- min(temp_th)

    temp_rows  <- which(feat_table$feat_table$Compound == x & feat_table$feat_table$NORM_pIC50 >=temp_th)
    temp_rows  <- setdiff(temp_rows, test_rows)

    return(temp_rows)
    })

  extra_rows   <- as.vector(unlist(extra_rows))
  non_rows     <- setdiff(1:length(feat_table$target), extra_rows)
  extra_rows   <- setdiff(extra_rows, test_rows)
  non_rows     <- setdiff(non_rows, test_rows)

  extra_rows   <- sample(extra_rows, length(non_rows), replace=T)

  all_rows     <- c(non_rows, extra_rows)

  set.seed(1234)
  train_rows   <- sample(all_rows, length(all_rows)*0.75)

  valid_rows   <- setdiff(all_rows, train_rows)
  print("act_rebalancetop_")

} else if (grepl("percent_all_", samples)==T){

  portion  <- as.numeric(gsub("percent_all_", "", samples))
  portion  <- portion/100
  print(portion)

  set.seed(1234)
  all_rows     <- sample(1:length(feat_table$target), length(feat_table$target)* portion )

  set.seed(1234)
  train_rows   <- sample(all_rows, length(all_rows)*0.7)

  testing_rows <- setdiff(all_rows, train_rows)
  set.seed(1234)
  valid_rows   <- sample(testing_rows, length(testing_rows)*0.5)
  test_rows    <- setdiff(testing_rows, valid_rows)

} else if (samples == "whole"){
  #AKA all data training

  set.seed(1234)
  train_rows       <- sample(1:length(feat_table$target), length(feat_table$target)*0.8)
  valid_rows       <- setdiff(1:length(feat_table$target), train_rows)
  test_rows        <- valid_rows

} else if (samples == "all_split"){
  print("all_split")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  train_drug       <- sample(all_drugs, length(all_drugs)*0.7)
  testing_drugs    <- setdiff(all_drugs, train_drug) #all but training
  set.seed(1234)
  valid_drug       <- sample(testing_drugs, length(testing_drugs)*0.5)
  test_drug        <- setdiff(testing_drugs, valid_drug)

  train_rows       <- which(feat_table$feat_table$Compound %in% train_drug)
  valid_rows       <- which(feat_table$feat_table$Compound %in% valid_drug)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else if (samples == "semi_split"){
  print ("semi_split")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  #test_drug        <- sample(all_drugs, 5)
  #test_drug        <- c("GDC0449", "MS-275", "PAC-1", "RDEA119", "TG101348")
  test_drug        <- sample(all_drugs, 0.2*length(all_drugs))
  # test_drug        <- c("17-AAG", "BHG712", "Bleomycin", "BX-912", "CEP-701", "CH5424802",
  #                       "CMK", "CP724714", "Docetaxel", "EHT 1864", "FR-180204", "FTI-277",
  #                       "GSK1070916", "GSK2126458", "IPA-3", "Ispinesib Mesylate", "JQ1",
  #                       "KU-55933", "Lapatinib", "LAQ824", "Midostaurin", "NSC-87877",
  #                       "NU-7441", "Nutlin-3a (-)", "Parthenolide", "PHA-793887", "Rapamycin",
  #                       "Sunitinib", "TAK-715", "Temozolomide", "Tipifarnib", "Trametinib",
  #                       "Tubastatin A", "YM201636")
  testing_drugs    <- setdiff(all_drugs, test_drug)

  testing_rows     <- which(feat_table$feat_table$Compound %in% testing_drugs)
  set.seed(1234)
  train_rows       <- sample(testing_rows, length(testing_rows)*0.80)
  valid_rows       <- setdiff(testing_rows, train_rows)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else if(samples == "all_split_balanced"){
  print("all_split_balanced")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  test_drug        <- sample(all_drugs, length(all_drugs)*0.15)
  testing_drugs    <- setdiff(all_drugs, test_drug)

  # Sample based on similarity to testing for balanaced validation
  drug_sim         <- drug_sim[Var1 %in% test_drug,][Var2 %in% testing_drugs,]
  drug_sim         <- drug_sim[,list(mean_sim = mean(value)), by="Var2"] # Mean similarity to target
  valid_drug       <- sample(drug_sim$Var2, length(testing_drugs)*0.25 , prob= Function.range.0.1(drug_sim$mean_sim))
  train_drug       <- setdiff(testing_drugs, valid_drug)

  train_rows       <- which(feat_table$feat_table$Compound %in% train_drug)
  valid_rows       <- which(feat_table$feat_table$Compound %in% valid_drug)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else if(grepl("zero_", samples)==T) {
  print ("zero morgan target compound")

  if (grepl("^", samples)==T){
    drug_name <- gsub("^", " ", samples, fixed=T)
    drug_name <- strsplit(drug_name, "zero_")[[1]][2]

  } else{
    drug_name <- strsplit(samples, "zero_")[[1]][2]
  }

  # Uses zero morgan features and splits dataset for specific compound
  temp_rows        <- which(feat_table$feat_table$Compound == drug_name)

  # Check if we need to rebalance that
  if (grepl("act_rebalance_", samples)==T){

    th_rebalance <- strsplit(samples, "_zero")[[1]][1]
    th_rebalance <- as.numeric(gsub("act_rebalance_", "", th_rebalance))

    extra_rows   <- which(feat_table$feat_table$Compound==drug_name &  (feat_table$target >= th_rebalance))
    non_rows     <- which(feat_table$feat_table$Compound==drug_name &  (feat_table$target < th_rebalance))

    if (length(non_rows) > length(extra_rows)){
      set.seed(1234)
      extra_rows   <- sample(extra_rows, length(non_rows), replace=T)
    } else{
      non_rows     <- sample(non_rows, length(extra_rows), replace=T)
    }

    temp_rows    <- c(non_rows, extra_rows) #NOTE: Effectively replacing "temp_rows" indexes obtained above
  }

  if (fold=="fold_early"){
    # NOTE:Train with all data using early stopping

    set.seed(1234)
    train_rows       <- sample(temp_rows, length(temp_rows)*0.75)
    valid_rows       <- setdiff(temp_rows, train_rows)
    test_rows        <- valid_rows # Sample as valid rows for comparisson

    Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                          file_name, scale_tables=F)

  } else if (fold=="fold_none"){
    # NOTE:Train with all data, literally, no split done for early stopping
    train_rows       <- temp_rows

    Function_write_train_tables(feat_table, train_rows, max_drugs, pca, out_folder,
                                file_name, scale_tables=F)

  } else if (fold=="fold_all"){
    for (split_fold in 1:10){
      print(fold)
      print(split_fold)
      splits           <- split(temp_rows, 1:10)
      test_rows        <- splits[[split_fold]]

      testing_rows     <- setdiff(temp_rows, test_rows)
      set.seed(1234)
      train_rows       <- sample(testing_rows, length(testing_rows)*0.66)
      valid_rows       <- setdiff(testing_rows, train_rows)

      print(file_name)
      print(gsub("fold_all", split_fold, file_name))
      Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                            gsub("fold_all", split_fold, file_name), scale_tables=F)
    }
  }

} else if (grepl("zeroall_", samples)==T) {
  print ("zero morgan target compound - all data")
  # Uses zero morgan features and splits dataset for specific compound, trains on all data
  temp_rows        <- which(feat_table$feat_table$Compound == gsub("zeroall_", "", samples))

  set.seed(1234)
  train_rows       <- sample(temp_rows, length(temp_rows)*0.66)
  valid_rows       <- setdiff(temp_rows, train_rows)
  test_rows        <- valid_rows

  Function_write_tables(feat_table, train_rows, valid_rows, test_rows, max_drugs, pca, out_folder,
                        file_name, scale_tables=F)

} else {

  test_rows        <- which(feat_table$feat_table$Compound == samples)
  testing_rows     <- which(feat_table$feat_table$Compound != samples)

  set.seed(1234)
  train_rows       <- sample(testing_rows, length(testing_rows)*0.8)
  valid_rows       <- setdiff(testing_rows, train_rows)
}
