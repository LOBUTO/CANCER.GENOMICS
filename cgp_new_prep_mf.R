#cgp_new_prep_mf.R

library(data.table)
library(reshape2)

# FUNCTIONS
Function_top_cell_morgan_bits_features_extracted_mf <- function(feats, exp_table, morgan_table, max_cells=10, max_bits=10, pic50_scaled=T,
                                                                pic50_class=F){

  # Extract most variable cell features
  cell_feat <- cor(cgp_exp, method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features
  if (max_bits > 0){
    top_bits  <- morgan_table[,list(VAR = var(value)), by="bit_pos"][order(-VAR)]$bit_pos[1:max_bits]
    drug_feat <- morgan_table[bit_pos %in% top_bits, ]
    drug_feat <- acast(drug_feat, Compound~bit_pos, value.var="value")
    drug_feat <- data.table(drug_feat, keep.rownames = T)

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  }

  # Construct feature table based on classification/regression
  if(pic50_class==F){

    if(pic50_scaled==T){
      feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
    } else {
      feat_table <- feats[, c("Compound", "cell_name", "pIC50"), with=F]
    }

  } else{
    feat_table <- feats[, c("Compound", "cell_name", "pic50_class"),with=F]
    feat_table <- feat_table[pic50_class!=2,]
  }

  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  # Merge to obtain total number of possible combinations
  if (max_bits > 0){
    feat_table <- merge(feat_table, drug_feat[, 1:2, with=F], by="Compound")
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name")
  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Extract indeces for combinations
  drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
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

# LOAD INPUT
args        <- commandArgs(trailingOnly = TRUE)

max_cells   <- as.numeric(args[1])
max_drugs   <- as.numeric(args[2])
file_name   <- args[3]
met_type    <- args[4]
class_mlp   <- as.logical(args[5])
samples     <- args[6]

# LOAD DATA
in_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress
out_folder  <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/"

DRUGS.MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
cgp_new           <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
cgp_exp           <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))
morgan_bits       <- fread("/tigress/zamalloa/MORGAN_FILES/morgan_bits.txt")[radius==16 & bits==2048,]
morgan_counts     <- fread("/tigress/zamalloa/MORGAN_FILES/morgan_counts.txt",
                       colClasses = c("character", "numeric", "character", "numeric", "numeric"))[radius==12,]

#################################################### EXECUTE ####################################################
if (met_type=="drug_cor"){
  feat_table <- Function_top_cell_drug_features_extracted_mf(cgp_new, cgp_exp, DRUGS.MET.PROFILE,
                                                          max_cells = max_cells, max_drugs = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T)
} else if (met_type=="morgan_bits"){
  feat_table <- Function_top_cell_morgan_bits_features_extracted_mf(cgp_new, cgp_exp, morgan_bits,
                                                          max_cells = max_cells, max_bits = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T)

} else if (met_type=="morgan_counts"){
  feat_table <- Function_top_cell_morgan_counts_features_extracted_mf(cgp_new, cgp_exp, morgan_counts,
                                                          max_cells = max_cells, max_counts = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T)
}

# Split tables depending on all or drug samples
if (samples == "all"){

  set.seed(1234)
  train_rows       <- sample(1:length(feat_table$target), length(feat_table$target)*0.7)
  testing_rows     <- setdiff(1:length(feat_table$target), train_rows)
  set.seed(1234)
  valid_rows       <- sample(testing_rows, length(testing_rows)*0.5)
  test_rows        <- setdiff(testing_rows, valid_rows)

  train_drug       <- unique(names(feat_table$drug_index[train_rows]))
  valid_drug       <- unique(names(feat_table$drug_index[valid_rows]))
  test_drug        <- unique(names(feat_table$drug_index[test_rows]))

  train_cell       <- unique(names(feat_table$cell_index[train_rows]))
  valid_cell       <- unique(names(feat_table$cell_index[valid_rows]))
  test_cell        <- unique(names(feat_table$cell_index[test_rows]))

  train_drug_table <- feat_table$drug_feat[Compound %in% train_drug,]
  valid_drug_table <- feat_table$drug_feat[Compound %in% valid_drug,]
  test_drug_table  <- feat_table$drug_feat[Compound %in% test_drug,]

  train_cell_table <- feat_table$cell_feat[Compound %in% train_cell,]
  valid_cell_table <- feat_table$cell_feat[Compound %in% valid_cell,]
  test_cell_table  <- feat_table$cell_feat[Compound %in% test_cell,]

  # Obtain new index
  train_feat_table <- feat_table$feat_table[train_rows,]
  valid_feat_table <- feat_table$feat_table[valid_rows,]
  test_feat_table  <- feat_table$feat_table[test_rows,]

  train_drug_index <- unlist(sapply(train_feat_table$Compound, function(x) which(x==train_drug_table$Compound)))
  valid_drug_index <- unlist(sapply(valid_feat_table$Compound, function(x) which(x==valid_drug_table$Compound)))
  test_drug_index  <- unlist(sapply(test_feat_table$Compound,  function(x) which(x==test_drug_table$Compound)))

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

} else {
  testing_table <- feat_table[Compound!=samples, ]
  test_table    <- feat_table[Compound==samples, ]

  set.seed(1234)
  train_rows    <- sample(1:nrow(testing_table), nrow(testing_table)*0.7)
  valid_rows    <- setdiff(1:nrow(testing_table), train_rows)

  train_table   <- testing_table[train_rows, ]
  valid_table   <- testing_table[valid_rows, ]
}


# Scale tables with respect to training set
if (met_type=="drug_cor"){
  not_scaling <- 1:3
  scaling     <- (max(not_scaling) + 1):ncol(feat_table)

} else if (met_type=="morgan_bits"){
  not_scaling <- 1
  scaling     <- 2:ncol(train_cell_table)

} else if (met_type=="morgan_counts"){
  not_scaling <- 1:(3 + max_drugs)
  scaling     <- (max(not_scaling) + 1):ncol(feat_table)
}

scale_train_cell  <- scale(train_cell_table[, scaling, with=F])
train_cell_table  <- cbind(train_cell_table[, not_scaling, with=F], scale_train)

train_cell_mean   <- attributes(scale_train_cell)$`scaled:center`
train_cell_sd     <- attributes(scale_train_cell)$`scaled:scale`

valid_cell_scale  <- sweep(valid_cell_table[, scaling, with=F], 2, train_cell_mean, "-")
valid_cell_scale  <- sweep(valid_cell_scale, 2, train_cell_sd, "/")
valid_cell_table  <- cbind(valid_cell_table[, not_scaling, with=F],
                     valid_cell_scale)

test_cell_scale  <- sweep(test_cell_table[, scaling, with=F], 2, train_cell_mean, "-")
test_cell_scale  <- sweep(test_cell_scale, 2, train_cell_sd, "/")
test_cell_table  <- cbind(test_cell_table[, not_scaling, with=F],
                     test_cell_scale)

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
