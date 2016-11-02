# cgp_new_prep.R

# Preps mixed cgp data for overall data prediction
# This script was written initially for classification purposes

library(data.table)
library(reshape2)

# FUNCTIONS
Function_top_cell_drug_features_extracted <- function(feats, exp_table, met_table, max_cells=10, max_drugs=10, met_scaled=F, pic50_scaled=T,
                                                      pic50_class=F){
  # Constructs feature table using drugs and cell correlations as features, but limiting to most variable features
  # Both max_cells and max_drugs are based in terms of decreasing variance
  # NOTE: Function can now be used for classification or regression
  #       If classification is chosen, pic50-based classes are chosen (1/0), while 2 is discarded

  # Extract most variable cell features
  #cell_feat <- cor(t(scale(log(t(exp_table)))), method = "pearson") #Temporary change
  cell_feat <- cor(cgp_exp, method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features - is met scaling necessary?
  if (met_scaled==T){
    met_table[, SD:=sd(TC), by ="METABOLITE"]
    met_table <- met_table[SD!=0,]
    met_table <- met_table[, 1:3, with=F]

    met_table <- met_table[, TC := scale(TC), by="METABOLITE"]
  }

  drug_feat <- cor(acast(met_table, METABOLITE~DRUG, value.var = "TC"), method = "pearson")
  drug_var  <- data.table(drugs = colnames(drug_feat),
                          VAR   = apply(drug_feat, 2, var))
  top_drugs <- drug_var[order(-VAR),]$drugs[1:max_drugs]

  drug_feat <- data.table(drug_feat[, top_drugs], keep.rownames = T)
  setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))

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

  feat_table <- merge(feat_table, cell_feat, by="cell_name")
  feat_table <- merge(feat_table, drug_feat, by="Compound")

  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

# LOAD INPUT
args        <- commandArgs(trailingOnly = TRUE)

max_cells   <- as.numeric(args[1])
max_drugs   <- as.numeric(args[2])
file_name   <- args[3]

# LOAD DATA
in_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress
out_folder  <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/"

DRUGS.MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
cgp_new     <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
cgp_exp     <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))

#################################################### EXECUTE ####################################################
feat_table <- Function_top_cell_drug_features_extracted(cgp_new, cgp_exp, DRUGS.MET.PROFILE,
                                                        max_cells = max_cells, max_drugs = max_drugs,
                                                        pic50_class = T)

# Split tables
train_rows   <- sample(1:nrow(feat_table), nrow(feat_table)*0.7)
testing_rows <- setdiff(1:nrow(feat_table), train_rows)
valid_rows   <- sample(testing_rows, length(testing_rows)*0.5)
test_rows    <- setdiff(testing_rows, valid_rows)

train_table  <- feat_table[train_rows,]
valid_table  <- feat_table[valid_rows,]
test_table   <- feat_table[test_rows,]

# Scale tables with respect to training set
scale_train  <- scale(train_table[, 4:ncol(train_table), with=F])
train_table  <- cbind(train_table[, 1:3, with=F], scale_train)

train_mean   <- attributes(scale_train)$`scaled:center`
train_sd     <- attributes(scale_train)$`scaled:scale`

valid_scale <- sweep(valid_table[, 4:ncol(valid_table), with=F], 2, train_mean, "-")
valid_scale <- sweep(valid_scale, 2, train_sd, "/")
valid_table <- cbind(valid_table[, 1:3, with=F],
                     valid_scale)

test_scale  <- sweep(test_table[, 4:ncol(test_table), with=F], 2, train_mean, "-")
test_scale  <- sweep(test_scale, 2, train_sd, "/")
test_table  <- cbind(test_table[, 1:3, with=F],
                     test_scale)

# Write out
write.table(train_table, paste0(out_folder, file_name, "_train"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(valid_table, paste0(out_folder, file_name, "_valid"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(test_table, paste0(out_folder, file_name, "_test"),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing sel tables")
