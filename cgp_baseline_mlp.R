# cgp_baseline_mlp.R
# Builds training and testing sets for modeling cgp drugs

library(data.table)
library(reshape2)

Function_top_cell_drug_features_extracted <- function(feats, exp_table, met_table, max_cells=10, max_drugs=10){
  # Constructs feature table using drugs and cell correlations as features, but limiting to most variable features
  # Both max_cells and max_drugs are based in terms of decreasing variance

  # Extract most variable cell features
  cell_feat <- cor(t(scale(log(t(exp_table)))), method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features
  drug_feat <- cor(acast(met_table, METABOLITE~DRUG, value.var = "TC"), method = "pearson")
  drug_var  <- data.table(drugs = colnames(drug_feat),
                          VAR   = apply(drug_feat, 2, var))
  top_drugs <- drug_var[order(-VAR),]$drugs[1:max_drugs]

  drug_feat <- data.table(drug_feat[, top_drugs], keep.rownames = T)
  setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))

  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  feat_table <- merge(feat_table, cell_feat, by="cell_name")
  feat_table <- merge(feat_table, drug_feat, by="Compound")

  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  feat_table <- feat_table[sample(nrow(feat_table)),] #Randomizes entries
  return(feat_table)
}

Function_randomize_cell_drug_draw <- function(feat_table, target_drug, random_cells=c(), random_compounds=c()){
  # Function to randomly sample either cell_name or Compound attributes and return filtered down feature table

  target_table <- feat_table[Compound==target_drug,]
  feat_table   <- feat_table[Compound!=target_drug,]

  all_cells    <- unique(feat_table$cell_name)
  all_comp     <- unique(feat_table$Compound)

  # Select random samples
  if (length(random_cells)>0){
    new_cells  <- sample(all_cells, random_cells)
  } else {
    new_cells  <- all_cells
  }

  if (length(random_compounds)>0){
    new_compounds <- sample(all_comp, random_compounds)
  } else {
    new_compounds <- all_comp
  }

  # Obtain new feature table out of random samples
  feat_table <- feat_table[Compound %in% new_compounds,]
  feat_table <- feat_table[cell_name %in% new_cells,]

  # Clean up and Return
  feat_table <- rbind(feat_table, target_table)
  return(feat_table)
}

######################################################################################################
# LOAD DATA
args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
cell_n_feat <- args[2]
drug_n_feat <- args[3]

in_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress
out_table   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/" #For tigress

MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
cgp_new     <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
cgp_exp     <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))

#####################################################################################################
# EXECUTE
feat_table  <- Function_top_cell_drug_features_extracted(cgp_new, cgp_exp, MET.PROFILE,
                                                         max_cells = cell_n_feat, max_drugs = drug_n_feat)

#Split for training 80/20
target      <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ") #To account for drug names such as Mitomycin_C

test_table  <- feat_table[Compound==target,]
temp_table  <- feat_table[Compound!=target,]

all_rows    <- 1:nrow(temp_table)
train_rows  <- sample(all_rows, length(all_rows)*0.8)
valid_rows  <- setdiff(all_rows, train_rows)

train_table <- temp_table[train_rows,]
valid_table <- temp_table[valid_rows,]

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, "TRAIN_",  target_drug, "_C_", cell_n_feat, "_D_", drug_n_feat),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, "VALID_", target_drug, "_C_", cell_n_feat, "_D_", drug_n_feat),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table, paste0(out_table,  "TEST_",  target_drug, "_C_", cell_n_feat, "_D_", drug_n_feat),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing sel tables")
