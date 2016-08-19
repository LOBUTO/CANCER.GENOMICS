# cgp_new_sel_feat.R
# Function to split cgp table

library(data.table)
library(reshape2)

######################################################################################################
# LOAD DATA
args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
target_drug <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ")
usage       <- args[2]
modifier    <- args[3]

in_folder   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/" #For lab
in_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress

MET.PROFILE <- readRDS(paste0(in_folder, "080716.DRUG.MET.PROFILE.rds"))
if (usage=="nci60"){
  feat_table <- readRDS(paste0(in_folder,"081016_cgp_new_feat_combat.rds"))
} else if (usage=="ccle"){
  feat_table <- readRDS(paste0(in_folder,"081616_cgp_new_feat_combat_ccle_based.rds"))
}

out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/" #For lab
out_table   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/" #For tigress
#####################################################################################################
# EXECUTE

test_table  <- feat_table[Compound==target_drug]
temp_table  <- feat_table[Compound!=target_drug]

if (modifier=="target_cells"){
  target_cells <- unique(test_table$cell_name)
  temp_table   <- temp_table[cell_name %in% target_cells,]

} else if (modifier=="target_drugs"){
  drug_met_cor <- cor( acast(MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
  drug_met_cor <- data.table(melt(drug_met_cor))
  drug_met_cor <- drug_met_cor[Var1==target_drug,][value>0,]

  target_drugs <- unique(drug_met_cor$Var2)
  temp_table   <- temp_table[Compound %in% target_drugs,]

}else {
  print("no modifier")
}

train_rows  <- sample(1:nrow(temp_table), 0.8*nrow(temp_table))
valid_rows  <- setdiff(1:nrow(temp_table), train_rows)

train_table <- temp_table[train_rows, ]
valid_table <- temp_table[valid_rows, ]

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, usage , ".", modifier, "_TRAIN_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, usage , ".", modifier, "_VALID_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table,  paste0(out_table, usage , ".", modifier, "_TEST_CGP_SEL.",  target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing sel tables")
