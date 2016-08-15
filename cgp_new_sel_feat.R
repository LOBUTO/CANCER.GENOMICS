# cgp_new_sel_feat.R
# Function to split cgp table

library(data.table)
library(reshape2)

######################################################################################################
# LOAD DATA

feat_table <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_cgp_new_feat_combat.rds") #For lab
feat_table <- readRDS("/tigress/zamalloa/CGP_FILES/081016_cgp_new_feat_combat.rds") #For tigress

args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
target_drug <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ")

out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/" #For lab
out_table   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/" #For tigress
#####################################################################################################
# EXECUTE

test_table  <- feat_table[Compound==target_drug]
temp_table  <- feat_table[Compound!=target_drug]

train_rows  <- sample(1:nrow(temp_table), 0.8*nrow(temp_table))
valid_rows  <- setdiff(1:nrow(temp_table), train_rows)

train_table <- temp_table[train_rows, ]
valid_table <- temp_table[valid_rows, ]

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, "TRAIN_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, "VALID_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table,  paste0(out_table, "TEST_CGP_SEL.",  target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing sel tables")
