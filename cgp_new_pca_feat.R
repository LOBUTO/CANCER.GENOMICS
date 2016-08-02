# cgp_new_pca_feat.R
# Function to split pca cgp table into training, validation and testing tables

library(data.table)
library(reshape2)

######################################################################################################
# LOAD DATA
feat_pca    <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080116_cgp_feat_table.rds")
target_drug <- "Erlotinib"
out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

#####################################################################################################
# EXECUTE
test_table  <- feat_pca[Compound==target_drug,]
temp_table  <- feat_pca[Compound!=target_drug,]

train_rows  <- sample(1:nrow(temp_table), 0.8*nrow(temp_table))
valid_rows  <- setdiff(1:nrow(temp_table), train_rows)

train_table <- temp_table[train_rows, ]
valid_table <- temp_table[valid_rows, ]

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, "TRAIN_PCA.", target_drug), quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, "VALID_PCA.", target_drug), quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table,  paste0(out_table, "TEST_PCA.", target_drug), quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing pca tables")
