# cgp_new_pca_feat.R
# Function to split pca cgp table into training, validation and testing tables

library(data.table)
library(reshape2)

######################################################################################################
# LOAD DATA
feat_pca    <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080116_cgp_feat_table.rds")

args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
m_pca       <- as.numeric(args[2])
g_pca       <- as.numeric(args[3])

out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

#####################################################################################################
# EXECUTE
m_features  <- sapply(1:m_pca, function(x)  paste0("PC", x, ".M"))
g_features  <- sapply(1:g_pca, function(x)  paste0("PC", x, ".E"))
features    <- c(colnames(feat_pca)[1:3], g_features, m_features)
feat_pca    <- feat_pca[ , features, with=F]

test_table  <- feat_pca[Compound==target_drug,]
temp_table  <- feat_pca[Compound!=target_drug,]

train_rows  <- sample(1:nrow(temp_table), 0.8*nrow(temp_table))
valid_rows  <- setdiff(1:nrow(temp_table), train_rows)

train_table <- temp_table[train_rows, ]
valid_table <- temp_table[valid_rows, ]

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, "TRAIN_CGP_PCA.", target_drug, "_", m_pca, "_", g_pca),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, "VALID_CGP_PCA.", target_drug, "_", m_pca, "_", g_pca),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table,  paste0(out_table, "TEST_CGP_PCA.",  target_drug, "_", m_pca, "_", g_pca),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing pca tables")
