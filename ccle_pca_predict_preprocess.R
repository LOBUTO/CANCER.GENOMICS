# ccle_pca_predict_preprocess.R
# Function to apply pca transformation from cgp to ccle

library(data.table)
library(reshape2)

##########################################################################################################
# Load files

ccle_data <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080216_ccle_data.rds")
pca_data  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080216_pca_data.rds")

args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
m_pca       <- as.numeric(args[2])
g_pca       <- as.numeric(args[3])

out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"
##########################################################################################################
# Execute

# PCA on ccle data
genes_used     <- pca_data[["exp_pcas"]][["cgp_feat_order"]]
ccle_exp_pca   <- ccle_data[["exp_table"]]
ccle_exp_pca   <- ccle_exp_pca[genes_used, ]
ccle_exp_pca   <- scale(t(ccle_exp_pca)) %*% pca_data[["exp_pcas"]][["cgp"]]$rotation
ccle_exp_pca   <- data.table(ccle_exp_pca, keep.rownames = T)
setnames(ccle_exp_pca, c("cell_name", colnames(ccle_exp_pca)[2:ncol(ccle_exp_pca)]))

mets_used      <- pca_data[["drug_pcas"]][["cgp_feat_order"]]
ccle_met_pca   <- ccle_data[["met_table"]][METABOLITE %in% mets_used,]
ccle_met_pca   <- scale(acast(ccle_met_pca, DRUG~METABOLITE, value.var="TC")) %*% pca_data[["drug_pcas"]][["cgp"]]$rotation
ccle_met_pca   <- data.table(ccle_met_pca, keep.rownames = T)
setnames(ccle_met_pca, c("Compound", colnames(ccle_met_pca)[2:ncol(ccle_met_pca)]))

# Build PCA feature table for specific drug
feat_pca_table <- ccle_data[["main_table"]][,c("cell_name", "Compound", "pic50_class"), with=F]
feat_pca_table <- feat_table[pic50_class!=2,]
feat_pca_table <- merge(feat_pca_table, ccle_exp_pca, by = "cell_name")
feat_pca_table <- merge(feat_pca_table, ccle_met_pca, by = "Compound")

# Filter based on given features and target drug
m_features     <- sapply(1:m_pca, function(x)  paste0("PC", x, ".M"))
g_features     <- sapply(1:g_pca, function(x)  paste0("PC", x, ".E"))

features       <- c(colnames(feat_pca_table)[1:3], g_features, m_features)
feat_pca_table <- feat_pca_table[ , features, with=F]

feat_pca_table <- feat_pca_table[Compound == target_drug,]

##########################################################################################################
# Write out
write.table(feat_pca_table, paste0(out_table, "PRE_CCLE_PCA.", target_drug, "_", m_pca, "_", g_pca),
            quote=F, sep="\t", row.names=F, col.names=T)

print ("Done pre-processing")
