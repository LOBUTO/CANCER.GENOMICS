# cgp_new_sel_nci_feat.R
# Function to predict nci60 drug activity using cgp_new on selected features

library(data.table)
library(reshape2)

######################################################################################################
# LOAD DATA

drug_hmdb   <- fread("TABLES/TCGA.TRAINING/NCI60.TC.HMDB.FP4",
                    header=T, colClasses=c("character", "character", "numeric"))

feat_table  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716_cgp_new_feat.rds")
MET.PROFILE <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.DRUG.MET.PROFILE.rds")
nci60.exp   <- readRDS("/home/zamalloa/Documents/FOLDER/OBJECTS/061916.NCI60.EXP.rds")
cgp_exp     <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716_cgp_new_exp.rds")
nci60.gi50  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci.gi50.rds")
nci_to_cgp  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci60_names_to_cgp.rds")

args        <- commandArgs(trailingOnly = TRUE)
target_drug <- args[1]
out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

#####################################################################################################
# EXECUTE

# Prep drug table first
target_nsc   <- as.character(nci_to_cgp[Compound==target_drug,]$NSC)
drug_hmdb    <- drug_hmdb[DRUG==target_nsc,]
drug_hmdb$DRUG <- target_drug
common_met   <- intersect(unique(drug_hmdb$METABOLITE), unique(MET.PROFILE$METABOLITE))
drug_table   <- rbind(drug_hmdb[METABOLITE %in% common_met,],
                    MET.PROFILE[METABOLITE %in% common_met,][!DRUG %in% target_drug,])

drug_table   <- cor( acast(drug_table, METABOLITE~DRUG, value.var = "TC")  , method="spearman")
drug_table   <- drug_table[target_drug, unique(MET.PROFILE$DRUG)]
drug_table   <- data.table(drug_table, keep.rownames=T)

setnames(drug_table, c("Compound", colnames(drug_table)[2:ncol(drug_table)]))

# Then expression table
nci60.exp    <- t(scale(t(nci60.exp)))
cgp_exp      <- t(scale(t(cgp_exp)))

common_genes <- intersect(rownames(nci60.exp), rownames(cgp_exp))
common_cells <- intersect(colnames(nci60.exp), colnames(cgp_exp))

nci60_cells  <- setdiff(colnames(nci60.exp), common_cells)

cell_table   <- cbind(nci60.exp[common_genes, nci60_cells],
                      cgp_exp[common_genes])
cell_table   <- cor(cell_table, method="spearman")
cell_table   <- cell_table[colnames(nci60.exp), colnames(cgp_exp)]
cell_table   <- data.table(cell_table, keep.rownames = T)
setnames(cell_table, c("cell_name", colnames(cell_table)[2:ncol(cell_table)]))

# Then merge it all
main_table   <- merge(nci.gi50, nci_to_cgp, by="NSC")
main_table   <- main_table[ ,c("cell_name", "Compound", "SCALE.ACT"),with=F]
setnames(main_table, c("cell_name", "Compound", "NORM_pIC50")) # To be in accordance with the pred.py script

main_table   <- merge(main_table, drug_table, by = "Compound")
main_table   <- merge(main_table, cell_table, by = "cell_name")

# Filter features and return
main_table   <- main_table[ , c("cell_name", "Compound", "NORM_pIC50", colnames(feat_table)[4:ncol(feat_table)] ), with= F]

write.table(main_table, paste0(out_table, "NCI_CGP_SEL_FEAT.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing nci60 features")
