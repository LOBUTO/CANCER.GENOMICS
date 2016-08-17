# cgp_new_sel_nci_pred.R
# Function to take predictions and graph them

library(data.table)
library(reshape2)
library(ggplot2)

Function.NRMSE <- function(pred, actual){

  NRMSE <- sqrt(mean((pred-actual)^2)) / diff(range(actual))
  return(NRMSE)
}

#####################################################################################
# LOAD DATA
args        <- commandArgs(trailingOnly = TRUE)
target_drug <- args[1]
target_drug <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ")
extra       <- args[2]
usage       <- args[3]

if (usage=="nci60"){
  nci60.gi50  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci.gi50.rds")
  cgp_table   <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_cgp_new_feat_combat.rds")
} else if (usage=="ccle"){
  nci60.gi50  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081616_ccle.rds")
  cgp_table   <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081616_cgp_new_feat_combat_ccle_based.rds")
}

out_folder  <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/"
in_folder   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_RESULTS/"

date_out    <- Sys.Date()
#####################################################################################
# EXECUTE

#Obtain comparable cgp identifiers (Depending on usage)
drug_table     <- fread(paste0(in_folder, "cgp_new_modeling_nci_", target_drug)) #Keeping same name for now

if (usage=="nci60"){
  setnames(drug_table, c("Compound", "CELL", "Actual", "Predicted"))
  drug_table     <- merge(drug_table, unique(nci60.gi50[,c("CELL", "cell_name"),with=F]), by="CELL")
  drug_table$CELL<-NULL
} else if (usage=="ccle"){
  drug_table     <- drug_table
}

# Predict for all
drug_all_pred  <- drug_table[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                   Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

# Predict for common
common_cells   <- intersect(unique(drug_table$cell_name),
                            unique(cgp_table[Compound==target_drug,]$cell_name))

drug_table_com <- drug_table[cell_name %in% common_cells,]
drug_comb_pred <- drug_table_com[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                       Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

# Predict for uncommon
drug_table_unc <- drug_table[!cell_name %in% common_cells,]
drug_unc_pred  <- drug_table_unc[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                       Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

# Predict for cgp self
drug_table_cgp <- fread(paste0(in_folder, "cgp_new_modeling_cgp_", target_drug)) #Keeping same name for now

drug_cgp_pred  <- drug_table_cgp[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                   Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

# Plot
if (nchar(extra)>1){
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", usage , "_" , target_drug, "_", extra ,".pdf"), width=12, height=8)
} else {
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", usage , "_" , target_drug, ".pdf"), width=12, height=8)
}

ggplot(drug_table_cgp, aes(Actual, Predicted)) + geom_point(size=1.5) +
  theme_bw() + stat_smooth(method="lm", se=F, color = "purple") +
  ggtitle(paste0("CGP based prediction on CGP Compound: ", target_drug,
                 "- Cor:", drug_cgp_pred$Cor, "\n", "All CGP cells"))

ggplot(drug_table, aes(Actual, Predicted)) + geom_point(size=1.5) +
  theme_bw() + stat_smooth(method="lm", se=F, color = "purple" ) +
  ggtitle(paste0("CGP based prediction on ", usage, " Compound: ", target_drug ,
                 "- Cor:", drug_all_pred$Cor, "\n", "All ", usage, " Cells" ))

ggplot(drug_table_com, aes(Actual, Predicted)) + geom_point(size=1.5) +
 theme_bw() + stat_smooth(method="lm", se =F, color = "purple" ) +
 ggtitle(paste0("CGP based prediction on ", usage ," Compound: ", target_drug ,
                "- Cor:", drug_comb_pred$Cor, "\n", "Common CGP/", usage, " Cells" ))

ggplot(drug_table_unc, aes(Actual, Predicted)) + geom_point(size=1.5) +
 theme_bw() + stat_smooth(method="lm", se =F, color = "purple" ) +
 ggtitle(paste0("CGP based prediction on ", usage, " Compound: ", target_drug ,
                "- Cor:", drug_unc_pred$Cor, "\n", "Non-CGP ", usage, " Cells" ))

dev.off()

print("Done plotting")
