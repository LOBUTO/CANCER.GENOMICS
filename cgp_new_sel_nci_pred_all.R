# cgp_new_sel_nci_pred_all.R
# Stat analysis of all nci cgp-based predictions

library(data.table)
library(reshape2)
library(ggplot2)

Function.NRMSE <- function(pred, actual){

  NRMSE <- sqrt(mean((pred-actual)^2)) / diff(range(actual))
  return(NRMSE)
}

#####################################################################################
# LOAD DATA

cgp_table   <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_cgp_new_feat_combat.rds")
nci60.gi50  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci.gi50.rds")
nci_to_cgp  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci60_names_to_cgp.rds")

out_folder  <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/"
in_folder   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_RESULTS/"

all_drugs   <- c("Sunitinib", "Mitomycin C", "Midostaurin", "Imatinib", "Cisplatin",
                 "Camptothecin", "17-AAG", "Dasatinib", "Gefitinib", "Nilotinib")
date_out    <- Sys.Date()

#####################################################################################
# EXECUTE

main_table <- lapply(all_drugs, function(target_drug) {

  drug_table     <- fread(paste0(in_folder, "cgp_new_modeling_nci_", target_drug))

  # Predict for all
  drug_all_pred  <- drug_table[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                     Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

  # Predict for common
  main_table     <- merge(nci60.gi50, nci_to_cgp, by="NSC")
  main_table     <- main_table[ ,c("cell_name", "Compound", "SCALE.ACT"),with=F]
  main_table     <- merge(main_table, cgp_table, by=c("Compound", "cell_name"))
  common_cells   <- main_table$cell_name

  drug_table_com <- drug_table[cell_name %in% common_cells,]
  drug_comb_pred <- drug_table_com[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                         Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

  # Predict for uncommon
  drug_table_unc <- drug_table[!cell_name %in% common_cells,]
  drug_unc_pred  <- drug_table_unc[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                        Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

  # Predict for cgp self
  drug_table_cgp <- fread(paste0(in_folder, "cgp_new_modeling_cgp_", target_drug))

  drug_cgp_pred  <- drug_table_cgp[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                    Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

  # Add to main table
  return(data.table(Compound = target_drug,
                    Prediction = c(drug_all_pred$Cor, drug_comb_pred$Cor,
                                   drug_unc_pred$Cor, drug_cgp_pred$Cor),
                    Count      = c(nrow(drug_table), nrow(drug_table_com),
                                   nrow(drug_table_unc), nrow(drug_table_cgp)),
                    Type       = c("All NCI-60", "Common NCI-60/CGP",
                                   "Non-CGP NCI-60", "Self CGP"),
                    Perc_common= c(nrow(drug_table_com) / nrow(drug_table)),
                    cgp_pred   = drug_cgp_pred$Cor
                    ))
})

main_table <- do.call(rbind, main_table)

# Plot
pdf(paste0(out_folder, date_out, "cgp_new_modeling_nci_all.pdf"), width=12, height=8)

ggplot(main_table, aes(Type, Prediction, colour=Type)) + geom_boxplot() + geom_jitter(size=0.4) +
  theme_bw() + scale_colour_brewer(palette="Set1") +
  ggtitle("CGP-based predictions comparisson on NCI-60") + xlab("Type of comparisson") +
    ylab("Accuracy in terms of correlation")

ggplot(main_table[Type=="All NCI-60",], aes(Perc_common * 100, Prediction)) + geom_point(size=1.5) +
  theme_bw() +
  ggtitle("CGP-based predictions on all NCI-60 cells") +
  xlab("Percent common cell count") + ylab("Accuracy in terms of correlation")

ggplot(main_table[Type=="All NCI-60",], aes(cgp_pred, Prediction)) + geom_point(size=1.5) +
  theme_bw() +
  ggtitle("Correlation between CGP-based model self-dataset accuracy and all NCI-60 cells") +
  xlab("CGP accuracy in terms of correlation") + ylab("NCI-60 accuracy in terms of correlation")

ggplot(main_table[Type=="Common NCI-60/CGP",], aes(Perc_common * 100, Prediction)) + geom_point(size=1.5) +
  theme_bw() +
  ggtitle("CGP-based predictions on common cells across CGP and NCI-60 datasets") +
  xlab("Percent common cell count") + ylab("Accuracy in terms of correlation")

ggplot(main_table[Type=="Common NCI-60/CGP",], aes(cgp_pred, Prediction)) + geom_point(size=1.5) +
  theme_bw() +
  ggtitle("Correlation between CGP-based model self-dataset accuracy and common NCI-60/CGP cells") +
  xlab("CGP accuracy in terms of correlation") + ylab("NCI-60 accuracy in terms of correlation")

dev.off()

print("Done plotting")
