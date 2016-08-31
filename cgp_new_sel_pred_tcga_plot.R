# cgp_new_sel_pred_tcga_plot.R
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
cancer      <- args[2]
extra       <- args[3]

in_folder   <- "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TCGA_NEW_RESULTS/"
out_folder  <- "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TCGA_NEW_RESULTS/"

date_out    <- Sys.Date()
#####################################################################################
# EXECUTE

self_table     <- fread(paste0(in_folder, "cgp_new_modeling_tcga_", target_drug))
tcga_table     <- fread(paste0(in_folder, "cgp_new_modeling_tcga_", cancer, "_", target_drug))

# Predict for self
self_pred      <- self_table[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                   Cor   = cor(Predicted, Actual, method="pearson")), by = "Compound"]

# Plot
if (nchar(extra)>1){
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", cancer , "_" , target_drug, "_", extra ,".pdf"), width=12, height=8)
} else {
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", cancer , "_" , target_drug, ".pdf"), width=12, height=8)
}

ggplot(self_table, aes(Actual, Predicted)) + geom_point(size=1.5) +
  theme_bw() + stat_smooth(method="lm", se=F, color = "purple") +
  ggtitle(paste0("CGP based prediction on CGP Compound: ", target_drug,
                 "- Cor:", self_pred$Cor, "\n", "All CGP cells"))

ggplot(tcga_table, aes(Actual, Predicted, colour=Actual)) + geom_boxplot() + gemo_jitter(size=0.4) +
  theme_bw() + xlab("Response") + ylab("Predicted response") +
  ggtitle(paste0(toupper(cancer), " - CGP based prediction on TCGA samples for Compound: ", target_drug, "\n", "N = ", nrow(tcga_table)))

dev.off()

print("Done plotting")
