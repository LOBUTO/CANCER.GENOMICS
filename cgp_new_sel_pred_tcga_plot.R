# cgp_new_sel_pred_tcga_plot.R
# Function to take predictions and graph them

library(data.table)
library(reshape2)
library(ggplot2)

Function.NRMSE <- function(pred, actual){

  NRMSE <- sqrt(mean((pred-actual)^2)) / diff(range(actual))
  return(NRMSE)
}

fun_length <- function(x){

  # Found @ http://stackoverflow.com/questions/23330279/ggplot2-annotate-labelling-geom-boxplot-with-position-dodge
  return(data.frame(y=median(x),label= paste0("n=", length(x))))
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

# Clean up for tcga
tcga_table$Actual <- factor(tcga_table$Actual, levels = c("Clinical Progressive Disease", "Stable Disease",
                                                              "Partial Response", "Complete Response"))

# Predict more stringently
tcga_table$Binary_response <- ifelse(as.character(tcga_table$Actual) %in% c("Clinical Progressive Disease", "Stable Disease"),
                                    "Uneffective", "Effective")

# Plot
if (nchar(extra)>1){
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", cancer , "_" , target_drug, "_", extra ,".pdf"), width=12, height=8)
} else {
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_", cancer , "_" , target_drug, ".pdf"), width=12, height=8)
}

ggplot(self_table, aes(Actual, Predicted)) + geom_point(size=1.5) +
  theme_bw() + stat_smooth(method="lm", se=F, color = "purple") +
  scale_colour_brewer(palette = "Set1") +
  ggtitle(paste0("CGP based prediction on CGP Compound: ", target_drug,
                 "- Cor:", self_pred$Cor, "\n", "All CGP cells"))

ggplot(tcga_table, aes(Actual, Predicted, colour=Actual)) + geom_boxplot() + geom_jitter(size=0.4) +
  stat_summary(aes(x = factor(Actual)), fun.data = fun_length, geom = "text", vjust = +0.5, size=4) +
  theme_bw() + scale_colour_brewer(palette = "Set1") +
  xlab("Clinical response") + ylab("Predicted response") +
  ggtitle(paste0(toupper(cancer), " - CGP based prediction for clinical TCGA drug response for Compound: ", target_drug, "\n", "N = ", nrow(tcga_table)))

ggplot(tcga_table, aes(Binary_response, Predicted, colour=Binary_response)) + geom_boxplot() + geom_jitter(size=0.4) +
  stat_summary(aes(x = factor(Binary_response)), fun.data = fun_length, geom = "text", vjust = +0.5, size=4) +
  theme_bw() + scale_colour_brewer(palette = "Set1") +
  xlab("Clinical binary response") + ylab("Predicted response") +
  ggtitle(paste0(toupper(cancer), " - CGP based prediction for binary clinical TCGA drug response for Compound: ", target_drug, "\n", "N = ", nrow(tcga_table)))

dev.off()

print("Done plotting")
