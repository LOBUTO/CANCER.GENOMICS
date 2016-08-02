# ccle_pca_predict.R
# Functions to plot predictions of ccle

library(ggplot2)
library(data.table)
library(reshape2)

##########################################################################################
# Load files

args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
m_pca       <- as.numeric(args[2])
g_pca       <- as.numeric(args[3])

file_in     <- paste0("/home/zamalloa/Documents/FOLDER/CGP_FILES/CCLE_PRED/cgp_pca_model_M", m_pca ,"_G_", g_pca, "_D_", target_drug)

##########################################################################################
# Execute

predictions <- fread(file_in, sep="\t", header=T)

print (mean(predictions$ACTUAL==predictions$PREDICTION))


print("Done prediction")
