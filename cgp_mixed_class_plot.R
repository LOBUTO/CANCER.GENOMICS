# cgp_mixed_class_plot.R
# Plots classes/regression mixed results

library(data.table)
library(ggplot2)
library(RColorBrewer)

# FUNCTIONS
Function_extract_max_acc <- function(epoch_file){

  epoch_table <- fread(epoch_file, sep="\t", header=T)

  if(nrow(epoch_table) > 1 ){

    max_epoch   <- max(epoch_table$EPOCH)
    epoch_table <- epoch_table[EPOCH == max_epoch,]
    mean_acc    <- mean(epoch_table$ACTUAL == epoch_table$PREDICTED)

  } else {
    mean_acc    <- "None"
  }

  return(mean_acc)
}

Function_extract_max_cor <- function(epoch_file){

  epoch_table <- fread(epoch_file, sep="\t", header=T)

  if(nrow(epoch_table) > 1 ){

    max_epoch   <- max(epoch_table$EPOCH)
    epoch_table <- epoch_table[EPOCH == max_epoch,]
    cor_acc    <- cor(epoch_table$ACTUAL, epoch_table$PREDICTED, method="pearson")

  } else {
    cor_acc    <- "None"
  }

  return(cor_acc)
}

#################################################### LOAD INPUT ##################################################
args      <- commandArgs(trailingOnly = TRUE)
mlp_type  <- args[1]
neural_type <- args[2]
file_out    <- args[3]
date      <- Sys.Date()

if (mlp_type == "classification") {
  Function_extract_max_epoch <- Function_extract_max_acc
  FOLDER                     <- "/tigress/zamalloa/CGP_FILES/CLASS_RESULTS"
  #FOLDER_IN                  <- "/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/combined_D_values.all_scaled_C_"
  FOLDER_IN                  <- "/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/combined_D_values.all_split_scaled_C_"
  file_out                   <- paste0("/tigress/zamalloa/FIGURES/CGP_CLASS/", date, "_", file_out, "_", neural_type ,".pdf")

  y_lim                      <- c(0.5, 0.9)

} else if (mlp_type == "regression"){
  Function_extract_max_epoch <- Function_extract_max_cor
  FOLDER                     <- "/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS"
  FOLDER_IN                  <- "/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS/combined_D_values.all_scaled_C_"
  file_out                   <- paste0("/tigress/zamalloa/FIGURES/CGP_REGRESSION/", date, "_" ,file_out, "_", neural_type ,".pdf")

  y_lim                      <- c(0, 0.65)

}

#################################################### LOAD DATA ###################################################
# Load drug_cor feature data
# drug_feat <- lapply(c(10, 20, 50, 100, 200, 500, 750), function(cell){
#
#   temp_table <- data.table()
#   for (drug in c(0, 10, 20, 50, 100, 200, 250)){
#
#     file_in <- paste0(FOLDER_IN, cell, "_D_", drug, ".txt")
#
#     if (file_in %in% list.files(FOLDER, full.names=T)){
#       acc     <- Function_extract_max_epoch(file_in)
#
#       temp_table <- rbind(temp_table, data.table(CELL = cell,
#                                                  DRUG = drug,
#                                                  ACC  = acc))
#     } else {
#       temp_table <- rbind(temp_table, data.table())
#     }
#   }
#   return(temp_table)
# })
#
# drug_feat <- do.call(rbind, drug_feat)

# Load drug morgan bit feature data
bit_feat  <- lapply(c(10, 50, 250), function(cell){

  temp_table <- data.table()
  for (bit in c(10, 50, 250)){

    file_in <- paste0(FOLDER_IN, cell, "_MB_", bit, "_mf_T_dn_250_cn_250_fn_250", ".txt")

    if (file_in %in% list.files(FOLDER, full.names=T)){
      acc     <- Function_extract_max_epoch(file_in)

      temp_table <- rbind(temp_table, data.table(CELL = cell,
                                                 DRUG = bit,
                                                 ACC  = acc))
    } else {
      temp_table <- rbind(temp_table, data.table())
    }
  }
  return(temp_table)
})

bit_feat  <- do.call(rbind, bit_feat)

# Load drug morgan count feature data
# count_feat <- lapply(c(10, 20, 50, 100, 200, 500, 750), function(cell){
#
#   temp_table <- data.table()
#   for (count in c(0, 10, 20, 50, 100, 200, 500, 1000, 2000)){
#
#     file_in <- paste0(FOLDER_IN, cell, "_MC_", count, ".txt")
#
#     if (file_in %in% list.files(FOLDER, full.names=T)){
#
#       acc     <- Function_extract_max_epoch(file_in)
#
#       temp_table <- rbind(temp_table, data.table(CELL = cell,
#                                                  DRUG = count,
#                                                  ACC  = acc))
#     } else{
#       temp_table <- rbind(temp_table, data.table())
#     }
#   }
#   return(temp_table)
# })
#
# count_feat <- do.call(rbind, count_feat)
#################################################### PROCESS DATA ###################################################
# drug_feat$type  <- "drug_cor"
# count_feat$type <- "morgan_count"
bit_feat$type   <- "morgan_bit"

# all_feat        <- do.call(rbind, list(drug_feat, bit_feat, count_feat))
all_feat        <- do.call(rbind, list(bit_feat))
all_feat        <- all_feat[ACC!="None",]
all_feat$ACC    <- as.numeric(all_feat$ACC)
print(all_feat)

# Calculate percent gains
# all_feat[, Baseline:= ACC[which(DRUG==0)], by=c("type", "CELL")]
# all_feat <- all_feat[!is.na(Baseline),]
#
# all_feat$Percent_Gain <- (all_feat$ACC - all_feat$Baseline) / all_feat$Baseline
# all_feat$Percent_Gain_per_drug_feture <- all_feat$Percent_Gain / all_feat$DRUG
#
# all_feat$Percent_Gain[is.na(all_feat$Percent_Gain)] <- 0
# all_feat$Percent_Gain_per_drug_feture[is.na(all_feat$Percent_Gain_per_drug_feture)] <- 0

#################################################### PLOT DATA ###################################################
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
print(all_feat)
pdf(file_out, width=12, height=10)
ggplot(all_feat, aes(factor(CELL), ACC, fill = factor(DRUG))) + geom_bar(stat="identity", position="dodge") +
  facet_wrap(~type) + theme_classic() +
  scale_fill_manual(values = getPalette( length(unique(all_feat$DRUG)) )) +
  coord_cartesian(ylim = y_lim)

# ggplot(all_feat, aes(factor(CELL), Percent_Gain, fill = factor(DRUG))) + geom_bar(stat="identity", position="dodge") +
#   facet_wrap(~type) + theme_classic() +
#   scale_fill_manual(values = getPalette( length(unique(all_feat$DRUG)) ))
#
# ggplot(all_feat, aes(DRUG, Percent_Gain, colour = factor(CELL))) + geom_point(size=0.8) + geom_line() +
#   facet_wrap(~type, scales="free_x") + theme_bw() +
#   scale_colour_brewer(palette = "Set3")
#
# ggplot(all_feat, aes(CELL, Percent_Gain_per_drug_feture, colour = factor(DRUG))) + geom_point(size=0.8) + geom_line() +
#   facet_wrap(~type) + theme_bw() +
#   scale_fill_manual(values = getPalette( length(unique(all_feat$DRUG)) ))

dev.off()

print("Done plotting")
