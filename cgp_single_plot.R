# cgp_single_plot.R

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
date      <- Sys.Date()

if (mlp_type == "classification") {
  Function_extract_max_epoch <- Function_extract_max_acc
  FOLDER                     <- "/tigress/zamalloa/CGP_FILES/CLASS_RESULTS"
  FOLDER_IN                  <- "/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/combined_D_values."
  file_out                   <- paste0("/tigress/zamalloa/FIGURES/CGP_CLASS/", date, "_mlp_class_single_cgp_", neural_type,".pdf")

  count_cells                <- c(500)
  bit_cells                  <- c(200)
  counts                     <- c(0, 500, 1000, 2000)
  bits                       <- c(10, 50 , 200)
  y_lim                      <- c(0.5, 0.9)

} else if (mlp_type == "regression"){
  Function_extract_max_epoch <- Function_extract_max_cor
  FOLDER                     <- "/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS"
  FOLDER_IN                  <- "/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS/combined_D_values."
  file_out                   <- paste0("/tigress/zamalloa/FIGURES/CGP_REGRESSION/", date, "_mlp_regression_single_cgp_", neural_type,".pdf")

  count_cells                <- c(200)
  counts                     <- c(0, 500, 1000)
  y_lim                      <- c(0, 0.65)

}

target_drugs <- gsub("\n", " ", "FK866 IPA-3 NSC-207895 UNC0638 CX-5461 Trametinib SNX-2112 OSI-027 QS11 AT-7519 PAC-1
SN-38 PI-103 I-BET-762 5-Fluorouracil PHA-793887 YM201636 LY317615 TAK-715 RDEA119 Gemcitabine
Bleomycin CHIR-99021 VX-11e EKB-569 GSK-650394 17-AAG AG-014699 GDC0941 Y-39983")
target_drugs <- strsplit(target_drugs, " ")[[1]]

#################################################### LOAD DATA ###################################################
# Load drug morgan count feature data
# count_feat <- lapply(count_cells, function(cell){
#
#   temp_table <- data.table()
#   for (count in counts){
#
#     for (drug in target_drugs){
#
#       file_in   <- paste0(FOLDER_IN, drug, "_scaled_C_" , cell, "_MC_", count, ".txt")
#
#       if (file_in %in% list.files(FOLDER, full.names=T)){
#
#         acc         <- Function_extract_max_epoch(file_in)
#
#         temp_table  <- rbind(temp_table, data.table(CELL = cell,
#                                                    DRUG = count,
#                                                    Compound = drug,
#                                                    ACC  = acc))
#       } else{
#         temp_table  <- rbind(temp_table, data.table())
#       }
#     }
#   }
#   return(temp_table)
# })
#
# count_feat <- do.call(rbind, count_feat)

# Load drug morgan bit feature data
bit_feat <- lapply(bit_cells, function(cell){

  temp_table <- data.table()
  for (bit in bits){

    for (drug in target_drugs){

      file_in   <- paste0(FOLDER_IN, drug, "_scaled_C_" , cell, "_MB_", bit, "_mf_T_dn_250_cn_250_fn_250", ".txt")

      if (file_in %in% list.files(FOLDER, full.names=T)){

        acc         <- Function_extract_max_epoch(file_in)

        temp_table  <- rbind(temp_table, data.table(CELL = cell,
                                                   DRUG = bit,
                                                   Compound = drug,
                                                   ACC  = acc))
      } else{
        temp_table  <- rbind(temp_table, data.table())
      }
    }
  }
  return(temp_table)
})

bit_feat <- do.call(rbind, bit_feat)

#################################################### PROCESS DATA ###################################################
# count_feat$type <- "morgan_count"
bit_feat$type <- "morgan_bit"

all_feat        <- do.call(rbind, list(bit_feat))
all_feat        <- all_feat[ACC!="None",]
all_feat$ACC    <- as.numeric(all_feat$ACC)

#################################################### PLOT DATA ###################################################
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

pdf(file_out, width=12, height=10)
ggplot(all_feat, aes(factor(CELL), ACC, fill = factor(DRUG))) + geom_boxplot() + geom_jitter(size=0.4) +
  facet_wrap(~type) + theme_classic() +
  scale_fill_manual(values = getPalette( length(unique(all_feat$DRUG)) )) +
  coord_cartesian(ylim = y_lim)

dev.off()
