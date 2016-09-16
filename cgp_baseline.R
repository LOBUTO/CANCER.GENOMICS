# cgp_baseline.R
# Function to calculate and plot baseline results

library(data.table)
library(ggplot2)
library(reshape2)

Function_cell_weights <- function(cell_exp) {
  # Function to obtain cell-cell weights
  # Expression with samples in rows and genes in columns

  w <- data.table(melt(cor(cell_exp, method = "pearson")))
  setnames(w, c("cell_1", "cell_2", "w"))

  w$cell_1 <- as.character(w$cell_1)
  w$cell_2 <- as.character(w$cell_2)

  setkey(w)
  w <- unique(w)

  return(w)

}

Function_drug_weights <- function(drug_feat) {
  # Function to obtain drug-drug weights
  # Tables such as DRUGS.MET.PROFILE

  w <- acast(drug_feat, METABOLITE~DRUG, value.var = "TC")
  w <- data.table(melt(cor(w, method="pearson")))
  setnames(w, c("drug_1", "drug_2", "w"))

  w$drug_1 <- as.character(w$drug_1)
  w$drug_2 <- as.character(w$drug_2)

  setkey(w)
  w <- unique(w)

  return(w)
}

Function_cgp_model_baseline <- function(cgp_table, target_drug, cell_w, drug_w, exponential=F){
  # Finds baseline prediction for target drug based on all the data and type of prediction
  # This will be in terms of weighted average

  # Make sure all data is present from weights
  cgp_table    <- cgp_table[cell_name %in% cell_w[[1]],]
  cgp_table    <- cgp_table[Compound %in% drug_w[[1]],]
  cgp_table    <- cgp_table[, c("Compound", "cell_name", "NORM_pIC50"), with=F]

  # Split into training and testing set
  train_table  <- cgp_table[Compound!=target_drug,]
  test_table   <- cgp_table[Compound==target_drug,]

  # Calculate combined weights
  setnames(drug_w, c("drug_1", "drug_2", "weight_d"))
  setnames(cell_w, c("cell_1", "cell_2", "weight_c"))

  drug_w       <- drug_w[drug_2 == target_drug,]
  cell_w       <- cell_w[cell_2  %in% test_table$cell_name,]

  weight_table <- merge(train_table,  drug_w, by.x="Compound",  by.y="drug_1", allow.cartesian=TRUE)
  weight_table <- merge(weight_table, cell_w, by.x="cell_name", by.y="cell_1", allow.cartesian=TRUE)
    #[Compound, cell_name, NORM_pIC50, drug_2, weight_d, cell_2, weight_c]

  setkey(weight_table)
  weight_table <- unique(weight_table)

  # Predicting each target cell's IC50 in other drugs
  if (exponential == T) {

    weight_table$weight_d <- Function.range.0.1(exp(weight_table$weight_d))
    weight_table$weight_c <- Function.range.0.1(exp(weight_table$weight_c))
  }
  print(weight_table)

  weight_table <- weight_table[, list(other_cell_pred = sum(weight_c * NORM_pIC50) / sum(weight_c),
                                      weight_d = weight_d),
                               by = c("cell_2", "Compound")]
  setkey(weight_table)
  weight_table <- unique(weight_table)
  print(weight_table)

  # Predicting each target cell's IC50 in target drug
  weight_table <- weight_table[, list(Prediction = sum(other_cell_pred * weight_d) / sum(weight_d) ),
                               by=c("cell_2")]

  # Outputing
  output       <- merge(test_table, weight_table, by.x="cell_name", by.y="cell_2")
  output$Cor   <- cor(output$NORM_pIC50, output$Prediction)
  output$NRMSE <- Function.NRMSE(output$Prediction, output$NORM_pIC50)

  # Return
  return(output)

}

Function.range.0.1 <- function(x){
  scale.1 <-  (x-min(x))/(max(x)-min(x))
  return(scale.1)
}

Function.NRMSE <- function(pred, actual){

  NRMSE <- sqrt(mean((pred-actual)^2)) / diff(range(actual))
  return(NRMSE)
}

############################################################################################################################################
# LOAD FILES

cgp_new_feat      <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080816_cgp_new_feat_combat.rds")
cgp_exp           <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716_cgp_new_exp.rds")
DRUGS.MET.PROFILE <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/082316.DRUG.MET.PROFILE.rds")

date_out    <- Sys.Date()
############################################################################################################################################
# EXECUTE
# "Vinblastine",  "Midostaurin", "Paclitaxel" , "Camptothecin",
# "Lapatinib"  ,  "Erlotinib" ,  "Vorinostat" , "Cyclopamine", "Cisplatin",

baseline <- data.table()
for (d in c("Gemcitabine", "Sunitinib",   "Doxorubicin", "Mitomycin C", "Vinorelbine",
            "Elesclomol" ,  "17-AAG"    ,  "ATRA",        "Gefitinib"  , "Parthenolide")){

  print(d)
  baseline <- rbind(baseline,
                         Function_cgp_model_baseline(cgp_new_feat, d,
                                               Function_cell_weights(cgp_exp),
                                               Function_drug_weights(DRUGS.MET.PROFILE),
                                               exponential = T)
  )
}

############################################################################################################################################
# WRITE and PLOT

out_folder <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/"

write.table(baseline, paste0(out_folder, date_out, ".cgp_baseline.txt"), quote=F, sep="\t", row.names=F, col.names=T)

pdf(paste0(out_folder, date_out, ".cgp_baseline.pdf"), width=10, height=8)

ggplot(baseline, aes(NORM_pIC50, Prediction)) + geom_point() +
  theme_classic() + ggtitle(paste0(unique(Compound), " - Cor: ", round(unique(Cor),2))) +
  stat_smooth(method="lm", se = F, colour="red")

ggplot(baseline[,c("Compound", "Cor"),with=F], aes(Compound, Cor)) + geom_bar(stat="identity") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  ggtitle("CGP Baseline")

dev.off()
