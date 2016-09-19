# cgp_baseline.R
# Function to calculate and plot baseline results

library(data.table)
library(ggplot2)
library(reshape2)
library(parallel)

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

Function_cgp_weight_separate <- function(cgp_table, target_drug, cell_w, drug_w, exponential=F) {

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

Function_cgp_weight_compounded <- function(cgp_table, target_drug, exponential=F) {
  # Uses whole vector to calculate distances

  # It is numerically faster to calculate the distance per each sample being analyzed
  train_table  <- cgp_table[Compound!=target_drug,]
  target_table <- cgp_table[Compound==target_drug,]

  predictions  <- apply(target_table, 1, function(x) {

    target_vector <- as.numeric(x[4:ncol(target_table)])

    weights       <- apply(train_table, 1, function(y) {

      train_vector  <- as.numeric(y[4:ncol(train_table)])
      return(cor(target_vector, train_vector, method="pearson"))
    })

    if (exponential==T){
      weights <- Function.range.0.1(exp(weights))
    }

    prediction <- sum(weights * train_table$NORM_pIC50) / sum(weights)
    return(prediction)
  })

  #Return formatted predictions
  output       <- data.table(cell_name  = target_table$cell_name,
                             Compound   = target_drug,
                             NORM_pIC50 = target_table$NORM_pIC50,
                             Prediction = prediction)

  output$Cor   <- cor(output$NORM_pIC50, output$Prediction)
  output$NRMSE <- Function.NRMSE(output$Prediction, output$NORM_pIC50)

  return(output)
}

Function_cgp_model_baseline <- function(cgp_table, target_drug, cell_w, drug_w, exponential=F, type="separate"){
  # Finds baseline prediction for target drug based on all the data and type of prediction
  # This will be in terms of weighted average
  # type can be:
    # "separate"   - using each cell_w and drug_w
    # "compounded" - using whole vector as distance in cgp_table rows [4:ncol]

  if (type == "separate"){

    output <- Function_cgp_weight_separate(cgp_table, target_drug, cell_w, drug_w, exponential=exponential)

  } else if (type == "compounded"){

    output <- Function_cgp_weight_compounded(cgp_table, target_drug, exponential=exponential)

  }

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

args        <- commandArgs(trailingOnly = TRUE)
date_out    <- Sys.Date()
file_name   <- args[1]

############################################################################################################################################
# EXECUTE
drugs         <- c("Gemcitabine", "Sunitinib",   "Doxorubicin", "Mitomycin C", "Vinorelbine",
                   "Vinblastine",  "Midostaurin", "Paclitaxel" , "Camptothecin", "Embelin",
                   "Bleomycin", "Axitinib", "Docetaxel", "Nilotinib", "Sorafenib", "Cytarabine",
                   "Shikonin", "Roscovitine", "Etoposide", "Pyrimethamine", "Methotrexate",
                   "PAC-1", "Temsirolimus", "Rapamycin", "Bortezomib", "Imatinib", "Pazopanib",
                   "Dasatinib",
                   "Lapatinib"  ,  "Erlotinib" ,  "Vorinostat" , "Cyclopamine", "Cisplatin",
                   "Elesclomol" ,  "17-AAG"    ,  "ATRA",        "Gefitinib"  , "Parthenolide")

# baseline <- data.table()
# for (d in drugs) {
#
#   print(d)
#
#   baseline <- rbind(baseline, Function_cgp_model_baseline(cgp_new_feat, d,
#                                                           Function_cell_weights(cgp_exp),
#                                                           Function_drug_weights(DRUGS.MET.PROFILE),
#                                                           exponential = T, type="compounded"))
#
# }

nodes<-detectCores()
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("as.data.table","data.table", "drugs", "cgp_new_feat", "cgp_exp", "DRUGS.MET.PROFILE",
                            "Function_cgp_model_baseline", "Function_cell_weights",
                            "Function_cgp_weight_separate", "Function_cgp_weight_compounded",
                            "Function_drug_weights", "Function.range.0.1", "Function.NRMSE",
                            "melt", "acast", "setnames", "setkey"),
                            envir=environment())
print ("Done exporting values")

baseline   <- parLapply(cl, drugs, function(d) {

  #Write to log to track
  write.table(d, "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/log", append=T)
  write.table("\t", "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/log", append=T)

  return(Function_cgp_model_baseline(cgp_new_feat, d,
                        Function_cell_weights(cgp_exp),
                        Function_drug_weights(DRUGS.MET.PROFILE),
                        exponential = T, type="compounded"))
  })

baseline <- do.call(rbind, baseline)

stopCluster(cl)

############################################################################################################################################
# WRITE and PLOT

out_folder <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_FIGURES/"

write.table(baseline, paste0(out_folder, date_out, file_name, ".cgp_baseline.txt"), quote=F, sep="\t", row.names=F, col.names=T)

#baseline <- fread(paste0(out_folder, date_out, file_name, ".cgp_baseline.txt"), header=T, sep="\t")
pdf(paste0(out_folder, date_out, ".cgp_baseline.pdf"), width=10, height=8)

ggplot(baseline, aes(NORM_pIC50, Prediction)) + geom_point(size=0.8) +
  theme_classic() +
  stat_smooth(method="lm", se = F, colour="red") + facet_wrap(~Compound)

ggplot(unique(baseline[,c("Compound", "Cor"),with=F]), aes(Compound, Cor)) + geom_bar(stat="identity") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) +
  ggtitle("CGP Baseline")

dev.off()
