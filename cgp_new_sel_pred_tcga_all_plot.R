# cgp_new_sel_pred_tcga_all_plot.R
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
extra       <- args[1]

tcga_resp   <- readRDS("/home/zamalloa/Documents/FOLDER/TCGA_FILES/090616_fireshose_all_response.rds")
tcga_exp    <- readRDS("/home/zamalloa/Documents/FOLDER/TCGA_FILES/090616_fireshose_all_exp.rds")
cgp_new     <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/082916_cgp_new.rds")

in_folder   <- "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TCGA_NEW_RESULTS/"
out_folder  <- "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TCGA_NEW_RESULTS/"

date_out    <- Sys.Date()
#####################################################################################
# EXECUTE

# Trim response
all_exp_samples  <- unlist(lapply(names(tcga_exp), function(z) colnames(tcga_exp[[z]][["tumor"]])))
response_sel     <- tcga_resp[sample %in% all_exp_samples,][,list(N=length(sample)), by=c("cancer", "Compound", "binary_response")]
response_sel     <- response_sel[N>5,]
response_sel[, N_length := length(N), by=c("cancer", "Compound")]
response_sel     <- response_sel[N_length==2,]
response_sel     <- response_sel[,1:2,with=F]
setkey(response_sel)
response_sel     <- unique(response_sel)

# Obtain data per cancer-drug pair
self_list        <- apply(response_sel, 1, function(x) {

  cancer       <- x[1]
  target_drug  <- x[2]

  self_table   <- fread(paste0(in_folder, "cgp_new_modeling_tcga_", target_drug))

  return(self_table)
  })

tcga_list        <- apply(response_sel, 1, function(x) {

  cancer       <- x[1]
  target_drug  <- x[2]

  tcga_table   <- fread(paste0(in_folder, "cgp_new_modeling_tcga_", cancer, "_", target_drug))

  return(tcga_table)
  })

self_list  <- do.call(rbind, self_list)
tcga_list  <- do.call(rbind, tcga_list)

# Classify for those that we actually have cgp data for
tcga_list$cgp  <- ifelse(tcga_list$Compound %in% unique(cgp_new$Compound), "CGP", "Non-CGP")

# Predict for self
self_pred      <- self_list[,list(NRMSE = Function.NRMSE(Predicted, Actual),
                                   Cor   = cor(Predicted, Actual, method="pearson")), by = c("Compound")]

# Predict for tcga binary response
tcga_list$Binary_response <- ifelse(as.character(tcga_list$Actual) %in% c("clinical progressive disease", "stable disease"),
                                    "Uneffective", "Effective")

tcga_list$Binary_response <- factor(tcga_list$Binary_response, levels = c("Uneffective", "Effective"))
tcga_list[, P_val:= wilcox.test(  Predicted[which(Binary_response == "Effective")],
                                  Predicted[which(Binary_response == "Uneffective")],
                                  paired=F, alternative="greater" )$p.val , by=c("Compound", "Cancer")]

tcga_pval      <- tcga_list[,c("Compound", "Cancer", "P_val", "cgp"), with=F]
tcga_pval$self <- sapply(tcga_pval$Compound, function(x) unique(self_pred[Compound==x,]$Cor) )
setkey(tcga_pval)
tcga_pval <- unique(tcga_pval)

# Plot
if (nchar(extra)>1){
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_all_cancers_and_drugs_", extra ,".pdf"), width=10, height=8)
} else {
  pdf(paste0(out_folder, date_out, "cgp_new_modeling_all_cancers_and_drugs.pdf"), width=10, height=8)
}

ggplot(tcga_list, aes(Binary_response, Predicted, colour=Binary_response)) + geom_boxplot() + geom_jitter(size=0.4) +
  stat_summary(aes(x = factor(Binary_response)), fun.data = fun_length, geom = "text", vjust = +2, size=4) +
  theme_bw() + scale_colour_brewer(palette = "Set3") +
  facet_wrap(Cancer~Compound) +
  xlab("Clinical binary response") + ylab("Predicted response") +
  ggtitle("CGP based prediction for binary clinical TCGA drug response for all Compounds")

ggplot(tcga_pval, aes(x = Compound, y = -log(P_val), fill=Cancer)) + geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~Cancer) + scale_fill_brewer(palette = "Set3") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  geom_hline(aes(yintercept = -log(0.1)), color="black", linetype="dashed") +
  ylab("-log(P_val)") +
  ggtitle("CGP based prediction for binary clinical TCGA drug response for all Compounds \n -log(P-values)")

# Evaluate those for which we have performance in cgp and see if we do better as model increases
ggplot(tcga_pval[cgp=="CGP",], aes(self, -log(P_val), colour=Compound, label=Cancer)) +
  geom_point(size=0.8) + geom_label() +
  scale_fill_brewer(palette = "Set3") +
  geom_hline(aes(yintercept = -log(0.1)), color="black", linetype="dashed") +
  xlab("Cross-validation CGP accuracy in terms of correlation") +
  ylab("-log(P_val)") +
  ggtitle("TCGA drug-cancer prediction accuracy vs Cross-validation CGP Accuracy")

dev.off()

print("Done plotting")
