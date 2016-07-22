#cgp_pic50_testing.R
library(data.table)
library(ggplot2)
library(reshape2)

#####################################################################################################################
#LOAD FILES
args <- commandArgs(trailingOnly = TRUE)
main_table <- fread(args[1], sep="\t") #Like CGP_FILES/NCI_RESULTS/nci60_prediction_ErlotinibTH_0.85
nci_cgp_cor <- fread(args[2], sep="\t") #Like CGP_FILES/nci60_cgp_cor.csv
target_drug <- args[3]
nci_thr <- args[4]

file_out <- paste0("FIGURES/CGP.NCI/", target_drug, "_", nci_thr, ".pdf")

#####################################################################################################################
#EXECUTE
main_table <- main_table[,list(PRED_COR=cor(ACTUAL, PREDICT)), by=c("MODEL_DRUG", "NSC")]
main_table <- merge(main_table, nci_cgp_cor, by ="NSC")

#Plot accuracy versus similarity
main_table <- main_table[order(-PRED_COR),]


pdf(file_out, width=12, height=8)

print(ggplot(main_table, aes(factor(COR), PRED_COR)) + geom_bar(position="dodge", stat="identity" ) +
          theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
      )

print(ggplot(main_table, aes(COR, PRED_COR)) + geom_point(size=0.5) +
          theme_classic()
      )

dev.off()
