# cgp_baseline_mlp_plots.R
# Plots baseline results

library(data.table)
library(ggplot2)

Function.NRMSE <- function(pred, actual){

  NRMSE <- sqrt(mean((pred-actual)^2)) / diff(range(actual))
  return(NRMSE)
}
date_out    <- Sys.Date()

##############################################################################################################
# EXECUTE
folder_in = "/tigress/zamalloa/CGP_FILES/CGP_BASELINE_RESULTS/"

baseline_table <- data.table()

for (c in c(10, 20, 50, 100, 200)){
  for (d in c(10, 20, 50, 100, 200)){
    for (drug in c("Dasatinib", "Camptothecin", "Mitomycin_C", "Imatinib", "Cisplatin", "Sunitinib",
                     "Midostaurin",  "17-AAG", "Gefitinib", "Nilotinib", "Doxorubicin", "Vorinostat",
                     "Gemcitabine", "Cytarabine", "Axitinib", "Vinorelbine", "Methotrexate", "Shikonin",
                     "Etoposide", "Paclitaxel", "Embelin", "Cyclopamine", "PAC-1", "Bleomycin",
                     "Docetaxel", "Rapamycin", "ATRA", "Sorafenib", "Erlotinib", "Temsirolimus",
                     "Parthenolide", "Lapatinib", "Pazopanib", "Vinblastine", "Bortezomib", "Pyrimethamine",
                     "Elesclomol", "Roscovitine")){

                       print(c(drug, c, d))

                       drug_file  <- paste0(folder_in, "combined_D_values.Baseline_", drug,
                                           "_C_", c, "_D_", d, ".txt")
                       drug_file  <- data.table(read.table(drug_file, sep="\t", header=T))
                       max_epoch  <- max(drug_file$EPOCH)
                       drug_file  <- drug_file[EPOCH==max_epoch,]
                       drug_cor   <- cor(drug_file$ACTUAL, drug_file$PREDICTED, method="pearson")
                       drug_NRMSE <- Function.NRMSE(drug_file$PREDICTED, drug_file$ACTUAL)

                       drug_file  <- data.table(cell_name  = "place_holder",
                                                Compound   = drug,
                                                NORM_pIC50 = drug_file$ACTUAL,
                                                Prediction = drug_file$PREDICTED,
                                                Cor        = drug_cor,
                                                NRMSE      = drug_NRMSE,
                                                C          = c,
                                                D          = d)

                       baseline_table <- rbind(baseline_table, drug_file)
                     }
  }
}

##############################################################################################################
# PLOT
folder_out <- "/tigress/zamalloa/CGP_FILES/CGP_BASELINE_FIGURES/"

pdf(paste0(folder_out, date_out, ".baseline_mlp_C_D.pdf"), width=10, height=10)

print(ggplot(unique(baseline_table[,c("Compound", "Cor", "C", "D"),with=F]), aes(x="Correlation", y= Cor)) + geom_boxplot() + geom_jitter(size=0.5) +
          facet_grid(C~D) + theme_bw() + xlab("Drug features used for Model") + ylab("Cell features used for Model") +
          ggtitle("DL CGP Drug vs Cells features in MLP \n MLP applied to entire set \n
                  Number of Drug and Cell features pre-filtered prior to model based on variance"))

print(ggplot(unique(baseline_table[,c("Compound", "Cor", "C", "D"),with=F])[,list(MEDIAN=median(Cor)), by=c("C","D")],
             aes(C, MEDIAN, fill=factor(D) )) +
        geom_bar(stat="identity", position="dodge") + scale_fill_brewer(palette = "Set1") + theme_bw() +
        ggtitle("DL CGP Drug vs Cells features in MLP \n MLP applied to entire set \n
                Number of Drug and Cell features pre-filtered prior to model based on variance"))

dev.off()
##############################################################################################################
# WRITE TABLE
saveRDS(baseline_table, paste0(folder_out, date_out, ".baseline_mlp_C_D.rds"))

print("Done with baseline analysis")
