#cgp.cor.AUC_class_pred_plot.R
#Plotting of predictions per drug compound for results obtained from cgp.cor.AUC_class_pred.py

library(data.table)
library(ggplot2)

IN_FOLDER <- "/tigress/zamalloa/RESULTS/TCGA.TRAINING/" #For tigress
IN_FOLDER <- "/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/" #For Lab

FIGURES <- "/tigress/zamalloa/FIGURES/TCGA.TRAINING/" #For tigress
FIGURES <- "/home/zamalloa/Documents/FOLDER/FIGURES/TCGA.TRAINING/" #For Lab

#############################################################################################
file_in <- paste0(IN_FOLDER, "cgp_auc_class_predction")

prediction <- fread(file_in, header=T)
prediction <- prediction[, list(Accuracy=mean(ACTUAL==PREDICTED), COUNT=length(ACTUAL)), by="DRUG"]
prediction$BIN <- cut(prediction$COUNT, 4)

#Plot across drugs
file.name <- paste0(FIGURES, "cgp.pca.nci60.based.predictions.pdf")
pdf(file.name, width=20, height=16)

ggplot(prediction, aes(DRUG, Accuracy, fill=BIN)) + geom_bar(stat="identity", position="dodge") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Set1") +
  ggtitle("CGP predictions using NCI60 Big Data MLP Model") + facet_grid(~BIN, scales="free_x")

dev.off()

#Boxplots of accuracies across bins
file.name <- paste0(FIGURES, "cgp.pca.nci60.based.predictions.boxplot.pdf")
pdf(file.name, width=16, height=12)

ggplot(prediction, aes(BIN, Accuracy, fill=BIN)) + geom_boxplot() +
  geom_jitter(colour="steelblue4", size=0.2) + scale_fill_brewer(palette="Set1") +
  theme_bw() +
  ggtitle("CGP predictions using NCI60 Big Data MLP Model") +
  geom_text(data=prediction[,list(MEAN=mean(Accuracy)), by="BIN"],
      aes(x=BIN, y=MEAN, label=MEAN), size=7)

dev.off()
#Done
print ("Done plotting!!")
