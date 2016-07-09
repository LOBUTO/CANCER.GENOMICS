###Script to plot predictions from nci60/hybrid models

library(data.table)
library(ggplot2)
library(reshape2)
library(survival)
library(GGally)

# args<-commandArgs(trailingOnly=T)
# if (length(args)>1){
#   target.name <- paste0(args[1:length(args)], collapse="_")
# } else {
#   target.name <- args[1]
# }
target.name <- "nci.60_based.pca700_model.700.700_drugs_"

IN_FOLDER <- "/tigress/zamalloa/RESULTS/TCGA.TRAINING/" #For tigress
IN_FOLDER <- "/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/" #For Lab

FIGURES <- "/tigress/zamalloa/FIGURES/TCGA.TRAINING/" #For tigress
FIGURES <- "/home/zamalloa/Documents/FOLDER/FIGURES/TCGA.TRAINING/" #For Lab

#Load prediction table
prediction <- fread(paste0(IN_FOLDER, "tcga_prediction_table.txt"), header=T)

#Load original clinical
master.clinical <- fread("/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.ALL.TCGA.CLINICAL.csv", header=T)

#Execute
drugs <- unique(prediction$DRUG)

prediction <- merge(prediction, master.clinical[,c("SAMPLE", "LIVED"),with=F] , by ="SAMPLE")

prediction[,MEAN:=mean(PREDICTED), by="DRUG"]
print(prediction)
prediction <- prediction[MEAN<=0.9,]
prediction$MEAN <- NULL

print(prediction)

P.VALS <- sapply(drugs, function(x) {
  print(x)
  p.value <- wilcox.test(prediction[DRUG==x,][PREDICTED==1,]$LIVED,
                         prediction[DRUG==x,][PREDICTED==0,]$LIVED,
                         paired=F, alternative="greater")$p.value
  return(p.value)
  } )
P.VALS <- data.table(DRUG = drugs, P.VAL = P.VALS)

#Box plot
file.name <- paste0(FIGURES, target.name, ".boxplot.pdf")
pdf(file.name, width=12, height=8)

print(prediction)
ggplot(prediction, aes(factor(DRUG), LIVED)) +
    geom_boxplot(aes(fill=factor(PREDICTED))) + geom_jitter(colour="steelblue4", size=0.2) +
    geom_text(data=P.VALS, aes(x=DRUG, y=4000, label=paste0("P-val=", round(P.VAL,3)) )) +
    scale_fill_brewer(palette="Set1")
dev.off()

#Survival plot
master.clinical$STATUS <- ifelse(master.clinical$DEATH=="[Not Applicable]", 0, 1)

master.clinical <- merge(master.clinical, prediction[,c("SAMPLE", "PREDICTED", "DRUG"),with=F], by="SAMPLE")
master.clinical$CASE <- ifelse(master.clinical$PREDICTED==1, "EFFECTIVE", "NOT_EFFECTIVE")

for (drug in drugs){

  cancer.clinical <- master.clinical[DRUG==drug,]

  test.survival<-survfit(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)
  SURV.DIFF <- survdiff(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)

  P.VAL <- pchisq(SURV.DIFF$chisq, length(SURV.DIFF$n)-1, lower.tail = FALSE)

  file.name <- paste0("/home/zamalloa/Documents/FOLDER/FIGURES/TCGA.TRAINING/",
                      target.name, "." , drug, ".survival.pdf")
  pdf(file.name, width=12, height=18)

  print(ggsurv(test.survival, surv.col=c("black", "darkviolet")) + theme(legend.position="bottom") +
          theme_classic() +
          geom_text(aes(mean(cancer.clinical$LIVED), 0.85, label= P.VAL), size=8.0))
  dev.off()

}


#DONE
print ("Done plotting!!")
