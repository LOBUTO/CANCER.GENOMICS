#cgp_model_tcga_plot.R

library(data.table)
library(ggplot2)
library(reshape2)
library(survival)
library(GGally)

Function.classify.lived.pred <- function(x, sd.multiplier=1, effective="NEG"){

  sd.factor <- sd(x) * sd.multiplier
  sd.mean <- mean(x)
  #sd.mean <- 0

  above.sd <- x[x > (sd.mean + sd.factor)]
  below.sd <- x[x < (sd.mean - sd.factor)]

  if (effective=="NEG"){
    main.table <- data.table(PREDICTION=below.sd, CASE="EFFECTIVE")
    main.table <- rbind(main.table,
                        data.table(PREDICTION=above.sd, CASE="NOT_EFFECTIVE"))
  } else {
    main.table <- data.table(PREDICTION=above.sd, CASE="EFFECTIVE")
    main.table <- rbind(main.table,
                        data.table(PREDICTION=below.sd, CASE="NOT_EFFECTIVE"))
  }

  #Return
  return (main.table)
}

################################################################################################################################

target.name <- "cgp_based_model_for_tcga_pcas_"

IN_FOLDER <- "/tigress/zamalloa/RESULTS/TCGA.TRAINING/" #For tigress
IN_FOLDER <- "/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/" #For Lab

FIGURES <- "/tigress/zamalloa/FIGURES/TCGA.TRAINING/" #For tigress
FIGURES <- "/home/zamalloa/Documents/FOLDER/FIGURES/TCGA.TRAINING/" #For Lab

#Load original clinical
master.clinical <- fread("/home/zamalloa/Documents/FOLDER/TABLES/TCGA.TRAINING/062116.ALL.TCGA.CLINICAL.csv", header=T)
master.clinical$STATUS <- ifelse(master.clinical$DEATH=="[Not Applicable]", 0, 1)

#Load prediction table per PCA used
for (pca in c(500, 800, 1000)){

  prediction <- fread(paste0(IN_FOLDER, "cgp_auc_tcga_prediction_", pca), header=T)

  #Do we need to filter?
  #prediction <- prediction[ACTUAL>100, ]

  prediction[,COUNT:=length(SAMPLE), by="CANCER"]
  prediction <- prediction[COUNT>50,]
  prediction$COUNT <- NULL

  #Execute
  cancers <- unique(prediction$CANCER)
  prediction <- merge(prediction, master.clinical[,c("SAMPLE", "LIVED", "STATUS"),with=F] , by ="SAMPLE")

  cancer.cors <- sapply(cancers, function(x) {

    y <- cor(prediction[CANCER==x,]$ACTUAL,
            prediction[CANCER==x,]$PREDICTED,
            method="pearson")
    return(y)

    })
  cancer.cors <- data.table(CANCER = cancers, COR = cancer.cors)

  #REGRESSION PLOT
  file.name <- paste0(FIGURES, target.name, pca ,"_regression.pdf")
  pdf(file.name, width=12, height=8)

  print(prediction)
  print(ggplot(prediction, aes(ACTUAL, PREDICTED)) + geom_point(colour="steelblue4", size=0.2) +
    geom_text(data=cancer.cors, aes(x=1000, y=0.5, label=paste0("Cor=", round(COR,3)) )) +
    scale_fill_brewer(palette="Set1") + theme_bw() + geom_smooth(method="lm", se=T, color = "black", size=0.3) +
    facet_wrap(~CANCER, scales="free"))
  dev.off()

  #DEFINED CLASSES PLOT
  cancers <- unique(prediction$CANCER)

  pred.classes <- lapply(cancers, function(x) {

    pred.temp <- Function.classify.lived.pred(prediction[CANCER==x,]$PREDICTED, sd.multiplier=0.3, effective="POS")
    pred.temp <- merge(prediction[CANCER==x,]$PREDICTED, pred.classes, by.x="PREDICTED", by.y="PREDICTION")

    return(pred.temp)
    } )
  prediction <- do.call(rbind, pred.classes)

  P.VALS <- sapply(cancers, function(x) {

    p.value <- wilcox.test(prediction[CANCER==x,][CASE=="EFFECTIVE",]$LIVED,
                           prediction[CANCER==x,][CASE=="NOT_EFFECTIVE",]$LIVED,
                           paired=F, alternative="greater")$p.value
    return(p.value)
    } )
  P.VALS <- data.table(CANCER = cancers, P.VAL = P.VALS)

  file.name <- paste0(FIGURES, target.name, pca ,"_est.classes.pdf")
  pdf(file.name, width=12, height=8)

  print( ggplot(prediction, aes(CANCER, LIVED)) + geom_boxplot(aes(fill=factor(CASE))) +
    geom_jitter(colour="steelblue4", size=0.2) +
    geom_text(data=P.VALS, aes(x=CANCER, y=4000, label=paste0("P-val=", round(P.VAL,3)) )) +
    scale_fill_brewer(palette="Set1") + theme_bw() )

  dev.off()

  #Survival plot
  # for (cancer in cancers){
  #
  #   cancer.clinical <- master.clinical[CANCER==cancer,]
  #
  #   test.survival<-survfit(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)
  #   SURV.DIFF <- survdiff(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)
  #
  #   P.VAL <- pchisq(SURV.DIFF$chisq, length(SURV.DIFF$n)-1, lower.tail = FALSE)
  #
  #   file.name <- paste0(FIGURES, target.name, "." , cancer, ".survival.pdf")
  #   pdf(file.name, width=12, height=18)
  #
  #   print(ggsurv(test.survival, surv.col=c("black", "darkviolet")) + theme(legend.position="bottom") +
  #           theme_classic() +
  #           geom_text(aes(mean(cancer.clinical$LIVED), 0.85, label= P.VAL), size=8.0))
  #   dev.off()
  #
  # }
  #

}

#DONE
print ("Done plotting!!")
