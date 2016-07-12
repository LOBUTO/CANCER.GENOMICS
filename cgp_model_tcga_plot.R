#cgp_model_tcga_plot.R

library(data.table)
library(ggplot2)
library(reshape2)
library(survival)
library(GGally)
library(grid)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  #plots <- c(list(...), plotlist)
  plots <- plotlist

  numPlots = length(plots)
  print(numPlots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]],
            vp = viewport(layout.pos.row = matchidx$row,
                          layout.pos.col = matchidx$col))
    }
  }
}

Function.classify.lived.pred <- function(x, sd.multiplier=1, effective="POS"){

  sd.factor <- sd(x) * sd.multiplier
  sd.mean <- mean(x)
  # sd.mean <- 0
  # sd.factor <- sd.multiplier

  above.sd <- x[x > (sd.mean + sd.factor)]
  below.sd <- x[x < (sd.mean - sd.factor)]

  if (effective=="NEG"){

    CASE  <- ifelse(x %in% below.sd , "EFFECTIVE",
                    ifelse(x %in% above.sd, "NOT_EFFECTIVE", "NO_CLASS"))

  } else {

    CASE  <- ifelse(x %in% below.sd , "NOT_EFFECTIVE",
                    ifelse(x %in% above.sd, "EFFECTIVE", "NO_CLASS"))

  }

  #Return
  return (CASE)
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
  prediction <- prediction[ACTUAL>50, ]

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

    pred.temp <- prediction[CANCER==x,]
    pred.temp$CASE <- Function.classify.lived.pred(pred.temp$PREDICTED, sd.multiplier=0.3, effective="POS")

    pred.temp <- pred.temp[CASE!="NO_CLASS",]

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

  #SURVIVAL PLOTS
  surv.plots <- lapply(cancers, function(cancer) {

    cancer.clinical <- prediction[CANCER==cancer,]

    test.survival<-survfit(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)
    SURV.DIFF <- survdiff(Surv(LIVED, STATUS)~CASE, data=cancer.clinical)

    P.VAL <- pchisq(SURV.DIFF$chisq, length(SURV.DIFF$n)-1, lower.tail = FALSE)

    # file.name <- paste0(FIGURES, target.name, pca,  "." , cancer, ".survival.pdf")
    # pdf(file.name, width=12, height=18)

    temp.plot <- ggsurv(test.survival, surv.col=c("black", "darkviolet")) + theme(legend.position="bottom") +
                  theme_classic() + ggtitle(paste0(cancer, " - P-val: ", round(P.VAL,3))) #+

    return(temp.plot)
    # dev.off()

  })

  file.name <- paste0(FIGURES, target.name, pca,  "." , ".survival.pdf")
  pdf(file.name, width=12, height=18)

  multiplot(plotlist = surv.plots, cols=3)

  dev.off()

}

#DONE
print ("Done plotting!!")
