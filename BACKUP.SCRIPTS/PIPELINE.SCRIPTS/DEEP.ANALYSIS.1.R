library(data.table)
library(reshape2)
library(ggplot2)

#MANUAL BINARY
ALL<-do.call(rbind, lapply(paste("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/",list.files("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/", pattern = "051215"), sep =""), function(x)  readRDS(x)) )

ggplot(ALL[METHOD=="Tanh",], aes(factor(FEATURES), TRAIN.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_wrap(~HIDDEN) + theme(strip.text.x = element_text(size = 14))
ggplot(ALL[METHOD=="Tanh",], aes(factor(FEATURES), TEST.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_wrap(~HIDDEN)+ theme(strip.text.x = element_text(size = 14))
ggplot(ALL[METHOD=="TanhWithDropout",], aes(factor(FEATURES), TRAIN.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_grid(HIDDEN~INPUT.DR) + theme(strip.text.x = element_text(size = 14))
ggplot(ALL[METHOD=="TanhWithDropout",], aes(factor(FEATURES), TEST.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_grid(INPUT.DR~HIDDEN) + theme(strip.text.x = element_text(size = 14))

ALL.MEAN<-ALL[,list(TRAIN.ACC=mean(TRAIN.ACC), TEST.ACC=mean(TEST.ACC)), by=c("FOLD", "FEATURES", "METHOD", "HIDDEN", "INPUT.DR", "HIDDEN.DR")]
ggplot(ALL.MEAN, aes(factor(FEATURES), TRAIN.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)+ theme(strip.text.x = element_text(size = 14))
ggplot(ALL.MEAN, aes(factor(FEATURES), TEST.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)+ theme(strip.text.x = element_text(size = 14))

#NEW BINARY
ALL.BIN<-do.call(rbind, lapply(paste("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/",list.files("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/", pattern = "051315"), sep =""), function(x)  readRDS(x)) )
ggplot(ALL.BIN, aes(factor(FEATURES), TRAIN.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD) + theme.format + theme(strip.text.x = element_text(size = 14))
ggplot(ALL.BIN, aes(factor(FEATURES), TEST.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD) + theme.format+ theme(strip.text.x = element_text(size = 14))
