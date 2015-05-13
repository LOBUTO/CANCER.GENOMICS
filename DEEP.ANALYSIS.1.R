library(data.table)
library(reshape2)
library(ggplot2)

#MANUAL BINARY
ALL<-do.call(rbind, lapply(paste("RESULTS/",list.files("RESULTS", pattern = "051215"), sep =""), function(x)  readRDS(x)) )

ggplot(ALL[METHOD=="Tanh",], aes(factor(FEATURES), TRAIN.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_wrap(~HIDDEN)
ggplot(ALL[METHOD=="Tanh",], aes(factor(FEATURES), TEST.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_wrap(~HIDDEN)
ggplot(ALL[METHOD=="TanhWithDropout",], aes(factor(FEATURES), TRAIN.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_grid(HIDDEN~INPUT.DR)
ggplot(ALL[METHOD=="TanhWithDropout",], aes(factor(FEATURES), TEST.ACC, colour=factor(FOLD) )) +geom_boxplot() + facet_grid(INPUT.DR~HIDDEN)

ALL.MEAN<-ALL[,list(TRAIN.ACC=mean(TRAIN.ACC), TEST.ACC=mean(TEST.ACC)), by=c("FOLD", "FEATURES", "METHOD", "HIDDEN", "INPUT.DR", "HIDDEN.DR")]
ggplot(ALL.MEAN, aes(factor(FEATURES), TRAIN.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)
ggplot(ALL.MEAN, aes(factor(FEATURES), TEST.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)

#NEW BINARY
ALL.BIN<-do.call(rbind, lapply(paste("RESULTS/",list.files("RESULTS", pattern = "051315"), sep =""), function(x)  readRDS(x)) )
ggplot(ALL.BIN, aes(factor(FEATURES), TRAIN.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)
ggplot(ALL.BIN, aes(factor(FEATURES), TEST.ACC, colour=HIDDEN)) + geom_boxplot() + facet_grid(INPUT.DR~METHOD)
