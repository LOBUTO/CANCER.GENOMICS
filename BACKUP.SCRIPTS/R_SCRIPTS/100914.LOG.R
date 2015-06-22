#100914.LOG.R

###Produce accurate measure of background mutation rate
brca.table.1<-readRDS("BRCA/100514.BRCA.Table.1.rds")
GBM.table.1<-readRDS("GBM/100814.GBM.Table.1.rds")
COAD.table.1<-readRDS("COAD/100914.COAD.Table.1.rds")
OV.table.1<-readRDS("OV/101014.OV.Table.1.rds")

###Correlation of missense and silent mutations to total mutations in brca - RESULT: YES
brca.table.1.count<-brca.table.1$table.1[,list(Patient.Missense=sum(Missense), 
    Patient.Silent=sum(Silent), Patient.TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
brca.table.1.melt<-melt(brca.table.1.count, id.vars=c("PATIENT", "Patient.TOTAL"))
ggplot(brca.table.1.melt, aes(Patient.TOTAL, value, colour=variable)) + geom_point() + stat_smooth(method="lm") + theme(legend.position="bottom")
ggplot(brca.table.1.count, aes(Patient.TOTAL, Patient.Silent)) + geom_point()  + stat_smooth(method="lm")

###Does this happen in other cancers? Try a few - RESULT: YES
GBM.table.1.count<-GBM.table.1$table.1[,list(Patient.Missense=sum(Missense), 
    Patient.Silent=sum(Silent), Patient.TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
GBM.table.1.melt<-melt(GBM.table.1.count, id.vars=c("PATIENT", "Patient.TOTAL"))
ggplot(GBM.table.1.melt, aes(Patient.TOTAL, value, colour=variable)) + geom_point()  + stat_smooth(method="lm") + theme(legend.position="bottom")
ggplot(GBM.table.1.count, aes(Patient.TOTAL, Patient.Silent)) + geom_point()  + stat_smooth(method="lm")

COAD.table.1.count<-COAD.table.1$table.1[,list(Patient.Missense=sum(Missense), 
    Patient.Silent=sum(Silent), Patient.TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
COAD.table.1.melt<-melt(COAD.table.1.count, id.vars=c("PATIENT", "Patient.TOTAL"))
ggplot(COAD.table.1.melt, aes(Patient.TOTAL, value, colour=variable)) + geom_point()  + stat_smooth(method="    ") + theme(legend.position="bottom")
ggplot(COAD.table.1.count, aes(Patient.TOTAL, Patient.Silent)) + geom_point()  + stat_smooth(method="lm")

OV.table.1.count<-OV.table.1$table.1[,list(Patient.Missense=sum(Missense), 
    Patient.Silent=sum(Silent), Patient.TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
OV.table.1.melt<-melt(OV.table.1.count, id.vars=c("PATIENT", "Patient.TOTAL"))
ggplot(OV.table.1.melt, aes(Patient.TOTAL, value, colour=variable)) + geom_point()  + stat_smooth(method="lm") + theme(legend.position="bottom")
    ggplot(OV.table.1.count, aes(Patient.TOTAL, Patient.Silent)) + geom_point()  + stat_smooth(method="lm")

###Can we model this at the gene level? - RESULT: MAYBE, NEED TO BORROW INFO FROM OTHER GENES

#Try with a few genes
ttn.table.1.count<-brca.table.1$table.1[Hugo_Symbol=="TTN",][,list(Missense=sum(Missense), Silent=sum(Silent),TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
ttn.table.1.melt<-melt(ttn.table.1.count, id.vars=c("PATIENT", "TOTAL"))
ggplot(ttn.table.1.melt, aes(TOTAL, value, colour=variable)) + geom_point() + stat_smooth(method="lm")
ggplot(ttn.table.1.melt, aes(variable, value, colour=variable)) + geom_boxplot() 

tp53.table.1.count<-brca.table.1$table.1[Hugo_Symbol=="TP53",][,list(Missense=sum(Missense), Silent=sum(Silent),TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]
tp53.table.1.melt<-melt(tp53.table.1.count, id.vars=c("PATIENT", "TOTAL"))
ggplot(tp53.table.1.melt, aes(TOTAL, value, colour=variable)) + geom_point() + stat_smooth(method="lm")

###But then, can we model it at the patient level?
ggplot(brca.table.1.count, aes(Patient.TOTAL, Patient.Silent)) + geom_point()  + stat_smooth(method="lm")
brca.table.1.count
brca.lm<-lm(Patient.Silent~Patient.TOTAL, data=brca.table.1.count)
coef(brca.lm)[[1]]
summary(brca.lm)    

count.exp<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST.1.5), as.data.frame(brca.table.1.count)))
count.exp.model<-lm(Patient.Missense~Patient.TOTAL+exp(DIST.TO.NORMAL), data=count.exp)
summary(count.exp.model)
ggplot(count.exp, aes(DIST.TO.NORMAL, Patient.Missense)) + geom_point() + scale_y_log10()
ggplot(count.exp, aes(Patient.Silent, DIST.TO.NORMAL)) + geom_point() + scale_x_log10()

summary(lm(Patient.Missense~Patient.Silent+DIST.TO.NORMAL, data=count.exp))

count.exp.melt<-melt(count.exp, id.vars=c("PATIENT","Patient.TOTAL"))
ggplot(count.exp.melt, aes(Patient.TOTAL, value, colour=))

count.exp.melt<-as.data.table(merge(as.data.frame(BRCA.CANCER.DIST.1.5), as.data.frame(brca.table.1.melt)))
ggplot(count.exp.melt, aes(Patient.TOTAL, DIST.TO.NORMAL)) + geom_point() + scale_x_log10()
ggplot(count.exp.melt, aes(value, DIST.TO.NORMAL, colour=variable)) + geom_point() + scale_x_log10()
