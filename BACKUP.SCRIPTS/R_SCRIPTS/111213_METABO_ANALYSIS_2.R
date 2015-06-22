#METABO ANALYSIS #2
#11/12/13 - FROM 153.py
#Comparisson of synonymous versus non-synonymous mutations

library(reshape)

#LOAD NON-SYNONYMOUS TABLE
NS_TABLE=read.table("NETWORK/R_ANALYSIS/111213_PATIENT_NS_MUTATIONS_PER_AA", header=TRUE, sep="\t")
NS_DATA=as.data.frame(NS_TABLE)
CANCERS.ONLY<-as.vector(unique(NS_DATA$CANCER)); CANCERS.ONLY

#LOAD SYNOYMOUS TABLE
S_TABLE=read.table("NETWORK/R_ANALYSIS/111213_PATIENT_S_MUTATIONS_PER_AA", header=TRUE, sep="\t")
S_DATA=as.data.frame(S_TABLE)

##PROCESS
#PLOT non-synonymous distribution
head(NS_DATA)
MELTED_NS_DATA<-melt(NS_DATA[,c(1,2,5,6)], id=c("CANCER", "PATIENT")) #MELTED got rid of number of genes columns
tail(MELTED_NS_DATA)

ggplot(MELTED_NS_DATA, aes(factor(CANCER), value)) + geom_boxplot(aes(fill=factor(variable))) + scale_y_log10() + 
  theme(legend.position="bottom", axis.text.x=element_text(size=rel(3)), 
        axis.text.y=element_text(size=rel(1.5)), plot.title = element_text(size = rel(2))) + 
  ylab("Mutations per amino acid") + labs(title="NON_SYNONYMOUS")

#MANN_WHITNEY for NON-SYNONIMOUS Boxplot
M.WHITNEY.NS<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in as.vector(unique(MELTED_NS_DATA$CANCER))) {
  SAMPLE<-MELTED_NS_DATA[MELTED_NS_DATA$CANCER==cancer,]
  NS.PVALUE<-wilcox.test(SAMPLE[SAMPLE$variable=="METABOLIC_MUTATIONS_PER_AA",]$value, 
                         SAMPLE[SAMPLE$variable=="NON_METABOLIC_MUTATIONS_PER_AA",]$value)$p.value
  M.WHITNEY.NS<-rbind(M.WHITNEY.NS, data.frame(CANCER=cancer, P.VALUE=NS.PVALUE))
}
M.WHITNEY.NS

#PLOT synonymous distribution
head(S_DATA)
MELTED_S_DATA<-melt(S_DATA[,c(1,2,5,6)], id=c("CANCER", "PATIENT"))

ggplot(MELTED_S_DATA, aes(factor(CANCER), value)) + geom_boxplot(aes(fill=factor(variable))) + scale_y_log10() + 
  theme(legend.position="bottom", axis.text.x=element_text(size=rel(3)), 
        axis.text.y=element_text(size=rel(1.5)), plot.title = element_text(size = rel(2)) ) + 
  ylab("Mutations per amino acid")   + labs(title="SYNONYMOUS")

#MANN_WHITNEY for SYNONIMOUS Boxplot
M.WHITNEY.S<-data.frame(CANCER=c(), P.VALUE=c())
for (cancer in as.vector(unique(MELTED_S_DATA$CANCER))) {
  SAMPLE<-MELTED_S_DATA[MELTED_S_DATA$CANCER==cancer,]
  S.PVALUE<-wilcox.test(SAMPLE[SAMPLE$variable=="METABOLIC_MUTATIONS_PER_AA",]$value, 
                         SAMPLE[SAMPLE$variable=="NON_METABOLIC_MUTATIONS_PER_AA",]$value)$p.value
  M.WHITNEY.S<-rbind(M.WHITNEY.S, data.frame(CANCER=cancer, P.VALUE=S.PVALUE))
}
M.WHITNEY.S

#PLOT NM vs M for synonymous mutated genes per patient
head(S_DATA)
ggplot(S_DATA, aes(x=METABOLIC_MUTATED_GENES, y=NON_METABOLIC_MUTATED_GENES, colour=CANCER)) + geom_point() + 
  scale_y_log10() + scale_x_log10() + facet_wrap(~CANCER)

#PEARSON for NM vs M for synonymous mutated genes per patient
head(S_DATA)
PEARSON.S_DATA<-data.frame(CANCER=c(), P.VALUE=c(), RHO=c())
for (cancer in CANCERS.ONLY) {
  SAMPLE<-S_DATA[S_DATA$CANCER==cancer,]
  CORRELATION<-cor.test(SAMPLE$METABOLIC_MUTATED_GENES, SAMPLE$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  PEARSON.S_DATA<-rbind(PEARSON.S_DATA, data.frame(CANCER=cancer, P.VALUE=format.pval(CORRELATION$p.value), RHO=CORRELATION$estimate))
}
PEARSON.S_DATA

#SLOPE for ABOVE pearson correlation - MANUALLY
head(S_DATA)
coef(lm(NON_METABOLIC_MUTATED_GENES~METABOLIC_MUTATED_GENES, S_DATA[S_DATA$CANCER=="COAD",]))

#HYPERGEOMETRIC TEST for synonymous mutations
head(S_DATA)
S_DATA$P.VALUE<-phyper(S_DATA$METABOLIC_MUTATED_GENES, 5113, 21000-5113, S_DATA$NON_METABOLIC_MUTATED_GENES, lower.tail=FALSE)
S_DATA$P.VALUE.CORRECTED<-p.adjust(S_DATA$P.VALUE, method="fdr")
S_DATA$SIGNIFICANT<-S_DATA$P.VALUE.CORRECTED<0.05

#PLOT significant vs non-significant from HYPERGEOMETRIC test results
S_DATA.SPLIT.LISTS<-list()

for (cancer in CANCERS.ONLY) {
  S_DATA.SPLIT.LISTS[[cancer]]$SIG<-subset(S_DATA, CANCER==cancer & SIGNIFICANT==TRUE)
  S_DATA.SPLIT.LISTS[[cancer]]$NON_SIG<-subset(S_DATA, CANCER==cancer & SIGNIFICANT==FALSE)
  
  S_DATA.SPLIT.LISTS[[cancer]]$SIG$TYPE<-"SIG"
  S_DATA.SPLIT.LISTS[[cancer]]$NON_SIG$TYPE<-"NON_SIG"  
}
head(S_DATA.SPLIT.LISTS$LUSC$NON_SIG)

S_DATA.SPLIT.BOX<-data.frame(CANCERS=c(), TOTAL=c(),PATIENTS=c(), TYPE=c())
for (cancer in CANCERS.ONLY) {
  SAMPLE.SIG<-data.frame(CANCERS=cancer, PATIENTS=nrow(S_DATA.SPLIT.LISTS[[cancer]]$SIG), 
                         TOTAL=nrow(S_DATA.SPLIT.LISTS[[cancer]]$SIG) + nrow(S_DATA.SPLIT.LISTS[[cancer]]$NON_SIG), TYPE="SIG")
  SAMPLE.NON_SIG<-data.frame(CANCERS=cancer, PATIENTS=nrow(S_DATA.SPLIT.LISTS[[cancer]]$NON_SIG), 
                          TOTAL=nrow(S_DATA.SPLIT.LISTS[[cancer]]$SIG) + nrow(S_DATA.SPLIT.LISTS[[cancer]]$NON_SIG), TYPE="NON_SIG")
  S_DATA.SPLIT.BOX<-rbind(S_DATA.SPLIT.BOX, SAMPLE.SIG, SAMPLE.NON_SIG)
}
S_DATA.SPLIT.BOX

ggplot(S_DATA.SPLIT.BOX, aes(CANCERS, PATIENTS, fill=TYPE)) +geom_bar(position="dodge") +
  geom_text(aes(label=paste(round(100*PATIENTS/TOTAL,2), "%")), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(title="Significant Number of Patients") + theme(axis.text.x=element_text(size=rel(1.5)))

#PLOT above along with NON-SYNONIMOUS (SPLIT.BOX.MELTED)
colnames(SPLIT.BOX.MELTED)<-c("CANCERS", "TOTAL", "TYPE", "PATIENTS")
S_DATA.SPLIT.BOX$ORIGIN<-"SYNONYMOUS"
SPLIT.BOX.MELTED$ORIGIN<-"NON_SYNONYMOUS"
SYN.NON_SYN.SIG.NON_SIG<-rbind(S_DATA.SPLIT.BOX, SPLIT.BOX.MELTED)
SYN.NON_SYN.SIG.NON_SIG

ggplot(SYN.NON_SYN.SIG.NON_SIG, aes(CANCERS, PATIENTS, fill=TYPE)) +geom_bar(position="dodge") +
  geom_text(aes(label=paste(round(100*PATIENTS/TOTAL,2), "%")), position=position_dodge(width=0.9), vjust=-0.25) + 
  labs(title="Significant Number of Patients") + theme(axis.text.x=element_text(size=rel(1.5))) + facet_grid(ORIGIN~.)

facet_grid(LENGTH.BINS~., scales="free")

#PLOT UNCOUPLED PEARSON for SYNONYMOUS, GET SLOPES!
head(S_DATA)
ggplot(S_DATA, aes(x=METABOLIC_MUTATED_GENES, y=NON_METABOLIC_MUTATED_GENES, col=SIGNIFICANT)) +
  geom_point() + scale_y_log10() + scale_x_log10() + facet_wrap(~CANCER)

S_DATA_UNCOUPLED_PEARSON<-data.frame(CANCER=c(), SIG_P.VALUE=c(), SIG_RHO=c(), NON_SIG_P.VALUE=c(), NON_SIG_RHO=c())
for (cancer in CANCERS.ONLY) {
  
  SAMPLE.SIG<-S_DATA[S_DATA$CANCER==cancer & S_DATA$SIGNIFICANT==TRUE,]
  SIG.TEST<-cor.test(SAMPLE.SIG$METABOLIC_MUTATED_GENES, SAMPLE.SIG$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  
  SAMPLE.NON_SIG<-S_DATA[S_DATA$CANCER==cancer & S_DATA$SIGNIFICANT==FALSE,]
  NON_SIG.TEST<-cor.test(SAMPLE.NON_SIG$METABOLIC_MUTATED_GENES, SAMPLE.NON_SIG$NON_METABOLIC_MUTATED_GENES, method=c("pearson"))
  
  S_DATA_UNCOUPLED_PEARSON<-rbind(S_DATA_UNCOUPLED_PEARSON, 
                                  data.frame(CANCER=cancer, SIG_P.VALUE=format.pval(SIG.TEST$p.value), SIG_RHO=SIG.TEST$estimate,
                                             NON_SIG_P.VALUE=format.pval(NON_SIG.TEST$p.value), NON_SIG_RHO=NON_SIG.TEST$estimate))
}
S_DATA_UNCOUPLED_PEARSON

#GET SLOPES FROM UNCOUPLED PEARSON - MANUALLY
head(S_DATA)
coef(lm(NON_METABOLIC_MUTATED_GENES~METABOLIC_MUTATED_GENES, S_DATA[S_DATA$CANCER=="OV" & S_DATA$SIGNIFICANT==TRUE,]))

#PLOT PEARSON IN TERMS OF MUTATIONS PER AMINO ACID - NON_METABOLIC VS METABOLIC
ggplot(S_DATA, aes(x=METABOLIC_MUTATIONS_PER_AA, y=NON_METABOLIC_MUTATIONS_PER_AA, col=SIGNIFICANT)) +
  geom_point() + scale_y_log10() + scale_x_log10() + facet_wrap(~CANCER)

#PEARSON for above
S_DATA_UNCOUPLED_PEARSON.PER_AA<-data.frame(CANCER=c(), SIG_P.VALUE=c(), SIG_RHO=c(), NON_SIG_P.VALUE=c(), NON_SIG_RHO=c())
for (cancer in CANCERS.ONLY) {
  
  SAMPLE.SIG<-S_DATA[S_DATA$CANCER==cancer & S_DATA$SIGNIFICANT==TRUE,]
  SIG.TEST<-cor.test(SAMPLE.SIG$METABOLIC_MUTATIONS_PER_AA, SAMPLE.SIG$NON_METABOLIC_MUTATIONS_PER_AA, method=c("pearson"))
  
  SAMPLE.NON_SIG<-S_DATA[S_DATA$CANCER==cancer & S_DATA$SIGNIFICANT==FALSE,]
  NON_SIG.TEST<-cor.test(SAMPLE.NON_SIG$METABOLIC_MUTATIONS_PER_AA, SAMPLE.NON_SIG$NON_METABOLIC_MUTATIONS_PER_AA, method=c("pearson"))
  
  S_DATA_UNCOUPLED_PEARSON.PER_AA<-rbind(S_DATA_UNCOUPLED_PEARSON.PER_AA, 
                                  data.frame(CANCER=cancer, SIG_P.VALUE=format.pval(SIG.TEST$p.value), SIG_RHO=SIG.TEST$estimate,
                                             NON_SIG_P.VALUE=format.pval(NON_SIG.TEST$p.value), NON_SIG_RHO=NON_SIG.TEST$estimate))
}
S_DATA_UNCOUPLED_PEARSON.PER_AA

#GET SLOPE for above - SEMI-MANUALLY
for (cancer in CANCERS.ONLY) {
  print (cancer)
  TEST<-coef(lm(NON_METABOLIC_MUTATIONS_PER_AA~METABOLIC_MUTATIONS_PER_AA, S_DATA[S_DATA$CANCER==cancer & S_DATA$SIGNIFICANT==FALSE,]))
  print (TEST)
}

