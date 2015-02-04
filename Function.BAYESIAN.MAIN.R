#####FUNCTION.MAIN.BAYESIAN######
#020315
#Will be used to build main bayesian functions
#NOTES:
#   We are using chen table mid 4
#   Using COSMIC tables containing somatic missense mutations only (~170 genes total as of 020315)

Function.THOUSAND.Prob.Joint<-function(THOUSAND.PHAST.TABLE,  CHEN.TABLE, FEATURES, exp.table,hic.table, biogrid, exon, noise=1, rounded=3, nt.cut=4, rep.cut=4){
  
  require(data.table)
  require(car)
  
  #Unique tables
  setkey(THOUSAND.PHAST.TABLE)
  THOUSAND.PHAST.TABLE<-unique(THOUSAND.PHAST.TABLE)
  
  #Filter thousand table
  main.table<-THOUSAND.PHAST.TABLE[EXON==TRUE,]
  
  ####INTRODUCE FEATURES#####
  #Introduce rep times and cuts
  if (length(intersect(FEATURES, c("REP.TIME","REP.CLASS")))>0){
    CHEN.TABLE$REP.CLASS<-cut(CHEN.TABLE$REP.TIME, quantile(CHEN.TABLE$REP.TIME, seq(0,1,1/rep.cut)), include.lowest=T, labels=letters[1:rep.cut])
    main.table<-merge(main.table, CHEN.TABLE, by="Hugo_Symbol")
    main.table$REP.TIME<-round(main.table$REP.TIME, rounded)  
  }
  
  #Prep classes - Recode ref.alt to reflect binomial chance on either strand
  if ("REF.ALT" %in% FEATURES){
    main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
    main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')  
  }
  
  #Introduce expression info
  #main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  #main.table$logFC<-abs(main.table$logFC)
  #main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce PPI info
  if (length(intersect(FEATURES, c("degree","DEGREE.CUT")))>0){
    main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  }
  
  #Introduce chromatin open state info
  if (length(intersect(FEATURES, c("hic","HIC.CUT")))>0){
    main.table<-merge(main.table, hic.table,by="Hugo_Symbol")
  }
  
  #Introduce exon info - Modified nt.cuts
  if (length(intersect(FEATURES, c("NT.CUT","nt.length")))>0){
    exon$NT.CUT<-cut(exon$nt.length, quantile(exon$nt.length, seq(0,1,1/nt.cut)), include.lowest=T, labels=letters[1:nt.cut])
    main.table<-merge(main.table, exon, by="Hugo_Symbol")  
  }
  
  #Modified joint
  n.maf<-sum(main.table$MAF)
  #main.table$PHAST.CUT<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
  #   main.table[,TYPE.P:=sum(MAF)/n.maf, by="TYPE"]
  #   main.table[,REF.ALT.P:=sum(MAF)/n.maf, by= "REF.ALT"]
  #   main.table[,REP.CLASS.P:=sum(MAF)/n.maf, by="REP.CLASS"]
  #   main.table[,REP.TIME.P:=sum(MAF)/n.maf, by="REP.TIME"]
  #   main.table[,Chrom.P:=sum(MAF)/n.maf, by="Chrom"]
  #   main.table[,DEGREE.CUT.P:=sum(MAF)/n.maf,by="DEGREE.CUT"]
  #   main.table[,Hugo.P:=sum(MAF)/n.maf, by="Hugo_Symbol"]
  #   main.table$THOUSAND.PROB<-main.table$TYPE.P * main.table$REF.ALT.P   *main.table$Hugo.P * main.table$REP.TIME.P
  
  main.table[,THOUSAND.PROB:=sum(MAF)/n.maf, by=FEATURES]
  
  #Add fudge factor
  main.table[,THOUSAND.FF:=sum(MAF)/n.maf,by="Hugo_Symbol"]
  
  #Clean up and return
  main.table<-main.table[,unique(c("Hugo_Symbol","THOUSAND.FF", "THOUSAND.PROB", FEATURES)),with=F]
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.TCGA.Prob.Joint<-function(TCGA.PHAST.TABLE, THOUSAND.PHAST.TABLE, FEATURES, CHEN.TABLE, exp.table, hic.table,biogrid, exon, noise=1, rounded=3, nt.cut=4, rep.cut=4){
  require(data.table)
  
  #Unique tables
  setkey(TCGA.PHAST.TABLE)
  setkey(THOUSAND.PHAST.TABLE)
  TCGA.PHAST.TABLE<-unique(TCGA.PHAST.TABLE)
  THOUSAND.PHAST.TABLE<-unique(THOUSAND.PHAST.TABLE)
  
  #Filter TCGA table - change labels for compatibility with thousand table
  main.tcga<-TCGA.PHAST.TABLE[TYPE %in% c("Missense_Mutation", "Silent"),]
  main.tcga$TYPE<-ifelse(main.tcga$TYPE=="Missense_Mutation", "nonsynonymous SNV", "synonymous SNV")
  
  #Filter THOUSAND table - change labels for compatibility with tcga table
  main.thousand<-THOUSAND.PHAST.TABLE[EXON==TRUE,]
  main.thousand$MUT.FREQ<-main.thousand$MAF*noise
  
  #"Merge" tables
  main.tcga<-main.tcga[,c("Chrom", "Position","Hugo_Symbol", "MUT.FREQ","REF","ALT","TYPE"),with=F]
  main.thousand<-main.thousand[,colnames(main.tcga),with=F]
  main.table<-rbind(main.tcga, main.thousand)
  main.table<-main.table[,list(C.MUT.FREQ=mean(MUT.FREQ)),by=c("Chrom","Position","Hugo_Symbol","REF","ALT","TYPE")]
  
  ####INTRODUCE FEATURES AS NEEDED####
  #Introduce rep times and cuts
  if (length(intersect(FEATURES, c("REP.TIME","REP.CLASS")))>0){
    CHEN.TABLE$REP.CLASS<-cut(CHEN.TABLE$REP.TIME, quantile(CHEN.TABLE$REP.TIME, seq(0,1,1/rep.cut)), include.lowest=T, labels=letters[1:rep.cut])
    main.table<-merge(main.table, CHEN.TABLE, by="Hugo_Symbol")
    main.table$REP.TIME<-round(main.table$REP.TIME, rounded)  
  }
  
  #Prep classes - Recode ref.alt to reflect binomial chance on either strand
  if ("REF.ALT" %in% FEATURES){
    main.table$REF.ALT<-paste(as.vector(main.table$REF), as.vector(main.table$ALT), sep="_")
    main.table$REF.ALT<-recode(main.table$REF.ALT, ' "T_G"="A_C"; "T_C"="A_G"; "T_A"="A_T"; "G_C"="C_G" ; "G_T"="C_A" ; "G_A"="C_T" ')  
  }
  
  #Introduce expression info
  #main.table<-merge(main.table, exp.table[,c("logFC","Hugo_Symbol"),with=F], by="Hugo_Symbol")
  #main.table$logFC<-abs(main.table$logFC)
  #main.table$EXP.CUTS<-cut(main.table$logFC, c(min(exp.table$logFC),quantile(exp.table$logFC, c(0.25,0.5,0.75)), max(exp.table$logFC)), include.lowest=T)
  
  #Introduce PPI info
  if (length(intersect(FEATURES, c("degree","DEGREE.CUT")))>0){
    main.table<-merge(main.table, biogrid, by="Hugo_Symbol")
  }
  
  #Introduce chromatin open state info
  if (length(intersect(FEATURES, c("hic","HIC.CUT")))>0){
    main.table<-merge(main.table, hic.table,by="Hugo_Symbol")
  }
  
  #Introduce exon info - Modified nt.cuts
  if (length(intersect(FEATURES, c("NT.CUT","nt.length")))>0){
    exon$NT.CUT<-cut(exon$nt.length, quantile(exon$nt.length, seq(0,1,1/nt.cut)), include.lowest=T, labels=letters[1:nt.cut])
    main.table<-merge(main.table, exon, by="Hugo_Symbol")  
  }
  
  #Modified joint
  #main.table$PHAST.CUT<-cut(main.table$PHAST.45, c(0,0.0005,0.005,0.02,0.1,0.3,0.8,0.98,0.995,0.999,1.0), include.lowest=T)
  n.maf<-sum(main.table$C.MUT.FREQ)
  
  #   main.table[,TYPE.P:=sum(C.MUT.FREQ)/n.maf, by="TYPE"]
  #   main.table[,REF.ALT.P:=sum(C.MUT.FREQ)/n.maf, by= "REF.ALT"]
  #   main.table[,REP.CLASS.P:=sum(C.MUT.FREQ)/n.maf, by="REP.CLASS"]
  #   main.table[,REP.TIME.P:=sum(C.MUT.FREQ)/n.maf, by="REP.TIME"]
  #   main.table[,Chrom.P:=sum(C.MUT.FREQ)/n.maf, by="Chrom"]
  #   main.table[,DEGREE.CUT.P:=sum(C.MUT.FREQ)/n.maf,by="DEGREE.CUT"]
  #   main.table[,Hugo.P:=sum(C.MUT.FREQ)/n.maf, by="Hugo_Symbol"]
  #   main.table$TCGA.PROB=main.table$TYPE.P * main.table$REF.ALT.P   *main.table$Hugo.P *main.table$REP.TIME.P
  
  main.table[,TCGA.PROB:=sum(C.MUT.FREQ)/n.maf, by=FEATURES]
  
  #Add fudge factor
  main.table[,TCGA.FF:=sum(C.MUT.FREQ)/n.maf,by="Hugo_Symbol"]
  
  #Clean up and return
  setkey(main.table)
  main.table<-unique(main.table)
  return(main.table)
}

Function.Main.Bayes.Joint<-function(thousand.prop.joint, tcga.prob.joint, cancer.prob=0.12, FF=F, FEATURES){
  #Merges cancer and non-cancer probabilities to calculate non-naive bayesian probability per site
  
  library(data.table)
  
  #Merge probabilites (TCGA.PROB, THOUSAND.PROB)
  main.table<-merge(tcga.prob.joint, thousand.prop.joint, by=unique(c("Hugo_Symbol", FEATURES)))
  
  #Calculate bayes prob per site based on features
  #FF ?
  if (FF==T){
    main.table$BAYES.PROB<-(cancer.prob*main.table$TCGA.PROB*main.table$TCGA.FF)/
      (cancer.prob*main.table$TCGA.PROB*main.table$TCGA.FF + (1-cancer.prob)*(main.table$THOUSAND.PROB*main.table$THOUSAND.FF))  
  } else {
    main.table$BAYES.PROB<-(cancer.prob*main.table$TCGA.PROB)/
      (cancer.prob*main.table$TCGA.PROB + (1-cancer.prob)*(main.table$THOUSAND.PROB))  
  }
  
  #Clean up and return
  main.table<-main.table[order(BAYES.PROB, decreasing=T),]
  return(main.table)
}

Function.Bayes.Test<-function(main.table, cosmic.table, cancer.terms){
  #Extract cancer genes from cosmic table related to genes
  
  require (data.table)
  
  #Extract genes of interest
  cancer.genes<-c()
  for (term in cancer.terms){
    cancer.genes<-c(cancer.genes, cosmic.table[grepl(term, TUMOR, ignore.case=T),]$Hugo_Symbol)
  }
  
  #Classify
  main.table$cancer<-main.table$Hugo_Symbol %in% cancer.genes
  main.table$all.cancers<-main.table$Hugo_Symbol %in% cosmic.table$Hugo_Symbol
  
  #Analyze...
  
  #Return
  return (main.table) 
}

Function.bayes.prob.plot<-function(test.bayes.plot) {
  bayes.box.test<-data.table()
  for (prob in c(0.99,0.95,0.9,0.85, 0.80, 0.75, 0.7, 0.6, 0.5,0.4,0.3,0.2,0)){
    n.all.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & all.cancers==TRUE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
    n.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & cancer==TRUE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
    n.non.cancer<-length(unique(as.vector(test.bayes.plot[BAYES.PROB>=prob & all.cancers==FALSE & TYPE=="nonsynonymous SNV",]$Hugo_Symbol)))
    
    bayes.box.test<-rbind(bayes.box.test, data.table(PROB=prob, ALL.CANCER=n.all.cancer, CANCER=n.cancer, NON.CANCER=n.non.cancer))
    print (c(prob, n.all.cancer,n.cancer, n.non.cancer))
  }
  ggplot(melt(bayes.box.test,id.vars="PROB"), aes(PROB, value, colour=variable)) + geom_bar(stat="identity", position="dodge") + theme.format +
    geom_text(aes(label=value), position=position_dodge(width=0.04), vjust=-0.25)   
  
}

Function.Main<-function(THOUSAND.PHAST.TABLE, TCGA.PHAST.TABLE, FEATURES, CHEN.TABLE, hic.table, biogrid, exon, rounded=3, nt.cut=4, rep.cut=4, cancer.prob=0.12, FF=F){
  
  #Construct prior
  prior<-Function.THOUSAND.Prob.Joint(THOUSAND.PHAST.TABLE, CHEN.TABLE, FEATURES, hic.table=hic.table, biogrid=biogrid, exon=exon, rounded=rounded)
  
  #Construct posterior
  posterior<-Function.TCGA.Prob.Joint(TCGA.PHAST.TABLE,THOUSAND.PHAST.TABLE, FEATURES, CHEN.TABLE, hic.table=hic.table, biogrid=biogrid, exon=exon, rounded=rounded)
  
  #Calculate bayesian probabilities
  bayes.table<-Function.Main.Bayes.Joint(prior, posterior, cancer.prob, FF, FEATURES)
  
  #Return
  return(bayes.table)
}

test.bayes<-Function.Main(THOUSAND.PHAST.45[MAF!=0,], TCGA.BRCA, FEATURES, CHEN.REP, hic, biogrid.degree, exon, rounded=2, cancer.prob=0.133, FF=T)
test.bayes.plot<-Function.Bayes.Test(test.bayes, COSMIC, c("UTERINE"))
Function.bayes.prob.plot(test.bayes.plot)


cancer.prob<-data.table(CANCERS=c("BRCA","GBM","SKCM","AML","OV","COAD","HNSC","KIRC","LUAD","LUSC","PRAD","READ","STAD","UCEC"),
                        PROBS=c(0.133,0.10,0.02,0.12,0.05,0.025,0.04,0.11,0.025,0.02,0.11,0.015,0.023,0.007),
                        CALLS=c("breast","gbm,glioblastoma","skcm,skin","aml","ovarian","colorectal","neck","cell renal", "lung","lung","prostate","colorectal","","uterine"))

CANCER.THRESHOLD<-0.8
cancer.main.table<-data.table()
for (x in 1:nrow(cancer.prob)) {
  tcga.table<-get(paste0("TCGA.",cancer.prob[x,1,with=F]))
  prob<-as.numeric(cancer.prob[x,2,with=F])
  cancer.calls<-unlist(strsplit(as.character(cancer.prob[x,3,with=F]),","))
  print (cancer.calls)
  
  test.bayes<-Function.Main(THOUSAND.PHAST.45[MAF!=0,], tcga.table, FEATURES, CHEN.REP, hic, biogrid.degree, exon, rounded=2, cancer.prob=prob, FF=T)
  test.bayes.plot<-Function.Bayes.Test(test.bayes, COSMIC, cancer.calls)
  test.bayes.plot<-test.bayes.plot[BAYES.PROB>=CANCER.THRESHOLD,]
  
  test.bayes.plot$TUMOR<-as.character(cancer.prob[x,1,with=F])
  cancer.main.table<-rbind(cancer.main.table, test.bayes.plot)
}


cancer.main.table<-cancer.main.table[,c("Hugo_Symbol","BAYES.PROB","cancer","TUMOR"),with=F]
setkey(cancer.main.table)
cancer.main.table<-unique(cancer.main.table)
cancer.main.table[order(BAYES.PROB,decreasing=T),]

biogrid.comb[A %in% cancer.main.table[TUMOR=="BRCA",]$Hugo_Symbol,]

