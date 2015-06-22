theme.format<-theme(axis.text.y=element_text(size=rel(2.5)), axis.text.x=element_text(size=rel(2.5)),
                    axis.title.y=element_text(size=22), axis.title.x=element_text(size=22),
                    legend.text = element_text(size = 22))

normalize.vector<-function(x){
  y=(x-min(x))/(max(x)- min(x))
  return(y)
}

Function.paired.t.tests.tailored<-function(vector.a, vector.b, ALT="two.sided", MU=0, var=F, paired=F, type="parametric"){
  #Understand CAREFULLY what you are inputting as MU!!
  #RETURNS LIST!!
  #Default for alternative hypothesis is "two.sided"
  if (type=="parametric"){
    #pair t.test
    if(  (length(unique(vector.a))>1) | (length(unique(vector.b))>1)  ) {
      dummy.test<-t.test(vector.a, vector.b, var.equal=var, alternative=ALT, mu=MU, paired=paired)
    }
  } else if (type=="non.parametric") {
    #Wilcoxon-signed rank test
    dummy.test<-wilcox.test(vector.a, vector.b, var.equal=var, alternative=ALT, mu=MU, paired=paired)
  }
  #Return results of all tests plus number of patients covered by SIGNIFICANCE vector (vector.a) - LENGTH OF VECTOR USED FOR COUNTING COVERAGE!!!!
  return(list(stat=dummy.test$statistic, p.val=dummy.test$p.value,
              PATIENTS.COVERED=length(vector.a)))
}

Function.p.rank.enrichment<-function(table.p, interest.columns, target.genes, normalize=F) {
  #Rank cummulative and non-cummulative enrichment for proteins p in order of decreasing SCORE
  #interest.columns assumes that first column is ID and second column is scores  ie. [KEGG_ID, wilcoxon.stat]
  #Calculates PRECISSION-RECALL values as well for each cummulative count
  require(data.table)
  require(base)
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  #Pre-process tables
  table.p<-table.p[,interest.columns, with=F]
  if (normalize==T) {
    table.p$SCORE<-normalize.vector(table.p[,2, with=F])
    table.p<-table.p[,c(1,3), with=F]
  }
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  #Order j.table by score to calculate cummulative enrichemnt [Hugo_Symbol, SCORE]
  table.p<-table.p[order(SCORE, decreasing=T),]
  #CALCULATIONS
  background.genes<-unique(as.vector(table.p$Hugo_Symbol)) #Genes that we can choose from
  target.in.p<-intersect(background.genes, target.genes) #Targets we could have chosen (White in urn)
  #Breaks for tables 0-1
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.0001)
  PR.BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  #FIRST CUMMULATIVE [BREAKS, CUM.SCORE.PVALS, CUM.SCORE.COVERAGE]
  #Calculate p.value
  CUM.P.VAL<-sapply(PR.BREAKS, function(x)
    phyper(q=length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))-1,
           m=length(target.in.p),
           n=length(background.genes)-length(target.in.p),
           k= length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))),
           lower.tail=F
    ))
  #Calculate coverage
  CUM.SCORE.COVERAGE<-sapply(PR.BREAKS, function (x)
    length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/length(target.in.p) )
  #Calcualte fdr for CUM
  CUM.TABLE<-as.data.table(cbind(PR.BREAKS, CUM.P.VAL, CUM.SCORE.COVERAGE))
  setnames(CUM.TABLE, colnames(CUM.TABLE), c("BREAKS", "CUM.P.VAL", "CUM.SCORE.COVERAGE"))
  CUM.TABLE$CUM.ADJ.P.VAL<-p.adjust(CUM.P.VAL, method="fdr")
  #THEN NON-CUMMULATIVE - [Hugo_Symbol, SCORE, NON.CUM.P.VAL, NON.CUM.SCORE.COVERAGE]
  table.p$NON.CUM.P.VAL<-sapply(as.vector(table.p$Hugo_Symbol), function(z) phyper(q=length(intersect(z, target.in.p))-1,
                                                                                   m=length(target.in.p),
                                                                                   n=length(background.genes)-length(target.in.p),
                                                                                   k=length(unique(z)), lower.tail=F))
  table.p$NON.CUM.SCORE.COVERAGE<-sapply(as.vector(table.p$Hugo_Symbol), function(z) length(intersect(z, target.in.p))/length(target.in.p))
  #Assign to Hugo_Symbols for NON.CUM
  table.p$NON.CUM.ADJ.P.VAL<-p.adjust(table.p$NON.CUM.P.VAL, method="fdr")
  #Lastly add CUM.RECALL and CUM.PRECISION columns to CUM
  CUM.TABLE$CUM.RECALL<-CUM.SCORE.COVERAGE
  CUM.TABLE$CUM.PRECISION<-sapply(PR.BREAKS, function (x)
    length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/ length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))))
  #Return as list
  return(list(CUM.TABLE=CUM.TABLE,NON.CUM.TABLE=table.p))
}

Function.PR.PATIENT.COVERAGE<-function(table.p, interest.columns, CANCER, c.table.1 ){
  require(data.table)
  require(parallel)
  #Process scoring table
  table.p<-table.p[,interest.columns, with=F]
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  #CALCULATIONS
  background.genes<-unique(as.vector(table.p$Hugo_Symbol))
  background.patients<-unique(as.vector(c.table.1$Tumor_Sample_Barcode)) #Patients that we can choose from
  target.patients<-unique(as.vector(c.table.1[TYPE==CANCER & Hugo_Symbol %in% background.genes,]$Tumor_Sample_Barcode)) #Targets we could have chosen (White in urn) - ALL POSITIVES
  #Breaks for tables 0-1
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.0001)
  PR.BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  c.table.1<-copy(c.table.1)
  clusterExport(cl, varlist=c("table.p", "target.patients", "c.table.1"),envir=environment())
  #PR
  SCORE.RECALL<-parSapply(cl, PR.BREAKS, function (x)
    length(unique(as.vector(c.table.1[c.table.1$TYPE==CANCER & c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) / length(target.patients))
  SCORE.PRECISION<-parSapply(cl, PR.BREAKS, function(x)
    length(unique(as.vector(c.table.1[c.table.1$TYPE==CANCER & c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) /
                               length(unique(as.vector(c.table.1[c.table.1$Hugo_Symbol %in% unique(as.vector(table.p[table.p$SCORE>=x,]$Hugo_Symbol)),]$Tumor_Sample_Barcode))) )
  #Stop parallelization
  stopCluster(cl)
  #Clean up and Return
  dummy.return<-as.data.table(do.call(cbind, list(PR.BREAKS, SCORE.RECALL, SCORE.PRECISION)))
  setnames(dummy.return, colnames(dummy.return), c("PR.BREAKS", "RECALL", "PRECISION"))
  return(dummy.return)
}

Function.composite.table.1<-function(tables.1, cnv.tables, cancers) {
  #Takes table.1 and cnv.table from GISTIC to create a combined table for all cancers supplied
  
  require(base)
  require(data.table)
  dummy.list<-list()
  for (type in 1:length(cancers)){
    cancer.cnv<-readRDS(cnv.tables[type])
    setnames(cancer.cnv, colnames(cancer.cnv), c("Tumor_Sample_Barcode", "Hugo_Symbol"))
    cancer.table.1<-readRDS(tables.1[type])
    cancer.table.1<-cancer.table.1$table.1[,1:2,with=F]
    cancer.table.1$Tumor_Sample_Barcode<-sapply(as.character(cancer.table.1$Tumor_Sample_Barcode), function(x) paste0(strsplit(x,"-")[[1]][1:4],collapse=".")  )
    cancer.table.1<-unique(rbind(cancer.table.1, cancer.cnv))
    cancer.table.1$TYPE<-cancers[type]
    dummy.list[[type]]<-cancer.table.1
  }
  dummy.return<-as.data.table(do.call(rbind, dummy.list))
  return(dummy.return)
}
#Example
# composite.tables.1<-Function.composite.table.1(c("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/072214_Table1_COAD.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/072314_Table1_UCEC.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/072214_Table1_GBM.rds"),
#                                               c("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/080714.BRCA.GISTIC.TH.2.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/080714.COAD.GISTIC.TH.2.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/080714.GBM.GISTIC.TH.2.rds",
#                                                 "PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/080714.UCEC.GISTIC.TH.2.rds"),
#                                               c("BRCA", "COAD","UCEC","GBM"))
# composite.tables.1

Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  require(plyr)
  require(limma)
  require(edgeR)
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  #Combine matrices
  #shared.genes = intersect(rownames(normal.matrix), rownames(cancer.matrix))
  #normal = normal.matrix[shared.genes, ]
  #tumor = cancer.matrix[shared.genes, ]
  #G.all = as.matrix(cbind(normal, tumor))
  #G.all=G.all[complete.cases(G.all),]
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 20 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  #Return topTable
  return(all.G.fit)
}
       
Function.path.rank.enrichment<-function(table.j, interest.columns, pathway.uniprot.file, uniprot.mapping.file, target.genes, normalize=F) {
  #Rank cummulative and non-cummulative enrichment for metabolites j in order of decreasing SCORE
  #interest.columns assumes that first column is KEGG_ID and second column is scores  ie. [KEGG_ID, wilcoxon.stat]
  require(data.table)
  require(base)
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  #Pre-process tables
  table.j<-table.j[,interest.columns, with=F]
  if (normalize==T) {
    table.j$SCORE<-normalize.vector(table.j[,2, with=F])
    table.j<-table.j[,c(1,3), with=F]
  }
  setnames(table.j, colnames(table.j), c("PATH", "SCORE"))
  #Process pathway files
  uniprot.mapping<-fread(uniprot.mapping.file, header=T, sep="\t", drop=c(2:4, 6:7)) #[UNIPROT, Hugo_Symbol]
  setnames(uniprot.mapping, colnames(uniprot.mapping), c("UNIPROT", "GENES"))
  uniprot.mapping$Hugo_Symbol<-sapply(uniprot.mapping$GENES, function(x) strsplit(x," ")[[1]][1])
  uniprot.mapping$GENES<-NULL
  pathway.uniprot<-as.data.table(read.csv(pathway.uniprot.file, header=F, sep="\t")) #[UNIPROT, PATH]
  pathway.uniprot<-pathway.uniprot[,c(1,4), with=F]
  setnames(pathway.uniprot, colnames(pathway.uniprot), c("UNIPROT", "PATH"))
  dummy.table.2<-as.data.table(merge(as.data.frame(uniprot.mapping), as.data.frame(pathway.uniprot), by="UNIPROT")) #[PATH, Hugo_Symbol]
  dummy.table.2<-unique(dummy.table.2[,c("PATH", "Hugo_Symbol"), with=F])
  #Order j.table by score to calculate cummulative enrichemnt [PATH, SCORE]
  table.j<-table.j[order(SCORE, decreasing=T),]
  #CALCULATIONS
  background.genes<-unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j$PATH),]$Hugo_Symbol)) #Genes that we can choose from (All in urn)
  target.in.j<-intersect(background.genes, target.genes) #Targets we could have chosen (White in urn)
  #FIRST CUMMULATIVE
  #Calculate p.value
  CUM.SCORE.PVALS<-sapply(1:nrow(table.j), function(x)
    phyper(q=length(intersect(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol)),target.in.j))-1,
           m=length(target.in.j),
           n=length(background.genes)-length(target.in.j),
           k= length(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol))),
           lower.tail=F
    ))
  #Calculate coverage
  CUM.SCORE.COVERAGE<-sapply(1:nrow(table.j), function (x)
    length(intersect(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol)),target.in.j))/length(target.in.j) )
  #THEN NON-CUMMULATIVE
  table.non.cum<-as.data.table(merge(as.data.frame(table.j), as.data.frame(dummy.table.2), by="PATH")) #[PATH, SCORE, Hugo_Symbol]
  table.non.cum<-table.non.cum[,list(NON.CUM.P.VAL=phyper(q=length(intersect(Hugo_Symbol, target.in.j))-1,
                                                          m=length(target.in.j),
                                                          n=length(background.genes)-length(target.in.j),
                                                          k=length(unique(Hugo_Symbol)), lower.tail=F),
                                     NON.CUM.SCORE.COVERAGE=length(intersect(Hugo_Symbol, target.in.j))/length(target.in.j)),
                               by=c("PATH", "SCORE")]
  #Assign to PATHs
  table.j$CUM.SCORE.COVERAGE<-CUM.SCORE.COVERAGE
  table.j$CUM.P.VAL<-CUM.SCORE.PVALS
  table.j$CUM.ADJ.P.VAL<-p.adjust(table.j$CUM.P.VAL, method="fdr")
  table.j<-as.data.table(merge(as.data.frame(table.j), as.data.frame(table.non.cum), by=c("PATH", "SCORE")))
  table.j$NON.CUM.ADJ.P.VAL<-p.adjust(table.j$NON.CUM.P.VAL, method="fdr")
  #Return
  return(table.j)
}         
         
Function.j.rank.enrich.plot<-function(Table.enrich.result){
  #Plots results of Function.j.rank.enrichment()
  require(ggplot2)
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # Multiple plot function
    #
    # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    # - cols:   Number of columns in layout
    # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
    #
    # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    # then plot 1 will go in the upper left, 2 will go in the upper right, and
    # 3 will go all the way across the bottom.
    #
    require(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
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
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  p1<-ggplot(Table.enrich.result$CUM.TABLE, aes(x=BREAKS, y=CUM.SCORE.COVERAGE)) +geom_point(size=4) +theme.format
  p3<-ggplot(Table.enrich.result$CUM.TABLE, aes(x=BREAKS, y=-log(CUM.ADJ.P.VAL))) + geom_point(size=4) +theme.format +
    geom_hline(yintercept=-log(0.05), colour="red")
  p2<-ggplot(Table.enrich.result$NON.CUM.TABLE, aes(x=SCORE, y=NON.CUM.SCORE.COVERAGE)) +geom_point(size=4) +theme.format
  p4<-ggplot(Table.enrich.result$NON.CUM.TABLE, aes(x=SCORE, y=-log(NON.CUM.ADJ.P.VAL))) + geom_point(size=4) +theme.format +
    geom_hline(yintercept=-log(0.05), colour="red")
  multiplot(p1,p2,p3,p4, cols=2)
}

Function.p.rank.patient.enrichment<-function(table.p, interest.columns, processed.table.1.rds, normalize=F) {
  #Ranks cummulative and non-cumulative coverage per protein p by decreasing SCORE
  #Interest columns assumes that first column is ID and second column is scores ie. [KEGG_ID, wilcoxon.stat]
  #Calculates PRECISSION-RECALL values as well for each cummulative count
  
  require(data.table)
  require(base)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-process tables
  table.p<-table.p[,interest.columns, with=F]
  if (normalize==T) {
    table.p$SCORE<-normalize.vector(table.p[,2, with=F])
    table.p<-table.p[,c(1,3), with=F]
  }
  
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  table.1<-readRDS(processed.table.1.rds)
  table.1<-table.1$table.1[, 1:2, with=F] #[Tumor_Sample_Barcode, Hugo_Symbol]
  
  #Order j.table by score to calculate cummulative enrichemnt [Hugo_Symbol, SCORE]
  table.p<-table.p[order(SCORE, decreasing=T),]
  
  #CALCULATIONS
  starting.patients<-unique(as.vector(table.1$Tumor_Sample_Barcode)) #Patients fed to model
  target.patients<-unique(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p$Hugo_Symbol),]$Tumor_Sample_Barcode)) #Patients we can actually cover
  
  #Calculate Recall
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.01)
  BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  CUM.RECALL<-sapply(BREAKS, function (x)
    length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[SCORE>=x,]$Hugo_Symbol),]$Tumor_Sample_Barcode),target.patients))/length(target.patients))
  
  #Calculate coverages
  CUM.PATIENT.COVERAGE<-sapply(1:nrow(table.p), function (x)
    length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[1:x,]$Hugo_Symbol),]$Tumor_Sample_Barcode), target.patients))/length(target.patients))
  NON.CUM.PATIENT.COVERAGE<-sapply(as.vector(table.p$Hugo_Symbol), function(z)
    length(intersect(as.vector(table.1[Hugo_Symbol==z,]$Tumor_Sample_Barcode), target.patients)) /length(target.patients))
  
  #Calculate precision - ONLY RECALL FOR NOW, MAY HAVE TO USE 1000G or OTHER CANCERS COMBINED FOR THIS
  #CUM.PRECISION<-sapply(BREAKS, function (x)
  #  length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[SCORE>=x,]$Hugo_Symbol)]$Tumor_Sample_Barcode),target.patients))/
  #  length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/ length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))))
  
  #Construct tables
  CUM.TABLE<-data.table(BREAKS=BREAKS, PATIENT.RECALL=CUM.RECALL)
  table.p$CUM.PATIENT.COVERAGE<-CUM.PATIENT.COVERAGE
  table.p$NON.CUM.PATIENT.COVERAGE<-NON.CUM.PATIENT.COVERAGE
  
  #Return as list
  return(list(PR.TABLE=CUM.TABLE,COVERAGE.TABLE=table.p, STARTING.PATIENTS=starting.patients, PATIENTS.IN.MODEL=target.patients))
}
#Example
#BRCA.v.p.filtered.2000.enrich.patients<-Function.p.rank.patient.enrichment(BRCA.v.p.filtered.2000, c("V1","V2"),
#"PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds", normalize=F)

Function.all.t.tests<-function(vector.a, vector.b, ALT="two.sided", MU=0){
  #Understand CAREFULLY what you are inputting as MU!!
  #RETURNS LIST!!
  #Default for alternative hypothesis is "two.sided"
  
  #student's t.test
  dummy.t.test<-t.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU)
  
  #pair t.test
  dummy.paired.t.test<-t.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU, paired=T)
  
  #Mann-whitney U-test
  dummy.mann.whitney<-wilcox.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU)
  
  #Wilcoxon-signed rank test
  dummy.wilcoxon<-wilcox.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU, paired=T)
  
  #Return results of all tests
  return(list(student.t.stat=dummy.t.test$statistic, student.p.val=dummy.t.test$p.value,
              paired.t.stat=dummy.paired.t.test$statistic, paired.t.p.val=dummy.paired.t.test$p.value,
              mann.whintey.u.stat=dummy.mann.whitney$statistic, mann.whitney.u.p.val=dummy.mann.whitney$p.value,
              wilcoxon.stat=dummy.wilcoxon$statistic, wilcoxon.p.val=dummy.wilcoxon$p.value))
}

Function.RNAseq.Matrices.Normalization<-function(normal.matrix, cancer.matrix, rm.batch.effect=T) {
  #Normalize RNAseq matrices for posterior differential expression analysis and returns fully normalized expression matrix 
  #+ list of normal and cancer patients
  
  require(plyr)
  require(limma)
  require(edgeR)
  require(sva)
  
  #Get patient cohorts
  G1.patients<-colnames(cancer.matrix)
  G0.patients<-colnames(normal.matrix)
  
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 3 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  
  if (rm.batch.effect==T) {
    print ("performing batch effect correction")
    
    #Obtain batch normalization parameters
    batch.design=data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
    batch.mod = model.matrix(~ G, batch.design) #What we have
    batch.mod0 = model.matrix(~ 1, batch.design) #Null
    batch.ss = svaseq(G.all, batch.mod, batch.mod0) #Apply function to get batch effect parameters
    print ("done correcting for batch effects")
  }
  
  #Quantile normalization
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Need to correct for batch effects?
  if (rm.batch.effect==T){
    #Obtain matrices corrected from batch effect
    f.ss = lmFit(G.all, model.matrix(~G + batch.ss$sv, batch.design)) #Just like regular fit but added batch effect normalization
    ss.removed = G.all$E - f.ss$coefficients[, -c(1, 2)] %*% t(f.ss$design[, -c(1, 2)])
  } else {
    ss.removed=G.all$E
  }
  
  #Return combined matrix fully normalized (quantile, log cpm converated + batch effect) GENESxPATIENTS + list of cancer and normal patients
  ss.removed<-list(combined.matrices=ss.removed, cancer.patients=G1.patients, normal.patients=G0.patients)
  return(ss.removed)
}
#Example
#batch.test<-Function.RNAseq.Matrices.Normalization(cancer.matrices$normal, cancer.matrices$tumor)         
#saveRDS(batch.test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082614.CANCER.MATRICES.NORMALIZED.OBJ.rds")

Function.RNAseq.Differential.Expression.V2<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  
  require(limma)
  require(edgeR)
  require(data.table)
  
  #Get target combined matrix
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  target.combined.matrix<-normalized.matrices.object$combined.matrices[,c(cancer.samples.in.matrix,normal.samples)]
  
  #Get design matrix
  G.design.matrix<-data.frame(G=c(rep("G1", length(cancer.samples.in.matrix)), rep("G0", length(normal.samples))))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Perform differential expression
  G.fit = lmFit(target.combined.matrix, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Convert to data.table
  if ("ID" %in% colnames(all.G.fit)) {
    all.G.fit<-as.data.table(all.G.fit)
  } else {
    all.G.fit$ID<-rownames(all.G.fit)
    all.G.fit<-as.data.table(all.G.fit)
  }
  
  #Return topTable
  return(all.G.fit)
}
         
Function.v.p.Version.3<-function(normalized.matrices.object, table1plus) {
  #Takes object from Function.RNAseq.Matrices.Normalization() and Table.1.Plus to obtain v(p) object
  #This will produce a diff.expression list per gene - Not to be confused with single table object per v(p) as past versions
  #Requires:
  #   Function.RNAseq.Differential.Expression.V2()
  #COMPUTATIONAL INTENSIVE!!!!
  
  require(parallel)
  require(reshape2)
  require(data.table)
  require(limma)
  require(edgeR)
  
  #Filter table 1 plus for genes that have greater than 5 patients
  patient.coverage<-table1plus[,list(size=length(PATIENT)), by="Hugo_Symbol"]
  patient.coverage<-patient.coverage[size>5,]
  table1plus<-table1plus[Hugo_Symbol %in% as.vector(unique(patient.coverage$Hugo_Symbol))]
  
  #Break table.1 into lists per gene
  table.1.lists<-split(table1plus, table1plus$Hugo_Symbol, drop=T)
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression.V2", "table.1.lists", "normalized.matrices.object"),envir=environment())
  
  #v(p) list object
  v.p.object<-parLapply(cl, names(table.1.lists),
                        function(x) Function.RNAseq.Differential.Expression.V2(normalized.matrices.object, as.vector(table.1.lists[[x]]$PATIENT) ) )
  names(v.p.object)<-names(table.1.lists)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(v.p.object)
}         
         
#Example:
#v.p.object.test<-Function.v.p.Version.3(batch.test, Table.1.PLUS[Hugo_Symbol %in% c("A1CF", "TP53")])
#saveRDS(v.p.object.test, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/082714.BRCA.v3.tables.list.rds")         

Function.v.p.Version.3.NULL<-function(normalized.matrices.object, null.patient.count.vector) {
  #Obtain Null lists for 100x random sampling from cancer patients in normalized matrices object
  #COMPUTATIONALLY INTENSIVE!!!!
  
  require(limma)
  require(data.table)
  require(data.table)
  require(reshape2)
  require(parallel)
  
  #Make sure patient null vector includes only samples greater than 5
  null.vector<-null.patient.count.vector[null.patient.count.vector>5]
  
  #All cancer patients to sample from
  cancer.patients<-normalized.matrices.object$cancer.patients
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.RNAseq.Differential.Expression.V2", "null.vector", "normalized.matrices.object", "cancer.patients"),envir=environment())
  
  #Build NULL distributions
  Null.dist<-parLapply(cl, null.vector, function(x) replicate(100, Function.RNAseq.Differential.Expression.V2(normalized.matrices.object, sample(cancer.patients,x)),simplify=F))
  names(Null.dist)<-null.vector
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(Null.dist)
}

Function.v.p.prediction<-function(v.p.object, cancer.diff.exp.table) {
  
  #Uses object from Function.v.p.Version.3() to obtain a coefficient per gene to determine how close it is to predicted diff gene for all thresholds in main cancer table
  require(data.table)
  require(reshape2)
  require(plyr)
  require(parallel)
  
  #Break cancer.diff.table into vectors of genes
  diff.1.5<-as.vector(cancer.diff.exp.table[abs(logFC)>1.5,]$ID)
  diff.2.0<-as.vector(cancer.diff.exp.table[abs(logFC)>2.0,]$ID)
  diff.2.5<-as.vector(cancer.diff.exp.table[abs(logFC)>2.5,]$ID)
  diff.3.0<-as.vector(cancer.diff.exp.table[abs(logFC)>3.0,]$ID)
  diff.3.5<-as.vector(cancer.diff.exp.table[abs(logFC)>3.5,]$ID)
  diff.4.0<-as.vector(cancer.diff.exp.table[abs(logFC)>4.0,]$ID)
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("v.p.object", "diff.1.5", "diff.2.0", "diff.2.5", "diff.3.0", "diff.3.5", "diff.4.0"),envir=environment())
  
  #Run
  v.p.pred.table<-parSapply(cl, v.p.object, function(x) {
    y.1.5<-as.vector(x[abs(x$logFC)>1.5,]$ID)
    y.2.0<-as.vector(x[abs(x$logFC)>2.0,]$ID)
    y.2.5<-as.vector(x[abs(x$logFC)>2.5,]$ID)
    y.3.0<-as.vector(x[abs(x$logFC)>3.0,]$ID)
    y.3.5<-as.vector(x[abs(x$logFC)>3.5,]$ID)
    y.4.0<-as.vector(x[abs(x$logFC)>4.0,]$ID)
    
    #How many of the ones we predict are right (prediction)
    p.1.5<-length(intersect(y.1.5, diff.1.5))/length(y.1.5)
    p.2.0<-length(intersect(y.2.0, diff.2.0))/length(y.2.0)
    p.2.5<-length(intersect(y.2.5, diff.2.5))/length(y.2.5)
    p.3.0<-length(intersect(y.3.0, diff.3.0))/length(y.3.0)
    p.3.5<-length(intersect(y.3.5, diff.3.5))/length(y.3.5)
    p.4.0<-length(intersect(y.4.0, diff.4.0))/length(y.4.0)
    
    #How many of the target we get right (coverage)
    c.1.5<-length(intersect(y.1.5, diff.1.5))/length(diff.1.5)
    c.2.0<-length(intersect(y.2.0, diff.2.0))/length(diff.2.0)
    c.2.5<-length(intersect(y.2.5, diff.2.5))/length(diff.2.5)
    c.3.0<-length(intersect(y.3.0, diff.3.0))/length(diff.3.0)
    c.3.5<-length(intersect(y.3.5, diff.3.5))/length(diff.3.5)
    c.4.0<-length(intersect(y.4.0, diff.4.0))/length(diff.4.0)
    
    #Return
    return(c(p.1.5, p.2.0, p.2.5, p.3.0, p.3.5, p.4.0, c.1.5, c.2.0, c.2.5, c.3.0,c.3.5,c.4.0))
    
  } ,USE.NAMES=T)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Clean up
  v.p.pred.table<-t(v.p.pred.table)
  v.p.pred.table<-as.data.table(v.p.pred.table, keep.rownames=T)
  setnames(v.p.pred.table, c("Hugo_Symbol", "p1.5","p2.0","p2.5","p3.0","p3.5","p4.0", "c1.5","c2.0","c2.5","c3.0","c3.5","c4.0"))
  
  #Return
  return(v.p.pred.table)
}

Function.v.p.corr.matrix.obj<-function(normalized.matrices.object, table1plus, hk.genes) {
  #Produces correlation matrix per causal gene in table1plus after removal of house-keeping genes (low noise) from expression matrix
  #Takes object from Function.RNAseq.Differential.Expression.V2 and Table.1.Plus to obtain v(p) correlation object
  #This will produce:
  #   A correlation matrix per v(p) gene - OR-
  #   If modified in code, will produce a mean variance of all gene acroos all patients
  
  require(parallel)
  require(reshape2)
  require(data.table)
  require(Hmisc)
  
  #Filter table 1 plus for genes that have greater than 5 patients and are in our gene expression matrix only
  table1plus<-table1plus[PATIENT %in% normalized.matrices.object$cancer.patients,]
  patient.coverage<-table1plus[,list(size=length(PATIENT)), by="Hugo_Symbol"]
  patient.coverage<-patient.coverage[size>5,]
  table1plus<-table1plus[Hugo_Symbol %in% as.vector(unique(patient.coverage$Hugo_Symbol)),]
  
  #Break table.1 into lists per gene
  table.1.lists<-split(table1plus, table1plus$Hugo_Symbol, drop=T)
  
  #Get normalized matrix for cancer patient only
  cancer.normalized.matrix<-normalized.matrices.object$combined.matrices[, normalized.matrices.object$cancer.patients]
  
  #Remove house keeping genes from cancer expression matrix
  matrix.genes<-rownames(cancer.normalized.matrix)
  non.hk.genes<-matrix.genes[!(matrix.genes %in% hk.genes)]
  cancer.normalized.matrix<-cancer.normalized.matrix[non.hk.genes, ]
  print(dim(cancer.normalized.matrix))
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("table.1.lists", "cancer.normalized.matrix", "rcorr"),envir=environment())
  
  #Obtain correlation matrix per gene - IF WANTED!
  #cancer.corr.matrix.object<-parLapply(cl, table.1.lists, function(x) rcorr(cancer.normalized.matrix[, as.vector(x$PATIENT)] , type="spearman")$r )
  #names(cancer.corr.matrix.object)<-names(table.1.lists)
  
  # -- OR--, mean variance across genes for all patietns affected by gene v(p)
  cancer.corr.matrix.object<-parSapply(cl, table.1.lists, function(x) mean(apply(cancer.normalized.matrix[,as.vector(x$PATIENT)], 1, var)) ,USE.NAMES=T)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(cancer.corr.matrix.object)
}

Function.v.p.table.V2<-function(v.p.list.object){
  #Takes output list object from Function.v.p.Version.3() and produces classical v(p) table per gene
  
  require(data.table)
  require(parallel)
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("v.p.list.object"),envir=environment())
  
  #Run
  dummy.return<-parSapply(cl, v.p.list.object, function(x) mean(as.vector(x$adj.P.Val)<0.05))
  dummy.return<-as.data.table(dummy.return, keep.rownames=T)
  setnames(dummy.return, c("Hugo_Symbol", "v.PROTEIN"))
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return
  return(dummy.return)
}

Function.v.p.prediction.NULL<-function(v.p.null.object, cancer.diff.exp.table) {
  
  #Uses object from Function.v.p.Version.3() to obtain a coefficient per gene to determine how close it is to predicted diff gene for all thresholds in main cancer table
  require(data.table)
  require(reshape2)
  require(plyr)
  #require(parallel)
  
  #Break cancer.diff.table into vectors of genes
  diff.1.5<-as.vector(cancer.diff.exp.table[abs(logFC)>1.5,]$ID)
  diff.2.0<-as.vector(cancer.diff.exp.table[abs(logFC)>2.0,]$ID)
  diff.2.5<-as.vector(cancer.diff.exp.table[abs(logFC)>2.5,]$ID)
  diff.3.0<-as.vector(cancer.diff.exp.table[abs(logFC)>3.0,]$ID)
  diff.3.5<-as.vector(cancer.diff.exp.table[abs(logFC)>3.5,]$ID)
  diff.4.0<-as.vector(cancer.diff.exp.table[abs(logFC)>4.0,]$ID)
  
  #Set up parallelization
  #nodes<-detectCores()
  #cl<-makeCluster(nodes)
  #setDefaultCluster(cl)
  #clusterExport(cl, varlist=c("v.p.null.object", "diff.1.5", "diff.2.0", "diff.2.5", "diff.3.0", "diff.3.5", "diff.4.0",
  #                            "as.data.table","setnames"),envir=environment())
  
  #Run
  v.p.pred.table<-lapply(v.p.null.object, function(x) {
    
    internal<-sapply(x, function(y)  {
    
      #Get thresholded genes per each distribution in sample size sampling
      y.1.5<-as.vector(y[abs(y$logFC)>1.5,]$ID)
      y.2.0<-as.vector(y[abs(y$logFC)>2.0,]$ID)
      y.2.5<-as.vector(y[abs(y$logFC)>2.5,]$ID)
      y.3.0<-as.vector(y[abs(y$logFC)>3.0,]$ID)
      y.3.5<-as.vector(y[abs(y$logFC)>3.5,]$ID)
      y.4.0<-as.vector(y[abs(y$logFC)>4.0,]$ID)
      
      #Calculate
      a1.5<-length(intersect(y.1.5, diff.1.5))/length(y.1.5) + length(intersect(y.1.5, diff.1.5))/length(diff.1.5)
      a2.0<-length(intersect(y.2.0, diff.2.0))/length(y.2.0) + length(intersect(y.2.0, diff.2.0))/length(diff.2.0)
      a2.5<-length(intersect(y.2.5, diff.2.5))/length(y.2.5) + length(intersect(y.2.5, diff.2.5))/length(diff.2.5)
      a3.0<-length(intersect(y.3.0, diff.3.0))/length(y.3.0) + length(intersect(y.3.0, diff.3.0))/length(diff.3.0)
      a3.5<-length(intersect(y.3.5, diff.3.5))/length(y.3.5) + length(intersect(y.3.5, diff.3.5))/length(diff.3.5)
      a4.0<-length(intersect(y.4.0, diff.4.0))/length(y.4.0) + length(intersect(y.4.0, diff.4.0))/length(diff.4.0)
      
      #Return
      return(c(a1.5,a2.0,a2.5,a3.0,a3.5,a4.0))
    })
    internal<-as.data.table(t(internal))
    setnames(internal, c("1.5","2.0","2.5","3.0","3.5","4.0"))
    
    #Return internal
    return(internal)
  })
  
  #Stop parallelization
  #stopCluster(cl)
  
  #Clean up
  names(v.p.pred.table)<-names(v.p.null.object)
  
  #Return function output
  return(v.p.pred.table)
}

Function.logFC.Median.Differential.Expression<-function(normalized.matrices.object, target.cancer.samples) {
  #Produces logFC.median table per gene between cancer and normal matrices found in matrices.object
  
  require(data.table)
  
  #Get target cancer and normal matrices
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  
  target.normal.matrix<-normalized.matrices.object$combined.matrices[,normal.samples]
  target.cancer.matrix<-normalized.matrices.object$combined.matrices[,cancer.samples.in.matrix]
  
  #Perform median logFC differential expression
  cancer.medians<-apply(target.cancer.matrix,1,median)
  normal.medians<-apply(target.normal.matrix,1,median)
  diff.medians<-cancer.medians-normal.medians
  
  #Clean up
  dummy.return<-data.frame(Hugo_Symbol=rownames(target.cancer.matrix), logFC.Median=diff.medians)
  dummy.return<-as.data.table(dummy.return)
  dummy.return<-dummy.return[order(abs(logFC.Median),decreasing=T),]
  
  #Return topTable
  return(dummy.return)
}

Function.v.p.Version.4.NULL<-function(normalized.matrices.object, patient.coverage.vector) {
  #Constructs null tables in the form of a data table for all size of PATIENT.COVERAGE using median.logFC function
  #Requires:
  #   Function.logFC.Median.Differential.Expression()
  
  require(data.table)
  require(parallel)
  
  #Make sure patient null vector includes only samples greater than 5
  null.vector<-unique(patient.coverage.vector[patient.coverage.vector>5])
  
  #All cancer patients to sample from
  cancer.patients<-normalized.matrices.object$cancer.patients
  
  #Get number of genes in normalized tables
  all.genes<-length(rownames(normalized.matrices.object$combined.matrices))
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.logFC.Median.Differential.Expression", "null.vector", "normalized.matrices.object", "cancer.patients","all.genes"),envir=environment())
  
  #Build NULL distributions
  Null.dist<-parLapply(cl, null.vector, function(x) {
    dummy.tables<-replicate(100, Function.logFC.Median.Differential.Expression(normalized.matrices.object, sample(cancer.patients,x)),simplify=F)
    
    #Combined all into single data table object [Hugo_Symbol, logFC.Median]
    single.table<-do.call(rbind,dummy.tables)
    
    #Add replication indices
    single.table$replication<-rep(1:100,each=all.genes)
    
    #Add PATIENT.COVERAGE indice
    single.table$PATIENT.COVERAGE<-x
    return(single.table)}
                       )
  
  #Stop parallelization
  stopCluster(cl)
  
  #Combined all into single data.table
  names(Null.dist)<-null.vector
  Null.dist<-do.call(rbind, Null.dist)
  
  #Return
  return(Null.dist)
}

Function.logFC.Median.count.TH.NULL<-function(v.p.Version.4.NULL) {
  #Takes data.table produced by Function.v.p.Version.4.NULL() to produce distributions per sample size(PATIENT.COVERAGE) at each replication for thresholded values
  
  require(data.table)
  require(reshape2)
  
  #Single step with data.table, no parallelization required
  dummy.return<-v.p.Version.4.NULL[,list(t1.5=sum(abs(logFC.Median)>1.5),
                                         t2.0=sum(abs(logFC.Median)>2.0),
                                         t2.5=sum(abs(logFC.Median)>2.5),
                                         t3.0=sum(abs(logFC.Median)>3.0),
                                         t3.5=sum(abs(logFC.Median)>3.5),
                                         t4.0=sum(abs(logFC.Median)>4.0)), 
                                   by=c("PATIENT.COVERAGE", "replication")]
  
  #Clean up and melt
  dummy.return$replication<-NULL #Don't need replication column anymore
  dummy.return<-melt(dummy.return, id="PATIENT.COVERAGE")
  
  #Return
  return(dummy.return)
  
}

Function.logFC.Median.count.TH.NULL.rem.out<-function(Function.logFC.Median.count.TH.NULL.table){
  #Takes output from Function.logFC.Median.count.TH.NULL() and removes first layer of outliers per PATIENT.COVERAGExvariable group
  #Based on 1.5*IQR limits
  
  require(data.table)
  
  dummy.function<-function(x) {
    stat=summary(x)
    top=stat[[5]]
    bottom=stat[[2]]
    iqr1.5=1.5*(top-bottom)
    clean=x[x>=(bottom-iqr1.5) & x<=(top+iqr1.5)]
    return (clean)
  }
  
  dummy.return<-Function.logFC.Median.count.TH.NULL.table[,dummy.function(value), by=c("PATIENT.COVERAGE", "variable")]
  setnames(dummy.return, c("PATIENT.COVERAGE", "variable", "value"))
  
  #Return
  return(dummy.return)
}

Function.v.p.Version.4<-function(normalized.matrices.object, table.1.plus) {
  #Constructs v.p.tables in the form of a data table for all v(p)  genes found in table.1.plus using median.logFC function
  #Adds patient.coverage to each gene
  #Requires:
  #   Function.logFC.Median.Differential.Expression()
  
  require(data.table)
  require(parallel)
  
  #Filter table.1.plus for patients present in expression matrix
  cancer.patients<-normalized.matrices.object$cancer.patients
  table.1.plus<-table.1.plus[PATIENT %in% cancer.patients,]
  
  #Make sure to only work with genes that occur in more than 5 patients
  patient.coverage.table<-table.1.plus[,list(PATIENT.COVERAGE=length(PATIENT)),by="Hugo_Symbol"]
  patient.coverage.table<-patient.coverage.table[PATIENT.COVERAGE>5,]
  table.1.plus<-table.1.plus[Hugo_Symbol %in% unique(as.vector(patient.coverage.table$Hugo_Symbol)),]
  
  #Do diff. expression function
  internal.function<-function(patients) {
    diff.exp<-Function.logFC.Median.Differential.Expression(normalized.matrices.object, patients)
    return(list(HUGOS=as.vector(diff.exp$Hugo_Symbol),
                logFC.Median=as.vector(diff.exp$logFC.Median)))
  }
  
  #Split genes into 8 portions
  all.genes<-unique(as.vector(table.1.plus$Hugo_Symbol))
  all.genes.8<-split(all.genes,factor(1:8))
  
  #Set up parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("Function.logFC.Median.Differential.Expression", "internal.function", "table.1.plus", "all.genes.8", "as.data.table"),envir=environment())
  
  #Execute parallelized differential expression in data.tables
  #diff.expression<-table.1.plus[,internal.function(PATIENT), by="Hugo_Symbol"]
  diff.expression<-parLapply(cl, all.genes.8, function(x)  as.data.table(table.1.plus[Hugo_Symbol %in% x,])[,internal.function(PATIENT), by="Hugo_Symbol"])
  
  #Stop parallelization
  stopCluster(cl)
                     
  #Clean up and add PATIENT.COVERAGE information
  diff.expression<-do.call(rbind, diff.expression)
  diff.expression<-as.data.table(merge(as.data.frame(diff.expression), as.data.frame(patient.coverage.table), by="Hugo_Symbol"))
  
  #Return
  return(diff.expression)
}

BRCA.Table.1.PATIENT.COVERAGE
dim(BRCA.NORM.MATRICES.OBJ$combined.matrices)
Table.1.PLUS

test<-Function.v.p.Version.4(BRCA.NORM.MATRICES.OBJ, Table.1.PLUS)
test
test[Hugo_Symbol=="ATF3",]

Function.v.p.logFC.Median.count.TH<-function(v.p.Version.4) {
  #Takes data.table produced by Function.v.p.Version.4()) to produce distributions per v(p) gene at each thresholded value
  
  require(data.table)
  require(reshape2)
  
  #Single step with data.table, no parallelization required
  dummy.return<-v.p.Version.4[,list(t1.5=sum(abs(logFC.Median)>1.5),
                                         t2.0=sum(abs(logFC.Median)>2.0),
                                         t2.5=sum(abs(logFC.Median)>2.5),
                                         t3.0=sum(abs(logFC.Median)>3.0),
                                         t3.5=sum(abs(logFC.Median)>3.5),
                                         t4.0=sum(abs(logFC.Median)>4.0)), 
                                   by=c("Hugo_Symbol","PATIENT.COVERAGE")]
  
  #Clean up and melt
  dummy.return<-melt(dummy.return, id=c("Hugo_Symbol","PATIENT.COVERAGE"))
  
  #Return
  return(dummy.return)
}

Function.heatmap.2.mut.color<-function(corr.matrix, table.1, target.gene){
    #Takes corr.matrix and produces a vector color for colnames based on given mutated gene in table.1
    #Patients with target gene will be colored in red while the rest will be in green

    GENE.COLOR<-as.character(colnames(corr.matrix) %in% as.vector(table.1[Hugo_Symbol==target.gene]$PATIENT))
    GENE.COLOR<-replace(GENE.COLOR, GENE.COLOR=="TRUE", "red")
    GENE.COLOR<-replace(GENE.COLOR, GENE.COLOR=="FALSE", "green")
    return(GENE.COLOR)
}

#To convert 082714.BRCA.v.p.object.NULL.rds to vectorize v(p) form
null.v.vector<-lapply(names(table.v.null), function(x) {
  target.list<-table.v.null[[x]]
  target.vector<-sapply(target.list, function(y) {
    all.genes<-nrow(y)
    sig.genes<-nrow(y[adj.P.Val<0.05,])
    v.p<-sig.genes/all.genes
    return(v.p)
  })
  return(target.vector)
})

table.v.p<-readRDS("BRCA/090214.BRCA.Table.v.p.rds")
table.1.obj$table.1
null.v.vector<-readRDS("BRCA/100814.BRCA.null.v.vectors.rds")

Function.v.pval<-function(table.v.p, table.1.obj, null.v.vector){
  #Calculates empirical p-values for table.v based on null distribution

  #Filter patients that have at least 5 missense mutations
  table.1.count<-table.1.obj$table.1[Missense!=0,]
  table.1.count<-table.1.count[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
  table.1.count<-table.1.count[N.PATIENTS>5,]

  #Apply count threshold to table.v.p
  table.v.p<-as.data.table(merge(as.data.frame(table.v.p), as.data.frame(table.1.count)))

  #Calculate p.values based on null vector
  table.v.p$P.VAL<-apply(as.matrix(table.v.p), 1, function(x) {
    y<-as.vector(c(x[2],as.numeric(as.character(x[3]))))
    P.VAL<-mean(as.numeric(null.v.vector[[y[2]]])>=as.numeric(y[1]))
    return(P.VAL)
  })

  table.v.p$P.VAL.ADJ<-p.adjust(table.v.p$P.VAL, method="fdr")

}

path<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/NCI.PATHWAYS/091614.path.gene.rds")

Function.path.enrichment<-function(HUGOS, path, universe.hugos){
  #Looks for enriched pathways based on hugos
  #path table in the form of 091614.path.gene.rds
  #universe.hugos are all potential hugos in the universe from where we started filtering our HUGOS set

  require(data.table)

  #Filter path table for universe.hugos
  path<-path[Hugo_Symbol %in% universe.hugos,]

  #Calculate enrichment
  HUGO.P.VAL<-path[,list(P.VAL=
        phyper(q=sum(HUGOS %in% Hugo_Symbol)-1,
            m=length(Hugo_Symbol),
            n=length(universe.hugos)-length(Hugo_Symbol),
            k=length(HUGOS), lower.tail=F)
        , HUGOS.IN.PATH=sum(HUGOS %in% Hugo_Symbol), ALL.IN.PATH=length(Hugo_Symbol)),by=c("Path")]

  #Adjust for multiple hypothesis correction
  HUGO.P.VAL$P.VAL.ADJ<-p.adjust(HUGO.P.VAL$P.VAL, method="fdr")

  #Clean up and return
  HUGO.P.VAL<-HUGO.P.VAL[HUGOS.IN.PATH>=2,]
  HUGO.P.VAL<-HUGO.P.VAL[order(P.VAL.ADJ),]
  return(HUGO.P.VAL)
}

Function.BRCA.SUBTYPE<-function(brca.normalized.obj, version=1){
  #Takes normalized expression object from breast cancer and returns subtype score and assigns subtype to each patient based on maximum score
  
  require(genefu)
  require(data.table)
  
  #Choose pam50 method
  if (version==1){
    data(pam50)
    pam<-copy(pam50)
  } else if(version==2){
    data(pam50.scale)
    pam<-copy(pam50.scale)
  } else if(version==3){
    data(pam50.robust)
    pam<-copy(pam50.robust)
  }
  
  #Obtain target pam gene in expression matrix
  target.genes<-intersect(rownames(pam$centroids), rownames(brca.normalized.obj$combined.matrices))
  
  #Simplify pam and expression matrices
  simplified.pam<-pam$centroids[target.genes,]
  simplified.exp<-brca.normalized.obj$combined.matrices[target.genes, brca.normalized.obj$cancer.patients]
  
  #Obtain predicted type based on correlation
  BRCA.SCORES<-cor(simplified.exp, simplified.pam, method="spearman")
  BRCA.SCORES<-as.data.frame(BRCA.SCORES)
  BRCA.SCORES$PATIENT<-rownames(BRCA.SCORES)
  BRCA.SCORES$TYPE<-colnames(BRCA.SCORES)[max.col(BRCA.SCORES[,1:5])]
  
  #Clean up and return
  BRCA.SCORES<-as.data.table(BRCA.SCORES)
  BRCA.SCORES<-BRCA.SCORES[order(TYPE),]
  return(BRCA.SCORES)
}
