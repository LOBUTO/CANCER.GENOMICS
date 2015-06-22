#Function.Clust.Gene.Enrich.R
#100414
#Function aims to find driver genes at the patient level in a multi-step process:
#   1. Performs differential expression and filters based on  log(FC) threshold
#   2. Perform hierarchical clustering based on euclidean distance of patient vectors for genes filtered on step 1
#   3. Uses table.1.pval (BMR method) on each cluster to find driver genes by wilcoxon based on background mutation rate - Use update table.1.pval with aa length information

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
  DIFF.EXP<-all.G.fit
  return(DIFF.EXP)
}

Function.Hclust.Groups<-function(DIFF.EXP, EXP.OBJ, diff.threshold){
    #Finds patients groups based on hclustering using thresholded differential expression results

    require(data.table)

    #TEST
    #DIFF.EXP<-copy(BRCA.DIFF.EXP)
    #EXP.OBJ<-copy(BRCA.EXP.OBJ)
    #diff.threshold<-3

    #Filter differentially expressed genes by threshold
    DIFF.EXP<-DIFF.EXP[logFC>=diff.threshold,]
    TARGET.GENES<-as.vector(DIFF.EXP$ID)

    #Use hierarchical clustering to build cluster tree based on target gene expression levels
    EXP.DIST<-dist(t(BRCA.EXP.OBJ$combined.matrices[TARGET.GENES, BRCA.EXP.OBJ$cancer.patients]), method="euclidean")
    EXP.HCLUST<-hclust(EXP.DIST, method="ward.D")
    
    #Obtain cluster groups for different number of clusters
    EXP.HCLUST.GROUPS<-lapply(5:10, function(x) {
        
        #Obtain cluster groups per k
        HCLUST.GROUP<-as.data.table(cutree(EXP.HCLUST, k=x), keep.rownames=T)
        setnames(HCLUST.GROUP, c("PATIENT", "CLUSTER.NUMBER"))
        HCLUST.GROUP<-HCLUST.GROUP[order(CLUSTER.NUMBER),]
        HCLUST.GROUP$K.CHOICE<-x

        #Return
        return(HCLUST.GROUP)
    })

    #Combine cluster groups per k into single table 
    HCLUST.TABLE<-do.call(rbind, EXP.HCLUST.GROUPS)

    #Return
    return(HCLUST.TABLE)
}

Function.Hclust.Enrichment<-function(HCLUST.TABLE, table.1.pval, count.filter, cnv.table){
    #TEST 1 - Caculates wilcoxon signed rank test for each gene with cluster based on gene mutation rates vs background mutation rates within clusters
    #TEST 2- Calculates hypergeometric per gene within cluster to ratio outside cluster

    require(data.table)
    require(reshape2)

    #Obtain total background length - Theoretical protein length of all possible sequenced genes
    #USE ALL PROTEINS FOR THIS!!, DON'T FILTER FOR FDR YET
    background.length<-unique(table.1.pval[,c("Hugo_Symbol", "Length"), with=F])
    background.length<-sum(as.vector(background.length$Length))   

    #Get total mutations per patients [PATIENT, BMR]
    #USE ALL PROTEINS FOR THIS!!, DON'T FILTER FOR FDR YET
    BM.table<-table.1.pval[,list(BM=sum(Missense)), by="PATIENT"]
    BM.table$BMR<-BM.table$BM/background.length
    BM.table$BM<-NULL

    #Filter table 1 for fdr<0.05
    table.1.pval<-table.1.pval[P.VAL.ADJ<0.05,]

    #Obtain mutation rate per gene in each patient
    GM.table<-table.1.pval[,list(GMR=Missense/Length), by=c("PATIENT", "Hugo_Symbol")]

    #Cast GM.table to convert to main.table to have x and 0s for BMR
    main.table<-acast(GM.table, PATIENT~Hugo_Symbol, value.var="GMR", fill=0)
    main.table<-as.data.table(melt(main.table))
    setnames(main.table, c("PATIENT", "Hugo_Symbol", "GMR"))

    #Merge with background mutation rate
    main.table<-as.data.table(merge(as.data.frame(main.table), as.data.frame(BM.table), by="PATIENT"))

    #TEST 1 - Perform wilcoxon
    internal.function<-function(samples){
        #Obtains wilcoxon per gene in each cluster
        wilcox.function<-function(gmr,bmr){
            w.test<-wilcox.test(gmr, bmr, paired=T, alternative="greater")
            return(list(W=w.test$statistic, P.VAL=w.test$p.value))
        }

        #Obtains patients in cluster
        target.table<-main.table[PATIENT %in% samples,]
        
        #Wilcox signed ranked test
        target.table<-target.table[,wilcox.function(GMR,BMR), by="Hugo_Symbol"]

        #Return
        return(list(Hugo_Symbol=as.vector(target.table$Hugo_Symbol),
            P.VAL=as.vector(target.table$P.VAL), W=as.vector(target.table$W)))
    }

    HCLUST.WILCOXON<-HCLUST.TABLE[, internal.function(PATIENT) , by=c("CLUSTER.NUMBER", "K.CHOICE")]
    HCLUST.WILCOXON$P.VAL.ADJ<-p.adjust(HCLUST.WILCOXON$P.VAL, method="fdr")

    #TEST 2 - Perform Hypergeometric using mutations and CNV
    #Since we are testing presence, this will only be binary, not presence or absence

    #TEST
    count.filter<-3
    cnv.table<-copy(brca.cnv.table)

    #Integrate mutation and CNV tables
    variation.table<-table.1.pval[,c("PATIENT","Hugo_Symbol", "P.VAL.ADJ"),with=F]
    variation.table$TYPE<-"MUTATION"

    cnv.table<-cnv.table[,c("Hugo_Symbol", "PATIENT", "P.VAL.ADJ"), with=F]
    cnv.table$TYPE<-"CNV"    

    variation.table<-rbind(variation.table, cnv.table)

    #Obtain total count of all significant variations across all patients in clusters
    #Maybe need to make one for mutation and one for CNV!!
    OVERALL.GENE.COUNT<-as.data.table(table(variation.table[PATIENT %in% as.vector(HCLUST.TABLE$PATIENT),]$Hugo_Symbol))
    setnames(OVERALL.GENE.COUNT, c("Hugo_Symbol", "Freq"))
    TOTAL.GENE.COUNT<-sum(as.vector(OVERALL.GENE.COUNT$Freq))

    #Build thresholded table 1
    #table.1.thresholded<-table.1.pval[,c("PATIENT","Hugo_Symbol"),with=F]

    #Merge with cluster table to obtain per cluster mutated genes and frequencies with filtering
    HCLUST.TABLE.THRESHOLDED<-as.data.table(merge(as.data.frame(HCLUST.TABLE), as.data.frame(variation.table), by="PATIENT"))
    HCLUST.TABLE.THRESHOLDED$COUNT<-1
    CLUSTER.COUNT<-HCLUST.TABLE.THRESHOLDED[,list(CLUSTER.COUNT=sum(COUNT)), by=c("K.CHOICE","CLUSTER.NUMBER")]

    #Apply filtering before merging since we have total counts already
    HCLUST.TABLE.THRESHOLDED<-HCLUST.TABLE.THRESHOLDED[P.VAL.ADJ<0.05, ]
    HCLUST.TABLE.THRESHOLDED<-as.data.table(merge(as.data.frame(HCLUST.TABLE.THRESHOLDED), as.data.frame(CLUSTER.COUNT), by=c("CLUSTER.NUMBER", "K.CHOICE")))

    #Merge to obtain total gene frequency across all clusters and threshold by count.filter
    OVERALL.GENE.COUNT<-OVERALL.GENE.COUNT[Freq>=count.filter,]
    HCLUST.TABLE.THRESHOLDED<-as.data.table(merge(as.data.frame(HCLUST.TABLE.THRESHOLDED), as.data.frame(OVERALL.GENE.COUNT), by="Hugo_Symbol"))
    #HCLUST.TABLE.THRESHOLDED<-HCLUST.TABLE.THRESHOLDED[Freq>=count.filter,]

    #Calculate hypergeometric
    HCLUST.HYPER<-HCLUST.TABLE.THRESHOLDED[,list(P.VAL=phyper(q=sum(COUNT)-1,
        m=Freq, n=TOTAL.GENE.COUNT-Freq, k=CLUSTER.COUNT, lower.tail=F)),
    by=c("CLUSTER.NUMBER", "K.CHOICE", "Hugo_Symbol")]

    HCLUST.HYPER$P.VAL.ADJ<-p.adjust(HCLUST.HYPER$P.VAL, method="fdr")

    unique(HCLUST.HYPER[P.VAL.ADJ<0.05,][,1:3,with=F])
    unique(HCLUST.HYPER[P.VAL.ADJ<0.05,][,2:3,with=F])[,list(GENE.COUNT=length(Hugo_Symbol)), by="K.CHOICE"]


}