#Function.Cummulative.Correlation - HALT, DOESN'T SEEM TO WORK
#103014
#Obtains pairwise correlations across genes bases on cummulative significant mutations related to cancer phenotype based on expresssion distance measure

Function.dist.matrix<-function(exp.matrix){
    #Pre-process exp matrix to produced median distance matrix to normal
    require(data.table)

    #First get pairwise distance matrix
    CORR.MATRIX<-cor(exp.matrix$combined.matrices[,c(exp.matrix$normal.patients,
    exp.matrix$cancer.patients)], method="spearman")
    DISS.MATRIX<-1-CORR.MATRIX

    #Get median distance of cancer patients to normal patients
    CANCER.DIST<-sapply(exp.matrix$cancer.patients, function(x) median(DISS.MATRIX[x, exp.matrix$normal.patients]))
    CANCER.DIST<-as.data.table(CANCER.DIST, keep.rownames=T)
    setnames(CANCER.DIST, c("PATIENT", "DIST.TO.NORMAL"))

    #Return
    return(CANCER.DIST)
}

Function.main.table<-function(table.1.pval, CANCER.DIST){
    #Pre-filters main table

    require(data.table)

    #Select significant mutations above background
    table.1.pval<-table.1.pval[P.VAL.ADJ<0.05,]

    #Filter  table 1 for patients above coverage threshold
    #table.1.coverage<-table.1.pval[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
    #table.1.coverage<-table.1.coverage[N.PATIENTS>=coverage.threshold,]
    #table.1.pval<-table.1.pval[Hugo_Symbol %in% unique(as.vector(table.1.coverage$Hugo_Symbol)),]

    #Assign phenotype distance to genes
    main.table<-as.data.table(merge(as.data.frame(table.1.pval), as.data.frame(CANCER.DIST)))

    #Return
    return(main.table)
}

Function.Cum.Corr<-function(main.table){
    #Obtain correlations for cumulative vectors based on cancer distance

    require(data.table)
    require(reshape2)
    require(Hmisc)

    #Convert main table to wide format
    main.table.melt<-acast(main.table, PATIENT~Hugo_Symbol,value.var="Missense" , fill=0,fun.aggregate=max)
    main.table.melt<-as.data.table(main.table.melt, keep.rownames=T)
    setnames(main.table.melt, c("PATIENT", 
        colnames(main.table.melt)[2: length(colnames(main.table.melt))]))

    #Get cummulative mutations vectors per gene based on phenotype distance
    main.table.melt<-as.data.table(merge(as.data.frame(main.table.melt),
        as.data.frame(unique(main.table[,c("PATIENT","DIST.TO.NORMAL"), with=F])),
        by="PATIENT"))
    main.table.melt<-main.table.melt[order(DIST.TO.NORMAL),]

    #Gene vectors
    genes<-colnames(main.table.melt)[2:(length(colnames(main.table.melt))-1)]
    cum.vectors<-apply(main.table.melt[,genes,with=F], 2, cumsum)

    #Obtain gene correlations
    main.corr<-rcorr(cum.vectors, type="spearman")

    #Clean up RHO and p-values
    rho.dist<-as.dist(main.corr$r)
    rho.corr<-as.data.table(data.frame(t(combn(rownames(main.corr$r), 2)), as.numeric(rho.dist)))
    setnames(rho.corr, c("Hugo.1", "Hugo.2", "RHO"))

    p.dist<-as.dist(main.corr$P)
    p.corr<-as.data.table(data.frame(t(combn(rownames(main.corr$P), 2)), as.numeric(p.dist)))
    setnames(p.corr, c("Hugo.1", "Hugo.2", "P.VAL"))

    #Correct for multiple hypothesis testing
    p.corr$P.VAL.ADJ<-p.adjust(p.corr$P.VAL, method="fdr")
    p.corr<-p.corr[P.VAL.ADJ<0.05,]

    #Merge RHO and corrected p-values
    main.corr<-join(p.corr, rho.corr, by=c("Hugo.1", "Hugo.2"))

}

