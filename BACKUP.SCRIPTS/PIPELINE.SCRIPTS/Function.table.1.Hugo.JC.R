#Function.table.1.Hugo.JC.R - HALTED
#092714
#Converts significant genes found in tables based on Function.table.1.bmr() to a table for cytoscape


Function.JC<-function(table.1.pval, coverage.threshold){

    require(data.table)
    require(reshape2)
    require(Hmisc)
    require(plyr)

    #Filter for significant values
    table.1.pval<-table.1.pval[P.VAL.ADJ<0.05,]

    #Filter  table 1 for patients above coverage threshold
    table.1.coverage<-table.1.pval[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
    table.1.coverage<-table.1.coverage[N.PATIENTS>=coverage.threshold,]
    table.1.pval<-table.1.pval[Hugo_Symbol %in% unique(as.vector(table.1.coverage$Hugo_Symbol)),]
    print ("done filtering for gene coverage")

    #Convert to wide format
    table.1.pval<-unique(table.1.pval[,c("PATIENT", "Hugo_Symbol","Missense"), with=F])
    table.1.cast<-acast(table.1.pval, PATIENT~Hugo_Symbol, fill=0,fun.aggregate=max)
    print ("done converting to wide format")

    #Obtain gene pairwise spearman correlation table of rhos and pvalues
    #This may take a while
    #NOTE - Keep in mind that this is a correlation done for all present known mutations while removing all not-present mutations across patients as missing values

    #Set up parallelization
    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    print ("done setting up parallelization")

    print("calculating correlations, this may take a while")
    gene.combos<-t(combn(colnames(table.1.cast), 2))
    clusterExport(cl, varlist=c("table.1.cast","cor.test", "gene.combos"),envir=environment())
    
    table.1.corr<-parSapply(cl, 1:nrow(gene.combos), function(x) {
        
        #Obtain pairwise table
        target.table<-table.1.cast[,c(gene.combos[x,1] , gene.combos[x,2])]
        
        #Remove missing zeroes corresponding to patients that have no mutations for both genes
        target.table<-target.table[apply(target.table, 1, sum)!=0,]

        #Correlation on true known values
        spearman<-cor.test(target.table[,1], target.table[,2], method="spearman")

        #Return
        return(c(spearman$p.value, spearman$estimate))

    }, USE.NAMES=T, simplify=T)

    #Stop parallelization
    print ("done building correlations")
    stopCluster(cl)

    #Clean up and assign correct Hugo identifiers
    table.1.corr<-t(table.1.corr)
    table.1.corr<-as.data.table(cbind(gene.combos, table.1.corr))
    setnames(table.1.corr, c("Hugo.1", "Hugo.2", "P.VAL", "RHO"))

    table.1.corr$P.VAL<-as.numeric(as.character(table.1.corr$P.VAL))
    table.1.corr$RHO<-as.numeric(as.character(table.1.corr$RHO))

    #Correct for multiple hypothesis    
    table.1.corr$P.VAL.ADJ<-p.adjust(table.1.corr$P.VAL, method="fdr")



    #ALTERNATIVE
    print ("Obtaining correlation matrix, this step may take a while")
    table.1.corr<-rcorr(table.1.cast, type="spearman")

    #Convert rho and p.values tables to short format 
    rho.dist<-as.dist(table.1.corr$r)
    rho.melt<-as.data.table(data.frame(t(combn(rownames(table.1.corr$r), 2)), as.numeric(rho.dist)))
    setnames(rho.melt, c("Hugo.1", "Hugo.2", "RHO"))

    p.dist<-as.dist(table.1.corr$P)
    p.melt<-as.data.table(data.frame(t(combn(rownames(table.1.corr$P), 2)), as.numeric(p.dist)))
    setnames(p.melt, c("Hugo.1", "Hugo.2", "P.VAL"))

    #Correct p-value table for multiple hypothesis testing
    print ("Correcting correlations for multiple hypothesis testing")
    p.melt$P.VAL.ADJ<-p.adjust(p.melt$P.VAL, method="fdr")
    p.melt<-p.melt[P.VAL.ADJ<0.05,]

    #Match to RHO values for those that pass multiple hypothesis correction
    print ("Obtaining matching significant RHOs")
    table.1.melt<-join(p.melt, rho.melt, by=c("Hugo.1", "Hugo.2"))    

    #Append thresholded patient coverage
    table.1.melt<-as.data.table(merge(as.data.frame(table.1.melt), as.data.frame(table.1.coverage), by.x="Hugo.1", by.y="Hugo_Symbol"))
    setnames(table.1.melt, c(colnames(table.1.melt)[1:5], "Hugo.1.N.PATIENTS"))

    table.1.melt<-as.data.table(merge(as.data.frame(table.1.melt), as.data.frame(table.1.coverage), by.x="Hugo.2", by.y="Hugo_Symbol"))
    setnames(table.1.melt, c(colnames(table.1.melt)[1:6], "Hugo.2.N.PATIENTS"))


}