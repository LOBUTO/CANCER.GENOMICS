#Function.j.clustering.null.ai.R
#091314
#Creates null distribution for Silhouette a(i) metric 
#To be applied to Function.j.clustering.R() a(i) reults
#This will help us support the hypothesis that, given a cluster, how significant it is that all its members were found together not by chance, that is, that they had something in common.

Function.j.diss.matrix<-function(exp.matrix, dissimilarity=1){
    #Builds dissimilarity matrix for all patients in expression matrix

    require(Hmisc)

    #obtain target matrix
    target.matrix<-exp.matrix$combined.matrices[,exp.matrix$cancer.patients]

    #Get correlation
    corr.matrix<-cor(target.matrix, method="spearman")

    #Get dissimilarity matrix based on choice
    if (dissimilarity==1){
        diss.matrix<-1-corr.matrix
    } else if (dissimilarity==2){
        diss.matrix<-sqrt(1-corr.matrix^2)
    }
    #Return dissimilarity matrix
    return(diss.matrix)
}

Function.j.null.pre.clustering<-function(table.1, table.2){
    #Obtains sampling size for a(i) null distribution based on KEGG coverage

    require(data.table)
    require(reshape2)

    #Setup tables
    table.2<-table.2[,2:3, with=F]
    setnames(table.2, c("KEGG_ID", "Hugo_Symbol"))

    #Combine to get a per metabolite patient information [KEGG_ID, PATIENT]
    main.table<-as.data.table(merge(as.data.frame(table.2), as.data.frame(table.1), by="Hugo_Symbol"))
    main.table<-unique(main.table[,2:3, with=F])

    #Drop unused levels
    main.table<-droplevels(main.table)
    patient.count<-main.table[,list(N.PATIENTS=length(PATIENT)), by="KEGG_ID"]
    patient.vector<-unique(as.vector(patient.count$N.PATIENTS))

    return(patient.vector)
}

Function.j.ai.null<-function(patient.vector, diss.matrix) {
    #Constructs null distribution for each sample size in patient.vector
    
    require(data.table)

    internal.null.function<-function(sample.size){
        #Gets dissimilarity median per cluster of patients given sample.size
        patient.vector<-sample(colnames(diss.matrix), sample.size)
        internal.matrix<-diss.matrix[patient.vector, patient.vector]

        #a(i) - Get median dissimilarity for each member of the cluster
        diag(internal.matrix)<-NA
        internal.a.i<-apply(internal.matrix, 1, function(x) median(x, na.rm=T))
        internal.a.i<-as.data.table(internal.a.i, keep.rownames=T)
        setnames(internal.a.i, c("PATIENT", "A.I"))
        internal.null.a.i<-median(as.vector(internal.a.i$A.I))

        #Return
        return(internal.null.a.i)
    }

    #Set up parallelization
    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    clusterExport(cl, varlist=c("patient.vector", "internal.null.function", "diss.matrix","as.data.table","setnames"),envir=environment())

    #Calculate a(i) for each patient in each cluster
    a.i.null.table<-parLapply(cl,patient.vector, function(x) {
        print(x)
        a.i.size.vector<-replicate(100,internal.null.function(x))
        a.i.size.table<-data.frame(SIZE=x, RANDOM.CLUSTER.A.I=a.i.size.vector)
        a.i.size.table<-as.data.table(a.i.size.table)
        return(a.i.size.table)
        } )

    #Stop parallelization
    stopCluster(cl)

    a.i.null.table<-do.call(rbind, a.i.null.table)
    print ("done building a.i")

    return(a.i.null.table)
}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
table.2<-readRDS(args[2])
exp.matrix<-readRDS(args[3])
output.file<-args[4]
print("opened files")

diss.matrix<-Function.j.diss.matrix(exp.matrix, 1)
print ("done building diss matrix")

patient.vector<-Function.j.null.pre.clustering(table.1, table.2)
print ("done pre-clustering")

s.i.results<-Function.j.ai.null(patient.vector,diss.matrix)
print ("done building s.i")

saveRDS(s.i.results, output.file)
print ("done saving")