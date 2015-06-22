#Function.p.clustering.null.ai.R
#091414
#Creates null distribution for Silhouette a(i) metric (x1000 per sample size)
#To be applied to Function.p.clustering.R() a(i) reults
#This will help us support the hypothesis that, given a cluster, how significant it is that all its members were found together not by chance, that is, that they had something in common.

Function.p.diss.matrix<-function(exp.matrix, dissimilarity=1){
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

Function.p.null.pre.clustering<-function(table.1){
    #Obtains sampling size for a(i) null distribution based on GENE coverage

    require(data.table)
    require(reshape2)

    #Remove genes that cover less than 5 patients
    table.1.coverage<-table.1[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
    table.1.coverage<-table.1.coverage[N.PATIENTS>4,]

    #Get patient vector
    patient.vector<-unique(as.vector(table.1.coverage$N.PATIENTS))

    return(patient.vector)
}

Function.p.ai.null<-function(patient.vector, diss.matrix) {
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
        
        a.i.size.vector<-replicate(100,internal.null.function(x))
        a.i.size.table<-data.frame(SIZE=x, RANDOM.CLUSTER.A.I=a.i.size.vector)
        a.i.size.table<-as.data.table(a.i.size.table)
        
        return(a.i.size.table)
        } )

    #Stop parallelization
    stopCluster(cl)

    a.i.null.table<-do.call(rbind, a.i.null.table)
    print ("done building a.i")

    #Return
    return(a.i.null.table)
}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
exp.matrix<-readRDS(args[2])
output.file<-args[3]
print("opened files")

diss.matrix<-Function.p.diss.matrix(exp.matrix, 1)
print ("done building diss matrix")

patient.vector<-Function.p.null.pre.clustering(table.1)
print ("done pre-clustering")

s.null.i.results<-Function.p.ai.null(patient.vector, diss.matrix)
print ("done building a.i null")

saveRDS(s.null.i.results, output.file)
print ("done saving")
