#Function.p.clustering.R
#091114
#Function to validate patients clusters based on their mutation similarity on gene expression

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

Function.p.pre.clustering<-function(table.1, diss.matrix){
    #Cleans table.1 for patients present in expression matrix

    require(data.table)

    #Cleans for those that we actually have expression data for
    main.table<-table.1[PATIENT %in% colnames(diss.matrix),]
    main.table<-droplevels(main.table)

    #Cleans for genes that are found in more than 5 patients
    patient.count<-main.table[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
    patient.count<-patient.count[N.PATIENTS>5,]
    main.table<-main.table[Hugo_Symbol %in% unique(as.vector(patient.count$Hugo_Symbol)),]

    #Return - [Hugo_Symbol, PATIENT]
    return(main.table)
}

Function.p.Silhouette<-function(diss.matrix, main.table){
    #S(j) - Get Modified Silhouette (S) score per Hugo_Symbol
    #This modified S(j) score calculates medians rather than average dissimilarities
    require(data.table)

    internal.function<-function(patient.vector){
        #Get dissimilarity matrix of clustered patients for gene p
        internal.matrix<-diss.matrix[patient.vector, patient.vector]

        #a(i) - Get median dissimilarity for each member of the cluster
        diag(internal.matrix)<-NA
        internal.a.i<-apply(internal.matrix, 1, function(x) median(x, na.rm=T))
        internal.a.i<-as.data.table(internal.a.i, keep.rownames=T)
        setnames(internal.a.i, c("PATIENT", "A.I"))
        
        #Return
        return(internal.a.i)
    }

    internal.function.2<-function(other.patients, self.patient) {
        internal.matrix<-diss.matrix[self.patient, other.patients]
        internal.median<-median(internal.matrix)
        return(internal.median)
    }

    #Calculate a(i) for each patient in each cluster [Hugo_Symbol, PATIENT, A.I]
    a.i.table.p<-main.table[,internal.function(PATIENT), by="Hugo_Symbol"]
    print ("done building a.i")

    #Calculate b(i) for each patient in each cluster with parallelization
    p.vector<-unique(as.vector(main.table$Hugo_Symbol))

    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    clusterExport(cl, varlist=c("p.vector", "main.table","as.data.table","internal.function.2", "diss.matrix"),envir=environment())

    b.i.table.p<-parLapply(cl,p.vector, function(x) {
        print (x)
        self.genes<-x
        self.patients<-as.vector(as.data.table(main.table)[Hugo_Symbol==self.genes,]$PATIENT)
        
        non.self.main.table<-as.data.table(main.table)[Hugo_Symbol!=self.genes,]
        
        #Get lowest median similarity b(i) to any other cluster b.i.c that is not self
        min.vector<-sapply(self.patients, function(y) {

            #Filter for clusters in wich self (y) is not a member
            y.clusters<-unique(as.vector(non.self.main.table[PATIENT==y,]$Hugo_Symbol))
            if (length(y.clusters)!=0) {
                non.y.table<-non.self.main.table[Hugo_Symbol %in% y.clusters,]    
            } else {
                non.y.table<-non.self.main.table
            }
            min.diss<-min(as.vector(non.y.table[,list(b.i.c=internal.function.2(PATIENT, y)), by="Hugo_Symbol"]$b.i.c))
            return(min.diss)
        })

        cluster.b.i<-data.frame(PATIENT=self.patients, B.I=min.vector)
        cluster.b.i$Hugo_Symbol<-self.genes
        cluster.b.i<-as.data.table(cluster.b.i)

        #Return
        return(cluster.b.i)
    })

    print ("done building b.i")
    names(b.i.table.p)<-p.vector
    stopCluster(cl)
    b.i.table.p<-do.call(rbind, b.i.table.p)
 
    #Combine a(i) and b(i) tables to calculate s(i) per patient in cluster
    s.i.table<-as.data.table(merge(as.data.frame(a.i.table.p), as.data.frame(b.i.table.p), by=c("Hugo_Symbol", "PATIENT")))
    s.i.table$S.I<-(s.i.table$B.I - s.i.table$A.I) /pmax(s.i.table$B.I, s.i.table$A.I)

    #Calculate Silhouette coeffient S(c) per cluster to measure its validity
    #This measures how appropriate the data has been clustered
    s.i.per.cluster<-s.i.table[,list(S.I.CLUSTER=median(S.I)), by="Hugo_Symbol"]

    #Return Cluster Silhouette S(c), a.i.table and b.i.table
    return (list(S.I.C.TABLE=s.i.per.cluster, A.I.TABLE=a.i.table.p, B.I.TABLE=b.i.table.p))

}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
exp.matrix<-readRDS(args[2])
output.file<-args[3]
print("opened files")

diss.matrix<-Function.p.diss.matrix(exp.matrix, 1)
print ("done building diss matrix")

main.table<-Function.p.pre.clustering(table.1, diss.matrix)
print ("done pre-clustering")

s.i.results<-Function.p.Silhouette(diss.matrix, main.table)
print ("done building s.i")

saveRDS(s.i.results, output.file)
print ("done saving")