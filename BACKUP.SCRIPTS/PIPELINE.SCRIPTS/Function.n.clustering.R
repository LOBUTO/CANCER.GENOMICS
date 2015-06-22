#Function.n.clustering.R
#091214
#Function to validate patient clusters from kmeans results

Function.n.diss.matrix<-function(exp.matrix, dissimilarity=1){
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

Function.n.pre.clustering<-function(cluster.table, diss.matrix){
    #Cleans up cluster.table obtained from kmeans result
    #cluster.table pre-processed as following
    #   cluster.table<-as.data.table(k.means.result$cluster, keep.rownames=T)
    #   setnames(cluster.table, c("PATIENT", "CLUSTER"))

    require(data.table)

    #Cleans for those that we actually have expression data for
    main.table<-cluster.table[PATIENT %in% colnames(diss.matrix),]
    main.table<-droplevels(main.table)

    #Filters for clusters that contain more than 5 patients
    patient.count<-main.table[,list(N.PATIENTS=length(PATIENT)), by="CLUSTER"]
    patient.count<-patient.count[N.PATIENTS>5,]
    main.table<-main.table[CLUSTER %in% unique(as.vector(patient.count$CLUSTER)),]

    #Return - [CLUSTER, PATIENT]
    return(main.table)
}

Function.n.Silhouette<-function(diss.matrix, main.table){
    #S(j) - Get Modified Silhouette (S) score per CLUSTER
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
    a.i.table<-main.table[,internal.function(PATIENT), by="CLUSTER"]
    print ("done building a.i")

    #Calculate b(i) for each patient in each cluster with parallelization
    n.vector<-unique(as.vector(main.table$CLUSTER))

    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    clusterExport(cl, varlist=c("n.vector", "main.table","as.data.table","internal.function.2", "diss.matrix"),envir=environment())

    b.i.table<-parLapply(cl,n.vector, function(x) {
        print (x)
        self.cluster<-x
        self.patients<-as.vector(as.data.table(main.table)[CLUSTER==self.cluster,]$PATIENT)
        
        non.self.main.table<-as.data.table(main.table)[CLUSTER!=self.cluster,]
        
        #Get lowest median similarity b(i) to any other cluster b.i.c that is not self
        min.vector<-sapply(self.patients, function(y) {

            #Filter for clusters in wich self (y) is not a member
            y.clusters<-unique(as.vector(non.self.main.table[PATIENT==y,]$CLUSTER))
            if (length(y.clusters)!=0) {
                non.y.table<-non.self.main.table[CLUSTER %in% y.clusters,]    
            } else {
                non.y.table<-non.self.main.table
            }
            min.diss<-min(as.vector(non.y.table[,list(b.i.c=internal.function.2(PATIENT, y)), by="CLUSTER"]$b.i.c))
            return(min.diss)
        })

        cluster.b.i<-data.frame(PATIENT=self.patients, B.I=min.vector)
        cluster.b.i$CLUSTER<-self.cluster
        cluster.b.i<-as.data.table(cluster.b.i)

        #Return
        return(cluster.b.i)
    })

    print ("done building b.i")
    names(b.i.table)<-n.vector
    stopCluster(cl)
    b.i.table<-do.call(rbind, b.i.table)
 
    #Combine a(i) and b(i) tables to calculate s(i) per patient in cluster
    s.i.table<-as.data.table(merge(as.data.frame(a.i.table), as.data.frame(b.i.table), by=c("CLUSTER", "PATIENT")))
    s.i.table$S.I<-(s.i.table$B.I - s.i.table$A.I) /pmax(s.i.table$B.I, s.i.table$A.I)

    #Calculate Silhouette coeffient S(c) per cluster to measure its validity
    #This measures how appropriate the data has been clustered
    s.i.per.cluster<-s.i.table[,list(S.I.CLUSTER=median(S.I)), by="CLUSTER"]

    #Return Cluster Silhouette S(c)
    return (s.i.per.cluster)
}

