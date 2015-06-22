#Function.j.clustering.R
#090914
#Function to validate patient cluster based on their metabolic mutation similarity on gene expression

Function.j.pre.clustering<-function(table.2, table.1){
    #Clusters patients by metabolite

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

    #Return clustered patients by metabolite
    return(main.table)
}

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

Function.j.Silhouette<-function(diss.matrix, main.table){
    #S(j) - Get Modified Silhouette (S) score per kegg_id (j metabolite)
    #This modified S(j) score calculates medians rather than average dissimilarities
    require(data.table)

    internal.function<-function(patient.vector){
        #Get dissimilarity matrix of clustered patients for metabolite j
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

    #Filter main.table for patient information that we actually have
    main.table<-main.table[PATIENT %in% colnames(diss.matrix),]

    #Calculate a(i) for each patient in each cluster
    a.i.table<-main.table[,internal.function(PATIENT), by="KEGG_ID"]
    print ("done building a.i")

    #Calculate b(i) for each patient in each cluster with parallelization
    j.vector<-unique(as.vector(main.table$KEGG_ID))

    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)
    clusterExport(cl, varlist=c("j.vector", "main.table","as.data.table","internal.function.2", "diss.matrix"),envir=environment())

    b.i.table<-parLapply(cl, j.vector, function(x) {
        print (x)
        self.kegg<-x
        self.patients<-as.vector(as.data.table(main.table)[KEGG_ID==self.kegg,]$PATIENT)
        
        non.self.main.table<-as.data.table(main.table)[KEGG_ID!=self.kegg,]
        
        #Get lowest median similarity b(i) to any other cluster b.i.c that is not self
        min.vector<-sapply(self.patients, function(y) {

            #Filter for clusters in wich self (y) is not a member
            y.clusters<-unique(as.vector(non.self.main.table[PATIENT==y,]$KEGG_ID))
            if (length(y.clusters)!=0) {
                non.y.table<-non.self.main.table[KEGG_ID %in% y.clusters,]    
            } else {
                non.y.table<-non.self.main.table
            }
            min.diss<-min(as.vector(non.y.table[,list(b.i.c=internal.function.2(PATIENT, y)), by="KEGG_ID"]$b.i.c))
            return(min.diss)
        })
        
        cluster.b.i<-data.frame(PATIENT=self.patients, B.I=min.vector)
        cluster.b.i$KEGG_ID<-self.kegg
        cluster.b.i<-as.data.table(cluster.b.i)

        #Return
        return(cluster.b.i)
    })

    print ("done building b.i")
    names(b.i.table)<-j.vector
    stopCluster(cl)
    b.i.table<-do.call(rbind, b.i.table)
 
    #Combine a(i) and b(i) tables to calculate s(i) per patient in cluster
    s.i.table<-as.data.table(merge(as.data.frame(a.i.table), as.data.frame(b.i.table), by=c("KEGG_ID", "PATIENT")))
    s.i.table$S.I<-(s.i.table$B.I - s.i.table$A.I) /pmax(s.i.table$B.I, s.i.table$A.I)

    #Calculate Silhouette coeffient S(c) per cluster to measure its validity
    #This measures how appropriate the data has been clustered
    s.i.per.cluster<-s.i.table[,list(S.I.CLUSTER=median(S.I)), by="KEGG_ID"]

    #Return Cluster Silhouette S(c) along with a(i) and b(i) tables
    return (list(S.I.C.TABLE=s.i.per.cluster, A.I.TABLE=a.i.table, B.I.TABLE=b.i.table))
}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
table.2<-readRDS(args[2])
exp.matrix<-readRDS(args[3])
output.file<-args[4]
print("opened files")

main.table<-Function.j.pre.clustering(table.2, table.1)
print ("done pre-clustering")

diss.matrix<-Function.j.diss.matrix(exp.matrix, 1)
print ("done building diss matrix")

s.i.results<-Function.j.Silhouette(diss.matrix, main.table)
print ("done building s.i")

saveRDS(s.i.results, output.file)
print ("done saving")