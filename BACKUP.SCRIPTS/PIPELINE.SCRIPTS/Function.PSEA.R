#Function.PSEA.R
#091614
#Obtains enrichment score based on ability of gene to cluster cancer patients

#table.1, exp.matrix, 

Function.PSEA.clust<-function(target.patients, dissimilarity=1){
    #Builds dissimilarity matrix for all patients in expression matrix

    print(count.var/total.count)
    
    #obtain target matrix
    target.matrix<-exp.matrix$combined.matrices[,c(target.patients, exp.matrix$normal.patients)]

    #Get correlation
    corr.matrix<-cor(target.matrix, method="spearman")

    #Get dissimilarity matrix based on choice
    if (dissimilarity==1){
        diss.matrix<-1-corr.matrix
    } else if (dissimilarity==2){
        diss.matrix<-sqrt(1-corr.matrix^2)
    }

    #Get hclust tree and clustered vector
    target.hclust<-hclust(as.dist(diss.matrix), method="complete")
    patient.clust.vector<-target.hclust$labels[target.hclust$order]

    #Classify for target patient vector
    patient.clust.vector<-as.character(patient.clust.vector %in% target.patients)
    patient.clust.vector<-replace(patient.clust.vector, patient.clust.vector=="TRUE", 1)
    patient.clust.vector<-replace(patient.clust.vector, patient.clust.vector=="FALSE", -1)

    #PSEA for pos and negative running vectors
    POS.PSEA<-max(cumsum(patient.clust.vector))
    NEG.PSEA<-max(cumsum(rev(patient.clust.vector)))

    #Maybe weight patient by number of mutations????

    #Update count
    count.var<<-count.var+1

    #Return maximum weighted PSEA - Weighted based on the maximum peaked it could have reached
    #   That is, if all patients lined up perfectly together
    PSEA<-max(POS.PSEA, NEG.PSEA)/length(target.patients)

    return(list(PSEA.SCORE=PSEA))
}

Function.p.pre.clustering<-function(table.1, exp.matrix){
    #Cleans table.1 for patients present in expression matrix

    require(data.table)

    #Cleans for those that we actually have expression data for
    main.table<-table.1[PATIENT %in% exp.matrix$cancer.patients,]
    main.table<-droplevels(main.table)

    #Cleans for genes that are found in more than 5 patients
    patient.count<-main.table[,list(N.PATIENTS=length(PATIENT)), by="Hugo_Symbol"]
    patient.count<<-patient.count[N.PATIENTS>5,]
    main.table<-main.table[Hugo_Symbol %in% unique(as.vector(patient.count$Hugo_Symbol)),]

    #Return - [Hugo_Symbol, PATIENT]
    return(main.table)
}

Function.PSEA<-function(main.table) {
    #Calculates PSEA score per gene    
    
    #Total count
    total.count<<-length(unique(as.vector(main.table$Hugo_Symbol)))

    #Set up parallelization
    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)

    hugo.portions<-split(unique(as.vector(main.table$Hugo_Symbol)), 1:nodes)

    clusterExport(cl, varlist=c("hugo.portions","main.table", "as.data.table", "Function.PSEA.clust", "count.var", "total.count", "exp.matrix","print"),envir=environment())

    #Run
    main.table<-as.data.table(main.table)
    PSEA.Table<-parLapply(cl, hugo.portions, function(x) {
        hugo.table<-main.table[Hugo_Symbol %in% x,]
        hugo.PSEA<-hugo.table[,Function.PSEA.clust(PATIENT) , by="Hugo_Symbol"]
        return (hugo.PSEA)
    })

    #Stop parallelization
    stopCluster(cl)
    print("done with PSEA calculation")

    #Clean up
    PSEA.Table<-do.call(rbind,PSEA.Table)
    PSEA.Table<-as.data.table(merge(as.data.frame(PSEA.Table), as.data.frame(patient.count),
        by="Hugo_Symbol"))
    PSEA.Table<-PSEA.Table[order(PSEA.SCORE, decreasing=T),]

    #Return
    return(PSEA.Table)
}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
exp.matrix<-readRDS(args[2])
output.file<-args[3]
print("opened files")

count.var<-0

main.table<-Function.p.pre.clustering(table.1, exp.matrix)
print ("done pre-clustering")

PSEA.scores<-Function.PSEA(main.table)
print ("done with PSEA")

saveRDS(PSEA.scores, output.file)
print ("done saving")