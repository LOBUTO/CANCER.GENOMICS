#Function.mut.dist.over.path
#091614
#Given an expression matrix, table.1 and path.gene table, it calculates the Sum of intra-patient correlation with a given mutation in table.1 in a correlation matrix where all patients are included. Correlation matrix is derived from expression matrix using 1-Corr.matrix

Function.mut.correlation<-function(path.genes, target.patients){
    cat(count.var/total.count)
    #Obtains a correlation score per path

    #obtain target matrix
    pre.target.matrix<-exp.matrix$combined.matrices[,exp.matrix$cancer.patients]
    target.genes<-intersect(rownames(pre.target.matrix), path.genes)
    target.matrix<-pre.target.matrix[target.genes,]

    #Get correlation
    corr.matrix<-cor(target.matrix, method="spearman")

    #Get correlation score for target patients
    max.score<-length(target.patients)^2
    corr.sum<-sum(corr.matrix[target.patients, target.patients])
    path.score<-corr.sum/max.score

    #Update count
    count.var<<-count.var+1

    #Return dissimilarity matrix
    return(list(Score=path.score))
}

Function.target.distances<-function(table.1, target.gene, exp.matrix, gene.path) {
    #Produces data.table of distances per pathway found in gene.path table for patients containing target.mutation

    library(data.table)

    #Gets patients of interested from table.1
    target.patients<-unique(as.vector(table.1[Hugo_Symbol==target.gene,]$PATIENT))

    #Keep those that are actually in diss.matrix
    target.patients<-target.patients[target.patients %in% exp.matrix$cancer.patients]

    #Filter gene.path table for paths that have at least 10 genes
    gene.path<-gene.path[Hugo_Symbol %in% rownames(exp.matrix$combined.matrices),]
    gene.path.count<-gene.path[,list(gene.count=length(Hugo_Symbol)), by="Path"]
    gene.path.count<-gene.path.count[gene.count>10,]
    gene.path<-gene.path[Path  %in% unique(as.vector(gene.path.count$Path)),]

    #Total paths to go through
    total.count<<-length(unique(as.vector(gene.path$Path)))

    #Set up parallelization
    library(parallel)
    nodes<-detectCores()
    cl<-makeCluster(nodes)
    setDefaultCluster(cl)

    path.portions<-split(unique(as.vector(gene.path$Path)), 1:nodes)

    clusterExport(cl, varlist=c("path.portions", "gene.path", "as.data.table", "target.patients", "Function.mut.correlation", "count.var", "total.count", "exp.matrix"),envir=environment())

    #Calculate
    path.results<-parLapply(cl,path.portions, function(x) {
        target.table<-as.data.table(gene.path)[Path %in% x,]
        table.portion<-target.table[ ,Function.mut.correlation(Hugo_Symbol, target.patients), by="Path"]
        return(table.portion)
    })
    stopCluster(cl)

    #Clean up and return [Path, Score]
    path.results<-do.call(rbind, path.results)
    path.results<-as.data.table(merge(as.data.frame(path.results), as.data.frame(gene.path.count), by="Path"))
    path.results<-path.results[order(Score, decreasing=T),]
    return(path.results)
}

args<-commandArgs(trailingOnly=T)
table.1<-readRDS(args[1])
exp.matrix<-readRDS(args[2])
gene.path<-readRDS(args[3])
target.gene<-as.character(args[4])
output.file<-args[5]
print("opened files")

count.var<-0

path.scores<-Function.target.distances(table.1, target.gene, exp.matrix, gene.path)
print ("done with path scores")

saveRDS(path.scores, output.file)
print ("done saving")