library(data.table)
library(reshape2)
library(parallel)
library(igraph)

Function.mut.bfs<-function(kegg.path.enzyme, kegg.path, kegg.edges, path.common, path=c(), n.layers=2, mut.table, cancer=T, met.edges="constant"){
  #Breadth-first search
  #NOTE: Keep in mind that there should only be one mutation class (SILENT, LOF, GOF) per gene
  #NOTE: met.edges can be "constant" or "degree", a degree weight will be based on all outward degrees of substrate/product/enzyme to each other, a
  #   constant weight will be based on a value of 0.8
  
  require(igraph)
  
  #Obtain mets and genes present in paths
  path.genes<-c()
  path.mets<-c()
  for (p in path){
    path.genes<-union(path.genes, unique(kegg.path.enzyme[grepl(p, PATH, ignore.case = T),]$Hugo_Symbol))
    path.mets<-union(path.mets, unique(kegg.path[grepl(p, DESCRIPTION, ignore.case = T),]$COMPOUND))
  }
  
  #Add edge weight depending on a constant value or degree
  if (met.edges=="degree"){
    kegg.edges[,H.WEIGHT:=1/length(SUBSTRATE), by="Hugo_Symbol"]
    kegg.edges[,S.WEIGHT:=1/length(PRODUCT), by="SUBSTRATE"]
    kegg.edges[,P.WEIGHT:=1/length(Hugo_Symbol), by="PRODUCT"]  
  } else if(met.edges=="constant"){
    kegg.edges[,H.WEIGHT:=0.8, by="Hugo_Symbol"]
    kegg.edges[,S.WEIGHT:=0.8, by="SUBSTRATE"]
    kegg.edges[,P.WEIGHT:=0.8, by="PRODUCT"]  
  }
  
  #Filter edge table for mets and genes 
  net.edges<-kegg.edges[Hugo_Symbol %in% path.genes,]
  net.edges<-net.edges[SUBSTRATE %in% path.mets | PRODUCT %in% path.mets, ]
  
  #Construct bipartite igraph (Integration of metabolite and genes into single graph)
  #NOTE: Adding concept of PRODUCT INHIBITION INTO EDGES: That is, a substrate will eventually inhibit the enzyme that produces it
  #NOTE: Order: Enzyme -> M1 -> M2 -| Enzyme
  #NOTE: Weight of enzyme on product will be based on its degree, analogously, substrate on product will be based on degree contribution of substrate on product
  net.edges<-as.matrix(net.edges)
  bip.edges<-data.table()
  for (i in 1:nrow(net.edges)){
    bip.edges<-rbind(bip.edges, data.table(FROM=net.edges[i,1], TO=net.edges[i,2], WEIGHT=as.numeric(net.edges[i,4])))
    bip.edges<-rbind(bip.edges, data.table(FROM=net.edges[i,2], TO=net.edges[i,3], WEIGHT=as.numeric(net.edges[i,5])))
    bip.edges<-rbind(bip.edges, data.table(FROM=net.edges[i,3], TO=net.edges[i,1], WEIGHT=-as.numeric(net.edges[i,6])))
  }
  
  setkey(bip.edges)
  bip.edges<-unique(bip.edges)
  
  #Add number of layers to pathway if necessary - This adds both genes as well as accompanying metabolites of genes
  path.common<-path.common[WEIGHT!=0,]
  path.common<-path.common[,c("PARENT", "CHILD", "WEIGHT"), with=F]
  setkey(path.common)
  path.common<-unique(path.common)
  setnames(path.common, c("FROM", "TO", "WEIGHT"))
  
  if (n.layers>0){
    for (l in 1:n.layers){
      
      #Add associated gene edges
      path.filtered<-path.common[FROM %in% union(bip.edges$FROM, bip.edges$TO) | TO %in% union(bip.edges$FROM, bip.edges$TO),]
      bip.edges<-rbind(bip.edges, path.filtered)
      
      #Add associated metabolites
      add.genes<-union(path.filtered$FROM, path.filtered$TO)
      add.mets<-kegg.edges[Hugo_Symbol %in% add.genes,]
      
      add.mets<-as.matrix(add.mets)
      path.bip.edges<-data.table()
      for (i in 1:nrow(add.mets)){
        path.bip.edges<-rbind(path.bip.edges, data.table(FROM=add.mets[i,1], TO=add.mets[i,2], WEIGHT=as.numeric(add.mets[i,4])))
        path.bip.edges<-rbind(path.bip.edges, data.table(FROM=add.mets[i,2], TO=add.mets[i,3], WEIGHT=as.numeric(add.mets[i,5])))
        path.bip.edges<-rbind(path.bip.edges, data.table(FROM=add.mets[i,3], TO=add.mets[i,1], WEIGHT=-as.numeric(add.mets[i,6])))
      }
      bip.edges<-rbind(bip.edges, path.bip.edges)
      setkey(bip.edges)
      bip.edges<-unique(bip.edges)
    }
  }
  setkey(bip.edges)
  bip.edges<-unique(bip.edges)
  
  #NOTE: Remove edges that are not both repressed and activate by at least 1 node - HOLD
  #   bip.edges[,N.ACT:=sum(WEIGHT>0), by="TO"]
  #   bip.edges[,N.REP:=sum(WEIGHT<0), by="TO"]
  #   bip.edges<-bip.edges[N.ACT>0 & N.REP>0,]
  #   bip.edges$N.ACT<-NULL
  #   bip.edges$N.REP<-NULL
  
  #Apply mutation status - Trimming of regulatory edges in cancer if we cannot diffuse through them
  mut.genes<-unique(mut.table$Hugo_Symbol)
  print (mut.genes)
  mut.genes<-mut.genes[mut.genes %in% unique(bip.edges$FROM)]
  print (mut.genes)
  
  if (cancer==T){
    
    for (mut in mut.genes){
      mut.stat<-unique(mut.table[Hugo_Symbol==mut,]$MUTATION)
      
      if (mut.stat=="GOF"){
        #Remove negative regulation onto gene - Add self zero for next function
        bip.edges<-bip.edges[!(TO==mut & WEIGHT<0),]
        
      } else if (mut.stat=="LOF"){
        #Remove positive regulation onto gene - Add self zero for next function
        bip.edges<-bip.edges[!(TO==mut & WEIGHT>0),]
      }
    }
  }
  
  #Create graph object
  main.graph<-graph.data.frame(bip.edges, directed = T)
  V(main.graph)$color<-ifelse(V(main.graph)$name %in% path.mets, "green",
                              ifelse(V(main.graph)$name %in% path.genes, "purple", "white"))
  V(main.graph)$shape<-ifelse(V(main.graph)$name %in% union(kegg.edges$SUBSTRATE, kegg.edges$PRODUCT), "csquare", "circle")
  
  ####Obtain BFS order list per mutated gene####### - RETURN AT THE END
  bfs.list<-lapply(mut.genes, function(x) {
    bfs.dist<-graph.bfs(main.graph, x, neimode = "out", unreachable = F, dist=T)[["dist"]]
    bfs.dist<-lapply(1:max(bfs.dist[!is.na(bfs.dist)]), function(y) V(main.graph)$name[which(bfs.dist==y)])
    return(bfs.dist)
  })
  names(bfs.list)<-mut.genes
  
  #Return
  return(list(BFS=bfs.list, GRAPH=main.graph, EDGES=bip.edges))
}

Function.boolnet.2<-function(bfs.list, net.edges, net.graph){
  #Apply BOOLEAN NETWORK Algorithm on output from Function.mut.bfs() using non-dysregular edges
  #NOTE: Keep in mind that the BFS layers should not contain the original mutated node as the first layer, since its parent nodes are not active yet,
  #   that is to say, that the mutated node is the first parent layer
  #NOTE: net.edges should have been filtered out from sink nodes (EXCEPT MUTATED NODE)
  
  require(reshape2)
  
  #Mix bfs lists into single ordered gene list
  #NOTE: MAKE SURE we actually have a bfs.list to go through
  
  if (length(bfs.list) > 0){
    
    mut.genes<-names(bfs.list)
    maximum.layers<-max(sapply(bfs.list, function(x) length(x)))
    ordered.bfs.genes<-c()
    for (m in 1:maximum.layers){
      for (g in mut.genes){
        if(length(bfs.list[[g]])>=m){
          ordered.bfs.genes<-c(ordered.bfs.genes, bfs.list[[g]][[m]])
        }
      }
    }
    ordered.bfs.genes<-unique(ordered.bfs.genes)
    
    #Complete edges tables with "0" weights to be able to construct transition matrix
    net.edges<-rbind(net.edges, data.table(FROM=ordered.bfs.genes, TO=ordered.bfs.genes, WEIGHT=0))
    net.edges<-net.edges[,list(WEIGHT=sum(WEIGHT)), by=c("FROM", "TO")]
    
    #Construct activator and repressor transition matrices
    act.table<-net.edges[WEIGHT>=0,]
    rep.table<-net.edges[WEIGHT<=0,]
    act.trans<-acast(act.table, TO~FROM, fill = 0, value.var = "WEIGHT")
    rep.trans<-abs(acast(rep.table, TO~FROM, fill = 0, value.var = "WEIGHT"))
    
    #Initialize boolnet.matrix (with 0.1)
    all.nodes<-union(net.edges$FROM, net.edges$TO)
    bool.matrix<-matrix(0, ncol=1, nrow=length(all.nodes), dimnames = list(all.nodes, "0" ))
    bool.matrix[,"0"]<-0.1
    
    #Iterate through layers for n.iteratinons until stable
    delta<-Inf
    iter<-1
    
    print ("optimizing...")
    while (delta>0.01){
      
      #Add entry in boolean matrix as previous iteration, non-updated genes will stay the same
      cnames<-c(colnames(bool.matrix), as.character(iter))
      bool.matrix<-cbind(bool.matrix, bool.matrix[,as.character(iter -1)])
      colnames(bool.matrix)<-cnames
      
      #Update entry 
      for (g in ordered.bfs.genes){
        parent.act<-act.table[TO==g,]$FROM
        parent.rep<-rep.table[TO==g,]$FROM
        act.score<- 1- prod(1 - (as.numeric(bool.matrix[parent.act ,as.character(iter-1)]) * as.numeric(act.trans[g, parent.act])))
        rep.score<-    prod(1 - (as.numeric(bool.matrix[parent.rep ,as.character(iter-1)]) * as.numeric(rep.trans[g, parent.rep])))
        bool.matrix[g, as.character(iter)]<-act.score * rep.score
      }
      
      #Update delta and iteration count (only for those that are actually changing)
      delta<-mean(abs(as.numeric(bool.matrix[ordered.bfs.genes, as.character(iter)]) - 
                        as.numeric(bool.matrix[ordered.bfs.genes, as.character(iter-1)])))
      print (c(iter, delta, colMeans(bool.matrix)))
      iter<-iter+1
    }
    
    #Return 
    return(list(BOOL.MATRIX=bool.matrix, GRAPH=net.graph))
    
  } else {
    return(list(BOOL.MATRIX=c(), GRAPH=net.graph))
  }
}

Function.master.boolnet.cancer<-function(tang.matrix, tcga.mut, paths=c(), layers=1){
  #Master of CANCER-boolnet.2: Create a boolnet.2 instance for each individual patient having GOF and LOF mutations only
  #NOTE: This will run bfs and boolnet on all WITH EXCLUSION of cancer edges, bfs will be run from the mutation starting points
  
  #Filter cancer samples from tang.matrix
  tcga.mut<-tcga.mut[SAMPLE %in% colnames(tang.matrix),]
  
  #Filter for GOF and LOF mutations only
  tcga.mut<-tcga.mut[MUTATION!="SILENT",]
  tcga.samples<-unique(tcga.mut$SAMPLE)
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table","tcga.mut", "kegg.path.enzyme","kegg.path", "kegg.edges","path.cancer.breast",
                              "paths", "layers",  "Function.mut.bfs", "Function.boolnet.2", "tcga.samples", "setkey", "setnames") ,envir=environment())
  print ("done exporting variables for parallelization")
  
  #Apply CANCER bfs and boolnet for each individual (parallelize if necessary)
  master.list<-parLapply(cl, tcga.samples, function(x) { 
    print (x)
    
    #Apply CANCER-BFS (path.cancer.breast) - Last used with constant edges
    mut.table<-tcga.mut[SAMPLE==x,]
    sample.bfs<-Function.mut.bfs(kegg.path.enzyme, kegg.path, kegg.edges, path.cancer.breast, path = paths , n.layers = layers, mut.table ,cancer = T, met.edges = "constant")
    
    #Apply boolnet.2
    sample.boolnet<-Function.boolnet.2(sample.bfs$BFS, sample.bfs$EDGES, sample.bfs$GRAPH)
    
    #Return
    return(sample.boolnet)
  })
  names(master.list)<-tcga.samples
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  #Return
  return(master.list)
}



#Arguments
args<-commandArgs(trailingOnly=T)
tang.matrix<-readRDS(args[1])
tcga.mut<-readRDS(args[2])
kegg.path.enzyme<-fread(args[3], header=T, sep="\t", stringsAsFactors = F)
kegg.path<-fread(args[4], header=T, sep="\t")
kegg.edges<-readRDS(args[5])
path.cancer.breast<-readRDS(args[6])

output.file<-args[7]
print ("done loading files")

MAIN.OBJ<-Function.master.boolnet.cancer(tang.matrix, tcga.mut,  c("Glycolysis", "TCA"), 2)

#Save to output
saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")