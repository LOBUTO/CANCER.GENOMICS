library(org.Hs.eg.db)
library(annotate)

Function.BRCA.SUBTYPE<-function(brca.normalized.obj, version=1){
  #Takes normalized expression object from breast cancer and returns subtype score and assigns subtype to each patient based on maximum score
  
  require(genefu)
  require(data.table)
  
  #Choose pam50 method
  if (version==1){
    data(pam50)
    pam<-copy(pam50)
  } else if(version==2){
    data(pam50.scale)
    pam<-copy(pam50.scale)
  } else if(version==3){
    data(pam50.robust)
    pam<-copy(pam50.robust)
  }
  
  #Obtain target pam gene in expression matrix
  target.genes<-intersect(rownames(pam$centroids), rownames(brca.normalized.obj$combined.matrices))
  
  #Simplify pam and expression matrices
  simplified.pam<-pam$centroids[target.genes,]
  simplified.exp<-brca.normalized.obj$combined.matrices[target.genes, brca.normalized.obj$cancer.patients]
  
  #Obtain predicted type based on correlation
  BRCA.SCORES<-cor(simplified.exp, simplified.pam, method="spearman")
  BRCA.SCORES<-as.data.frame(BRCA.SCORES)
  BRCA.SCORES$PATIENT<-rownames(BRCA.SCORES)
  BRCA.SCORES$TYPE<-colnames(BRCA.SCORES)[max.col(BRCA.SCORES[,1:5])]
  
  #Clean up and return
  BRCA.SCORES<-as.data.table(BRCA.SCORES)
  BRCA.SCORES<-BRCA.SCORES[order(TYPE),]
  return(BRCA.SCORES)
}
hist(BRCA.EXP$combined.matrices[,BRCA.EXP$cancer.patients])

#Build Claudin low predictor
Function.CL.CENT<-function(UNC.file, unprocessed.cl.list, type="unc", filter=c()){
  
  #Process cl gene list
  cl.list<-fread(unprocessed.cl.list, header=F, stringsAsFactors=F)
  cl.list$V2<-sapply(cl.list$V2, function(x) unlist(strsplit(x, " "))[1])
  cl.list<-cl.list[!(is.na(V2)),]
  cl.list<-unique(cl.list)
  
  #Filter?
  if (length(filter)>0){
    cl.list<-cl.list[V1 %in% filter,]
  }
  print(dim(cl.list))
  
  if (type=="unc"){
    
    setnames(cl.list, c("subtype", "Hugo_Symbol"))
    #Load UNC337 matrix and filter for CL list
    #UNC.file<-"DATABASES/CANCER_DATA/PEROU/DifferentiationPredictor/UNC337arraydata_imputedCollapsed.txt"
    UNC337<-fread(UNC.file, header=T, sep="\t", stringsAsFactors=F, skip=4)
    UNC337<-data.frame(UNC337)
    UNC337.MATRIX<-data.matrix(UNC337[,2:ncol(UNC337)])
    UNC337<-cbind(UNC337[,1,drop=F], data.frame(scale(UNC337.MATRIX - apply(UNC337.MATRIX, 1, median))))
    UNC337<-merge(UNC337, cl.list, by="subtype")    
    
  } else if (type=="9cell") {
    
    setnames(cl.list, c("Group", "Hugo_Symbol"))
    #Load 9CELL matrix and filter for CL list
    UNC337<-fread(UNC.file, skip=1, drop=2)
    UNC337<-data.frame(UNC337)
    UNC337.MATRIX<-data.matrix(UNC337[,2:ncol(UNC337)])
    UNC337<-cbind(UNC337[,1,drop=F], data.frame(scale(UNC337.MATRIX - apply(UNC337.MATRIX, 1, median))))
    UNC337<-merge(UNC337, cl.list, by="Group")
  }
  
  UNC337<-data.frame(UNC337,row.names=as.vector(UNC337$Hugo_Symbol))
  UNC337$subtype<-NULL
  UNC337$Group<-NULL
  UNC337$Hugo_Symbol<-NULL
  UNC337<-data.matrix(UNC337)  
  
  #Apply scaling according to paper
  UNC337<-UNC337-apply(UNC337, 1, median)
  
  #Obtain centroids
  claudin.samples<-colnames(UNC337)[grepl("Claudin", colnames(UNC337))]
  other.samples<-setdiff(colnames(UNC337), claudin.samples)
  
  centroids.table<-data.table(Hugo_Symbol=rownames(UNC337), 
                              CLAUDIN.LOW=apply(UNC337, 1, function(x) mean(x[claudin.samples])),
                              OTHERS=apply(UNC337, 1, function(x) mean(x[other.samples])))
  
  #Return
  return(centroids.table)
}

CL.CENTROIDS<-Function.CL.CENT("DATABASES/CANCER_DATA/PEROU/DifferentiationPredictor/UNC337arraydata_imputedCollapsed.txt",
                               "DATABASES/CANCER_DATA/PEROU/CLAUDIN.LOW.LIST", type="unc")
CL.CENTROIDS<-Function.CL.CENT("DATABASES/CANCER_DATA/PEROU/T.E.9CELL-LINE.txt",
                               "DATABASES/CANCER_DATA/PEROU/CLAUDIN.LOW.LIST", type="9cell")

CELL.9<-data.frame(fread("DATABASES/CANCER_DATA/PEROU/T.E.9CELL-LINE.txt", skip=1, drop=2))
rownames(CELL.9)<-CELL.9$Group
CELL.9$Group<-NULL
CELL.9<-as.matrix(CELL.9)
CELL.9<-normalizeBetweenArrays(CELL.9, method="quantile")

claudin.samples<-colnames(CELL.9)[grepl("Claudin", colnames(CELL.9))]
other.samples<-setdiff(colnames(CELL.9), claudin.samples)
CELL.9.DIFF<-data.table(EntrezGene.ID=rownames(CELL.9),
                        P.VAL=apply(CELL.9, 1, function(x) wilcox.test(x[claudin.samples], x[other.samples], paired=F)$p.value))
CELL.9.DIFF$P.VAL.ADJ<-p.adjust(CELL.9.DIFF$P.VAL, method="fdr")
CELL.9.DIFF<-CELL.9.DIFF[order(P.VAL.ADJ),]
CELL.9.DIFF[P.VAL.ADJ<0.05,]

Function.CL.PRED<-function(cl.centroids, exp.matrix, scale=T){
  
  #Scale testing matrix?
  if (scale==T){
    exp.matrix<-exp.matrix-apply(exp.matrix,1,median)
    exp.matrix<-scale(exp.matrix) 
  }
  
  #Filter for cl gene predictors
  cl.genes<-intersect(rownames(exp.matrix), unique(cl.centroids$Hugo_Symbol))
  exp.matrix<-exp.matrix[cl.genes,]
  
  #Obtain centroid matrix
  cl.matrix<-data.frame(cl.centroids, row.names=as.vector(cl.centroids$Hugo_Symbol))
  cl.matrix$Hugo_Symbol<-NULL
  cl.matrix<-as.matrix(cl.matrix)
  cl.matrix<-cl.matrix[cl.genes,]
    
  #Obtain eucledian distance
  cl.dist<-apply(exp.matrix, 2, function(x) sqrt(sum((x-as.vector(cl.matrix[,"CLAUDIN.LOW"]))^2)) )
  others.dist<-apply(exp.matrix, 2, function(x) sqrt(sum((x-as.vector(cl.matrix[,"OTHERS"]))^2)) )
  
  #cl.dist<-apply(exp.matrix, 2, function(x) cor(x,as.vector(cl.matrix[,"CLAUDIN.LOW"]), method = "spearman") ) 
  #others.dist<-apply(exp.matrix, 2, function(x) cor(x,as.vector(cl.matrix[,"OTHERS"]), method = "spearman") )
  
  cl.scores<-data.table(SAMPLE=colnames(exp.matrix), CL.DIST=cl.dist, OTHERS.DIST=others.dist)
  
  #Apply score based on highest eucledian distance
  cl.scores$PRED<-ifelse( cl.scores$CL.DIST<cl.scores$OTHERS.DIST , "CLAUDIN.LOW", "OTHER")
  
  #Clean up and return
  cl.scores<-cl.scores[order(PRED),]
  return(cl.scores)
}

mouse.cl.pred<-Function.CL.PRED(CL.CENTROIDS, mouse , scale=T)
mouse.cl.pred[grepl("claudin", SAMPLE, ignore.case=T),]
mouse.cl.pred[PRED=="CLAUDIN.LOW",][order(CL.DIST, decreasing=T),]

BRCA.TCGA.PRED<-Function.CL.PRED(CL.CENTROIDS, brca.test)
BRCA.TCGA.PRED[PRED=="CLAUDIN.LOW",]
hist(BRCA.TCGA.PRED$CL.DIST/BRCA.TCGA.PRED$OTHERS.DIST)
BRCA.TCGA.PRED[SAMPLE=="TCGA.A8.A08H.01A",]

BRCA.TCGA.PAM50<-Function.BRCA.SUBTYPE(BRCA.EXP)
BRCA.TCGA.PAM50[TYPE %in% c("Normal"),]

BRCA.TCGA.PRED[PRED=="CLAUDIN.LOW",]$SAMPLE %in% BRCA.TCGA.PAM50[TYPE %in% c("Basal"),]$PATIENT

mouse<-fread("DATABASES/CANCER_DATA/PEROU/DifferentiationPredictor/Jason-122-Mouse.txt", header=T, sep="\t", stringsAsFactors=F, skip=1)
mouse<-mouse[3:nrow(mouse),]
mouse<-data.frame(mouse)
setnames(mouse, c("subtype", colnames(mouse)[2:ncol(mouse)]))
mouse$subtype<-getSYMBOL(mouse$subtype, data = 'org.Hs.eg')
mouse<-mouse[!is.na(mouse$subtype),]
mouse<-data.frame(mouse,row.names = 1)
mouse<-data.matrix(mouse)

cl.list<-fread("DATABASES/CANCER_DATA/PEROU/CLAUDIN.LOW.LIST.ALL", header=F, stringsAsFactors=F)
cl.list$V2<-sapply(cl.list$V2, function(x) unlist(strsplit(x, " "))[1])
cl.list<-cl.list[!(is.na(V2)),]
cl.list<-unique(cl.list)
setnames(cl.list, c("subtype", "Hugo_Symbol"))

mouse<-merge(mouse, cl.list, by="subtype")
mouse<-data.frame(mouse, row.names=as.vector(mouse$Hugo_Symbol))
mouse$Hugo_Symbol<-NULL
mouse$subtype<-NULL
mouse<-data.matrix(mouse)

dim(mouse)
mouse[1:3,1:3]

