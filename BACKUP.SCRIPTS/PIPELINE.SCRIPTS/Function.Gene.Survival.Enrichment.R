#Function.Gene.Survival.Enrichment.R - NO RESULTS
#100314
#Function to find genes correlated to disease progression based on gene expression patterns on pathways

Function.RNAseq.Differential.Expression.V2<-function(normalized.matrices.object, target.cancer.samples) {
  #Takes object from function Function.RNAseq.Matrices.Normalization() and target cancer samples to perform differential expression
  
  require(limma)
  require(edgeR)
  require(data.table)
  
  #Get target combined matrix
  cancer.samples.in.matrix<-intersect(target.cancer.samples, normalized.matrices.object$cancer.patients)
  normal.samples<-normalized.matrices.object$normal.patients
  target.combined.matrix<-normalized.matrices.object$combined.matrices[,c(cancer.samples.in.matrix,normal.samples)]
  
  #Get design matrix
  G.design.matrix<-data.frame(G=c(rep("G1", length(cancer.samples.in.matrix)), rep("G0", length(normal.samples))))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Perform differential expression
  G.fit = lmFit(target.combined.matrix, G.design.matrix) #fitting data
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Convert to data.table
  if ("ID" %in% colnames(all.G.fit)) {
    all.G.fit<-as.data.table(all.G.fit)
  } else {
    all.G.fit$ID<-rownames(all.G.fit)
    all.G.fit<-as.data.table(all.G.fit)
  }
  
  #Return topTable
  DIFF.EXP<-all.G.fit
  return(DIFF.EXP)
}

Function.Thresholded.Distance.Matrix<-function(DIFF.EXP, EXP.OBJ, PATH.TABLE,exp.threshold, path.threshold){
    #Obtains thresholded distance table per pathway based on differential expression analysis, gene pathways and normalized matrix object per gene

    require(data.table)

    #Obtain target genes based on exp.threshold
    DIFF.EXP<-DIFF.EXP[abs(logFC)>=exp.threshold, ]
    TARGET.GENES<-as.vector(DIFF.EXP$ID)

    #TEST
    TARGET.GENES<-as.vector(BRCA.DIFF.EXP[abs(logFC)>=1.0,]$ID)
    path.threshold<-5
    EXP.OBJ<-copy(BRCA.EXP.OBJ)
    PATH.TABLE<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/NCI.PATHWAYS/091614.path.gene.rds")

    #Filter path.table [Path, Hugo_Symbol] based on path.threshold and available genes from differential expression analysis
    PATH.TABLE<-PATH.TABLE[Hugo_Symbol %in% TARGET.GENES,]
    PATH.COUNT<-PATH.TABLE[,list(N.HUGO=length(Hugo_Symbol)), by="Path"]
    PATH.COUNT<-PATH.COUNT[N.HUGO>=path.threshold,]
    PATH.TABLE<-PATH.TABLE[Path %in% as.vector(PATH.COUNT$Path),]

    #Obtain matrix containing per pathway patient distances to normal
    CANCER.DIST.PER.PATH<-function(hugos) {

        #Build correlation matrix based on chosen gene
        CORR.MATRIX<-cor(EXP.OBJ$combined.matrices[hugos, c(EXP.OBJ$normal.patients, EXP.OBJ$cancer.patients)], method="spearman")

        #Get distance to normal based on genes
        DISS.MATRIX<-1-CORR.MATRIX
        CANCER.DIST<-sapply(EXP.OBJ$cancer.patients, function(y) median(DISS.MATRIX[y, EXP.OBJ$normal.patients]))
        CANCER.DIST<-as.data.table(CANCER.DIST, keep.rownames=T)
        setnames(CANCER.DIST, c("PATIENT", "DIST.TO.NORMAL"))

        #Return
        return(list(PATIENT=CANCER.DIST$PATIENT, DIST.TO.NORMAL=CANCER.DIST$DIST.TO.NORMAL))
    }

    PATH.DISTANCES<-PATH.TABLE[,CANCER.DIST.PER.PATH(Hugo_Symbol),by=Path]

    #Return
    return(PATH.DISTANCES)
}

Function.Survival.Procurement.Correlation<-function(PATH.DISTANCES. CLINICAL.TABLE) {
    #Produces Correlation of pathways distances to Overall survival and Procurement

    require(data.table)

    #TEST
    CLINICAL.TABLE<-copy(CLINICAL)

    #Filters clinical table
    CLINICAL.TABLE<-CLINICAL.TABLE[,c("PATIENT", "VITAL.STATUS", "OVERALL.SURVIVAL", 
        "DAYS.TO.PROCUREMENT"), with=F]
    CLINICAL.TABLE<-CLINICAL.TABLE[OVERALL.SURVIVAL>=10 & DAYS.TO.PROCUREMENT>=0, ]

    #Introduce clinical survival and procurement information into distances table
    PATH.DISTANCES<-as.data.table(merge(as.data.frame(PATH.DISTANCES), as.data.frame(CLINICAL.TABLE), by="PATIENT"))

    #Obtain distance to clinical correlation
    internal.correlation<-function(distances, survival, proc) {

        SURVIVAL.CORR<-cor.test(distances, survival, method="spearman")
        PROCUREMENT.CORR<-cor.test(distances, proc, method="spearman")

        return(list(SURVIVAL.RHO=SURVIVAL.CORR$estimate, SURVIVAL.P.VAL=SURVIVAL.CORR$p.value,
            PROCUREMENT.RHO=PROCUREMENT.CORR$estimate, PROCUREMENT.P.VAL=PROCUREMENT.CORR$p.value))
    }

    PATH.CORR<-PATH.DISTANCES[, internal.correlation(DIST.TO.NORMAL, OVERALL.SURVIVAL, DAYS.TO.PROCUREMENT) ,by=c("Path", "VITAL.STATUS")]

    #Correct for multiple hypothesis testing
    PATH.CORR$SURVIVAL.P.VAL.ADJ<-p.adjust(PATH.CORR$SURVIVAL.P.VAL, method="fdr")
    PATH.CORR$PROCUREMENT.P.VAL.ADJ<-p.adjust(PATH.CORR$PROCUREMENT.P.VAL, method="fdr")

    #Return
    return(PATH.CORR)
}