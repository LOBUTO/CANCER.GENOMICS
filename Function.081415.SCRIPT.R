#Met.scoring.3

Function.met.score.3<-function(tcga.mut, exp.obj, tcga.clinical, pval.th=0.1, fold.th=1, met.sample.th=5){
  #Calculate metabolite score based on BOTH mutation and expression data. CONTRIBUTION of mutation to phenotype.
  #   NOTE: Contribution is determined based on number of differentially expression genes for samples that have metabolite mutations vs rest of samples in cancer
  #Initial score will not discriminate between GOF/LOF  mutations
  #Scores based on ER-/+ subpopulations
  #NOTE: exp.obj is object containing both normal and cancer matrices
  
  require(parallel)
  require(data.table)
  require(reshape2)
  
  #Set up parallelization
  print ("Importing values for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("tcga.mut", "exp.obj", "tcga.clinical", "pval.th", "fold.th", "met.sample.th", 
                              "as.data.table", "data.table"),envir=environment())
  
  internal.scores<-function(exp.samples){
    
    #Consolidate mutation data
    kegg<-unique(rbind(data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$SUBSTRATE),
                       data.table(Hugo_Symbol=kegg.edges$Hugo_Symbol, MET=kegg.edges$PRODUCT)))
    
    ###Obtain score for each metabolite based on phenotype contribution to exp.samples cohort###
    #Given that there are actually metabolites with associated mutation
    merged.kegg<-merge(kegg, tcga.mut, by="Hugo_Symbol")
    merged.kegg<-merged.kegg[SAMPLE %in% exp.samples,] #Based on SAMPLE COHORT
    
    #Filter for those metabolites that have samples in more than met.sample.th
    merged.kegg[,COUNT.TH:=length(unique(SAMPLE)), by="MET"]
    merged.kegg<-merged.kegg[COUNT.TH>=met.sample.th,]
    all.mets<-unique(merged.kegg$MET)
    
    count<-0
    met.scores<-sapply(all.mets, function(x) {
      
      #Separate those samples affected by metabolite associated mutations and the rest
      met.samples<-unique(merged.kegg[MET==x,]$SAMPLE)
      rest.samples<-setdiff(exp.samples, met.samples)
      
      #Calculate number of differentially expressed genes between the 2 groups (Phenotype contribution)
      gene.pvals<-parApply(cl, exp.obj$tumor, 1, function(y) wilcox.test(y[met.samples], y[rest.samples])$p.value)
      gene.lfc<-parApply(cl, exp.obj$tumor, 1, function(y) log2(mean(y[met.samples])/mean(y[rest.samples])))
      met.table<-data.table(PVAL=gene.pvals, LFC=gene.lfc)
      met.table$PVAL.ADJ<-p.adjust(met.table$PVAL, method="fdr")
      
      #Return score
      met.table<-met.table[!(is.na(LFC)),][!(is.na(PVAL.ADJ)),] #NOTE: Temporary clean up to account for inconsistencies!
      met.score<-mean(met.table$PVAL.ADJ<pval.th & abs(met.table$LFC)>=fold.th)
      
      count<<-count+1
      print(count/length(all.mets))
      return(met.score)
    })
    
    #Clean up and return
    main.table<-data.table(MET=all.mets, MET.SCORE=met.scores)
    return(main.table)
  }
  
  #Only consider silent mutations, for the time being, treat GOF or LOF the same
  tcga.mut<-tcga.mut[MUTATION!="SILENT",]
  
  #Pre-filter expression matrices for non-variable genes (sd==0)
  exp.obj$tumor<-data.matrix(exp.obj$tumor)
  exp.obj$tumor<-exp.obj$tumor[apply(exp.obj$tumor, 1, sd)!=0,]
  
  #Do for both er+/- For those that have both expression and mutation data
  common.samples<-intersect(unique(tcga.mut$SAMPLE), colnames(exp.obj$tumor))
  pos.samples<-intersect(unique(tcga.clinical[ER.STATUS=="Positive",]$SAMPLE), common.samples)
  neg.samples<-intersect(unique(tcga.clinical[ER.STATUS=="Negative",]$SAMPLE), common.samples)
  
  exp.pos<-internal.scores(pos.samples)
  exp.neg<-internal.scores(neg.samples)
  
  #Stop parallelization
  stopCluster(cl)
  
  #Return both
  return(list(EXP.POS=exp.pos, EXP.NEG=exp.neg)) 
}

#Arguments
args<-commandArgs(trailingOnly=T)
tcga.mut<-readRDS(args[1])
exp.obj<-readRDS(args[2])
tcga.clinical<-readRDS(args[3])
output.file<-args[4]
print ("done loading files")

MAIN.OBJ<-Function.met.score.3(tcga.mut, exp.obj, tcga.clinical, pval.th = 0.1, fold.th = 1, met.sample.th = 10)

#Save to output
saveRDS(object = MAIN.OBJ, file = output.file)
print ("Done writing output")