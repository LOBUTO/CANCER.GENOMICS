#DECAY MERGING FUNCTION - PLUS OTHER FUNCTIONS
#We have one parameter to model: The decay rate lambda

Function.Maf.Coverage.Plot<-function(gene.list, exons.table, thousand.table) {
  require(data.table)
  
  #Filter gene.list for those that we actually have data for
  print (length(gene.list))
  gene.list<-gene.list[gene.list %in% unique(exons.table$Hugo_Symbol)]
  print (length(gene.list))
  
  gene.table<-list()
  no.gene.table<-c()
  
  for (record in gene.list){
    print (record)
    
    #Create target matrix
    gene.exons.table<-exons.table[Hugo_Symbol==record,]
    start<-unique(gene.exons.table$START)
    end<-unique(gene.exons.table$END)
    chrom<-unique(gene.exons.table$Chrom)
    target.matrix<-thousand.table[Chrom==chrom & Position>start & Position<end,]
    target.matrix<-target.matrix[order(Position),]
    
    #Check if we have actualy records for gene in 1000G data
    if (nrow(target.matrix)>1){
      
      #Clean target matrix for position duplicates and aggregate(sum) MAFs contained in duplicates
      target.matrix<-aggregate(MAF~Position,data=target.matrix,FUN=sum)
      
      #Estimate data coverage
      coverage<-Function.Maf.Exon.Coverage(gene.exons.table, target.matrix, gene=record)
      
      #Add to table
      gene.table<-c(gene.table, list(c(record, coverage$gene.coverage, coverage$exon.coverage)))  
      
    } else if (nrow(target.matrix)==1) {
      
      #Estimate data coverage
      coverage<-Function.Maf.Exon.Coverage(gene.exons.table, target.matrix, gene=record)
      
      #Add to table
      gene.table<-c(gene.table, list(c(record, coverage$gene.coverage, coverage$exon.coverage)))  
      
    } else {
       no.gene.table<-c(no.gene.table, record)
    }
  }
  
  gene.table<-as.data.table(do.call(rbind, gene.table))
  setnames(gene.table, c("Hugo_Symbol", "gene.coverage", "exon.coverage"))
  gene.table$gene.coverage<-as.numeric(gene.table$gene.coverage)
  gene.table$exon.coverage<-as.numeric(gene.table$exon.coverage)
  
  #Print not found genes
  if (length(no.gene.table)>0){
    cat ("No MAF data found for the following genes:", no.gene.table )
  }
  
  #Return
  return(gene.table)
}

GENE.COVERAGE<-Function.Maf.Coverage.Plot(c("TP53", "TTN", "PIK3CA", "BRCA1", "BRCA2", "MUC4", "CDH1", sample(unique(EXONS$Hugo_Symbol), 100) ),EXONS,THOUSAND)

ggplot(GENE.COVERAGE, aes(gene.coverage, exon.coverage)) + geom_point() + theme.format

Function.Maf.Exon.Coverage<-function(exons.table, thousand.processed.table, gene){
  #Gives coverage of maf found in thousand genome data with respect to all annoted gene and all exons
  
  require(data.table)
  require(IRanges)
  
  #Calculate total exon nt count
  gene.exon<-exons.table[Hugo_Symbol==gene,]
  gene.ranges<-IRanges(gene.exon$FEAT_START, gene.exon$FEAT_END)
  gene.ranges<-reduce(gene.ranges)
  gene.exon.nt<-sum(end(gene.ranges) - start(gene.ranges))
  
  #Calculate total gene nt count
  gene.total.nt<-unique(gene.exon$END)-unique(gene.exon$START)
  
  #Calculate coverage
  mafs.in.exon<-apply(gene.exon[,c("FEAT_START", "FEAT_END"), with=F], 1, 
                      function(x) as.vector(thousand.processed.table$Position) %in% x[1]:x[2])
  
  mafs.in.exon<-apply(mafs.in.exon, 1, sum)
  mafs.in.exon<-sum(mafs.in.exon>0)
  
  exon.coverage<-mafs.in.exon/gene.exon.nt
  gene.coverage<-nrow(thousand.processed.table)/gene.total.nt
  
  #Return as list
  return(list(exon.coverage=exon.coverage, gene.coverage=gene.coverage))
}

Function.Exp.Decay<-function(No, expansion, lambda=NULL, phastcon=NULL, fixed.lambda=T, a=2){
  #"a" is used only when a phastcons score is provided and refers to the curvature in the inverse relationship of lambda and phastcon.
  #   This relationship is modeled by lambda=a/(a+phastcon)
  #"No" is the minor allele frequency obtained from the thousand geneome data and that we are trying to use to propagate maf data to neighbor sites with no data
  #"expansion.l" and "expansion.r" are the length of the positions we are trying to "expand" our data to. left(l) and right(r)
  
  require(data.table)
  
  #First use if fixed lambda is provided, otherwise use phastcon score
  if (fixed.lambda==T){
    decay<-sapply(1:expansion, function(x) No*exp(-lambda*x))
    #decay.l<-sapply(1:expansion.l, function(x) No*exp(-lambda*x))
  } else{
    decay<-sapply(1:expansion, function(x) No*exp(-x*(phastcon^a/(1+phastcon^a))))
    #decay.l<-sapply(1:expansion.l, function(x) No*exp(-x*(a/(2*a+phastcon))))
  } 
  
  #Return
  return(decay)
}

Function.Pairwise.Decay<-function(thousand.table, exons.table, gene, lambda, phastcons.table=NULL, a.set=2, min.threshold=T, filter=F) {
  
  #Implemented weighted means for opposite decaying functions
  #Implemented multiple maf factors (need to remove from final version)
  
  require(data.table)
  
  #Create target matrix
  exons.table<-exons.table[Hugo_Symbol==gene,]
  start<-unique(exons.table$START)
  end<-unique(exons.table$END)
  chrom<-unique(exons.table$Chrom)
  target.matrix<-thousand.table[Chrom==chrom & Position>start & Position<end,]
  
  #Filter for low count MAFs?
  if (filter==T){
    pre.target.matrix<-pre.target.matrix[MAF>=0.000199682,]
  }
  
  #Check if we have data in 1000G data, if not, stop function - (12.14.14)
  if (nrow(target.matrix)<1){
    return (paste(gene, "not found in 1000G data"))
  }
  
  target.matrix<-target.matrix[order(Position),]
  
  #Clean target matrix for position duplicates and aggregate(sum) MAFs contained in duplicates
  target.matrix<-aggregate(MAF~Position,data=target.matrix,FUN=sum)
  print ("Done cleaning up target matrix")
  
  #Estimate data coverage - TEMPORAL!!!! (12.13.14)
  #coverage<-Function.Maf.Exon.Coverage(exons.table, target.matrix, gene=gene)
  #print (coverage)
  
  #Get minimum standard maf for region based on available data (has to be done before modeling)
  #min.maf<-min(target.matrix$MAF)
  #print (min.maf)
  
  ######PHASTCONS###### - DYNAMIC
  if (!is.null(phastcons.table)){
    print ("Calculating")
    
    phastcons.table<-phastcons.table[Chrom==as.character(chrom) & Position>start & Position<end,]
    print ("Done processing phastcons")
    target.matrix<-merge(target.matrix, phastcons.table, by="Position")
    print ("Done integrating target and phastcons")
    
    #Prep storage vectors
    POS<-as.vector(target.matrix$Position)
    MAF<-as.vector(target.matrix$MAF)
    PHAST<-as.vector(target.matrix$Score)
    EXPANSION.MEAN<-c()
    
    #Loop through a.factors -REMOVE, DONE FOR TESTING MULTIPLE a.factors only
    a.list<-list()
    for (a.factor in a.set) {
      
      #Construct left decay first (before first maf data)
      decay.vector<-rev(Function.Exp.Decay(MAF[1], POS[1]-start, fixed.lambda=F, phastcon=PHAST[1], a=a.factor))
      
      #Check that we actually have more than 1 position for middle loop
      if (length(POS)>1){
        #Then middle
        for (n in 1:(length(POS)-1)){
          
          #If they are separated by 1 or not, either way store maf value
          decay.vector<-c(decay.vector, MAF[n])
          
          if (POS[n]!=(POS[n+1]-1)){
            
            #If we are not separated by 1, then apply extension - IMPLEMENTED WEIGHTED MEANS (121014)
            decay.expansion<-POS[n+1]-POS[n]-1
            EXPANSION.MEAN<-c(EXPANSION.MEAN, decay.expansion)
            
            right.decay<-Function.Exp.Decay(MAF[n], decay.expansion, fixed.lambda=F, phastcon=PHAST[n], a=a.factor)
            weighted.right.decay<-right.decay*(decay.expansion:1)/decay.expansion
            
            left.decay<-rev(Function.Exp.Decay(MAF[n+1], decay.expansion, fixed.lambda=F, phastcon=PHAST[n], a=a.factor))
            weighted.left.decay<-left.decay*(1:decay.expansion)/decay.expansion
            
            decay.vector<-c(decay.vector, colSums(rbind(weighted.right.decay, weighted.left.decay)))
          } 
        }
      } 
 
      #Then right tail
      decay.vector<-c(decay.vector, tail(MAF,1))
      last.vector<-Function.Exp.Decay(tail(MAF,1), end-tail(POS,1), fixed.lambda=F, phastcon=tail(PHAST,1), a=a.factor)
      decay.vector<-c(decay.vector,last.vector)
      
      #Then append to a.list
      a.list<-c(a.list,list(decay.vector))
    }
    
    names(a.list)<-as.character(a.set)
    
    #POST-PROCESS
    #Loop through a.factors to assign separate lists
    multiple.decay.table<-list()
    for (a.factor in a.set){
      decay.vector<-a.list[[as.character(a.factor)]]
      
      #Replace all values in decay function that are equal to zero by the minimum value in the decay vector - FIXES BRCA1 problem of low comparisson
      decay.vector[decay.vector==0]<-min(decay.vector[decay.vector>0])
      
      #Match up to actual positions
      position.vector<-start:end
      decay.table<-data.table(Chrom=chrom, Position=position.vector, MAF=decay.vector)
      
      #Append to decay table
      decay.table$a.factor<-a.factor
      multiple.decay.table[[as.character(a.factor)]]<-decay.table
    }
    
    multiple.decay.table<-do.call(rbind,multiple.decay.table)
    multiple.decay.table<-as.data.table(multiple.decay.table)
    
    #Add median expansion as reference
    multiple.decay.table$expansion.median<-median(EXPANSION.MEAN)
    
    #Return
    print (median(EXPANSION.MEAN))
    return(multiple.decay.table)
    
    ########### NON-PHASTCON ######## - STATIC
  } else {
    print ("Calculating")
    
    #Prep storage vectors
    POS<-as.vector(target.matrix$Position)
    MAF<-as.vector(target.matrix$MAF)
    EXPANSION.MEAN<-c()
    
    #Construct left decay first (before first maf data)
    decay.vector<-rev(Function.Exp.Decay(MAF[1], POS[1]-start, lambda=lambda))
    
    #Check that we actually have data for more than one position before doing loop
    if (length(POS)>1){
      #Then middle
      for (n in 1:(length(POS)-1)){
        
        #If they are separated by 1 or not, either way store maf value
        decay.vector<-c(decay.vector, MAF[n])
        
        if (POS[n]!=(POS[n+1]-1)){
          
          #If we are not separated by 1, then apply extension
          decay.expansion<-POS[n+1]-POS[n]-1
          EXPANSION.MEAN<-c(EXPANSION.MEAN, decay.expansion)
          
          right.decay<-Function.Exp.Decay(MAF[n], decay.expansion, lambda=lambda)
          weighted.right.decay<-right.decay*(decay.expansion:1)/decay.expansion
          
          left.decay<-rev(Function.Exp.Decay(MAF[n+1], decay.expansion, lambda=lambda))
          weighted.left.decay<-left.decay*(1:decay.expansion)/decay.expansion
          
          decay.vector<-c(decay.vector, colSums(rbind(weighted.right.decay, weighted.left.decay)))
        } 
      }
    } 
    
    #Then right tail
    decay.vector<-c(decay.vector, tail(MAF,1))
    last.vector<-Function.Exp.Decay(tail(MAF,1), end-tail(POS,1), lambda=lambda)
    decay.vector<-c(decay.vector,last.vector)
    
    #POST - PROCESS
    #Replace all values in decay function that are equal to zero by the minimum value in the decay vector - FIXES BRCA1 problem of low comparisson
    decay.vector[decay.vector==0]<-min(decay.vector[decay.vector>0])
    
    #Match up to actual positions
    position.vector<-start:end
    decay.table<-data.table(Chrom=chrom, Position=position.vector, MAF=decay.vector)
    
    #Add median expansion as reference
    decay.table$expansion.median<-median(EXPANSION.MEAN)
    
    #Return
    print (median(EXPANSION.MEAN))
    return(decay.table) 
  }
}

Function.maf.bayes<-function(thousand.amaf, tcga.maf, gene,cancer.prob=0.12) {
  #Calculate bayesian probability per site using background artifical maf (amaf) applied to tcga samples
  #The default cancer probability is for breast cancer and was obtained from http://www.breastcancer.org/symptoms/understand_bc/risk/understanding
  #This a per gene calculation, make sure tcga maf chrom and position coordinates match to entered thousand.amaf!!
  
  require(data.table)
  
  #Set up tables
  tcga.maf<-tcga.maf[Hugo_Symbol==gene, c("Start_Position", "Chrom", "MUT.FREQ", "Sample", "Hugo_Symbol"), with=F]
  setnames(tcga.maf, c("Position", "Chrom", "MUT.FREQ", "Sample", "Hugo_Symbol"))
  
  target.table<-merge(thousand.amaf, tcga.maf, by=c("Position","Chrom"))
  
  #Calculate the conditional probability for cancer - Probability of having variation given that you have cancer, which is position maf given by MUT.FREQ
  #Calculate the conditional probability for non-cacner - Probability of having variation given that you don't have cancer, which is AMAF given by MAF in thousand 
  
  #Apply bayesian method per position on each patient to obtain posterior probabilities
  target.table$PROB<-(target.table$MUT.FREQ*cancer.prob)/(target.table$MUT.FREQ*cancer.prob  + target.table$MAF*(1-cancer.prob))
  
  #Clean up and return
  target.table<-target.table[order(PROB, decreasing=T),]
  return(target.table)
}

Function.Pairwise.Decay.Exon<-function(thousand.table, exons.table, gene, lambda, phastcons.table=NULL, a.set, filter=F){
  #Calculates decay based on exon information only
  #   NOTE: May need to change in the future to depend on the amount of information (1000G) we have per exon (MAF count per exon length)
  #Implemented weighted means for opposite decaying functions
  #Implemented multiple maf factors (need to remove from final version
  
  require(data.table)
  require(IRanges)
  
  #####CREATE TARGET MATRIX BASED ON NON-OVERLAPPING EXONS######
  exons.table<-exons.table[Hugo_Symbol==gene,]
  start<-min(unique(exons.table$FEAT_START))
  end<-max(unique(exons.table$FEAT_END))
  chrom<-unique(exons.table$Chrom)
  
  #Reduce to non-overlapping exons
  exon.ranges<-IRanges(exons.table$FEAT_START, exons.table$FEAT_END)
  exon.ranges<-reduce(exon.ranges)
  exon.ranges<-data.table(exon.start=start(exon.ranges), exon.end=end(exon.ranges))
  
  #Filter gene matrix for exon MAFs only
  pre.target.matrix<-thousand.table[Chrom==chrom & Position>start & Position<end,]
  
  #Filter for low count MAFs?
  if (filter==T){
    pre.target.matrix<-pre.target.matrix[MAF>=0.000199682,]
  }
  
  exon.mafs<-sapply(as.vector(pre.target.matrix$Position), function(x) {
    any(apply(exon.ranges, 1, function(y) x>=y[1] & x<=y[2]))
  })

  target.matrix<-pre.target.matrix[exon.mafs,]
  
  #Check if we have data in 1000G data, if not, stop function - (12.14.14)
  if (nrow(target.matrix)<1){
    return (paste(gene, "not found in 1000G data"))
  }
  
  target.matrix<-target.matrix[order(Position),]
  
  #Clean target matrix for position duplicates and aggregate(sum) MAFs contained in duplicates
  target.matrix<-aggregate(MAF~Position,data=target.matrix,FUN=sum)
  print ("Done cleaning up target matrix")
  
  #Estimate data coverage - TEMPORAL!!!! (12.13.14)
  #coverage<-Function.Maf.Exon.Coverage(exons.table, target.matrix, gene=gene)
  #print (coverage)
  
  #Get minimum standard maf for region based on available data (has to be done before modeling)
  #min.maf<-min(target.matrix$MAF)
  #print (min.maf)
  
  ######PHASTCONS###### - DYNAMIC
  if (!is.null(phastcons.table)){
    print ("Calculating")
    
    phastcons.table<-phastcons.table[Chrom==as.character(chrom) & Position>start & Position<end,]
    print ("Done processing phastcons")
    target.matrix<-merge(target.matrix, phastcons.table, by="Position")
    print ("Done integrating target and phastcons")
    
    #Prep storage vectors
    POS<-as.vector(target.matrix$Position)
    MAF<-as.vector(target.matrix$MAF)
    PHAST<-as.vector(target.matrix$Score)
    EXPANSION.MEAN<-c()
    
    #Loop through a.factors -REMOVE, DONE FOR TESTING MULTIPLE a.factors only
    a.list<-list()
    for (a.factor in a.set) {
      #Construct left decay first (before first maf data)
      decay.vector<-rev(Function.Exp.Decay(MAF[1], POS[1]-start, fixed.lambda=F, phastcon=PHAST[1], a=a.factor))
      
      #Check that we actually have more than 1 position for middle loop
      if (length(POS)>1){
        #Then middle
        for (n in 1:(length(POS)-1)){
          
          #If they are separated by 1 or not, either way store maf value
          decay.vector<-c(decay.vector, MAF[n])
          
          if (POS[n]!=(POS[n+1]-1)){
            
            #If we are not separated by 1, then apply extension - IMPLEMENTED WEIGHTED MEANS (121014)
            decay.expansion<-POS[n+1]-POS[n]-1
            EXPANSION.MEAN<-c(EXPANSION.MEAN, decay.expansion)
            
            right.decay<-Function.Exp.Decay(MAF[n], decay.expansion, fixed.lambda=F, phastcon=PHAST[n], a=a.factor)
            weighted.right.decay<-right.decay*(decay.expansion:1)/decay.expansion
            
            left.decay<-rev(Function.Exp.Decay(MAF[n+1], decay.expansion, fixed.lambda=F, phastcon=PHAST[n], a=a.factor))
            weighted.left.decay<-left.decay*(1:decay.expansion)/decay.expansion
            
            decay.vector<-c(decay.vector, colSums(rbind(weighted.right.decay, weighted.left.decay)))
          } 
        }
      } 
      
      #Then right tail
      decay.vector<-c(decay.vector, tail(MAF,1))
      last.vector<-Function.Exp.Decay(tail(MAF,1), end-tail(POS,1), fixed.lambda=F, phastcon=tail(PHAST,1), a=a.factor)
      decay.vector<-c(decay.vector,last.vector)
      
      #Then append to a.list
      a.list<-c(a.list,list(decay.vector))
    }
    
    names(a.list)<-as.character(a.set)
    
    #POST-PROCESS
    #Loop through a.factors to assign separate lists
    multiple.decay.table<-list()
    for (a.factor in a.set){
      decay.vector<-a.list[[as.character(a.factor)]]
      
      #Replace all values in decay function that are equal to zero by the minimum value in the decay vector - FIXES BRCA1 problem of low comparisson
      decay.vector[decay.vector==0]<-min(decay.vector[decay.vector>0])
      
      #Match up to actual positions
      position.vector<-start:end
      decay.table<-data.table(Chrom=chrom, Position=position.vector, MAF=decay.vector)
      
      #Append to decay table
      decay.table$a.factor<-a.factor
      multiple.decay.table[[as.character(a.factor)]]<-decay.table
    }
    
    multiple.decay.table<-do.call(rbind,multiple.decay.table)
    multiple.decay.table<-as.data.table(multiple.decay.table)
    
    #Add median expansion as reference
    multiple.decay.table$expansion.median<-median(EXPANSION.MEAN)
    
    #Return
    print (median(EXPANSION.MEAN))
    return(multiple.decay.table)
    
    ########### NON-PHASTCON ######## - STATIC
  } else {
    print ("Calculating - NON - DYNAMIC")
    
    #Prep storage vectors
    POS<-as.vector(target.matrix$Position)
    MAF<-as.vector(target.matrix$MAF)
    EXPANSION.MEAN<-c()
    
    #Construct left decay first (before first maf data)
    decay.vector<-rev(Function.Exp.Decay(MAF[1], POS[1]-start, lambda=lambda))
    
    #Check that we actually have data for more than one position before doing loop
    if (length(POS)>1){
      #Then middle
      for (n in 1:(length(POS)-1)){
        
        #If they are separated by 1 or not, either way store maf value
        decay.vector<-c(decay.vector, MAF[n])
        
        if (POS[n]!=(POS[n+1]-1)){
          
          #If we are not separated by 1, then apply extension
          decay.expansion<-POS[n+1]-POS[n]-1
          EXPANSION.MEAN<-c(EXPANSION.MEAN, decay.expansion)
          
          right.decay<-Function.Exp.Decay(MAF[n], decay.expansion, lambda=lambda)
          weighted.right.decay<-right.decay*(decay.expansion:1)/decay.expansion
          
          left.decay<-rev(Function.Exp.Decay(MAF[n+1], decay.expansion, lambda=lambda))
          weighted.left.decay<-left.decay*(1:decay.expansion)/decay.expansion
          
          decay.vector<-c(decay.vector, colSums(rbind(weighted.right.decay, weighted.left.decay)))
        } 
      }
    } 
    
    #Then right tail
    decay.vector<-c(decay.vector, tail(MAF,1))
    last.vector<-Function.Exp.Decay(tail(MAF,1), end-tail(POS,1), lambda=lambda)
    decay.vector<-c(decay.vector,last.vector)
    
    #POST - PROCESS
    #Replace all values in decay function that are equal to zero by the minimum value in the decay vector - FIXES BRCA1 problem of low comparisson
    decay.vector[decay.vector==0]<-min(decay.vector[decay.vector>0])
    
    #Match up to actual positions
    position.vector<-start:end
    decay.table<-data.table(Chrom=chrom, Position=position.vector, MAF=decay.vector)
    
    #Add median expansion as reference
    decay.table$expansion.median<-median(EXPANSION.MEAN)
    
    #Return
    print (median(EXPANSION.MEAN))
    return(decay.table) 
  }
}

#Load phastcon scores - USED PHAST45 on 121714
PHAST<-fread("DATABASES/PHASTCONS/120314.CHRM.FILTERED.EXONS.45.PARALLEL", header=F, sep="\t", stringsAsFactors=F)
setnames(PHAST, c("Chrom", "Position", "Score"))

#Test gene
BRCA.COSMIC.GENES<-c("AKT1", "BAP1", "BRCA1", "BRCA2", "BRIP1", "CDH1", "EP300", "ERBB2",
                     "FOXA1", "MAP2K4","PALB2","PBRM1", "PIK3CA", "RB1","TP53")

gene.testing.plot<-function(gene, a.factor, static=F, lambda.factors=c(), EXON=F, filter=F) {
  #Testing probabilities across different parameters.
  #static=F implies that we are not using phastcon scores to predict decay rate and therefore pre-determined decay rates are used
  
  target.matrix<-data.frame()
  
  if (EXON==F){
    if(static==F){
      gene.decay<-Function.Pairwise.Decay(THOUSAND, EXONS, gene, phastcons.table=PHAST, a.set=a.factor, filter=filter)
      
      decay.split<-split(gene.decay, gene.decay$a.factor)
      
      for (decay in decay.split){
        gene.bayes<-Function.maf.bayes(decay[,!c("a.factor"),with=F], TCGA.BRCA.maf, gene)
        gene.bayes$a.factor<-unique(as.vector(decay$a.factor))
        target.matrix<-rbind(target.matrix, gene.bayes)   
      }
    } else { 
      for (lambda in lambda.factors){
        gene.decay<-Function.Pairwise.Decay(THOUSAND, EXONS, gene, lambda=lambda, filter=filter)
        gene.bayes<-Function.maf.bayes(gene.decay, TCGA.BRCA.maf, gene)
        gene.bayes$lambda<-lambda
        target.matrix<-rbind(target.matrix, gene.bayes)
      }
    }
  } else {
    if(static==F){
      gene.decay<-Function.Pairwise.Decay.Exon(THOUSAND, EXONS, gene, phastcons.table=PHAST, a.set=a.factor, filter=filter)
      
      decay.split<-split(gene.decay, gene.decay$a.factor)
      
      for (decay in decay.split){
        gene.bayes<-Function.maf.bayes(decay[,!c("a.factor"),with=F], TCGA.BRCA.maf, gene)
        gene.bayes$a.factor<-unique(as.vector(decay$a.factor))
        target.matrix<-rbind(target.matrix, gene.bayes)   
      }
    } else { 
      for (lambda in lambda.factors){
        gene.decay<-Function.Pairwise.Decay.Exon(THOUSAND, EXONS, gene, lambda=lambda,filter=filter)
        gene.bayes<-Function.maf.bayes(gene.decay, TCGA.BRCA.maf, gene)
        gene.bayes$lambda<-lambda
        target.matrix<-rbind(target.matrix, gene.bayes)
      }
    }
  }
  
  #Return
  target.matrix<-as.data.table(target.matrix)
  return(target.matrix)
}

#Test with PHASTCONS - DYNAMICCS - NEED TO RE-RUN TO OBTAIN GRAPHICS DUE TO BUG
#COSMIC
for (cancer in BRCA.COSMIC.GENES){
  print (cancer)
  gene.bayes<-gene.testing.plot(cancer, a.factor=seq(1,9,2))
  assign(paste0(cancer,".BAYES"), gene.bayes)
}

TP53.BAYES<-gene.testing.plot("TP53", a.factor=seq(1,9,2))
TTN.BAYES<-gene.testing.plot("TTN", a.factor=seq(1,9,2)) #NON.CANCER
ggplot(TTN.BAYES, aes(PROB)) +geom_histogram() + theme.format + facet_grid(~a.factor)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 16))

count=1
for (gene in intersect(sample(unique(as.vector(EXONS$Hugo_Symbol)),10),unique(as.vector(TCGA.BRCA.maf$Hugo_Symbol)))) {
  if (!(paste0(as.character(gene),".BAYES") %in% BAYES.SET)){
    print (gene)
    gene.test<-gene.testing.plot(gene, a.factor=seq(1,9,2))
    assign(paste0(gene,".BAYES"), gene.test, envir = .GlobalEnv) 
    print (count)
    count=count+1
  }
}

BAYES.SET<-setdiff(setdiff(setdiff(ls(pattern=".BAYES"), ls(pattern=".BAYES.STATIC")), ls(pattern=".BAYES.PHAST")), ls(pattern="BAYES.EXON"))
sum(sapply(BAYES.SET, function(x) strsplit(x,".BA")[[1]][1]) %in% BRCA.COSMIC.GENES)

#Visualize probabilities in phastcons samples
BAYES.DISTRIBUTION<-data.frame()
for (gene.bayes in BAYES.SET) {
  bayes.table<-get(gene.bayes)
  bayes.table<-bayes.table[,c("Hugo_Symbol", "PROB", "a.factor"), with=F]
  bayes.table$CANCER<-bayes.table$Hugo_Symbol %in% BRCA.COSMIC.GENES
  BAYES.DISTRIBUTION<-rbind(BAYES.DISTRIBUTION, bayes.table)
}
ggplot(BAYES.DISTRIBUTION, aes(factor(a.factor), PROB, colour=CANCER)) + geom_boxplot() + geom_jitter() + theme.format +
  facet_grid(~CANCER) + theme(strip.text.x = element_text(size = 16))
length(unique(BAYES.DISTRIBUTION[PROB>0.9 & a.factor>=9,]$Hugo_Symbol))

######TEST with STATIC lambda - No Phastcons
TP53.BAYES.STATIC<-gene.testing.plot("TP53", static=T, lambda.factors=seq(0.01,0.09,0.02))
TTN.BAYES.STATIC<-gene.testing.plot("TTN", static=T, lambda.factors=seq(0.01,0.09,0.02))
TTN.BAYES.STATIC[lambda==0.03 & PROB>=0.9,]

ggplot(TTN.BAYES.STATIC, aes(PROB)) +geom_histogram() + theme.format + facet_grid(~lambda)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 16))

for (gene in setdiff(BRCA.COSMIC.GENES, sapply(BAYES.SET.STATIC, function(x) strsplit(x,".BAYES")[[1]][1])) ){
  print (gene)
  gene.test<-gene.testing.plot(gene, static=T, lambda.factors=seq(0.01,0.09,0.02))
  assign(paste0(gene,".BAYES.STATIC"), gene.test, envir = .GlobalEnv) 
}

#Visualize probabilities in NON-DYNAMIC samples
BAYES.SET.STATIC<-setdiff(ls(pattern=".BAYES.STATIC"), ls(pattern=".BAYES.STATIC.EXON"))
BAYES.DISTRIBUTION.STATIC<-data.frame()
for (gene.bayes in BAYES.SET.STATIC) {
  bayes.table<-get(gene.bayes)
  bayes.table<-bayes.table[,c("Hugo_Symbol", "PROB", "lambda"), with=F]
  bayes.table$CANCER<-bayes.table$Hugo_Symbol %in% BRCA.COSMIC.GENES
  BAYES.DISTRIBUTION.STATIC<-rbind(BAYES.DISTRIBUTION.STATIC, bayes.table)
}
ggplot(BAYES.DISTRIBUTION.STATIC, aes(factor(lambda), PROB, colour=CANCER)) + geom_boxplot() + geom_jitter() + theme.format +
  facet_grid(~CANCER) + theme(strip.text.x = element_text(size = 16))
BAYES.DISTRIBUTION.STATIC[lambda==0.03 & PROB>0.9,]

BAYES.STATIC.CHART<-BAYES.DISTRIBUTION.STATIC[,internal.1(Hugo_Symbol, PROB),by=c("lambda","CANCER")]
BAYES.STATIC.CHART<-melt(BAYES.STATIC.CHART, id.vars=c("lambda","CANCER"))
setnames(BAYES.STATIC.CHART, c("lambda","CANCER","PROB","COUNT"))
ggplot(BAYES.STATIC.CHART, aes(lambda,COUNT, fill=CANCER)) + geom_histogram(stat="identity",position="dodge") + theme.format +
  facet_wrap(~PROB)+ theme(strip.text.x = element_text(size = 16))

######TEST with STATIC-EXON lambda - No Phastcons
TP53.BAYES.STATIC.EXON<-gene.testing.plot("TP53", static=T, lambda.factors=seq(0.005, 0.03,0.005),EXON=T)
TTN.BAYES.STATIC.EXON<-gene.testing.plot("TTN", static=T, lambda.factors=seq(0.005, 0.03,0.005),EXON=T, filter=T)

TTN.BAYES.STATIC.EXON[lambda==0.01 & PROB>=0.9,]
TCGA.BRCA[Hugo_Symbol=="TTN" & Start_Position %in% as.vector(TTN.BAYES.STATIC.EXON[lambda==0.01 & PROB>=0.9,]$Position),]

ggplot(TTN.BAYES.STATIC.EXON, aes(PROB)) +geom_histogram() + theme.format + facet_grid(~lambda)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 16))

count=1
for (gene in intersect(sample(unique(as.vector(EXONS$Hugo_Symbol)),100),unique(as.vector(TCGA.BRCA.maf$Hugo_Symbol)))) {
  if (!(paste0(as.character(gene),".BAYES.STATIC.EXON") %in% ls(pattern=".BAYES.STATIC.EXON") )){
    print (gene)
    gene.test<-gene.testing.plot(gene, static=T, lambda.factors=seq(0.005, 0.03,0.005), EXON=T)
    assign(paste0(gene,".BAYES.STATIC.EXON"), gene.test, envir = .GlobalEnv) 
    print (count)
    count=count+1
  }
}

#Visualize probabilities in NON-DYNAMIC samples
BAYES.SET.STATIC.EXON<-ls(pattern=".BAYES.STATIC.EXON")
BAYES.DISTRIBUTION.STATIC.EXON<-data.frame()
for (gene.bayes in BAYES.SET.STATIC.EXON) {
  bayes.table<-get(gene.bayes)
  bayes.table<-bayes.table[,c("Hugo_Symbol", "PROB", "lambda", "expansion.median"), with=F][lambda %in% seq(0.005, 0.03,0.005), ]
  bayes.table$CANCER<-bayes.table$Hugo_Symbol %in% BRCA.COSMIC.GENES
  BAYES.DISTRIBUTION.STATIC.EXON<-rbind(BAYES.DISTRIBUTION.STATIC.EXON, bayes.table)
}
ggplot(BAYES.DISTRIBUTION.STATIC.EXON, aes(expansion.median, PROB, colour=lambda)) + geom_point() + theme.format + facet_grid(~lambda) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(strip.text.x = element_text(size = 16)) + scale_x_log10()

UBUNTU.STATIC.EXON<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/121714.BREAST.STATIC.EXON.800.SAMPLES.DISTRIBUTION.rds")
BAYES.DISTRIBUTION.STATIC.EXON<-unique(rbind(BAYES.DISTRIBUTION.STATIC.EXON, UBUNTU.STATIC.EXON))

ggplot(BAYES.DISTRIBUTION.STATIC.EXON, aes(factor(lambda), PROB, colour=CANCER)) + geom_boxplot() + geom_jitter(size=1) + theme.format +
  facet_grid(~CANCER) + theme(strip.text.x = element_text(size = 16))

internal.1<-function(genes, probs){
  test.table<-data.table(GENE=genes, PROB=probs)
  return(list(P.70=length(unique(test.table[PROB>=0.7,]$GENE)),
              P.80=length(unique(test.table[PROB>=0.8,]$GENE)),
              P.90=length(unique(test.table[PROB>=0.9,]$GENE)),
              P.100=length(unique(test.table[PROB>=1,]$GENE))))
}
BAYES.STATIC.EXON.CHART<-BAYES.DISTRIBUTION.STATIC.EXON[,internal.1(Hugo_Symbol, PROB),by=c("lambda","CANCER")]
BAYES.STATIC.EXON.CHART<-melt(BAYES.STATIC.EXON.CHART, id.vars=c("lambda","CANCER"))
setnames(BAYES.STATIC.EXON.CHART, c("lambda","CANCER","PROB","COUNT"))
ggplot(BAYES.STATIC.EXON.CHART, aes(lambda,COUNT, fill=CANCER)) + geom_histogram(stat="identity",position="dodge") + theme.format +
  facet_wrap(~PROB)+ theme(strip.text.x = element_text(size = 16))

#####TEST with DYNAMIC-EXON - PHASTCONS
TP53.BAYES.EXON<-gene.testing.plot("TP53", a.factor=seq(1,9,2), EXON=T)
TTN.BAYES.EXON<-gene.testing.plot("TTN", a.factor=seq(1,9,2), EXON=T)
BRCA1.BAYES.EXON<-gene.testing.plot("BRCA1", a.factor=seq(1,9,2), EXON=T)
ZZZ3.BAYES.EXON<-gene.testing.plot("ZZZ3", a.factor=seq(1,9,2), EXON=T)
PIK3CA.BAYES.EXON<-gene.testing.plot("PIK3CA", a.factor=seq(1,9,2), EXON=T)
CCND1.BAYES.EXON<-gene.testing.plot("CCND1", a.factor=seq(1,9,2), EXON=T)

ggplot(TTN.BAYES.EXON, aes(PROB)) +geom_histogram() + theme.format + facet_grid(~a.factor)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 16))

#######TEST with DYNAMICCS -PHASTCONS 45
#COSMIC
for (cancer in BRCA.COSMIC.GENES){
  print (cancer)
  gene.bayes<-gene.testing.plot(cancer, a.factor=seq(1,9,2))
  assign(paste0(cancer,".BAYES"), gene.bayes)
}

TP53.BAYES.DYNAMIC.45<-gene.testing.plot("TP53", a.factor=seq(1,9,2))
TTN.BAYES.DYNAMIC.45<-gene.testing.plot("TTN", a.factor=seq(1,9,2)) #NON.CANCER
ggplot(TTN.BAYES.DYNAMIC.45, aes(PROB)) +geom_histogram() + theme.format + facet_grid(~a.factor)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.x = element_text(size = 16))


################################################################################################
#Test with phastcons-dependent lambda
TP53.DECAY.PHAST<-Function.Pairwise.Decay(THOUSAND, EXONS, "TP53", phastcons.table=PHAST)
TP53.BAYES.PHAST<-Function.maf.bayes(TP53.DECAY.PHAST, TCGA.BRCA.maf, "TP53")

for (gene in c("TP53", "TTN", "PIK3CA", "BRCA1","CDH1")){
  name<-paste0(gene,".DECAY.PHAST.a.1")
  name.function<-Function.Pairwise.Decay(THOUSAND, EXONS, gene, phastcons.table=PHAST, a=1)
  assign(name, name.function, envir = .GlobalEnv)
}
for (gene in c("TP53", "TTN", "PIK3CA", "BRCA1", "CDH1")){
  name<-paste0(gene,".BAYES.PHAST.a.1")
  name.function<-Function.maf.bayes(get(paste0(gene,".DECAY.PHAST.a.1")),TCGA.BRCA.maf, gene )
  assign(name, name.function, envir = .GlobalEnv)
}
hist(TTN.BAYES.PHAST.a.1$PROB)

for (gene in c("TP53", "TTN", "PIK3CA", "BRCA1", "BRCA2", "MUC4","CDH1")){
  filename<-paste0("FIGURES/121114.",gene, ".PHAST.100.DECAY.a.2.jpeg")
  jpeg(filename,width=1080, height=1080, quality=100, res=100)
  object<-get(paste0(gene, ".BAYES.PHAST.a.2"))
  hist(object[["PROB"]])
  dev.off()
}

  
hist(PIK3CA.BAYES.PHAST.a.2$PROB)

system.time(Function.Exp.Decay(0.5,30000,20000,fixed.lambda=F,phastcon=1))

#Potential decay functions
plot(1:150, Function.Exp.Decay(0.5,150,lambda=0.00))

x=seq(0,1,0.05)
plot(x, 1/(100^x))
plot(x, 0.005/(0.01+x))
plot(x, 1-(1/(1+x^2)))
plot(x, x^9/(1+x^9))