#Function.u.R
#071714
#Compares the background mutation rate and the metabolic mutation rate per metabolite across patients
#It calculates significant mutations associated with metabolites using 4 tests: student-t-test, paired-t-test, mann-whitney-u and wilcoxon signed-rank test
#It returns a list of:
#   STATS - Statistics calculates
#   STARTING.PATIENTS - Number of patients fed to the script
#   PATIENTS.IN.MODEL - Number of patients actually used by the script
#   PATIENTS.AFFECTED - List per test of number of patients that have mutations in genes calculated as significant by the test, after fdr

############################################################################################################################################################

Function.all.t.tests<-function(vector.a, vector.b, ALT="two.sided", MU=0){
  #Understand CAREFULLY what you are inputting as MU!!
  #RETURNS LIST!!
  #Default for alternative hypothesis is "two.sided"
  
  #student's t.test
  dummy.t.test<-t.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU)
  
  #pair t.test
  dummy.paired.t.test<-t.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU, paired=T)
  
  #Mann-whitney U-test
  dummy.mann.whitney<-wilcox.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU)
  
  #Wilcoxon-signed rank test
  dummy.wilcoxon<-wilcox.test(vector.a, vector.b, var.equal=T, alternative=ALT, mu=MU, paired=T)
  
  #Return results of all tests
  return(list(student.t.stat=dummy.t.test$statistic, student.p.val=dummy.t.test$p.value,
              paired.t.stat=dummy.paired.t.test$statistic, paired.t.p.val=dummy.paired.t.test$p.value,
              mann.whintey.u.stat=dummy.mann.whitney$statistic, mann.whitney.u.p.val=dummy.mann.whitney$p.value,
              wilcoxon.stat=dummy.wilcoxon$statistic, wilcoxon.p.val=dummy.wilcoxon$p.value))
}

Function.Outliers.Table.1<-function(table.1) {
  
  #Get mutations counts per patient
  table.1.count<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)),by = "Tumor_Sample_Barcode"]
  
  #Calculate non-outlier range
  summary.stats<-summary(table.1.count$TOTAL.MUTATIONS)
  IR.RANGE<-summary.stats[[5]]-summary.stats[[2]]
  LOW.1.5<-summary.stats[[2]]-1.5*IR.RANGE
  HIGH.1.5<-summary.stats[[5]]+1.5*IR.RANGE
  
  #Remove outliers
  table.1.count<-table.1.count[TOTAL.MUTATIONS>=LOW.1.5,]
  table.1.count<-table.1.count[TOTAL.MUTATIONS<=HIGH.1.5,]
  table.1<-table.1[Tumor_Sample_Barcode %in% table.1.count$Tumor_Sample_Barcode,]
  
  #Return filtered table.1
  return(table.1)
}

Function.BMR<-function(main.table, background.length) {
    #Calculates mutation significance of each gene based on the background mutation rate per patient
    #This is done for missense mutations only
    
    require(data.table)
    
    #Get total mutations per patients [PATIENT, BM]
    BM.table<-main.table[,list(BM=sum(Missense)), by="PATIENT"]
    
    #Merge to main table
    main.table<-as.data.table(merge(as.data.frame(main.table),
    as.data.frame(BM.table), by="PATIENT"))
    
    #Calculate significance
    sig.table<-main.table[,list(P.VAL=phyper(q=Missense-1, m=Length,
    n=background.length-Length, k=BM, lower.tail=F)),
    by=c("PATIENT", "Hugo_Symbol", "Missense", "Silent", "Length")]
    sig.table$P.VAL.ADJ<-p.adjust(sig.table$P.VAL, method="fdr")
    
    #Clean up and Return
    sig.table<-sig.table[order(P.VAL.ADJ),]
    return(sig.table)
}

Function.pre.process<-function(table.1, background.sequenced.genes, length.file, bmr=TRUE){
    #Pre process table.1 to give appropriate protein length and obtain total background amino acid length
    #If bmr=TRUE, will apply background mutation function before wilcoxon analysis
    
    require(data.table)
    
    #Use uniprot gene names and synonyms to map lengths for table.1 [Hugo_Symbol, Length]
    uniprot.table<-as.data.table(read.csv(length.file, header=T, sep="\t", stringsAsFactors=F))
    uniprot.table<-uniprot.table[,c("Gene.names", "Length"),with=F]
    uniprot.table$ID<-1:nrow(uniprot.table)
    uniprot.table<-uniprot.table[,list(Hugo_Symbol= sub(";", "",strsplit(Gene.names," ")[[1]]), Length=Length ), by="ID"]
    uniprot.table$ID<-NULL
    
    #Filter for duplicated gene records based on largest gene size
    uniprot.table<-uniprot.table[order(Length, decreasing=T),]
    uniprot.table<-uniprot.table[!duplicated(Hugo_Symbol),]
    
    #Merge with table.1 to get length information
    #[Hugo_Symbol, Missense, Silent, PATIENT, Length]
    #main.table<-table.1.pval[P.VAL.ADJ<0.05,][,c("Hugo_Symbol", "Missense","Silent", "PATIENT"), with=F]
    #main.table<-as.data.table(merge(as.data.frame(main.table), as.data.frame(uniprot.table), by="Hugo_Symbol"))
    main.table<-as.data.table(merge(as.data.frame(table.1),as.data.frame(uniprot.table), by="Hugo_Symbol"))
    
    #Get total hypothetical background mutation length
    background.length<-uniprot.table[Hugo_Symbol %in% as.vector(background.sequenced.genes),]
    print (background.length)
    background.length<-sum(as.vector(background.length$Length))
    print (background.length)
    #background.length<-9000000
    #background.length<-10000000
    #background.length<-11317118
    
    #Are we pre-filtering by bmr.function()?
    if (bmr==TRUE){
        main.table<-Function.BMR(main.table, background.length)
        main.table<-main.table[,c("Hugo_Symbol", "Missense","Silent", "PATIENT","Length","P.VAL.ADJ"), with=F]
    }
    
    #Return
    return(main.table)
}

Function.u.t.bmr<-function(table.1, table.2, length.table, background.sequenced.genes, remove.outliers=F, bmr.test=TRUE) {
  #Calculate u(j) weights per metabolite j 
  #This compares the bmr per patient of proteins associated with j versus a constant sequenced proteome bmr
  #This function assumes that all genes in table.2 are covered in length.table
  
  require(data.table)
  require(reshape)
  
  ####QUICK FIX FOR NEW table.1 version#####101514
  #This includes optional corrections for bmr using hypergeometric test
  if (bmr.test==TRUE){
      table.1<-Function.pre.process(table.1, background.sequenced.genes)
  } else {
      table.1<-table.1[Missense!=0,]
      table.1$Tumor_Sample_Barcode<-table.1$PATIENT
      table.1$N.MUTATIONS<-table.1$Missense
  }
  
  
  
  ####PREP####
  
  #Get total count of patients we are starting with
  TOTAL.PATIENT.COUNT<-length(unique(as.vector(table.1$Tumor_Sample_Barcode)))
  
  #Remove outliers if necessary
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }
  
  #Rename tables
  setnames(length.table, colnames(length.table), c("Hugo_Symbol", "Length"))
  length.table<-as.data.table(length.table)
  
  #Remove duplicates from table.2 [KEGG_ID, Hugo_Symbol]
  table.2<-unique(table.2[,c(2,3),with=F])
  setnames(table.2, colnames(table.2), c("KEGG_ID", "Hugo_Symbol"))
  
  #Get total background amino acid length (those that could be measured) - USED FOR TOTAL.BMR
  background.length<-length.table[Hugo_Symbol %in% background.sequenced.genes,]
  background.length<-sum(as.vector(background.length$Length)) #Theoretical length of all translatable amino acids added up
  
  #Get a per metabolite length table [KEGG_ID, METABOLIC.LENGTH]
  metabolite.length<-merge(table.2, length.table, by="Hugo_Symbol")
  metabolite.length<-metabolite.length[,list(METABOLIC.LENGTH=sum(Length)), by="KEGG_ID"]
  
  #Get total count of mutations per patient (meatabolic and non-metabolic) [Tumor_Sample_Barcode, TOTAL.MUTATIONS]
  patient.all.mutations<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)), by="Tumor_Sample_Barcode"]
  
  ####PROCESS####
  
  #Construct matrix of patients(table.1) x metabolites(table.2) with values as total affected metabolic mutations per patient per metabolite
  patients.x.metabolites<-merge(as.data.frame(table.1), as.data.frame(table.2), by="Hugo_Symbol")
  patients.x.metabolites<-patients.x.metabolites[,2:4] #Don't need gene column anymore
  patients.x.metabolites<-cast(patients.x.metabolites, Tumor_Sample_Barcode~KEGG_ID, value="N.MUTATIONS", sum, fill=0) #MATRIX
  
  #Prep data.table for per meatbolite across patient statistic calculations
  patients.t.metabolites<-melt(patients.x.metabolites, id="Tumor_Sample_Barcode") #[Tumor_Sample_Barcode, value, KEGG_ID]
  setnames(patients.t.metabolites, colnames(patients.t.metabolites), c("Tumor_Sample_Barcode", "METABOLIC.MUTATIONS", "KEGG_ID"))
  patients.t.metabolites<-merge(patients.t.metabolites, as.data.frame(patient.all.mutations), by="Tumor_Sample_Barcode") #To get background-affected per patient
  
  #[KEGG_ID, Tumor_Sample_Barcode, METABOLIC.MUTATIONS, KEGG_ID, TOTAL.MUTATIONS, METABOLIC.LENGTH]
  patients.t.metabolites<-merge(patients.t.metabolites, as.data.frame(metabolite.length), by="KEGG_ID")
  patients.t.metabolites<-as.data.table(patients.t.metabolites)
  
  #Get patients that are used in the model (number of patients that we actually use and contain metabolic information)
  PATIENT.FOR.MODEL.COUNT<-length(unique(as.vector(patients.t.metabolites$Tumor_Sample_Barcode)))
  
  #Obtain rates
  patients.t.metabolites$METABOLIC.RATE<-patients.t.metabolites$METABOLIC.MUTATIONS/patients.t.metabolites$METABOLIC.LENGTH
  patients.t.metabolites$PATIENT.BMR<-patients.t.metabolites$TOTAL.MUTATIONS/background.length
  
  #STATS - Calculate student's t.test, paired t.test, mann whitney u and wilcoxon-signed ranked test per metabolite across patients
  KEGG.STATS<-patients.t.metabolites[, Function.all.t.tests(METABOLIC.RATE, PATIENT.BMR, ALT="greater", MU=0) , by="KEGG_ID"]
  
  #Correct for multiple hypothesis testing
  KEGG.STATS$STUDENT.P.ADJ<-p.adjust(KEGG.STATS$student.p.val, method="fdr")
  KEGG.STATS$PAIRED.P.ADJ<-p.adjust(KEGG.STATS$paired.t.p.val, method="fdr")
  KEGG.STATS$MANN.U.P.ADJ<-p.adjust(KEGG.STATS$mann.whitney.u.p.val, method="fdr")
  KEGG.STATS$WILCOXON.P.ADJ<-p.adjust(KEGG.STATS$wilcoxon.p.val, method="fdr")
  
  #Order by WILCOXON.SIGNED.RANK.TEST (NON-PARAMETRIC PAIRED T.TEST) (fdr corrected)
  KEGG.STATS<-KEGG.STATS[order(WILCOXON.P.ADJ),]
  
  #Get affected patients (patients that have mutations related to enriched metabolites)
  STUDENT.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[KEGG_ID %in% KEGG.STATS[STUDENT.P.ADJ<0.05,]$KEGG_ID ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  PAIRED.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[KEGG_ID %in% KEGG.STATS[PAIRED.P.ADJ<0.05,]$KEGG_ID ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  MANN.U.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[KEGG_ID %in% KEGG.STATS[MANN.U.P.ADJ<0.05,]$KEGG_ID ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  WILCOXON.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[KEGG_ID %in% KEGG.STATS[WILCOXON.P.ADJ<0.05,]$KEGG_ID ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  
  PATIENT.AFFECTED<-list(STUDENT.AFFECTED=STUDENT.AFFECTED,
                         PAIRED.AFFECTED=PAIRED.AFFECTED,
                         MANN.U.AFFECTED=MANN.U.AFFECTED,
                         WILCOXON.AFFECTED=WILCOXON.AFFECTED)
  
  #Return list including STATS, count of initial patients, patients used in model and patients per test affected (have mutations related to enriched j)
  dummy.return=list(STATS=KEGG.STATS, STARTING.PATIENTS=TOTAL.PATIENT.COUNT, PATIENTS.IN.MODEL=PATIENT.FOR.MODEL.COUNT, PATIENT.AFFECTED=PATIENT.AFFECTED)
  return(dummy.return)
  
}

#PROCESS ENTRIES
require(base)
require(data.table)

args<-commandArgs(trailingOnly=T)

input.table.1.rds<-args[1]
input.table.2.rds<-args[2]
length.file<-args[3] #File that contains a table of 2 columns (Hugo_Symbo, Length in amino acids)
output.file<-args[4]

table.1.input<-readRDS(input.table.1.rds)
table.2.input<-readRDS(input.table.2.rds)
length.file.input<-as.data.table(read.csv(length.file, header=T, sep="\t"))

#EXECUTE
Run<-Function.u.t.bmr(table.1.input$table.1, table.2.input, length.file.input, table.1.input$background.genes, remove.outliers=F)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)


