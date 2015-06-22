#Function.u.p.R
#072014
#Compares the background mutation rate and the per gene mutation rate  across patients
#It calculates significant mutations associated with a protein using 4 tests: student-t-test, paired-t-test, mann-whitney-u and wilcoxon signed-rank test
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

Function.u.t.bmr.p<-function(table.1, length.table, background.sequenced.genes, remove.outliers=F) {
  #Calculate u(j) weights per gene j 
  #This compares the bmr per patient of protein p versus a constant sequenced proteome bmr
  #This function will filter all proposed genes for only those where there is length information available in length.table
  
  require(reshape2)
  require(data.table)
  
  ####PREP####
  
  #Get total count of patients we are starting with
  TOTAL.PATIENT.COUNT<-length(unique(as.vector(table.1$Tumor_Sample_Barcode)))
  
  #Remove outliers if necessary
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }
  
  #Rename tables
  setnames(length.table, colnames(length.table), c("Hugo_Symbol", "Length"))
  length.table<-as.data.table(length.table) #[Hugo_Symbol, Length]
  
  #Filter table.1 for genes that contain length information in length.table
  table.1<-table.1[Hugo_Symbol %in% unique(as.vector(length.table$Hugo_Symbol)), ]
  
  #Get total background amino acid length (those that could be measured) - USED FOR TOTAL.BMR
  background.length<-length.table[Hugo_Symbol %in% background.sequenced.genes,]
  background.length<-sum(as.vector(background.length$Length)) #Theoretical length of all translatable amino acids added up
  
  #Get total count of mutations per patient (meatabolic and non-metabolic) [Tumor_Sample_Barcode, TOTAL.MUTATIONS]
  patient.all.mutations<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)), by="Tumor_Sample_Barcode"]
  
  ####PROCESS####
  
  #Construct matrix of patients (table.1) x Hugo_Symbol with values as total mutations per gene per patient
  patients.x.genes<-dcast(table.1, Hugo_Symbol~Tumor_Sample_Barcode, value.var="N.MUTATIONS", fun.aggregate=sum, fill=0)
  
  #Prep data.table for per gene across patient statistic calculation [Tumor_Sample_Barcode, GENE.MUTATIONS, Hugo_Symbol, TOTAL.MUTATIONS]
  patients.t.genes<-melt(patients.x.genes, id="Hugo_Symbol") #[Hugo_Symbol, variable, value]
  setnames(patients.t.genes, colnames(patients.t.genes), c("Hugo_Symbol","Tumor_Sample_Barcode","GENE.MUTATIONS"))
  patients.t.genes<-merge(patients.t.genes, as.data.frame(patient.all.mutations), by="Tumor_Sample_Barcode") #To get background-affected per patient
  
  #[Hugo_Symbol, Tumor_Sample_Barcode, GENE.MUTATIONS, TOTAL.MUTATIONS, Length]
  patients.t.genes<-merge(patients.t.genes, as.data.frame(length.table), by="Hugo_Symbol")
  patients.t.genes<-as.data.table(patients.t.genes)
  
  #Get patients that are used in the model (number of patients that we actually use and contain metabolic information)
  PATIENT.FOR.MODEL.COUNT<-length(unique(as.vector(patients.t.genes$Tumor_Sample_Barcode)))
  
  #Obtain rates
  patients.t.genes$GENE.RATE<-patients.t.genes$GENE.MUTATIONS/patients.t.genes$Length
  patients.t.genes$PATIENT.BMR<-patients.t.genes$TOTAL.MUTATIONS/background.length
  
  #STATS - Calculate student's t.test, paired t.test, mann whitney u and wilcoxon-signed ranked test per metabolite across patients
  GENE.STATS<-patients.t.genes[, Function.all.t.tests(GENE.RATE, PATIENT.BMR, ALT="greater", MU=0) , by="Hugo_Symbol"]
  
  #Correct for multiple hypothesis testing
  GENE.STATS$STUDENT.P.ADJ<-p.adjust(GENE.STATS$student.p.val, method="fdr")
  GENE.STATS$PAIRED.P.ADJ<-p.adjust(GENE.STATS$paired.t.p.val, method="fdr")
  GENE.STATS$MANN.U.P.ADJ<-p.adjust(GENE.STATS$mann.whitney.u.p.val, method="fdr")
  GENE.STATS$WILCOXON.P.ADJ<-p.adjust(GENE.STATS$wilcoxon.p.val, method="fdr")
  
  #Order by WILCOXON.SIGNED.RANK.TEST (NON-PARAMETRIC PAIRED T.TEST) (fdr corrected)
  GENE.STATS<-GENE.STATS[order(WILCOXON.P.ADJ),]
  
  #Get affected patients (patients that have mutations related to enriched metabolites)
  STUDENT.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in% GENE.STATS[STUDENT.P.ADJ<0.05,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  PAIRED.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  GENE.STATS[PAIRED.P.ADJ<0.05,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  MANN.U.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  GENE.STATS[MANN.U.P.ADJ<0.05,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  WILCOXON.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in% GENE.STATS[WILCOXON.P.ADJ<0.05,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  
  PATIENT.AFFECTED<-list(STUDENT.AFFECTED=STUDENT.AFFECTED,
                         PAIRED.AFFECTED=PAIRED.AFFECTED,
                         MANN.U.AFFECTED=MANN.U.AFFECTED,
                         WILCOXON.AFFECTED=WILCOXON.AFFECTED)
  
  #Return list including STATS, count of initial patients, patients used in model and patients per test affected (have mutations related to enriched j)
  dummy.return=list(STATS=GENE.STATS, STARTING.PATIENTS=TOTAL.PATIENT.COUNT, PATIENTS.IN.MODEL=PATIENT.FOR.MODEL.COUNT, PATIENT.AFFECTED=PATIENT.AFFECTED)
  return(dummy.return)
  
}

#PROCESS ENTRIES
require(base)
require(data.table)

args<-commandArgs(trailingOnly=T)

input.table.1.rds<-args[1]
length.file<-args[2] #File that contains a table of 2 columns (Hugo_Symbo, Length in amino acids)
output.file<-args[3]

table.1.input<-readRDS(input.table.1.rds)
length.file.input<-as.data.table(read.csv(length.file, header=T, sep="\t"))

#EXECUTE
Run<-Function.u.t.bmr.p(table.1.input$table.1, length.file.input, table.1.input$background.genes, remove.outliers=F)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)

####MODIFICATIONS####
