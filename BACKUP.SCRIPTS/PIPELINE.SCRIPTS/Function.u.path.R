#Function.u.path.R
#072214
#Compares the background mutation rate and the pathway mutation rate per pathway across patients
#It calculates significant mutations associated with pathways using 4 tests: student-t-test, paired-t-test, mann-whitney-u and wilcoxon signed-rank test
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

Function.u.t.bmr.path<-function(table.1, pathway.uniprot.file, uniprot.mapping.file, length.table, background.sequenced.genes, remove.outliers=F) {
  #Calculate u(j) weights per metabolite j 
  #This compares the bmr per patient of proteins associated with j versus a constant sequenced proteome bmr
  #This function assumes that all genes in table.2 are covered in length.table
  
  require(data.table)
  require(reshape)
  
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
  
  #Process pathway files
  uniprot.mapping<-fread(uniprot.mapping.file, header=T, sep="\t", drop=c(2:4, 6:7)) #[UNIPROT, Hugo_Symbol]  
  setnames(uniprot.mapping, colnames(uniprot.mapping), c("UNIPROT", "GENES"))
  uniprot.mapping$Hugo_Symbol<-sapply(uniprot.mapping$GENES, function(x) strsplit(x," ")[[1]][1])
  uniprot.mapping$GENES<-NULL
  
  pathway.uniprot<-as.data.table(read.csv(pathway.uniprot.file, header=F, sep="\t")) #[UNIPROT, PATH]
  pathway.uniprot<-pathway.uniprot[,c(1,4), with=F]
  setnames(pathway.uniprot, colnames(pathway.uniprot), c("UNIPROT", "PATH"))
  
  table.2<-as.data.table(merge(as.data.frame(uniprot.mapping), as.data.frame(pathway.uniprot), by="UNIPROT")) #[PATH, Hugo_Symbol]
  table.2<-unique(table.2[,c("PATH", "Hugo_Symbol"), with=F])
  
  #Get total background amino acid length (those that could be measured) - USED FOR TOTAL.BMR
  background.length<-length.table[Hugo_Symbol %in% background.sequenced.genes,]
  background.length<-sum(as.vector(background.length$Length)) #Theoretical length of all translatable amino acids added up
  
  #Get a per pathway length table [PATH, PATHWAY.LENGTH]
  pathway.length<-as.data.table(merge(as.data.frame(table.2), as.data.frame(length.table), by="Hugo_Symbol"))
  pathway.length<-pathway.length[,list(PATHWAY.LENGTH=sum(Length)), by="PATH"]
  
  #Get total count of mutations per patient (pathway and non-pathway) [Tumor_Sample_Barcode, TOTAL.MUTATIONS]
  patient.all.mutations<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)), by="Tumor_Sample_Barcode"]
  
  ####PROCESS####
  
  #Construct matrix of patients(table.1) x pathways(table.2) with values as total affected pathways mutations per patient per pathway
  patients.x.pathways<-merge(as.data.frame(table.1), as.data.frame(table.2), by="Hugo_Symbol")
  patients.x.pathways<-patients.x.pathways[,2:4] #Don't need gene column anymore
  patients.x.pathways<-cast(patients.x.pathways, Tumor_Sample_Barcode~PATH, value="N.MUTATIONS", sum, fill=0) #MATRIX
  
  #Prep data.table for per pathway across patient statistic calculations
  patients.t.pathways<-melt(patients.x.pathways, id="Tumor_Sample_Barcode") #[Tumor_Sample_Barcode, value, PATH]
  setnames(patients.t.pathways, colnames(patients.t.pathways), c("Tumor_Sample_Barcode", "PATHWAY.MUTATIONS", "PATH"))
  patients.t.pathways<-merge(patients.t.pathways, as.data.frame(patient.all.mutations), by="Tumor_Sample_Barcode") #To get background-affected per patient
  
  #[PATH, Tumor_Sample_Barcode, PATHWAY.MUTATIONS, PATH, TOTAL.MUTATIONS, PATHWAY.LENGTH]
  patients.t.pathways<-merge(patients.t.pathways, as.data.frame(pathway.length), by="PATH")
  patients.t.pathways<-as.data.table(patients.t.pathways)
  
  #Get patients that are used in the model (number of patients that we actually use and contain metabolic information)
  PATIENT.FOR.MODEL.COUNT<-length(unique(as.vector(patients.t.pathways$Tumor_Sample_Barcode)))
  
  #Obtain rates
  patients.t.pathways$PATHWAY.RATE<-patients.t.pathways$PATHWAY.MUTATIONS/patients.t.pathways$PATHWAY.LENGTH
  patients.t.pathways$PATIENT.BMR<-patients.t.pathways$TOTAL.MUTATIONS/background.length
  
  #STATS - Calculate student's t.test, paired t.test, mann whitney u and wilcoxon-signed ranked test per metabolite across patients
  PATH.STATS<-patients.t.pathways[, Function.all.t.tests(PATHWAY.RATE, PATIENT.BMR, ALT="greater", MU=0) , by="PATH"]
  
  #Correct for multiple hypothesis testing
  PATH.STATS$STUDENT.P.ADJ<-p.adjust(PATH.STATS$student.p.val, method="fdr")
  PATH.STATS$PAIRED.P.ADJ<-p.adjust(PATH.STATS$paired.t.p.val, method="fdr")
  PATH.STATS$MANN.U.P.ADJ<-p.adjust(PATH.STATS$mann.whitney.u.p.val, method="fdr")
  PATH.STATS$WILCOXON.P.ADJ<-p.adjust(PATH.STATS$wilcoxon.p.val, method="fdr")
  
  #Order by WILCOXON.SIGNED.RANK.TEST (NON-PARAMETRIC PAIRED T.TEST) (fdr corrected)
  PATH.STATS<-PATH.STATS[order(WILCOXON.P.ADJ),]
  
  #Get affected patients (patients that have mutations related to enriched pathways)
  STUDENT.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[PATH %in% PATH.STATS[STUDENT.P.ADJ<0.05,]$PATH ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  PAIRED.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[PATH %in% PATH.STATS[PAIRED.P.ADJ<0.05,]$PATH ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  MANN.U.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[PATH %in% PATH.STATS[MANN.U.P.ADJ<0.05,]$PATH ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  WILCOXON.AFFECTED<-unique(as.vector(table.1[Hugo_Symbol %in%  table.2[PATH %in% PATH.STATS[WILCOXON.P.ADJ<0.05,]$PATH ,]$Hugo_Symbol, ]$Tumor_Sample_Barcode))
  
  PATIENT.AFFECTED<-list(STUDENT.AFFECTED=STUDENT.AFFECTED,
                         PAIRED.AFFECTED=PAIRED.AFFECTED,
                         MANN.U.AFFECTED=MANN.U.AFFECTED,
                         WILCOXON.AFFECTED=WILCOXON.AFFECTED)
  
  #Return list including STATS, count of initial patients, patients used in model and patients per test affected (have mutations related to enriched j)
  dummy.return=list(STATS=PATH.STATS, STARTING.PATIENTS=TOTAL.PATIENT.COUNT, PATIENTS.IN.MODEL=PATIENT.FOR.MODEL.COUNT, PATIENT.AFFECTED=PATIENT.AFFECTED)
  return(dummy.return)
  
}

#PROCESS ENTRIES
require(base)
require(data.table)

args<-commandArgs(trailingOnly=T)

input.table.1.rds<-args[1]
input.pathway.uniprot.file<-args[2]
input.uniprot.mapping.file<-args[3]
length.file<-args[4] #File that contains a table of 2 columns (Hugo_Symbo, Length in amino acids)
output.file<-args[5]

table.1.input<-readRDS(input.table.1.rds)
length.file.input<-as.data.table(read.csv(length.file, header=T, sep="\t"))

#EXECUTE
Run<-Function.u.t.bmr.path(table.1.input$table.1, input.pathway.uniprot.file, input.uniprot.mapping.file, length.file.input, table.1.input$background.genes, remove.outliers=F)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)

