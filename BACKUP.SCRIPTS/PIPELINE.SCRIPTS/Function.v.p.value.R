#Function.v.p.value.R
#071814
#Calculates an empirical p.value for each v(j) using a pre-populated v.pval null distribution matrix
#This is based on the chance of getting an equal or greater ratio of differentially expressed genes if we were to pick the same number of samples that were used to
# to calculate v(j)
#This assaumes that the v.pval null distribution was constructed using Table.v
#Returns a table with p.values and fdr corrected pvalues per metabolite j including original v(j) feature

################################################################################################################################################################

Function.v.p.value<-function(processed.table.1.rds, processed.table.3.rds, v.table.rds, v.pval.null.rds) {
  #Per metabolite j:
  #   Find number of samples that have at least 1 mutation associated with KEGG(j)
  #   Sample from b$tumor number of such samples 100 times
  #   Do differential expression and construct null distribution per j
  #   Calculate p-value
  #Assumes that p.val.null contains all distributions for KEGG_IDs in v.table that had a mutation in at least 2 patients
  
  require(data.table)
  require(base)
  
  #Pre-Process tables
  table.1<-readRDS(processed.table.1.rds)
  table.1<-table.1$table.1[,1:2,with=F]
  table.1$PATIENT<-substr(table.1$Tumor_Sample_Barcode, 1,16)
  table.1$Tumor_Sample_Barcode<-NULL
  
  table.3<-readRDS(processed.table.3.rds)
  setnames(table.3, colnames(table.3), c("Hugo_Symbol", "KEGG_ID"))
  
  v.pval.null<-readRDS(v.pval.null.rds) #MATRIX!!! with column names as N.SAMPLES and rows as null distribution per N.SAMPLE
  
  v.table<-readRDS(v.table.rds)
  
  #Get number patients per metabolite [KEGG_ID, N.SAMPLES]
  patient.count.per.j<-as.data.table(merge(as.data.frame(table.1), as.data.frame(table.3), by="Hugo_Symbol"))
  patient.count.per.j<-patient.count.per.j[,list(N.SAMPLES=length(unique(PATIENT))), by="KEGG_ID"]
  
  #Filter out those that don't have at least one matched patient
  patient.count.per.j.0<-patient.count.per.j[N.SAMPLES<2,]
  patient.count.per.j.1<-patient.count.per.j[N.SAMPLES>1,]
  
  #Merge ">1" with v.table by KEGG_ID [KEGG_ID, v.METABOLITE, N.SAMPLES]
  v.table.1<-as.data.table(merge(as.data.frame(v.table), as.data.table(patient.count.per.j.1), by="KEGG_ID"))
  
  #Merge rest for later [KEGG_ID, v.METABOLITE, N.SAMPLES]
  v.table.0<-as.data.table(merge(as.data.frame(v.table), as.data.frame(patient.count.per.j.0), by="KEGG_ID"))
  
  #Calculate p.value
  sampling.number<-nrow(v.pval.null)
  v.table.1<-v.table.1[, list(v.METABOLITE=v.METABOLITE,
                              P.VAL=sum(v.pval.null[,as.character(N.SAMPLES)]>=v.METABOLITE)/sampling.number), by="KEGG_ID"]
  
  #Do multiple hypothesis testing [KEGG_ID, v.METABOLITE, P.VAL, ADJ.P.VAL]
  v.table.1$ADJ.P.VAL<-p.adjust(v.table.1$P.VAL, method="fdr")
  
  #Merge with "less than 2" KEGG_IDs
  v.table.0$P.VAL<-1
  v.table.0$ADJ.P.VAL<-1
  v.table.0$N.SAMPLES<-NULL
  dummy.result<-rbind(v.table.1, v.table.0)
  
  #Clean up and return
  dummy.result<-dummy.result[order(ADJ.P.VAL),]
  return(dummy.result)
}

#PROCESS ENTRIES
args<-commandArgs(trailingOnly=T)

input.table.1.rds<-args[1]
input.table.3.rds<-args[2]
input.table.v<-args[3]
input.table.p.vall.null<-args[4]
output.file<-args[5]

#EXECUTE
Run<-Function.v.p.value(input.table.1.rds, input.table.3.rds, input.table.v, input.table.p.vall.null)

#WRITE TO OUTPUT
saveRDS(Run, file=output.file)