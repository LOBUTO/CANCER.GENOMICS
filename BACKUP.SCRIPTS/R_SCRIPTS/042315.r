#042315.R
library(data.table)
library(pheatmap)

#Parse and polish recon datasets
Function.process.recon.simple<-function(gene.file, metabolite.file, reaction.file, filter=60){
  #Merge recon files to their simplest forms.
  #This will remove info on:
  # Compartment (location of molecule)
  # Get rid of all metabolite IDs but KEGG_ID
  # Wether a reaction is reversible or not
  # EC of a reaction enzyme
  
  require(data.table)
  
  #Process gene file [MODIFIERES, Hugo_Symbol]
  recon.genes<-fread(gene.file, header=T, sep="\t", stringsAsFactors=F)
  recon.genes<-recon.genes[Hugo_Symbol!="null",] #Filter out "null genes"
  recon.genes$COMPARTMENT<-NULL
  setkey(recon.genes)
  recon.genes<-unique(recon.genes)
  
  #Process metabolite file [RECON.ID, NAME, KEGG_ID]
  recon.metabolites<-fread(metabolite.file, header=T, sep="\t", stringsAsFactors=F)  
  #recon.metabolites$COMPARTMENT<-NULL
  #recon.metabolites$HMDB<-NULL #MODIFY
  #recon.metabolites$CHEBI<-NULL
  recon.metabolites$EHMN<-NULL
  setkey(recon.metabolites)
  recon.metabolites<-unique(recon.metabolites)
  
  #Apply size filtering
  recon.metabolites.a<-recon.metabolites[WEIGHT!="NONE",]
  recon.metabolites.b<-recon.metabolites[WEIGHT=="NONE",]
  recon.metabolites.a$WEIGHT<-as.numeric(recon.metabolites.a$WEIGHT)
  recon.metabolites.a<-recon.metabolites.a[WEIGHT>=filter,]
  recon.metabolites.a$WEIGHT<-NULL
  recon.metabolites.b$WEIGHT<-NULL
  recon.metabolites<-rbind(recon.metabolites.a, recon.metabolites.b)
  
  #Process reactions file
  recon.reactions<-fread(reaction.file, header=T, sep="\t", stringsAsFactors=F)
  recon.reactions$REVERSIBLE<-NULL
  recon.reactions$EC<-NULL
  
  #Merge first to obtain hugo symbol
  main.table<-merge(recon.reactions, recon.genes, by="MODIFIERS", allow.cartesian=T)
  main.table$MODIFIERS<-NULL
  main.table<-setkey(main.table)
  main.table<-unique(main.table)
  
  #Split into substrates table and product table
  substrate.table<-main.table[,setdiff(colnames(main.table), "PRODUCTS"), with=F]
  product.table<-main.table[,setdiff(colnames(main.table), "SUBSTRATES"), with=F]
  setkey(substrate.table)
  setkey(product.table)
  substrate.table<-unique(substrate.table)
  product.table<-unique(product.table)
  
  #Split recon molecule identifiers
  product.table<-product.table[,list(RECON.ID=unlist(strsplit(PRODUCTS,"|", fixed=T))), by=setdiff(colnames(product.table), "PRODUCTS")]
  substrate.table<-substrate.table[,list(RECON.ID=unlist(strsplit(SUBSTRATES,"|", fixed=T))), by=setdiff(colnames(substrate.table), "SUBSTRATES")]
  
  #Assign molecule names and kegg identifiers for product and substrate tables
  product.table<-merge(product.table, recon.metabolites, by="RECON.ID")
  substrate.table<-merge(substrate.table, recon.metabolites, by="RECON.ID")
  
  #Clean up and return
  setkey(product.table)
  setkey(substrate.table)
  product.table<-unique(product.table)
  substrate.table<-unique(substrate.table)
  return(list(PRODUCT=product.table, SUBSTRATE=substrate.table))
}

recon.table<-Function.process.recon.simple("DATABASES/RECON/042215.PROCESSED.GENE", "DATABASES/RECON/042215.PROCESSED.METABOLITES",
                                           "DATABASES/RECON/042215.PROCESSED.REACTIONS", filter=60)
recon.table$PRODUCT
recon.table$SUBSTRATE
recon.table$PRODUCT[,c("Hugo_Symbol","NAME" ,"KEGG_ID", "HMDB"),with=F]

test<-fread("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt", skip=3, header=F, stringsAsFactors=F, drop=c(1,3:43,45:109))
test<-test[V44 %in% c("Negative", "Positive"),]
test$V2<-gsub("-",".",test$V2)
test<-rbind(test[,list(sample=paste(V2, "01A", sep=".")  ), by="V44"], test[,list(sample=paste(V2, "01B", sep=".")  ), by="V44"])
test
dim(brca.agilent$tumor[,colnames(brca.agilent$tumor) %in% test$sample])

#Apply algorithm
Function.met.score<-function(sub.table, prod.table, agilent.obj, clinical.file){
  #Scores based on output from Function.process.recon.simple()
  
  #Process clinical file
  clinical<-fread(clinical.file, skip=3, header=F, stringsAsFactors=F, drop=c(1,3:43,45:109))
  setnames(clinical, c("SAMPLE", "ER"))
  clinical<-clinical[ER %in% c("Negative", "Positive"),]
  clinical$SAMPLE<-gsub("-",".",clinical$SAMPLE)
  clinical<-rbind(clinical[,list(SAMPLE=paste(SAMPLE, "01A", sep=".")  ), by="ER"], clinical[,list(SAMPLE=paste(SAMPLE, "01B", sep=".")  ), by="ER"])
  
  #Filter agilent tumor samples by er+/er- samples
  agilent.obj$tumor<-agilent.obj$tumor[,colnames(agilent.obj$tumor) %in% clinical$SAMPLE ]
  
  #Process agilent samples
  combined.matrices<-cbind(agilent.obj$normal, agilent.obj$tumor[rownames(agilent.obj$normal),])
  combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
  combined.matrices<-data.matrix(combined.matrices)
  
  #Proceed with differential expression
  design.matrix = data.frame(type=factor(c(rep("normal", ncol(agilent.obj$normal)), rep("cancer", ncol(agilent.obj$tumor)))) )
  design.matrix = model.matrix(~type, design.matrix)
  colnames(design.matrix)<- c("normal", "tumor")
  
  combined.matrices = normalizeBetweenArrays(combined.matrices, method="quantile")
  
  fit<-lmFit(combined.matrices, design.matrix)
  eb = eBayes(fit)
  exp.table<-data.table(topTable(eb, coef=2, n=Inf),keep.rownames=T)
  setnames(exp.table, c("Hugo_Symbol", setdiff(colnames(exp.table), "rn")))
  
  ####Do without prior filtering first####
  #Apply diff exp to ER- to ER+
  agilent.obj$er.plus<-agilent.obj$tumor[,colnames(agilent.obj$tumor) %in% clinical[ER=="Positive",]$SAMPLE]
  agilent.obj$er.neg<-agilent.obj$tumor[,colnames(agilent.obj$tumor) %in% clinical[ER=="Negative",]$SAMPLE]
  print (dim(agilent.obj$er.plus))
  print (dim(agilent.obj$er.neg))
  
  combined.matrices<-cbind(agilent.obj$er.plus, agilent.obj$er.neg[rownames(agilent.obj$er.plus),])
  combined.matrices<-combined.matrices[complete.cases(combined.matrices),]
  combined.matrices<-data.matrix(combined.matrices)
  print (dim(combined.matrices))
  combined.matrices<-combined.matrices[rownames(combined.matrices) %in% exp.table[adj.P.Val<0.1]$Hugo_Symbol,]#Filter for differentiated genes
  print (dim(combined.matrices))
  
  design.matrix = data.frame(type=factor(c(rep("plus", ncol(agilent.obj$er.plus)), rep("neg", ncol(agilent.obj$er.neg)))) )
  design.matrix = model.matrix(~type, design.matrix)
  colnames(design.matrix)<- c("plus", "neg")
  
  combined.matrices = normalizeBetweenArrays(combined.matrices, method="quantile")
  
  fit<-lmFit(combined.matrices, design.matrix)
  eb = eBayes(fit)
  exp.table<-data.table(topTable(eb, coef=2, n=Inf),keep.rownames=T)
  setnames(exp.table, c("Hugo_Symbol", setdiff(colnames(exp.table), "rn")))
  
  #Remove extracellular input and output so we don't account for metabolite levels outside of the cell
  sub.table<-sub.table[COMPARTMENT!="e",]
  prod.table<-prod.table[COMPARTMENT!="e",]
  sub.table$COMPARTMENT<-NULL
  prod.table$COMPARTMENT<-NULL
  setkey(sub.table)
  setkey(prod.table)
  sub.table<-unique(sub.table)
  prod.table<-unique(prod.table)
  
  #Filter recon tables by genes actually found in expression matrices
  prod.table<-prod.table[Hugo_Symbol %in% unique(exp.table$Hugo_Symbol),]
  sub.table<-sub.table[Hugo_Symbol %in% unique(exp.table$Hugo_Symbol),]
  
  #Integrate log folds
  prod.table<-merge(prod.table, exp.table[,c("Hugo_Symbol", "logFC"), with=F], by="Hugo_Symbol")
  sub.table<-merge(sub.table, exp.table[,c("Hugo_Symbol", "logFC"), with=F], by="Hugo_Symbol")
  
  #Convert to main table by median of each substrate and product enzymes categories
  prod.table<-prod.table[,list(PROD.LOG=median(logFC)), by="NAME"]
  sub.table<-sub.table[,list(SUB.LOG=median(logFC)), by="NAME"]
  main.table<-merge(prod.table, sub.table, by="NAME",all=T)
  main.table[is.na(main.table)]<-0
  
  #Obtain total fold
  main.table$PROD.LOG<-2^main.table$PROD.LOG
  main.table$SUB.LOG<-2^main.table$SUB.LOG
  main.table$MET.LOGFC<-log2(main.table$PROD.LOG / main.table$SUB.LOG)
  
  ##Clean up and Return
  main.table<-main.table[order(abs(MET.LOGFC), decreasing=T),]
  return(main.table)
}

met.score<-Function.met.score(recon.table$SUBSTRATE, recon.table$PRODUCT, brca.agilent, 
                              "DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/031315/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
met.score[1:20,]
hist(met.score$MET.LOGFC)
ER.GS
met.score[grepl("methionine", NAME, ignore.case=T) ,]

recon.table$PRODUCT[NAME=="2-methylbutyrylglycine",]
exp.table[Hugo_Symbol=="GLYAT",]


###Introduce gold standard between ER+ and ER-##
#####From Terunuma#######
GS.TER<-read.csv("DATABASES/METABOLOMICS/TERUNUMA.2014/cleaned.met.csv", header=T, skip=1, stringsAsFactors=F)
rownames(GS.TER)<-GS.TER$STATUS
GS.TER$STATUS<-NULL
GS.TER<-as.matrix(GS.TER)
GS.TER<-GS.TER[apply(GS.TER, 1, sd )!=0, ] #remove constant variables
GS.TER<-GS.TER[!grepl("X - ", rownames(GS.TER)),]

#First filter for normally diff exp metabolites (normal vs cancer)
TER.METABOLITES<-rownames(GS.TER)
NORMAL<-colnames(GS.TER)[grepl("NORMAL", colnames(GS.TER))]
CANCER<-setdiff(colnames(GS.TER), NORMAL)
GS.TER.DIFF<-data.table(MET=rownames(GS.TER),
                        P.VAL=apply(GS.TER, 1, function(x) wilcox.test(x[NORMAL], x[CANCER], paired=F)$p.value))
GS.TER.DIFF$P.VAL.ADJ<-p.adjust(GS.TER.DIFF$P.VAL, method="fdr")
GS.TER<-GS.TER[GS.TER.DIFF[P.VAL.ADJ<0.1,]$MET,] #at 0.1 fdr

#Then do ER metabolite marker
ER.NEG<-colnames(GS.TER)[grepl("NEG", colnames(GS.TER))]
ER.POS<-colnames(GS.TER)[grepl("POS", colnames(GS.TER))]
      #heatplot(log(GS.TER), scale="none")
GS.TER<-data.table(metabolite=rownames(GS.TER), 
                   P.VAL=apply(GS.TER, 1, function(x) wilcox.test(x[ER.NEG], x[ER.POS], paired=F, alternative="greater")$p.value),
                   FC=  apply(GS.TER, 1, function(x) median(x[ER.NEG])/median(x[ER.POS])))
GS.TER$P.VAL.ADJ<-p.adjust(GS.TER$P.VAL, method="fdr")
GS.TER<-GS.TER[order(P.VAL.ADJ),]
GS.TER[P.VAL.ADJ<0.1,]
GS.TER[grepl("glutarate", metabolite),]

#######From Tang######
GS.TANG<-fread("DATABASES/METABOLOMICS/TANG.2014/clean.met.er.csv", header=T, skip=1)
setnames(GS.TANG, as.vector(sapply(colnames(GS.TANG), function(x) unlist(strsplit(x,"/"))[1])))
colnames(GS.TANG)<-ifelse(colnames(GS.TANG)=="ER-", "ER.NEG", 
                          ifelse(colnames(GS.TANG)=="ER+", "ER.POS", colnames(GS.TANG))) 
GS.TANG<-data.frame(GS.TANG)
rownames(GS.TANG)<-GS.TANG$TYPE
GS.TANG$TYPE<-NULL
GS.TANG<-as.matrix(GS.TANG)
GS.TANG<-GS.TANG[apply(GS.TANG, 1, sd )!=0, ] #remove constant variables

#First filter for normally diff exp metabolites (normal vs cancer)
NORMAL<-colnames(GS.TANG)[grepl("Normal", colnames(GS.TANG))]
CANCER<-setdiff(colnames(GS.TANG), NORMAL)
GS.TANG.DIFF<-data.table(MET=rownames(GS.TANG),
                        P.VAL=apply(GS.TANG, 1, function(x) wilcox.test(x[NORMAL], x[CANCER], paired=F)$p.value))
GS.TANG.DIFF$P.VAL.ADJ<-p.adjust(GS.TANG.DIFF$P.VAL, method="fdr")
GS.TANG<-GS.TANG[GS.TANG.DIFF[P.VAL.ADJ<0.1,]$MET,] #at 0.1 fdr

#Then do ER metabolite marker
ER.NEG<-colnames(GS.TANG)[grepl("NEG",colnames(GS.TANG))]
ER.POS<-colnames(GS.TANG)[grepl("POS",colnames(GS.TANG))]
      #heatplot(log(GS.TANG), scale="none")
GS.TANG<-data.table(metabolite=rownames(GS.TANG), 
                    P.VAL=apply(GS.TANG, 1, function(x) wilcox.test(x[ER.NEG], x[ER.POS], paired=F, alternative="greater")$p.value ),
                    FC=  apply(GS.TANG, 1, function(x) median(x[ER.NEG])/median(x[ER.POS])))
GS.TANG<-GS.TANG[!is.na(P.VAL),]
GS.TANG$P.VAL.ADJ<-p.adjust(GS.TANG$P.VAL, method="fdr")
GS.TANG[P.VAL.ADJ<0.05,]
dim(GS.TANG)

intersect(GS.TER[P.VAL.ADJ<0.1,]$metabolite, GS.TANG[P.VAL.ADJ<0.1,]$metabolite) #22, 5
setdiff(GS.TANG[P.VAL.ADJ<0.1,]$metabolite, GS.TER[P.VAL.ADJ<0.1,]$metabolite)
setdiff(GS.TANG[P.VAL.ADJ<0.1,]$metabolite, GS.TER$metabolite)
hist(GS.TANG[metabolite %in% GS.TER[P.VAL.ADJ<0.1]$metabolite,]$P.VAL.ADJ)

COMPARE.TANG<-GS.TANG[metabolite %in% TER.METABOLITES,]
COMPARE.TANG$TER.SIG<-COMPARE.TANG$metabolite %in% GS.TER[P.VAL.ADJ<0.1,]$metabolite
COMPARE.TANG$COMMON<-COMPARE.TANG$metabolite %in% intersect(GS.TER[P.VAL.ADJ<0.1,]$metabolite, GS.TANG[P.VAL.ADJ<0.1,]$metabolite)

wilcox.test(-log(COMPARE.TANG[COMMON==T,]$P.VAL.ADJ),-log(COMPARE.TANG[COMMON==F,]$P.VAL.ADJ), alternative="greater")
wilcox.test(COMPARE.TANG[COMMON==T,]$FC,COMPARE.TANG[COMMON==F,]$FC, alternative="greater")

COMPARE.TANG$NEG.LOG.PVAL<--log(COMPARE.TANG$P.VAL.ADJ)
COMPARE.TANG<-melt(COMPARE.TANG[,c(3,5:7),with=F], id.vars=c("COMMON", "TER.SIG"))

ggplot(COMPARE.TANG, aes(TER.SIG, value, colour=TER.SIG)) + geom_boxplot() + geom_jitter() + theme.format + facet_wrap(~variable)
ggplot(COMPARE.TANG, aes(COMMON, value, colour=COMMON)) + geom_boxplot() + geom_jitter() + theme.format + facet_wrap(~variable)

####Finally setttle on ER+/ER- Gold standard####
ER.GS<-intersect(GS.TER[P.VAL.ADJ<0.1,]$metabolite, GS.TANG[P.VAL.ADJ<0.1,]$metabolite) #22 for now


######Predicting 2HG patients######

