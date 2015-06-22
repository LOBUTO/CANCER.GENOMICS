####FIND SIGNIFICANCE OF METABOLIC NODE WEIGHTS#####
##This will be initially approached by using clinical data

#########FUNCTIONS#########
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Function.BRCA.clinical.type<-function(file) {
  #Take in clinical BRCA file of the type "nationwidechildrens.org_clinical_patient_brca.txt"
  #Returns data.table of patients classified by subtype and plot of counts per subtype
  
  #Load and process file
  a<-read.csv(file, header=T,sep="\t",stringsAsFactors=FALSE)
  a<-a[-c(1,2),]
  a<-a[,c("bcr_patient_barcode", "er_status_by_ihc","pr_status_by_ihc","her2_status_by_ihc")]
  a<-as.data.table(a)
  setnames(a, colnames(a), c(colnames(a)[1], as.vector(sapply(colnames(a)[2:4], function(x) strsplit(x,"_")[[1]][1]))))
  
  #Remove non-conclusive lines
  filter=c("Positive","Negative")
  a<-a[er %in% filter & pr %in% filter & her2 %in% filter,]
  
  #Classify samples by combined subtype
  a[a=='Negative']<-"-"
  a[a=='Positive']<-"+"
  a<-a[,list(er=paste0("er",er,""),
             pr=paste0("pr",pr,""),
             her2=paste0("her2",her2,"")), by="bcr_patient_barcode"]
  a$subtype<-paste0(a$er,"_",a$pr,"_",a$her2,"")
  
  #Create plot
  d<-a[,list(count=length(bcr_patient_barcode)), by="subtype"]
  d$SINGLE<-sapply(d$subtype, function(x)  length(which(strsplit(x,"")[[1]]=="+"))==1 )
  PLOT.1<-ggplot(d, aes(x=subtype, y=count, fill=SINGLE)) + geom_bar(stat="identity") + theme.format+
    geom_text(aes(label=count), vjust=-0.25, size=14)
  
  return(list(C.TABLE=a,PLOT=PLOT.1))
}

Function.read.RNAseq.files<- function(folder, processed.map.matrix, cancer.sep=T) {
  #Produces matrices of genesxpatient.samples based on folder locations of RNAseq files and processed.map.matrix obtained from
  # Function.process.RNAseq.map.files()
  # Could produce 2 matrices corresponding to normal 11A and cancer 01A or simply a single matrix 
  #processed.map.matrix has first column as filename and second column as patient.sample
  
  Internal.Function.1<-function(input.matrix, folder){
    #Merges RNAseq files in "folder" depending on instructions from input matrix
    
    #Read files into list
    matrices<-lapply(1:nrow(input.matrix), function(f) read.csv(paste0(folder,"/", input.matrix[f,1],""), header=T, sep="\t",
                                                                        col.names=c("gene_id", input.matrix[f,2])))
    
    #Merge vector matrices to create expression table
    dummy.expression.table<-join_all(matrices, by="gene_id", type="inner")
    dummy.expression.table$gene_id<-as.character(dummy.expression.table$gene_id)
    
    #Minor fix to account for ANNOTATION ERROR
    dummy.expression.table[dummy.expression.table=="SLC35E2|728661"]<-"SLC35E2B|728661"
    
    #Remove "?" genes
    dummy.expression.table$gene<-sapply(dummy.expression.table$gene_id, function(x) strsplit(x, "[|]")[[1]][1])
    dummy.expression.table<-dummy.expression.table[dummy.expression.table$gene!="?",]
    
    #Assign genes as rownames and clean up
    rownames(dummy.expression.table)<-dummy.expression.table$gene
    dummy.expression.table$gene<-NULL
    dummy.expression.table$gene_id<-NULL
    dummy.expression.table<-dummy.expression.table[complete.cases(dummy.expression.table),]
    
    #Return
    return(dummy.expression.table)
  }
  
  require(plyr)
  
  if (cancer.sep==T) {
    processed.map.matrix<-as.data.frame(processed.map.matrix, stringAsFactors=F)
    
    #Separate based on normal or cancer
    processed.map.matrix$type<-substr(processed.map.matrix$patient.sample, 14,15)
    processed.map.matrix.normal<-processed.map.matrix[processed.map.matrix$type=="11",]
    processed.map.matrix.cancer<-processed.map.matrix[processed.map.matrix$type=="01",]
    
    dummy.normal<-Internal.Function.1(as.matrix(processed.map.matrix.normal), folder)
    dummy.cancer<-Internal.Function.1(as.matrix(processed.map.matrix.cancer), folder)
    
    dummy.return<-list(dummy.normal, dummy.cancer)
    names(dummy.return)<-c("normal", "tumor")
    
  } else if (cancer.sep==F) {
    dummy.return<-Internal.Function.1(processed.map.matrix, folder)
  } else
    print ("Please choose correct cancer.sep")
  
  #Return
  return(dummy.return)
}

Function.process.RNAseq.map.files<-function(map.file) {  
  #Processes map.file from RNAseq data
  #Filters for gene name containing files and shortens barcode to represent patient.sample
  #Produces a 2-column matrix of "filename" and "patient.sample" columns
  
  dummy.map<-read.csv(map.file, header=T, sep="\t")
 
  #Keep only files that are for processed genes
  dummy.map<-dummy.map[grepl("rsem.genes.normalized_results", dummy.map$filename),]
 
  #Process sample name from barcode
  dummy.map$patient.sample<-substr(dummy.map$barcode.s.,1,16)
  dummy.map$barcode.s.<-NULL
 
  #Return
  return(as.matrix(dummy.map))
}

Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  
  require(plyr)
  require(limma)
  require(edgeR)
  
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 20 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
                                # upper quantile normalized (TCGA)
  
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data 
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Return topTable
  return(all.G.fit)
}

Function.Expression.Heatmap<-function(normal.matrix, cancer.matrix) {

}

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

#######060414######
#Catalog patients based on subtype
e<-Function.BRCA.clinical.type("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/99d90670-467c-4c6b-b870-027712c19014/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
unique(as.vector(e[[1]]$subtype))

#See if subtypes can be achieved by microarray information only
a<-read.csv("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
head(a)
a$patient<-substr(a$barcode.s., 1, 12)
length(unique(as.vector(a$patient)))

target.subtypes<-c("er-_pr-_her2-", "er-_pr-_her2+", "er+_pr-_her2-", "er+_pr+_her2+")
f<-e[[1]]
f<-f[subtype %in% target.subtypes,]
f<-f[bcr_patient_barcode %in% a$patient,]
ggplot(f, aes(subtype, fill=factor(subtype))) + geom_bar() + theme.format

#######060614#######
d<-Function.read.RNAseq.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/UCEC/OS_CONTROL/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3",cancer.sep=F,
                              Function.process.RNAseq.map.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/UCEC/OS_CONTROL/FILE_SAMPLE_MAP.txt"))
b<-Function.read.RNAseq.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3",cancer.sep=T,
                              Function.process.RNAseq.map.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/FILE_SAMPLE_MAP.txt"))
dim(b$tumor)
head(as.matrix(b$normal)[,1:3])

#######060914########
g<-Function.BRCA.clinical.type("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/060914_BRCA_CLINICAL/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
g$PLOT
h<-Function.RNAseq.Differential.Expression(b$normal, b$tumor)
hist(h$adj.P.Val)
head(h)

b$normal$rn<-rownames(b$normal)
b$tumor$rn<-rownames(b$tumor)

#Normalize and discard low cpm before doing it
i<-join_all(list(b$normal, b$tumor), by="rn", type="inner")
rownames(i)<-i$rn
i$rn<-NULL

i.isexpr<-rowSums(cpm(i)>1)>=20
i<-i[i.isexpr,]
i<-normalizeBetweenArrays(as.matrix(i), method="quantile")

ncol(b$normal), ncol(b$tumor)

heatmap.2(i, trace="none",
          ColSideColors=c(rep("green",ncol(b$normal)-1), rep("blue", ncol(b$tumor)-1)))

i.annotation<-data.frame(samples=colnames(i))
i.annotation$TYPE<-substr(i.annotation$samples,14,15)
i.annotation[i.annotation=="11"]<-"Normal"
i.annotation[i.annotation=="01"]<-"Cancer"
rownames(i.annotation)<-i.annotation$samples
i.annotation$samples<-NULL
i.annotation$TYPE<-as.factor(i.annotation$TYPE)
pheatmap(i,annotation=i.annotation)

#######061014#######

Function.process.maf.files<-function(main.file=c(), additional.files=c()) {
  #Process maf.files from TCGA to 2 objects:
  #First is a table containing:
  #   Tumor_Sample_Barcode
  #   Hugo_Symbol
  #   N.MUTATIONS per gene per patient sample
  #Second is a vector containing all theoretically sequenced genes
  #   This is regardless of type ("silent"), since this are all possible sequenced genes
  
  require(data.table)
  
  #Process main.file
  dummy.a<-read.csv(main.file, header=T,  sep="\t")
  dummy.a$Line_Number<-NULL #To remove duplicates
  dummy.a<-unique(dummy.a)
  SEQUENCED.GENES<-unique(as.vector(dummy.a$Hugo_Symbol)) #ALL POTENTIALLY SEQUENCED GENES
  dummy.a<-dummy.a[dummy.a$Variant_Classification!="Silent",] #Remove silent mutations
  
  #Filter for info we want
  dummy.a<-dummy.a[,c("Tumor_Sample_Barcode", "Hugo_Symbol")] #Only want sample names and mutations
  dummy.a<-as.data.table(dummy.a)
  dummy.a$fill<-1 #To count number of mutations per gene
  dummy.a<-dummy.a[,list(N.MUTATIONS=length(fill)), by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]
  
  #Process additional files:
  for (files in additional.files) {
    a.1<-read.csv(files, header=T, sep="\t")
    
    #Add to background genes
    SEQUENCED.GENES<-unique(c(SEQUENCED.GENES, unique(as.vector(a.1$Hugo_Symbol))))
    
    #Only keep samples not found in main.file
    a.1<-a.1[!(a.1$Tumor_Sample_Barcode %in% unique(as.vector(dummy.a$Tumor_Sample_Barcode))),]
    
    #Continue processing
    a.1$Line_Number<-NULL
    a.1<-unique(a.1)
    a.1<-a.1[a.1$Variant_Classification!="Silent",]
    a.1<-a.1[,c("Tumor_Sample_Barcode", "Hugo_Symbol")] #Only want sample names and mutations
    a.1<-as.data.table(a.1)
    a.1$fill<-1 #To count number of mutations per gene
    a.1<-a.1[,list(N.MUTATIONS=length(fill)), by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]  
    
    #Add to main file
    dummy.a<-rbind(dummy.a,a.1)
  }
  
  #Return
  dummy.list<-list(table.1=dummy.a, background.genes=SEQUENCED.GENES)
  return (dummy.list)
}

Function.post.process.table.2<-function(table.2, KEGG.IDS) {
  #Filters table.2 for KEGG.IDS, preferrably from table.3
  #Table.2 of the form BRCA.table.2 (data.table with KEGG_ID column)
  
  require (data.table)
  table.2<-table.2[KEGG_ID %in% KEGG.IDS,]
  return(table.2)
}

Function.post.process.table.3<-function(table.3, output.file) {
  #Small fine tuning
  #   unique() to filter for duplicates obtained from 
  #   Filters out "Product" column from table.3
  #   Filters out water and H2O
  #Table.3 of the form BRCA.Table.3 (data.table with 3 columns)
  
  table.3<-unique(table.3[,c(1,2),with=F])
  table.3<-table.3[table.3$Product!="H2O",]
  table.3<-table.3[table.3$Product!="Water",]
  
  return(table.3)
}

#######061614######
Function.u.bg<-function(table.1, table.2, length.table, background.sequenced.genes, remove.outliers=F) {
  #Calculates significance of mutations related to metabolite using amino acid background mutation rate
  #Returns adjusted p.value of hypergeometric test per metabolite
  
  require(data.table)
  require(reshape)
  
  hyper.function<-function(metabolic.mutations, total.metabolic.sequence,
                           all.mutations, all.sequences){
    dummy.hyper<-phyper(metabolic.mutations-1, background.mutations,
                        all.sequences*3-all.mutations, total.metabolic.sequence*3, lower.tail=F)
    return(list(p.hyper=dummy.hyper))
  }
  
  #Remove outliers if necessary
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }
  
  #Prep tables
  length.table<-as.data.table(length.table)
  setnames(length.table, colnames(length.table), c("Hugo_Symbol", "Length"))
  
  table.2<-unique(table.2[,2:3, with=F]) #Remove Name column to remove duplicates
  setnames(table.2, colnames(table.2), c("KEGG_ID", "Hugo_Symbol"))
  
  #Filter for those where we have length information (Theoretically we are filtering out pseudogenes, though may not be necessarily correct)
  table.1<-table.1[Hugo_Symbol %in% length.table$Hugo_Symbol,]
  
  #Get background mutation rate - This only accounts for protein lengths we can tell, if we can't then we don't count it in
  background.merge<-merge(as.data.frame(table.1), as.data.frame(length.table), by="Hugo_Symbol")
  background.mutations<-sum(background.merge$N.MUTATIONS) #TOTAL MUTATIONS ACROSS ALL TUMORS
  
  background.total.length<-sum(as.vector(length.table[Hugo_Symbol %in% background.sequenced.genes,]$Length))
  sequenced.samples<-length(unique(as.vector(table.1$Tumor_Sample_Barcode))) #Number of sequence samples
  background.total.length<-as.numeric(background.total.length)*as.numeric(sequenced.samples) #TOTAL LENGTH OF SEQUENCED SAMPLES
  
  #Get sequences associated with metabolite - This assumes that all genes in table.2 are covered by length.table
  metabolite.sequences<-as.data.table(merge(as.data.frame(table.2), as.data.frame(length.table), by="Hugo_Symbol")) #Length per protein ass. with metabolite
  metabolite.sequences<-metabolite.sequences[,list(METABOLIC.SEQUENCE=sum(Length)),by="KEGG_ID"] #Total aa length per metabolite
  metabolite.sequences$METABOLIC.SEQUENCE<-metabolite.sequences$METABOLIC.SEQUENCE*sequenced.samples #TOTAL LENGTH OF SEQUENCED SAMPLES ASSOCIATED WITH METABOLITES
  
  #Get sum of mutations associated with metabolites across patients
  patient.x.metabolite<-merge(as.data.frame(table.1), as.data.frame(table.2), by="Hugo_Symbol")
  patient.x.metabolite<-patient.x.metabolite[,2:4] #Don't need gene column anymore
  patient.x.metabolite<-cast(patient.x.metabolite, Tumor_Sample_Barcode~KEGG_ID, value="N.MUTATIONS", sum) #MATRIX
  patient.t.metabolite<-melt(patient.x.metabolite, id="Tumor_Sample_Barcode") #TABLE PER PATIENT ACROSS METABOLITE MUTATIONS
  setnames(patient.t.metabolite, colnames(patient.t.metabolite), c("Tumor_Sample_Barcode", "METABOLIC.MUTATIONS", "KEGG_ID"))
  
  #Merge to get mutations seen and total sequence per metabolite
  patient.t.metabolite<-as.data.table(merge(as.data.frame(patient.t.metabolite), as.data.frame(metabolite.sequences), by="KEGG_ID"))
  
  #Calculate
  metabolite.hyper<-patient.t.metabolite[,hyper.function(sum(METABOLIC.MUTATIONS),mean(METABOLIC.SEQUENCE),
                                                         background.mutations, background.total.length),by="KEGG_ID"]
  #Correct for multiple hypothesis
  metabolite.hyper$p.hyper.adj<-p.adjust(metabolite.hyper$p.hyper, method="fdr")
  
  #Clean up
  metabolite.hyper<-metabolite.hyper[order(p.hyper.adj),]
  
  #Return
  return(metabolite.hyper)
}

dummy.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds")
dummy.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
dummy.table.length<-as.data.table(read.csv("DATABASES/UNIPROT/042314_GENE_LENGTH", header=T, sep="\t"))

test<-Function.u.bg(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes,remove.outliers=T)
hist(test$p.hyper.adj)
test[p.hyper.adj<0.05,]

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

#######061714#######
dummy.clinical<-Function.BRCA.clinical.type("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/060914_BRCA_CLINICAL/Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt")
test.1<-split(dummy.clinical$C.TABLE, dummy.clinical$C.TABLE$subtype)
test.1$`er-_pr-_her2-`

Function.u.bg.per.subtype<-function(table.1, table.2, length.table, background.sequenced.genes, clinical.data.table, remove.outliers=F) {
  #Does Function.u.bg() calculation for each subtype found in clinical.data.table
  #Clinical.data.table is obtained from Function.BRCA.clinical.type() function
  #   We assume at the moment that BRCA is the only cancer that has clinical subtype documented in TCGA
  
  require(data.table)
  
  #Prep clinical table - Only care about bcr_patient_code and subtype columnds
  clinical.data.table<-clinical.data.table[,c(1,5),with=F]
  
  #Classify and filter patients by subtype 
  table.1$bcr_patient_barcode<-substr(table.1$Tumor_Sample_Barcode,1,12)
  table.1<-merge(as.data.frame(table.1), as.data.frame(clinical.data.table), by="bcr_patient_barcode")
  
  #Apply Function.u.bg()
  table.1.subtypes.list<-split(table.1, table.1$subtype) #list of data.frames
  table.1.u.by.subtype<-lapply(table.1.subtypes.list, function(x) 
                               Function.u.bg(as.data.table(x[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "N.MUTATIONS")]), 
                                             table.2, length.table, background.sequenced.genes, remove.outliers)
                               )
  
  #Return
  return(table.1.u.by.subtype)
  
}

test.clinical<-Function.u.bg.per.subtype(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, dummy.clinical$C.TABLE,remove.outliers=T)
test.clinical[[1]]
test.clinical$`er+_pr-_her2-`

Function.u.ratio<-function(table.1, table.2, length.table, background.sequenced.genes, remove.outliers=F) {
  #Calculate significance per metabolite based on probability of finding mutations related to metabolite give sequence length ratio
  
  require(data.table)
  
  #Remove outliers if necessary
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }
  
  #Prep tables
  table.2<-unique(table.2[,c(2,3),with=F])
  setnames(table.2, colnames(table.2), c("KEGG_ID", "Hugo_Symbol"))
  table.1<-table.1[Hugo_Symbol %in% length.table$Hugo_Symbol,]
  
  #Get background total amino acid length per metabolite <- 2 COLUMNS [KEGG_ID, metabolic.aa.length]
  metabolite.lengths<-as.data.table(merge(as.data.frame(table.2), as.data.frame(length.table), by="Hugo_Symbol"))
  metabolite.lengths<-metabolite.lengths[,list(metabolic.aa.length=sum(Length)), by="KEGG_ID"]
  
  #Get background total amino acid length for sequenceable proteome
  proteome.aa.length<-length.table[Hugo_Symbol %in% background.sequenced.genes,]
  proteome.aa.length<-sum(proteome.aa.length$Length)
  
  #Get total mutations across patient (of those proteins we have length information for) 
  table.1.count<-table.1[,list(TOTAL.MUTATIONS=sum(N.MUTATIONS)), by="Tumor_Sample_Barcode"]
  table.1.mutations<-sum(table.1.count$TOTAL.MUTATIONS)
  
  #Get total mutations per metabolite across patients - [KEGG_ID, METABOLIC.MUTATIONS]
  table.1.metabolite.count<-as.data.table(merge(as.data.frame(table.1), as.data.frame(table.2), by="Hugo_Symbol"))
  table.1.metabolite.count<-table.1.metabolite.count[,list(METABOLIC.MUTATIONS=sum(N.MUTATIONS)), by="KEGG_ID"]
  
  #Set up table for hyper [KEGG_ID, METABOLIC.MUTATIONS, metabolic.aa.length]
  hyper.table<-as.data.table(merge(as.data.frame(table.1.metabolite.count), as.data.frame(metabolite.lengths), by="KEGG_ID"))
  hyper.table$p.hyper<-phyper(hyper.table$METABOLIC.MUTATIONS-1, hyper.table$metabolic.aa.length,
                              proteome.aa.length-hyper.table$metabolic.aa.length, table.1.mutations, lower.tail=F)
  hyper.table$p.hyper.adj<-p.adjust(hyper.table$p.hyper, method="fdr")
  hyper.table<-hyper.table[order(p.hyper.adj),]
  
  #Return
  return(hyper.table)
  
}

dummy.u.ratio<-Function.u.ratio(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=T)
ggplot(dummy.u.ratio, aes(x=p.hyper.adj)) + geom_histogram()

Function.u.ratio.per.subtype<-function(table.1, table.2, length.table, background.sequenced.genes, clinical.data.table, remove.outliers=F) {
  #Does Function.u.ratio() calculation for each subtype found in clinical.data.table
  #Clinical.data.table is obtained from Function.BRCA.clinical.type() function
  #   We assume at the moment that BRCA is the only cancer that has clinical subtype documented in TCGA
  
  require(data.table)
  
  #Prep clinical table - Only care about bcr_patient_code and subtype columnds
  clinical.data.table<-clinical.data.table[,c(1,5),with=F]
  
  #Classify and filter patients by subtype 
  table.1$bcr_patient_barcode<-substr(table.1$Tumor_Sample_Barcode,1,12)
  table.1<-merge(as.data.frame(table.1), as.data.frame(clinical.data.table), by="bcr_patient_barcode")
  
  #Apply Function.u.bg()
  table.1.subtypes.list<-split(table.1, table.1$subtype) #list of data.frames
  table.1.u.by.subtype<-lapply(table.1.subtypes.list, function(x) 
    Function.u.ratio(as.data.table(x[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "N.MUTATIONS")]), 
                  table.2, length.table, background.sequenced.genes, remove.outliers)
  )
  
  #Return
  return(table.1.u.by.subtype)
}

dummy.u.ratio.subtype<-Function.u.ratio.per.subtype(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, dummy.clinical$C.TABLE,remove.outliers=T)
dummy.u.ratio.subtype[[7]]

Function.u<-function(table.1, table.2, length.table, background.sequenced.genes, method="mean",remove.outliers=F) {
  ###
  
  require(reshape)
  require(data.table)
  
  #Remove outliers if necessary
  if (remove.outliers==T) {
    table.1<-Function.Outliers.Table.1(table.1)
  }
  
  #Prep tables
  table.2<-unique(table.2[,c(2,3), with=F])
  setnames(table.2, colnames(table.2), c("KEGG_ID", "Hugo_Symbol"))
  length.table<-as.data.table(length.table)
  setnames(length.table, colnames(length.table), c("Hugo_Symbol", "Length"))
  
  #Get length of metabolic genes [KEGG_ID, METABOLIC.LENGTH]
  table.2.length<-as.data.table(merge(as.data.frame(table.2), as.data.frame(length.table), by="Hugo_Symbol"))
  table.2.length<-table.2.length[,list(METABOLIC.LENGTH=sum(Length)), by="KEGG_ID"]
  
  #Build patientxsample matrix based on metabolic mutated genes total length (common length)
  patients.x.metabolites<-merge(as.data.frame(table.1), as.data.frame(table.2), by="Hugo_Symbol")  
  patients.x.metabolites<-as.data.table(merge(patients.x.metabolites, as.data.frame(length.table), by="Hugo_Symbol")) #add length to calculate common metabolic length
  patients.x.metabolites<-patients.x.metabolites[,list(COMMON.LENGTH=sum(Length)), by=c("Tumor_Sample_Barcode", "KEGG_ID")] #3 columns, incomplete
  patients.x.metabolites<-cast(patients.x.metabolites, Tumor_Sample_Barcode~KEGG_ID, value="COMMON.LENGTH",fill=0) #MATRIX,complete
  
  #Prep data.table for u(j) calculation [Tumor_Sample_Barcode, COMMON.LENGTH, KEGG_ID]
  patients.t.metabolites<-melt(patients.x.metabolites, id="Tumor_Sample_Barcode")
  setnames(patients.t.metabolites, colnames(patients.t.metabolites), c("Tumor_Sample_Barcode", "COMMON.LENGTH", "KEGG_ID"))
  
  #Add METABOLIC LENGTH [Tumor_Sample_Barcode, COMMON.LENGTH, KEGG_ID, METABOLIC.LENGTH]
  patients.t.metabolites<-as.data.table(merge(as.data.frame(patients.t.metabolites), as.data.frame(table.2.length), by="KEGG_ID"))
  
  #Calculate u(j) based on method of choice (mean/median) [KEGG_ID, METABOLITE.U]
  if (method=="mean") {
    metabolites.u<-patients.t.metabolites[,list(METABOLITE.U= mean(COMMON.LENGTH/METABOLIC.LENGTH)), by="KEGG_ID"]
  } else if (method=="median") {
    metabolites.u<-patients.t.metabolites[,list(METABOLITE.U= median(COMMON.LENGTH/METABOLIC.LENGTH)), by="KEGG_ID"]
  }
  
  #ASSIGN SIGNIFICANCE 
  #This is based on random sampling of mutated genes across all patients
  Function.u.significance<-function(n.metabolic.genes, background.genes, sig.table, table.length,KEGG) {
    print(KEGG)
    total.patients.count<-length(unique(as.vector(sig.table$Tumor_Sample_Barcode))) #Total number of patients we started with
    sample.genes<-sample(background.genes, n.metabolic.genes,replace=T) #Sampled genes for test
    sample.genes.length<-sum(table.length[Hugo_Symbol %in% sample.genes, ]$Length) #Total length of sampled genes
    
    #Measure mock u(j) across patients
    common.patients.calc<-sig.table[Hugo_Symbol %in% sample.genes,]
    common.patients.calc<-common.patients.calc[,list(COMMON.MEASURE=sum(Length)/sample.genes.length), by="Tumor_Sample_Barcode"]
    common.patients.calc<-as.vector(common.patients.calc$COMMON.MEASURE)
    common.patients.calc<-c(common.patients.calc, rep(0, total.patients.count-length(common.patients.calc))) #To account for those that did not have sample genes
    
    #Return vector
    return(common.patients.calc)
  }
  
  #Get number of genes per metabolite and metabolic.lenght [KEGG_ID, N.METABOLIC.GENES, METABOLITE.U]
  table.2.count<-table.2[,list(N.METABOLIC.GENES=length(Hugo_Symbol)), by="KEGG_ID"]
  table.2.sig<-as.data.table(merge(as.data.frame(metabolites.u), as.data.frame(table.2.count), by="KEGG_ID"))
  
  #Sampling pool will be from all background.genes we can get the length of (sequenced.background.genes)
  background.genes.sig<-background.sequenced.genes[background.sequenced.genes %in% length.table$Hugo_Symbol]
  table.1.sig<-as.data.table(merge(as.data.frame(table.1), as.data.frame(length.table), by="Hugo_Symbol"))
  table.1.sig<-table.1.sig[Tumor_Sample_Barcode %in% patients.t.metabolites$Tumor_Sample_Barcode,] #Use only patients analyzed to calculate u(j)
  
  if (method=="mean"){
    table.2.sig<-table.2.sig[,
                             list(P.VAL=sum(replicate(100,mean(Function.u.significance(N.METABOLIC.GENES, background.genes.sig, table.1.sig, length.table,KEGG_ID)))>=METABOLITE.U)/100), 
                             by="KEGG_ID"]
  } else if (method=="median") {
    table.2.sig<-table.2.sig[,
                             list(P.VAL=sum(replicate(100,median(Function.u.significance(N.METABOLIC.GENES, background.genes.sig, table.1.sig, length.table)))>=METABOLITE.U)/100), 
                             by="KEGG_ID"]
  }
  
  #Assign significance table
  metabolites.u<-as.data.table(merge(as.data.frame(metabolites.u), as.data.frame(table.2.sig), by="KEGG_ID"))
  metabolites.u$P.VAL.ADJ<-p.adjust(metabolites.u$P.VAL, method="fdr")
  metabolites.u<-metabolites.u[order(P.VAL.ADJ, decreasing=T),]
  
  #Return
  return(metabolites.u)
}

dummy.u<-Function.u(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, method="mean",remove.outliers=T)
dummy.u

#######061814#######
#Test on "Silent" mutations to prove that test is significant
dummy.1.bg<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061414_Table1_BRCA.bg.rds")

dummy.u.bg.bg<-Function.u.bg(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes,remove.outliers=T)
dummy.u.bg.ratio<-Function.u.ratio(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes,remove.outliers=T)
dummy.u.bg.u<-Function.u(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes, method="mean",remove.outliers=T)

dummy.u.bg.bg.types<-Function.u.bg.per.subtype(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes, dummy.clinical$C.TABLE,remove.outliers=T)
dummy.u.bg.ratio.types<-Function.u.ratio.per.subtype(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes, dummy.clinical$C.TABLE,remove.outliers=T)

#Per subtype plot function
Function.per.subtype.plot<-function(function.u.subtype.result, threshold, test.columns=c(1,3), p.val.column=NULL, table.2) {
  #Plots heatmap from result of u(j) function (list)
  #Asks for "KEGG_ID" and "P.VALUE" columns (test.columns)
  #Threshold is the significance level to be considered
  #Table.2 comes straight from Function.post.process.table.2.R()
  #Assumes 8 subtypes (MAY CHANGE THIS LATER ON)
  
  require(pheatmap)
  
  #x[as.vector(x[,p.val.column,with=F]<threshold), test.columns , with=F]
  
  test.clinical.plot<-lapply(function.u.subtype.result, function(x) as.data.frame(x[as.vector(x[,p.val.column,with=F]<threshold), test.columns , with=F]))
  test.clinical.plot.types<-names(test.clinical.plot)
  test.clinical.plot<-Reduce(function(x,y) merge(x,y, all=T, by="KEGG_ID"), test.clinical.plot )
  colnames(test.clinical.plot)<-c("KEGG_ID", test.clinical.plot.types)
  test.clinical.plot[,2:9]<--log(test.clinical.plot[,2:9]) #Based on number of subtypes
  test.clinical.plot[is.na(test.clinical.plot)]<-0
  
  #Get KEGG Metabolite names from table.2
  test.clinical.plot<-merge(test.clinical.plot, unique(as.data.frame(table.2[,1:2,with=F])), by="KEGG_ID")
  rownames(test.clinical.plot)<-test.clinical.plot$METABOLITE
  test.clinical.plot$KEGG_ID<-NULL
  test.clinical.plot$METABOLITE<-NULL
  
  #Plot heatmap
  PLOT<-pheatmap(as.matrix(test.clinical.plot),scale="none",fontsize=22)
  
  #Return
  return(PLOT)
}

Function.per.subtype.plot(dummy.u.ratio.subtype, 0.05, c(1,5), 5, dummy.2)

########061914######

#Try background mutation rate per patient and do t.test

Function.u.t.bmr<-function(table.1, table.2, length.table, background.sequenced.genes, remove.outliers=F) {
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

dummy.u.t<-Function.u.t.bmr(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
dummy.u.t$STARTING.PATIENTS
dummy.u.t$PATIENTS.IN.MODEL
length(dummy.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED)
dummy.u.t$STATS[WILCOXON.P.ADJ<0.05,] #3 we wanted found

dummy.bg.t<-Function.u.t.bmr(dummy.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes)
dummy.bg.t$STATS[order(WILCOXON.P.ADJ),] #nothing found, as desired

Function.u.t.bmr.per.subtype<-function (table.1, table.2, length.table, background.sequenced.genes, clinical.data.table, remove.outliers=F) {
  #Does Function.u.t.bmr() calculation for each subtype found in clinical.data.table
  #Clinical.data.table is obtained from Function.BRCA.clinical.type() function
  #   We assume at the moment that BRCA is the only cancer that has clinical subtype documented in TCGA
  
  require(data.table)
  
  #Prep clinical table - Only care about bcr_patient_code and subtype columnds
  clinical.data.table<-clinical.data.table[,c(1,5),with=F]
  
  #Classify and filter patients by subtype 
  table.1$bcr_patient_barcode<-substr(table.1$Tumor_Sample_Barcode,1,12)
  table.1<-merge(as.data.frame(table.1), as.data.frame(clinical.data.table), by="bcr_patient_barcode")
  
  #Apply Function.u.bg()
  table.1.subtypes.list<-split(table.1, table.1$subtype) #list of data.frames
  table.1.u.by.subtype<-lapply(table.1.subtypes.list, function(x) 
    Function.u.t.bmr(as.data.table(x[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "N.MUTATIONS")]), 
                     table.2, length.table, background.sequenced.genes, remove.outliers)
  )
  
  #Return
  return(table.1.u.by.subtype)
}

dummy.u.t.per.subtypes<-Function.u.t.bmr.per.subtype(dummy.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, dummy.clinical$C.TABLE,remove.outliers=F)
dummy.u.t.per.subtypes
Function.per.subtype.plot(dummy.u.t.per.subtypes, 0.05, c(1,13),13, dummy.2) #Only one subtype

########062014#######
#Find significant metabolites by u(j) in other cancers

#OV
OV.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_OV.rds")
OV.u.t<-Function.u.t.bmr(OV.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
OV.u.t$STATS #Nothing in OV (when looking at wilcoxon)
OV.u.t$STATS[PAIRED.P.ADJ<0.05,] #2 found

#GBM
GBM.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_GBM.rds")
GBM.u.t<-Function.u.t.bmr(GBM.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
GBM.u.t$STATS # Same 3 as in BRCA

GBM.1.bg<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_GBM.bg.rds")
GBM.bg.u.t<-Function.u.t.bmr(GBM.1.bg$table.1, dummy.2, dummy.table.length, dummy.1.bg$background.genes, remove.outliers=F)
GBM.bg.u.t #Nothing in the GBM background

#KIRC
KIRC.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_KIRC.rds")
KIRC.u.t<-Function.u.t.bmr(KIRC.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
KIRC.u.t #Nothing in KIRC

#LUAD
LUAD.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_LUAD.rds")
LUAD.u.t<-Function.u.t.bmr(LUAD.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
LUAD.u.t$STATS #One found which has been suggested previously in LUNG but not BREAST

LUAD.1.bg<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_LUAD.bg.rds")
LUAD.bg.u.t<-Function.u.t.bmr(LUAD.1.bg$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
LUAD.bg.u.t#NONE in background

#LUSC
LUSC.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_LUSC.rds")
LUSC.u.t<-Function.u.t.bmr(LUSC.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
LUSC.u.t$STATS #Found 2, same as LUAD plus cAMP

LUSC.1.bg<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_LUSC.bg.rds")
LUSC.bg.u.t<-Function.u.t.bmr(LUSC.1.bg$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
LUSC.bg.u.t#NONE in background

#COAD
COAD.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062014_Table1_COAD.rds")
COAD.u.t<-Function.u.t.bmr(COAD.1$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
COAD.u.t$STATS #Nothing on COAD (when looking at wilcoxon)
COAD.u.t$STATS[PAIRED.P.ADJ<0.05,] #One found

#LGG
LGG.1<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062114_Table1_LGG.rds")
LGG.u.t<-Function.u.t.bmr(LGG.1$table.1, dummy.2, dummy.table.length, LGG.1$background.genes, remove.outliers=F)
head(LGG.u.t$STATS[,c("KEGG_ID", "WILCOXON.P.ADJ"), with=F],10) #Found 6, specific to LGG (Brain Lower Grade Glioma)

LGG.1.bg<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/062114_Table1_LGG.bg.rds")
LGG.bg.u.t<-Function.u.t.bmr(LGG.1.bg$table.1, dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=T)
LGG.bg.u.t$STATS #4 in background

########062314#########
#Classify BRCA stages
BRCA.CLINICAL<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/CLINICAL/BRCA/060914_BRCA_CLINICAL/Clinical/Biotab//nationwidechildrens.org_clinical_patient_brca.txt",
                        header=T, sep="\t",skip=1))
BRCA.CLINICAL<-BRCA.CLINICAL[-1,]

ggplot(BRCA.CLINICAL, aes(x=factor(number_of_lymphnodes_positive_by_he))) + geom_bar()
ggplot(BRCA.CLINICAL, aes(x=pathologic_T)) + geom_bar()
ggplot(BRCA.CLINICAL, aes(x=pathologic_N)) + geom_bar()
ggplot(BRCA.CLINICAL, aes(x=pathologic_M)) + geom_bar()
ggplot(BRCA.CLINICAL, aes(x=pathologic_stage)) + geom_bar() #May not be necessarily helpful unless we know what we are looking for

#######062414##########
re.dummy<-Function.u.t.bmr(dummy.1$table.1[!(Tumor_Sample_Barcode %in%  dummy.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED),],
                           dummy.2, dummy.table.length, dummy.1$background.genes, remove.outliers=F)
re.dummy$STATS #Nothing found on left over co-hort
#Find what metabolites (re-calculate u(j)) using cohort not covered here (NOTHING FOUND)
#Are these phenotypically different than the rest of breast cancer patients? (NOT WORTH IT SINCE NOTHING FOUND)

#Explore Budczies meatabolites
BUD.METABOLITES<-as.data.table(read.csv("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/062414_BUDCZIES_METABOLITES", header=T, sep="\t"))

BUD.METABOLITES[KEGG_ID %in% dummy.u.t$STATS[WILCOXON.P.ADJ<0.05,]$KEGG_ID,] #BRCA significant u(j) not found here

#Explore significant u(j) genes
heatmap.u.cancers<-list(BRCA=as.data.frame(dummy.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13),with=F]),
                        OV=as.data.frame(OV.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        GBM=as.data.frame(GBM.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        KIRC=as.data.frame(KIRC.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        LUAD=as.data.frame(LUAD.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        LUSC=as.data.frame(LUSC.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        COAD=as.data.frame(COAD.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F]),
                        LGG=as.data.frame(LGG.u.t$STATS[WILCOXON.P.ADJ<0.05,c(1,13), with=F])
)
heatmap.u.cancers.names<-names(heatmap.u.cancers)
heatmap.u.cancers<-Reduce(function(x,y) merge(x,y, all=T, by="KEGG_ID",), heatmap.u.cancers)
colnames(heatmap.u.cancers)<-c("KEGG_ID", heatmap.u.cancers.names)
heatmap.u.cancers<-merge(heatmap.u.cancers, unique(as.data.frame(dummy.2[,1:2,with=F])), by="KEGG_ID",)
heatmap.u.cancers$KEGG_ID<-NULL
rownames(heatmap.u.cancers)<-heatmap.u.cancers$METABOLITE
heatmap.u.cancers$METABOLITE<-NULL
heatmap.u.cancers<--log(heatmap.u.cancers)
heatmap.u.cancers[is.na(heatmap.u.cancers)]<-0

library(RColorBrewer)
pheatmap(as.matrix(heatmap.u.cancers),scale="none",fontsize=22, col=colorRampPalette(c("white", "blue"))(100) )

cancers.u.patient.count<-data.frame(CANCERS=c("BRCA", "OV", "GBM","KIRC", "LUAD", "LUSC", "COAD", "LGG"),
                                    ALL.PATIENTS=c(dummy.u.t$STARTING.PATIENTS, OV.u.t$STARTING.PATIENTS, GBM.u.t$STARTING.PATIENTS, KIRC.u.t$STARTING.PATIENTS,
                                                   LUAD.u.t$STARTING.PATIENTS, LUSC.u.t$STARTING.PATIENTS, COAD.u.t$STARTING.PATIENTS, LGG.u.t$STARTING.PATIENTS),
                                    IN.MODEL=c(dummy.u.t$PATIENTS.IN.MODEL, OV.u.t$PATIENTS.IN.MODEL, GBM.u.t$PATIENTS.IN.MODEL, KIRC.u.t$PATIENTS.IN.MODEL,
                                               LUAD.u.t$PATIENTS.IN.MODEL, LUSC.u.t$PATIENTS.IN.MODEL, COAD.u.t$PATIENTS.IN.MODEL, LGG.u.t$PATIENTS.IN.MODEL),
                                    AFFECTED=c(length(dummy.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), length(OV.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), 
                                               length(GBM.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), length(KIRC.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED),
                                               length(LUAD.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), length(LUSC.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), 
                                               length(COAD.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED), length(LGG.u.t$PATIENT.AFFECTED$WILCOXON.AFFECTED))
                                                   )
ggplot(melt(cancers.u.patient.count), aes(x=CANCERS,y=value,fill=variable)) + geom_bar(position="dodge",stat="identity") +theme.format +
  ylab("PATIENTS")

########062514##########
ggplot(dummy.u.t$STATS, aes(wilcoxon.stat)) + geom_histogram() + theme.format
ggplot(dummy.u.t$STATS, aes(x=wilcoxon.stat, y=-log(wilcoxon.p.val))) + geom_point() + theme.format 

#######062714##########
dummy.3<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.3.rds")
dim(dummy.3)

b<-Function.read.RNAseq.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3",cancer.sep=T,
                              Function.process.RNAseq.map.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/FILE_SAMPLE_MAP.txt"))
dim(b$normal)
head(b$normal[,1:3])
log2(rowMeans(b$tumor["PAMR1",]) / rowMeans(b$normal["PAMR1",]))

Function.v.dysregulation<-function(normal.matrix, cancer.matrix, table.3, method="mean") {
  #Calculates v(j) based on average(or median, max) LogFC expression of enzymes that produce metabolite j
  #This differential expression is based on all cancer vs all normal
  #Calculates also empirical p-value given how many same number of enzymes would produce an average logFC expression level similar or greater than for j
  
  require(data.table)
  
  ##Prep tables
  setnames(table.3, colnames(table.3), c("Hugo_Symbol", "KEGG_ID"))
  
  ##Calculate mean(logFC) per metabolite
  
  #Do differential expression given normal and cancer matrices and obtain logFC table per gene [Hugo_Symbol, logFC]
  dummy.diff.exp<-Function.RNAseq.Differential.Expression(normal.matrix, cancer.matrix)
  dummy.diff.exp$Hugo_Symbol<-rownames(dummy.diff.exp)
  dummy.diff.exp<-dummy.diff.exp[,c("Hugo_Symbol", "logFC")]
  
  #Break logFC by metabolite [Hugo_Symbol, logFC, KEGG_ID]
  dummy.diff.exp.met<-as.data.table(merge(dummy.diff.exp, table.3, by="Hugo_Symbol"))
  
  #Get average logFC per metabolite [KEGG_ID, LOG.FC.MEAN] <- choice of "mean", "median", "max", "sum"
  if (method=="mean"){
    dummy.avg.logfc.met<-dummy.diff.exp.met[,list(LOG.FC.MEAN=mean(logFC)), by="KEGG_ID"]
  } else if (method=="median") {
    dummy.avg.logfc.met<-dummy.diff.exp.met[,list(LOG.FC.MEAN=median(logFC)), by="KEGG_ID"]
  } else if (method=="max") {
    dummy.avg.logfc.met<-dummy.diff.exp.met[,list(LOG.FC.MEAN=max(logFC)), by="KEGG_ID"]
  } else if (method=="sum") {
    dummy.avg.logfc.met<-dummy.diff.exp.met[,list(LOG.FC.MEAN=sum(logFC)), by="KEGG_ID"]
  }
  
  #Get count of enzymes per metabolite [KEGG_ID, ENZYME.COUNT]
  table.3.enzyme.count<-table.3[,list(ENZYME.COUNT=length(Hugo_Symbol)), by="KEGG_ID"]
  
  #Merge to get expected value (logfc mean) and count of enzymes that provided such value for each metabolite [KEGG_ID, LOG.FC.MEAN, ENZYME.COUNT]
  dummy.avg.logfc.met<-as.data.table(merge(as.data.frame(dummy.avg.logfc.met), as.data.frame(table.3.enzyme.count), by="KEGG_ID"))
  
  ##Do empirical test
  
  #Get bacground logFC population
  dummy.background.logfc<-as.vector(dummy.diff.exp$logFC)
  
  #Do test - choice of "mean", "median", "sum" and "max"
  #This is a two.tailed test, since we are not certain what to expect and value could be in either extreme (i.e. don't know if it is over- or under- expressed)
  if (method=="mean"){
    dummy.logfc.test<-dummy.avg.logfc.met[,list(P.VAL=sum(replicate(1000, abs(mean(sample(dummy.background.logfc, ENZYME.COUNT))))>=
                                                            abs(LOG.FC.MEAN))/1000), by="KEGG_ID"]
  } else if (method=="median"){
    dummy.logfc.test<-dummy.avg.logfc.met[,list(P.VAL=sum(replicate(1000, abs(median(sample(dummy.background.logfc, ENZYME.COUNT))))>=
                                                            abs(LOG.FC.MEAN))/1000), by="KEGG_ID"]
  } else if (method=="max"){
    dummy.logfc.test<-dummy.avg.logfc.met[,list(P.VAL=sum(replicate(1000, abs(max(sample(dummy.background.logfc, ENZYME.COUNT))))>=
                                                            abs(LOG.FC.MEAN))/1000), by="KEGG_ID"]
  } else if (method=="sum"){
    dummy.logfc.test<-dummy.avg.logfc.met[,list(P.VAL=sum(replicate(1000, abs(sum(sample(dummy.background.logfc, ENZYME.COUNT))))>=
                                                            abs(LOG.FC.MEAN))/1000), by="KEGG_ID"]
  }
  
  dummy.avg.logfc.met<-as.data.table(merge(as.data.frame(dummy.avg.logfc.met), as.data.frame(dummy.logfc.test), by="KEGG_ID"))
  
  #Do multiple hypothesis test
  dummy.avg.logfc.met$ADJ.P.VAL<-p.adjust(dummy.avg.logfc.met$P.VAL,method="fdr")
  
  #Clean up
  if (method=="median") {
    setnames(dummy.avg.logfc.met, colnames(dummy.avg.logfc.met), c("KEGG_ID", "LOG.FC.MEDIAN", colnames(dummy.avg.logfc.met)[3:5]))
  } else if (method=="max"){
    setnames(dummy.avg.logfc.met, colnames(dummy.avg.logfc.met), c("KEGG_ID", "LOG.FC.MAX", colnames(dummy.avg.logfc.met)[3:5]))
  } else if (method=="sum"){
    setnames(dummy.avg.logfc.met, colnames(dummy.avg.logfc.met), c("KEGG_ID", "LOG.FC.SUM", colnames(dummy.avg.logfc.met)[3:5]))
  }
  
  dummy.avg.logfc.met<-dummy.avg.logfc.met[order(ADJ.P.VAL),]
  
  #Return
  return(dummy.avg.logfc.met)
}

BRCA.v.dys<-Function.v.dysregulation(b$normal, b$tumor, dummy.3,method="sum")
BRCA.v.dys$BUDCZIES<-BRCA.v.dys$KEGG_ID %in% BUD.METABOLITES$KEGG_ID
BRCA.v.dys[KEGG_ID %in% BUD.METABOLITES$KEGG_ID, ]
ggplot(BRCA.v.dys, aes(ADJ.P.VAL)) +geom_histogram() + theme.format
ggplot(BRCA.v.dys, aes(x=ADJ.P.VAL, y=LOG.FC.SUM, colour=BUDCZIES)) + geom_point(size=6) + geom_vline(xintercept=0.05,colour="blue") +theme.format

Function.j.expression<-function(normal.matrix, cancer.matrix, table.3, method="mean"){
  #Performs differential expression on predicted expression values of metabolites
  #Obtains a topTable for all metabolites
  
  require(data.table)
  
  ##Prep tables
  setnames(table.3, colnames(table.3), c("Hugo_Symbol", "KEGG_ID"))
  
  ##Convert gene expression matrices to metabolic expression matrices
  normal.matrix$Hugo_Symbol<-rownames(normal.matrix)
  cancer.matrix$Hugo_Symbol<-rownames(cancer.matrix)
  
  
  dummy.normal.matrix<-as.data.table(merge(as.data.frame(normal.matrix) , as.data.frame(table.3), by="Hugo_Symbol"))
  dummy.cancer.matrix<-as.data.table(merge(as.data.frame(cancer.matrix) , as.data.frame(table.3), by="Hugo_Symbol"))
  
  dummy.normal.matrix<-dummy.normal.matrix[,-1,with=F] #Removes Hugo_Symbol column
  dummy.cancer.matrix<-dummy.cancer.matrix[,-1,with=F] #Removes Hugo_Symbol column
  
  dummy.normal.matrix$KEGG_ID<-as.character(dummy.normal.matrix$KEGG_ID)
  dummy.cancer.matrix$KEGG_ID<-as.character(dummy.cancer.matrix$KEGG_ID)
  
  dummy.normal.lists<-split(dummy.normal.matrix, dummy.normal.matrix$KEGG_ID) #Splits into each KEGG across patients
  dummy.cancer.lists<-split(dummy.cancer.matrix, dummy.cancer.matrix$KEGG_ID) #to perform colMeans per KEGG
  
  if (method=="mean"){
    dummy.normal.kegg.matrix<-sapply(dummy.normal.lists, function(x) colMeans(x[,!"KEGG_ID",with=F]) ,USE.NAMES=T ) #Drops KEGG column
    dummy.cancer.kegg.matrix<-sapply(dummy.cancer.lists, function(x) colMeans(x[,!"KEGG_ID",with=F]) ,USE.NAMES=T ) #and performs colMeans
  } else if (method=="median"){
    dummy.normal.kegg.matrix<-sapply(dummy.normal.lists, function(x) apply(x[,!"KEGG_ID",with=F],2,median) ,USE.NAMES=T ) #Drops KEGG column
    dummy.cancer.kegg.matrix<-sapply(dummy.cancer.lists, function(x) apply(x[,!"KEGG_ID",with=F],2,median) ,USE.NAMES=T ) #and performs colMedians
  } else if (method=="sum") {
    dummy.normal.kegg.matrix<-sapply(dummy.normal.lists, function(x) colSums(x[,!"KEGG_ID",with=F]) ,USE.NAMES=T ) #Drops KEGG column
    dummy.cancer.kegg.matrix<-sapply(dummy.cancer.lists, function(x) colSums(x[,!"KEGG_ID",with=F]) ,USE.NAMES=T ) #and performs colSums
  }
  
  dummy.normal.kegg.matrix<-as.data.frame(t(dummy.normal.kegg.matrix)) #Need to transpose
  dummy.cancer.kegg.matrix<-as.data.frame(t(dummy.cancer.kegg.matrix)) #because of sapply output
  
  ##Do KEGG metabolite differential expression
  dummy.kegg.diff.exp<-Function.RNAseq.Differential.Expression(dummy.normal.kegg.matrix, dummy.cancer.kegg.matrix)
  
  #Return
  return(dummy.kegg.diff.exp)  
}

BRCA.j.exp<-Function.j.expression(b$normal, b$tumor, dummy.3,method="median")
BRCA.j.exp$ADJ.P.VAL<-p.adjust(BRCA.j.exp$P.Value, method="fdr")
BRCA.j.exp$KEGG_ID<-rownames(BRCA.j.exp)
BRCA.j.exp$BUDCZIES<-BRCA.j.exp$KEGG_ID %in% BUD.METABOLITES$KEGG_ID
dim(BRCA.j.exp)
head(BRCA.j.exp)
hist(BRCA.j.exp$logFC)
dim(BRCA.j.exp[BRCA.j.exp$ADJ.P.VAL<0.05,])

#1446, 1414, 1371
dim(BRCA.j.exp[BRCA.j.exp$KEGG_ID %in% BUD.METABOLITES$KEGG_ID & BRCA.j.exp$adj.P.Val<0.05, ])
ggplot(BRCA.j.exp, aes(x=ADJ.P.VAL, y=logFC, colour=BUDCZIES)) + geom_point(size=3) + geom_vline(xintercept=0.05,colour="blue") + theme.format


test<-merge(as.data.frame(dummy.1$table.1[,1:2,with=F]), as.data.frame(dummy.3), by.x="Hugo_Symbol", by.y="Enzyme" )
dim(test)
head(test[,1:3])

#######062714#########

##Look at significant u(j) in BRCA
Function.u.patient.x.genes<-function(table.1, significant.keggs, table.2) {
  #Builds matrix of patients vs genes depending on fed u(j) metabolites
  #Returns list that includes matrix, significant u(j) associated genes and heatmap
  
  require(data.table)
  require(pheatmap)
  require(reshape)
  require(gplots)
  
  #Prep tables
  table.2$METABOLITE<-NULL
  setnames(table.2, colnames(table.2), c("KEGG_ID", "Hugo_Symbol"))
  
  #Get significant u(j) genes
  dummy.sig.genes<-unique(as.vector(table.2[KEGG_ID %in% significant.keggs, ]$Hugo_Symbol))
  
  #Pick patients with significant genes
  dummy.table.1<-table.1[Hugo_Symbol %in% dummy.sig.genes,]
  
  #Build Patient x gene matrix
  dummy.matrix.1<-cast(dummy.table.1, Tumor_Sample_Barcode~Hugo_Symbol, value="N.MUTATIONS", sum)
  rownames(dummy.matrix.1)<-dummy.matrix.1$Tumor_Sample_Barcode
  dummy.matrix.1$Tumor_Sample_Barcode<-NULL
  
  #Build plot
  PLOT.1<-heatmap.2(as.matrix(dummy.matrix.1),trace="none",scale="none", col=hmcol)
  
  #Return
  return(list(matrix=as.data.frame(dummy.matrix.1), heatmap=PLOT.1, sig.genes=dummy.sig.genes))
}

dummy.1$table.1
as.vector(dummy.u.t$STATS[WILCOXON.P.ADJ<0.05,]$KEGG_ID)
BRCA.u.patient.genes<-Function.u.patient.x.genes(dummy.1$table.1, as.vector(dummy.u.t$STATS[WILCOXON.P.ADJ<0.05,]$KEGG_ID), dummy.2) 
#Not much found with mutations

##Look at significant mutations u(j) compared to 1000G
THOUSAND.G<-as.data.table(read.csv("DATABASES/CANCER_DATA/1000GENOME/062714_1000G_EXON", header=T, sep="\t"))

rect<-data.frame(xmin=5.0e+07,xmax=1.0e+08,ymin=-Inf,ymax=Inf)
ggplot(THOUSAND.G[B.Chrom %in% 1:4,], aes(x=B.Position, y=ALT.ALLELES, colour=Variant_Class)) + geom_point() + theme.format +
  facet_wrap(~B.Chrom,ncol=1) +
  geom_rect(data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), alpha=0.5, color="grey20", inherit.aes=F)


Function.maf.maps<-function(main.maf, extra.mafs=c()){
  #Process maf files to return mutation count and location
  #This contains silent mutations
  
  require(data.table)
  
  #Process main.file
  dummy.a<-read.csv(main.maf, header=T,  sep="\t")
  dummy.a$Line_Number<-NULL #To remove duplicates
  dummy.a<-unique(dummy.a)
  
  #Filter for info we want
  dummy.a<-dummy.a[,c("Hugo_Symbol", "Chrom", "Start_Position", "Reference_Allele",
                      "Tumor_Seq_Allele2","Variant_Classification","Tumor_Sample_Barcode")]
  setnames(dummy.a, colnames(dummy.a), c(colnames(dummy.a)[1:2], "Position",  colnames(dummy.a)[4:5],
                                         "Variant_Class", "Tumor_Sample_Barcode"))
  dummy.a<-as.data.table(dummy.a)
  
  #Process additional files
  for (maf in extra.mafs) {
    a.1<-read.csv(maf, header=T, sep="\t")
    a.1$Line.Number<-NULL
    a.1<-unique(a.1)
    
    a.1<-a.1[,c("Hugo_Symbol", "Chrom", "Start_Position", "Reference_Allele",
                "Tumor_Seq_Allele2","Variant_Classification","Tumor_Sample_Barcode")]
    setnames(a.1, colnames(a.1), c(colnames(a.1)[1:2], "Position",  colnames(a.1)[4:5],
                                   "Variant_Class", "Tumor_Sample_Barcode"))
    a.1<-as.data.table(a.1)
    
    #Add to main file
    dummy.a<-rbind(dummy.a,a.1)
  }
  
  #Count total number of samples
  total.samples<-length(unique(dummy.a$Tumor_Sample_Barcode))
  
  #Count of all specific types of mutations
  dummy.factors<-colnames(dummy.a)[1:6]
  dummy.count<-dummy.a[, list(ALT.MUTATIONS=length(Tumor_Sample_Barcode)), by=dummy.factors]
  dummy.count$TOTAL.SAMPLES<-total.samples
  
  #Return
  return(dummy.count)
}

BRCA.maf.maps<-Function.maf.maps("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/061014/Somatic_Mutations/WUSM__IlluminaGA_DNASeq_curated/Level_2/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", c("DATABASES/CANCER_DATA/TCGA/SOMATIC_MUTATIONS/BRCA/061014/Somatic_Mutations/WUSM__IlluminaGA_DNASeq/Level_2/genome.wustl.edu__Illumina_Genome_Analyzer_DNA_Sequencing_level2.maf"))

table(BRCA.maf.maps$ALT.MUTATIONS)

########060214#######
#Function to look at 1000G frequencies vs u(j) gene mutation frequencies over gene coordinates

BRCA.u.patient.genes$sig.genes
THOUSAND.G
BRCA.maf.maps

#Load gene mRNA coordinates
mRNA.coor<-as.data.table(read.csv("DATABASES/REFSEQ/070114_GENE_COORDIANTES",sep="\t",header=T))

Function.gene.mut.analysis<-function (gene, thousand.table, maf.table, mrna.coordinates) {
  
  require(ggplot2)
  require(data.table)
  
  #Filter all tables for gene of interest
  thousand.table<-thousand.table[Hugo_Symbol==gene,]
  maf.table<-maf.table[Hugo_Symbol==gene,]
  mrna.coordinates<-mrna.coordinates[Hugo_Symbol==gene,]
  
  #Filter for lines that we are interested in
  thousand.table<-thousand.table[,c("Hugo_Symbol", "B.Chrom", "B.Position", "TOTAL.ALLELES", "ALT.ALLELES","Variant_Class"),with=F]
  maf.table<-maf.table[,c("Hugo_Symbol","Chrom","Position", "TOTAL.SAMPLES","ALT.MUTATIONS", "Variant_Class"), with=F]
  
  #Classify based only on "synonymous" or "nonsynonymous"
  maf.table$Variant_Class<-as.character(maf.table$Variant_Class)
  maf.table$Variant_Class[maf.table$Variant_Class=="Silent"]<-"synonymous"
  maf.table$Variant_Class[maf.table$Variant_Class!="synonymous"]<-"nonsynonymous"
  thousand.table$Variant_Class[thousand.table$Variant_Class!="synonymous"]<-"nonsynonymous"
  
  #Join 1000G and cancer.maf information
  setnames(thousand.table, colnames(thousand.table), colnames(maf.table))
  thousand.table$TYPE<-"1000G"
  maf.table$TYPE<-"CANCER"
  dummy.table<-rbind(maf.table, thousand.table)
  
  #Get gene coordinates to plot in graph
  START=unique(as.vector(mrna.coordinates$START))
  END=unique(as.vector(mrna.coordinates$END))
  
  #Initialize plot
  dummy.plot<-ggplot(dummy.table, aes(x=Position, y=ALT.MUTATIONS/TOTAL.SAMPLES, colour=TYPE)) + geom_point(size=3) +
    facet_wrap(~Variant_Class, ncol=1) + xlim(START, END) +
    theme(strip.text.x = element_text(size = 24))
  
  for (line in 1:nrow(mrna.coordinates)){
    EXON_START=as.vector(mrna.coordinates[line,]$FEAT_START)
    EXON_END=as.vector(mrna.coordinates[line,]$FEAT_END)
    
    rect.frame<-data.frame(xmin=EXON_START, xmax=EXON_END, ymin=-Inf, ymax=Inf)
    
    dummy.plot<-dummy.plot + geom_rect(data=rect.frame, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.3, color="grey20", inherit.aes=F)    
  }
  
  #Return
  return(dummy.plot)
}

Function.gene.mut.analysis("PLD1", THOUSAND.G, BRCA.maf.maps, mRNA.coor) + theme.format 

#######070714########
#Function to look at 1000G frequencies vs u(j) gene mutation frequencies over amino acid mutations in pfams
THOUSAND.G.AA<-fread("DATABASES/CANCER_DATA/1000GENOME/070714_1000G_AA", header=T, sep="\t")
PFAM.AA<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/070714_HUGO_PFAM.rds")

Function.post.TCGA.Dario<-function(TCGA.Dario.files=c()){
  #Post processes Dario's output of TCGA mutations to obtain frequencies per gene across the specific cancer population
  #Keep in mind that this would only account for mutation in the first isoform .001
  
  require(data.table)
  require(plyr)
  
  #Open files
  dummy.files<-lapply(TCGA.Dario.files, function(x) read.csv(x, sep="\t",header=F, stringsAsFactors=F))
  
  #Combine tables and eliminate replicates
  dummy.combined<-Reduce(function(x,y) rbind.fill(x,y), dummy.files)
  dummy.combined<-unique(dummy.combined)
  dummy.combined<-as.data.table(dummy.combined)
  
  #Keep only first isoform
  setnames(dummy.combined, colnames(dummy.combined), c("Isoform", "Patient", "Position", "REF.AA", "ALT.AA"))
  dummy.combined$First<-grepl(".001", dummy.combined$Isoform)
  dummy.combined<-dummy.combined[First==TRUE,]
  dummy.combined$First<-NULL
  dummy.combined$Hugo_Symbol<-sapply(dummy.combined$Isoform, function(y) strsplit(y, ".001")[[1]][1])
  dummy.combined$Isoform<-NULL
  
  #Classify mutations (syn vs nonsyn) [Hugo_Symbol, Patient, Position, REF.AA, ALT.AA, Variant_Class]
  dummy.combined$Variant_Class<-dummy.combined$REF.AA==dummy.combined$ALT.AA
  dummy.combined$Variant_Class[dummy.combined$Variant_Class==TRUE]<-"synonymous"
  dummy.combined$Variant_Class[dummy.combined$Variant_Class==FALSE]<-"nonsynonymous"
  
  #Get total patient count 
  TOTAL.PATIENT.COUNT<-length(unique(as.vector(dummy.combined$Patient)))
  
  #Get patient counts per mutation [Hugo_Symbol, Position, REF.AA, ALT.AA, Variant_Class, ALT.ALLELES, TOTAL.ALLELES]
  dummy.combined<-dummy.combined[, list(ALT.ALLELES=length(Patient)), by=c("Hugo_Symbol", "Position", "REF.AA", "ALT.AA", "Variant_Class")]
  dummy.combined$TOTAL.ALLELES<-TOTAL.PATIENT.COUNT
  
  #Return table
  return(dummy.combined)
}

TCGA.BRCA.AA.FREQ<-Function.post.TCGA.Dario(c("DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/BRCA/070714_BRCA_CURATED_NON_SYN.mut",
                                              "DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/BRCA/070714_BRCA_CURATED_SYN.mut",
                                              "DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/BRCA/070714_BRCA_NON_SYN.mut",
                                              "DATABASES/CANCER_DATA/TCGA/PROCESSED_MUTATIONS/BRCA/070714_BRCA_SYN.mut"))

Function.AA.mut.analysis<-function(gene, thousand.aa.table, tcga.aa.table, pfam.table) {
  #Plots mutations frequency of gene of interest at the amino acid comparing 1000G vs TCGA on PFAMs
  
  require(data.table)
  require(ggplot2)
  
  #Filter all tables for gene of interest
  thousand.aa.table<-thousand.aa.table[Hugo_Symbol==gene,]
  tcga.aa.table<-tcga.aa.table[Hugo_Symbol==gene,]
  pfam.table<-pfam.table[Hugo_Symbol==gene,]
  
  #Classify based only on "synonymous" or "nonsynonymous"
  thousand.aa.table$Variant_Class[thousand.aa.table$Variant_Class!="synonymous"]<-"nonsynonymous"
  
  #Join 1000G and cancer.maf information
  thousand.aa.table<-thousand.aa.table[,c(1,3:5,8,7,6), with=F]
  setnames(thousand.aa.table, colnames(thousand.aa.table), colnames(tcga.aa.table))
  thousand.aa.table$TYPE<-"1000G"
  tcga.aa.table$TYPE<-"CANCER"
  dummy.table<-rbind(tcga.aa.table, thousand.aa.table)
  
  #Get protein limits to plot in graph
  START=unique(as.vector(pfam.table$START))
  END=unique(as.vector(pfam.table$END))
  
  #Initialize plot
  dummy.table$Position<-as.numeric(dummy.table$Position)
  dummy.plot<-ggplot(dummy.table, aes(x=Position, y=ALT.ALLELES/TOTAL.ALLELES, colour=TYPE)) + geom_point(size=3) +
    facet_wrap(~Variant_Class, ncol=1) + xlim(START, END) +
    theme(strip.text.x = element_text(size = 24))
  
  for (line in 1:nrow(pfam.table)){
    PFAM_START=as.vector(pfam.table[line,]$FEAT.START)
    PFAM_END=as.vector(pfam.table[line,]$FEAT.END)
    
    rect.frame<-data.frame(xmin=PFAM_START, xmax=PFAM_END, ymin=-Inf, ymax=Inf)
    
    dummy.plot<-dummy.plot + geom_rect(data=rect.frame, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha=0.3, color="grey20", inherit.aes=F)    
  }
  
  #Return
  return(dummy.plot)
  
}

Function.AA.mut.analysis("ABCD4", THOUSAND.G.AA, TCGA.BRCA.AA.FREQ, PFAM.AA) + theme.format

########070914########
#Function to calculate enrichment of mutations in PFAMs per gene
length(BRCA.u.patient.genes$sig.genes)
THOUSAND.G.AA
TCGA.BRCA.AA.FREQ
PFAM.AA

Function.pfam.mut.enrichment<-function(gene.list, thousand.aa.table, tcga.aa.table, pfam.table) {
  #Calculates enrichments of tcga mutations vs 1000G mutations of selected genes in PFMA domains
  #This calculation is done by a hypergeometric calculation on each genes, this does not take into account the number of samples covered by mutation
  
  require(data.table)
  
  #Filter tables for gene of interest
  thousand.aa.table<-thousand.aa.table[Hugo_Symbol %in% gene.list,]
  tcga.aa.table<-tcga.aa.table[Hugo_Symbol %in% gene.list,]
  pfam.table<-pfam.table[Hugo_Symbol %in% gene.list,]
  
  #Prep table
  #Classify based only on "synonymous" or "nonsynonymous"
  thousand.aa.table$Variant_Class[thousand.aa.table$Variant_Class!="synonymous"]<-"nonsynonymous"
  
  #First calculate for tcga genes [Hugo_Symbol, Variant_Class, PFAM, P.VALUE]
  dummy.tcga<-as.data.table(merge(as.data.frame(tcga.aa.table), as.data.frame(pfam.table), by="Hugo_Symbol"))
  dummy.tcga$Position<-as.numeric(dummy.tcga$Position)
  dummy.tcga$IN.PFAM<-dummy.tcga$Position >dummy.tcga$FEAT.START & dummy.tcga$Position <dummy.tcga$FEAT.END
  dummy.tcga<-dummy.tcga[,list(P.VAL=phyper(sum(IN.PFAM)-1,
                                            unique(as.vector(FEAT.END))-unique(as.vector(FEAT.START)),
                                            unique(as.vector(END))-(unique(as.vector(FEAT.END))-unique(as.vector(FEAT.START))),
                                            length(IN.PFAM), lower.tail=F)) ,
                         by=c("Hugo_Symbol", "Variant_Class", "PFAM")]
  dummy.tcga$TYPE<-"CANCER"
  
  #Then for 1000G [Hugo_Symbol, Variant_Class, PFAM, P.VALUE]
  dummy.thousand<-as.data.table(merge(as.data.frame(thousand.aa.table), as.data.frame(pfam.table), by="Hugo_Symbol"))
  dummy.thousand$B.Position<-as.numeric(dummy.thousand$B.Position)
  dummy.thousand$IN.PFAM<-dummy.thousand$B.Position>dummy.thousand$FEAT.START & dummy.thousand$B.Position<dummy.thousand$FEAT.END
  dummy.thousand<-dummy.thousand[,list(P.VAL=phyper(sum(IN.PFAM)-1,
                                              unique(as.vector(FEAT.END))-unique(as.vector(FEAT.START)),
                                              unique(as.vector(END))-(unique(as.vector(FEAT.END))-unique(as.vector(FEAT.START))),
                                              length(IN.PFAM), lower.tail=F)) ,
                         by=c("Hugo_Symbol", "Variant_Class", "PFAM")]
  dummy.thousand$TYPE<-"1000G"
  
  #Merge into single table
  dummy.result<-rbind(dummy.tcga, dummy.thousand)
  
  #Multiple hypothesis correction
  dummy.result$P.ADJ<-p.adjust(dummy.result$P.VAL,method="fdr")
  
  #Clean up and return
  dummy.result<-dummy.result[order(P.ADJ),]
  return(dummy.result)
}

BRCA.PFAM.u<-Function.pfam.mut.enrichment(BRCA.u.patient.genes$sig.genes, THOUSAND.G.AA, TCGA.BRCA.AA.FREQ, PFAM.AA)
head(BRCA.PFAM.u,20 ) #Nothing significant found

ggplot(unique(BRCA.PFAM.u[,c(1:2,4:6),with=F]), aes(x=P.VAL, fill=TYPE)) + 
  geom_histogram(position="dodge") +theme.format + 
  facet_wrap(~Variant_Class,ncol=1) + geom_vline(xintercept=0.05,colour="blue") +
  theme(strip.text.x = element_text(size = 20))

#Function to calculate enrichment based on conservation score - LEAVE FOR LATER

########071114########
#Integrate network
dummy.u.t$STATS #Function.u.t.bmr

Function.v<-function(RNA.SEQ.FOLDER, cancer.map.file, processed.table.1.rds, processed.table.3.rds, Table.u.KEGG_IDs) {
  #Calculates v(j) based on degree of influence of KEGG(j) mutations on cancer phenotype when compared to normal
  #Need functions pre-loaded:
  #   Function.read.RNAseq.files()
  #   Function.process.RNAseq.map.files()
  require(base)
  require(data.table)
  
  #Load RDS files
  dummy.table.1<-readRDS(processed.table.1.rds)
  dummy.table.1<-dummy.table.1$table.1[,1:2, with=F]
  dummy.table.3<-readRDS(processed.table.3.rds)
  
  #Create tumor and cancer matrices using map file on RNA-seq folder
  RNASEQ.MATRICES<-Function.read.RNAseq.files(RNA.SEQ.FOLDER, cancer.sep=T,
                                              Function.process.RNAseq.map.files(cancer.map.file))
  
  #Get cancer samples in processed.table.1 that have expression information in cancer matrix [Hugo_Symbol, PATIENT]
  cancer.rnaseq.patients<-colnames(RNASEQ.MATRICES$tumor)
  dummy.table.1$PATIENT<-sapply(as.character(dummy.table.1$Tumor_Sample_Barcode), function(x) paste(strsplit(x, "-")[[1]][1:4] , collapse="."))
  dummy.table.1<-dummy.table.1[PATIENT %in% cancer.rnaseq.patients, ]
  dummy.table.1$Tumor_Sample_Barcode<-NULL
  
  #Prep storeage matrix for v(j)
  dummy.table.v<-data.frame(KEGG_ID=c(), v.METABOLITE=c())
  
  #To keep count
  x=0
  
  #Process each metabolite for v(j)
  for (metabolite in Table.u.KEGG_IDs) {
    
    #Get enzymes that produce metabolite
    dummy.metabolite.enzymes<-unique(as.vector(dummy.table.3[KEGG_ID==metabolite,]$Enzyme))
    
    #Obtain cancer samples in dummy.table.1 that contain mutations in metabolite's enzymes
    dummy.G1.samples<-unique(as.vector(dummy.table.1[Hugo_Symbol %in% dummy.metabolite.enzymes, ]$PATIENT))
    
    #Account for fact that there may not be more than 1 sample with mutation, cannot do diff expression in that case
    if (length(dummy.G1.samples)>1) {
      
      #Get expression matrix for patients affected with mutations associated with metabolic enzymes
      dummy.G1.expression<-RNASEQ.MATRICES$tumor[,dummy.G1.samples]
      
      #Do differential expression against normal
      dummy.diff.exp<-Function.RNAseq.Differential.Expression(RNASEQ.MATRICES$normal, dummy.G1.expression)
      
      #Get v(j) by dividing that were differentially expressed over all genes tested
      dummy.G.diff.exp.genes<-nrow(dummy.diff.exp[dummy.diff.exp$adj.P.Val<0.05,]) #Value set at p<0.05 for bonferroni corrected p-values
      dummy.v.metabolite<-dummy.G.diff.exp.genes/nrow(dummy.diff.exp)
      
    } else #If 1 or no sample have mutation, then its influence will be equal to zero
      dummy.v.metabolite<-0
    
    #Add v(metabolite) to table
    dummy.table.v<-rbind(dummy.table.v, data.frame(KEGG_ID=metabolite, v.METABOLITE=dummy.v.metabolite))
    
    #To count
    x=x+1
    print (x/length(Table.u.KEGG_IDs))
  }
  #Return v(j) table
  dummy.table.v<-as.data.table(dummy.table.v)
  return(dummy.table.v)
}

Table.v.BRCA<-Function.v("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3", 
                         "DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/FILE_SAMPLE_MAP.txt",
                         "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds",
                         "PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.3.rds",
                         as.vector(dummy.u.t$STATS$KEGG_ID))
saveRDS(Table.v.BRCA, file="PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071214.BRCA.Table.v.rds")
ggplot(Table.v.BRCA, aes(v.METABOLITE)) + geom_histogram() + theme.format

dummy.u.t$STATS
Table.v.BRCA

#Is there a correlation between u(j) and v(j)

test<-as.data.table(merge(as.data.frame(dummy.u.t$STATS[,c("KEGG_ID", "wilcoxon.stat", "WILCOXON.P.ADJ"), with=F]), as.data.frame(Table.v.BRCA), by="KEGG_ID" ))
test$wilcoxon.stat.norm<-normalize.vector(test$wilcoxon.stat)
test$u.sig<-test$WILCOXON.P.ADJ<0.05
ggplot(test, aes(x=v.METABOLITE, y=wilcoxon.stat.norm)) + geom_point(aes(colour=u.sig)) +theme.format
cor.test(test$v.METABOLITE, test$wilcoxon.stat.norm, method="spearman") #YES, STRONG

########071214##########
#Significant test for v(j)
#For number of samples that produced v(j) for j, randomly sample the same number from the cancer distribution and calculate how many would produce the same or better ratio of differentially expressed genes by chance v(j)

BRCA.CANCER.MATRICES<-Function.read.RNAseq.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3",cancer.sep=T,
                              Function.process.RNAseq.map.files("DATABASES/CANCER_DATA/TCGA/RNA_SEQ/BRCA/060614_BRCA_RNASEQ/FILE_SAMPLE_MAP.txt"))
saveRDS(BRCA.CANCER.MATRICES, file="PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES.rds")

#Function.v.p.value needs to be modified and broken so that it takes *.v.pval.null.rds object - DONE
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

Table.v.BRCA.PVAL<-Function.v.p.value("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds",
                                      "PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.3.rds",
                                      "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071214.Table.v.BRCA.rds",
                                      "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.v.pval.null.rds")
dim(Table.v.BRCA.PVAL)
head(Table.v.BRCA.PVAL,20)

#Test interaction again
test<-as.data.table(merge(as.data.frame(dummy.u.t$STATS[,c("KEGG_ID", "wilcoxon.stat", "WILCOXON.P.ADJ"), with=F]), as.data.frame(Table.v.BRCA.PVAL), by="KEGG_ID" ))
test$wilcoxon.stat.norm<-normalize.vector(test$wilcoxon.stat)
test$u.sig<-test$WILCOXON.P.ADJ<0.05
test$v.sig<-test$ADJ.P.VAL<0.05
test$sig<-ifelse(test$u.sig==T & test$v.sig==T, "Both",
                 ifelse(test$u.sig==T & test$v.sig==F, "u.Sig",
                        ifelse(test$u.sig==F & test$v.sig==T, "v.Sig",
                               "None")))
test$sig.overall<-ifelse(test$sig!="None","Significant", 
                         "Not Significant")
  
test.sig<-test[sig!="None",][order(wilcoxon.stat.norm),]
cat(as.vector(test[sig!="None",][order(wilcoxon.stat.norm),]$KEGG_ID))

library(stats)
ggplot(test, aes(x=v.METABOLITE, y=wilcoxon.stat.norm)) + theme.format + geom_point(size=1.7)+
  geom_point(aes(colour=sig, size=sig.overall))+
  scale_color_manual(values=c("red","green4","yellow1","black"))  +
  geom_smooth(data=subset(test, sig.overall=="Significant"), method="lm",se=F, formula=y~exp(x*11)) +
  scale_colour_discrete(name="Significance\nType") + scale_size_discrete(name="Significance\nStatus") +
  ylab(label="u(j)")+xlab(label="v(j)") + 
  ggtitle("Classification of Driver and Passanger Metabolites\nby PMN Metabolic Features") + 
  theme(plot.title = element_text(lineheight=1, face="bold",size=24),
        legend.title=element_text(size=14))

cor.test(test.sig$v.METABOLITE, test.sig$wilcoxon.stat.norm, method="spearman") #YES, STRONG


########071414######
#TO DO:
#   Modify Function.v() so that it incorporates null file to add p.val
#   Write official null.pval and Function.v() scripts. State that null.pval has to be done on a powerful machine with parallelization=T

#Construct network for driver metabolites based on significant metabolites (obtained from significance table i.e., test.sig)
dim(test.sig)
test[v.METABOLITE==0,]$KEGG_ID
hist(test.sig$v.METABOLITE)

#Calculate V - WRITE APPROPRITE SCRIPT !!!!
Function.V.construct<-function(feature.sig.table, beta.range) {
  #Constructs multiple bound tables per beta in beta range using KEGG_ID of significant.feature table
  #Assumes columns wilcoxon.stat.norm (u(j)) and v.METABOLITE (v(j)) are present
  #Assumes all values are between 0-1!!!
  
  require(data.table)
  require(base)
  
  #Filter for columns of interest
  V.Table<-feature.sig.table[,c("KEGG_ID", "wilcoxon.stat.norm", "v.METABOLITE"),with=F]

  #APPLY MAX NORMALIZATION before calculating [KEGG_ID, v.METABOLITE, u.METABOLITE]
  OMEGA<-max(V.Table$v.METABOLITE)/max(V.Table$wilcoxon.stat.norm)
  V.Table$u.METABOLITE<-V.Table$wilcoxon.stat.norm * OMEGA
  V.Table$wilcoxon.stat.norm<-NULL
  
  #Calculate V
  dummy.result<-data.frame(KEGG_ID=c(), v.METABOLITE=c(), u.METABOLITE=c(), V.METABOLITE=c(), BETA=c())
  for (beta in beta.range) {
   
    temp.dummy<-copy(V.Table)
    temp.dummy$V.METABOLITE<-temp.dummy$u.METABOLITE*beta + temp.dummy$v.METABOLITE*(1-beta)
    temp.dummy$BETA<-beta
    
    dummy.result<-rbind(dummy.result, temp.dummy)
  }
  
  #Clean up and return
  dummy.result<-as.data.table(dummy.result)
  return(dummy.result)
  
}

BRCA.Table.V<-Function.V.construct(test.sig, seq(0,1,0.1))
saveRDS(BRCA.Table.V, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.V.rds")

ggplot(BRCA.Table.V, aes(x=factor(BETA), y=V.METABOLITE)) +geom_boxplot()

#Calculate W - WRITE APPRORPIATE SCRIPT !!!!!
Function.W.construct<-function(expression.matrices.rds, processed.table.3.rds, Table.V.rds, differential.threshold=0.1){
  #This function is a modified version of calculation the network edges W(j,j')
  #   It is bi-modal, only assign a non-zero edge if enzymes producing metabolites are differentially expressed
  #   It accounts for the fact that normal levels of metabolic enzymes seen in cancer may be stable across two metabolites
  #   This will cause for there to be a meaningul correlation between two metabolites where there is none 
  #   Check previous way of calculating Table.W for comparisson
  #Normal and cancer.expression.samples are from processed agilent data
  #   These are in the form expression.samples (an object of class "matrix)
  #differential.threshold tells the function at least how many differential genes do we want to be differentially expressed to be considered an edge
  #   The values could be between 0-1.0
  #   Default is at least 10% (0.1)
  
  require(data.table)
  require(base)
  require(reshape)
  
  #Load tables
  CANCER.MATRICES<-readRDS(expression.matrices.rds)
  dummy.table.3<-readRDS(processed.table.3.rds)
  Table.V<-readRDS(Table.V.rds)
  
  #STEP 0 - RNAseq differential expression analysis to get differentially expressed genes
  DIFF.EXP<-Function.RNAseq.Differential.Expression(CANCER.MATRICES$normal, CANCER.MATRICES$tumor)
  print (head(DIFF.EXP))
  if ("ID" %in% colnames(DIFF.EXP)){
    DIFF.EXP.GENES<-as.vector(DIFF.EXP[DIFF.EXP$adj.P.Val<0.05,]$ID)
  } else {
    DIFF.EXP$ID<-rownames(DIFF.EXP)
    DIFF.EXP.GENES<-as.vector(DIFF.EXP[DIFF.EXP$adj.P.Val<0.05,]$ID)
  }
  print ((length(DIFF.EXP.GENES)))
  
  #STEP 1 - Use network.kegg.ids to filter table.3 (enzyme table) - THIS WILL BE THE TOTAL NUMBER OF NODES IN NETWORK
  network.KEGG.IDs<-unique(as.vector(Table.V$KEGG_ID))
  dummy.table.3$Enzyme<-as.character(dummy.table.3$Enzyme)
  dummy.table.3<-dummy.table.3[KEGG_ID %in% network.KEGG.IDs,]
  print ("STEP 1 - DONE")
  
  #STEP 2 - Construct background network (ZERO NETWORK)
  dummy.zero.network<-as.data.table(t(combn(unique(as.vector(dummy.table.3$KEGG_ID)),2)))
  setnames(dummy.zero.network, c("V1","V2"),c("KEGG.ID.1","KEGG.ID.2"))
  dummy.zero.network$W.METABOLITE<-0
  print ("STEP 2 - DONE")
  print (dim(dummy.zero.network))
  
  #STEP 3 Internal Function
  
  #STEP 3 - Use genes that are differentially expressed to drop metabolites that have less differentially expressed genes than threshold [Enzyme, KEGG_ID]
  dummy.sign.metabolites.table<-dummy.table.3[,list(DIFF.EXP=(length(intersect(Enzyme, DIFF.EXP.GENES))>=
                                                                length(Enzyme)*differential.threshold)), 
                                              by="KEGG_ID"]
  
  dummy.sign.metabolites<-unique(as.vector(dummy.sign.metabolites.table[DIFF.EXP==TRUE,]$KEGG_ID)) #Metabolites that pass threshold
  dummy.table.3<-dummy.table.3[KEGG_ID %in% dummy.sign.metabolites,] #Only metabolites that have above threshold diff exp genes and all of their enzymes,
  #including those that are not differentially expressed (since rest passed threshold)
  print ("STEP 3 - DONE")
  
  #STEP 4 - Reduce cancer.expression matrix to contain only enzyme values for metabolites that pass threshold (dummy.table.3)
  dummy.cancer.sign.expression<-as.data.frame(CANCER.MATRICES$tumor,keep.rownames=T)
  dummy.cancer.sign.expression$Enzyme<-rownames(dummy.cancer.sign.expression) #To merge with sign metabolites later
  dummy.cancer.sign.expression<-as.data.table(dummy.cancer.sign.expression)
  
  dummy.cancer.sign.expression<-merge(dummy.cancer.sign.expression, dummy.table.3, by="Enzyme") #To filter out non-enzymes and enzymes from non-diff exp metabolites
  dummy.cancer.sign.expression[,Enzyme:=NULL] #Delete enzyme column to do mean calculation
  print ("STEP 4 - DONE")
  
  #STEP 5 - Calculate the means of compound enzymes per patient (to create vector for all sign metabolites across vectors)
  dummy.cancer.sign.expression.mean<-cast(melt(dummy.cancer.sign.expression, id=c("KEGG_ID")), KEGG_ID~variable, mean)
  dummy.cancer.sign.expression.mean<-as.data.frame(dummy.cancer.sign.expression.mean[complete.cases(dummy.cancer.sign.expression.mean),])
  print ("STEP 5 - DONE")
  
  #STEP 6 - Perform ABSOLUTE spearman correlation across all metabolite vectors - This is W(j,j')
  rownames(dummy.cancer.sign.expression.mean)<-as.character(dummy.cancer.sign.expression.mean$KEGG_ID) #To build correlation matrix
  dummy.cancer.sign.expression.mean$KEGG_ID<-NULL
  dummy.pre.W.table<-cor(t(dummy.cancer.sign.expression.mean), method="spearman")
  dummy.pre.W.table<-melt(dummy.pre.W.table) #To convert it to a pairwise table
  dummy.pre.W.table<-dummy.pre.W.table[dummy.pre.W.table$X1!=dummy.pre.W.table$X2,] #To delete diagonal correlation (node to itself)
  colnames(dummy.pre.W.table)<-c("KEGG.ID.1","KEGG.ID.2","W.METABOLITE")
  dummy.pre.W.table$W.METABOLITE<-abs(dummy.pre.W.table$W.METABOLITE)
  print ("STEP 6 - DONE")
  
  #STEP 7 - Integrate with background table to obtain final W.table (to have info of 0 edge weight)
  dummy.W.table<-as.data.frame(rbind(dummy.pre.W.table, dummy.zero.network)) #Mix background and sign metabolic tables
  print (dim(dummy.W.table))
  dummy.W.table<-dummy.W.table[!duplicated(dummy.W.table[,1:2]),]#Removes zeroes from background for diff metabolites
  #This works because pairwise table above has two-sided combinations
  print (dim(dummy.W.table))
  dummy.W.table.combn<-as.data.frame(dummy.zero.network)[,1:2]#To remove two sided duplicates - this has all the edge permutations we need
  
  dummy.W.table<-as.data.table(merge(dummy.W.table, dummy.W.table.combn))
  print (dim(dummy.W.table))
  print ("STEP 7 - DONE")
  
  #STEP 8 - Clean table
  dummy.W.table[is.na(dummy.W.table)]<-0
  
  #END
  return (dummy.W.table)
}

BRCA.Table.W<-Function.W.construct("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061314.BRCA.CANCER.MATRICES",
                                   "PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.3.rds",
                                   "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.V.rds",
                                   0.5)
saveRDS(BRCA.Table.W, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.W.rds")

ggplot(BRCA.Table.W, aes(W.METABOLITE)) + geom_histogram()

#Integrate V and W to do clustering analysis - WRITE APPROPRIATE SCRIPT!!!!!!
Function.W.V.construct<-function(Table.V.rds, Table.W.rds) {
  
  require(base)
  require(data.table)
  
  #Load tables
  V.table.beta<-readRDS(Table.V.rds)
  W.table<-readRDS(Table.W.rds)
  
  #Filter for columns we need
  V.table.beta<-V.table.beta[,c("KEGG_ID", "V.METABOLITE", "BETA"), with=F]
  
  #Create lists of V.table per beta
  V.table.beta.list<-split(V.table.beta, V.table.beta$BETA)
  
  #Obtain integrated V.W edges
  W.V.TABLE<-data.frame(KEGG.ID.1=c(), KEGG.ID.2=c(), INTEGRATED.W=c(), BETA=c())
  for (beta in names(V.table.beta.list)) {
    
    #[KEGG_ID, V.METABOLITE]
    dummy<-V.table.beta.list[[beta]]
    dummy$BETA<-NULL
    
    #Do first node assignment 
    setnames(dummy, colnames(dummy), c("KEGG.ID.1", "NODE.WEIGHT.1"))
    v.w.table<-merge(x=dummy, y=W.table, by="KEGG.ID.1")
    
    #Do second node assignment
    setnames(dummy, colnames(dummy), c("KEGG.ID.2", "NODE.WEIGHT.2"))
    v.w.table<-merge(x=v.w.table, y=dummy, by="KEGG.ID.2")
    
    #Get integrated edge weight with normalized omega parameters
    OMEGA.NODE.1<-max(v.w.table$W.METABOLITE)/max(v.w.table$NODE.WEIGHT.1)
    OMEGA.NODE.2<-max(v.w.table$W.METABOLITE)/max(v.w.table$NODE.WEIGHT.2)
    v.w.table$INTEGRATED.W<-ifelse(v.w.table$W.METABOLITE!=0, v.w.table$NODE.WEIGHT.1*OMEGA.NODE.1 + v.w.table$NODE.WEIGHT.2*OMEGA.NODE.2 + v.w.table$W.METABOLITE,
                                   0)
    
    #Normalize to 0-1 for later cluster analysis
    v.w.table$INTEGRATED.W<-normalize.vector(v.w.table$INTEGRATED.W)
    
    #Add to table
    v.w.table<-as.data.table(v.w.table)
    v.w.table$BETA<-beta
    v.w.table<-v.w.table[,c("KEGG.ID.1","KEGG.ID.2", "INTEGRATED.W", "BETA"), with=F]
    W.V.TABLE<-rbind(W.V.TABLE, v.w.table)
  }  
  
  #Return
  return(W.V.TABLE)
}

BRCA.Table.W.V<-Function.W.V.construct("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.V.rds",
                                     "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.W.rds")
saveRDS(BRCA.Table.W.V, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071514.BRCA.Table.W.V.rds")
write.table(file="PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/071514.BRCA.Table.W.V", BRCA.Table.W.V, quote=F, row.names=F)

ggplot(BRCA.Table.W.V, aes(x=factor(BETA), y=INTEGRATED.W)) +geom_boxplot()
ggplot(BRCA.Table.W.V, aes(INTEGRATED.W)) + geom_histogram() + facet_wrap(~BETA)

#Find KEGG pathway enirchment in clusters - WRITE FORMAL FUNCTIONS!!
test.ccsf<-common.cluster.spici.function("PIPELINES/METABOLIC.DRIVERS/CLUSTER.ANALYSIS/BRCA/SPICI.RESULTS", "")
test.ccsf$FILENAME<-as.character(test.ccsf$FILENAME)
test.ccsf$BETA<-sapply(test.ccsf$FILENAME, function(x)  strsplit(strsplit(x, "_")[[1]][2],".cluster" )[[1]][1] )

test.capvfv3<-cluster.analysis.p.value.function.V3(test.ccsf, BRCA.Table.W.V)
test.cckpe<-cluster.cpd.kegg.pathways.enrichment(test.capvfv3, "DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND", 
                                                 unique(c(as.vector(BRCA.Table.W.V$KEGG.ID.1) ,as.vector(BRCA.Table.W.V$KEGG.ID.2))),0.8,0.8)
head(test.cckpe)
test.capvfv3

Table.enrichment.plot.1.function(test.cckpe, unique(as.vector(test.cckpe$Beta)))
Table.enrichment.plot.2.function(test.cckpe, unique(as.vector(test.cckpe$Beta)))

Function.kegg.path.enrich<-function(capvfv3, pathway.file, kegg.in.network, cluster.p.cutoff, kegg.p.cutoff) {
  
  require(data.table)
  
  #Load pathway file
  path.table<-as.data.table(read.csv(pathway.file, header=T, sep="\t", colClasses="character"))
  path.table$DESCRIPTION<-gsub(" \\(human)", "", path.table$DESCRIPTION)
  
  #Filter cluster.table for p.val threshold [Beta, CLUSTER.SIZE, P.VALUE, CLUSTER, ADJ.P.VALUE]
  capvfv3<-as.data.table(capvfv3)
  cluster.table<-capvfv3[ADJ.P.VALUE<cluster.p.cutoff,]
  
  #Filter pathway table for metabolites that we are only analyzing in cluster (to save time) [PATHWAY, COMPOUND, DESCRIPTION]
  cluster.metabolites<-unique(unlist(strsplit(as.vector(cluster.table$CLUSTER),"_")))
  path.table<-path.table[COMPOUND %in% cluster.metabolites, ]
  
  #Get a count for how many network metabolites are in each pathway in pathway.table [PATHWAY, COMPOUND, DESCRIPTION, IN.NETWORK]
  path.count<-copy(path.table)
  path.count<-path.count[,list(IN.NETWORK=length(intersect(COMPOUND,kegg.in.network))),by=c("PATHWAY", "DESCRIPTION")]
  path.table<-as.data.table(merge(as.data.frame(path.table), as.data.frame(path.count), by=c("PATHWAY", "DESCRIPTION")))
  
  #Loop through clusters to get pathway p.val
  dummy.result<-data.frame(CLUSTER=c(), PATHWAY=c(), PATHWAY.DESCRIPTION=c(), PATHWAY.PVALUE=c(), PATHWAY.COVERAGE=c())
  for (cluster in as.vector(cluster.table$CLUSTER)) {
    
    #Add cluster to be analyzed (met1,met2,...)
    cluster.metabolites<-strsplit(cluster,"_")[[1]]
    
    #Calculate p.value per pathway for cluster
    path.table.pval<-path.table[,list(PATHWAY.PVALUE=phyper(q=length(intersect(cluster.metabolites, COMPOUND))-1,
                                                            m=unique(as.vector(IN.NETWORK)),
                                                            n=length(kegg.in.network) - unique(as.vector(IN.NETWORK)),
                                                            k=length(cluster.metabolites), lower.tail=F),
                                      PATHWAY.COVERAGE=length(intersect(cluster.metabolites, COMPOUND))/length(cluster.metabolites)
                                      ), 
                                by=c("PATHWAY", "DESCRIPTION")]
    
    #Do multiple hypothesis correction
    path.table.pval$PATHWAY.ADJ.PVALUE<-p.adjust(path.table.pval$PATHWAY.PVALUE, method="fdr")
    
    #Add cluster information
    path.table.pval$CLUSTER<-cluster
    
    #Add to result
    dummy.result<-rbind(dummy.result, path.table.pval)
    
  }
  
  #Merge result with cluster.table
  dummy.result<-as.data.table(merge(as.data.frame(dummy.result), as.data.frame(cluster.table), by="CLUSTER"))
  
  #Filter for p.val threshold
  dummy.result<-dummy.result[PATHWAY.PVALUE<kegg.p.cutoff,]
  
  #Return
  return(dummy.result)
}

BRCA.KEGG.ENRICH<-Function.kegg.path.enrich(test.capvfv3,"DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND", 
                                            unique(c(as.vector(BRCA.Table.W.V$KEGG.ID.1) ,as.vector(BRCA.Table.W.V$KEGG.ID.2))),0.05,0.05)
BRCA.KEGG.ENRICH
Table.enrichment.plot.1.function(as.data.frame(BRCA.KEGG.ENRICH), unique(as.vector(BRCA.KEGG.ENRICH$Beta)))
Table.enrichment.plot.2.function(as.data.frame(BRCA.KEGG.ENRICH), unique(as.vector(BRCA.KEGG.ENRICH$Beta)))
head(BRCA.KEGG.ENRICH)

#Validation
test.sig
BUD.METABOLITES$KEGG_ID<-as.character(BUD.METABOLITES$KEGG_ID)
BUD.METABOLITES[KEGG_ID %in% as.character(test.sig$KEGG_ID),] #FOUND NONE

COSMIC<-as.data.table(read.csv("DATABASES/CANCER_DATA/COSMIC/cancer_gene_census.csv", header=T,sep=","))
COSMIC.BRCA<-COSMIC[grepl("breast",Tumour.Types...Somatic.Mutations.,ignore.case=T) ,] #17 BRCA cosmic genes

test.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")
length(unique(as.vector(test.2[KEGG_ID %in% test.sig$KEGG_ID,]$GENE))) #1295
COSMIC.BRCA[Symbol %in% unique(as.vector(test.2[KEGG_ID %in% test.sig$KEGG_ID,]$GENE)),]
sum(as.vector(COSMIC.BRCA$Symbol) %in% unique(as.vector(test.2[KEGG_ID %in% test.sig$KEGG_ID,]$GENE))) #5 found in sig.test out of possible 6
sum(as.vector(COSMIC.BRCA$Symbol) %in% unique(as.vector(test.3[KEGG_ID %in% test.sig$KEGG_ID,]$Enzyme)))
sum(as.vector(COSMIC$Symbol) %in% unique(as.vector(test.2[KEGG_ID %in% test.sig$KEGG_ID,]$GENE)))

sum(as.vector(COSMIC.BRCA$Symbol) %in% unique(as.vector(test.2$GENE)))
sum(as.vector(COSMIC.BRCA$Symbol) %in% unique(as.vector(test.3$Enzyme)))
length(unique(as.vector(test.2$GENE))) #3897
phyper(5-1,6,3897-6,1295, lower.tail=F )

#####
setnames(test.3, colnames(test.3), c("Hugo_Symbol", "KEGG_ID"))
sample.test<-as.data.table(merge(as.data.frame(test.3), as.data.frame(test.1$table.1), by="Hugo_Symbol"))
sample.patient<-sample.test[,list(N.PATIENT=length(unique(Tumor_Sample_Barcode))), by="KEGG_ID"]
sample.patient[N.PATIENT<2,]

LOG.1<-read.csv("LOGS/LOG1.csv")
LOG.1$PROCESSED<-substr(LOG.1$X.1...C00032., 6,11)
test[KEGG_ID %in% LOG.1$PROCESSED,]
test[v.METABOLITE==0,]

##########071814#########

#Rank Validation
Table.u.BRCA<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/071814.BRCA.Table.u.rds")
Table.v.BRCA.PVAL
table.2<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/061214_Table.2.rds")

Function.j.rank.enrichment<-function(table.j, interest.columns, dummy.table.2, target.genes, normalize=F) {
  #Rank cummulative and non-cummulative enrichment for metabolites j in order of decreasing SCORE
  #interest.columns assumes that first column is KEGG_ID and second column is scores  ie. [KEGG_ID, wilcoxon.stat]
  
  require(data.table)
  require(base)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-process tables 
  table.j<-table.j[,interest.columns, with=F] 
  if (normalize==T) {
    table.j$SCORE<-normalize.vector(table.j[,2, with=F])
    table.j<-table.j[,c(1,3), with=F]
  }
  setnames(table.j, colnames(table.j), c("KEGG_ID", "SCORE"))
  
  dummy.table.2<-unique(dummy.table.2[,2:3, with=F])
  setnames(dummy.table.2, colnames(dummy.table.2), c("KEGG_ID", "Hugo_Symbol")) #[KEGG_ID, Hugo_Symbol]
  
  #Order j.table by score to calculate cummulative enrichemnt [KEGG_ID, SCORE]
  table.j<-table.j[order(SCORE, decreasing=T),]
  
  #CALCULATIONS
  background.genes<-unique(as.vector(dummy.table.2[KEGG_ID %in% as.vector(table.j$KEGG_ID),]$Hugo_Symbol)) #Genes that we can choose from (All in urn)
  target.in.j<-intersect(background.genes, target.genes) #Targets we could have chosen (White in urn)
  
  #FIRST CUMMULATIVE
  #Calculate p.value
  CUM.SCORE.PVALS<-sapply(1:nrow(table.j), function(x)  
    phyper(q=length(intersect(unique(as.vector(dummy.table.2[KEGG_ID %in% as.vector(table.j[1:x]$KEGG_ID),]$Hugo_Symbol)),target.in.j))-1,
           m=length(target.in.j),
           n=length(background.genes)-length(target.in.j),
           k= length(unique(as.vector(dummy.table.2[KEGG_ID %in% as.vector(table.j[1:x]$KEGG_ID),]$Hugo_Symbol))),
           lower.tail=F
              ))  
  
  #Calculate coverage
  CUM.SCORE.COVERAGE<-sapply(1:nrow(table.j), function (x) 
    length(intersect(unique(as.vector(dummy.table.2[KEGG_ID %in% as.vector(table.j[1:x]$KEGG_ID),]$Hugo_Symbol)),target.in.j))/length(target.in.j) )
  
  #THEN NON-CUMMULATIVE
  table.non.cum<-as.data.table(merge(as.data.frame(table.j), as.data.frame(dummy.table.2), by="KEGG_ID")) #[KEGG_ID, SCORE, Hugo_Symbol]
  table.non.cum<-table.non.cum[,list(NON.CUM.P.VAL=phyper(q=length(intersect(Hugo_Symbol, target.in.j))-1,
                                                                m=length(target.in.j),
                                                                n=length(background.genes)-length(target.in.j),
                                                                k=length(unique(Hugo_Symbol)), lower.tail=F),
                                     NON.CUM.SCORE.COVERAGE=length(intersect(Hugo_Symbol, target.in.j))/length(target.in.j)),
                               by=c("KEGG_ID", "SCORE")]
  
  #Assign to KEGG_IDs
  table.j$CUM.SCORE.COVERAGE<-CUM.SCORE.COVERAGE
  table.j$CUM.P.VAL<-CUM.SCORE.PVALS
  table.j$CUM.ADJ.P.VAL<-p.adjust(table.j$CUM.P.VAL, method="fdr")
  table.j<-as.data.table(merge(as.data.frame(table.j), as.data.frame(table.non.cum), by=c("KEGG_ID", "SCORE")))
  table.j$NON.CUM.ADJ.P.VAL<-p.adjust(table.j$NON.CUM.P.VAL, method="fdr")
  
  #Return
  return(table.j)
}
Function.j.rank.enrich.plot<-function(Table.enrich.result){
  #Plots results of Function.j.rank.enrichment()
  
  require(ggplot2)
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # Multiple plot function
    #
    # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    # - cols:   Number of columns in layout
    # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
    #
    # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    # then plot 1 will go in the upper left, 2 will go in the upper right, and
    # 3 will go all the way across the bottom.
    #
    
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  p1<-ggplot(Table.enrich.result$CUM.TABLE, aes(x=BREAKS, y=CUM.SCORE.COVERAGE)) +geom_point(size=4) +theme.format
  p3<-ggplot(Table.enrich.result$CUM.TABLE, aes(x=BREAKS, y=-log(CUM.ADJ.P.VAL))) + geom_point(size=4) +theme.format + 
    geom_hline(yintercept=-log(0.05), colour="red")
  p2<-ggplot(Table.enrich.result$NON.CUM.TABLE, aes(x=SCORE, y=NON.CUM.SCORE.COVERAGE)) +geom_point(size=4) +theme.format
  p4<-ggplot(Table.enrich.result$NON.CUM.TABLE, aes(x=SCORE, y=-log(NON.CUM.ADJ.P.VAL))) + geom_point(size=4) +theme.format + 
    geom_hline(yintercept=-log(0.05), colour="red")
  
  multiplot(p1,p2,p3,p4, cols=2)
}

#With COSMIC
COSMIC.BRCA
Table.BRCA.u.rank.enrich<-Function.j.rank.enrichment(Table.u.BRCA$STATS, c("KEGG_ID", "wilcoxon.stat"), table.2, unique(as.vector(COSMIC.BRCA$Symbol)), normalize=T)
Function.j.rank.enrich.plot(Table.BRCA.u.rank.enrich)

Table.BRCA.v.rank.enrich<-Function.j.rank.enrichment(Table.v.BRCA.PVAL, c("KEGG_ID", "v.METABOLITE"), table.2, as.vector(COSMIC.BRCA$Symbol))
Function.j.rank.enrich.plot(Table.BRCA.v.rank.enrich)

#With Curated Mutsig
MUTSIG.BRCA.CURATED<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072014.BRCA.CURATED.RESULTS.sig_genes.txt",header=T,sep="\t"))
MUTSIG.BRCA.CURATED<-MUTSIG.BRCA.CURATED[q<0.05,]
MUTSIG.BRCA<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072014.BRCA.RESULTS.sig_genes.txt", header=T, sep="\t"))
MUTSIG.BRCA<-MUTSIG.BRCA[q<0.05,]

Table.BRCA.u.rank.enrich<-Function.j.rank.enrichment(Table.u.BRCA$STATS, c("KEGG_ID", "wilcoxon.stat"), table.2, unique(as.vector(MUTSIG.BRCA.CURATED$gene)), normalize=T)
Function.j.rank.enrich.plot(Table.BRCA.u.rank.enrich)

Table.BRCA.v.rank.enrich<-Function.j.rank.enrichment(Table.v.BRCA.PVAL, c("KEGG_ID", "v.METABOLITE"), table.2, as.vector(MUTSIG.BRCA.CURATED$gene))
Function.j.rank.enrich.plot(Table.BRCA.v.rank.enrich)

###########072014############
#Enrichment based on u(p) and v(p)
Table.u.BRCA.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/072014.BRCA.Table.u.p.rds")
Table.u.BRCA.p$STATS
Table.v.BRCA.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/072014.BRCA.Table.v.p.rds")

Function.p.rank.enrichment<-function(table.p, interest.columns, target.genes, normalize=F) {
  #Rank cummulative and non-cummulative enrichment for proteins p in order of decreasing SCORE
  #interest.columns assumes that first column is ID and second column is scores  ie. [KEGG_ID, wilcoxon.stat]
  #Calculates PRECISSION-RECALL values as well for each cummulative count
  
  require(data.table)
  require(base)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-process tables 
  table.p<-table.p[,interest.columns, with=F] 
  if (normalize==T) {
    table.p$SCORE<-normalize.vector(table.p[,2, with=F])
    table.p<-table.p[,c(1,3), with=F]
  }
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  
  #Order j.table by score to calculate cummulative enrichemnt [Hugo_Symbol, SCORE]
  table.p<-table.p[order(SCORE, decreasing=T),]
  
  #CALCULATIONS
  background.genes<-unique(as.vector(table.p$Hugo_Symbol)) #Genes that we can choose from
  target.in.p<-intersect(background.genes, target.genes) #Targets we could have chosen (White in urn)
  
  #Breaks for tables 0-1
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.0001)
  PR.BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  
  #FIRST CUMMULATIVE [BREAKS, CUM.SCORE.PVALS, CUM.SCORE.COVERAGE]
  #Calculate p.value
  CUM.P.VAL<-sapply(PR.BREAKS, function(x)  
    phyper(q=length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))-1,
           m=length(target.in.p),
           n=length(background.genes)-length(target.in.p),
           k= length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))),
           lower.tail=F
    ))  
  
  #Calculate coverage
  CUM.SCORE.COVERAGE<-sapply(PR.BREAKS, function (x) 
    length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/length(target.in.p) )
  
  #Calcualte fdr for CUM
  CUM.TABLE<-as.data.table(cbind(PR.BREAKS, CUM.P.VAL, CUM.SCORE.COVERAGE))
  setnames(CUM.TABLE, colnames(CUM.TABLE), c("BREAKS", "CUM.P.VAL", "CUM.SCORE.COVERAGE"))
  CUM.TABLE$CUM.ADJ.P.VAL<-p.adjust(CUM.P.VAL, method="fdr")
  
  #THEN NON-CUMMULATIVE - [Hugo_Symbol, SCORE, NON.CUM.P.VAL, NON.CUM.SCORE.COVERAGE]
  table.p$NON.CUM.P.VAL<-sapply(as.vector(table.p$Hugo_Symbol), function(z) phyper(q=length(intersect(z, target.in.p))-1,
                                                                                  m=length(target.in.p),
                                                                                  n=length(background.genes)-length(target.in.p),
                                                                                  k=length(unique(z)), lower.tail=F))
  table.p$NON.CUM.SCORE.COVERAGE<-sapply(as.vector(table.p$Hugo_Symbol), function(z) length(intersect(z, target.in.p))/length(target.in.p))
  
  #Assign to Hugo_Symbols for NON.CUM
  table.p$NON.CUM.ADJ.P.VAL<-p.adjust(table.p$NON.CUM.P.VAL, method="fdr")
  
  #Lastly add CUM.RECALL and CUM.PRECISION columns to CUM
  CUM.TABLE$CUM.RECALL<-CUM.SCORE.COVERAGE
  CUM.TABLE$CUM.PRECISION<-sapply(PR.BREAKS, function (x) 
    length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/ length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))))
  
  #Return as list
  return(list(CUM.TABLE=CUM.TABLE,NON.CUM.TABLE=table.p))
}

Table.BRCA.u.rank.enrich.p<-Function.p.rank.enrichment(Table.u.BRCA.p$STATS, c("Hugo_Symbol", "paired.t.stat"), as.vector(COSMIC.BRCA$Symbol), normalize=T)
Function.j.rank.enrich.plot(Table.BRCA.u.rank.enrich.p)

Table.BRCA.v.rank.enrich.p<-Function.p.rank.enrichment(Table.v.BRCA.p, c("Hugo_Symbol", "v.PROTEIN"), as.vector(COSMIC.BRCA$Symbol), normalize=F)
Function.j.rank.enrich.plot(Table.BRCA.v.rank.enrich.p)

#Enrichment based on u(path) and v(path)
Table.u.BRCA.path<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/072114.BRCA.Table.u.path.rds")
Table.u.BRCA.path$STATS[,]

Function.path.rank.enrichment<-function(table.j, interest.columns, pathway.uniprot.file, uniprot.mapping.file, target.genes, normalize=F) {
  #Rank cummulative and non-cummulative enrichment for metabolites j in order of decreasing SCORE
  #interest.columns assumes that first column is KEGG_ID and second column is scores  ie. [KEGG_ID, wilcoxon.stat]
  
  require(data.table)
  require(base)
  
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-process tables 
  table.j<-table.j[,interest.columns, with=F] 
  if (normalize==T) {
    table.j$SCORE<-normalize.vector(table.j[,2, with=F])
    table.j<-table.j[,c(1,3), with=F]
  }
  setnames(table.j, colnames(table.j), c("PATH", "SCORE"))
  
  #Process pathway files
  uniprot.mapping<-fread(uniprot.mapping.file, header=T, sep="\t", drop=c(2:4, 6:7)) #[UNIPROT, Hugo_Symbol]  
  setnames(uniprot.mapping, colnames(uniprot.mapping), c("UNIPROT", "GENES"))
  uniprot.mapping$Hugo_Symbol<-sapply(uniprot.mapping$GENES, function(x) strsplit(x," ")[[1]][1])
  uniprot.mapping$GENES<-NULL
  
  pathway.uniprot<-as.data.table(read.csv(pathway.uniprot.file, header=F, sep="\t")) #[UNIPROT, PATH]
  pathway.uniprot<-pathway.uniprot[,c(1,4), with=F]
  setnames(pathway.uniprot, colnames(pathway.uniprot), c("UNIPROT", "PATH"))
  
  dummy.table.2<-as.data.table(merge(as.data.frame(uniprot.mapping), as.data.frame(pathway.uniprot), by="UNIPROT")) #[PATH, Hugo_Symbol]
  dummy.table.2<-unique(dummy.table.2[,c("PATH", "Hugo_Symbol"), with=F])
  
  #Order j.table by score to calculate cummulative enrichemnt [PATH, SCORE]
  table.j<-table.j[order(SCORE, decreasing=T),]
  
  #CALCULATIONS
  background.genes<-unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j$PATH),]$Hugo_Symbol)) #Genes that we can choose from (All in urn)
  target.in.j<-intersect(background.genes, target.genes) #Targets we could have chosen (White in urn)
  
  #FIRST CUMMULATIVE
  #Calculate p.value
  CUM.SCORE.PVALS<-sapply(1:nrow(table.j), function(x)  
    phyper(q=length(intersect(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol)),target.in.j))-1,
           m=length(target.in.j),
           n=length(background.genes)-length(target.in.j),
           k= length(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol))),
           lower.tail=F
    ))  
  
  #Calculate coverage
  CUM.SCORE.COVERAGE<-sapply(1:nrow(table.j), function (x) 
    length(intersect(unique(as.vector(dummy.table.2[PATH %in% as.vector(table.j[1:x]$PATH),]$Hugo_Symbol)),target.in.j))/length(target.in.j) )
  
  #THEN NON-CUMMULATIVE
  table.non.cum<-as.data.table(merge(as.data.frame(table.j), as.data.frame(dummy.table.2), by="PATH")) #[PATH, SCORE, Hugo_Symbol]
  table.non.cum<-table.non.cum[,list(NON.CUM.P.VAL=phyper(q=length(intersect(Hugo_Symbol, target.in.j))-1,
                                                          m=length(target.in.j),
                                                          n=length(background.genes)-length(target.in.j),
                                                          k=length(unique(Hugo_Symbol)), lower.tail=F),
                                     NON.CUM.SCORE.COVERAGE=length(intersect(Hugo_Symbol, target.in.j))/length(target.in.j)),
                               by=c("PATH", "SCORE")]
  
  #Assign to PATHs
  table.j$CUM.SCORE.COVERAGE<-CUM.SCORE.COVERAGE
  table.j$CUM.P.VAL<-CUM.SCORE.PVALS
  table.j$CUM.ADJ.P.VAL<-p.adjust(table.j$CUM.P.VAL, method="fdr")
  table.j<-as.data.table(merge(as.data.frame(table.j), as.data.frame(table.non.cum), by=c("PATH", "SCORE")))
  table.j$NON.CUM.ADJ.P.VAL<-p.adjust(table.j$NON.CUM.P.VAL, method="fdr")
  
  #Return
  return(table.j)
}

Table.BRCA.u.rank.enrich.path<-Function.path.rank.enrichment(Table.u.BRCA.path$STATS, c("PATH", "wilcoxon.stat"), 
                                                             "DATABASES/NCI.PATHWAYS/uniprot.tab","DATABASES/UNIPROT/070714_UNIPROT_MAPPING.tab",
                                                             as.vector(COSMIC.BRCA$Symbol), normalize=T)
Function.j.rank.enrich.plot(Table.BRCA.u.rank.enrich.path)

#Enrichement based on u(MUTSIG)
MUTSIG.BRCA<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072014.BRCA.CURATED.RESULTS.sig_genes.txt", header=T, sep="\t"))
MUTSIG.BRCA$NEG.LOG.q<--log(MUTSIG.BRCA$q)
MUTSIG.BRCA$NEG.LOG.q[MUTSIG.BRCA$NEG.LOG.q==Inf]<-60
MUTSIG.BRCA$RANK<-rank(-MUTSIG.BRCA$q)
head(MUTSIG.BRCA)
Table.MUTSIG.BRCA.u.rank.enrich.p<-Function.p.rank.enrichment(MUTSIG.BRCA, c("gene", "RANK"), as.vector(COSMIC.BRCA$Symbol), normalize=T)
Function.j.rank.enrich.plot(Table.MUTSIG.BRCA.u.rank.enrich.p)

sum(COSMIC.COAD$Symbol %in% MUTSIG.COAD$gene)
sum(COSMIC.COAD$Symbol %in% Table.u.COAD.p$STATS$Hugo_Symbol)
sum(COSMIC.BRCA$Symbol %in% Table.v.BRCA.p$Hugo_Symbol)

log.1<-Table.BRCA.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.1$METHOD<-"BMR"
log.2<-Table.MUTSIG.BRCA.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.2$METHOD<-"MUTSIG"
log.3<-Table.BRCA.v.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.3$METHOD<-"DIFF.EXP.RATIO"
log<-rbind(log.1, log.2,log.3)

p1<-ggplot(log, aes(x=BREAKS, y=CUM.SCORE.COVERAGE, colour=METHOD)) + geom_point(size=4) + theme.format + geom_line()
p2<-ggplot(log, aes(x=BREAKS, y=-log(CUM.ADJ.P.VAL), colour=METHOD)) + geom_point(size=4) + theme.format + geom_line()+
  geom_hline(yintercept=-log(0.05), colour="red")
multiplot(p1,p2, cols=2)
ggplot(log, aes(x=CUM.RECALL, y=CUM.PRECISION, colour=METHOD)) + geom_line() + theme.format + geom_point(size=4)

#Enrichment for other cancers and PR curves
COSMIC.GBM<-COSMIC[grepl("glioblastoma",Tumour.Types...Somatic.Mutations.,ignore.case=T) ,] #11 GBM cosmic genes 
COSMIC.COAD<-COSMIC[grepl("colorectal",Tumour.Types...Somatic.Mutations.,ignore.case=T) ,] #21 COAD cosmic genes
COSMIC.UCEC<-COSMIC[grepl("endometrial",Tumour.Types...Somatic.Mutations.,ignore.case=T) ,] #10 COAD cosmic genes
COSMIC[grepl("endometrial",Tumour.Types...Somatic.Mutations.,ignore.case=T) ,]

MUTSIG.GBM<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072314.GBM.RESULTS.sig_genes.txt", header=T, sep="\t"))
MUTSIG.GBM$NEG.LOG.q<--log(MUTSIG.GBM$q)
MUTSIG.GBM$NEG.LOG.q[MUTSIG.GBM$NEG.LOG.q==Inf]<-50
Table.MUTSIG.GBM.u.rank.enrich.p<-Function.p.rank.enrichment(MUTSIG.GBM, c("gene", "NEG.LOG.q"), as.vector(COSMIC.GBM$Symbol), normalize=T)

MUTSIG.COAD<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072314.COAD.RESULTS.sig_genes.txt", header=T, sep="\t"))
MUTSIG.COAD$NEG.LOG.q<--log(MUTSIG.COAD$q)
MUTSIG.COAD$NEG.LOG.q[MUTSIG.COAD$NEG.LOG.q==Inf]<-52
Table.MUTSIG.COAD.u.rank.enrich.p<-Function.p.rank.enrichment(MUTSIG.COAD, c("gene", "NEG.LOG.q"), as.vector(COSMIC.COAD$Symbol), normalize=T)

MUTSIG.UCEC<-as.data.table(read.csv("SOFTWARE/MUTSIG/MutSigCV_1.4/072314.UCEC.RESULTS.sig_genes.txt", header=T, sep="\t"))
MUTSIG.UCEC$NEG.LOG.q<--log(MUTSIG.UCEC$q)
MUTSIG.UCEC$NEG.LOG.q[MUTSIG.UCEC$NEG.LOG.q==Inf]<-54
Table.MUTSIG.UCEC.u.rank.enrich.p<-Function.p.rank.enrichment(MUTSIG.UCEC, c("gene", "NEG.LOG.q"), as.vector(COSMIC.UCEC$Symbol), normalize=T)

Table.u.GBM.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/GBM/072214.GBM.Table.u.p.rds")
Table.GBM.u.rank.enrich.p<-Function.p.rank.enrichment(Table.u.GBM.p$STATS, c("Hugo_Symbol", "paired.t.stat"), as.vector(COSMIC.GBM$Symbol), normalize=T)

Table.u.COAD.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/COAD/072214.COAD.Table.u.p.rds")
Table.COAD.u.rank.enrich.p<-Function.p.rank.enrichment(Table.u.COAD.p$STATS, c("Hugo_Symbol", "paired.t.stat"), as.vector(COSMIC.COAD$Symbol), normalize=T)

Table.u.UCEC.p<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/UCEC/072314.UCEC.Table.u.p.rds")
Table.UCEC.u.rank.enrich.p<-Function.p.rank.enrichment(Table.u.UCEC.p$STATS, c("Hugo_Symbol", "paired.t.stat"), as.vector(COSMIC.UCEC$Symbol), normalize=T)

log.1<-Table.COAD.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.1$CANCER<-"COAD"
log.1$METHOD<-"BMR"
log.2<-Table.MUTSIG.COAD.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.2$CANCER<-"COAD"
log.2$METHOD<-"MUTSIG"
log.3<-Table.GBM.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.3$CANCER<-"GBM"
log.3$METHOD<-"BMR"
log.4<-Table.MUTSIG.GBM.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.4$CANCER<-"GBM"
log.4$METHOD<-"MUTSIG"
log.5<-Table.UCEC.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.5$CANCER<-"UCEC"
log.5$METHOD<-"BMR"
log.6<-Table.MUTSIG.UCEC.u.rank.enrich.p$CUM.TABLE[,c(1,3:6),with=F]
log.6$CANCER<-"UCEC"
log.6$METHOD<-"MUTSIG"

log<-rbind(log.1, log.2, log.3,log.4, log.5, log.6)
ggplot(log, aes(x=CUM.RECALL, y=CUM.PRECISION, colour=METHOD)) + geom_line() + theme.format + geom_point(size=4) +
         facet_wrap(~CANCER)

#PR on BRCA u(p) for breast cancer metabolic census genes vs MUTSIG

COSMIC.BRCA.MET<-as.vector(COSMIC.BRCA[Symbol %in% table.2[KEGG_ID %in% as.vector(Table.u.BRCA$STATS$KEGG_ID),]$GENE, ]$Symbol)
Table.BRCA.u.rank.enrich.on.met.p<-Function.p.rank.enrichment(Table.u.BRCA.p$STATS, c("Hugo_Symbol", "paired.t.stat"), COSMIC.BRCA.MET, normalize=T)
Table.MUTSIG.BRCA.u.rank.enrich.on.met.p<-Function.p.rank.enrichment(MUTSIG.BRCA, c("gene", "NEG.LOG.q"), COSMIC.BRCA.MET, normalize=T)

log.7<-Table.BRCA.u.rank.enrich.on.met.p$CUM.TABLE[,c(1,3:6),with=F]
log.7$CANCER<-"BRCA"
log.7$METHOD<-"BMR"
log.8<-Table.MUTSIG.BRCA.u.rank.enrich.on.met.p$CUM.TABLE[,c(1,3:6),with=F]
log.8$CANCER<-"BRCA"
log.8$METHOD<-"MUTSIG"
log2<-rbind(log.7,log.8)
ggplot(log2, aes(x=CUM.RECALL, y=CUM.PRECISION, colour=METHOD)) + geom_line() + theme.format + geom_point(size=4)

##############072314##############
#POST GM
Table.u.BRCA.p$STATS
Table.v.BRCA.p
Table.u.v.BRCA.p<-as.data.table(merge(as.data.frame(Table.u.BRCA.p$STATS[,c("Hugo_Symbol", "paired.t.stat"),with=F]), as.data.frame(Table.v.BRCA.p), by="Hugo_Symbol"))
ggplot(Table.u.v.BRCA.p, aes(x=v.PROTEIN, y=normalize.vector(paired.t.stat))) + geom_point() + theme.format

dummy.1$table.1

Function.p.rank.patient.enrichment<-function(table.p, interest.columns, processed.table.1.rds, normalize=F) {
  #Ranks cummulative and non-cumulative coverage per protein p by decreasing SCORE
  #Interest columns assumes that first column is ID and second column is scores ie. [KEGG_ID, wilcoxon.stat]
  #Calculates PRECISSION-RECALL values as well for each cummulative count
  
  require(data.table)
  require(base)
  normalize.vector<-function(x){
    y=(x-min(x))/(max(x)- min(x))
    return(y)
  }
  
  #Pre-process tables 
  table.p<-table.p[,interest.columns, with=F] 
  if (normalize==T) {
    table.p$SCORE<-normalize.vector(table.p[,2, with=F])
    table.p<-table.p[,c(1,3), with=F]
  }
  setnames(table.p, colnames(table.p), c("Hugo_Symbol", "SCORE")) #[Hugo_Symbol, SCORE]
  table.p$SCORE<-round(table.p$SCORE,digits=10)
  
  table.1<-readRDS(processed.table.1.rds)
  table.1<-table.1$table.1[, 1:2, with=F] #[Tumor_Sample_Barcode, Hugo_Symbol]
  
  #Order j.table by score to calculate cummulative enrichemnt [Hugo_Symbol, SCORE]
  table.p<-table.p[order(SCORE, decreasing=T),]
  
  #CALCULATIONS
  starting.patients<-unique(as.vector(table.1$Tumor_Sample_Barcode)) #Patients fed to model
  target.patients<-unique(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p$Hugo_Symbol),]$Tumor_Sample_Barcode)) #Patients we can actually cover

  #Calculate Recall
  #BREAKS<-seq(max(as.vector(table.p$SCORE)), min(as.vector(table.p$SCORE)), -0.01)
  BREAKS<-sort(unique(as.vector(table.p$SCORE)), decreasing=T)
  
  CUM.RECALL<-sapply(BREAKS, function (x) 
      length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[SCORE>=x,]$Hugo_Symbol),]$Tumor_Sample_Barcode),target.patients))/length(target.patients))
  
  #Calculate coverages
  CUM.PATIENT.COVERAGE<-sapply(1:nrow(table.p), function (x)
    length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[1:x,]$Hugo_Symbol),]$Tumor_Sample_Barcode), target.patients))/length(target.patients))
    
  NON.CUM.PATIENT.COVERAGE<-sapply(as.vector(table.p$Hugo_Symbol), function(z) 
    length(intersect(as.vector(table.1[Hugo_Symbol==z,]$Tumor_Sample_Barcode), target.patients)) /length(target.patients))
  
  #Calculate precision - ONLY RECALL FOR NOW, MAY HAVE TO USE 1000G or OTHER CANCERS COMBINED FOR THIS
  #CUM.PRECISION<-sapply(BREAKS, function (x) 
  #  length(intersect(as.vector(table.1[Hugo_Symbol %in% as.vector(table.p[SCORE>=x,]$Hugo_Symbol)]$Tumor_Sample_Barcode),target.patients))/ 
  #  length(intersect(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol)),target.in.p))/ length(unique(as.vector(table.p[SCORE>=x,]$Hugo_Symbol))))
  
  #Construct tables
  CUM.TABLE<-data.table(BREAKS=BREAKS, PATIENT.RECALL=CUM.RECALL)
  table.p$CUM.PATIENT.COVERAGE<-CUM.PATIENT.COVERAGE
  table.p$NON.CUM.PATIENT.COVERAGE<-NON.CUM.PATIENT.COVERAGE
  
  #Return as list
  return(list(PR.TABLE=CUM.TABLE,COVERAGE.TABLE=table.p, STARTING.PATIENTS=starting.patients, PATIENTS.IN.MODEL=target.patients))
}

Table.BRCA.v.rank.patient.enrich.p<-Function.p.rank.patient.enrichment(Table.v.BRCA.p, c("Hugo_Symbol", "v.PROTEIN"), 
                                                                       "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds", normalize=F)

ggplot(Table.BRCA.v.rank.patient.enrich.p$PR.TABLE, aes(BREAKS, PATIENT.RECALL)) + geom_point() +theme.format
multiplot(ggplot(Table.BRCA.v.rank.patient.enrich.p$COVERAGE.TABLE, aes(SCORE, CUM.PATIENT.COVERAGE)) + geom_point() + theme.format,
          ggplot(Table.BRCA.v.rank.patient.enrich.p$COVERAGE.TABLE, aes(SCORE, NON.CUM.PATIENT.COVERAGE)) + geom_point() + theme.format, 
          cols=1)

Table.BRCA.u.rank.patient.enrich.p<-Function.p.rank.patient.enrichment(Table.u.BRCA.p$STATS, c("Hugo_Symbol", "paired.t.stat"), 
                                                                       "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds", normalize=T)
ggplot(Table.BRCA.u.rank.patient.enrich.p$PR.TABLE, aes(BREAKS, PATIENT.RECALL)) + geom_point() +theme.format
multiplot(ggplot(Table.BRCA.u.rank.patient.enrich.p$COVERAGE.TABLE, aes(SCORE, CUM.PATIENT.COVERAGE)) + geom_point() + theme.format,
          ggplot(Table.BRCA.u.rank.patient.enrich.p$COVERAGE.TABLE, aes(SCORE, NON.CUM.PATIENT.COVERAGE)) + geom_point() + theme.format, 
          cols=1)

Table.BRCA.u.rank.patient.enrich.p.2<-Function.p.rank.patient.enrichment(Table.u.BRCA.p$STATS, c("Hugo_Symbol", "wilcoxon.stat"), 
                                                                       "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061414_Table1_BRCA.rds", normalize=T)
ggplot(Table.BRCA.u.rank.patient.enrich.p.2$PR.TABLE, aes(BREAKS, PATIENT.RECALL)) + geom_point() +theme.format
multiplot(ggplot(Table.BRCA.u.rank.patient.enrich.p.2$COVERAGE.TABLE, aes(SCORE, CUM.PATIENT.COVERAGE)) + geom_point() + theme.format,
          ggplot(Table.BRCA.u.rank.patient.enrich.p.2$COVERAGE.TABLE, aes(SCORE, NON.CUM.PATIENT.COVERAGE)) + geom_point() + theme.format, 
          cols=1)

Table.BRCA.u.rank.patient.enrich.p.2$COVERAGE.TABLE[CUM.PATIENT.COVERAGE<1,]
Table.BRCA.u.rank.patient.enrich.p$COVERAGE.TABLE[CUM.PATIENT.COVERAGE<1,]
head(Table.u.BRCA.p$STATS[order(wilcoxon.stat,decreasing=T),],10)
head(Table.u.BRCA.p$STATS[order(paired.t.stat,decreasing=T),],20)