#######042214#######
#Pipeline for Sub-Aim 2 to build network using processed tables

library(data.table)
library(limma)
library(reshape2)
library(reshape)
library(ggplot2)
library(gplots)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(Hmisc)

normalize.vector<-function(x){
  y=(x-min(x))/(max(x)- min(x))
  return(y)
}

#Load processed tables - Built in 032614_Tables.R
load("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042214_MTables.RData")

##Obtain u(j) using mutation data##
#This node metric will be obtained using processed Tables 2 and 1 
BRCA.table.2.processed
BRCA.table.1.processed

colnames(BRCA.table.2.processed)<-c("METABOLITE", "KEGG_ID", "Hugo_Symbol") #To match table.1 #DIFFERENT
Table.2.split<-split(BRCA.table.2.processed, BRCA.table.2.processed$METABOLITE,drop=T) #Make sure to drop unused levels or droplevels() in whole data frame
Table.1.split<-split(BRCA.table.1.processed, BRCA.table.1.processed$Tumor_Sample_Barcode, drop=T)

#Loop through metabolites and get median coverage across patient sample
#Use this loop to create a function for calculating Table.u on 042414
Table.u<-data.frame(METABOLITE=c(), KEGG_ID=c(), U.METABOLITE=c() ) #To store u(j)
Table.u.average<-data.frame(METABOLITE=c(), KEGG_ID=c(), U.METABOLITE=c() )

for (metabolite in names(Table.2.split)) {
  metabolite.genes<-length(as.vector(Table.2.split[[metabolite]]$Hugo_Symbol))
  u.metabolite<-c()
  
  for (sample in names(Table.1.split)) {
    
    #Find if sample has mutations in genes associated with metabolite
    Common.genes<-length(intersect(as.vector(Table.2.split[[metabolite]]$Hugo_Symbol) , as.vector(Table.1.split[[sample]]$Hugo_Symbol)))
    
    #Get sample coverage for metabolite
    sample.coverage<-Common.genes/metabolite.genes
    
    #Append for u(j) calculation
    u.metabolite<-c(u.metabolite, sample.coverage)
    
  }
  #Add u(metabolite) to table
  u.metabolite.median<-mean(u.metabolite)
  
  Table.u.average<-rbind(Table.u.average, data.frame(METABOLITE=metabolite, 
                                     KEGG_ID=unique(Table.2.split[[metabolite]]$KEGG_ID),
                                     U.METABOLITE=u.metabolite.median))
}

##Plot u(j) values 
Table.u<-as.data.table(Table.u)
Table.u<-Table.u[order(U.METABOLITE, decreasing=T),]
head(Table.u,10) #With medium calculation only 7 metabolites have a value greater than 0
ggplot(Table.u, aes(x=KEGG_ID, y=U.METABOLITE)) + geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="u(j) - Median") + xlab("Metabolites") + ylab("Median across samples")

Table.u.average<-as.data.table(Table.u.average)
Table.u.average<-Table.u.average[order(U.METABOLITE, decreasing=T),]
head(Table.u.average,20) 
nrow(Table.u.average[U.METABOLITE!=0,]) #With average calculation 765 metabolites have a value greater than zero
Table.u.average$KEGG_ID<-factor(x=Table.u.average$KEGG_ID,levels=Table.u.average$KEGG_ID,ordered=T)
ggplot(Table.u.average, aes(x=KEGG_ID, y=U.METABOLITE)) + geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="u(j) - Average") + xlab("Metabolites") + ylab("Average across samples")

save(Table.u, file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042414_Table.u.RData")
load("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042414_Table.u.RData")

##Obtain v(j) using expression data##

#Get expression samples
expression.samples<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/FILE_SAMPLE_MAP.txt", header=T, sep="\t"))
head(expression.samples)
dim(expression.samples)
expression.samples$type<-substr(expression.samples$barcode.s., 14,15)
expression.samples<-expression.samples[type=="01",] #Remove normal and metastatic, keep tumor samples only, left with 526 tumor samples
expression.samples$Sample<-substr(expression.samples$barcode.s., 1,16)
length(unique(expression.samples$Tumor_Sample_Barcode))

#Get samples that have both mutation and expression data
BRCA.table.1.processed$Sample<-substr(BRCA.table.1.processed$Tumor_Sample_Barcode, 1,16)
v.samples<-unique(intersect(BRCA.table.1.processed$Sample, expression.samples$Sample)) #This is just to know how many
length(v.samples) #505 BRCA samples with both mutation and expression data

#Get expression filenames that correspond to these samples
expression.samples<-expression.samples[Sample %in% v.samples,] #Not needed
BRCA.table.1.v<-BRCA.table.1.processed[Sample %in% v.samples,] #Not needed
dim(BRCA.table.1.v)
BRCA.table.1.v<-merge(BRCA.table.1.v, expression.samples[,c(1,4),with=F], by="Sample")

#Load Table.3 with KEGG enzymes to product
BRCA.table.3.processed

#Calculate by looping through each metabolite in Table.u by KEGG_IDs - Keep in mind that there are duplicate KEGG_IDs, but we need them to get enzymes
read.expression.files = function(path="DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3", files) {
  #There are "null" in table that function reads as NA instead to deal with easier
  tables = lapply(files, function(f) read.csv(paste0(path, "/", strsplit(f,"[$]")[[1]][2]), header=FALSE, sep="\t", skip=2, na.strings=c("null"),
                                              col.names=c("gene", strsplit(f,"[$]")[[1]][1]   )))
  
  m = tables[[1]]
  for (i in 2:length(tables)) {
    m = merge(m, tables[[i]], all=FALSE) #So that only common genes in tables are merged
  }
  rownames(m) = m$gene
  m$gene = NULL
  
  m<-as.matrix(m)
  m<-m[complete.cases(m),]
  
  return (m)
}

BRCA.table.1.v2<-BRCA.table.1.v[,c(1,4),with=F] #To introduce patient identifier into function as column name (look at function above)
BRCA.table.1.v2$code<-paste(BRCA.table.1.v2$Sample, BRCA.table.1.v2$filename, sep="$")
Table.expression<-read.expression.files(files=unique(BRCA.table.1.v2$code))
save(Table.expression, file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042414_Table.expression.RData") #SAVED EXPRESSION DATA
load("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042414_Table.expression.RData")
dim(Table.expression)
Table.expression[1:3,1:3]

Table.v<-data.frame(KEGG_ID=c(), V.METABOLITE=c()) #To store v(j)
for (metabolite in unique(as.vector(Table.u$KEGG_ID))){
  print (metabolite)
  
  #Get enzymes that produce metabolite
  metabolite.enzymes<-unique(as.vector(BRCA.table.3.processed[KEGG_ID==metabolite,]$Enzyme))
  
  #Separate cancer samples depending on wether they have mutated enzyme or not
  G.groups<-copy(BRCA.table.1.v)
  G.groups$G<-BRCA.table.1.v$Hugo_Symbol %in% metabolite.enzymes
  G1.samples<-unique(as.vector(G.groups[G==TRUE,]$Sample)) #Samples that contain mutation in enzyme
  G0.samples<-setdiff(unique(as.vector(G.groups$Sample)), G1.samples) #Rest of files that do not have a mutation for enzyme
  
  #Have to account for fact that there might be only one sample with mutation, differential expression does not work for only one sample
  if (length(G1.samples)>1) {
  
    #Transform IDs so they can be read by Table.expression syntax ("." separated rather than "-" separated)
    G1.samples<-unlist(lapply(G1.samples, function(g) paste(strsplit(g,"-")[[1]],collapse="." ) ) )
    G0.samples<-unlist(lapply(G0.samples, function(g) paste(strsplit(g,"-")[[1]],collapse="." ) ) )
    
    #Do differential expression using already read files in Table.expression
    G1.expression<-Table.expression[, G1.samples]
    G0.expression<-Table.expression[, G0.samples]
    
    colnames(G1.expression) = paste0(colnames(G1.expression), ".G1")
    colnames(G0.expression) = paste0(colnames(G0.expression), ".G0")
    
    G.all = cbind(G1.expression, G0.expression[rownames(G1.expression), ]) #combined to single matrix
    G1.n.samples<-length(colnames(G1.expression))
    G0.n.samples<-length(colnames(G0.expression))
    G.design.matrix = data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
    
    G.all<-G.all[complete.cases(G.all),] #Remove NAs
    G.all<-normalizeBetweenArrays(G.all, method="quantile") #Normalize across arrays
    
    G.fit = lmFit(G.all, model.matrix(~ G, G.design.matrix)) #fitting data 
    G.eb = eBayes(G.fit)
    
    #Get v(j) by dividing that were differentially expressed over all genes tested
    all.G.fit<-topTable(G.eb, coef=2, n=Inf)
    G.diff.exp.genes<-nrow(all.G.fit[all.G.fit$adj.P.Val<0.05,]) #Value set at p<0.05 for bonferroni corrected p-values
    
    v.metabolite<-G.diff.exp.genes/nrow(all.G.fit)
  
  } else #Otherwise if only one sample has mutation, then phenotypic effect is 0
    v.metabolite<-0
  
  #Add v(metabolite) to table
  Table.v<-rbind(Table.v, data.frame(KEGG_ID=metabolite, V.METABOLITE=v.metabolite))
}

##Plot v(j) values
Table.v<-as.data.table(Table.v)
Table.v<-Table.v[order(V.METABOLITE, decreasing=T),]
save(file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042414_Table.v.RData", Table.v)

dim(Table.v[V.METABOLITE>0,]) #380 have v(j) greater than 0
Table.v$KEGG_ID<-factor(x=Table.v$KEGG_ID, levels=Table.v$KEGG_ID,ordered=T)
ggplot(Table.v, aes(x=KEGG_ID, y=V.METABOLITE)) + geom_histogram() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="v(j)") + xlab("Metabolites") + ylab("Proportion of Differentially Expressed Genes")

ggplot(Table.v, aes(x="METABOLITES", y=V.METABOLITE)) + geom_boxplot() + scale_y_log10() +
  theme(title=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(2.5)),
        axis.title.x=element_text(size=18), axis.title.y=element_text(size=18)) +
  labs(title="v(j)") + xlab("Metabolites") + ylab("Proportion of Differentially Expressed Genes")

##Obtain W(j) using expression data##
#Run on PC
Table.expression[1:3,1:3]

Network.Kegg.IDs<-unique(as.vector(Table.u$KEGG_ID))
Table.W<-data.frame(KEGG.ID.1=c(), KEGG.ID.2=c(), W.METABOLITE=c())
for (metabolite in Network.Kegg.IDs){
  print (metabolite)
  
  #Remove current metabolite
  Network.Kegg.IDs<-Network.Kegg.IDs[Network.Kegg.IDs!=metabolite]
  
  #Obtain Average enzyme expression values across patient samples for metabolite
  Current.enzymes<-unique(as.vector(BRCA.table.3.processed[KEGG_ID==metabolite,]$Enzyme))
  
  #Check that current.enzymes are present in expression matrix  
  if (length(intersect(Current.enzymes, rownames(Table.expression)))>0) {
  
    Current.averages<-as.vector(colMeans(as.data.frame(Table.expression)[Current.enzymes,],na.rm=T))
    
    #Now that it is remove loop through others and obtain value per pair
    for (other.metabolite in Network.Kegg.IDs) {
      
      #Obtain average enzyme expression values across patient samples for other metabolite
      Other.enzymes<-unique(as.vector(BRCA.table.3.processed[KEGG_ID==other.metabolite,]$Enzyme))
      
      #Check that enzymes for other.metabolite are in expression matrix
      if (length(intersect(Other.enzymes, rownames(Table.expression)))>0) {
      
        Other.averages<-as.vector(colMeans(as.data.frame(Table.expression)[Other.enzymes,],na.rm=T))
        
        #Get correlation as absolute value of spearman - MAY CHANGE TO PEARSON LATER!
        spearman.met<-abs(cor(Current.averages, Other.averages,method="spearman"))
        
        #Store in table
        Table.W<-rbind(Table.W, data.frame(KEGG.ID.1=metabolite, KEGG.ID.2=other.metabolite, W.METABOLITE=spearman.met))
      
      } else #If other is not present then the spearman.met is 0, other.metabolite will not be connected to network (Unless node weight contributes to it)
        Table.W<-rbind(Table.W, data.frame(KEGG.ID.1=metabolite, KEGG.ID.2=other.metabolite, W.METABOLITE=0))
    }
    
  } else #If enzymes not present in expression data then metabolite will not be connected to network, give edges of 0 to it (Unless node weight contributes to it)
    Table.W<-rbind(Table.W, data.frame(KEGG.ID.1=rep(metabolite, length(Network.Kegg.IDs)),
                                       KEGG.ID.2=Network.Kegg.IDs,
                                       W.METABOLITE=rep(0, length(Network.Kegg.IDs))))
}

######042314########

##Modify u(j) calculation to be normalized by protein length, quantile approach##

#Load table of lengths
Table.length<-as.data.table(read.csv("DATABASES/UNIPROT/042314_GENE_LENGTH", header=T, sep="\t"))

#Get protein lengths of genes we are interested on
BRCA.table.2.processed
length(unique(BRCA.table.2.processed$Hugo_Symbol)) #3896 gene in Table.2
length(intersect(unique(Table.length$GENE), unique(BRCA.table.2.processed$Hugo_Symbol))) #We have gene length information for all of them
BRCA.table.2.length<-Table.length[GENE %in% unique(BRCA.table.2.processed$Hugo_Symbol),]

#Plot distribution of amino acid lengths for proteins in  HMDB(j) 
ggplot(BRCA.table.2.length, aes(x="GENE",y=LENGTH)) + geom_boxplot() + geom_jitter(size=1.0,aes(colour=Q.Rank)) + coord_flip() +
  theme(axis.text.x=element_text(size=rel(2.5)),
        axis.title.x=element_text(size=22))

ggplot(BRCA.table.2.length, aes(x=LENGTH, colour=as.factor(Q.Rank))) + geom_histogram(binwidth=40) +
  theme(axis.text.y=element_text(size=rel(2.5)), axis.text.x=element_text(size=rel(2.5)),
        axis.title.y=element_text(size=22), axis.title.x=element_text(size=22))

#Get quartiles of HMDB(j) protein lengths and assign quartile rank to gene lengths
quantile(BRCA.table.2.length$LENGTH)
quantile(BRCA.table.2.length$LENGTH, c(0,0.333,0.666,1))

BRCA.table.2.length$Q.Rank<-as.numeric(cut(BRCA.table.2.length$LENGTH, breaks=quantile(BRCA.table.2.length$LENGTH, c(0,0.333,0.666,1)), include.lowest=T,labels=1:3))


#Created function to re-calculate u(j) with and without ranks (as before)
Table.u.function<-function (Table.2.s, Table.1.s, method="mean", Q.rank.table=NULL) {
  #Methods can be: mean, median, mean.rank, median.rank. Default is mean
  #To use*.rank a Q.rank.table must be provided such as BRCA.table.2.length
  
  #To keep track
  current.count<-0
  max.count<-length(Table.2.s)
  
  #To store result
  Table.u.result<-data.frame(METABOLITE=c(), KEGG_ID=c(), U.METABOLITE=c() ) #To store u(j)
  
  #Loop through metabolite and patients in split tables
  for (metabolite in names(Table.2.s)) {
    print (current.count/max.count)
    
    metabolite.genes<-as.vector(Table.2.s[[metabolite]]$Hugo_Symbol)
    u.metabolite<-c()
    
    for (sample in names(Table.1.s)) {
      
      #Find if sample has mutations in genes associated with metabolite
      #This step is the same whether we are doing ranked or unranked calculation
      Common.genes<-intersect(as.vector(Table.2.s[[metabolite]]$Hugo_Symbol) , as.vector(Table.1.s[[sample]]$Hugo_Symbol))
      
      #Depending on Ranked or Unranked calculation
      if ( method %in% c("mean", "median")) {
        
        #Get sample coverage for metabolite
        sample.coverage<-length(Common.genes)/length(metabolite.genes) #If there are no common.genes then this will be zero
        
        #Append for u(j) calculation
        u.metabolite<-c(u.metabolite, sample.coverage)
        
      } else if (method %in% c("mean.rank", "median.rank")) {
        
        #Use quantile ranks to get normalized sample coverage
        gene.qranks<-Q.rank.table[Q.rank.table$GENE %in% Common.genes,]$Q.Rank #We get as many qranks as there are genes
        qrank.calc<-sum(rep(1,length(Common.genes))/gene.qranks) #Normalize each mutated gene count (1) by its rank weight and sum over them to get counts
                                                                  #This calculation is analogus to length (Common genes), but normalized by rank
        sample.coverage<-qrank.calc/length(metabolite.genes) #the length of metabolite genes should be normalized as well
        
        #Append for u(j) calculation
        u.metabolite<-c(u.metabolite, sample.coverage)
        
      } else
        print ("Error, please choose appropriate method")
      
    }
    
    #Add u(metabolite) to table depending on method
    if (method %in% c("mean", "mean.rank")) {
      
      u.metabolite.calc<-mean(u.metabolite)
    } else if (method %in% c("median", "median.rank")) {
      
      u.metabolite.calc<-median(u.metabolite)
    } else
      print ("Error with method")
    
    Table.u.result<-rbind(Table.u.result, data.frame(METABOLITE=metabolite, 
                                                     KEGG_ID=unique(Table.2.s[[metabolite]]$KEGG_ID),
                                                     U.METABOLITE=u.metabolite.calc))
    
    current.count<-current.count+1
  }
  return(Table.u.result)
}

#Test function, works as before on unranked situation
Table.u.median<-Table.u.function(Table.2.split, Table.1.split ,method="median")
head(Table.u.median[order(Table.u.median$U.METABOLITE,decreasing=T),],10)

#Use to obtain Table.u using quantile ranked
Table.u.mean.ranked<-Table.u.function(Table.2.split, Table.1.split, method="mean.rank", Q.rank.table=BRCA.table.2.length)
Table.u.mean.ranked<-as.data.table(Table.u.mean.ranked)
Table.u.mean.ranked<-Table.u.mean.ranked[order(U.METABOLITE,decreasing=T),]
Table.u.mean.ranked
save(file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042614_Table.u.mean.ranked.RData", Table.u.mean.ranked)
dim(Table.u.mean.ranked)
dim(Table.u.mean.ranked[U.METABOLITE!=0,])

Table.u.median.ranked<-Table.u.function(Table.2.split, Table.1.split, method="median.rank", Q.rank.table=BRCA.table.2.length)
Table.u.median.ranked<-as.data.table(Table.u.median.ranked)
head(Table.u.median.ranked[order(U.METABOLITE, decreasing=T),],8)

########0042514########

##Compare effect of quartile normalized on u(j)##
theme.format<-theme(axis.text.y=element_text(size=rel(2.5)), axis.text.x=element_text(size=rel(2.5)),
                    axis.title.y=element_text(size=22), axis.title.x=element_text(size=22),
                    legend.text = element_text(size = 22))

ggplot(Table.u.mean.ranked, aes(x=U.METABOLITE)) + geom_histogram() + theme.format
ggplot(Table.u.average, aes(x=U.METABOLITE)) + geom_histogram() + theme.format
ggplot(melt(merge(Table.u.average,Table.u.mean.ranked, by=c("KEGG_ID","METABOLITE")), measure.vars=c("U.METABOLITE.x", "U.METABOLITE.y")),
            aes(x=value, fill=variable)) + geom_histogram() + theme.format

#Merge and normalized between 0-1 for fair comparisson
Table.u.mean.compare<-merge(Table.u.average, Table.u.mean.ranked, by=c("METABOLITE", "KEGG_ID"))
colnames(Table.u.mean.compare)<-c("METABOLITE", "KEGG_ID", "U.METABOLITE.MEAN", "U.METABOLITE.MEAN.RANKED")

Table.u.mean.compare$U.METABOLITE.MEAN<-normalize.vector(Table.u.mean.compare$U.METABOLITE.MEAN)
Table.u.mean.compare$U.METABOLITE.MEAN.RANKED<-normalize.vector(Table.u.mean.compare$U.METABOLITE.MEAN.RANKED)

#Plot Table.u.mean vs Table.u.mean.rank normalized for comparisson
ggplot(melt(Table.u.mean.compare, measure.vars=c("U.METABOLITE.MEAN","U.METABOLITE.MEAN.RANKED"))
       , aes(x=value, fill=variable)) + geom_histogram(binwidth=0.05) + theme.format

head(melt(Table.u.mean.compare,measure.vars=c("U.METABOLITE.MEAN","U.METABOLITE.MEAN.RANKED")))
t.test(Table.u.mean.compare$U.METABOLITE.MEAN, Table.u.mean.compare$U.METABOLITE.MEAN.RANKED, paired=T) #P<0.05, ranked makes significant difference in u(j)
t.test(merge(Table.u.average,Table.u.mean.ranked, by=c("KEGG_ID","METABOLITE"))$U.METABOLITE.x, 
       merge(Table.u.average,Table.u.mean.ranked, by=c("KEGG_ID","METABOLITE"))$U.METABOLITE.y,
       paired=T)

##Obtain V(j) using v(j) and u(j)##

Table.u.mean.ranked[duplicated(KEGG_ID),]
Table.u.mean.ranked[KEGG_ID %in% Table.u.mean.ranked[duplicated(KEGG_ID),]$KEGG_ID,][order(KEGG_ID),]
Table.v
dim(Table.v[V.METABOLITE!=0,])
ggplot(Table.v, aes(x=normalize.vector(V.METABOLITE))) + geom_histogram(binwidth=0.002)

#Make function to get weight node V(j)
Table.V.function<-function (u.table, v.table, beta){
  
  #Merge them by KEGG_ID
  Table.V<-merge(u.table, v.table, by="KEGG_ID")
  
  #Normalize each u(j) and v(j) column before calculation
  Table.V$U.METABOLITE<-normalize.vector(Table.V$U.METABOLITE)
  Table.V$V.METABOLITE<-normalize.vector(Table.V$V.METABOLITE)
  
  #Assign a weight dependency factor based on comparisson of u(j) vs v(j)
  Table.V$WEIGHT.FACTOR<-ifelse(Table.V$U.METABOLITE>Table.V$V.METABOLITE, "u(j)",
                                ifelse(Table.V$U.METABOLITE<Table.V$V.METABOLITE, "v(j)",
                                       "None"))
  
  #Merge by beta - this will be in the form V=beta*u+(1-beta)*v
  Table.V$NODE.WEIGHT<-Table.V$U.METABOLITE*beta + Table.V$V.METABOLITE*(1-beta)
  
  #Return values
  return(Table.V)
}

#Use function to obtain V(j) at different betas - THIS IS USING THE u.mean.ranked
Table.V.0.0<-Table.V.function(Table.u.mean.ranked, Table.v, 0.0)
Table.V.0.1<-Table.V.function(Table.u.mean.ranked, Table.v, 0.1)
Table.V.0.2<-Table.V.function(Table.u.mean.ranked, Table.v, 0.2)
Table.V.0.3<-Table.V.function(Table.u.mean.ranked, Table.v, 0.3)
Table.V.0.4<-Table.V.function(Table.u.mean.ranked, Table.v, 0.4)
Table.V.0.5<-Table.V.function(Table.u.mean.ranked, Table.v, 0.5)
Table.V.0.6<-Table.V.function(Table.u.mean.ranked, Table.v, 0.6)
Table.V.0.7<-Table.V.function(Table.u.mean.ranked, Table.v, 0.7)
Table.V.0.8<-Table.V.function(Table.u.mean.ranked, Table.v, 0.8)
Table.V.0.9<-Table.V.function(Table.u.mean.ranked, Table.v, 0.9)
Table.V.1.0<-Table.V.function(Table.u.mean.ranked, Table.v, 1.0)

#Plot V(j) distributions as function of beta
Beta.V.plot.1<-data.frame(Beta=c(), NODE.WEIGHT=c(), WEIGHT.FACTOR=c())
for (beta in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  V.TABLE<-get(paste0("Table.V.", beta))
  Beta.V.plot.1<-rbind(Beta.V.plot.1, data.frame(Beta=beta, NODE.WEIGHT=V.TABLE$NODE.WEIGHT,
                                                 WEIGHT.FACTOR=V.TABLE$WEIGHT.FACTOR))
}
head(Beta.V.plot.1)

ggplot(Beta.V.plot.1, aes(x=factor(Beta), y=NODE.WEIGHT)) + geom_boxplot() + theme.format + geom_jitter(size=1.5, aes(colour=WEIGHT.FACTOR)) +
  scale_colour_brewer(palette="Spectral")

##Analyze W(j) obtained from PC##
load("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042514_Table.W.RData")
head(Table.W)
dim(Table.W)
Table.W[is.na(Table.W)]<-0 #HAD TO TAKE CARE OF NA BECAUSE SOME OF THE VALUES IN THE EXPRESSION MATRIX HAD NAs in THEM

ggplot(Table.W, aes(x=W.METABOLITE)) +geom_histogram() + theme.format
ggplot(Table.W, aes(x="METABOLITES",y=W.METABOLITE)) +geom_boxplot() + theme.format

########042614########

##Integrate Node Weigths V(j) into Edge Weights W(j)##

#Write function for integration
network.integration.V.W<-function(V.table.beta, w.table) {
  #Function takes V.table.beta (don't confuse with v.table) table that have V(j) beta assigend node weights
  
  dummy<-V.table.beta[,c("KEGG_ID", "NODE.WEIGHT"), with=F]
  
  #Remove duplicates in V(j) if necessary
  dummy<-dummy[!duplicated(dummy$KEGG_ID)]
  
  #Do first node assignment
  colnames(dummy)<-c("KEGG.ID.1", "NODE.WEIGHT.1")
  v.w.table<-merge(x=dummy, y=w.table, by="KEGG.ID.1")
  
  #Do second node assignment
  colnames(dummy)<-c("KEGG.ID.2", "NODE.WEIGHT.2")
  v.w.table<-merge(x=v.w.table, y=dummy, by="KEGG.ID.2")
  
  #Get integrated edge weight
  v.w.table$INTEGRATED.W<-v.w.table$NODE.WEIGHT.1 + v.w.table$NODE.WEIGHT.2 + v.w.table$W.METABOLITE
  
  #Return
  return(v.w.table)
}

#Get integrate w.nodes files for different Betas
Table.W.V.00<-network.integration.V.W(Table.V.0.0, Table.W)
Table.W.V.01<-network.integration.V.W(Table.V.0.1, Table.W)
Table.W.V.02<-network.integration.V.W(Table.V.0.2, Table.W)
Table.W.V.03<-network.integration.V.W(Table.V.0.3, Table.W)
Table.W.V.04<-network.integration.V.W(Table.V.0.4, Table.W)
Table.W.V.05<-network.integration.V.W(Table.V.0.5, Table.W)
Table.W.V.06<-network.integration.V.W(Table.V.0.6, Table.W)
Table.W.V.07<-network.integration.V.W(Table.V.0.7, Table.W)
Table.W.V.08<-network.integration.V.W(Table.V.0.8, Table.W)
Table.W.V.09<-network.integration.V.W(Table.V.0.9, Table.W)
Table.W.V.10<-network.integration.V.W(Table.V.1.0, Table.W)

#Save integrated w.node files - keep in mind INTEGRATED.W here is unnormalized
save(list=c(paste0("Table.W.V.0", 0:9), "Table.W.V.10"), file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042814_Tables.W.V.RData")

#Plot W.V(j) distributions (integrated edge weight) as function of beta
Beta.W.V.plot.1<-data.frame(Beta=c(), INTEGRATED.W=c(), WEIGHT.FACTOR=c())
for (beta in c(paste0("0",seq(0,9)),"10")) {
  W.V.TABLE<-get(paste0("Table.W.V.", beta))
  W.V.TABLE$WEIGHT.FACTOR<-ifelse( (W.V.TABLE$NODE.WEIGHT.1 + W.V.TABLE$NODE.WEIGHT.2)>W.V.TABLE$W.METABOLITE, "V(j)",
                                   ifelse( (W.V.TABLE$NODE.WEIGHT.1 + W.V.TABLE$NODE.WEIGHT.2)<W.V.TABLE$W.METABOLITE,"W(j)",
                                           ifelse( (W.V.TABLE$NODE.WEIGHT.1 + W.V.TABLE$NODE.WEIGHT.2)==W.V.TABLE$W.METABOLITE, "Both",
                                                   "None")))
  
  Beta.W.V.plot.1<-rbind(Beta.W.V.plot.1, data.frame(Beta=beta, INTEGRATED.W=W.V.TABLE$INTEGRATED.W,
                                                 WEIGHT.FACTOR=W.V.TABLE$WEIGHT.FACTOR))
}
head(Beta.W.V.plot.1) #I believe there is no much use of knowing the weight factor(if V is greater than W)
Beta.W.V.plot.1<-as.data.table(Beta.W.V.plot.1)
Beta.W.V.plot.1$SPLIT.BETA<-paste0(Beta.W.V.plot.1$Beta, Beta.W.V.plot.1$WEIGHT.FACTOR)
Beta.W.V.plot.1$Beta<- factor(paste0(substr(Beta.W.V.plot.1$Beta,1,1 ),".",substr(Beta.W.V.plot.1$Beta,2,2))  )

ggplot(Beta.W.V.plot.1, aes(x=Beta, y=INTEGRATED.W)) + geom_boxplot() + theme.format + 
  scale_colour_brewer(palette="Spectral") #See the effect of V(j) on Total edge weight, since W(j) does not depend on Beta but V(j) does

ggplot(Beta.W.V.plot.1, aes(x=factor(Beta), y=INTEGRATED.W)) + geom_violin() + theme.format + 
  scale_colour_brewer(palette="Spectral") + geom_jitter()

ggplot(Beta.W.V.plot.1, aes(x=INTEGRATED.W)) + geom_histogram() + theme.format +
  facet_wrap(~Beta) + theme(strip.text.x = element_text(size = 18, colour = "black"))

#Write tables.tsf for cytoscape analysis
for (i in c(paste0("0",seq(0,9)),"10")) {
  file.name=paste0("NETWORK/TSF_FILES/042614_WV",i,".csv")
  name.table<-get(paste0("Table.W.V.",i))
  write.csv(file=file.name, name.table, sep="\t", quote=F, row.names=F)
}

#Write tables for SPICi cluster analysis
  #Unnormalized
for (i in 1:9) {
  file.name=paste0("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/042614_NETWORK_MEAN_RANKED_0",i)
  name.table<-get(paste0("Table.W.V.0",i))[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  write.table(file=file.name, name.table, sep="\t", quote=F, row.names=F, col.names=F)
}

  #Normalized
for (i in c(paste0("0",seq(0,9)),"10")) {
  file.name=paste0("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/042614_NETWORK_MEAN_RANKED_",i,"N")
  name.table<-get(paste0("Table.W.V.",i))[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  name.table$INTEGRATED.W<-normalize.vector(name.table$INTEGRATED.W) #Normalized step vs above unnormalized
  write.table(file=file.name, name.table, sep="\t", quote=F, row.names=F, col.names=F)
}

#######042714#######

##Analyze obtained clusters##
cluster.analysis.function<-function(folder, cluster.name.pattern){
  
  #Get cluster files
  cluster.files<-list.files(folder, cluster.name.pattern)
  
  #Create table to store cluster counts per file
  cluster.counts<-data.frame(FILENAME=c(), CLUSTER.SIZE=c())
  
  #Get clusters for each file
  for (c.file in cluster.files) {
    dummy.cluster.file<-read.csv(paste0(folder,c.file), sep="\t", header=F,stringsAsFactors=F)
    
    #Count number of metabolites per cluster
    for (cluster in 1:nrow(dummy.cluster.file)) {
      dummy.vector<-as.character(as.vector(dummy.cluster.file[cluster,]))
      dummy.vector<-dummy.vector[dummy.vector!=""]
      dummy.vector.size<-length(dummy.vector)
      
      #Append to cluster.table
      cluster.counts<-rbind(cluster.counts, 
                            data.frame(FILENAME=c.file, CLUSTER.SIZE=dummy.vector.size))
    }
  }
  #Assign beta
  cluster.counts<-as.data.table(cluster.counts)
  cluster.counts$BETA<-sapply(cluster.counts$FILENAME, function(d) strsplit(as.character(d),"_")[[1]][5])
  cluster.counts<-cluster.counts[order(BETA,CLUSTER.SIZE),]
  
  #Return
  return(cluster.counts)
}

#SPICi analysis of unnormalized clusters
cluster.analysis.5C<-cluster.analysis.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/", "5C")
cluster.analysis.5C<-as.data.table(cluster.analysis.5C)
cluster.analysis.5C$BETA<-sapply(cluster.analysis.5C$FILENAME, function(d) strsplit(as.character(d),"_")[[1]][5])

ggplot(cluster.analysis.5C, aes(x=BETA, y=CLUSTER.SIZE)) + geom_boxplot(position="identity") + 
  theme.format + geom_jitter()  +scale_y_log10()
ggplot(cluster.analysis.5C, aes(x=BETA,y=CLUSTER.SIZE, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity") + 
  theme.format #+scale_y_log10()

#SPICI analysis of normalized clusters
cluster.analysis.N_5C<-cluster.analysis.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/", "N_5C")
cluster.analysis.N_5C<-as.data.table(cluster.analysis.N_5C)
cluster.analysis.N_5C$BETA<-sapply(cluster.analysis.N_5C$FILENAME, function(d) strsplit(as.character(d),"_")[[1]][5])
cluster.analysis.N_5C<-cluster.analysis.N_5C[order(BETA,CLUSTER.SIZE),]

ggplot(cluster.analysis.N_5C, aes(x=BETA, y=CLUSTER.SIZE)) + geom_boxplot(position="identity") + 
  theme.format + geom_jitter(aes(colour=CLUSTER.SIZE))
ggplot(cluster.analysis.N_5C, aes(x=BETA,y=CLUSTER.SIZE, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity") + 
  theme.format 
ggplot(cluster.analysis.N_5C, aes(x=BETA,y=length(CLUSTER.SIZE)/264, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity", colour="black") + 
  theme.format +
  ylab("Number of Clusters")

#########042814##########

##Do significance test on clusters##
#Analyze clusters that were obtained from SPICi on normalized integreated edge weights

#Created function for it
cluster.analysis.p.value.function<-function(filename, BETA, background.network) {

  #Load file
  dummy.cluster.file<-read.csv(filename, sep="\t", header=F, stringsAsFactors=F)
  
  #Prep background network - KEEP IN MIND THAT THIS SHOULD BE THE NETWORK THAT PRODUCED CLUSTERS!!
  #Background network are of the form: i.e. Table.W.V.03
  dummy.background<-copy(background.network)
  dummy.background<-dummy.background[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  dummy.background<-rbind(dummy.background, dummy.background[,c(2,1,3),with=F], use.names=F)
  
  #Prep vector to store p.values
  dummy.cluster.scores<-data.frame(FILENAME=c(), Beta=c(), CLUSTER.SIZE=c(), P.VALUE=c())

  for (cluster in 1:nrow(dummy.cluster.file)) {
    dummy.vector<-as.character(as.vector(dummy.cluster.file[cluster,]))
    dummy.vector<-dummy.vector[dummy.vector!=""]
    dummy.vector.size<-length(dummy.vector)
    
    #Create table of edges for cluster to get scores
    cluster.table<-as.data.frame(t(combn(dummy.vector,2))) #No repeats
    colnames(cluster.table)<-colnames(dummy.background)[1:2] #To merge by KEGG_IDs
    
    #Merge with background to get scores
    cluster.table<-merge(cluster.table, dummy.background, by=colnames(dummy.background)[1:2])
    cluster.score<-sum(cluster.table$INTEGRATED.W)
    
    #Permutation test per score, this depends on cluster size (x1000)
    permutation.vector<-c()
    for (i in 1:1000) {
      #Get sample x times depending on cluster size
      permutation.cluster<-as.vector(sample(dummy.background[[colnames(dummy.background)[1]]],dummy.vector.size))
      
      #Get cluster score for permutation vector
      permutation.table<-as.data.frame(t(combn(permutation.cluster,2)))
      colnames(permutation.table)<-colnames(dummy.background)[1:2]
      permutation.table<-merge(permutation.table, dummy.background, by=colnames(dummy.background)[1:2])
      permutation.score<-sum(permutation.table$INTEGRATED.W)
      
      #Store permutated score
      permutation.vector<-c(permutation.vector, permutation.score)
    }
    
    #Get p-value per cluster
    dummy.p.value<-length(permutation.vector[permutation.vector>=cluster.score])/length(permutation.vector)
    
    #Store cluster information
    dummy.cluster.scores<-rbind(dummy.cluster.scores, data.frame(FILENAME=filename, Beta=BETA,
                                                                 CLUSTER.SIZE=dummy.vector.size, P.VALUE=dummy.p.value))
  }
  #Do multiple hypothesis correction across p.values
  dummy.cluster.scores$ADJ.P.VALUE<-p.adjust(dummy.cluster.scores$P.VALUE,method="fdr")
  
  #Return table
  return (dummy.cluster.scores)
}

cluster.analysis.p.value.function.V2<-function(filename, BETA, background.network) {
  
  #To keep count
  x=1
  
  #Load file
  dummy.cluster.file<-read.csv(filename, sep="\t", header=F, stringsAsFactors=F)
  
  #Prep background network - KEEP IN MIND THAT THIS SHOULD BE THE NETWORK THAT PRODUCED CLUSTERS!!
  #Background network are of the form: i.e. Table.W.V.03
  dummy.background<-copy(background.network)
  dummy.background<-dummy.background[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  
  #Double for dummy cluster, regular for permutation
  dummy.background.double<-rbind(dummy.background, dummy.background[,c(2,1,3),with=F], use.names=F)
  
  #Prep vector to store p.values
  dummy.cluster.scores<-data.frame(FILENAME=c(), Beta=c(), CLUSTER.SIZE=c(), P.VALUE=c())
  
  for (cluster in 1:nrow(dummy.cluster.file)) {
    dummy.vector<-as.character(as.vector(dummy.cluster.file[cluster,]))
    dummy.vector<-dummy.vector[dummy.vector!=""]
    dummy.vector.size<-length(dummy.vector) #Number of nodes in cluster
    
    #Create table of edges for cluster to get scores
    cluster.table<-as.data.frame(t(combn(dummy.vector,2))) #No repeats
    colnames(cluster.table)<-colnames(dummy.background.double)[1:2] #To merge by KEGG_IDs
    
    #Merge with background to get scores
    cluster.table<-merge(cluster.table, dummy.background.double, by=colnames(dummy.background.double)[1:2])
    cluster.score<-sum(cluster.table$INTEGRATED.W)
    
    #Permutation test per score, this depends on cluster size (x1000)
    dummy.vector.edges<-dummy.vector.size*(dummy.vector.size-1)/2  #Number of edges in cluster
    permutation.vector<-replicate(10000, sum(sample(dummy.background$INTEGRATED.W, dummy.vector.edges)))
    
    #Get p-value per cluster
    dummy.p.value<-length(permutation.vector[permutation.vector>=cluster.score])/length(permutation.vector)
    
    #Store cluster information
    dummy.cluster.scores<-rbind(dummy.cluster.scores, data.frame(FILENAME=filename, Beta=BETA,
                                                                 CLUSTER.SIZE=dummy.vector.size, P.VALUE=dummy.p.value))
    #To keep count
    print (x/nrow(dummy.cluster.file))
    x=x+1
    
  }
  #Do multiple hypothesis correction across p.values
  dummy.cluster.scores$ADJ.P.VALUE<-p.adjust(dummy.cluster.scores$P.VALUE,method="fdr")
  
  #Return table
  return (dummy.cluster.scores)
}

#Applied function in PC - Version 2 is much faster and uses 10000x permutations
dummy.test<-cluster.analysis.p.value.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/042614_NETWORK_MEAN_RANKED_03N_5C.cluster",
                                              BETA="0.3", Table.W.V.03)

#########042914########

##First find most common cluster across all Beta SPICi outputs##
common.cluster.spici.function<-function(folder, PATTERN) {
  
  #Prep vector to store clusters
  dummy.cluster.collapsed.table<-data.frame(FILENAME=c(), CLUSTER.SIZE=c(), CLUSTER=c())
  
  #Extract files
  dummy.files=list.files(folder,pattern=PATTERN)
  
  #Extract for each file
  for (filename in dummy.files) {
  
    #Load file
    dummy.cluster.file<-read.csv(paste0(folder,"/",filename), sep="\t", header=F, stringsAsFactors=F)
    
    for (cluster in 1:nrow(dummy.cluster.file)) {
      dummy.vector<-as.character(as.vector(dummy.cluster.file[cluster,]))
      dummy.vector<-dummy.vector[dummy.vector!=""]
      dummy.vector.size<-length(dummy.vector)
      
      #Linearized cluster elements into a string, make sure to order it first for fair comparisson
      dummy.vector<-dummy.vector[order(dummy.vector)]
      dummy.vector.collapsed<-paste(dummy.vector, collapse="_")
      
      #Append to table of Beta clusters
      dummy.cluster.collapsed.table<-rbind(dummy.cluster.collapsed.table,
                                           data.frame(FILENAME=filename, 
                                                      CLUSTER.SIZE=dummy.vector.size,
                                                      CLUSTER=dummy.vector.collapsed))
      
    }
  }
  return(dummy.cluster.collapsed.table)
}

common.clusters<-common.cluster.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "N_5C")
common.clusters$BETA<-sapply(common.clusters$FILENAME, function(x) substr(strsplit(as.character(x) ,"_")[[1]][5], 1,2))
dim(common.clusters)
colnames(common.clusters)
common.clusters.freq<-as.data.frame(table(common.clusters$CLUSTER))
common.clusters.freq<-common.clusters.freq[order(common.clusters.freq$Freq,decreasing=T),]
head(common.clusters.freq,10)
dim(common.clusters.freq)

#Plot clusters found across betas 
common.clusters.plot.1<-common.clusters[,c(3,4)]
common.clusters.plot.1$PRESENCE<-1
common.clusters.plot.1<-cast(common.clusters.plot.1, CLUSTER~BETA, value="PRESENCE",fill=0)
head(as.matrix(common.clusters.plot.1))
heatmap.2(as.matrix(common.clusters.plot.1),trace="none")

##Second find most common metabolites across all Beta SPICi outputs##
common.kegg.spici.function<-function(folder, PATTERN) {
  
  #Prep vector to store clusters
  dummy.kegg.collapsed.table<-data.frame(FILENAME=c(), KEGG=c())
  
  #Extract files
  dummy.files=list.files(folder,pattern=PATTERN)
  
  #Extract for each file
  for (filename in dummy.files) {
    
    #Load file
    dummy.cluster.file<-read.csv(paste0(folder,"/",filename), sep="\t", header=F, stringsAsFactors=F)
    
    for (cluster in 1:nrow(dummy.cluster.file)) {
      dummy.vector<-as.character(as.vector(dummy.cluster.file[cluster,]))
      dummy.vector<-dummy.vector[dummy.vector!=""]
      
      #Append to table of Beta KEGGs
      dummy.kegg.collapsed.table<-rbind(dummy.kegg.collapsed.table,
                                           data.frame(FILENAME=filename, 
                                                      KEGG=dummy.vector))
      
    }
  }
  return(dummy.kegg.collapsed.table)
}

common.kegg<-common.kegg.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "N_5C")
head(common.kegg)
common.kegg$BETA<-sapply(common.kegg$FILENAME, function(x) substr(strsplit(as.character(x) ,"_")[[1]][5], 1,2))
common.kegg.freq<-as.data.frame(table(common.kegg$KEGG))
common.kegg.freq<-common.kegg.freq[order(common.kegg.freq$Freq, decreasing=T),]
tail(common.kegg.freq,20)
dim(common.kegg.freq)

#Plot KEGGs found across betas
common.kegg.plot.1<-common.kegg[,c(2,3)]
common.kegg.plot.1$PRESENCE<-1
common.kegg.plot.1<-cast(common.kegg.plot.1, KEGG~BETA, value="PRESENCE",fill=0)
heatmap.2(as.matrix(common.kegg.plot.1), trace="none")

#####Re-calculate Integrated Node and Edge weights to account for v(j) lack of presence#####

##Modify V.funciton to use Omega factor (median(u)/median(v))
#Will not use max like in the paper because v(j) has very big outliers
#Will change this to .V3 if I ever decide to calculate v(j) in terms of cancer vs normal gene expression
median(Table.v$V.METABOLITE)#Mean is heaviliy influenced by outlier, so use median instead
median(Table.u.mean.ranked$U.METABOLITE) #Very similar to mean
Table.V.function.V2<-function (u.table, v.table, BETA){
  
  #Merge them by KEGG_ID
  Table.V<-merge(u.table, v.table, by="KEGG_ID")
  
  #APPLY MAX NORMALIZATION
  OMEGA<-median(Table.V$U.METABOLITE)/median(Table.V$V.METABOLITE)
  Table.V$V.METABOLITE<-Table.V$V.METABOLITE * OMEGA
  
  #Normalize u and v
  #This is still done due to major gap between u(j),v(j) and W(j)
  #Table.V$U.METABOLITE<-normalize.vector(Table.V$U.METABOLITE)
  #Table.V$V.METABOLITE<-normalize.vector(Table.V$V.METABOLITE)
  
  #Assign a weight dependency factor based on comparisson of u(j) vs v(j)
  Table.V$WEIGHT.FACTOR<-ifelse(Table.V$U.METABOLITE>Table.V$V.METABOLITE, "u(j)",
                                ifelse(Table.V$U.METABOLITE<Table.V$V.METABOLITE, "v(j)",
                                       "None"))
  
  #Merge by beta - this will be in the form V=beta*u+(1-beta)*v
  Table.V$NODE.WEIGHT<-Table.V$U.METABOLITE*BETA + Table.V$V.METABOLITE*(1-BETA)
  
  #Log10 normalize and remove -Inf values
  Table.V$NODE.WEIGHT<-log10(Table.V$NODE.WEIGHT)
  Table.V<-Table.V[is.finite(NODE.WEIGHT),]
  
  #Back normalize after removing -Inf value (original zero value node weight)
  Table.V$NODE.WEIGHT<-normalize.vector(Table.V$NODE.WEIGHT)
  
  #Return values
  return(Table.V)
}

#Apply to obtain new V(j)
for (beta in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  dummy.beta<-Table.V.function.V2(Table.u.mean.ranked, Table.v, BETA=as.numeric(beta))
  assign(paste0("Table.V.",beta,".V2",collapse="") , dummy.beta)
}

#Plot V(j).V2 distributions as function of beta
Beta.V2.plot.1<-data.frame(Beta=c(), NODE.WEIGHT=c(), WEIGHT.FACTOR=c())
for (beta in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  V2.TABLE<-get(paste0("Table.V.", beta,".V2"))
  Beta.V2.plot.1<-rbind(Beta.V2.plot.1, data.frame(Beta=beta, NODE.WEIGHT=V2.TABLE$NODE.WEIGHT,
                                                 WEIGHT.FACTOR=V2.TABLE$WEIGHT.FACTOR))
}
Beta.V2.plot.1<-as.data.table(Beta.V2.plot.1)
ggplot(Beta.V2.plot.1, aes(x=factor(Beta), y=NODE.WEIGHT)) + geom_boxplot() + theme.format + geom_jitter(size=1.5, aes(colour=WEIGHT.FACTOR))  # +
  scale_colour_brewer(palette="Spectral")

##Integrate V(j)V2 into W(J)##
#Modify network.integration.V.W function for Omega parameter
network.integration.V.W.V2<-function(V.table.beta, w.table) {
  #Function takes V.table.beta (don't confuse with v.table) table that have V(j) beta assigend node weights
  
  dummy<-V.table.beta[,c("KEGG_ID", "NODE.WEIGHT"), with=F]
  
  #Remove duplicates in V(j) if necessary
  dummy<-dummy[!duplicated(dummy$KEGG_ID)]
  
  #Do first node assignment
  colnames(dummy)<-c("KEGG.ID.1", "NODE.WEIGHT.1")
  v.w.table<-merge(x=dummy, y=w.table, by="KEGG.ID.1")
  
  #Do second node assignment
  colnames(dummy)<-c("KEGG.ID.2", "NODE.WEIGHT.2")
  v.w.table<-merge(x=v.w.table, y=dummy, by="KEGG.ID.2")
  
  #Get integrated edge weight with normalized omega parameters
  OMEGA.NODE.1<-max(v.w.table$W.METABOLITE)/max(v.w.table$NODE.WEIGHT.1)
  OMEGA.NODE.2<-max(v.w.table$W.METABOLITE)/max(v.w.table$NODE.WEIGHT.2)
  v.w.table$INTEGRATED.W<-v.w.table$NODE.WEIGHT.1*OMEGA.NODE.1 + v.w.table$NODE.WEIGHT.2*OMEGA.NODE.2 + v.w.table$W.METABOLITE
  
  #Return
  return(v.w.table)
}

for (beta in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  dummy.beta<-network.integration.V.W.V2(get(paste0("Table.V.",beta,".V2",collapse="")),Table.W)
  assign(paste0("Table.W.V.",beta,".V2",collapse="") , dummy.beta)
}

#Save new V2 integrated w.node files - keep in mind INTEGRATED.W here is unnormalized
save(list=paste0(paste0("Table.W.V.", c("0.0",seq(0.1,0.9,0.1),"1.0"),""),".V2",""), file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042814_Tables.W.V.V2.RData")

##Plot W.V(j) distributions (integrated edge weight) as function of beta - V2 VERSION##
Beta.W.V.V2.plot.1<-data.frame(Beta=c(), INTEGRATED.W=c(), WEIGHT.FACTOR=c())
for (beta in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  W.V.V2.TABLE<-get(paste0("Table.W.V.", beta,".V2",collapse=""))
  W.V.V2.TABLE$WEIGHT.FACTOR<-ifelse( (W.V.V2.TABLE$NODE.WEIGHT.1 + W.V.V2.TABLE$NODE.WEIGHT.2)>W.V.V2.TABLE$W.METABOLITE, "V(j)",
                                   ifelse( (W.V.V2.TABLE$NODE.WEIGHT.1 + W.V.V2.TABLE$NODE.WEIGHT.2)<W.V.V2.TABLE$W.METABOLITE,"W(j)",
                                           ifelse( (W.V.V2.TABLE$NODE.WEIGHT.1 + W.V.V2.TABLE$NODE.WEIGHT.2)==W.V.V2.TABLE$W.METABOLITE, "Both",
                                                   "None")))
  
  Beta.W.V.V2.plot.1<-rbind(Beta.W.V.V2.plot.1, data.frame(Beta=beta, INTEGRATED.W=W.V.V2.TABLE$INTEGRATED.W,
                                                     WEIGHT.FACTOR=W.V.V2.TABLE$WEIGHT.FACTOR))
}
head(Beta.W.V.V2.plot.1) #I believe there is no much use of knowing the weight factor(if V is greater than W)
Beta.W.V.V2.plot.1<-as.data.table(Beta.W.V.V2.plot.1)
Beta.W.V.V2.plot.1$SPLIT.BETA<-paste0(Beta.W.V.V2.plot.1$Beta, Beta.W.V.V2.plot.1$WEIGHT.FACTOR)
Beta.W.V.V2.plot.1$Beta<- factor(Beta.W.V.V2.plot.1$Beta )

ggplot(Beta.W.V.V2.plot.1, aes(x=Beta, y=INTEGRATED.W)) + geom_boxplot() + theme.format + 
  scale_colour_brewer(palette="Spectral") #See the effect of V(j) on Total edge weight, since W(j) does not depend on Beta but V(j) does

ggplot(Beta.W.V.V2.plot.1, aes(x=INTEGRATED.W)) + geom_histogram() + theme.format +
  facet_wrap(~Beta) + theme(strip.text.x = element_text(size = 18, colour = "black"))

##Write tables for SPICi cluster analysis for -  V2 VERSION##
#Normalized
for (i in c("0.0",seq(0.1,0.9,0.1),"1.0")) {
  file.name=paste0("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/043014_NETWORK_MEAN_RANKED_",i,"NV2")
  name.table<-get(paste0("Table.W.V.",i,".V2",collapse=""))[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  name.table$INTEGRATED.W<-normalize.vector(name.table$INTEGRATED.W) 
  write.table(file=file.name, name.table, sep="\t", quote=F, row.names=F, col.names=F)
}

##Analyze obtained cluster from V2 VERSION##
#SPICI analysis of normalized clusters
cluster.analysis.NV2_5C<-cluster.analysis.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/", "NV2_5C")
cluster.analysis.NV2_5C<-as.data.table(cluster.analysis.NV2_5C)
cluster.analysis.NV2_5C$BETA<-sapply(cluster.analysis.NV2_5C$FILENAME, function(d) strsplit(as.character(d),"_")[[1]][5])
cluster.analysis.NV2_5C<-cluster.analysis.NV2_5C[order(BETA,CLUSTER.SIZE),]

ggplot(cluster.analysis.NV2_5C, aes(x=BETA, y=CLUSTER.SIZE)) + geom_boxplot(position="identity") + 
  theme.format + geom_jitter(aes(colour=CLUSTER.SIZE))
ggplot(cluster.analysis.NV2_5C, aes(x=BETA,y=CLUSTER.SIZE, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity") + 
  theme.format 
ggplot(cluster.analysis.NV2_5C, aes(x=BETA,y=length(CLUSTER.SIZE)/219, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity", colour="black") + 
  theme.format +
  ylab("Number of Clusters")

##Do significance test on cluster from V2 VERSION##
#Applied function in PC - Version is much faster and uses 10000x permutations
dummy.test<-cluster.analysis.p.value.function.V2("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/043014_NETWORK_MEAN_RANKED_1.0NV2_5C.cluster",
                                              BETA="1.0", Table.W.V.1.0.V2)

####Analysze Enzymes from Top 10 KEGG clusters####
#Get KEGG ids from top 10 clusters
TOP.KEGG<-unlist(sapply(head(common.clusters.freq,10)$Var1, function(ss) strsplit(as.character(ss),"_")[[1]]))

#Get enzymes for these KEGG ids
TOP.KEGG.ENZYMES<-as.vector(BRCA.table.3.processed[KEGG_ID %in% TOP.KEGG,]$Enzyme)

#Get gene expression values for enzymes
TOP.KEGG.EXPRESSION<-as.data.frame(Table.expression)[TOP.KEGG.ENZYMES,]
dim(TOP.KEGG.EXPRESSION)
TOP.KEGG.EXPRESSION<-TOP.KEGG.EXPRESSION[complete.cases(TOP.KEGG.EXPRESSION),]

#Draw heatmap
heatmap.2(as.matrix(TOP.KEGG.EXPRESSION),trace="none",scale="column")

##Do PCA on TOP KEGG expression
TOP.KEGG.PCA<-prcomp(t(TOP.KEGG.EXPRESSION),scale.=T,center=T)
plot(TOP.KEGG.PCA)
head(TOP.KEGG.PCA$x)
barplot(TOP.KEGG.PCA$sdev/TOP.KEGG.PCA$sdev[1])

#Plot PCAs
TOP.KEGG.PCA.SCORES<-as.data.frame(TOP.KEGG.PCA$x)
head(TOP.KEGG.PCA.SCORES)
ggplot(TOP.KEGG.PCA.SCORES, aes(x=PC1, y=PC3 )) + geom_point() +theme.format

#Plot variables (enzymes) as correlations matrices
t(TOP.KEGG.EXPRESSION)[1:3,1:3]
TOP.KEGG.PCA.CORRELATIONS<-as.data.frame(cor(t(TOP.KEGG.EXPRESSION), TOP.KEGG.PCA$x))
TOP.KEGG.PCA.CORRELATIONS[1:4,1:4]
dim(TOP.KEGG.PCA.CORRELATIONS)
corrplot(as.matrix(TOP.KEGG.PCA.CORRELATIONS[140:174,1:6], order="hclust"),tl.cex=0.5,cl.pos="r") 

#CIRCLE PLOT
# function to create a circle
circle <- function(center=c(0,0), npoints=100)
{
  r = 1
  tt = seq(0, 2*pi, length=npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0,0), npoints = 100)

# data frame with arrows coordinates
TOP.KEGG.PCA.FILTERED<-TOP.KEGG.PCA.CORRELATIONS[,c(1,3)]
TOP.KEGG.PCA.FILTERED<-TOP.KEGG.PCA.FILTERED[sqrt(TOP.KEGG.PCA.FILTERED$PC1^2 + TOP.KEGG.PCA.FILTERED$PC3^2)>0.3, ]
head(TOP.KEGG.PCA.FILTERED)

arrows = data.frame(x1=rep(0,nrow(TOP.KEGG.PCA.FILTERED)), y1=rep(0,nrow(TOP.KEGG.PCA.FILTERED)),
                    x2=TOP.KEGG.PCA.FILTERED$PC1, y2=TOP.KEGG.PCA.FILTERED$PC3)

# geom_path will do open circles
ggplot() +
  geom_path(data=corcir, aes(x=x, y=y), colour="gray65") +
  geom_segment(data=arrows, aes(x=x1, y=y1, xend=x2, yend=y2), colour="gray65") +
  geom_text(data=TOP.KEGG.PCA.FILTERED, aes(x=PC1, y=PC3, label=rownames(TOP.KEGG.PCA.FILTERED))) +
  geom_hline(yintercept=0, colour="gray65") +
  geom_vline(xintercept=0, colour="gray65") +
  xlim(-1.1,1.1) + ylim(-1.1,1.1) +
  labs(x="pc1 aixs", y="pc3 axis") #+opts(title="Circle of correlations")

########043014#########
##After meeting with Mona decided that I need to find a way to normalize Table.v##

#Find suitable normalization function
Table.v
ggplot(Table.v, aes(V.METABOLITE)) + geom_histogram(binwidth=0.001)

dummy<-copy(Table.v)

#The inverse hyper sine normalizatoin is suitable for distribution, won't drop zeros
#into -inf values, however we need to choose a suitable theta for it
hyper.sine.normalization<-function(distribution, THETA) {
  hyper.sine<-ifelse(distribution==0,0,
                     log(THETA*distribution + sqrt(distribution^2 * THETA^2  +1))/THETA)
  return(hyper.sine)
}

dummy$hyper.sine<-hyper.sine.normalization(dummy$V.METABOLITE,THETA=10000)
nrow(dummy[normalize.vector(hyper.sine)>0.25,])
ggplot(dummy, aes(normalize.vector(hyper.sine))) + geom_histogram()
ggplot(dummy, aes(hyper.sine)) + geom_histogram()
ggplot(dummy, aes(x=V.METABOLITE, y=normalize.vector(hyper.sine))) + geom_point()

#Apply to v(j) to see effect - if it works we can apply to both before integrating into V(j)
ggplot(Table.u.mean.ranked, aes(x=U.METABOLITE)) + geom_histogram()
ggplot(Table.u.mean.ranked, aes(x=normalize.vector(U.METABOLITE))) + geom_histogram()
ggplot(Table.u.mean.ranked, aes(x=hyper.sine.normalization(U.METABOLITE,1000))) + geom_histogram()

##Apply new normalization scheme into new version of function get V(j) - Version V3##
Table.V.function.V3<-function (u.table, v.table, BETA, THETA, TYPE){
  
  #Merge them by KEGG_ID
  Table.V<-merge(u.table, v.table, by="KEGG_ID")
  
  #Use inverse hypersine normalization on one,both or neither.
  #If you don't want to apply the hyper.sine.normalization function then 
  #you can look at the other Table.V function
  if (TYPE==1) {
    Table.V$U.METABOLITE<-hyper.sine.normalization(Table.V$U.METABOLITE, THETA)
  } else if (TYPE==2) {
    Table.V$V.METABOLITE<-hyper.sine.normalization(Table.V$V.METABOLITE, THETA)
  } else if (TYPE==3) {
    Table.V$U.METABOLITE<-hyper.sine.normalization(Table.V$U.METABOLITE, THETA)
    Table.V$V.METABOLITE<-hyper.sine.normalization(Table.V$V.METABOLITE, THETA)
  } else if (TYPE==0) {
    Table.V$U.METABOLITE<-Table.V$U.METABOLITE
    Table.V$V.METABOLITE<-Table.V$V.METABOLITE
  } else if (TYPE==5){
    Table.V$U.METABOLITE<-normalize.vector(Table.V$U.METABOLITE)
    Table.V$V.METABOLITE<-normalize.vector(Table.V$V.METABOLITE)
  }else
    print ("ERROR!!!, choose correct type: 0,1,2,3")
  
  #APPLY MAX NORMALIZATION
  #If TYPE==5, then this step does not have any effect
  OMEGA<-max(Table.V$V.METABOLITE)/max(Table.V$U.METABOLITE)
  Table.V$U.METABOLITE<-Table.V$U.METABOLITE * OMEGA
  
  #Assign a weight dependency factor based on comparisson of u(j) vs v(j)
  Table.V$WEIGHT.FACTOR<-ifelse(Table.V$U.METABOLITE>Table.V$V.METABOLITE, "u(j)",
                                ifelse(Table.V$U.METABOLITE<Table.V$V.METABOLITE, "v(j)",
                                       "None"))
  
  #Merge by beta - this will be in the form V=beta*u+(1-beta)*v
  Table.V$NODE.WEIGHT<-Table.V$U.METABOLITE*BETA + Table.V$V.METABOLITE*(1-BETA)
  
  #Finally, normalize V(j) between 0-1
  Table.V$NODE.WEIGHT<-normalize.vector(Table.V$NODE.WEIGHT)
  
  #Return values
  return(Table.V)
}

#Apply to obtain new V(j)
for (beta in c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0", seq(1.1,1.5,0.1))) {
  dummy.beta<-Table.V.function.V3(Table.u.mean.ranked, Table.v, BETA=as.numeric(beta), 10000,2)
  assign(paste0("Table.V.",beta,".V3",collapse="") , dummy.beta)
}

#Plot V(j).V3 distributions as function of beta
Beta.V3.plot.1<-data.frame(Beta=c(), NODE.WEIGHT=c(), WEIGHT.FACTOR=c())
for (beta in c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1))) {
  V3.TABLE<-get(paste0("Table.V.", beta,".V3"))
  Beta.V3.plot.1<-rbind(Beta.V3.plot.1, data.frame(Beta=beta, NODE.WEIGHT=V3.TABLE$NODE.WEIGHT,
                                                   WEIGHT.FACTOR=V3.TABLE$WEIGHT.FACTOR))
}
Beta.V3.plot.1<-as.data.table(Beta.V3.plot.1)
ggplot(Beta.V3.plot.1, aes(x=factor(Beta), y=NODE.WEIGHT)) + geom_boxplot() + theme.format + geom_jitter(size=1.5, aes(colour=WEIGHT.FACTOR)) 

#######050214#########

##Integrate V(j).V3 into W(j) for Version 3#
for (beta in c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1))) {
  dummy.beta<-network.integration.V.W.V2(get(paste0("Table.V.",beta,".V3",collapse="")),Table.W)
  assign(paste0("Table.W.V.",beta,".V3",collapse="") , dummy.beta)
}

#Save new V3 integrated w.node files - keep in mind INTEGRATED.W here is unnormalized
save(list=paste0(paste0("Table.W.V.", c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1)),""),".V3",""), 
     file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042814_Tables.W.V.V3.RData")

##Plot W.V(j) distributions (integrated edge weight) as function of beta - V3 VERSION##
Beta.W.V.V3.plot.1<-data.frame(Beta=c(), INTEGRATED.W=c(), WEIGHT.FACTOR=c())
for (beta in c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1))) {
  W.V.V3.TABLE<-get(paste0("Table.W.V.", beta,".V3",collapse=""))
  W.V.V3.TABLE$WEIGHT.FACTOR<-ifelse( (W.V.V3.TABLE$NODE.WEIGHT.1 + W.V.V3.TABLE$NODE.WEIGHT.2)>W.V.V3.TABLE$W.METABOLITE, "V(j)",
                                      ifelse( (W.V.V3.TABLE$NODE.WEIGHT.1 + W.V.V3.TABLE$NODE.WEIGHT.2)<W.V.V3.TABLE$W.METABOLITE,"W(j)",
                                              ifelse( (W.V.V3.TABLE$NODE.WEIGHT.1 + W.V.V3.TABLE$NODE.WEIGHT.2)==W.V.V3.TABLE$W.METABOLITE, "Both",
                                                      "None")))
  
  Beta.W.V.V3.plot.1<-rbind(Beta.W.V.V3.plot.1, data.frame(Beta=beta, INTEGRATED.W=W.V.V3.TABLE$INTEGRATED.W,
                                                           WEIGHT.FACTOR=W.V.V3.TABLE$WEIGHT.FACTOR))
}
head(Beta.W.V.V3.plot.1) #I believe there is no much use of knowing the weight factor(if V is greater than W)
Beta.W.V.V3.plot.1<-as.data.table(Beta.W.V.V3.plot.1)
Beta.W.V.V3.plot.1$SPLIT.BETA<-paste0(Beta.W.V.V3.plot.1$Beta, Beta.W.V.V3.plot.1$WEIGHT.FACTOR)
Beta.W.V.V3.plot.1$Beta<- factor(Beta.W.V.V3.plot.1$Beta )

ggplot(Beta.W.V.V3.plot.1, aes(x=Beta, y=INTEGRATED.W)) + geom_boxplot() + theme.format + 
  scale_colour_brewer(palette="Spectral") 

ggplot(Beta.W.V.V3.plot.1, aes(x=INTEGRATED.W)) + geom_histogram() + theme.format +
  facet_wrap(~Beta) + theme(strip.text.x = element_text(size = 18, colour = "black"))

##Write tables for SPICi cluster analysis for -  V3 VERSION##
#Normalized
for (i in c(seq(-0.5,-0.1,0.1),"0.0",seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1))) {
  file.name=paste0("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/043014_NETWORK_MEAN_RANKED_",i,"NV3")
  name.table<-get(paste0("Table.W.V.",i,".V3",collapse=""))[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
  name.table$INTEGRATED.W<-normalize.vector(name.table$INTEGRATED.W) 
  write.table(file=file.name, name.table, sep="\t", quote=F, row.names=F, col.names=F)
}

##Analyze obtained cluster from V3 VERSION##
#SPICI analysis of normalized clusters
cluster.analysis.NV3_5C<-cluster.analysis.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/", "NV3_5C")
cluster.analysis.NV3_5C<-as.data.table(cluster.analysis.NV3_5C)
cluster.analysis.NV3_5C$BETA<-sapply(cluster.analysis.NV3_5C$FILENAME, function(d) as.numeric(strsplit(strsplit(as.character(d),"_")[[1]][5],"N")[[1]][1]))
cluster.analysis.NV3_5C<-cluster.analysis.NV3_5C[order(BETA,CLUSTER.SIZE),]

head(cluster.analysis.NV3_5C)
ggplot(cluster.analysis.NV3_5C, aes(x=factor(BETA), y=CLUSTER.SIZE)) + geom_boxplot(position="identity") + 
  theme.format + geom_jitter(aes(colour=CLUSTER.SIZE)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) 

ggplot(cluster.analysis.NV3_5C, aes(x=BETA,y=CLUSTER.SIZE, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity") + 
  theme.format + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  guides(fill = guide_legend(keywidth = 3, keyheight = 1))
ggplot(cluster.analysis.NV3_5C[BETA %in% round(seq(0,1,0.1),1),], aes(x=factor(BETA),y=length(CLUSTER.SIZE)/393, fill=factor(CLUSTER.SIZE))) + 
  geom_bar(position="stack",stat="identity", colour="black") + 
  theme.format + ylab("Number of Clusters") +  xlab("Beta")+
  guides(fill = guide_legend(title="Cluster Size", title.theme=element_text(angle=0,size=20),keywidth = 2, keyheight = 1.5)) 


##Do significance test on cluster from V2 VERSION##
#Applied function in PC - Version is much faster and uses 10000x permutations
dummy.test<-cluster.analysis.p.value.function.V2("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/043014_NETWORK_MEAN_RANKED_0.7NV3_5C.cluster",
                                                  BETA="0.7", Table.W.V.0.7.V3)
dummy.test

##First find most common cluster across all Beta SPICi outputs##
common.clusters.V3<-common.cluster.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "NV3_5C")
common.clusters.V3$BETA<-sapply(common.clusters.V3$FILENAME, function(x) strsplit(strsplit(as.character(x) ,"_")[[1]][5], "N")[[1]][1]  )
length(unique(common.clusters.V3$CLUSTER)) #135
common.clusters.freq.V3<-as.data.frame(table(common.clusters.V3$CLUSTER))
common.clusters.freq.V3<-common.clusters.freq.V3[order(common.clusters.freq.V3$Freq,decreasing=T),]
head(common.clusters.freq.V3,10)

#Plot clusters found across betas 
common.clusters.plot.1.function<-function(common.cluster.spici.function.output.beta, beta.range) {
  #Takes output from common.cluster.spici.function(), MAKE SURE IT HAS A BETA AT THE FOURTH COLUMN!!
  #Produces heatmap of clusters vs Betas using cast()
  
  dummy.plot<-common.cluster.spici.function.output.beta[,c(3,4)]
  dummy.plot<-dummy.plot[dummy.plot$BETA %in% beta.range,]
  
  dummy.plot$BETA<-as.numeric(dummy.plot$BETA)
  dummy.plot$PRESENCE<-1
  dummy.plot<-cast(dummy.plot,CLUSTER~BETA, value="PRESENCE",fill=0)
  
  dummy.heatmap<-pheatmap(as.matrix(dummy.plot),scale="none", color=hmcol, cluster_col=F,show_rownames=F,legend=F,
                          treeheight_row=120, fontsize=18,
                          main="Cluster Presence", cellwidth=45,fontsize_row=14, fontsize_col=18)    
  return(dummy.heatmap)
}

common.clusters.plot.1.function(common.clusters.V3, c("0.0",seq(0.1,0.9,0.1),"1.0"))

##Second find most common metabolites across all Beta SPICi outputs##
common.kegg.V3<-common.kegg.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "NV3_5C")
common.kegg.V3$BETA<-sapply(common.kegg.V3$FILENAME, function(x) substr(strsplit(as.character(x) ,"_")[[1]][5], 1,3))
head(common.kegg.V3)
length(unique(common.kegg.V3$KEGG)) #345
common.kegg.freq.V3<-as.data.frame(table(common.kegg.V3$KEGG))
common.kegg.freq.V3<-common.kegg.freq.V3[order(common.kegg.freq.V3$Freq, decreasing=T),]
head(common.kegg.freq.V3,10)
dim(common.kegg.freq.V3)

#Plot KEGGs found across betas
common.kegg.plot.1.function<-function(common.kegg.spici.function.output.beta) {
  #Takes output from common.kegg.spici.functino() MAKE SURE IT HAS BETA ON THIRD COLUMN
  #Produces heatmap of presence of individual metabolites across Betas
  
  dummy.plot<-common.kegg.spici.function.output.beta[,c(2,3)]
  dummy.plot$BETA<-as.numeric(dummy.plot$BETA)
  dummy.plot$PRESENCE<-1
  dummy.plot<-cast(dummy.plot, KEGG~BETA, value="PRESENCE",fill=0)
  
  dummy.heatmap<-pheatmap(as.matrix(dummy.plot),scale="none", color=hmcol, cluster_col=F,show_rownames=F,
                          treeheight_row=120, fontsize=18, legend=F,
                          main="Metabolite Presence", cellwidth=45,fontsize_row=14, fontsize_col=18)
  return(dummy.heatmap)
}

common.kegg.plot.1.function(common.kegg.V3)

####Analysze Enzymes from Top 10 KEGG clusters####
#Get KEGG ids from top 10 clusters
TOP.KEGG.V3<-unlist(sapply(head(common.clusters.freq.V3,10)$Var1, function(ss) strsplit(as.character(ss),"_")[[1]]))

#Get enzymes for these KEGG ids
TOP.KEGG.ENZYMES.V3<-as.vector(BRCA.table.3.processed[KEGG_ID %in% TOP.KEGG.V3,]$Enzyme)

#Get gene expression values for enzymes
TOP.KEGG.EXPRESSION.V3<-as.data.frame(Table.expression)[TOP.KEGG.ENZYMES.V3,]
dim(TOP.KEGG.EXPRESSION.V3)
TOP.KEGG.EXPRESSION.V3<-TOP.KEGG.EXPRESSION.V3[complete.cases(TOP.KEGG.EXPRESSION.V3),]

#Draw heatmap
heatmap.2(as.matrix(TOP.KEGG.EXPRESSION.V3),trace="none",scale="column")

########050314#########

##Draw new function to take in output from common.cluster.spici.function and calculate p-values for all clusters obtained##
#Make sure you have added a "BETA" column to it, similar to common.clusters.V3
#Make sure the total background network has a "BETA" column to be used for each beta cluster group, 

cluster.analysis.p.value.function.V3<-function(ccsf.output, background.network.all.beta) {
  
  #To keep count
  x=1
  
  #Prep vector to store p.values
  dummy.cluster.scores<-data.frame(Beta=c(), CLUSTER.SIZE=c(), P.VALUE=c(), CLUSTER=c())
  
  #Get betas
  dummy.betas<-as.vector(unique(ccsf.output$BETA)) #This will be numeric
  
  #Loop for each beta, have to loop since each background is different
  for (beta in dummy.betas) {
    
    #Get clusters for specified beta
    dummy.cluster.file<-ccsf.output[ccsf.output$BETA==beta,]
    
    #Get background for specified beta and extract desired columns
    dummy.background.beta<-background.network.all.beta[BETA==beta,]
    dummy.background.beta<-dummy.background.beta[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
    
    #Double for dummy cluster, regular for permutation
    dummy.background.beta.double<-rbind(dummy.background.beta, dummy.background.beta[,c(2,1,3),with=F], use.names=F)
    
    #Loop through clusters
    for (cluster in as.vector(dummy.cluster.file$CLUSTER)) {
      
      dummy.vector<-strsplit(cluster,"_")[[1]]
      dummy.vector.size<-length(dummy.vector)
      
      #Create table of edges for cluster to get scores
      cluster.table<-as.data.frame(t(combn(dummy.vector,2))) #No repeats
      colnames(cluster.table)<-colnames(dummy.background.beta.double)[1:2] #To merge by KEGG_IDs
      
      #Merge with background to get scores
      cluster.table<-merge(cluster.table, dummy.background.beta.double, by=colnames(dummy.background.beta.double)[1:2])
      cluster.score<-sum(cluster.table$INTEGRATED.W)      
      
      #Permutation test per score, this depends on cluster size (x10000)
      dummy.vector.edges<-dummy.vector.size*(dummy.vector.size-1)/2  #Number of edges in cluster
      
      permutation.vector<-replicate(10000, sum(sample(dummy.background.beta$INTEGRATED.W, dummy.vector.edges)))
      
      #Get p-value per cluster
      dummy.p.value<-length(permutation.vector[permutation.vector>=cluster.score])/length(permutation.vector)
      
      #Store cluster information
      dummy.cluster.scores<-rbind(dummy.cluster.scores, data.frame(Beta=beta, CLUSTER.SIZE=dummy.vector.size, 
                                                                   P.VALUE=dummy.p.value, CLUSTER=cluster))
    }
    x=x+1
    print (x/length(dummy.betas))
  }
  #Do multiple hypothesis correction across p.values
  dummy.cluster.scores$ADJ.P.VALUE<-p.adjust(dummy.cluster.scores$P.VALUE,method="fdr")
  
  #Return table
  return (dummy.cluster.scores)
}

#Get a combined Table.W.V. .V3 table for all BETAs
Table.W.V.ALL.V3<-data.frame()
for (beta in c(seq(-0.5,-0.1,0.1),"0.0", seq(0.1,0.9,0.1),"1.0",seq(1.1,1.5,0.1))) {
  dummy.table<-get(paste0("Table.W.V.", beta,".V3",collapse=""))
  dummy.table$BETA<-beta
  Table.W.V.ALL.V3<-rbind(Table.W.V.ALL.V3, dummy.table)
}
Table.W.V.ALL.V3

#Run function to calculate p.values for all betas
Table.NETWORK.CLUSTER.P.VALUE<-cluster.analysis.p.value.function.V3(common.clusters.V3, Table.W.V.ALL.V3)
tail(Table.NETWORK.CLUSTER.P.VALUE)
dim(Table.NETWORK.CLUSTER.P.VALUE)
min(Table.NETWORK.CLUSTER.P.VALUE$ADJ.P.VALUE)

##Get function to calculate p-values for KEGG pathways enrichment for each cluster##
#Allow user to have parameter for adjusted p-value cut-off
#This function can take the output straight from cluster.analysis.p.value.function.V3()
#The background pathway from KEGG has the format as in 050314_PATHWAY_TO_COMPOUND
#It will result in a matrix of similar format as input but containing only clusters that have enriched pathway!!!
#Need to include background network compounds

cluster.cpd.kegg.pathways.enrichment<-function(capvfv3, background.pathways.file, background.cpd.pool,cluster.p.cutoff, kegg.p.cutoff) {
    
  #Process background file
  dummy.background<-read.csv(background.pathways.file, sep="\t", header=T, ,colClasses="character")
  colnames(dummy.background)[3]<-"PATHWAY.DESCRIPTION"

  #Apply network backgrond pool to KEGG bacgrkound patwhays
  dummy.background<-dummy.background[dummy.background$COMPOUND %in% background.cpd.pool,]
  dummy.background.compounds<-length(unique(dummy.background$COMPOUND))

  #Process clusters so only those with ADJ.P.VALUE<threshold are kept
  dummy.cluster.matrix<-capvfv3[capvfv3$ADJ.P.VALUE<cluster.p.cutoff,]
  
  #Extract clusters for calculation
  dummy.clusters<-unique(as.vector(dummy.cluster.matrix$CLUSTER))
  
  #To store enriched pathways
  dummy.cluster.enrichment<-data.frame(CLUSTER=c(), PATHWAY=c(), PATHWAY.DESCRIPTION=c(), 
                                       PATHWAY.PVALUE=c(), PATHWAY.COVERAGE=c(), q=c(), m=c(), n=c(), k=c())
  
  #Find enriched pathways through clusters
  for (cluster in dummy.clusters) {
    
    #Unstring cluster
    cluster.vector<-strsplit(cluster,"_")[[1]]
    
    #Reduce background to only look in pathways were there is at least a compound from cluster
    reduced.pathways<-unique(as.vector(dummy.background[dummy.background$COMPOUND %in% cluster.vector,]$PATHWAY))
    reduced.background<-as.data.table(dummy.background[dummy.background$PATHWAY %in% reduced.pathways,])
    
    #Do enrichment analysis
    #Keep in mind that pathway coverage represents how much of the possible pathway is present in the cluster
    pathway.enrichment<-reduced.background[,list(PATHWAY.PVALUE=phyper(q=length(intersect(cluster.vector,COMPOUND))-1,
                                                   m=length(COMPOUND),
                                                   n=dummy.background.compounds-length(COMPOUND),
                                                   k=length(cluster.vector), lower.tail=F) ,
                             PATHWAY.COVERAGE=length(intersect(cluster.vector,COMPOUND))/length(cluster.vector)), 
                       by=c("PATHWAY", "PATHWAY.DESCRIPTION")]
    
    #Store
    dummy.cluster.enrichment<-rbind(dummy.cluster.enrichment,
                                    data.frame(CLUSTER=cluster, pathway.enrichment))  
  }
  
  #Do multiple hypothesis correction
  dummy.cluster.enrichment$PATHWAY.ADJ.PVALUE<-p.adjust(dummy.cluster.enrichment$PATHWAY.PVALUE, method="fdr")
  
  #Apply pathway p.value threshold
  dummy.cluster.enrichment<-dummy.cluster.enrichment[dummy.cluster.enrichment$PATHWAY.ADJ.PVALUE<kegg.p.cutoff,]
  
  #Apply to input to filter for enriched cluster
  dummy.cluster.enrichment<-merge(dummy.cluster.enrichment, dummy.cluster.matrix, by="CLUSTER")
  
  #Return updated cluster table with pathway enrichment analysis
  return(dummy.cluster.enrichment)
}

#Do enrichment analysis using function above
head(Table.NETWORK.CLUSTER.P.VALUE)
Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V3<-cluster.cpd.kegg.pathways.enrichment(Table.NETWORK.CLUSTER.P.VALUE,
                                                                               "DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND",
                                                                               as.vector(Table.v$KEGG_ID),
                                                                               0.05,0.05)

#######050414#######
#Plot enriched KEGG pathways in clusters in terms of BETAs
Table.enrichment.plot.1.function<-function(cluster.cpd.kegg.pathways.enrichment.output, beta.range) {
  #Plots the coverage of KEGG compounds with respect to cluster
  
  dummy.plot<-cluster.cpd.kegg.pathways.enrichment.output[cluster.cpd.kegg.pathways.enrichment.output$Beta %in% beta.range,]
  
  dummy.plot$Beta<-as.numeric(as.vector(dummy.plot$Beta))
  dummy.plot<-cast(dummy.plot[,c(3,5,7)],PATHWAY.DESCRIPTION~Beta,value="PATHWAY.COVERAGE",fill=0,mean)
  dummy.heatmap<-pheatmap(as.matrix(dummy.plot), scale="none", color=hmcol,cluster_cols=F,treeheight_row=120, fontsize=12,
                          main="Cluster Coverage",cellwidth=45,fontsize_row=14, fontsize_col=14)
  return(dummy.heatmap)
}

Table.enrichment.plot.2.function<-function(cluster.cpd.kegg.pathways.enrichment.output, beta.range) {
  #Measure of cluster significance across Betas  
  dummy.plot<-cluster.cpd.kegg.pathways.enrichment.output[,c(3,6,7)]
  dummy.plot<-dummy.plot[dummy.plot$Beta %in% beta.range,]
  dummy.plot$Beta<-as.numeric(as.vector(dummy.plot$Beta))
  dummy.plot<-as.data.table(dummy.plot)
  dummy.plot$LOG.ADJ.PVALUE<--log(dummy.plot$PATHWAY.ADJ.PVALUE)
  
  dummy.plot<-dummy.plot[order(LOG.ADJ.PVALUE, decreasing=T),]
  dummy.plot.casted<-cast(dummy.plot[,c(1,3,4),with=F], PATHWAY.DESCRIPTION~Beta, value="LOG.ADJ.PVALUE", fill=0,mean)
  
  hmcol<-brewer.pal(9, "Blues")
  dummy.heatmap<-pheatmap(as.matrix(dummy.plot.casted),scale="none",color=hmcol,cluster_cols=F,treeheight_row=120, fontsize=12,
                          main="Cluster Significance (-log P-value) vs Beta", cellwidth=45,fontsize_row=14, fontsize_col=14)
  
  return(dummy.heatmap)
}

Table.enrichment.plot.1.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V3, c("0.0",seq(0.1,0.9,0.1),"1.0"))
Table.enrichment.plot.2.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V3, c("0.0",seq(0.1,0.9,0.1),"1.0"))

#######050514#######
##Re-calculate v(j) using differential expression of cancer vs patient##

#Develop function for it
Table.v.function<-function(Table.u.KEGG_IDs, processed.table.3, processed.table.1, cancer.map.file, method="intra", 
                           expression.table.1, expression.table.0=NULL) {
  #Need KEGG_IDs from V(j) counterpart Table.u
  #Need the processed table of enzymes that produced KEGG_IDs (i.e. BRCA.table.3.processed), this is found in 042214_MTables.RData
  #Need the processed table of samples and corresponding mutations (i.e. BRCA.table.1.processed), this is found in 042214_MTables.RData
  #Need the TCGA cancer type map file that comes along with the expression files
  #"method" default will be intra-cancer, change to compare against normal samples
  #expression.table.1 has to be the output of read.expression.file()
  #   - The column names of this table have to be in the form TCGA.XX.XXXX.XXX
  #   - This should only contain cancer samples "01A" or "O1B" as the last 3 digits
  #expression.table.0 has the be the output of read.expression.file()
  #   - This should be supplied if method="inter"
  #   - This should only contain normal samples (i.e. "11A", "11B")
  
  ##Use the map file to process table.1 for samples that have expression data only##
  #Get samples that are cancer (not normal) in map file
  dummy.map<-as.data.table(read.csv(cancer.map.file, header=T, sep="\t"))
  dummy.map$type<-substr(dummy.map$barcode.s., 14,15)
  dummy.map<-dummy.map[type=="01",] #Keep cancer samples only
  dummy.map$Sample<-substr(dummy.map$barcode.s., 1, 16) #Get sample name to reduce duplicates
  
  #Get samples in table.1 that match map file cancer samples
  dummy.table.1<-copy(processed.table.1)
  dummy.table.1$Sample<-substr(dummy.table.1$Tumor_Sample_Barcode,1,16)
  dummy.table.1.v<-merge(dummy.table.1, dummy.map[,c(1,4),with=F], by="Sample") #Table.1 with correspoding expression files
  
  #Prep sample name in expression.table.1 column format
  dummy.table.1.v$Sample<-unlist(lapply(dummy.table.1.v$Sample, function(g) paste(strsplit(g,"-")[[1]],collapse="." ) ) )
  
  #Prep storeage matrix for v(j)
  dummy.table.v<-data.frame(KEGG_ID=c(), V.METABOLITE=c())
  
  #To keep count
  x=0
  
  #Process each metabolite for v(j)
  for (metabolite in Table.u.KEGG_IDs) {
    
    #Get enzymes that produce metabolite
    dummy.metabolite.enzymes<-unique(as.vector(processed.table.3[KEGG_ID==metabolite,]$Enzyme))
    
    #Separate cancer samples in dummy.table.1.v into groups depending on wether they have a mutation or not
    #G0 will not be used if method="inter"
    dummy.groups<-copy(dummy.table.1.v)
    dummy.groups$G<-dummy.table.1.v$Hugo_Symbol %in% dummy.metabolite.enzymes
    dummy.G1.samples<-unique(as.vector(dummy.groups[G==TRUE,]$Sample)) #Samples that have mutation in enzymes
    dummy.G0.samples<-setdiff(unique(as.vector(dummy.groups$Sample)),dummy.G1.samples) #Samples that do not
    
    #Account for fact that there may not be more than 1 sample with mutation, cannot do diff expression in that case
    if (length(dummy.G1.samples)>1) {
      
      #Get expression matrix for each group
      dummy.G1.expression<-expression.table.1[,dummy.G1.samples]
      
      #depending on method G0 matrix will be different
      if (method=="intra") {
        dummy.G0.expression<-expression.table.1[,dummy.G0.samples]
      } else if (method=="inter") {
        dummy.G0.expression<-expression.table.0
      } else {
        dummy.G0.expression<-NULL
        print ("Choose correct method")
      }
      
      ##Do differential expression
      dummy.G.all<-merge(dummy.G1.expression, dummy.G0.expression, by="row.names") #combine into single matrix
      rownames(dummy.G.all)<-dummy.G.all$Row.names
      dummy.G.all$Row.names<-NULL
      dummy.G.all<-as.matrix(dummy.G.all)
      
      #Get design matrix
      dummy.G1.n.samples<-length(colnames(dummy.G1.expression))
      dummy.G0.n.samples<-length(colnames(dummy.G0.expression))
      dummy.design.matrix<-data.frame(G=c(rep("G1", dummy.G1.n.samples), rep("G0", dummy.G0.n.samples)))
      
      #Clean and normalize
      dummy.G.all<-dummy.G.all[complete.cases(dummy.G.all),] #Remove NAs
      dummy.G.all<-normalizeBetweenArrays(dummy.G.all, method="quantile",)
      
      #Fit data
      dummy.fit<-lmFit(dummy.G.all, model.matrix(~G, dummy.design.matrix))
      dummy.eb<-eBayes(dummy.fit)
      
      #Get v(j) by dividing that were differentially expressed over all genes tested
      dummy.all.G.fit<-topTable(dummy.eb, coef=2, n=Inf)
      dummy.G.diff.exp.genes<-nrow(dummy.all.G.fit[dummy.all.G.fit$adj.P.Val<0.05,]) #Value set at p<0.05 for bonferroni corrected p-values
      
      dummy.v.metabolite<-dummy.G.diff.exp.genes/nrow(dummy.all.G.fit)
      
    } else #If 1 or no sample have mutation, then its influence will be equal to zero
      dummy.v.metabolite<-0
    
    #Add v(metabolite) to table
    dummy.table.v<-rbind(dummy.table.v, data.frame(KEGG_ID=metabolite, V.METABOLITE=dummy.v.metabolite))
    
    #To count
    x=x+1
    print (x/length(Table.u.KEGG_IDs))
  }
  #Return v(j) table
  return(dummy.table.v)
}

#Obtain expression data for normal patient
expression.samples.0<-as.data.table(read.csv("DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/FILE_SAMPLE_MAP.txt", header=T, sep="\t"))
expression.samples.0$type<-substr(expression.samples.0$barcode.s., 14,15)
expression.samples.0<-expression.samples.0[expression.samples.0$type=="11",]
expression.samples.0$code<-paste(expression.samples.0$sample, expression.samples.0$filename, sep="$")
head(expression.samples.0)
Table.expression.0<-read.expression.files(files=unique(expression.samples.0$code))
save(Table.expression.0, file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/050514_Table.expression.0.RData")

#Apply function to obtain new weigths for v(j) based on cancer vs normal
Table.v.2<-Table.v.function(unique(as.vector(Table.u.mean.ranked$KEGG_ID)), BRCA.table.3.processed, BRCA.table.1.processed,
                            "DATABASES/CANCER_DATA/TCGA/EXP_GENE/BRCA/040214_LEVEL3/FILE_SAMPLE_MAP.txt", method="inter",
                            Table.expression, expression.table.0=Table.expression.0)
Table.v.2<-as.data.table(Table.v.2)
head(Table.v.2)
dim(Table.v.2[V.METABOLITE!=0,]) #More metabolites with data!
save(Table.v.2, file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/050614_Table.v.2.RData")

#Plot new v(j) distribution
ggplot(Table.v.2, aes(x=V.METABOLITE)) + geom_histogram() + theme.format

##Obtain new V(j) with re-calculated v(j) - Call it V4##
Table.V.assign.version.function<-function(U.TABLE, V.TABLE, beta.range, version,sine.method=0) {
  for (beta in beta.range) {
    dummy.beta<-Table.V.function.V3(U.TABLE, V.TABLE,BETA=as.numeric(beta), THETA=10000, sine.method)
    assign(paste0("Table.V.",beta,".",version,collapse=""),dummy.beta, envir=.GlobalEnv)
  }
}

Table.V.assign.version.function(Table.u.mean.ranked, Table.v.2, c("0.0",seq(0.1,0.9,0.1),"1.0"), "V4",0)

#Plot V(j).V4 as function of beta
V.plot.1.function<-function(version, beta.range) {
  dummy.plot<-data.frame()
  
  #First mix tables
  for (beta in beta.range) {
    dummy.table<-get(paste0("Table.V.",beta,".",version,collapse=""))
    dummy.plot<-rbind(dummy.plot, data.frame(Beta=beta, dummy.table))
  }
  #Plot
  ggplot(dummy.plot, aes(x=factor(Beta), y=NODE.WEIGHT)) + geom_boxplot() + theme.format +
    geom_jitter(size=1.5, aes(colour=WEIGHT.FACTOR))
}

V.plot.1.function("V4",c("0.0",seq(0.1,0.9,0.1),"1.0") )

##Obtain new W(j,j').V4 for Version 4##
Table.W.assign.version.function<-function(W.TABLE, table.V.beta.range, version) {
  for (beta in table.V.beta.range) {
    dummy.beta<-network.integration.V.W.V2(get(paste0("Table.V.", beta,".",version,collapse="")), W.TABLE)
    assign(paste0("Table.W.V.",beta,".",version,collapse=""),dummy.beta, envir=.GlobalEnv)
  }
}
Table.W.assign.version.function(Table.W,c("0.0",seq(0.1,0.9,0.1),"1.0"), "V4" )

#Save new V4 integrated w.node files - keep in mind INTEGRATED.W here is unnormalized
save(list=paste0(paste0("Table.W.V.", c("0.0",seq(0.1,0.9,0.1),"1.0"),""),".V4",""), file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042814_Tables.W.V.V4.RData")

#Plot W.V4 distributions as functions of beta 
W.plot.1.function<-function(table.W.beta.range, version) {
  dummy.plot<-data.frame()
  
  for (beta in table.W.beta.range) {
    dummy.table<-get(paste0("Table.W.V.",beta,".",version,collapse=""))
    dummy.plot<-rbind(dummy.plot, data.frame(Beta=beta, dummy.table))
  }
  
  ggplot(dummy.plot, aes(x=factor(Beta), y=INTEGRATED.W)) + geom_boxplot() + theme.format
}
W.plot.2.function<-function(table.W.beta.range, version) {
  dummy.plot<-data.frame()
  
  for (beta in table.W.beta.range) {
    dummy.table<-get(paste0("Table.W.V.",beta,".",version,collapse=""))
    dummy.plot<-rbind(dummy.plot, data.frame(Beta=beta, dummy.table))
  }
  
  ggplot(dummy.plot, aes(x=INTEGRATED.W)) + geom_histogram() + theme.format +
    facet_wrap(~Beta) + theme(strip.text.x=element_text(size=18, colour="black"))
}

W.plot.1.function(c("0.0",seq(0.1,0.9,0.1),"1.0"), "V4")
W.plot.2.function(c("0.0",seq(0.1,0.9,0.1),"1.0"), "V4")

##Write tables (normalized) for SPICi cluster analysis for - V4 VERSION##
write.for.spici.normalized.function<-function(integrated.beta.range, folder, name, version) {
  for (beta in integrated.beta.range) {
    dummy.file.name<-paste0(folder,"/",name,beta,"N",version,collapse="")
    dummy.table<-get(paste0("Table.W.V.",beta,".",version,collapse=""))[,c("KEGG.ID.2", "KEGG.ID.1","INTEGRATED.W"),with=F]
    dummy.table$INTEGRATED.W<-normalize.vector(dummy.table$INTEGRATED.W) 
    write.table(file=dummy.file.name, dummy.table, sep="\t", quote=F, row.names=F, col.names=F)
  }
}
write.for.spici.normalized.function(c("0.0",seq(0.1,0.9,0.1),"1.0"),
                                    "DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS","050614_NETWORK_MEAN_RANKED_","V4")


##Analyze obtained cluster from V4 VERSION##
#SPICI analysis of normalized clusters
cluster.analysis.NV4_5C<-cluster.analysis.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS/", "NV4_5C")

#Plot clusters
CLUSTER.plot.1.function<-function(cluster.analysis.output) {
  ggplot(cluster.analysis.output, aes(x=BETA, y=CLUSTER.SIZE)) + geom_boxplot(position="identity") + 
    theme.format + geom_jitter(aes(colour=CLUSTER.SIZE)) + 
    stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) 
}
CLUSTER.plot.1.function(cluster.analysis.NV4_5C) + scale_y_log10()

CLUSTER.plot.2.function<-function(cluster.analysis.output) {
  ggplot(cluster.analysis.output, aes(x=BETA,y=CLUSTER.SIZE, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity") + 
    theme.format   
}
CLUSTER.plot.2.function(cluster.analysis.NV4_5C) 

ggplot(cluster.analysis.NV4_5C, aes(x=BETA,y=length(CLUSTER.SIZE)/160, fill=factor(CLUSTER.SIZE))) + geom_bar(position="stack",stat="identity", colour="black") + 
  theme.format + ylab("Number of Clusters")

##Plot common clusters across Betas##
common.clusters.V4<-common.cluster.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "NV4_5C")
common.clusters.V4$BETA<-sapply(common.clusters.V4$FILENAME, function(x) substr(strsplit(as.character(x) ,"_")[[1]][5], 1,3))
common.clusters.plot.1.function(common.clusters.V4)

##Plot common KEGG cpds across Betas##
common.kegg.V4<-common.kegg.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "NV4_5C")
common.kegg.V4$BETA<-sapply(common.kegg.V4$FILENAME, function(x) substr(strsplit(as.character(x) ,"_")[[1]][5], 1,3))
common.kegg.plot.1.function(common.kegg.V4)

##Do significance test for clusters in Version 4##
#Get a combined Table.W.V. .V4 table for all BETAs
Background.Table.W.V.function<-function(beta.range, version) {
  dummy.background<-data.frame()
  
  for(beta in beta.range) {
    dummy.table<-get(paste0("Table.W.V.",beta,".",version, collapse=""))
    dummy.table$BETA<-beta
    dummy.background<-rbind(dummy.background, dummy.table)
  }
  return(dummy.background)
}
Table.W.V.ALL.V4<-Background.Table.W.V.function(c("0.0", seq(0.1,0.9,0.1),"1.0"), "V4")

Table.NETWORK.CLUSTER.P.VALUE.V4<-cluster.analysis.p.value.function.V3(common.clusters.V4, Table.W.V.ALL.V4)

##Do KEGG pathway enrichemtn for Version 4##
Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V4<-cluster.cpd.kegg.pathways.enrichment(Table.NETWORK.CLUSTER.P.VALUE.V4,
                                                                               "DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND",
                                                                               as.vector(Table.v.2$KEGG_ID),
                                                                               0.05,0.05)

#Plot enriched pathways
Table.enrichment.plot.1.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V4)
Table.enrichment.plot.2.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.V4)

######050614#######
##Try different ways#
Table.V.assign.version.function(Table.u.mean.ranked, Table.v.2, c("0.0",seq(0.1,0.9,0.1),"1.0"), "V5",5)

V.plot.1.function("V5",c("0.0",seq(0.1,0.9,0.1),"1.0"))

Table.W.assign.version.function(Table.W,c("0.0",seq(0.1,0.9,0.1),"1.0"), "V5" )
save(list=paste0(paste0("Table.W.V.", c("0.0",seq(0.1,0.9,0.1),"1.0"),""),".V5",""), file="DATABASES/CANCER_DATA/METABOLOMICS/TABLES/042814_Tables.W.V.V5.RData")

W.plot.1.function(c("0.0",seq(0.1,0.9,0.1),"1.0"), "V5")
W.plot.2.function(c("0.0",seq(0.1,0.9,0.1),"1.0"), "V5")

write.for.spici.normalized.function(c("0.0",seq(0.1,0.9,0.1),"1.0"),
                                    "DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS","050614_NETWORK_MEAN_RANKED_","V5")

######052114########
##Design new edge weight calculation W(j,j') based on conditional differential gene expression##
#This calculation will result in a zero edge weight if none of the genes producing metabolites are differentially expressed
#This serves a different purpose than before, since instead of looking for concordance of levels, we are looking for only dysregulated metabolite concordance

Table.W.function<-function(normal.expression.samples,cancer.expression.samples, table.3.processed, network.KEGG.IDs,
                           differential.threshold=0.1){
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
  
  require(Hmisc)
  
  #STEP 0 - Find differentially expressed genes first between both groups
  dummy.exp.all<-merge(cancer.expression.samples, normal.expression.samples, by="row.names") #combine into single matrix
  rownames(dummy.exp.all)<-dummy.exp.all$Row.names
  dummy.exp.all$Row.names<-NULL
  dummy.exp.all<-as.matrix(dummy.exp.all)
  
  #Get design matrix
  dummy.cancer.n.samples<-length(colnames(cancer.expression.samples))
  dummy.normal.n.samples<-length(colnames(normal.expression.samples))
  dummy.design.matrix<-data.frame(EXP=c(rep("cancer", dummy.cancer.n.samples), rep("normal", dummy.normal.n.samples)))
  
  #Clean and normalize
  dummy.exp.all<-dummy.exp.all[complete.cases(dummy.exp.all),] #Remove NAs
  dummy.exp.all<-normalizeBetweenArrays(dummy.exp.all, method="quantile")
  
  #Fit data
  dummy.fit<-lmFit(dummy.exp.all, model.matrix(~EXP, dummy.design.matrix))
  dummy.eb<-eBayes(dummy.fit)  
  
  #Get set of differentially expressed genes between cancer and normal - CHECK LIMMA "ID" LIMITATION
  dummy.diff.exp.genes<-topTable(dummy.eb, coef=2, n=Inf)
  
  if("ID" %in% colnames(dummy.diff.exp.genes)) {
    dummy.diff.exp.genes<-as.vector(dummy.diff.exp.genes[dummy.diff.exp.genes$adj.P.Val<0.05,]$ID)  
  }else
    dummy.diff.exp.genes<-unique(rownames(dummy.diff.exp.genes[dummy.diff.exp.genes$adj.P.Val<0.05,]))
  
  print (length(dummy.diff.exp.genes)) #CHANGE!!
  
  #STEP 1 - Use network.kegg.ids to filter table.3 (enzyme table) - THIS WILL BE THE TOTAL NUMBER OF NODES IN NETWORK
  dummy.table.3<-as.data.table(table.3.processed)
  dummy.table.3$Enzyme<-as.character(dummy.table.3$Enzyme)
  dummy.table.3<-dummy.table.3[KEGG_ID %in% network.KEGG.IDs,]
  print ("STEP 1 - DONE")
  
  #STEP 2 - Construct background network (ZERO NETWORK)
  dummy.zero.network<-as.data.table(t(combn(unique(as.vector(dummy.table.3$KEGG_ID)),2)))
  setnames(dummy.zero.network, c("V1","V2"),c("KEGG.ID.1","KEGG.ID.2"))
  dummy.zero.network$W.METABOLITE<-0
  dummy.zero.network$P.VALUE<-NA #BEWARE!! - FOR CORRELATIONS WHEN BOTH ARE NOT DYSREGULATED
  print ("STEP 2 - DONE")
  print (dim(dummy.zero.network))
  
  #STEP 3 - Use genes that are differentially expressed to drop metabolites that have less differentially expressed genes than threshold
  dummy.sign.metabolites.table<-dummy.table.3[,list(DIFF.EXP=(length(intersect(Enzyme, dummy.diff.exp.genes))>=
                                                      length(Enzyme)*differential.threshold)), 
                                              by="KEGG_ID"]
  
  dummy.sign.metabolites<-unique(as.vector(dummy.sign.metabolites.table[DIFF.EXP==TRUE,]$KEGG_ID)) #Metabolites that pass threshold
  dummy.table.3<-dummy.table.3[KEGG_ID %in% dummy.sign.metabolites,] #Only metabolites that have above threshold diff exp genes and all of their enzymes,
                                                                      #including those that are not differentially expressed (since rest passed threshold)
  print ("STEP 3 - DONE")
  
  #STEP 4 - Reduce cancer.expression matrix to contain only enzyme values for metabolites that pass threshold (dummy.table.3)
  dummy.cancer.sign.expression<-as.data.frame(cancer.expression.samples,keep.rownames=T)
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
  
  dummy.pre.W.table.RHO<-rcorr(as.matrix(t(dummy.cancer.sign.expression.mean)), type="spearman")$r #To calcculate RHO
  dummy.pre.W.table.PVALUE<-rcorr(as.matrix(t(dummy.cancer.sign.expression.mean)), type="spearman")$P #To calculate p-value
  
  dummy.pre.W.table.RHO<-melt(dummy.pre.W.table.RHO) #To convert to pairwise table
  dummy.pre.W.table.PVALUE<-melt(dummy.pre.W.table.PVALUE) #To convert to pairwise table
  setnames(dummy.pre.W.table.RHO, c("X1", "X2", "value"), c("KEGG.ID.1", "KEGG.ID.2", "W.METABOLITE"))
  setnames(dummy.pre.W.table.PVALUE, c("X1", "X2", "value"), c("KEGG.ID.1", "KEGG.ID.2", "P.VALUE"))
  
  dummy.pre.W.table<-merge(dummy.pre.W.table.RHO, dummy.pre.W.table.PVALUE) #Merge to obtain p.vaues and rhos
  dummy.pre.W.table<-dummy.pre.W.table[dummy.pre.W.table$KEGG.ID.1!=dummy.pre.W.table$KEGG.ID.2,] #To delete diagonal correlation (node to itself)
  #dummy.pre.W.table<-cor(t(dummy.cancer.sign.expression.mean), method="spearman")
  #dummy.pre.W.table<-melt(dummy.pre.W.table) #To convert it to a pairwise table
  #dummy.pre.W.table<-dummy.pre.W.table[dummy.pre.W.table$X1!=dummy.pre.W.table$X2,] #To delete diagonal correlation (node to itself)
  #colnames(dummy.pre.W.table)<-c("KEGG.ID.1","KEGG.ID.2","W.METABOLITE")
  dummy.pre.W.table$W.METABOLITE<-abs(dummy.pre.W.table$W.METABOLITE) #ABSOLUTE VALUE OF SPEARMAN
  print ("STEP 6 - DONE")
  
  #STEP 7 - Integrate with background table to obtain final W.table (to have info of 0 edge weight)
  dummy.W.table<-as.data.frame(rbind(dummy.pre.W.table, dummy.zero.network)) #Mix background and sign metabolic tables
  dummy.W.table<-dummy.W.table[!duplicated(dummy.W.table[,1:2]),]#Removes zeroes from background for diff metabolites
                                                                  #This works because pairwise table above has two-sided combinations
  dummy.W.table.combn<-as.data.frame(dummy.zero.network)[,1:2]#To remove two sided duplicates - this has all the edge permutations we need
    
  dummy.W.table<-as.data.table(merge(dummy.W.table, dummy.W.table.combn))
  print ("STEP 7 - DONE")
  
  #STEP 8 - Clean table - ONLY USE WHEN NOT DOING P.VALUES!!!!
  #dummy.W.table[is.na(dummy.W.table)]<-0
  
  #STEP 9 - Do multiple hypothesis testing
  dummy.W.table$P.ADJUSTED<-p.adjust(dummy.W.table$P.VALUE, method="fdr")
  
  #END
  return (dummy.W.table)
}

Table.W.assign.function<-function(normal.expression.samples, cancer.expression.samples, table.3.processed, network.KEGG.IDs,
                                  threshold.range) {
  #Function calculates W(j,j') using Table.W.function() for different given thresholds betwen 0-1
  for (i in threshold.range) {
    dummy.W.table<-Table.W.function(normal.expression.samples, cancer.expression.samples, table.3.processed, network.KEGG.IDs, i)
    assign(paste0("Table.W.T.",i), dummy.W.table, envir=.GlobalEnv)
  }
}

Table.W.assign.function(Table.expression.0, Table.expression, BRCA.table.3.processed, unique(as.vector(Table.u.mean.ranked$KEGG_ID)),
                        threshold.range=c("0.0",seq(0.1,1,0.1)))

Table.W.figure.1<-function(table.name, known.threshold.range) {
  #This function draws threshold.range vs number of significant edges
  #This information is drawn from the Table.W.assin.function
  
  #To store
  dummy.all.table<-data.frame(Differential.Threshold=c(), Significant.W=c())
  
  #Loop through W threshold tables to get significant edges
  for (i in known.threshold.range) {
    dummy.table<-get(paste0(table.name,i))
    dummy.significant.W<-nrow(dummy.table[P.ADJUSTED<0.05,])
    dummy.all.table<-rbind(dummy.all.table, 
                           data.frame(Differential.Threshold=i, Significant.W=dummy.significant.W))
  }
  
  #Plot threshold vs number of significant W
  PLOT.1<-ggplot(dummy.all.table, aes(x=Differential.Threshold, y=Significant.W)) + geom_point(aes(size=2)) +
    theme.format + geom_line(aes(group=1))
  
  return(PLOT.1)
}

Table.W.Best.T<-function (table.name, known.threshold.range) {
  #This function calculats produces differential scores for all thresholds (Higher score means better choice of threshold)
  #This information is drawn from the Table.W.assin.function
  
  #To store
  dummy.all.table<-data.frame(Differential.Threshold=c(), Significant.W=c())
  
  #Loop through W threshold tables to get significant edges
  for (i in known.threshold.range) {
    dummy.table<-get(paste0(table.name,i))
    dummy.significant.W<-nrow(dummy.table[P.ADJUSTED<0.05,])
    dummy.all.table<-rbind(dummy.all.table, 
                           data.frame(Differential.Threshold=i, Significant.W=dummy.significant.W))
  }
  
  dummy.all.table$Differential.Threshold<-as.numeric(as.character(dummy.all.table$Differential.Threshold))
  
  #MAKE SURE Thresholds are ordered from smallest to largest!!
  dummy.all.table<-dummy.all.table[order(dummy.all.table$Differential.Threshold),]
  
  ##Calculation
  #This calculation weights the highest threhold levels against highest percent Significant edges kept 
  dummy.all.table$Threshold.Score<-sapply(rownames(dummy.all.table), 
                            function (x) prod(dummy.all.table[x,])*10* 
                              (dummy.all.table[x,"Significant.W"]/dummy.all.table[as.numeric(x)-1,"Significant.W"])
                            )

  #Replace NA score given to initial value with zero
  dummy.all.table$Threshold.Score<-as.numeric(dummy.all.table$Threshold.Score)
  dummy.all.table[is.na(dummy.all.table)]<-0
  
  #Normalize scores between 0-1 and order them in table
  dummy.all.table$Threshold.Score<-normalize.vector(dummy.all.table$Threshold.Score)
  dummy.all.table<-dummy.all.table[order(dummy.all.table$Threshold.Score, decreasing=T),]
  
  #Return table and scores
  return(dummy.all.table)
}

Table.W.figure.1("Table.W.T.", c("0.0",seq(0.1,1.0,0.1))) + 
  theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18))
Table.W.Best.T("Table.W.T.", c("0.0", seq(0.1,1.0,0.1))) #Found best Differential Threshold score to be 0.5 (Gamma)

#######053114#######
##Design function to compare v(j) and u(j) across all Differential Scores to justify Best Differential Threshold Score (Gamma) found
#   - Ideally we will find that at this Gamma we have minimal loss of high v(j) and u(j) while removing most low u(j) and low v(j)

##Write function
Tables.Gamma.Comparisson.function<-function(u.table, v.table, w.table.name, w.range) {
  #These are:
  #   - u.table in the format Table.u.mean.ranked
  #   - v.table in the format Table.v.2
  #   - w.table.name and w.range depend on what was used in Table.W.assign.function
  
  #Merge u.table and v.table to have node weight in single table
  dummy.u.v.table<-merge(u.table, v.table, by="KEGG_ID")
  
  #Table to store
  dummy.classified.table<-data.frame(METABOLITE=c(), KEGG_ID=c(), U.METABOLITE=c(),
                                     V.METABOLITE=c(), GAMMA.SIGNIFICANCE=c(), GAMMA=c())
  
  #Loop through w.table to obtain significant metabolites per Gamma
  for (gamma in w.range) {
    
    #Draw Table.W per Gamma
    dummy.w.table<-get(paste0("Table.W.T.", gamma))
    
    #Get significant metabolites
    dummy.significant.metabolites<-unique(union(as.vector(dummy.w.table[P.ADJUSTED<0.05,]$KEGG.ID.1) , 
                                                as.vector(dummy.w.table[P.ADJUSTED<0.05,]$KEGG.ID.2)))
    
    #Classify node weights in merge table based on Gamma sifnificance
    dummy.gamma.table<-copy(dummy.u.v.table)
    dummy.gamma.table$GAMMA.SIGNIFICANCE<-dummy.gamma.table$KEGG_ID %in% dummy.significant.metabolites
    
    #Assign Gamma column
    dummy.gamma.table$GAMMA<-gamma
    
    #Store
    dummy.classified.table<-rbind(dummy.classified.table, dummy.gamma.table)
  }
  
  dummy.classified.table<-as.data.table(dummy.classified.table)
  
  #Return
  return(dummy.classified.table)
}

#We will test the function on Table.v.2 first, although this may not necessarily be accurate due to the inverse hyper sine normalization used
A<-Tables.Gamma.Comparisson.function(Table.u.mean.ranked, Table.v.2, "Table.W.T.", c("0.0", seq(0.1, 1.0, 0.1)))
ggplot(A, aes(x=GAMMA, y=V.METABOLITE)) + geom_boxplot() + theme.format +
  geom_jitter(size=2.0, aes(colour=GAMMA.SIGNIFICANCE)) 
ggplot(A, aes(x=GAMMA, y=V.METABOLITE,fill=factor(GAMMA.SIGNIFICANCE))) + geom_boxplot() + theme.format +
  theme(legend.position="bottom")

ggplot(A, aes(x=GAMMA, y=U.METABOLITE,fill=factor(GAMMA.SIGNIFICANCE))) + geom_violin() + theme.format 
ggplot(A, aes(x=GAMMA, y=U.METABOLITE,fill=factor(GAMMA.SIGNIFICANCE))) + geom_boxplot() + theme.format

ggplot(A, aes(x=V.METABOLITE, fill=factor(GAMMA.SIGNIFICANCE))) + geom_histogram() + theme.format +
  facet_wrap(~GAMMA) + theme(strip.text.x = element_text(size = 30))

########060414#########
##Use new gamma (0.5) to filter out metabolites in u(j) and v(j)
#   - Will use Table.v.2, since this is what had differences, and it is more reliable

Function.Omega.Integration<-function(W.OMEGA.TABLE, TABLE.U, TABLE.V, beta.range, LABEL,
                                     DATA.FOLDER, SPICI.FOLDER, FILE.DATE) {
  #This function will build the whole network given an omega.table input to filter out individual node weights
  #It will produce a table per edge indicating the edge weight
  #LABEL such as "V4", "OMEGA.1"
  
  #Filter omega table for significant edges only, REMOVES ALL OTHER EDGES FROM NETWORK!!
  dummy.w.omega.table<-W.OMEGA.TABLE[P.ADJUSTED<0.05,]
  
  #Obtain omega significant metabolites
  dummy.omega.metabolites<-unique(union(as.vector(W.OMEGA.TABLE[P.ADJUSTED<0.05,]$KEGG.ID.1), 
                                        as.vector(W.OMEGA.TABLE[P.ADJUSTED<0.05,]$KEGG.ID.2)))
    
  #Filter out u(j) and v(j) by significant omega metabolites
  dummy.table.u<-TABLE.U[KEGG_ID %in% dummy.omega.metabolites, ]
  dummy.table.v<-TABLE.V[KEGG_ID %in% dummy.omega.metabolites, ]
  
  #Integrate using previous function
  #Type "0" doesn't use hyper sine normalization and instead uses a max ratio
  Table.V.assign.version.function(dummy.table.u, dummy.table.v, beta.range, LABEL,0)
  
  #Integrate to edge weight using previous function
  #This works since the limiting dataset will be whatever we feed it for w.omega.table, that is, originally zero nodes in these network will
  #not count as they did in previous calculations when V(j) and V(j') were added regardless if W(j,j') was zero.
  Table.W.assign.version.function(dummy.w.omega.table, beta.range, LABEL)
 
  #Save Integrated W tables as data files
  save(list=paste0(paste0("Table.W.V.", beta.range,""),".",LABEL,""), file=paste0(DATA.FOLDER,"/",FILE.DATE,"_", "Tables.W.V.",LABEL,".Rdata",""))
  
  #Save table for SPICI analysis
  write.for.spici.normalized.function(beta.range, SPICI.FOLDER, paste0(FILE.DATE,"_NETWORK_MEAN_RANKED_",""),LABEL)
}

#Run function with W.OMEGA=0.5
Function.Omega.Integration(Table.W.T.0.5, Table.u.mean.ranked, Table.v.2, c("0.0",seq(0.1,1,0.1)), "OMEGA.T.0.5",
                           "DATABASES/CANCER_DATA/METABOLOMICS/TABLES",
                           "DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "060414")

#Plot
V.plot.1.function("OMEGA.T.0.5",c("0.0",seq(0.1,1,0.1)) )
W.plot.1.function(c("0.0",seq(0.1,1.0,0.1)), "OMEGA.T.0.5")
W.plot.2.function(c("0.0",seq(0.1,1.0,0.1)), "OMEGA.T.0.5")

#Analyze SPICi clusters
##Plot common clusters across Betas##
common.clusters.OMEGA.T.0.5<-common.cluster.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "T.0.5_5C")
common.clusters.OMEGA.T.0.5$BETA<-sapply(common.clusters.OMEGA.T.0.5$FILENAME, function(x) strsplit(strsplit(as.character(x) ,"_")[[1]][5],"N")[[1]][1])
common.clusters.plot.1.function(common.clusters.OMEGA.T.0.5, c("0.0", seq(0.1,1,0.1)))

##Plot common KEGG cpds across Betas## -FIX!!
common.kegg.OMEGA.T.0.5<-common.kegg.spici.function("DATABASES/CANCER_DATA/METABOLOMICS/TABLES/CLUSTER_ANALYSIS", "T.0.5_5C")
common.kegg.OMEGA.T.0.5$BETA<-sapply(common.kegg.OMEGA.T.0.5$FILENAME, function(x) strsplit(strsplit(as.character(x) ,"_")[[1]][5],"N")[[1]][1])
common.kegg.plot.1.function(common.kegg.OMEGA.T.0.5)

##Do significance test for clusters in OMEGA.T.0.5##
Table.W.V.ALL.OMEGA.T.0.5<-Background.Table.W.V.function(c("0.0", seq(0.1,1,0.1)), "OMEGA.T.0.5")

Table.NETWORK.CLUSTER.P.VALUE.OMEGA.T.0.5<-cluster.analysis.p.value.function.V3(common.clusters.OMEGA.T.0.5, Table.W.V.ALL.OMEGA.T.0.5)

##Do KEGG pathway enrichemtn for Version 4##
Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.OMEGA.T.0.5<-cluster.cpd.kegg.pathways.enrichment(Table.NETWORK.CLUSTER.P.VALUE.OMEGA.T.0.5,
                                                                                  "DATABASES/KEGG/050314_PATHWAY_TO_COMPOUND",
                                                                                  as.vector(Table.v.2$KEGG_ID),
                                                                                  0.05,0.05)

#Plot enriched pathways
Table.enrichment.plot.1.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.OMEGA.T.0.5, c("0.0",seq(0.1,1,0.1)))
Table.enrichment.plot.2.function(Table.NETWORK.CLUSTER.P.VALUE.ENRICHMENT.OMEGA.T.0.5,c("0.0",seq(0.1,1,0.1)))
