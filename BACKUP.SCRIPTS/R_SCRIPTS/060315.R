#Metabolic Hub enrichment

library(data.table)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(made4)
library(pheatmap)
library(gplots)
library(network)
library(igraph)

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

#Load recon datasets and keep only metabolites that have either kegg of hmdb ids
recon.table<-Function.process.recon.simple("DATABASES/RECON/042215.PROCESSED.GENE", "DATABASES/RECON/042215.PROCESSED.METABOLITES",
                                           "DATABASES/RECON/042215.PROCESSED.REACTIONS", filter=60)
recon.table<-do.call(rbind, recon.table)
recon.table<-recon.table[CHEBI!="NONE",]
recon.table<-recon.table[,c("Hugo_Symbol", "NAME", "KEGG_ID", "HMDB"), with =F]
setkey(recon.table)
recon.table<-unique(recon.table)
recon.table<-recon.table[!grepl("_", Hugo_Symbol),]
recon.table<-recon.table[!(KEGG_ID=="NONE" & HMDB=="NONE"),]

saveRDS(recon.table, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/060415.recon.table.simple.rds")

#Load expression data
brca.ma<-readRDS("PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/041815.BRCA.MATRICES.AGILENT.rds")

#Calculate metabolic hub score #FIX FUNCTION TO AT LEAST 3 INSTEAD OF 4!!!!!
Function.MH.Exp.Score<-function(recon.table, brca.exp, mh.hmdb.kegg, HMDB.table, alpha=0.3){
  #This will use pre-filtered recon.table that contains either a viable KEGG or HMDB identifier
  
  #Double check and clean up expression datasets
  brca.exp$normal<-brca.exp$normal[complete.cases(brca.exp$normal),]
  brca.exp$tumor<-brca.exp$tumor[complete.cases(brca.exp$tumor),]
  
  #Get unique identifiers for recon table (Use primarily kegg id, only use hmdb id if kegg id not present)
  kegg.ids<-recon.table[KEGG_ID!="NONE",]
  hmdb.ids<-recon.table[KEGG_ID=="NONE" & HMDB!="NONE",]
  
  kegg.ids<-kegg.ids[,c("Hugo_Symbol", "KEGG_ID"),with=F]
  setnames(kegg.ids, c("Hugo_Symbol", "ID"))
  hmdb.ids<-hmdb.ids[,c("Hugo_Symbol", "HMDB"), with=F]
  setnames(hmdb.ids, c("Hugo_Symbol", "ID"))
  
  recon.table<-unique(rbind(kegg.ids, hmdb.ids))
  print ("Processed recon table")
  
  #Regularize some HMDB IDs to KEGG IDs
  id.change<-fread(mh.hmdb.kegg, header=T)
  id.change<-id.change[KEGG!="Not Found",]
  z<-id.change$KEGG
  names(z)<-id.change$HMDB
  
  recon.table$ID<-ifelse(recon.table$ID %in% names(z), z[recon.table$ID], recon.table$ID)
  
  #Add HMDB processed KEGG targeted table
  hmdb.source<-fread(HMDB.table, header=T, sep="\t", stringsAsFactors = F)
  hmdb.source<-hmdb.source[,c("KEGG_ID", "Hugo_Symbol"),with=F]
  setnames(hmdb.source, c("ID", "Hugo_Symbol"))
  recon.table<-unique(rbind(recon.table, hmdb.source[,c("Hugo_Symbol","ID"),with=F]))
  
  #Filter recon.table and brca.exp for common genes
  common.genes<-intersect(unique(recon.table$Hugo_Symbol), intersect(rownames(brca.exp$tumor),rownames(brca.exp$normal)) )
  recon.table<-recon.table[Hugo_Symbol %in% common.genes,]
  brca.exp$tumor<-brca.exp$tumor[common.genes,]
  brca.exp$normal<-brca.exp$normal[common.genes,]
  
  #Filter recon.table for mets that have at least 3 genes assigned to them
  print (length(unique(recon.table$ID)))
  recon.count<-recon.table[,list(N.GENES=length(Hugo_Symbol)), by ="ID"]
  recon.count<-recon.count[N.GENES>=3, ]
  recon.table<-recon.table[ID %in% recon.count$ID,]
  print (length(unique(recon.table$ID)))
  
  #Calculate correlation matrix in cancer and normal
  cor.normal<-rcorr(t(data.matrix(brca.exp$normal)), type = "spearman")
  cor.tumor<-rcorr(t(data.matrix(brca.exp$tumor)), type = "spearman")
  print ("Done with correlation tables")
  
  #Correct correlation p-values for multiple hypothesis testing and obtain their signficance status (P-val<0.05)
  cor.normal$P<-matrix(p.adjust(cor.normal$P, method = "fdr"), ncol=ncol(cor.normal$P), 
                              dimnames = list(rownames(cor.normal$P), colnames(cor.normal$P)))
  cor.tumor$P<-matrix(p.adjust(cor.tumor$P, method = "fdr"), ncol=ncol(cor.tumor$P), 
                              dimnames = list(rownames(cor.tumor$P), colnames(cor.tumor$P)))
  cor.normal$P<-cor.normal$P<0.05
  cor.tumor$P<-cor.tumor$P<0.05
  diag(cor.normal$P)<-FALSE 
  diag(cor.tumor$P)<-FALSE
  
  #Use corrrelation value only if it passes multiple hypothesis testing (Self correlation value will be zero)
  cor.normal$r<-cor.normal$r * cor.normal$P
  cor.tumor$r<-cor.tumor$r * cor.tumor$P
  
  #Obtain gene counts
  gene.counts<-recon.table[,list(N.GENES=length(unique(Hugo_Symbol))), by="ID"]
  
  #Get score per metabolite (pval of paired changes in edges)
  internal.function<-function(Hugo_Symbol){
    norm.vec<-cor.normal$r[Hugo_Symbol, Hugo_Symbol][upper.tri(cor.normal$r[Hugo_Symbol, Hugo_Symbol])]
    cancer.vec<-cor.tumor$r[Hugo_Symbol, Hugo_Symbol][upper.tri(cor.tumor$r[Hugo_Symbol, Hugo_Symbol])]
    
    SCORE=mean(abs(cancer.vec - norm.vec)>alpha)
    
    return(list(SCORE=SCORE))
  }
  
  met.scores<-recon.table[,internal.function(Hugo_Symbol), by="ID"]
  
  #Assign gene counts
  met.scores<-merge(met.scores, gene.counts, by="ID")
  
  #Clean up and return
  met.scores<-met.scores[order(SCORE, decreasing = T),]
  return(met.scores)
}

MH.Scores<-Function.MH.Exp.Score(recon.table, brca.ma,"PIPELINES/METABOLIC.DRIVERS/TABLES/MH.HMDB.KEGG.csv",
                                 "PIPELINES/METABOLIC.DRIVERS/TABLES/060515.HMDB.KEGG.RESULTS",alpha = 0.2)
MH.Scores[1:20,]
hist(MH.Scores$SCORE)
for (i in seq(0.1,0.9,0.1)){
  print (i)
  MH.Scores<-Function.MH.Exp.Score(recon.table, brca.ma,"PIPELINES/METABOLIC.DRIVERS/TABLES/MH.HMDB.KEGG.csv",
                                   "PIPELINES/METABOLIC.DRIVERS/TABLES/060515.HMDB.KEGG.RESULTS",alpha = i)
  i.name<-paste("MH.Scores", i, sep=".")
  assign(i.name, MH.Scores)
  
}

#Obtain p-values for scores
Function.MH.Exp.Score.PVAL<-function(MH.Table, MH.Null){
  #Make sure both table and null were calculate using the same ALPHA and same EXPRESSION matrix!!
  
  #Load null
  null<-readRDS(MH.Null)
  
  #Calculate Empirical P-values
  MH.Table[,PVAL:= mean(null[[as.character(N.GENES)]]>=SCORE) , by="ID"]
  
  #Adjust for fdr
  MH.Table$PVAL.ADJ<-p.adjust(MH.Table$PVAL, method="fdr")
  
  #Clean up and return
  MH.Table<-MH.Table[order(PVAL.ADJ),]
  return(MH.Table)
}
MH.Scores<-Function.MH.Exp.Score.PVAL(MH.Scores.0.2, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/060415.BRCA.MA.MET.SCORE.NULL.0.2.rds")
MH.Scores[PVAL.ADJ<0.05,]
hist(MH.Scores$PVAL.ADJ)
ggplot(MH.Scores, aes(N.GENES, -log(PVAL.ADJ+0.001))) + geom_point() + theme.format

#Look for enrichment in pathway dataset
kegg.path<-fread("DATABASES/KEGG/060415_PATHWAY_TO_COMPOUND", header=T, sep="\t")
kegg.path

Function.kegg.mh.map<-function(mh.scores, kegg.path){
  
  #Filter for significant metabolites
  mh.scores<-mh.scores[PVAL.ADJ<0.05,]
  
  #Remove "- Homo sapiens" from path description
  kegg.path$DESCRIPTION<-sapply(kegg.path$DESCRIPTION, function(x) unlist(strsplit(x, " - Homo"))[1])
  
  #Remove "Metabolic pathways" from kegg.path
  kegg.path<-kegg.path[!grepl("metabolic", DESCRIPTION, ignore.case = T),]
  
  #Merge tables
  setnames(kegg.path, c("PATHWAY", "ID", "DESCRIPTION"))
  main.table<-merge(kegg.path[,c("ID", "DESCRIPTION"),with=F], mh.scores[,c("ID", "SCORE"),with=F], by = "ID")
  
  #Convert to heatmap matrix
  main.matrix<-acast(main.table, ID~DESCRIPTION, fill = 0)
  
  #Return
  return(main.matrix)
}
  
mh.path.map<-Function.kegg.mh.map(MH.Scores, kegg.path)
pheatmap(mh.path.map, scale = "none")

Function.mh.kegg.path.enrich<-function(recon.table, brca.exp, mh.hmdb.kegg, HMDB.table ,MH.Scores, kegg.path){
  #Function to find enriched kegg pathways based on significant metabolic hubs
  
  ##############First find all the possible metabolites that we could withdraw from#############
  #Double check and clean up expression datasets
  brca.exp$normal<-brca.exp$normal[complete.cases(brca.exp$normal),]
  brca.exp$tumor<-brca.exp$tumor[complete.cases(brca.exp$tumor),]
  
  #Get unique identifiers for recon table (Use primarily kegg id, only use hmdb id if kegg id not present)
  kegg.ids<-recon.table[KEGG_ID!="NONE",]
  hmdb.ids<-recon.table[KEGG_ID=="NONE" & HMDB!="NONE",]
  
  kegg.ids<-kegg.ids[,c("Hugo_Symbol", "KEGG_ID"),with=F]
  setnames(kegg.ids, c("Hugo_Symbol", "ID"))
  hmdb.ids<-hmdb.ids[,c("Hugo_Symbol", "HMDB"), with=F]
  setnames(hmdb.ids, c("Hugo_Symbol", "ID"))
  
  recon.table<-unique(rbind(kegg.ids, hmdb.ids))
  print ("Processed recon table")
  
  #Regularize some HMDB IDs to KEGG IDs
  id.change<-fread(mh.hmdb.kegg, header=T)
  id.change<-id.change[KEGG!="Not Found",]
  z<-id.change$KEGG
  names(z)<-id.change$HMDB
  
  recon.table$ID<-ifelse(recon.table$ID %in% names(z), z[recon.table$ID], recon.table$ID)
  
  #Add HMDB processed KEGG targeted table
  hmdb.source<-fread(HMDB.table, header=T, sep="\t", stringsAsFactors = F)
  hmdb.source<-hmdb.source[,c("KEGG_ID", "Hugo_Symbol"),with=F]
  setnames(hmdb.source, c("ID", "Hugo_Symbol"))
  recon.table<-unique(rbind(recon.table, hmdb.source[,c("Hugo_Symbol","ID"),with=F]))
  
  #Filter recon.table and brca.exp for common genes
  common.genes<-intersect(unique(recon.table$Hugo_Symbol), intersect(rownames(brca.exp$tumor),rownames(brca.exp$normal)))
  recon.table<-recon.table[Hugo_Symbol %in% common.genes,]
  brca.exp$tumor<-brca.exp$tumor[common.genes,]
  brca.exp$normal<-brca.exp$normal[common.genes,]
  
  #Filter recon.table for mets that have at least 3 genes assigned to them
  print (length(unique(recon.table$ID)))
  recon.count<-recon.table[,list(N.GENES=length(Hugo_Symbol)), by ="ID"]
  recon.count<-recon.count[N.GENES>=3, ]
  recon.table<-recon.table[ID %in% recon.count$ID,]
  print (length(unique(recon.table$ID)))
  ################################################################################################
  
  ##########Process kegg pathway table, filter for available compounds in test########
  kegg.path<-kegg.path[COMPOUND %in% unique(recon.table$ID),]
  kegg.path[,N.KEGGS:=length(unique(COMPOUND)), by=c("PATHWAY", "DESCRIPTION")]
  kegg.cpd<-unique(kegg.path$COMPOUND)
  
  #######Calculates p-values#######
  MH.Scores<-MH.Scores[PVAL.ADJ<0.05,]
  MH.MET<-unique(MH.Scores$ID)
  MH.MET<-MH.MET[MH.MET %in% kegg.path$COMPOUND]
  
  main.table<-kegg.path[,list(PVAL=phyper(q = length(intersect(MH.MET, COMPOUND))-1, 
                              m =  length(COMPOUND), 
                              n = length(kegg.cpd)-length(COMPOUND),
                              k = length(unique(MH.MET)), lower.tail=F )  ), by="DESCRIPTION"]
  main.table$PVAL.ADJ<-p.adjust(main.table$PVAL, method="fdr")
  
  ####Clean up and Return####
  main.table<-main.table[order(PVAL.ADJ),]
  return(main.table)
}

MH.Enrich<-Function.mh.kegg.path.enrich(recon.table, brca.ma, "PIPELINES/METABOLIC.DRIVERS/TABLES/MH.HMDB.KEGG.csv",
                                        "PIPELINES/METABOLIC.DRIVERS/TABLES/060515.HMDB.KEGG.RESULTS", MH.Scores, kegg.path )
View(MH.Enrich[PVAL.ADJ<0.05,])
hist(MH.Enrich$PVAL.ADJ)

#########Inter Functions########
Function.Clean.Recon<-function(recon.table, mh.hmdb.kegg, HMDB.table){
  
  #Pre-clean up
  recon.table<-do.call(rbind, recon.table)
  recon.table<-recon.table[CHEBI!="NONE",]
  recon.table<-recon.table[,c("Hugo_Symbol", "NAME", "KEGG_ID", "HMDB"), with =F]
  setkey(recon.table)
  recon.table<-unique(recon.table)
  recon.table<-recon.table[!grepl("_", Hugo_Symbol),]
  recon.table<-recon.table[!(KEGG_ID=="NONE" & HMDB=="NONE"),] 
  
  #Get unique identifiers for recon table (Use primarily kegg id, only use hmdb id if kegg id not present)
  kegg.ids<-recon.table[KEGG_ID!="NONE",]
  hmdb.ids<-recon.table[KEGG_ID=="NONE" & HMDB!="NONE",]
  
  kegg.ids<-kegg.ids[,c("Hugo_Symbol", "KEGG_ID"),with=F]
  setnames(kegg.ids, c("Hugo_Symbol", "ID"))
  hmdb.ids<-hmdb.ids[,c("Hugo_Symbol", "HMDB"), with=F]
  setnames(hmdb.ids, c("Hugo_Symbol", "ID"))
  
  recon.table<-unique(rbind(kegg.ids, hmdb.ids))
  print ("Processed recon table")
  
  #Regularize some HMDB IDs to KEGG IDs
  id.change<-fread(mh.hmdb.kegg, header=T)
  id.change<-id.change[KEGG!="Not Found",]
  z<-id.change$KEGG
  names(z)<-id.change$HMDB
  
  recon.table$ID<-ifelse(recon.table$ID %in% names(z), z[recon.table$ID], recon.table$ID)
  
  #Add HMDB processed KEGG targeted table
  hmdb.source<-fread(HMDB.table, header=T, sep="\t", stringsAsFactors = F)
  hmdb.source<-hmdb.source[,c("KEGG_ID", "Hugo_Symbol"),with=F]
  setnames(hmdb.source, c("ID", "Hugo_Symbol"))
  recon.table<-unique(rbind(recon.table, hmdb.source[,c("Hugo_Symbol","ID"),with=F]))
  
  #Return
  return(recon.table)
}

cleaned.recon<-Function.Clean.Recon(recon.table, "PIPELINES/METABOLIC.DRIVERS/TABLES/MH.HMDB.KEGG.csv", 
                                    "PIPELINES/METABOLIC.DRIVERS/TABLES/060515.HMDB.KEGG.RESULTS")
saveRDS(cleaned.recon, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061215.BRCA.MA.CLEANED.RECON.rds")

Function.Clean.Corr<-function(cleaned.recon, brca.exp, p.sig=0.05){
  
  #Calculate spearman correlation matrices at significance p-value threshold after fdr-correction
  #Correlation matrices calculated for genes present in cleaned.recon only
  
  #Double check and clean up expression datasets
  brca.exp$normal<-brca.exp$normal[complete.cases(brca.exp$normal),]
  brca.exp$tumor<-brca.exp$tumor[complete.cases(brca.exp$tumor),]
  
  #Filter recon table and expression matrices for common genes
  common.genes<-intersect(cleaned.recon$Hugo_Symbol, intersect(rownames(brca.exp$normal), rownames(brca.exp$tumor)))
  recon.table<-cleaned.recon[Hugo_Symbol %in% common.genes,]
  brca.exp$normal<-brca.exp$normal[common.genes, ]
  brca.exp$tumor<-brca.exp$tumor[common.genes, ]
  
  #Calculate correlation matrix in cancer and normal
  cor.normal<-rcorr(t(data.matrix(brca.exp$normal)), type = "spearman")
  cor.tumor<-rcorr(t(data.matrix(brca.exp$tumor)), type = "spearman")
  print ("Done with correlation tables")
  
  #Correct correlation p-values for multiple hypothesis testing and obtain their signficance status (P-val<0.05)
  cor.normal$P<-matrix(p.adjust(cor.normal$P, method = "fdr"), ncol=ncol(cor.normal$P), 
                       dimnames = list(rownames(cor.normal$P), colnames(cor.normal$P)))
  cor.tumor$P<-matrix(p.adjust(cor.tumor$P, method = "fdr"), ncol=ncol(cor.tumor$P), 
                      dimnames = list(rownames(cor.tumor$P), colnames(cor.tumor$P)))
  cor.normal$P<-cor.normal$P < p.sig
  cor.tumor$P<-cor.tumor$P < p.sig
  diag(cor.normal$P)<-FALSE 
  diag(cor.tumor$P)<-FALSE
  
  #Use corrrelation value only if it passes multiple hypothesis testing (Self correlation value will be zero)
  cor.normal$r<-cor.normal$r * cor.normal$P
  cor.tumor$r<-cor.tumor$r * cor.tumor$P
  
  #Clean up and return
  return(list(COR.CANCER=cor.tumor$r, COR.NORMAL=cor.normal$r))
}

cleaned.cor<-Function.Clean.Corr(cleaned.recon, brca.ma, 0.05)
saveRDS(cleaned.cor, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061215.BRCA.MA.CLEANED.COR.rds")

#Calculates disturbance in inter eucledian distance from MH to MH and obtain score per MH dysregulation
Function.Inter.MH.Score<-function(cleaned.recon, cleaned.cor, beta){
  
  require(parallel)
  
  #Filter recon table and correlation matrices for common genes
  common.genes<-intersect(cleaned.recon$Hugo_Symbol, intersect(rownames(cleaned.cor$COR.NORMAL), rownames(cleaned.cor$COR.CANCER)))
  recon.table<-cleaned.recon[Hugo_Symbol %in% common.genes,]
  cor.normal<-cleaned.cor$COR.NORMAL[common.genes, common.genes]
  cor.tumor<-cleaned.cor$COR.CANCER[common.genes, common.genes]
  
  #Filter recon.table for mets that have at least 3 genes assigned to them
  print (length(unique(recon.table$ID)))
  recon.count<-recon.table[,list(N.GENES=length(Hugo_Symbol)), by ="ID"]
  recon.count<-recon.count[N.GENES>=3, ]
  recon.table<-recon.table[ID %in% recon.count$ID,]
  print (length(unique(recon.table$ID)))
  
  #Prep parallelization
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "cor.normal", "cor.tumor",
                              "recon.table", "beta") ,envir=environment())
  
  #Obtain Inter-dysregulation score per metabolite
  print ("Calculating scores")
  recon.mets<-unique(recon.table$ID)
  
  main.scores<-parSapply(cl, recon.mets, function(x) {
    met<-x
    Hugo.target<-unique(recon.table[ID==met, ]$Hugo_Symbol)
    other.hugos<-setdiff(unique(recon.table$Hugo_Symbol), Hugo.target)
    
    #Calculate dysregulation per MH node enzyme to whole metabolic network
    hub.scores<-abs(cor.tumor[Hugo.target, other.hugos] - cor.normal[Hugo.target, other.hugos])
    hub.scores<-apply(hub.scores, 1, median)
    
    INTER.SCORE<-mean(hub.scores>beta)
    return(INTER.SCORE)
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelization")
  
  main.table<-data.table(ID=recon.mets, INTER.SCORE=main.scores)
  
  #Add gene count info
  main.table<-merge(main.table, recon.count, by ="ID")
  
  #Clean up and return
  main.table<-main.table[order(INTER.SCORE, decreasing=T),]
  return(main.table)
}

MH.Inter.Scores<-Function.Inter.MH.Score(cleaned.recon, cleaned.cor, 0.15)
MH.Inter.Scores[1:20,]
hist(MH.Inter.Scores$INTER.SCORE)
ggplot(MH.Inter.Scores, aes(N.GENES, INTER.SCORE)) + geom_point() + theme.format + scale_x_log10() 

s.inter<-data.table()
for (s in c(0.1, 0.12, 0.15, 0.17, 0.2)){
  inter.mh<-Function.Inter.MH.Score(cleaned.recon, cleaned.cor, s)
  inter.mh$beta<-s
  s.inter<-rbind(s.inter, inter.mh)
}
ggplot(s.inter, aes(N.GENES, INTER.SCORE, colour=beta)) + geom_point() + theme.format + scale_x_log10() 

Function.Inter.Scores.PVAL<-function(MH.Inter.Scores, Inter.Null){
  #Make sure both table and null were calculate using the same ALPHA and same EXPRESSION matrix!!
  
  #Load null
  null<-readRDS(Inter.Null)
  
  #Calculate Empirical P-values
  MH.Inter.Scores[,PVAL:= mean(null[[as.character(N.GENES)]]>=INTER.SCORE) , by="ID"]
  
  #Adjust for fdr
  MH.Inter.Scores$PVAL.ADJ<-p.adjust(MH.Inter.Scores$PVAL, method="fdr")
  
  #Clean up and return
  MH.Inter.Scores<-MH.Inter.Scores[order(PVAL.ADJ),]
  return(MH.Inter.Scores)
}
MH.Inter.Scores.PVAL<-Function.Inter.Scores.PVAL(MH.Inter.Scores, "PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/061215.BRCA.MA.MET.INTER.NULL.0.15.rds")
hist(MH.Inter.Scores.PVAL$PVAL.ADJ)
MH.Inter.Scores.PVAL[PVAL.ADJ<0.05,]
ggplot(MH.Inter.Scores.PVAL, aes(N.GENES, -log(PVAL.ADJ+0.01))) + geom_point() + scale_x_log10() + theme.format

#Distribution of cancer compounds in cancer
kegg.cancer<-kegg.path[grepl("cancer", DESCRIPTION,ignore.case = T ) ,]
kegg.cancer$DESCRIPTION<-sapply(kegg.cancer$DESCRIPTION, function(x) strsplit(x, " - ")[[1]][1])
kegg.cancer$PRESENCE<-1
kegg.cancer<-acast(kegg.cancer, DESCRIPTION~COMPOUND, value.var = "PRESENCE", fill = 0)
pheatmap(kegg.cancer, dualScale = F, color=colorRampPalette(c("white","green","green4","violet","purple"))(100),fontsize = 14)

unique(kegg.path[grepl("cancer", DESCRIPTION,ignore.case = T ) ,]$DESCRIPTION)
unique(kegg.path[grepl("cancer", DESCRIPTION,ignore.case = T ) ,]$COMPOUND)
kegg.path[DESCRIPTION=="Pancreatic cancer - Homo sapiens (human)",]
kegg.path[DESCRIPTION=="Prostate cancer - Homo sapiens (human)",]
kegg.path[DESCRIPTION=="Pathways in cancer - Homo sapiens (human)",]


recon.table$SUBSTRATE[Hugo_Symbol %in% c("AKT1", "ARID1A", "ARID1B", "BAP1", "BRCA1", "BRCA2", "CASP8", "CCND1", "CDH1",
                                       "CDKN1B", "EP300", "ERBB2", "ETV6", "FOXA1", "GATA3", "MAP2K4", "MAP3K1", "MAP3K13", 
                                       "NCOR1", "NOTCH1", "NTRK3", "PBRM1", "PIK3CA", "RB1", "SMARCD1", "TBX3", "TP53"),]

###DEVELOP GOLD STANDARD#####
cleaned.recon
options( expressions = 5e5 )
Function.met.gold.std<-function(target.genes, met.table){
  #Function to build gold standard metabolite scoring metric
  #met.table preferably has 2 columns with the headers c("Hugo_Symbol", "ID")
  
  internal.function<-function(mets, wave, bag, main.met){
    
    #Store met waves in table
    current.wave<-wave+1
    sep.matrix[mets, main.met]<<-as.numeric(current.wave)
    
    #Remove them by assigning them to bag
    wave.bag<-c(bag, mets)
    
    print (c(current.wave, length(mets), length(wave.bag)))
    
    #Obtain wave's metabolites children if not already accounted for
    wave.mets.pre<-unique(met.uni[ID.1 %in% mets,]$ID.2)
    wave.mets.filt<-setdiff(wave.mets.pre, wave.bag)
    
    #Go to next wave to assign next degree of separation if there are still waves to go
    if (length(wave.mets.filt)>0){
      internal.function(wave.mets.filt, current.wave, wave.bag, main.met)
  }
  
  #Obtain cancer node metabolites based on their association to target.genes
  target.metabolites<-unique(met.table[Hugo_Symbol %in% target.genes,]$ID)
  all.metabolites<-unique(met.table$ID)
  
  #Assign initial score to target metabolites based on their number of associations 
  #This will penalize metabolites that while associated with cancer genes, do associated with many other genes as well
  #This is to take into account that metabolites that associate with many genes do not necessarily exert all their influence
  #   on cancer genes (i.e. ADP, NADH)
  n.met.assoc<-met.table[,list(N.GENES=length(unique(Hugo_Symbol))), by="ID"]
  score.matrix<-matrix(ncol=length(target.metabolites), nrow=length(all.metabolites), 
                       dimnames =  list(all.metabolites, target.metabolites))  
  
  #Assign degree of separation per target metabolite, how far in term of waves are all metabolites from cancer metabolite
  #First convert to unipartite network
  met.uni<-data.table()
  for (d in unique(met.table$ID)){
    
    hugos<-unique(met.table[ID==d,]$Hugo_Symbol)
    neighbor.met<-unique(met.table[Hugo_Symbol %in% hugos,]$ID)
    neighbor.met<-setdiff(neighbor.met, d)
    if (length(neighbor.met)>0){
      met.uni<-rbind(met.uni, data.table(ID.1=d, ID.2=neighbor.met))  
    } else {
      print (d) #Not all metabolites are connected to network
    }
  }
  
  #Then assign metabolite wave degree to separation matrix based on unipartite network
  sep.matrix<-matrix(ncol=length(target.metabolites), nrow=length(all.metabolites), 
                       dimnames =  list(all.metabolites, target.metabolites))  
  for (cancer.met in target.metabolites){
    print (cancer.met)
    
    #Assign wave number to cancer met
    sep.matrix[cancer.met, cancer.met]<-1
    
    #Obtain first propagation metabolites from parent
    met.parent.mets<-unique(met.uni[ID.1==cancer.met,]$ID.2)
    
    #Apply iteratively to next waves
    wave<-1
    wave.bag<-c(cancer.met) #store met
    internal.function(met.parent.mets, wave, wave.bag, cancer.met)
  }
  
  #Return score matrix
  return(sep.matrix)
}

met.gold.std<-Function.met.gold.std(c("AKT1", "ARID1A", "ARID1B", "BAP1", "BRCA1", "BRCA2", "CASP8", "CCND1", "CDH1",
                                      "CDKN1B", "EP300", "ERBB2", "ETV6", "FOXA1", "GATA3", "MAP2K4", "MAP3K1", "MAP3K13", 
                                      "NCOR1", "NOTCH1", "NTRK3", "PBRM1", "PIK3CA", "RB1", "SMARCD1", "TBX3", "TP53"), cleaned.recon)

met.gold.std[is.na(met.gold.std)]<-4
length(met.gold.std[is.na(met.gold.std)])/9
head(met.gold.std)
met.gold.std[c("C04185", "C00292"),]
data.table(met.gold.std,keep.rownames = T)[is.na(C04637),]
cleaned.recon.uni[ID.1=="C02355",]$ID.2
cleaned.recon.uni[ID.1 %in% cleaned.recon.uni[ID.1=="C02355",]$ID.2,]

cleaned.net<-graph.data.frame(cleaned.recon.uni)
par(mai=c(0,0,1,0))
plot(cleaned.net, layout=layout.auto)
dev.off()

hist(met.gold.std[,"C05981"])
pheatmap(met.gold.std[order(apply(met.gold.std, 1, median)),], scale="none", cluster_rows = F, fontsize = 14,
         color=colorRampPalette(c("white", "green","green4","violet","purple"))(100))

x<-matrix(ncol=2, nrow=5, dimnames = list(letters[1:5], c("A", "B")))
x
for (i in letters[3:4]){
  x[i,]<-i
}
length(x["a","A"])
x[c("b", "c"),"A"]<-2

z
plot(z$N.GENES, z$PROP)

cleaned.recon
cleaned.recon.uni<-data.table()
for (d in unique(cleaned.recon$ID)){
  
  hugos<-unique(cleaned.recon[ID==d,]$Hugo_Symbol)
  neighbor.met<-unique(cleaned.recon[Hugo_Symbol %in% hugos,]$ID)
  neighbor.met<-setdiff(neighbor.met, d)
  if (length(neighbor.met)>0){
    cleaned.recon.uni<-rbind(cleaned.recon.uni, 
                             data.table(ID.1=d, ID.2=neighbor.met))  
  } else {
    print (d)
  }
}

Function.jaccard<-function(a,b){
  j<-length(intersect(unique(a), unique(b))) / length(union(unique(a), unique(b)))
  return(j)
}
cleaned.recon.uni$JACCARD<-apply(as.matrix(cleaned.recon.uni), 1, function(x) 
  Function.jaccard(cleaned.recon[ID==x[1],]$Hugo_Symbol, cleaned.recon[ID==x[2],]$Hugo_Symbol))
hist(cleaned.recon.uni$JACCARD)
cleaned.recon.uni[ID.1=="C02355",]

w<-cleaned.recon.uni[,1:3,with=F]
w<-acast(w, ID.1~ID.2,fill = 0, value.var = "JACCARD")
heatplot(w, dualScale = F)
pheatmap(w, scale="none",color=colorRampPalette(c("white", "green","green4","violet","purple"))(100))

g<-random.graph.game(5, 5/5, directed=F)
V(g)$weight<-c(0,0,0,0,2)
V(g)$weight
plot(g)
page.rank(g)$vector

cbind(V(g)$weight, page.rank(g)$vector)


m<-matrix(c(0,0,1,0.5, 1/3, 0, 0, 0, 1/3, 0.5, 0, 0.5, 1/3, 0.5, 0, 0), ncol=4, byrow = T)
v<-c(0.25, 0.25, 0.25, 0.25) #weights of network

B<-matrix(ncol=ncol(m), nrow=nrow(m))
B[is.na(B)]<-1
A.1<-(m * (0.85) + (1-0.85)*(1/4)*(B))%*%v
delta<-Inf
while (delta>0.00000001){
  delta<-mean(abs(A.1 - (m * (0.85) + (1-0.85)*(1/4)*(B))%*%A.1))
  A.1<-(m * (0.85) + (1-0.85)*(1/4)*(B))%*%A.1
}
eigen(m)


#######CONSTRUCT DIFFUSION ALGORITHM THAT PREDICTS 
kegg.path
cleaned.recon.uni

t.m<-cleaned.recon.uni[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
t.m<-merge(t.m, cleaned.recon.uni[,c("ID.1", "ID.2"), with=F], by = "ID.1")
t.m$trans<-1/t.m$DEGREE
t.m$DEGREE<-NULL
t.m<-rbind(t.m, data.table(ID.1=unique(t.m$ID.1), ID.2=unique(t.m$ID.1), trans=0))
t.m<-acast(t.m, ID.2~ID.1, fill = 0)
heatplot(t.m, ="none")
pheatmap(t.m, scale="none", color=colorRampPalette(c("white", "green","green4","violet","purple"))(100))

Function.met.scoring<-function(network.table, kegg.path, target.path) {
  
  #network.table<-recon.directed
  #target.path<-random.path
  #Create transition matrix
  #t.m<-network.table[,list(DEGREE=length(unique(ID.2))), by="ID.1"]
  t.m<-network.table[,list(DEGREE=length(unique(PRODUCT))), by="SUBSTRATE"]
  #t.m<-merge(t.m, network.table[,c("ID.1", "ID.2"), with=F], by = "ID.1")
  t.m<-merge(t.m, network.table[,c("SUBSTRATE", "PRODUCT"), with=F], by = "SUBSTRATE")
  t.m$trans<-1/t.m$DEGREE #may modify depending on weighted jaccard edges
  #degree.info<-unique(data.table(ID=t.m$ID.1, DEGREE=t.m$DEGREE))
  degree.info<-unique(data.table(SUBSTRATE=t.m$SUBSTRATE, DEGREE=t.m$DEGREE))
  t.m$DEGREE<-NULL
  #t.m<-rbind(t.m, data.table(ID.1=unique(t.m$ID.1), ID.2=unique(t.m$ID.1), trans=0)) #To itself "0" transition
  t.m<-rbind(t.m, data.table(SUBSTRATE=unique(t.m$SUBSTRATE), PRODUCT=unique(t.m$SUBSTRATE), trans=0))
  #t.m<-acast(t.m, ID.2~ID.1, fill = 0)
  all.mets<-unique(c(t.m$SUBSTRATE, t.m$PRODUCT))
  adj.t.m<-matrix(ncol=length(all.mets), nrow=length(all.mets), dimnames = list(all.mets, all.mets))
  adj.t.m[is.na(adj.t.m)]<-0
  print (head(as.matrix(t.m)))
  apply(as.matrix(t.m), 1, function(x) {
    adj.t.m[x[2], x[1]]<<-as.numeric(x[3])
  })
  print (apply(adj.t.m, 2, sum))
  print (apply(adj.t.m, 1, sum))
  
  #t.m<-acast(t.m, PRODUCT~SUBSTRATE, fill = 0)
  print (dim(adj.t.m))
  
  #Obtain v vector (seeds) for path
  seeds<-unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  seeds<-ifelse(colnames(adj.t.m) %in% seeds, 1, 0)
  
  #####Iterate to minimize mean change per row and converge to stable network####
  #Generate google matrix
  B<-matrix(ncol=ncol(adj.t.m), nrow=nrow(adj.t.m))
  B[is.na(B)]<-1
  #g.m<-(t.m*0.85) + (1-0.85)*(1/nrow(t.m))*B
  
  #Obtain first iteration page rank
  #A.1<-g.m%*%seeds
  A.1<-(adj.t.m*0.85)%*%seeds + (1-0.85)*(1/nrow(adj.t.m))*rep(1, nrow(adj.t.m)) #Implementing personalized vector
  
  #Iterate to converge on optimal page rank
  delta<-Inf
  while(delta>0.0000000001){
    #delta<-mean(abs(A.1 - g.m%*%A.1))
    delta<-mean(abs(A.1 - ((adj.t.m*0.85)%*%A.1 + (1-0.85)*(1/nrow(adj.t.m))*rep(1, nrow(adj.t.m)))))
    A.1<-((adj.t.m*0.85)%*%A.1 + (1-0.85)*(1/nrow(adj.t.m))*rep(1, nrow(adj.t.m)))
    #A.1<-g.m%*%A.1
  }
  
  #Clean up and return
  main.table<-data.table(SUBSTRATE=rownames(adj.t.m), PAGE.RANK=as.vector(A.1))
  #main.table<-merge(main.table, degree.info, by="ID")
  print (main.table)
  print (degree.info)
  main.table<-merge(main.table, degree.info, by="SUBSTRATE")
  #main.table$TARGET.PATH<-main.table$ID %in% unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  main.table$TARGET.PATH<-main.table$SUBSTRATE %in% unique(kegg.path[DESCRIPTION==target.path,]$COMPOUND)
  main.table<-main.table[order(PAGE.RANK, decreasing = T),]
  return(main.table)
}

page.rank.9<-data.table()
for (random.path in sample(unique(kegg.path$DESCRIPTION), 9)){
  page.rank.1<-Function.met.scoring(cleaned.recon.uni, kegg.path, random.path)
  page.rank.1$PATH<-unlist(strsplit(random.path, " - "))[1]
  page.rank.9<-rbind(page.rank.9, page.rank.1)
}
page.rank.9
ggplot(page.rank.9, aes(DEGREE, PAGE.RANK)) + geom_point() + theme.format + facet_wrap(~PATH, scales = "free_y") +theme(strip.text.x = element_text(size = 12))
ggplot(page.rank.9, aes(TARGET.PATH, PAGE.RANK)) + geom_boxplot() + theme.format + facet_wrap(~PATH, scales = "free_y") + scale_y_log10() +
  theme(strip.text.x = element_text(size = 12))


Funcion.directed.network<-function(recon.table) {
  
  #Pre-clean up
  recon.product<-recon.table$PRODUCT[,c("Hugo_Symbol", "KEGG_ID"), with=F]
  recon.product<-recon.product[!grepl("_", Hugo_Symbol),]
  recon.product<-recon.product[KEGG_ID!="NONE",]
  recon.product<-unique(recon.product)
  
  recon.substrate<-recon.table$SUBSTRATE[,c("Hugo_Symbol", "KEGG_ID"), with=F]
  recon.substrate<-recon.substrate[!grepl("_", Hugo_Symbol),]
  recon.substrate<-recon.substrate[KEGG_ID!="NONE",]
  recon.substrate<-unique(recon.substrate)
  
  #Form edge network
  edge.net<-data.table()
  for (d in unique(recon.substrate$KEGG_ID)){
    hugos<-unique(recon.substrate[KEGG_ID==d,]$Hugo_Symbol)
    product.mets<-recon.product[Hugo_Symbol %in% hugos,]$KEGG_ID
    product.mets<-setdiff(product.mets, d)
    
    if (length(product.mets)>0){
      edge.net<-rbind(edge.net, data.table(SUBSTRATE=d, PRODUCT=product.mets))
    }
  }
  
  #Calculate jaccard coefficient for each edge (first column=substrate, second column=product)
  Function.jaccard<-function(a,b){
    j<-length(intersect(unique(a), unique(b))) / length(union(unique(a), unique(b)))
    return(j)
  }
  edge.net$JACCARD<-apply(as.matrix(edge.net), 1, function(x) 
    Function.jaccard(recon.substrate[KEGG_ID==x[1],]$Hugo_Symbol, recon.substrate[KEGG_ID==x[2],]$Hugo_Symbol))
  
  #Clean up and return
  edge.net<-edge.net[order(JACCARD,decreasing = T),]
  return(edge.net)
}

recon.directed<-Funcion.directed.network(recon.table)
hist(recon.directed$JACCARD)

recon.directed.graph<-graph.data.frame(recon.directed,directed = T)
par(mai=c(0,0,1,0))
plot(recon.directed.graph)
dev.off()

page.rank.9<-data.table()
for (random.path in sample(unique(kegg.path$DESCRIPTION), 9)){
  page.rank.1<-Function.met.scoring(recon.directed, kegg.path, random.path)
  page.rank.1$PATH<-unlist(strsplit(random.path, " - "))[1]
  page.rank.9<-rbind(page.rank.9, page.rank.1)
}
ggplot(page.rank.9, aes(DEGREE, PAGE.RANK)) + geom_point() + theme.format + facet_wrap(~PATH, scales = "free_y") +theme(strip.text.x = element_text(size = 12)) +
  scale_y_log10() + scale_x_log10()
ggplot(page.rank.9, aes(TARGET.PATH, PAGE.RANK)) + geom_boxplot() + theme.format + facet_wrap(~PATH, scales = "free_y") + scale_y_log10() +
  theme(strip.text.x = element_text(size = 12))

page.rank.9[PATH=="Melanogenesis",]
cor.test(page.rank.9[PATH=="Melanogenesis",]$PAGE.RANK, page.rank.9[PATH=="Melanogenesis",]$DEGREE, method="spearman")
