#Paper 1 Scripts
library(data.table)
library(reshape2)
library(parallel)
library(grid)
library(gridExtra)

Function_cgp_process <- function(exp_file, act_file){
	# exp_file as in sanger1018_brainarray_ensemblgene_rma.txt
	# act_file as in v17.3_fitted_dose_response.csv
 
	require(data.table)

	# DO ACTIVITY TABLE
	act_table       <- fread(act_file, header=T)[,c("COSMIC_ID", "CELL_LINE_NAME", "DRUG_NAME", "LN_IC50", "AUC"), with=F]
	setnames(act_table, c("COSMIC_ID", "cell_name", "Compound", "LN_IC50", "AUC"))

	act_table$IC50  <- exp(act_table$LN_IC50)
	act_table$pIC50 <- -log(act_table$IC50)

	# Remove duplicates activities (mean is calculated for now)
	act_table  <- act_table[,list(IC50      = mean(IC50),
	                           	 pIC50      = mean(pIC50),
	                             AUC        = mean(AUC),
	                             COSMIC_ID  = COSMIC_ID), 
							by=c("cell_name", "Compound")]

	# Clean up
	setkey(act_table)
	act_table <- unique(act_table)

	# DO EXPRESSION TABLE
	exp_table      <- fread(exp_file, header=T)
	exp_samples    <- exp_table$ensembl_gene
	exp_table      <- as.matrix(exp_table[,2:ncol(exp_table),with=F]) 
	rownames(exp_table) <- exp_samples

	common_samples <- unique(act_table[,c("cell_name","COSMIC_ID"),with=F])
	common_samples <- common_samples[COSMIC_ID %in% colnames(exp_table),]
	exp_table      <- exp_table[,as.character(common_samples$COSMIC_ID)]
	colnames(exp_table) <- common_samples$cell_name

	genes          <- Function_ensembl_to_hugo(rownames(exp_table))
	exp_table      <- exp_table[genes$ensembl_gene_id,]
	rownames(exp_table) <- genes$hgnc_symbol

	# FILTER
	act_table      <- act_table[cell_name %in% colnames(exp_table)]

	# Return
	return(list(act_table=act_table, exp_table=exp_table))
}

Function_cgp_new <- function(cgp_file, cgp_cell, cgp_drug, norm_resp="Compound") {
	# cgp_file as in "v17_fitted_dose_response.csv"
	# cgp_cell as in "Cell_Lines_Details.csv"
	# cgp_drug as in "Screened_Compounds.csv"
	# Files can be directly obtained from http://www.cancerrxgene.org/downloads

	act_table  <- fread(cgp_file, header = T)[,c("COSMIC_ID", "DRUG_ID", "LN_IC50", "AUC"), with=F]
	cell_table <- fread(cgp_cell, header = T)[,c("Line", "COSMIC_ID", "Site"),with=F]
	drug_table <- fread(cgp_drug, header = T)[,c("DRUG ID", "DRUG NAME"),with=F]
	setnames(drug_table, c("DRUG_ID", "DRUG NAME"))

	setkey(act_table)
	act_table  <- unique(act_table)
	act_table  <- merge(act_table, cell_table, by = c("COSMIC_ID"))
	act_table  <- merge(act_table, drug_table, by = c("DRUG_ID"))

	act_table$IC50  <- exp(act_table$LN_IC50)
	act_table$pIC50 <- -log(act_table$IC50)

	act_table  <- act_table[ ,c("Line", "DRUG NAME", "IC50", "pIC50", "AUC", "Site"), with=F]
	setnames(act_table, c("cell_name", "Compound", "IC50", "pIC50", "AUC", "Site"))

	# Remove duplicates activities (mean is calculated for now)
	act_table  <- act_table[,list(IC50  = mean(IC50),
	                           	 pIC50 = mean(pIC50),
	                             AUC   = mean(AUC),
	                             Site  = Site), 
							by=c("cell_name", "Compound")]

	# Clean up and Return
	setkey(act_table)
	act_table <- unique(act_table)
	return(act_table)
} 

Function_cgp_new_exp <- function(cgp_file, cgp_cell) {
	# cgp_file: RMA normalized CGP expression dataset
	# cgp_cell: Cell line details. Second panel of Cell_Lines_Details
	# Files can be directly obtained from http://www.cancerrxgene.org/downloads

	require(data.table)

	x <- fread(cgp_file, header=T, drop=2)
	x <- x[GENE_SYMBOLS!="",]

	x_labels <- x$GENE_SYMBOLS
	x$GENE_SYMBOLS <- NULL

	x <- as.matrix(x)
	rownames(x) <- x_labels

	#Replace cell identifiers with names
	cell_names  <- sapply(colnames(x), function(z) strsplit(z, "DATA.")[[1]][2])
	colnames(x) <- cell_names

	cell_table  <- fread(cgp_cell, header = T)[,c("Line", "COSMIC_ID"),with=F]
	cell_table$COSMIC_ID <- as.character(cell_table$COSMIC_ID) 
	cell_names  <- intersect(cell_names, cell_table$COSMIC_ID)

	x <- x[ ,cell_names]

	colnames(x) <- sapply(colnames(x), function(z)  unique(cell_table[COSMIC_ID==z,]$Line))

	#Return
	return(x)
}

Function_ensembl_to_hugo <- function(genes){
	require(biomaRt)
	require(data.table)

	# Get Hugo Identifiers
	mart     <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	g_list   <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
					  values=genes, mart= mart)
	g_list   <- data.table(g_list)
	g_list   <- g_list[hgnc_symbol!=""]

	# Verified duplicated
	g_list   <- g_list[hgnc_symbol!="LINC00856"]

	return(g_list)
}

Function_process_drug_5.0 <-function(fit_data_file, gdsc_out_file){
	# fit_data_file: as in gdsc_drug_sensitivity_fitted_data_w5.csv
	# gdsc_out_file: as in gdsc_en_output_w5.csv
	require(data.table)

	x <- fread(fit_data_file, sep=",", select=c("cell_line_name", "drug_id", "ic_50_est", "auc"))
	y <- fread(gdsc_out_file, sep=",", select=c("DRUG ID", "DRUG NAME"))
	y <- unique(y)

	x <- merge(x, y, by.x="drug_id", by.y="DRUG ID")
	x <- x[,c("cell_line_name", "DRUG NAME", "ic_50_est", "auc"),with=F]
	setnames(x, c("cell_name", "Compound", "IC50", "AUC")) # This values is probably the natural log of IC50, and not just IC50
	return(x)
}

Function_process_drug_2016 <- function(fit_file, cell_line_file, compound_file){
	# fit_file: as in v17_fitted_dose_response.csv
	# cell_line_file: as in Cell_Lines_Details.csv
	# compound_file: as in Screened_Compounds.csv
	require(data.table)

	x <- fread(fit_file, sep=",", select=c("COSMIC_ID", "DRUG_ID", "LN_IC50", "AUC"))
	y <- fread(cell_line_file, sep=",", select=c("Line", "COSMIC_ID", "Site"))
	z <- fread(compound_file, sep=",", select = c("Drug ID", "Drug Name"))

	x <- merge(x, y, by="COSMIC_ID")
	x <- merge(x, z, by.x="DRUG_ID", by.y="Drug ID")

	x <- x[,c("LN_IC50", "AUC", "Line", "Site", "Drug Name"),with=F]
	setnames(x, c("LN_IC50", "AUC", "cell_name", "Site", "Compound"))

	x$IC50 <- exp(x$LN_IC50)
	return(x)
}

Function_rma_normalize_affy <- function(affy_object){
	# RMA Normalizes AffyBatch object form CGP
	# CGP expression object: loaded E-MTAB-783.eSet.r into environment
	# Files obtained from: https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-783
	require(affy)
	require(hgu133a.db) 

	# RMA normalization
	rma_obj <- rma(affy_object)
	rma_exp <- exprs(rma_obj)

	#Assign identifiers and HUGO Symbols
	sample_names      <- data.table(rma_obj@phenoData@data[,1:2], keep.rownames = T)

	Symbols           <- unlist(mget(row.names(rma_exp), hgu133aSYMBOL, ifnotfound=NA))
	rma_exp           <- rma_exp[!is.na(Symbols),]
  	rownames(rma_exp) <- Symbols[!is.na(Symbols)]

  	rma_exp 		  <- data.table(rma_exp,keep.rownames = T)
	rma_exp 		  <- rma_exp[,lapply(.SD, mean), by=rn]
	rma_exp 		  <- data.frame(rma_exp,row.names = 1)
	rma_exp 		  <- data.matrix(rma_exp)

	colnames(rma_exp) <- sample_names$Source.Name
	return(rma_exp)
}

Function_process_geeleher_brain_array <- function(loaded_object, sdrf_file){
	# Script to process gene expression object obtained from Geeleher et al. to correct identifiers
	# loaded_object: loaded from gdsc_brainarray_syms.RData as gdsc_brainarray_syms
	# sdrf_file    : as in E-MTAB-783.sdrf.txt

	require(data.table)
	
	x <- fread(sdrf_file, sep="\t", select = c("Array Data File","Factor Value[CELL_LINE]"))
	setnames(x, c("cel_id", "cell_name"))
	
	common_cells  <- intersect(colnames(loaded_object), x$cel_id)
	x 			  <- x[cel_id %in% common_cells,]
	loaded_object <- loaded_object[,x$cel_id]
	colnames(loaded_object) <- x$cell_name

	# Return
	return(loaded_object)
}

Function_std_gene <- function(x, f=T) {
	# genes have to be in rows
	if (f==T){
		cat("applying gene normalization\n")
		std <- apply(x, 1, sd)
		mu  <- apply(x, 1, mean)
		x   <- (x-mu)/std
	}
	return(x)
}

Function_std_sample <- function(x, f=T) {
	# samples have to be in columns
	if (f==T){
		cat("applying sample normalization\n")
		std <- apply(x, 2, sd)
		mu  <- apply(x, 2, mean)

		x   <- sweep(x, 2, mu, "-")
		x   <- sweep(x, 2, std, "/")
	}
	return(x)
}

Function_common_genes <- function(x,y){
	common_genes <- intersect(rownames(x), rownames(y))
	return(list(x[common_genes,], y[common_genes,]))
}

Function_load_gsea <- function(gsea, gene_th, in_folder){
	# Gene sets obtained from: http://software.broadinstitute.org/gsea/downloads.jsp
	# Store them @ GSEA_FILES/ folder
	g_set      <- strsplit(gsea, "_")[[1]]

	gsea_all   <- data.table()
	for (g in g_set){
		cat(paste0(g, "\n"))
		gsea_set   <- fread(paste0(in_folder, "GSEA_FILES/", g, "_sets"))
		gsea_all   <- rbind(gsea_all, data.table(gsea_set, source=g))
	}

	gsea_all[,N:=length(unique(genes)), by="Gene_set"]
	gsea_all   <- gsea_all[N >= gene_th,]
	gsea_all$N <- NULL

	return(gsea_all)
}

Function_common_gsea <- function(gsea){
	require(parallel)
	# Gets common gene stats between two gene sets
	gene_sets  <- unique(gsea$Gene_set)
	gene_sets  <- data.table(t(combn(gene_sets,2)))

	# Setup parallelization
	nodes          <- detectCores()
	cat(paste0("Number of available nodes detected: ", nodes, "\n"))
	cl   		   <- makeCluster(nodes)
	setDefaultCluster(cl)
	clusterExport(cl, varlist=c("data.table", "as.data.table", "gene_sets", "gsea"),
				  envir=environment())
	cat("Done exporting to parallel\n")

	# Run parallel
	sets_inter <- parApply(cl, gene_sets, 1, function(x)  length(intersect(gsea[Gene_set==x[1]]$genes, gsea[Gene_set==x[2]]$genes)))
	stopCluster(cl)

	# Finish up with stats
	gene_sets$common <- sets_inter

	gsea[,n:=length(genes), by="Gene_set"]
	gene_sets  <- merge(gene_sets, unique(gsea[,c("Gene_set", "source", "n"),with=F]), by.x="V1", by.y="Gene_set")
	gene_sets  <- merge(gene_sets, unique(gsea[,c("Gene_set", "source", "n"),with=F]), by.x="V2", by.y="Gene_set")
	setnames(gene_sets, c("Gene_set_1", "Gene_set_2", "n_common", "source_1", "n_1", "source_2", "n_2"))

	gene_sets$tanimoto <- gene_sets$n_common / (gene_sets$n_1 + gene_sets$n_2 - gene_sets$n_common)

	# Output known at /tigress/zamalloa/GSEA_FILES/gene_sets_c2_common_count.rds
	return(gene_sets)
}

Function_combine_th <- function(gsea, common_table){
	# common_table obtained from Function_common_gsea()
	# removes gene_sets and specified tanimoto threshold and calculates lefover number of genes
	# Obtain coverage

	tanimoto_th <- c(seq(0.01, 0.09, 0.01), seq(0.1, 1, 0.05))

	genes_n  <- lapply(tanimoto_th, function(x) {
		print(x)
		common_sets   <- with(common_table[tanimoto>x], unique(c(Gene_set_1, Gene_set_2)))
		leftover_sets <- with(gsea, setdiff(Gene_set, common_sets))
		covered_genes <- length(unique(gsea[!Gene_set %in% common_sets]$genes))
		return(data.table(Tanimoto_th=x, Genes_n=covered_genes, Sets_n=length(leftover_sets) )  )
	})
	genes_n  <- do.call(rbind, genes_n)
	return(genes_n)
}

Function_combine_gsea <- function(gsea, common_table, th){
	# common_table obtained from Function_common_gsea()
	# th obtained from Function_combine_th() by visual inspection (maximizing number of genes)
	# NOTE: Need to fix

	common_table   <- common_table[tanimoto<th,]
	common_merge   <- merge(common_table[,list(N=sum(n_common!=0)),by="Gene_set_1"], 
						   common_table[,list(N=sum(n_common!=0)),by="Gene_set_2"], 
						   by.x="Gene_set_1", by.y="Gene_set_2", all.x=T, all.y=T)
	common_merge$N.x[is.na(common_merge$N.x)] <- 0
	common_merge$N.y[is.na(common_merge$N.y)] <- 0
	common_merge$N <- with(common_merge, N.x + N.y)
	common_merge   <- common_merge[,c(1,4),with=F]
	setnames(common_merge, c("Gene_set", "N"))

	common_table   <- merge(common_table, common_merge, by.x="Gene_set_1", by.y="Gene_set")
	common_table   <- merge(common_table, common_merge, by.x="Gene_set_2", by.y="Gene_set")
	common_table$N <- with(common_table, N.x + N.y)

	# common_table[,N_1:=sum(tanimoto!=0),by="Gene_set_1"]
	# common_table[,N_2:=sum(tanimoto!=0),by="Gene_set_2"]
	# common_table$N <- with(common_table, N_1+N_2)
	# common_sets   <- with(common_table[tanimoto>th], unique(c(Gene_set_1, Gene_set_2))) # Sets that we should exclude
	# gsea          <- gsea[!Gene_set %in% common_sets,]
	all_genes     <- length(unique(gsea$genes))
	source        <- unique(gsea$source) #HAS TO BE SINGLE SOURCE!!

	# combine to single set those which share genes
	# common_table  <- common_table[!Gene_set_1 %in% common_sets,][!Gene_set_2 %in% common_sets,] # Exlude sets
	print(common_table)
	# th_table      <- common_table[tanimoto!=0,]
	th_table      <- common_table
	# print(th_table)	

	# Filtering sets for tanimoto>0
	remove        <- c()
	n_combined    <- 1
	new_sets      <- data.table()

	# th_table$acc  <- with(th_table, n_1+n_2-n_common)
	# th_table[,N:=length(Gene_set_2), by="Gene_set_1"]
	th_table      <- th_table[order(N)]
	print(th_table)

	i       <- 1
	covered <- c()
	while( (i!=nrow(th_table))){ # && (length(covered)!=all_genes) ){
		
		set_1 <- th_table[i,]$Gene_set_1
		set_2 <- th_table[i,]$Gene_set_2

		if ((!set_1 %in% remove) && (!set_2 %in% remove)){
			cat(paste0( round(i*100/nrow(th_table),2), "%","\n"))
			new_genes    <- unique(gsea[Gene_set %in% c(set_1, set_2),]$genes)
			set_partners <- c(th_table[tanimoto>0][Gene_set_1 %in% c(set_1, set_2),]$Gene_set_2,
							  th_table[tanimoto>0][Gene_set_2 %in% c(set_1, set_2),]$Gene_set_1)
			new_sets     <- rbind(new_sets,
				data.table(Gene_set = paste0("Set_", n_combined), genes=new_genes, n=length(new_genes))
								)

			# Update
			covered    <- c(covered, unique(new_genes))
			cat(paste0("Number of covered genes: ", length(covered), "\n"))
			remove     <- c(remove, set_1, set_2, set_partners)
			n_combined <- n_combined+1
		}
		i <- i+1
	}
	new_sets$source <- source
	print(new_sets)

	# Removing combined sets from original gsea
	# gsea    <- gsea[!Gene_set %in% remove,]
	# print(gsea)

	# Return filtered table
	return(new_sets)
	# return(rbind(new_sets, gsea))
}

Function_gsea_setcover <- function(gsea, common_table, th){
	# common_table obtained from Function_common_gsea()

	common_table <- common_table[tanimoto < th]
	gsea         <- gsea[Gene_set %in% unique(c(common_table$Gene_set_1, common_table$Gene_set_2)),]
	all_genes    <- unique(gsea$genes)
	all_sets     <- unique(gsea$Gene_set)

	print(length(all_sets))
	print(length(all_genes))
	covered      <- c() #genes
	cover        <- c()	#sets

	while (length(covered)!=length(all_genes)){
		cat(paste0("Sets leftover: ", length(all_sets), "\n"))
		subset   <- sapply(all_sets, function(x)  length(setdiff(unique(gsea[Gene_set==x]$genes), covered)))
		gs       <- all_sets[which.max(subset)] # Gene set with best coverage

		cover    <- c(cover, gs)
		covered  <- unique(c(covered, unique(gsea[Gene_set==gs]$genes)))

		# Update
		all_sets <- setdiff(all_sets, gs) # So we don't calculate it next time
		cat(paste0("Number of covered genes: ",length(covered),"\n"))
	}

	cat(paste0("Total number of gene sets needed to cover it all: ", length(cover),"\n"))
	gsea  <- gsea[Gene_set %in% cover,]
	return(gsea)
}

Function_setcover_nonoverlap <- function(gsea, common_table) {
	# Greedely find sets that are non-overlapping
	# Taylor for gene_set inputs

	source    <- unique(gsea$source) #HAS TO BE SINGLE SOURCE!!
	th_table  <- common_table[tanimoto!=0,]

	# Start with sets that cover the most
	sets      <- common_table[,list(N=sum(tanimoto!=0)),by="Gene_set_1"][order(N)]$Gene_set_1
	# sets      <- unique(gsea[,c("Gene_set", "n"),with=F])[order(-n),]$Gene_set
	print(paste0("Number of leftover sets: ", length(sets)))
	all_genes <- length(unique(gsea$genes))

	cover     <- data.table()
	covered   <- c()

	count     <- 1
	while((length(covered)!=all_genes) && (count<=length(sets))){
		print(count)
		gene_set     <- sets[count]
		print(gene_set)

		# Remove its partners from main set list
		set_partners <- unique(c(th_table[Gene_set_1==gene_set]$Gene_set_2, th_table[Gene_set_2==gene_set]$Gene_set_1))
		sets         <- setdiff(sets, set_partners)
		print(paste0("Number of leftover sets: ", length(sets)))

		# Get genes that is covering
		genes        <- gsea[Gene_set==gene_set,]$genes
		covered      <- c(covered, genes)
		print(paste0("Number of covered genes: ", length(covered)))

		# Update
		count        <- count+1
		cover        <- rbind(cover, data.table(Gene_set=gene_set, genes=genes, n=length(genes)))
	}

	cover$source <- source

	#Return
	return(cover)
}

Function_load_gsea_updated <- function(gsea){
	# Gene sets obtained from: http://software.broadinstitute.org/gsea/downloads.jsp
	# gsea: "c4", "cancer", "c4_cancer"

	if (gsea=="l1000"){
		cat("Using l1000 instead of GSEA\n")
		l1000    <- Function_load_l1000("GSEA_FILES/OTHER_MODELS/L1000_2017.csv")
		gsea_all <- data.table(Gene_set="l1000", genes=l1000) 

	} else if (gsea=="cs"){
		cat("Using Compressed sensing instead of GSEA\n")
		gsea_all <- data.table()

	} else {
		cat("Using gseas\n")
		g_set      <- strsplit(gsea, "_")[[1]]

		if (gsea=="c2"){
			cat("Curated MSigDB collection (C2) gene sets\n")
			cat("Using a tanimoto similarity threshold filter of <0.1 for overlapping gene sets\n")
			gsea_set <- fread(paste0(in_folder, "GSEA_FILES/c2_sets"))
			gsea_ths <- readRDS(paste0(in_folder, "GSEA_FILES/gene_sets_c2_common_count.rds"))
			gsea_all <- gsea_set[!Gene_set %in% with(gsea_ths[tanimoto>0.1], unique(c(Gene_set_1, Gene_set_2)))]
		} else {
			gsea_all   <- data.table()
			for (g in g_set){
				cat(paste0(g, "\n"))
				gsea_set   <- fread(paste0(in_folder, "GSEA_FILES/", g, "_sets"))
				gsea_all   <- rbind(gsea_all, gsea_set)
			}
		}
	}
	return(gsea_all)
}

Function_gene_filter <- function(train_exp, test_exp, gsea_data, gene_th){
	# Filter gene sets, train and test expression dataset according to gene_th -OR- var_feat choice

	common_genes <- intersect(rownames(train_exp), rownames(test_exp))
	train_exp    <- train_exp[common_genes,]
	test_exp     <- test_exp[common_genes,]

	if (nrow(gsea_data)>0){
		gsea_data    <- gsea_data[genes %in% common_genes,]
		gsea_data[,N:=length(unique(genes)), by="Gene_set"]
		# gsea_data    <- gsea_data[N >= gene_th,]
		gsea_data    <- gsea_data[N >= gene_th,]
		gsea_data$N  <- NULL

		common_genes <- unique(gsea_data$genes)
		train_exp    <- train_exp[common_genes,]
		test_exp     <- test_exp[common_genes,]
	}
	
	return(list(train_exp, test_exp, gsea_data))
}

Function_gene_feat <- function(in_exp, gsea_data, method="median", var_feat, n_feat, features=c(), train=F, train_exp=matrix()){
	# Median/mean/cor by method based
	# NOTE: When "cor" is used, gsea_data is ignored
	# 		However, when *_cor_gsea is used, gsea_data is needed
	require(data.table)
	require(parallel)

	exp_data  <- data.table(melt(in_exp))
	exp_data  <- merge(exp_data, gsea_data, by.x="Var1", by.y="genes", allow.cartesian=TRUE)

	if (method=="median"){
		cat("featurization method used: median\n")
		exp_data  <- exp_data[,list(m=median(value)), by=c("Gene_set", "Var2")]
	} else if (method=="mean"){
		cat("featurization method used: mean\n")
		exp_data  <- exp_data[,list(m=mean(value)), by=c("Gene_set", "Var2")]
	} else if (method=="sum"){
		cat("featurization method used: sum\nKeep in mind that row standarization is inherently applied\n")
		in_exp    <- t(scale(t(in_exp)))
		in_exp    <- (in_exp < -0.8)*-1  + (in_exp > 0.8) # Binarize gene expression values
		exp_data  <- data.table(melt(in_exp))
		exp_data  <- merge(exp_data, gsea_data, by.x="Var1", by.y="genes", allow.cartesian=TRUE)
		exp_data  <- exp_data[,list(m=sum(value)), by=c("Gene_set", "Var2")]

	} else if (method=="cor"){
		cat("featurization method used: correlation\n")
		if (train==T){
			exp_data     <- cor(in_exp, method="spearman")
		} else{
			common_genes <- intersect(rownames(in_exp), rownames(train_exp))
			exp_data     <- cbind(in_exp[common_genes,], train_exp[common_genes,])
			exp_data     <- cor(exp_data, method="spearman")
			exp_data     <- exp_data[colnames(in_exp), colnames(train_exp)]
		}
		exp_data  <- data.table(melt(exp_data)) #Var1=samples, Var2=features
		setnames(exp_data, c("Var2", "Gene_set", "m"))

	} else if (method=="mean_cor_gsea"){
		cat("featurization method used: mean_cor_gsea\n")
		if (train==T){
			exp_data  <- exp_data[,list(m=mean(value)), by=c("Gene_set", "Var2")]
			exp_data  <- acast(exp_data, Gene_set~Var2, value.var="m")
			exp_data  <- cor(exp_data, method="spearman")
		} else{
			in_data     <- exp_data[,list(m=mean(value)), by=c("Gene_set", "Var2")]
			in_data     <- acast(in_data, Gene_set~Var2, value.var="m")

			train_data  <- data.table(melt(train_exp))
			train_data  <- merge(train_data, gsea_data, by.x="Var1", by.y="genes", allow.cartesian=TRUE)
			train_data  <- train_data[,list(m=mean(value)), by=c("Gene_set", "Var2")]
			train_data  <- acast(train_data, Gene_set~Var2, value.var="m")

			common_feat <- intersect(rownames(in_data), rownames(train_data))
			exp_data    <- cbind(in_data[common_feat,], train_data[common_feat,])
			exp_data    <- cor(exp_data, method="spearman")
			exp_data    <- exp_data[colnames(in_exp), colnames(train_exp)]

		}
		exp_data  <- data.table(melt(exp_data)) #Var1=samples, Var2=features
		setnames(exp_data, c("Var2", "Gene_set", "m"))

	} else if (method=="cor_gsea"){
		cat("featurization method used: cor_gsea\n")

		# Filter gseas to 15<n_genes<=50
		exp_data # Var1(genes), Var2(samples), Gene_set
		exp_data[,N:=length(unique(Var1)), by="Gene_set"]
		gseas     <- unique(exp_data[N>50,][N<=100,]$Gene_set)
		cat(paste0("Number of gene sets used: ", length(gseas), "\n"))
		in_exp    <- in_exp[unique(exp_data[Gene_set %in% gseas,]$Var1),]

		# Setup parallelization
		nodes          <- detectCores()
		cat(paste0("Number of available nodes detected: ", nodes, "\n"))
		cl   		   <- makeCluster(nodes)
		setDefaultCluster(cl)
		clusterExport(cl, varlist=c("data.table", "as.data.table", "acast",
									"train", "gseas", "exp_data", "in_exp", "train_exp"),
					  envir=environment())
		cat("Done exporting to parallel\n")

		if (train==T){
			train_cor <- parLapply(cl, gseas, function(g){
				genes     <- exp_data[Gene_set==g]$Var1
				gsea_cor  <- cor(in_exp[genes,], method="spearman")
				gsea_cor  <- data.table(melt(gsea_cor))
				gsea_cor$Var2 <- paste0(gsea_cor$Var2, g) # To add specific cell_cor_gsea identifier
				return(gsea_cor)
			})
		} else{
			train_exp    <- train_exp[unique(exp_data[Gene_set %in% gseas,]$Var1),]
			train_cor <- parLapply(cl, gseas, function(g){
				genes     <- exp_data[Gene_set==g]$Var1
				gsea_cor  <- cor(cbind(train_exp[genes,], in_exp[genes,]), method="spearman")
				gsea_cor  <- gsea_cor[colnames(in_exp), colnames(train_exp)] # To separate test samples from features
				gsea_cor  <- data.table(melt(gsea_cor))
				gsea_cor$Var2 <- paste0(gsea_cor$Var2, g) # To add specific cell_cor_gsea identifier
				train_cor <- rbind(train_cor, gsea_cor)
			})
		}
		train_cor <- do.call(rbind, train_cor)

		# Stop parallelization
		stopCluster(cl)
		exp_data <- train_cor
		setnames(exp_data, c("Var2", "Gene_set", "m"))
		print(exp_data)
	}
	
	if (var_feat==T){
		cat("gsea features chosen by top-n variance\n")
		if (length(features)>0){
			exp_data <- exp_data[Gene_set %in% features,]
		} else{
			#NOTE: Small modifications to scale variance 0-1
			exp_data[,scale_m:=Function_range_0_1(m), by="Gene_set"]
			top_var  <- exp_data[,list(VAR=var(m),
									   SCALED_VAR=var(scale_m)),
								 by=c("Gene_set")] # Variance per gene_set across samples
			top_var  <- top_var[order(-VAR),]
			top_var$CUM_VAR_PERC        <- cumsum(top_var$VAR)/sum(top_var$VAR)
			top_var  <- top_var[order(-SCALED_VAR),]
			top_var$CUM_SCALED_VAR_PERC <- cumsum(top_var$SCALED_VAR)/sum(top_var$SCALED_VAR)

			
			top_var  <- top_var[order(-SCALED_VAR)]$Gene_set[1:n_feat] # Top n_feat chosen VAR/SCALED_VAR
			# top_var  <- top_var[CUM_SCALED_VAR_PERC <= n_feat]$Gene_set # In this scenario, n_feat is cummulative variance
																		# explained by the features ordered by variance
																		# In this scenario, this number is in the range 0-1

			exp_data <- exp_data[Gene_set %in% top_var,]
		}
	}

	exp_data  <- acast(exp_data, Gene_set~Var2, value.var = "m") #MODIFIED CAST ORDER!!!
	# print(dim(exp_data))
	return(exp_data)
}

Function_target_data <- function(drug){
	# Datasets used in paper
	# drug choices are: Bortezomib, Docetaxel, Cisplatin, Erlotinib (for the independent datasets),
	#					tcga_Cisplatin, tcga_Gemcitabine, tcga_5-Fluorouracil and tcga_Temozolomide (tcga datasets).

	if (drug=="Bortezomib"){

	}
}

Function_exp_combat <- function(exp.matrix.1, exp.matrix.2) {
  #Function to normalize batch effects across two expression matrices
  #INPUT: 2 expression matrices with columns as samples and rows as genes
  #OUPOUT: 2 expression matrices, HOWEVER, genes that are not common across expression datasets will be removed
  
  require(sva)
  require(pamr)
  require(limma)
  
  # Remove null standard gene deviations from each matrix first
  exp.matrix.1 <- exp.matrix.1[apply(exp.matrix.1, 1, sd)!=0, ]
  exp.matrix.1 <- exp.matrix.1[!is.na(rownames(exp.matrix.1)),]
  exp.matrix.2 <- exp.matrix.2[apply(exp.matrix.2, 1, sd)!=0, ]
  exp.matrix.2 <- exp.matrix.2[!is.na(rownames(exp.matrix.2)),]

  common.genes <- intersect( unique(rownames(exp.matrix.1)) , 
                             unique(rownames(exp.matrix.2)))
  cat("common genes for batch normalization\n")
  print(length(common.genes))

  col.1  <- colnames(exp.matrix.1)
  col.2  <- colnames(exp.matrix.2)
  m.1    <- exp.matrix.1[common.genes,]
  m.2    <- exp.matrix.2[common.genes,]
  m.both <- cbind(m.1, m.2)
  
  # Filter out zero sd
  m.both <- m.both[apply(m.both, 1, sd)!=0, ]
  
  # Continue
  batch <- c( rep("m1", length(col.1)) , 
              rep("m2", length(col.2)) )
  
  modcombat <- model.matrix(~1, data.frame(1:length(batch)) ) # Intercept term
  
  combat.m <- ComBat(dat=m.both, batch=batch, mod=modcombat, par.prior = T, prior.plots = F)
  
  exp.matrix.1 <- combat.m[,col.1]
  exp.matrix.2 <- combat.m[,col.2]
  
  return(list(EXP.1=exp.matrix.1 , EXP.2=exp.matrix.2))
}

Function_remove_low_var <- function(bn_exp, gee_exp=F, geeproc_exp=F, gsea="F"){
	# Removes the lowest 20% varying genes in the training sets

	cgp_exp  <- bn_exp[[1]]
	test_exp <- bn_exp[[2]]

	if (gsea!="F"){
		cat("gsea method found, not removing lowest 20%\n")

	} else if (gee_exp==T | geeproc_exp==T){
		cat("removing lowest 20% gene features by variance\n")
		var_genes <- apply(cgp_exp, 1, var)
		top_80    <- sort(var_genes, decreasing = T)[1:(length(var_genes)*0.8)]
		
		print(dim(cgp_exp))
		cgp_exp   <- cgp_exp[var_genes %in% top_80,]
		print(dim(cgp_exp))
		text_exp  <- test_exp[rownames(cgp_exp),] # Homogenize external dataset
	}
	return(list(cgp_exp, test_exp))
}

Function_wrap_normalization <- function(exp_1, exp_2, norm, target_1, tissue) {
	print(dim(exp_1))
	print(dim(exp_2))
	if (norm=="batch"){
		cat("Applying batch normalization\n")
		bn_exp      <- Function_exp_combat(exp_1, exp_2)
	} else if (norm=="gene"){
		cat("Applying gene std\n")
		bn_exp      <- Function_common_genes(
			Function_std_gene(exp_1), Function_std_gene(exp_2)
			)
	} else if(norm=="sample"){
		cat("Applying sample std\n")
		# bn_exp      <- Function_common_genes(
		# 	Function_std_sample(exp_1), Function_std_sample(exp_2)
		# 	)
		bn_exp      <- Function_common_genes(exp_1, exp_2)
		bn_exp      <- list(Function_std_sample(bn_exp[[1]]),
							Function_std_sample(bn_exp[[2]]))

	} else if(norm=="gene_sample"){
		cat("Applying gene and then sample std\n")
		bn_exp      <- Function_common_genes(
			Function_std_gene(exp_1), Function_std_gene(exp_2)
			)
		bn_exp      <- Function_common_genes(
			Function_std_sample(bn_exp[[1]]), Function_std_sample(bn_exp[[2]])
			)
	} else if(norm=="batch_gene"){
		cat("Applying batch and then gene std\n")
		bn_exp      <- Function_exp_combat(exp_1, exp_2)
		bn_exp      <- Function_common_genes(
			Function_std_gene(bn_exp[[1]]), Function_std_gene(bn_exp[[2]])
			)
	} else if(norm=="batch_sample"){
		cat("Applying batch and then sample std\n")
		bn_exp      <- Function_exp_combat(exp_1, exp_2)
		bn_exp      <- Function_common_genes(
			Function_std_sample(bn_exp[[1]]), Function_std_sample(bn_exp[[2]])
			)
	} else if(norm=="batch_gene_sample"){
		cat("Applying batch and then gene and then sample std\n")
		bn_exp      <- Function_exp_combat(exp_1, exp_2)
		bn_exp      <- Function_common_genes(
			Function_std_gene(bn_exp[[1]]), Function_std_gene(bn_exp[[2]])
			)
		bn_exp      <- Function_common_genes(
			Function_std_sample(bn_exp[[1]]), Function_std_sample(bn_exp[[2]])
			)
	} else if (norm=="tissue_gene"){
		cat("Applying tissue-specific gene normalization\n")
		exp_1       <- Function_tissue_norm(exp_1, target_1, tissue)
		bn_exp      <- Function_common_genes(
			exp_1, Function_std_gene(exp_2)
			)
	} else if (norm=="tissue_gene_sample"){
		cat("Applying tissue-specific gene-sample normalization\n")
		exp_1       <- Function_tissue_norm(exp_1, target_1, tissue)
		bn_exp      <- Function_common_genes(
			exp_1, Function_std_gene(exp_2)
			)
		bn_exp      <- Function_common_genes(
			Function_std_sample(bn_exp[[1]]), Function_std_sample(bn_exp[[2]])
			)
	} else if (norm=="gtex"){
		cat("Applying gtex binary normalization\n")
		bn_exp      <- Function_common_genes(Function_gtex_norm(exp_1),
											 Function_gtex_norm(exp_2))
	} else{
		cat("No gene normalization applied\n")
		bn_exp      <- Function_common_genes(exp_1, exp_2)
	}
	return(bn_exp)
}

Function_prep_data <- function(train_exp, test_exp, train_table, test_table, drug){
	# train_exp   : as produced by Function_rma_normalize_affy() or Function_cgp_new_exp()
	#             : alternatively, matrix with genes as rows and samples as columns
	# train_table : as produced by Funcion_cgp_new()
	# test_exp    : as produced by Function_target_data()
	#             : alternatively, matrix with genes as rows and samples as columns
	# test_table  : as produced by Function_target_data()
	# 			  : alternatively, a data.table() object with column names: "Compound", "cell_name" and "target"

	train_table  <- train_table[Compound==drug,]
	test_table   <- test_table[Compound==drug,]

	# Fixing common samples
	train_common <- intersect(colnames(train_exp), train_table$cell_name)
	test_common  <- intersect(colnames(test_exp), test_table$cell_name)
	train_table  <- train_table[cell_name %in% train_common,]
	test_table   <- test_table[cell_name %in% test_common,]

	# Filtering common genes and gsea
	gsea         <- Function_load_gsea("c4_cancer", 20)
	filtered     <- Function_gene_filter(train_exp, test_exp, gsea, 20)
	train_exp    <- filtered[[1]]
	test_exp     <- filtered[[2]]
	gsea         <- filtered[[3]]

	# Standarizing expression datasets (gene standardization)
	train_exp    <- Function_std_gene(train_exp, f=T)
	test_exp     <- Function_std_gene(test_exp, f=T)

	# Obtaining genomic fingerprints
	train_exp    <- Function_gene_feat(train_exp, gsea) #Inverting axis
	test_exp     <- Function_gene_feat(test_exp, gsea) #Inverting axis
	train_exp    <- train_exp[as.character(train_table$cell_name),]
	test_exp     <- test_exp[as.character(test_table$cell_name),]

	# Fixing commom feature space
	common_feat    <- intersect(colnames(train_exp), colnames(test_exp))
	cat(paste0(length(common_feat), " common features used\n"))
	train_exp      <- train_exp[,common_feat]
	test_exp       <- test_exp[,common_feat]

	# Return
	return(list(train_exp, test_exp, train_table, test_table))
}

Function_cgp_name <- function(name){

	if ((name=="bortezomib_a") | (name=="bortezomib_b")){
		name <- "Bortezomib"
	} else if (grepl("tcga", name)){
		name <- gsub("tcga_", "", name)
		name <- strsplit(name, "_")[[1]][1]

	} else{
		substr(name, 1, 1) <- toupper(substr(name, 1, 1))
	}
	return(name)
}

Function_load_train <- function(target, gee_exp=F, gee_target=F, geeproc_exp=F){
	drug       <- Function_cgp_name(target)
	
	loc_folder <- "" #MODIFIED
	# loc_folder <- "~/Documents/Rotation/"

	# Target related options
	if (gee_target==T){
		cat("using geeleher sensitivity definitions - pre 2014\n")
		cgp_new <- do.call(rbind, list(
			fread(paste0(loc_folder,"DATABASES/GEELEHER/paper/Data/docetaxelData/sensitivity_data_for_drug_1007.csv")),
			fread(paste0(loc_folder,"DATABASES/GEELEHER/paper/Data/bortezomibData/sensitivity_data_for_drug_104.csv")),
			fread(paste0(loc_folder,"DATABASES/GEELEHER/paper/Data/cisplatinData/sensitivity_data_for_drug_1005.csv")),
			fread(paste0(loc_folder,"DATABASES/GEELEHER/paper/Data/Erlotinib/sensitivity_data_for_drug_1.csv"))
			))
		cgp_new <- cgp_new[,c("Drug Name", "Cell Line Name", "IC 50", "Cancer Type"),with=F]
		setnames(cgp_new, c("Compound", "cell_name", "IC50", "Site")) # NOTE: Added Tissue site information

		cgp_new <- cgp_new[Compound==drug,]
		cgp_new$AUC  <- cgp_new$IC50 # For preprocessing purposes
		# cgp_new$IC50 <- exp(cgp_new$IC50) # Because Geeleher et al. processed data outputs Ln(IC50)
		
	} else{
		cat("using most recent sensitivity definitions\n")
		cgp_new <- readRDS(paste0(in_folder,"CGP_FILES/082916_cgp_new.rds"))[Compound==drug,]
		cgp_new$IC50 <- log(cgp_new$IC50) #MODIFIED!!!
	}

	# Expression related options
	if (gee_exp==T){
		# UPDATED TO REMOVE DUPLICATED SAMPLES (Both copies)!!!
		cat("Using expression data as per Geeleher et al. (processed by us)\n")
		cgp_exp <- readRDS(paste0(in_folder,"CGP_FILES/cgp_exp_gee.rds"))
		cgp_exp <- cgp_exp[!is.na(rownames(cgp_exp)),]
		rep_samples <- colnames(cgp_exp)[duplicated(colnames(cgp_exp))] # Removing duplicated (including original)
		uni_samples <- setdiff(colnames(cgp_exp), rep_samples)
		cgp_exp     <- cgp_exp[,uni_samples]

	} else if(geeproc_exp==T){
		cat("Using pre-processed expression data by Geeleher et al.\n")
		cgp_exp <- readRDS(paste0(in_folder,"CGP_FILES/cgp_exp_geeproc.rds"))
		
		# Remove samples for which we have multiple array of information (duplicated)
		# May need to include this step earlier on
		cgp_exp <- cgp_exp[!is.na(rownames(cgp_exp)),]
		cat("Removing samples\n")
		rep_samples <- colnames(cgp_exp)[duplicated(colnames(cgp_exp))] # Removing duplicated (Both copies)
		uni_samples <- setdiff(colnames(cgp_exp), rep_samples)
		cgp_exp     <- cgp_exp[,uni_samples]

	} else{
		cat("Using latest expression array\n")
		cgp_exp <- readRDS(paste0(in_folder,"CGP_FILES/083016_cgp_exp.rds"))
	}
	
	common_samples <- intersect(cgp_new$cell_name, colnames(cgp_exp))
	cgp_new <- cgp_new[cell_name %in% common_samples,]
	cgp_exp <- cgp_exp[,common_samples]

	return(list(cgp_new, cgp_exp))
}

Function_load_test <- function(target, gee_target_proc_exp=F){
	loc_folder <- "" #MODIFIED
	# loc_folder <- "~/Documents/Rotation/"

	if (target=="bortezomib_a"){
		in_table <- readRDS(paste0(in_folder,"OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
		target_table <- in_table[["feat_table_a"]][Compound=="Bortezomib"]
		exp_table    <- in_table[["exp_table_a"]]
	} else if (target=="bortezomib_b"){
		in_table <- readRDS(paste0(in_folder,"OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
		target_table <- in_table[["feat_table_b"]][Compound=="Bortezomib"]
		exp_table    <- in_table[["exp_table_b"]]
	} else if (target=="docetaxel"){
		in_table <- readRDS(paste0(in_folder,"OBJECTS/030217_GEE_DOCETAXEL.rds"))
		target_table <- in_table[["feat_table"]][Compound=="Docetaxel"]
		
		if (gee_target_proc_exp==T){
			cat("Using target expression arrays processed by Geeleher et al. 2014\n")
			load(paste0(loc_folder, "DATABASES/GEELEHER/paper/Data/docetaxelData/doce_rma_syms_brainArray.RData"))
			exp_table      <- doceVivoNorm_syms
			colnames(exp_table) <- sapply(colnames(exp_table), function(x) as.character(strsplit(x, "[.]")[[1]][1]))	
		} else{
			exp_table    <- in_table[["exp_table"]]	
		}
		
	} else if (target=="cisplatin"){
		in_table <- readRDS(paste0(in_folder,"OBJECTS/030217_GEE_CISPLATIN.rds"))
		target_table <- in_table[["feat_table"]][Compound=="Cisplatin"]

		if (gee_target_proc_exp==T){
			cat("Using target expression arrays processed by Geeleher et al. 2014\n")
			load(paste0(loc_folder, "DATABASES/GEELEHER/paper/Data/cisplatinData/cisplatinBreast.RData"))
			exp_table    <- cisVivoNorm_syms
			colnames(exp_table) <- sapply(colnames(exp_table), function(x) strsplit(x, "[.]")[[1]][1])	
		} else{
			exp_table    <- in_table[["exp_table"]]	
		}
		
	} else if (target=="erlotinib"){
		in_table <- readRDS(paste0(in_folder,"OBJECTS/030417_GEE_ERLOTINIB.rds"))
		target_table <- in_table[["feat_table"]][Compound=="Erlotinib"]

		if (gee_target_proc_exp==T){
			cat("Using target expression arrays processed by Geeleher et al. 2014\n")
			load(paste0(loc_folder,"DATABASES/GEELEHER/paper/Data/Erlotinib/erlotinibExpr.RData"))
			exp_table    <- erlotProc_syms
			colnames(exp_table) <- sapply(colnames(exp_table), function(x) strsplit(x, "[.]")[[1]][1])
			colnames(exp_table) <- sapply(colnames(exp_table), function(x) strsplit(x, "_")[[1]][1])
		} else{
			exp_table    <- in_table[["exp_table"]]
		}

	} else if (grepl("tcga_", target)){
		target_table <- fread(paste0(in_folder, "CGP_FILES/", target, "_target"))
		exp_table    <- fread(paste0(in_folder, "CGP_FILES/", target, "_exp"))
		exp_table    <- data.matrix(data.frame(exp_table, row.names = 1))
	}

	target_table$cell_name <- as.character(target_table$cell_name)

	common_samples <- intersect(target_table$cell_name, colnames(exp_table))
	target_table   <- target_table[cell_name %in% common_samples,]
	exp_table      <- exp_table[,common_samples]

	return(list(target_table, exp_table))
}

Function_ridge_predict <- function(train_exp, test_exp, train_table, test_table){
	require(glmnet)
	require(car)
	require(reshape2)
	require(methods)
	# Input is output from Function_prep_data()

	# Powertransforming
	alpha          <- powerTransform(train_table$IC50)[[6]]
	train_table$pt <- train_table$IC50 ^ alpha

	# Modeling
	ridge_model   <- cv.glmnet(as(train_exp, "dgCMatrix"),
							   train_table$pt, 
							   alpha=0, lambda=10^seq(6, -10, by = -.1), standardize=T)

	ridge_predict <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.min, 
                     		 newx = test_exp)

	# Un-Powertransform
	ridge_predict <- (abs(ridge_predict) ^ (1/alpha)) * sign(ridge_predict)

	# Output
	predictions   <- merge(test_table, data.table(ridge_predict, keep.rownames=T), by.x="cell_name", by.y="rn")
	setnames(predictions, c("cell_name", "Compound", "target", "predicted"))

	return(predictions)
}

Function_svm_gee_model_cv <- function(train_table, test_table, train_exp, test_exp, gee_target=F){
	# Applying crossvalidation for svm
	# Models require finetuning so it might take some time to obtain output
	# Creating multiple models for each cross-validation fold to predict on both heldout-fold and external dataset
	require(caret)
	require(e1071)
	cat("Predicting\n")

	offset_used <- F
	# Modifications as per compute_phenotype_function.R found in Geeleher et al. (2014)
	if (min(train_table$IC50)<0){
		offset_used      <- T # Needs offset
		offset           <- -min(train_table$IC50) + 1 # Keep in mind that this handling Ln(IC50), may not be appropriate
		train_table$IC50 <- train_table$IC50 + offset	
	}
	# Power transform as usual
	alpha          <- powerTransform(train_table$IC50)[[6]]
	train_table$pt <- train_table$IC50 ^ alpha

	if (gee_target==T){
		target_choices <- c("pt")
	} else{
		target_choices <- c("pt", "AUC")
	}

	main_table     <- data.table()
	cv_table       <- data.table()
	int_table      <- data.table()
	train_frame    <- data.frame(t(train_exp[,train_table$cell_name]))
	test_frame     <- data.frame(t(test_exp[, test_table$cell_name]))

	# Applying per fold
	set.seed(1234)
	split_indeces  <- createFolds(1:nrow(train_frame), 5, list=T)
	split_names    <- names(split_indeces)

	for (n in split_names){
		cat(paste0(n,"\n"))
		test_index  <- split_indeces[[n]]
		train_index <- setdiff(1:nrow(train_frame), test_index)

		for (i in target_choices){
			
			# Tuning svm
			svm_model     <- tune(svm, target~., data=cbind(train_frame, 
								  target=train_table[[i]])[train_index,],
								  ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
			train_predict <- predict(svm_model$best.model, newdata=train_frame[train_index,]) # Predicting on train set
			svm_cv        <- predict(svm_model$best.model, newdata=train_frame[test_index,]) # Predicting for cv-error
			svm_predict   <- predict(svm_model$best.model, newdata=test_frame) # Predicting on external set

			if (i=="pt"){
				svm_predict <- (abs(svm_predict) ^ (1/alpha)) * sign(svm_predict)
				if(offset_used==T){
					cat("using offset\n")
					svm_predict <- svm_predict - offset
				}
			}
			print( data.table(svm_predict, keep.rownames=T))
			print(test_table)

			predictions <- data.table(test_table, predicted=svm_predict, fold=n)
			internal    <- data.table(train_table[test_index,], predicted=svm_cv, fold=n)
			train_error <- data.table(train_table[train_index,], predicted=train_predict, fold=n)
			main_table  <- rbind(main_table, 
								data.table(predictions, variable=i))
			cv_table    <- rbind(cv_table, 
								data.table(internal, variable=i))
			int_table   <- rbind(int_table,
								 data.table(train_error, variable=i))
		}
	}
	return(list(external=main_table, internal=cv_table, train=int_table))
}

Function_erlo_lasso_cv <- function(train_table, test_table, train_exp, test_exp, gsea="F"){
	# Applying crossvalidation for erlo on a select subset of train cell lines
	# Cell lines are split for 15 lowest IC50 (sensitive) and 55 highest IC50 (resistant)
	# This is a classification problem where powerTransform is not applied
	# Each fold model require finetuning so it might take some time to obtain output
	# Creating multiple models for each cross-validation fold to predict on both heldout-fold and external dataset
	# NOTE: If gsea==F, the script assumes that input features are gene expression
	require(caret)
	require(ridge)
	cat("Predicting\n")

	# Split into sensitive (15 cells) and resistant (55 cells) classes
	train_table    <- train_table[order(IC50)]
	sens_cells     <- train_table$cell_name[1:15]
	res_cells      <- train_table$cell_name[(nrow(train_table)-54):nrow(train_table)]
	train_table    <- train_table[cell_name %in% c(sens_cells, res_cells),]
	train_table$B  <- ifelse(train_table$cell_name %in% sens_cells, 1, 0)
	target_choices <- c("B")
	print(dim(train_table))

	# Select top-1000 differentiable genes between the two classes if gene expression used
	if (gsea=="F"){
		cat("not using gsea, selecting top-1000 differetiably genes\n")
		# p_vals    <- apply(train_exp, 1, function(x) t.test(x[sens_cells], x[res_cells])$p.value)
		# top_genes <- sort(p_vals)[1:1000]
		p_vals    <- rowttests(train_exp[,c(sens_cells, res_cells)], 
								as.factor(c(rep(1, length(sens_cells)),rep(0,length(res_cells)))))
		p_vals    <- data.table(p_vals, gene=rownames(train_exp))[order(p.value),]
		top_genes <- p_vals$gene[1:1000]

		train_exp <- train_exp[top_genes,]
	}

	# No powerTransform assumes that IC50 values are in the log-scale (ln(IC50))
	main_table     <- data.table()
	cv_table       <- data.table()
	int_table      <- data.table()
	train_frame    <- data.frame(t(train_exp[,train_table$cell_name]))
	test_frame     <- data.frame(t(test_exp[, test_table$cell_name]))

	# Applying per fold (Reduced fold-split from 10 to 5 due to smaller sample size)
	set.seed(1234)
	split_indeces  <- createFolds(1:nrow(train_frame), 5, list=T)
	split_names    <- names(split_indeces)

	for (n in split_names){
		cat(paste0(n,"\n"))
		test_index  <- split_indeces[[n]]
		train_index <- setdiff(1:nrow(train_frame), test_index)

		for (i in target_choices){
			
			ridge_model   <- logisticRidge(target~., data=cbind(train_frame, 
										 target=train_table[[i]])[train_index,])
			train_predict <- predict(ridge_model, newdata=train_frame[train_index,]) # Predicting on train set
			ridge_cv      <- predict(ridge_model, newdata=train_frame[test_index,]) # Predicting for cv-error
			ridge_predict <- predict(ridge_model, newdata=test_frame) # Predicting on external set

			print( data.table(ridge_predict, keep.rownames=T))
			print(test_table)

			predictions <- data.table(test_table, predicted=ridge_predict, fold=n)
			internal    <- data.table(train_table[test_index,], predicted=ridge_cv, fold=n)
			train_error <- data.table(train_table[train_index,], predicted=train_predict, fold=n)
			main_table  <- rbind(main_table, 
								data.table(predictions, variable=i))
			cv_table    <- rbind(cv_table, 
								data.table(internal, variable=i))
			int_table   <- rbind(int_table,
								 data.table(train_error, variable=i))
		}
	}
	return(list(external=main_table, internal=cv_table, train=int_table))
}

Function_gene_filter_count <- function(train_exp, test_exp, gsea_data){
	# Counts genes per input gsea

	common_genes <- intersect(rownames(train_exp), rownames(test_exp))
	train_exp    <- train_exp[common_genes,]
	test_exp     <- test_exp[common_genes,]

	gsea_data[,N_TOTAL:=length(unique(genes)), by="Gene_set"]
	gsea_data    <- gsea_data[genes %in% common_genes,]
	gsea_data[,N:=length(unique(genes)), by="Gene_set"]
	
	return(list(train_exp, test_exp, gsea_data))
}

Function_prep_gsea_n <- function(feat_matrix, train_target, gsea_data){
	# Obtain gene_feat x ordered sample (by IC50) to plot heatmap of train_feat values
	# Ideally gsea (N) can sort gene_feat on the y-axis

	feat_table <- data.table(melt(feat_matrix))
	setnames(feat_table, c("Gene_set", "cell_name", "feat"))

	feat_table <- merge(feat_table, train_target, by="cell_name")
	feat_table <- merge(feat_table, 
						unique(gsea_data[,c("Gene_set", "N"),with=F]),
						by="Gene_set")
	return(feat_table)
}

Function_range_0_1 <- function(x){
	normalized = (x-min(x))/(max(x)-min(x))
	return(normalized)
}

Function_se <- function(x) {
	# Calculates the standard error of the mean
	sem<-sd(x)/sqrt(length(x))
	return(sem)
}

Function_load_l1000 <- function(file_in){
	# file_in as in L1000_2017.csv
	require(data.table)

	x <- fread(file_in, skip=2, select = c(1,3))
	x <- x[!grepl("INV-", V1)]

	cat(paste0("Returning ", length(x$V3), " L1000 genes\n"))
	return(x$V3)
}

Function_gee_baseline <- function(){
	# Produces Geeleher baseline table of AUCs
	# lowercase name compounds as per our nomenclature

	# Gene expression
	x <- data.table(Compound=c("docetaxel", "cisplatin", "bortezomib_a", "bortezomib_b"),
					AUC=c(0.79, 0.63, 0.61, 0.68),
					type="Gene expression")

	# Top-var 300 mean
	y <- data.table(Compound=c("docetaxel", "cisplatin", "bortezomib_a", "bortezomib_b"),
					AUC=c(0.76, 0.67, 0.62, 0.65),
					type="Top-var 300")

	x <- rbind(x,y)
	return(x)
}

Function_rmse <- function(a,b){
  return(sqrt(mean((a-b)^2)))
}

Function_tissue <- function(target, gee_target){
	# Obtains tissue type depending target choice
	if (target=="docetaxel" | target=="cisplatin"){
		tissue <- "breast"
	} else if (target=="erlotinib"){
		tissue <- "lung"
	} else if (grepl("bortezomib", target)){
		tissue <- ifelse(gee_target==T, "blood", "haematopoietic_and_lymphoid_tissue")
	} else {
		tissue <- ""
	}

	cat(paste0("Tissue choice: ", tissue, "\n"))
	return(tissue)
}

Function_pairwise <- function(x, type="NRMSE"){
	# Samples in columns
	n_samples <- ncol(x)
	comb      <- combn(n_samples, 2)

	if (type=="NRMSE"){
		print("NRMSE")
		max_x  <- max(x)
		min_x  <- min(x)
		x_diff <- max_x-min_x

		nrmse     <- apply(comb, 2, function(y) Function_rmse(x[,y[1]], x[,y[2]])) #RMSE
		nrmse     <- nrmse/x_diff #NRMSE
		samples_1 <- comb[1,1:ncol(comb)]
		samples_2 <- comb[2,1:ncol(comb)]

		calc      <- data.table(Var1=samples_1, Var2=samples_2, nrmse=nrmse)
	} else if (type=="euclidean"){
		print("euclidean")
		
		calc <- data.table(melt(as.matrix(dist(t(x), method="euclidean"))))
		setnames(calc, c("Var1", "Var2", "euclidean"))
	}
	return(calc)
}

Function_prep_tcga_exp <- function(x){
	# x: tcga expression object as in 090616_fireshose_all_exp.rds
	# Returns:
	#   main_table: concatenation of all matrices with common ordered genes
	#   cancer_table: table of sample per cancer
	
	# Clean up
	repeated     <- c("COADREAD", "KIPAN", "GBMLGG", "STES", "PANGI")
	for (i in repeated){
		x[[i]] <- NULL
	}

	# Process
	genes        <- lapply(names(x), function(i) rownames(x[[i]][["tumor"]]))
	genes_common <- Reduce(intersect, genes)
	cat(paste0("Common genes: ", length(genes_common), "\n"))

	main_table   <- lapply(names(x), function(i) {
		return(x[[i]][["tumor"]][genes_common,])
	})
	main_table   <- do.call(cbind, main_table)

	target_table <- lapply(names(x), function(i) {
		return(data.table(cancer = i, samples = colnames(x[[i]][["tumor"]]) ))
	})
	target_table <- do.call(rbind, target_table)

	return(list(main_table, target_table))
}

Function_tcga_firehose_clinical <- function(main_folder){
	# main_folder where all cancer data is stored
	# main_folder such as stddata__2016_07_15

	# Obtain all folders
	folders <- list.dirs(main_folder, recursive=F)

	#Unzip necessary files in each folder:
	main_table  <- data.table()
	for (f in folders){
		print(f)
		cancer   <- strsplit(f, "__")[[1]][2]
		cancer   <- strsplit(cancer, "/")[[1]][2]

		internal <- list.dirs(f, recursive=F)
		files    <- list.files(internal, full.names=T, pattern="Level_4")
		files    <- files[grep("tar.gz", files)]
		files    <- files[!grepl("md5", files)]
		
		if (length(files)>10){
			
			# Untar files to directory
			untar(files, exdir=internal)
			folder   <- list.dirs(internal, recursive = F)
			folder   <- folder[grepl("Level_4", folder)]
			file_in  <- list.files(folder, full.names=T, pattern="All_CDE")

			# Read and process file
			x        <- fread(file_in)
			y        <- data.table(t(x), keep.rownames=T)
			setnames(y, as.character(y[1,]))
			y        <- y[2:nrow(y),]
			y        <- y[,c(1,3:6),with=F]
			print(y)

		} else{
			files    <- list.files(internal, full.names=T, pattern="Clinical.Level_1")
			files    <- files[grep("tar.gz", files)]
			files    <- files[!grepl("md5", files)]
			files    <- files[!grepl("FFPE", files)]

			# Untar files to directory
			untar(files, exdir=internal)
			folder   <- list.dirs(internal, recursive = F)
			folder   <- folder[grepl("Clinical.Level_1", folder)]
			file_in  <- list.files(folder, full.names=T, pattern="clin.merged")

			# Read and process file
			x        <- fread(file_in, skip=1)
			y        <- data.table(t(x), keep.rownames=T)
			setnames(y, as.character(y[1,]))
			y        <- y[2:nrow(y),]

			# Combine information per patient for all follow-ups
			vital    <- y[,grepl("vital_status", colnames(y)),with=F]
			vital    <- apply(vital, 1, function(x) ifelse("dead" %in% x, "dead", ifelse("alive" %in% x, "alive", NA)))

			death    <- y[,grepl("days_to_death", colnames(y)),with=F]
			death    <- apply(death, 1, function(x) ifelse( sum(!is.na(x)) == 0, NA, max(  as.numeric(x[!is.na(x)]) ) ))

			last     <- y[,grepl("days_to_last_followup", colnames(y)),with=F]
			last     <- apply(last, 1, function(x) ifelse( sum(!is.na(x)) == 0, NA, max( as.numeric(x[!is.na(x)]) ) ))

			y        <- data.table(y[,c("patient.bcr_patient_barcode", "patient.patient_id"),with=F],
									vital, death, last)
			setnames(y, colnames(main_table)[1:5])
			print(y)
		}

		main_table <- rbind(main_table, 
							data.table(y, cancer=cancer))
	}

	#Remove repeated cancer samples in table
	repeated     <- c("COADREAD", "KIPAN", "GBMLGG", "STES", "PANGI", "FPPP")
	main_table   <- main_table[!cancer %in% repeated,]
	return(main_table)
}

Filter_results_target_filter <- function(x, params, target, th, th_variable, x_target, roc_table){
	# x is table with "Compound" column
	# params are column names that are considered as parameters
	# x_target are columns that you want to retain

	target_list   <- lapply(x_target, function(l) x[[l]])
	names(target_list) <- x_target

	x             <- x[,c("Compound", params),with=F]
	x$param       <- do.call(paste0, lapply(colnames(x)[2:ncol(x)], function(j) x[[j]]))
	
	for (l in x_target){
		x[[l]] <- target_list[[l]]
	}

	y         <- x[Compound==target,]
	y         <- y[y[[th_variable]] > th]
	x         <- x[param %in% y$param,]
	
	# Merge to get best ROC result from table
	roc_table <- merge(roc_table, x[,-c("mean_auc"),with=F])

	return(list(x, roc_table))
}

Function_se <- function(x) {
	# Calculates the standard error of the mean
	sem<-sd(x)/sqrt(length(x))
	return(sem)
}
