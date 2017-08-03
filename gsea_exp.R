# gsea_exp.R
# Obtains gene expression binary scores based on pre-processed gene-set modules
# NOTE: This can be for either wilcoxon or ttest
# NOTE: No p.adjust is done!

library(data.table)
library(reshape2)
library(parallel)

# Functions
Function_exp_combat <- function(exp.matrix.1, exp.matrix.2) {
  #Function to normalize batch effects across two expression matrices
  #INPUT: 2 expression matrices with columns as samples and rows as genes
  #OUPOUT: 2 expression matrices, HOWEVER, genes that are not common across expression datasets will be removed
  
  require(sva)
  require(pamr)
  require(limma)
  
  # Remove null standard gene deviations from each matrix first
  exp.matrix.1 <- exp.matrix.1[apply(exp.matrix.1, 1, sd)!=0, ]
  exp.matrix.2 <- exp.matrix.2[apply(exp.matrix.2, 1, sd)!=0, ]

  common.genes <- intersect( unique(rownames(exp.matrix.1)) , 
                             unique(rownames(exp.matrix.2)))
  print(length(common.genes))

  col.1 <- colnames(exp.matrix.1)
  col.2 <- colnames(exp.matrix.2)
  
  m.1 <- exp.matrix.1[common.genes,]
  m.2 <- exp.matrix.2[common.genes,]
  m.both <- cbind(m.1, m.2)
  
  # Filter out zero sd
  m.both <- m.both[apply(m.both, 1, sd)!=0, ]
  
  # Continue
  
  batch <- c( rep("m1", length(col.1)) , 
              rep("m2", length(col.2)) )
  
  modcombat <- model.matrix(~1, data.frame(1:length(batch)) )
  
  combat.m <- ComBat(dat=m.both, batch=batch, mod=modcombat, par.prior = T, prior.plots = F)
  
  exp.matrix.1 <- combat.m[,col.1]
  exp.matrix.2 <- combat.m[,col.2]
  
  return(list(EXP.1=exp.matrix.1 , EXP.2=exp.matrix.2))
}

Function_process_ttest <- function(){
	nodes<-detectCores()
	print(nodes)
	cl<-makeCluster(nodes)
	setDefaultCluster(cl)
	clusterExport(cl, varlist=c("data.table", "as.data.table", "exp_samples", "gene_sets", "in_exp", "gsea"),envir=environment())

	main_table <- parApply(cl, gene_sets, 1, function(g) {

		name    <- paste0(g[1], "$", g[2])
		genes_1 <- gsea[Gene_set==g[1],]$genes
		genes_2 <- gsea[Gene_set==g[2],]$genes

		pvals <- sapply(exp_samples, function(s){

			genes_1_exp <- as.vector(in_exp[genes_1, s])
			genes_2_exp <- as.vector(in_exp[genes_2, s])

			return(t.test(genes_1_exp, genes_2_exp, alternative="greater")$p.value)
		})

		return(data.table(gs = name, sample = exp_samples, pvals))
	})

	stopCluster(cl)

	main_table <- do.call(rbind, main_table)
	return(main_table)
}

Function_process_wilcoxontest <- function(){
	nodes<-detectCores()
	print(nodes)
	cl<-makeCluster(nodes)
	setDefaultCluster(cl)
	clusterExport(cl, varlist=c("data.table", "as.data.table", "exp_samples", "gene_sets", "in_exp", "gsea"),envir=environment())

	main_table <- parApply(cl, gene_sets, 1, function(g) {

		name    <- paste0(g[1], "$", g[2])
		genes_1 <- gsea[Gene_set==g[1],]$genes
		genes_2 <- gsea[Gene_set==g[2],]$genes

		pvals <- sapply(exp_samples, function(s){

			genes_1_exp <- as.vector(in_exp[genes_1, s])
			genes_2_exp <- as.vector(in_exp[genes_2, s])

			return(wilcox.test(genes_1_exp, genes_2_exp, alternative="greater")$p.value)
		})

		return(data.table(gs = name, sample = exp_samples, pvals))
	})

	stopCluster(cl)

	main_table <- do.call(rbind, main_table)
	return(main_table)
}

# Load arguments
args      <- commandArgs(trailingOnly = TRUE)
exp_type  <- args[1]
base      <- args[2]
gsea_type <- args[3]
both      <- args[4]
test      <- args[5] # wilcoxon or ttest
bn        <- args[6] # against which other dataset is being normalize to (or None)

# Load files
in_folder <- "/tigress/zamalloa/"
if (exp_type=="cgp"){
	in_exp <- readRDS(paste0(in_folder,"CGP_FILES/083016_cgp_exp.rds"))
} else if (exp_type=="bortezomib_a"){
	in_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))[["exp_table_a"]]
} else if (exp_type=="bortezomib_b"){
	in_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))[["exp_table_b"]]
} else if (exp_type=="docetaxel"){
	in_exp <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_DOCETAXEL.rds"))[["exp_table"]]
} else if (exp_type=="erlotinib"){
	in_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_ERLOTINIB.rds"))[["exp_table"]]
} else if (exp_type=="cisplatin"){
	in_exp <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_CISPLATIN.rds"))[["exp_table"]]
}

# Do we need to apply batch normalization
if (bn!="None"){
	if (bn=="cgp"){
		bn_exp <- readRDS(paste0(in_folder,"CGP_FILES/083016_cgp_exp.rds"))
	} else if (bn=="bortezomib_a"){
		bn_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))[["exp_table_a"]]
	} else if (bn=="bortezomib_b"){
		bn_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))[["exp_table_b"]]
	} else if (bn=="docetaxel"){
		bn_exp <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_DOCETAXEL.rds"))[["exp_table"]]
	} else if (bn=="erlotinib"){
		bn_exp <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_ERLOTINIB.rds"))[["exp_table"]]
	} else if (bn=="cisplatin"){
		bn_exp <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_CISPLATIN.rds"))[["exp_table"]]
	}

	print("Batch normalizing")
	in_exp <- Function_exp_combat(in_exp, bn_exp)[["EXP.1"]]
}


# Load gseas
gsea    <- lapply(strsplit(gsea_type, "_")[[1]], function(x) {
	return(fread(paste0("GSEA_FILES/",x, "_sets")))
})
gsea    <- do.call(rbind, gsea)

# Pre-process
exp_samples <- colnames(in_exp)
exp_genes   <- rownames(in_exp)
print(nrow(gsea))
gsea        <- gsea[genes %in% exp_genes,]
print(nrow(gsea))

# Do we filter by base
if (base=="None"){
	gsea[,N:=length(unique(genes)), by="Gene_set"]
	gene_sets <- unique(gsea[N >= 100,]$Gene_set)
	gene_sets <- t(combn(gene_sets, 2))

	if (both == "T"){
		print("both one-tailed used")
		gene_sets <- rbind(gene_sets, gene_sets[,c(2,1)])
	}
	
} else {
	base      <- fread(paste0(in_folder, "GSEA_FILES/", gsea_type, "_gsea_", base, "_both_", both, "_pvals"))
	gene_sets <- unique(base$gs)
	gene_sets <- t(sapply(gene_sets, function(x) strsplit(x, "$", fixed = T)[[1]]))
}

print(c(gsea_type, nrow(gene_sets)))

# Process depending on test
if (test=="ttest"){
	print("performing ttest")
	main_table <- Function_process_ttest()
} else if (test=="wilcoxon"){
	print("performing wilcoxon")
	main_table <- Function_process_wilcoxontest()
}

print("Done calculating")

# saveRDS(main_table, "/tigress/zamalloa/GSEA_FILES/all_gsea_cgp_pvals.rds")
#saveRDS(main_table, paste0("/tigress/zamalloa/GSEA_FILES/",g,"_gsea_cgp_pvals.rds"))
file_out <- paste0(in_folder, "GSEA_FILES/", gsea_type, "_gsea_", exp_type, "_", test,"_both_", both, "_bn_",bn,"_pvals")
write.table(main_table, file_out, quote = F, sep="\t", row.names = F, col.names = T)
print("Done writing")