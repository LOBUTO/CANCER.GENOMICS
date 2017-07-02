# gsea_pick_bin.R
# Based on CGP data only (for the time being), takes pre-processed pvals obtained from gsea_exp.py and 
# mutation data from Rotation/DATABASES/CANCERRXGENE/mutations.txt to extract binarized GSEA features
# that are significantly associated with mutations

library(data.table)
library(reshape2)
library(ggplot2)
library(caret)

# Functions
Function_get_vertex_names <- function(v_name) {
	first  <- strsplit(v_name, ".", fixed=T)[[1]][1]
	second <- strsplit(first, "_")[[1]][1]

	return(second)
}

Function_tanimoto <- function(v1, v2, binary=T){
  # Computes the tanimoto coefficient between two vectors
  # Use names==T for binary vectors
  
  if (binary==T){
    v_intersect <- sum((v1+v2)==2)
  } else {
    v_intersect <- length(intersect(v1, v2))
  }
  
  tanimoto    <- v_intersect / (length(v1) + length(v2) - v_intersect)
  
  return(tanimoto)
}

Function_gsea_gene_sim <- function(gsea_table){
	# Based on tanimoto similarity across gene sets

	gene_sets  <- unique(gsea_table$Gene_set)

	print(length(gene_sets))
	count      <- 1
	main_table <- data.table()
	for (i in gene_sets){
		
		for (j in gene_sets){
			gs_1 <- gsea_table[Gene_set==i,]$genes
			gs_2 <- gsea_table[Gene_set==j,]$genes
			main_table <- rbind(main_table, 
								data.table(gs_1=i, gs_2=j,
											tanimoto =Function_tanimoto(gs_1, gs_2, binary=F)
										   )
								)
		}
		count <- count+1
	}
	return(main_table)
}

# Load arguments
args      <- commandArgs(trailingOnly = TRUE)
pval_th   <- as.numeric(args[1])
gsea_type <- args[2]
mut_th    <- as.numeric(args[3]) # Soft mutation threshold
hyper_th  <- as.numeric(args[4]) # Soft fdr-phyper() threshold

# Load files
cgp_mut <- fread("/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/mutations.txt", select=c(1,4,7))
gsea    <- fread(paste0("~/Documents/Rotation/PIPELINES/GSEA_FILES/",gsea_type ,"_gsea_cgp_both_T_pvals"))

# Process gsea
gsea[,p_vals:=p.adjust(pvals, method="fdr"), by="sample"]
gsea$binary <- as.numeric(gsea$p_vals < pval_th)
gsea[,var:=var(binary), by="gs"]
gsea        <- gsea[var!=0,]
all_edges   <- data.table(t(sapply(unique(gsea$gs), function(x) strsplit(x, "$", fixed=T)[[1]] )))

# Remove redundant gene sets (gene sets representing same bits across samples)
gsea_sum    <- gsea[,list(s=sum(binary)), by="gs"]
gsea_sum    <- gsea_sum[s > 20 & s < 900]
gsea        <- gsea[gs %in% gsea_sum$gs,]

# Obtain distance across gene sets based on tanimoto similarity (of genes within sets)
cancer_gsea      <- fread(paste0("~/Documents/Rotation/PIPELINES/GSEA_FILES/",gsea_type,"_sets"))
cancer_gsea_sim  <- Function_gsea_gene_sim(cancer_gsea)
cancer_gsea_sim  <- merge(cancer_gsea_sim, all_edges, by.x=c("gs_1", "gs_2"),by.y=c("V1", "V2"))

ggplot(cancer_gsea_sim, aes(gs_1, gs_2, fill=tanimoto)) + geom_raster() + theme_bw() + scale_fill_gradient(low="white", high="red") #SAVE PLOT/DONE

gs_1_dist        <- as.dist(1 - acast(cancer_gsea_sim, 
							gs_1~gs_2, value.var = "tanimoto")) # 1 - tanimoto to get distance
gs_1_clust       <- hclust(gs_1_dist)[["order"]]
gs_1_order		 <- rownames(acast(cancer_gsea_sim, gs_1~gs_2, value.var = "tanimoto"))[gs_1_clust]

gs_2_dist        <- as.dist(1 - acast(cancer_gsea_sim, 
							gs_2~gs_1, value.var = "tanimoto")) # 1 - tanimoto to get distance
gs_2_clust       <- hclust(gs_2_dist)[["order"]]
gs_2_order		 <- rownames(acast(cancer_gsea_sim, gs_2~gs_1, value.var = "tanimoto"))[gs_2_clust]

#  Found that ~150 gs can describe an individual
sample_1    <- gsea[sample=="BT-483"]

# Obtain sample map using binary scores
sample_1_all_edges        <- data.table(t(sapply(sample_1$gs, function(x) strsplit(x, "$", fixed=T)[[1]] )))
sample_1_all_edges$binary <- sample_1$binary

ggplot(sample_1_all_edges, aes(V1, V2, fill=binary)) + geom_raster() + theme_bw() + 
	scale_fill_gradient(low="white", high="black") +
	theme(axis.text  = element_text(size=5, colour="grey20"))

# Order binary scores based on gs gene similarity
sample_1_all_edges_ordered <- data.table(melt(acast(sample_1_all_edges, V1~V2, value.var = "binary")[gs_1_order, gs_2_order]))
sample_1_all_edges_ordered$value[is.na(sample_1_all_edges_ordered$value)] <- 0 #CHECK!!! P-val should be theoretically insignificant, so binary is 0

ggplot(sample_1_all_edges_ordered, aes(Var1, Var2, fill=value)) + geom_raster() + theme_bw() + 
	scale_fill_gradient(low="white", high="black") +
	theme(axis.text  = element_text(size=5, colour="grey20"))


# Process mutation data
cgp_mut$mutation <- paste0(cgp_mut$Gene, "_", cgp_mut$AA)
cgp_mut[,n:=length(unique(SAMPLE)), by="mutation"]
cgp_mut <- cgp_mut[n>=mut_th]

# Combine
gsea    <- gsea[,c("sample", "gs", "binary"),with=F] #Check that all gs have non-zero values across samples. Checked.
setkey(gsea)
gsea    <- unique(gsea)
cgp_mut <- unique(cgp_mut[,c("SAMPLE", "mutation"),with=F])

cgp_mut <- merge(gsea, cgp_mut, by.x="sample", by.y="SAMPLE", allow.cartesian=T)

# Calculate phyper
all_samples <- length(unique(cgp_mut$sample))
cgp_mut[,n_samples_mut:=length(unique(sample)), by="mutation"]
cgp_mut[,n_samples_gs:= length(unique(sample[which(binary==1)])) , by="gs"]

print("Calculating P-values")

# Hypothesis 1: white balls is you have the mutation, black is you don't have the mutation
# if you grab the "gs" group, how likely is it that you will have the mutation
# NOTE : If q==0/1 (number of samples with gs_j==1 at mutation_i), then by default p-val will not be significant
#	     So we filter those cases out to reduce number of calculations and not inflate multiple hypothesis correction
cgp_mut[,q_count:=sum(binary), by=c("mutation", "gs")]
cgp_mut  <- cgp_mut[q_count>2,]
cgp_pval <- cgp_mut[,list(p_hyper = phyper(q = sum(binary) - 1,
					  				       m = unique(n_samples_mut),
									       n = all_samples - unique(n_samples_mut),
									       k = unique(n_samples_gs),
									       lower.tail=F)), 
		 			by=c("mutation", "gs")]
cgp_pval[,p_fdr:=p.adjust(p_hyper, method="fdr"), by="mutation"]
cgp_pval$binary <- as.numeric(cgp_pval$p_fdr < hyper_th)

# Hypothesis 2: white balls is you the gs, black is you don't have the gs
# if you grab the "mutation" group, how likely is it that you will have the gs
cgp_gs_pval <- cgp_mut[,list(p_hyper = phyper(q = sum(binary) - 1,
											  m = unique(n_samples_gs),
											  n = all_samples - unique(n_samples_gs),
											  k = unique(n_samples_mut),
											  lower.tail = F)),
					    by=c("mutation", "gs")]
cgp_gs_pval$p_fdr <- p.adjust(cgp_gs_pval$p_hyper, method="fdr")
cgp_gs_pval[,p_fdr:=p.adjust(p_hyper, method="fdr"), by="gs"]
cgp_gs_pval$binary <- as.numeric(cgp_gs_pval$p_fdr < hyper_th)

# Obtained classified samples and score their gs across significantly validated mutations
cgp_pval[,n_gs:=sum(binary), by="mutation"]
cgp_pval   <- cgp_pval[n_gs>0,]
cgp_scores <- merge(cgp_mut[,c("sample", "gs", "mutation", "binary"),with=F],
					with(cgp_pval, data.table(mutation=mutation, gs=gs, score=binary)),
					by=c("mutation", "gs"))

cgp_final  <- cgp_scores[,list(precision = sum(binary[which(score==1)] == score[which(score==1)]) / (sum(binary[which(score==1)] == score[which(score==1)]) + ) ,
							   recall    = sum(binary[which(score==1)] == score[which(score==1)]) / sum(score) ),
						  by = c("sample", "mutation")]

cgp_final  <- cgp_scores[,list(precision = posPredValue(factor(binary), factor(score), positive="1"),
							   recall    = sensitivity(factor(binary), factor(score), positive="1")),
						  by = c("sample", "mutation")]

# Plot some visualizations
ccg <- fread("~/Documents/Rotation/DATABASES/CANCER_DATA/COSMIC/ccg_062417.tsv", sep="\t", select=c(1, 6, 8, 13, 14))
ccg <- ccg[Somatic=="yes",]

# Plot 1
cgp_pval_acast <- acast(cgp_pval[p_fdr<hyper_th,], mutation~gs, value.var="binary", fill=0)
mut_order      <- hclust(dist(cgp_pval_acast))[["order"]]
cgp_pval_acast <- cgp_pval_acast[mut_order,]
cgp_pval_acast <- data.table(melt(cgp_pval_acast))
cgp_pval_acast$gene <- sapply(as.character(cgp_pval_acast$Var1), function(x) strsplit(x, "_")[[1]][1])
cgp_pval_acast$ccg  <- cgp_pval_acast$gene %in% ccg[["Gene Symbol"]]

ggplot(cgp_pval_acast, aes(Var2, Var1, fill=value)) + geom_raster() +
  theme_bw() + 
  theme(axis.text.y = element_text(colour=ifelse(cgp_pval_acast$ccg==T, "red", "black"))) +
  theme(axis.text=element_text(size=6))

# Plot 2
cgp_gs_pval_acast <- acast(cgp_gs_pval[p_fdr<hyper_th,], mutation~gs, value.var="binary", fill=0)
mut_order         <- hclust(dist(cgp_gs_pval_acast))[["order"]]
cgp_gs_pval_acast <- cgp_gs_pval_acast[mut_order,]
cgp_gs_pval_acast <- data.table(melt(cgp_gs_pval_acast))
cgp_gs_pval_acast$gene <- sapply(as.character(cgp_gs_pval_acast$Var1), function(x) strsplit(x, "_")[[1]][1])
cgp_gs_pval_acast$ccg  <- cgp_gs_pval_acast$gene %in% ccg[["Gene Symbol"]]

ggplot(cgp_gs_pval_acast, aes(Var2, Var1, fill=value)) + geom_raster() +
  theme_bw() + 
  theme(axis.text.y = element_text(colour=ifelse(cgp_gs_pval_acast$ccg==T, "red", "black"))) +
  theme(axis.text=element_text(size=6))
