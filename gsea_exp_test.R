# gsea_exp_test.R
# Tests the efficiency of binary expression features on external datasets

library(data.table)
library(reshape2)
library(parallel)

Function_prep_for_pca <- function(internal_cast, external_cast){

	# Order features (to make sure PCA is applied correctly to external data)
	common_feat   <- intersect(colnames(internal_cast), colnames(external_cast))
	print(c("total number of features used:" , length(common_feat)))
	external_cast <- external_cast[, common_feat]
	internal_cast <- internal_cast[, common_feat] 

	# Scaling
	internal_cast <- scale(internal_cast)
	internal_mean <- attributes(internal_cast)$`scaled:center`
	internal_sd   <- attributes(internal_cast)$`scaled:scale`

	# Apply scaling to external
	external_cast <- sweep(external_cast, 2, internal_mean, "-")
	external_cast <- sweep(external_cast, 2, internal_sd, "/")

	# PCA
	internal_pca  <- prcomp(internal_cast, center = F, scale. = F)
	external_pca  <- external_cast %*% internal_pca$rotation

	return(list(internal=internal_pca$x, external=external_pca))
}

Function_target_data <- function(drug){
	# Loads external data for a particular drug (geeleher datasets) or ccle

	if(drug=="Bortezomib"){
		target_data  <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
		target_exp   <- target_data[["exp_table_a"]]
		target_data  <- target_data[["feat_table_a"]][Compound==target_drug]
		target_exp   <- target_exp[,unique(target_data$cell_name)]

	} else if (drug=="Docetaxel"){
		target_data  <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_DOCETAXEL.rds"))
		target_exp   <- target_data[["exp_table"]]
		target_data  <- target_data[["feat_table"]][Compound==target_drug]
		target_exp   <- target_exp[,unique(target_data$cell_name)]
	} else if (drug =="Cisplatin"){
		target_data  <- readRDS(paste0(in_folder, "OBJECTS/030217_GEE_CISPLATIN.rds"))
		target_exp   <- target_data[["exp_table"]]
		target_data  <- target_data[["feat_table"]][Compound==target_drug]
		target_exp   <- target_exp[,unique(target_data$cell_name)]
	} else if (drug=="Erlotinib"){
		target_data  <- readRDS(paste0(in_folder, "OBJECTS/030417_GEE_ERLOTINIB.rds"))
		target_exp   <- target_data[["exp_table"]]
		target_data  <- target_data[["feat_table"]][Compound==target_drug]
		target_exp   <- target_exp[,unique(target_data$cell_name)]
	} else if(drug=="ccle"){
		target_exp   <- readRDS(paste0(in_folder,"OBJECTS/121116_ccle_exp.rds"))
		target_data  <- readRDS(paste0(in_folder,"OBJECTS/121116_ccle.rds"))
		target_data  <- target_data[,c("cell_name", "Compound", "ActArea"),with=F]
		setnames(target_data, c("cell_name", "Compound", "target"))
	}

	return(list(target_data, target_exp))
}

# Parameters
args          <- commandArgs(trailingOnly = TRUE)
fdrs          <- c(0.01, 0.05, 0.1, 0.0001, 0.001)
# gsea_type     <- c("h", "c1", "c2.cp.biocarta", "c3", "cancer") 
gsea_type     <- c("cancer")
target_drug   <- "ccle" #"Bortezomib"
metric        <- "logIC50"
metric        <- ifelse(metric=="logIC50", metric, "")
multi_type    <- as.logical(args[1]) # Only call True when all prev fdrs have been calculated

mac_folder    <- "/Users/jzamalloa/Documents/Rotation/PIPELINES/"
lab_folder    <- "/home/zamalloa/Documents/Rotation/PIPELINES/"
gpu_folder    <- "/tigress/zamalloa/"
in_folder     <- gpu_folder

# Load common files
cgp_new      <- readRDS(paste0(in_folder, "CGP_FILES/082916_cgp_new.rds"))
target_data  <- Function_target_data(target_drug)
target_exp   <- target_data[[2]]
target_data  <- target_data[[1]]

# target_drug  <- ifelse(target_drug == "Bortezomib", "bortezo", target_drug )

# Process per fdr threshold
for (fdr_th in fdrs){
	print(c(target_drug, fdr_th))

	out_gsea_file <- paste0(in_folder, "GSEA_FILES/external_sets/", target_drug,"_a_cgp_",
							paste0(gsea_type, collapse = "_"), "_fdr_", fdr_th, ".rds")

	# Are we dealing with multiple gseas?
	if (multi_type==F){

		# Obtain binary internal features
		cgp_gsea_file   <- paste0(in_folder, "GSEA_FILES/", gsea_type, "_gsea_cgp_pvals.rds")
		cgp_gsea        <- readRDS(cgp_gsea_file)
		cgp_gsea$binary <- ifelse(cgp_gsea$p_vals < fdr_th, 1, 0)
		cgp_gsea[,var:=var(binary), by="gs"]
		cgp_gsea        <- cgp_gsea[var!=0,]

		# Do we need to calculate external features
		if (file.exists(out_gsea_file)){
			print("Using existing file")
			samples_gsea <- readRDS(out_gsea_file)

		} else {

			gsea_genes   <- fread(paste0(in_folder, "GSEA_FILES/", gsea_type, "_sets"))

			# Obtain external features if necessary
			samples         <- colnames(target_exp)
			gsea_genes      <- gsea_genes[genes %in% rownames(target_exp)]

			print(c("Variable (Non-zero variance) Number of GSEA features: ", length(unique(cgp_gsea$gs)) ))
			gsea_names      <- unique(cgp_gsea$gs)
			gsea_needed     <- matrix(sapply(gsea_names, function(s) strsplit(s, "$", fixed=T)[[1]]), ncol=2, byrow=T)

			count           <- 0
			all             <- length(gsea_names)

			# Setup parallelization
			nodes<-detectCores()
			print(nodes)
			cl<-makeCluster(nodes)
			setDefaultCluster(cl)
			clusterExport(cl, varlist=c("data.table", "as.data.table", "count", "all" ,"gsea_needed", "gsea_genes", "target_exp", "samples"),envir=environment())

			# Apply parallelization
			samples_gsea    <- parApply(cl, gsea_needed, 1, function(x){

				genes_1 <- gsea_genes[Gene_set==x[1]]$genes
				genes_2 <- gsea_genes[Gene_set==x[2]]$genes

				pvals <- sapply(samples, function(y){
					return(wilcox.test(target_exp[genes_1, y], target_exp[genes_2, y], alternative="greater")$p.value)
				})

				count <<- count+1
				print(count/all)
				return(data.table(sample = samples, gs = paste0(x[1], "$", x[2]), p_vals=pvals)) # Correct later per sample

			})

			# Close parallelization
			stopCluster(cl)
			samples_gsea    <- do.call(rbind, samples_gsea)

			# Correct P-values for FDR and binarize features
			samples_gsea    <- samples_gsea[,p_vals:=p.adjust(p_vals, method="fdr"), by="sample"]
			samples_gsea$binary <- ifelse(samples_gsea$p_vals < fdr_th, 1, 0)
			saveRDS(samples_gsea, out_gsea_file)
			print("Done writing gsea features on external dataset")
		}

	} else {
		# If multiple gsea are desired, this implies that fdr for all gseas have been calculated at the establisehd fdrs
		print("Using pre-calculated multiple gsea sets")

		# Load external
		samples_gsea    <- lapply(gsea_type, function(x) {
			return(readRDS(paste0(in_folder, "GSEA_FILES/external_sets/", target_drug,"_a_cgp_",
								  x, "_fdr_", fdr_th, ".rds")))
		})
		samples_gsea    <- do.call(rbind, samples_gsea)
		samples_gsea    <- samples_gsea[sample %in% colnames(target_exp),]

		# Load internal
		cgp_gsea        <- lapply(gsea_type, function(x) {
			return(readRDS(paste0(in_folder, "GSEA_FILES/", x, "_gsea_cgp_pvals.rds")))
		})
		cgp_gsea        <- do.call(rbind, cgp_gsea)
		cgp_gsea$binary <- ifelse(cgp_gsea$p_vals < fdr_th, 1, 0)
		cgp_gsea[,var:=var(binary), by="gs"]
		cgp_gsea        <- cgp_gsea[var!=0,]

	}

	# Build model
	# cgp_gsea     <- acast(cgp_gsea, sample~gs, value.var = "binary")
	# samples_gsea <- acast(samples_gsea, sample~gs, value.var = "binary")

	# pca_applied  <- Function_prep_for_pca(cgp_gsea, samples_gsea)

	# feat_table   <- cgp_new[Compound==target_drug]
	# if (metric =="logIC50"){
	# 	print("Using IC50")
	# 	feat_table     <- feat_table[,c("cell_name", "IC50"),with=F]
	# 	feat_table$AUC <- log(feat_table$IC50)
	# }
	
	# feat_table   <- feat_table[,c("cell_name", "AUC"),with=F]		
	# target_table <- target_data[Compound==target_drug][,c("cell_name", "target"),with=F]

	# feat_table   <- merge(feat_table, data.table(pca_applied[["internal"]], keep.rownames = T), 
	# 					  by.x="cell_name", by.y="rn")[,-c("cell_name"),with=F]
	# target_table <- merge(target_table, data.table(pca_applied[["external"]], keep.rownames = T), 
	# 					  by.x="cell_name", by.y="rn")[,-c("cell_name"),with=F]

	# # Predict
	# model_output <- data.table()
	# n_feat       <- c(2:10, 15, 20, 30, 40, 50, 60, 80, 100, 125, 150, 175, 200, 300, 400, 500, 600)
	# max_feat     <- ncol(feat_table) - 1
	# n_feat       <- n_feat[n_feat <= max_feat]

	# for (f in n_feat){
	# 	print(f)
	# 	pca_model     <- lm(AUC~., data=feat_table[, 1:(f+1), with=F])
	# 	pred_features <- setdiff(colnames(feat_table[, 1:(f+1), with=F]), "AUC")

	# 	predictions   <- predict.lm(object = pca_model, newdata = target_table[ ,pred_features, with=F])

	# 	model_output  <- rbind(model_output, data.table(n_feat = f, ACTUAL=target_table$target, PREDICTION=predictions, target=target_drug))
	# }

	# # Output
	# saveRDS(model_output, paste0(in_folder, "GSEA_FILES/cgp_", target_drug,"_",paste0(gsea_type, collapse = "_"), "_fdr_",fdr_th,"_pred",metric,".rds"))
	# print("Done writing predictions")
}

print("Done writing files")