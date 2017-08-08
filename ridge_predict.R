# ridge_predict.R
library(data.table)
library(reshape2)
library(parallel)
library(ROCR)

# Functions
Function_lasso_predict <- function(test_file, test_target, train_table, train_target, pval_th, sens=15, res=55){
  require(glmnet)
  require(car)
  
  target        <- fread(test_file)
  # target[,p_vals:=p.adjust(pvals, method="fdr"), by="sample"]
  target[,p_vals:=p.adjust(pvals), by="sample"]
  target$binary <- ifelse(target$p_vals < pval_th, 1, 0)
  
  # train_table[,p_vals:=p.adjust(pvals, method = "fdr"), by="sample"]
  train_table[,p_vals:=p.adjust(pvals), by="sample"]
  train_table$binary   <- ifelse(train_table$p_vals < pval_th, 1, 0)
  
  common_gs      <- intersect(unique(target$gs), unique(train_table$gs))
  target         <- target[gs %in% common_gs,]
  train_table    <- train_table[gs %in% common_gs,]
  
  train_target   <- train_target[cell_name %in% unique(train_table$sample),]
  target_samples <- unique(target$sample)
  
  # Separate sensitive and resistant
  sens_cells     <- train_target[order(IC50),]$cell_name[1:sens]
  res_cells      <- train_target[order(-IC50),]$cell_name[1:res]
  train_target   <- train_target[cell_name %in% c(sens_cells, res_cells),]
  train_target$class <- ifelse(train_target$cell_name %in% sens_cells, 1, 0)

  # Model
  lambdas       <- 10^seq(3, -2, by = -.1)
  ridge_model   <- cv.glmnet(acast(train_table, sample~gs, value.var="binary", fill=0)[train_target$cell_name, common_gs],
  							 type.measure = "auc",
                             as.factor(train_target$class), alpha=1, lambda = lambdas, family="binomial", nfolds=5, standardize=F)
  
  # Predict
  ridge_predict <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.min, 
                         newx = acast(target, sample~gs, value.var="binary", fill=0)[target_samples, common_gs],
                         type="response")
  
  # Output predictions and model
  predictions   <- data.table(sample = target_samples, prediction = ridge_predict)
  predictions   <- merge(predictions, test_target, by.x="sample", by.y="cell_name")
  setnames(predictions, c("sample", "prediction", "Compound", "target"))
  
  return(list(predictions=predictions, ridge_model=ridge_model)) 
}

Function_lasso_predict_wrap <- function(test_file, test_target, train_table, train_target, sens=c(), res=c()){
  
  p_s   <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  nodes <- detectCores()
  nodes <- ifelse(length(p_s) <= nodes, length(p_s), nodes)
  print(nodes)
  cl    <- makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("data.table", "as.data.table", "fread", "acast", "setnames",
  							  "Function_lasso_predict", "p_s",
  							  "test_file", "test_target", "train_table", "train_target", "sens", "res"), 
  				envir=environment())

  pred_table <- parLapply(cl, p_s, function(p){
  	print(p)

  	temp_table <- data.table()

  	for (s in sens){
  		for (r in res){
  			temp_ridge  <- Function_lasso_predict(test_file, test_target, 
                                         train_table, train_target,
                                         p, sens=s, res=r)
  			temp_table  <- rbind(temp_table, 
  								 data.table(temp_ridge$predictions, s, r, cv_acu = max(temp_ridge$ridge_model$cvm)))

  		}
  	}
  	return(data.table(temp_table, p_val_th = p))
  })

  pred_table <- do.call(rbind, pred_table)
  stopCluster(cl)
  
  return(pred_table)
}

Function_ridge_predict <- function(target, test_target, train_table, train_target, pval_th, powertrans=F, response="IC50"){
  require(glmnet)
  require(car)
  
  # target[,p_vals:=p.adjust(pvals, method="fdr"), by="sample"]
  target[,p_vals:=p.adjust(pvals, method="fdr"), by="sample"]
  target$binary <- ifelse(target$p_vals < pval_th, 1, 0)
  
  # train_table[,p_vals:=p.adjust(pvals, method = "fdr"), by="sample"]
  train_table[,p_vals:=p.adjust(pvals, method="fdr"), by="sample"]
  train_table$binary   <- ifelse(train_table$p_vals < pval_th, 1, 0)
  
  # Remove features with zero variance
  train_table   <- train_table[,var=var(binary), by="gs"]
  train_table   <- train_table[var!=0,]
  train_table$var <- NULL

  common_gs      <- intersect(unique(target$gs), unique(train_table$gs))
  target         <- target[gs %in% common_gs,]
  train_table    <- train_table[gs %in% common_gs,]
  
  train_target   <- train_target[cell_name %in% unique(train_table$sample),]
  target_samples <- unique(target$sample)
  
  # Use IC50 or AUC
  if (response=="AUC"){
  	train_target$IC50 <- train_target$AUC
  }

  # Powertransform?
  if (powertrans==T){
    alpha             <- powerTransform(train_target$IC50)[[6]]
    train_target$IC50 <- train_target$IC50 ^ alpha
  } else {
    train_target$IC50 <- log(train_target$IC50)
  }
  
  # Model
  lambdas       <- 10^seq(3, -2, by = -.1)
  ridge_model   <- cv.glmnet(acast(train_table, sample~gs, value.var="binary", fill=0)[train_target$cell_name, common_gs],
                             train_target$IC50, alpha=0, lambda = lambdas, standardize=F)
  
  # Predict
  ridge_predict <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.min, 
                         newx = acast(target, sample~gs, value.var="binary", fill=0)[target_samples, common_gs])
  
  if (powertrans==T){
    ridge_predict <- ridge_predict ^ (1/alpha)
  }
  
  # Output predictions and model
  predictions   <- data.table(sample = target_samples, prediction = ridge_predict)
  predictions   <- merge(predictions, test_target, by.x="sample", by.y="cell_name")
  setnames(predictions, c("sample", "prediction", "Compound", "target"))
  
  return(list(predictions=predictions, ridge_model=ridge_model)) 
}

Function_ridge_predict_wrap <- function(test_table, test_target, train_table, train_target, powertrans=F, response="IC50"){
  
  p_s   <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  nodes <- detectCores()
  nodes <- ifelse(length(p_s) <= nodes, length(p_s), nodes)
  print(nodes)
  cl    <- makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("data.table", "as.data.table", "fread", "acast", "setnames",
  							  "Function_ridge_predict", "p_s",
  							  "test_table", "test_target", "train_table", "train_target", "powertrans", "response"), 
  				envir=environment())

  pred_table <- parLapply(cl, p_s, function(p){
  	print(p)
  	temp_ridge  <- Function_ridge_predict(test_table, test_target, 
                                         train_table, train_target,
                                         p, powertrans, response)

  	return(data.table(temp_ridge$predictions, p_val_th = p, cv_error = min(temp_ridge$ridge_model$cvm)))

  })
  pred_table <- do.call(rbind, pred_table)
  stopCluster(cl)
  
  return(pred_table)
}

Function_gmv_fill <- function(name, fill="T"){
	if (fill=="T"){
		name <- paste0(name, "_gmv")
	}
	return(name)
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

Function_target_table <- function(target){
	if (target %in% c("tcga_5-Fluorouracil", "tcga_5-Fluorouracil_gmv","tcga_5-Fluorouracil_site_gmv")){
		target_table <- fread("CGP_FILES/tcga_5-Fluorouracil_target")		
	} else if (target=="bortezomib_a"){
		target_table <- readRDS("OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds")[["feat_table_a"]][Compound=="Bortezomib"]
	} else if (target=="bortezomib_b"){
		target_table <- readRDS("OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds")[["feat_table_b"]][Compound=="Bortezomib"]
	} else if (target=="docetaxel"){
		target_table <- readRDS("OBJECTS/030217_GEE_DOCETAXEL.rds")[["feat_table"]][Compound=="Docetaxel"]
	} else if (target=="cisplatin"){
		target_table <- readRDS("OBJECTS/030217_GEE_CISPLATIN.rds")[["feat_table"]][Compound=="Cisplatin"]
	} else if (target=="erlotinib"){
		target_table <- readRDS("OBJECTS/030417_GEE_ERLOTINIB.rds")[["feat_table"]][Compound=="Erlotinib"]
	}
	return(target_table)
}

Function_load_train <- function(exp_train, train_table){
	# Loading multiple exp_train if present
	exp_train <- strsplit(exp_train, "*", fixed = T)[[1]]

	count     <- 1

	main_pval  <- data.table()
	main_table <- data.table()

	for (e in exp_train){
		print(paste0("Used ", count, " train dataset(s)"))

		train_file   <- paste0("GSEA_FILES/", gsea_type,"_gsea_", e,"_both_T_ext_gmv_", gmv_train, "_pvals") # May need to modify this
		train_pvals  <- fread(train_file)
		
		train_pvals$sample    <- paste0(train_pvals$sample, "_", letters[count])
		train_table$cell_name <- paste0(train_table$cell_name, "_", letters[count])

		main_pval    <- rbind(main_pval,  train_pvals)
		main_table   <- rbind(main_table, train_table)

		count <- count + 1
	}

	return(list(main_pval, main_table))
}

# Load arguments
args         <- commandArgs(trailingOnly = TRUE)
exp_train    <- args[1] #cgp/ cgp_site_gmv
gmv_train    <- args[2] #T/F gmv normalized train expression
gmv_test     <- args[3] #T/F gmv normalized test expression
gsea_type    <- args[4] # cancer, c4..
target       <- args[5] #Done bortezomib_a, docetaxel, cisplatin, tcga_Cisplatin, tcga_Gemcitabine_gmv, tcga_Temozolomide_gmv

# train_file   <- paste0("GSEA_FILES/", gsea_type,"_gsea_", exp_train,"_both_T_ext_gmv_", gmv_train, "_pvals")
train_table  <- readRDS("CGP_FILES/082916_cgp_new.rds")[Compound==Function_cgp_name(target),]

target_file  <- paste0("GSEA_FILES/", gsea_type,"_gsea_", target ,"_both_T_ext_gmv_", gmv_test, "_pvals")
target_table <- Function_target_table(target)

response     <- "IC50"
# out_file     <- paste0("GSEA_FILES/RESULTS/",target,"_gmv_gee_lasso_pt_",response,"_3.rds")
out_file     <- paste0("GSEA_FILES/RESULTS/", gsea_type, "_", exp_train ,"_",target,"_gee_ridge_pt_",response,
						"_gmv_train_", gmv_train, "_gmv_test_", gmv_test , ".rds")
print(exp_train)
print(out_file)

# Execute
train        <- Function_load_train(exp_train, train_table) #exp_train separator -> "*"
main_table   <- Function_ridge_predict_wrap(fread(target_file), 
										    target_table, 
											train[[1]], 
											train[[2]],
											powertrans = T,
											response=response)

# main_table   <- Function_lasso_predict_wrap(target_file, 
# 										    target_table, 
# 											fread(train_file), 
# 											train_table,
# 											sens=seq(51,75,2),
# 											res=seq(50, 65,2))

# Return
saveRDS(main_table, out_file)
print("Done")