# cgp_pred_ridge.R
library(data.table)
library(glmnet)
library(reshape2)
library(car)
library(methods)

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
  
  combat.m <- ComBat(dat=m.both, batch=batch, mod=modcombat, par.prior = T, prior.plots = T)
  
  exp.matrix.1 <- combat.m[,col.1]
  exp.matrix.2 <- combat.m[,col.2]
  
  return(list(EXP.1=exp.matrix.1 , EXP.2=exp.matrix.2))
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

# Load all expression files
cgp_exp    <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/121216_cgp_exp.rds")
gee_bor    <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds")
gee_doc    <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030217_GEE_DOCETAXEL.rds")
gee_erlo   <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030417_GEE_ERLOTINIB.rds")
gee_cis    <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030217_GEE_CISPLATIN.rds")

# Load target file
cgp_new    <- readRDS("~/Documents/Rotation/PIPELINES/CGP_FILES/082916_cgp_new.rds")
cgp_new    <- cgp_new[cell_name %in% colnames(cgp_exp),]

target_drugs <- c("bortezomib_a", "bortezomib_b", "docetaxel", "cisplatin")

pred_table <- data.table()
for (d in target_drugs){
	cat(paste0(d, "\n"))

	target_table  <- cgp_new[Compound==Function_cgp_name(d),]
	train_samples <- target_table$cell_name

	# Apply power transform on target
	alpha             <- powerTransform(target_table$IC50)[[6]]
	target_table$IC50 <- target_table$IC50 ^ alpha

	#Load train
	if (d=="bortezomib_a"){
		test_exp   <- gee_bor[["exp_table_a"]]
		test_table <- gee_bor[["feat_table_a"]][Compound=="Bortezomib",]
	} else if(d=="bortezomib_b"){
		test_exp   <- gee_bor[["exp_table_b"]]
		test_table <- gee_bor[["feat_table_b"]][Compound=="Bortezomib",]
	} else if(d=="docetaxel"){
		test_exp   <- gee_doc[["exp_table"]]
		test_table <- gee_doc[["feat_table"]][Compound=="Docetaxel",]
	} else if(d=="cisplatin"){
		test_exp   <- gee_cis[["exp_table"]]
		test_table <- gee_cis[["feat_table"]][Compound=="Cisplatin",]
	}
	test_samples  <- test_table$cell_name

	common_genes  <- intersect(rownames(cgp_exp), rownames(test_exp))
	train_exp     <- cgp_exp[common_genes, train_samples]
	test_exp      <- test_exp[common_genes, test_samples]

	# Model
	lambdas       <- 10^seq(3, -2, by = -.1)

	ridge_model   <- cv.glmnet( t(train_exp), target_table$IC50, 
								alpha=0, lambda=lambdas, standardize=T)
	ridge_predict_min <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.min,
							 newx = t(test_exp))
	ridge_predict_1se <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.1se,
							 newx = t(test_exp))

	# Combat model
	combat_exp    <- Function_exp_combat(train_exp, test_exp)

	combat_model  <- cv.glmnet( t(combat_exp$EXP.1), target_table$IC50, 
								alpha=0, lambda=lambdas, standardize=T)
	combat_predict_min <- predict(combat_model$glmnet.fit, s = combat_model$lambda.min,
							 newx = t(combat_exp$EXP.2))
	combat_predict_1se <- predict(combat_model$glmnet.fit, s = combat_model$lambda.1se,
							 newx = t(combat_exp$EXP.2))

	# gm model
	train_mean    <- apply(train_exp, 1, mean)
	test_mean     <- apply(test_exp, 1, mean)
	train_var     <- apply(train_exp, 1, var)
	test_var      <- apply(test_exp, 1, var)
	train_exp     <- train_exp - train_mean
	test_exp      <- test_exp - test_mean

	gm_ridge_model   <- cv.glmnet( t(train_exp), target_table$IC50, 
									alpha=0, lambda=lambdas, standardize=T)
	gm_ridge_predict_min <- predict(gm_ridge_model$glmnet.fit, s = gm_ridge_model$lambda.min,
							 	 	newx = t(test_exp))
	gm_ridge_predict_1se <- predict(gm_ridge_model$glmnet.fit, s = gm_ridge_model$lambda.1se,
							 	 	newx = t(test_exp))

	# gmv model
	train_exp     <- train_exp/train_var
	test_exp      <- test_exp/test_var

	gmv_ridge_model   <- cv.glmnet( t(train_exp), target_table$IC50, 
									alpha=0, lambda=lambdas, standardize=T)
	gmv_ridge_predict_min <- predict(gmv_ridge_model$glmnet.fit, s = gmv_ridge_model$lambda.min,
							 	 	newx = t(test_exp))
	gmv_ridge_predict_1se <- predict(gmv_ridge_model$glmnet.fit, s = gmv_ridge_model$lambda.1se,
							 	 	newx = t(test_exp))

	# Unpower transform
	ridge_predict_min     <- (abs(ridge_predict_min) ^ (1/alpha)) * sign(ridge_predict_min)
	ridge_predict_1se     <- (abs(ridge_predict_1se) ^ (1/alpha)) * sign(ridge_predict_1se)
	combat_predict_min    <- (abs(combat_predict_min) ^ (1/alpha)) * sign(combat_predict_min)
	combat_predict_1se    <- (abs(combat_predict_1se) ^ (1/alpha)) * sign(combat_predict_1se)
	gm_ridge_predict_min  <- (abs(gm_ridge_predict_min) ^ (1/alpha)) * sign(gm_ridge_predict_min)
	gm_ridge_predict_1se  <- (abs(gm_ridge_predict_1se) ^ (1/alpha)) * sign(gm_ridge_predict_1se)
	gmv_ridge_predict_min <- (abs(gmv_ridge_predict_min) ^ (1/alpha)) * sign(gmv_ridge_predict_min)
	gmv_ridge_predict_1se <- (abs(gmv_ridge_predict_1se) ^ (1/alpha)) * sign(gmv_ridge_predict_1se)

	# Return
	predictions <- data.table(sample = rep(test_samples, 8),
							  prediction = c(ridge_predict_min, ridge_predict_1se,
							  				 combat_predict_min, combat_predict_1se,
							  				 gm_ridge_predict_min, gm_ridge_predict_1se,
							  				 gmv_ridge_predict_min, gmv_ridge_predict_1se),
							  type       = c(rep("none", length(test_samples)*2),
							  			   rep("combat", length(test_samples)*2),
							  			   rep("gm", length(test_samples)*2),
							  			   rep("gmv", length(test_samples)*2)),
							  lambda     = c(rep("min", length(test_samples)),
							  			   rep("1se", length(test_samples)),
							  			   rep("min", length(test_samples)),
							  			   rep("1se", length(test_samples)),
							  			   rep("min", length(test_samples)),
							  			   rep("1se", length(test_samples)),
							  			   rep("min", length(test_samples)),
							  			   rep("1se", length(test_samples)))
							  )
	print(predictions[is.na(prediction)])

	predictions <- merge(predictions, test_table, 
						 by.x="sample", by.y="cell_name")
	
	predictions$Compound <- d

	pred_table  <- rbind(pred_table, predictions)
}

# Store
saveRDS(pred_table, "/Users/jzamalloa/Documents/Rotation/PIPELINES/CGP_FILES/081417_cgp_gee_ridge_pred.rds")
cat("Done")