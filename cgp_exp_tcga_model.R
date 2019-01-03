# cgp_exp_tcga_model.R
library(data.table)
library(reshape2)
library(glmnet)
library(PRROC)
library(ROCR)
library(car)
library(ggplot2)
library(cowplot)
library(viridis)

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
  
  return(list(EXP_1=exp.matrix.1 , EXP_2=exp.matrix.2))
}


# Load files
cat("Loading files\n")
tcga_exp   <- fread("~/Documents/Rotation/PIPELINES/CGP_FILES/tcga_exp")
cgp_new    <- readRDS("~/Documents/Rotation/PIPELINES/CGP_FILES/082916_cgp_new.rds")
cgp_exp    <- fread("~/Documents/Rotation/PIPELINES/CGP_FILES/cgp_exp")

drugs      <- c("5-Fluorouracil", "Cisplatin", "Gemcitabine", "Temozolomide")

# Prep expression files
tcga_exp   <- as.matrix(data.frame(tcga_exp, row.names = "rn"))
cgp_exp    <- as.matrix(data.frame(cgp_exp, row.names = "rn"))
common     <- intersect(rownames(tcga_exp), rownames(cgp_exp))
tcga_exp   <- tcga_exp[common,]
cgp_exp    <- cgp_exp[common,]

# Prep bacth normalization of expression data
exp_norm   <- Function_exp_combat(cgp_exp, tcga_exp)

main_table <- data.table()
# Process
for (d in drugs){
	cat(paste0(d,"\n"))

	# Setup data training data per compound
	train_target      <- cgp_new[Compound==d,]
	train_cells       <- intersect(train_target$cell_name, colnames(cgp_exp))

	train_target      <- train_target[cell_name %in% train_cells,]
	train_feat        <- t(cgp_exp[,train_target$cell_name])
	train_bn_feat     <- t(exp_norm$EXP_1[,train_target$cell_name])

	# Setup testing data per compound
    test_target       <- fread(paste0("~/Documents/Rotation/PIPELINES/CGP_FILES/tcga_",d,"_target"))
    test_target       <- test_target[Compound==d,]
    tcga_samples      <- intersect(test_target$cell_name, colnames(tcga_exp))

    test_target       <- test_target[cell_name %in% tcga_samples,]
    test_feat         <- t(tcga_exp[,test_target$cell_name])
    test_bn_feat      <- t(exp_norm$EXP_2[,test_target$cell_name])

    # Model
    alpha             <- powerTransform(train_target$IC50)[[6]]
    train_target$IC50 <- train_target$IC50 ^ alpha
    
    lambdas           <- 10^seq(3, -2, by = -.1)
    ridge_model       <- cv.glmnet(train_feat, train_target$IC50,
    								alpha=0, lambda = lambdas, standardize=T)
    ridge_bn_model    <- cv.glmnet(train_bn_feat, train_target$IC50, 
    								alpha=0, lambda = lambdas, standardize=T)

    # Predict and unpower transform prediction
    ridge_predict     <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.1se, 
    							 newx = test_feat)
    ridge_bn_predict  <- predict(ridge_bn_model$glmnet.fit, s = ridge_model$lambda.1se, 
    							 newx = test_bn_feat)
    
    ridge_predict     <- (abs(ridge_predict) ^ (1/alpha)) * sign(ridge_predict)
    ridge_bn_predict  <- (abs(ridge_bn_predict) ^ (1/alpha)) * sign(ridge_bn_predict)

    # Output predictions
    predictions       <- rbind(data.table(sample=test_target$cell_name, prediction = ridge_predict, type="raw expression"),
    						   data.table(sample=test_target$cell_name, prediction = ridge_bn_predict, type="Bn expression"))
    predictions       <- merge(predictions, test_target, 
    						   by.x="sample", by.y="cell_name")
    setnames(predictions, c("sample", "prediction", "type", "Compound", "target"))
    print(predictions)

    print(predictions)
    main_table        <- rbind(main_table, predictions)
}

# Prep for plotting
main_auprc <- main_table[,list(auprc=pr.curve(-prediction[which(target==1)], -prediction[which(target==0)], curve=T)$auc.integral, 
							   baseline=sum(target)/length(target)),
						  by=c("Compound", "type")]
print(main_auprc)
main_auprc$log2_auprc <- with(main_auprc, log2(auprc/baseline))

main_curve <- main_table[,list(pr=performance(prediction(-prediction, target), "prec", "rec")@y.values[[1]],
                               rec=performance(prediction(-prediction, target), "prec", "rec")@x.values[[1]],
                               baseline=sum(target)/length(target)), 
						  by=c("Compound", "type")]
print(main_curve)

# Store predictions
saveRDS(main_auprc, "~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_expression_auprc_log2.rds")
saveRDS(main_curve, "~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_expression_auprc_curve.rds")

# Plot
pdf("~/Documents/FOLDER/LAB/GM/073017/cgp_exp/tcga_predictions.pdf", width=12, height=8)
print(
	ggplot(main_auprc, aes(Compound, log2_auprc, fill=type)) + 
		geom_bar(stat="identity", position="dodge") +
		scale_fill_viridis(option="A", discrete = T) +
		xlab("Compound") + ylab("log2(AUPRC/Baseline)")
	)
print(
	ggplot(main_curve, aes(rec, pr, colour=type)) +
		geom_line() +
		geom_hline(aes(yintercept=baseline), linetype="dashed", colour="black") +
		facet_wrap(~Compound) +
		xlab("Recall") + ylab("Precision") +
		scale_color_brewer(palette = "Set1")
	)
dev.off()

