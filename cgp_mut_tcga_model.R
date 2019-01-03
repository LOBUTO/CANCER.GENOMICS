# cgp_mut_tcga_model.R
# Predictive ridge regression model from cgp applied on tcga data
# Model uses mutation data only

library(data.table)
library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)
library(ROCR)
library(PRROC)
library(glmnet)
library(car)

# Load files
cgp_mut   <- fread("/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/mutations.txt")
cgp_new   <- readRDS("~/Documents/Rotation/PIPELINES/CGP_FILES/082916_cgp_new.rds")
tcga_mut  <- readRDS("~/Documents/Rotation/PIPELINES/OBJECTS/020617_tcga_maf.rds")
drugs     <- c("5-Fluorouracil", "Cisplatin", "Gemcitabine", "Temozolomide")

# Process
cgp_cells <- intersect(unique(cgp_new$cell_name), unique(cgp_mut$SAMPLE))
cgp_new   <- cgp_new[cell_name %in% cgp_cells,]
cgp_mut   <- cgp_mut[SAMPLE %in% cgp_cells,] 

cgp_mut   <- cgp_mut[,c("SAMPLE", "Gene", "AA"),with=F]
setkey(cgp_mut)
cgp_mut   <- unique(cgp_mut)

tcga_mut  <- tcga_mut[Variant_Classification!="Silent"]
tcga_mut  <- tcga_mut[,c("Hugo_Symbol", "sample"),with=F]
setkey(tcga_mut)
tcga_mut  <- unique(tcga_mut)

cgp_mut[,mut_count:=length(unique(SAMPLE)), by="Gene"]
cgp_mut <- cgp_mut[mut_count>=4]
cgp_mut[,sample_gene_mut:=length(AA), by=c("SAMPLE", "Gene")]
cgp_mut[,sample_n_mut:=length(AA), by="SAMPLE"]
cgp_mut$count <- with(cgp_mut, sample_gene_mut/sample_n_mut)

tcga_mut[,sample_gene_mut:=length(Hugo_Symbol), by=c("sample", "Hugo_Symbol")]
tcga_mut[,sample_n_mut:=length(Hugo_Symbol), by="sample"]
tcga_mut$count <- with(tcga_mut, sample_gene_mut/sample_n_mut)

cgp_mut  <- cgp_mut[,c("SAMPLE", "Gene", "count"),with=F]
setkey(cgp_mut)
cgp_mut  <- unique(cgp_mut)
setkey(tcga_mut)
tcga_mut <- unique(tcga_mut) 

# Add single count to model for presence only
cgp_mut$bin_count  <- 1
tcga_mut$bin_count <- 1


# Model per compound
main_table <- data.table()
for (d in drugs){
	print(d)

	# Setup data training data per compound
	train_target      <- cgp_new[Compound==d,]
	cgp_samples       <- intersect(unique(train_target$cell_name), unique(cgp_mut$SAMPLE))
	
	train_table       <- cgp_mut[SAMPLE %in% cgp_samples,]
	train_target      <- train_target[cell_name %in% cgp_samples,]

	# Setup testing data per compound
    test_target       <- fread(paste0("~/Documents/Rotation/PIPELINES/CGP_FILES/tcga_",d,"_target"))
    test_target       <- test_target[Compound==d,]
    tcga_samples      <- intersect(test_target$cell_name, unique(tcga_mut$sample))

    test_target       <- test_target[cell_name %in% tcga_samples,]
    test_table        <- tcga_mut[sample %in% tcga_samples,]

    # Common mutations
    common_features   <- intersect(unique(test_table$Hugo_Symbol), unique(train_table$Gene))
    train_table       <- train_table[Gene %in% common_features,]
    test_table        <- test_table[Hugo_Symbol %in% common_features,]

    tcga_samples      <- intersect(test_target$cell_name, unique(tcga_mut$sample))
    test_target       <- test_target[cell_name %in% tcga_samples,]
    test_table        <- tcga_mut[sample %in% tcga_samples,]

    # Cast bin tables to test
    train_bin_table   <- acast(train_table, SAMPLE~Gene, value.var="count", fill=0)[,common_features]
    test_bin_table    <- acast(test_table, sample~Hugo_Symbol, value.var = "count", fill=0)[,common_features]

    # Cast frequency tables to test
    train_table       <- acast(train_table, SAMPLE~Gene, value.var="count", fill=0)[,common_features]
    test_table        <- acast(test_table, sample~Hugo_Symbol, value.var = "count", fill=0)[,common_features]
    test_samples      <- test_target$cell_name

    print(length(common_features))
    print(dim(train_table))
    print(dim(test_table))

    # Model
    alpha             <- powerTransform(train_target$IC50)[[6]]
    train_target$IC50 <- train_target$IC50 ^ alpha
    
    lambdas           <- 10^seq(3, -2, by = -.1)
    ridge_model       <- cv.glmnet(train_table[train_target$cell_name,], train_target$IC50, 
    								alpha=0, lambda = lambdas, standardize=T)
    ridge_bin_model   <- cv.glmnet(train_bin_table[train_target$cell_name,], train_target$IC50, 
    								alpha=0, lambda = lambdas, standardize=F)

    # Predict and unpower transform prediction
    ridge_predict     <- predict(ridge_model$glmnet.fit, s = ridge_model$lambda.1se, 
    							 newx = test_table[test_samples,])
    ridge_bin_predict <- predict(ridge_bin_model$glmnet.fit, s = ridge_model$lambda.1se, 
    							 newx = test_bin_table[test_samples,])
    
    ridge_predict     <- (abs(ridge_predict) ^ (1/alpha)) * sign(ridge_predict)
    ridge_bin_predict <- (abs(ridge_bin_predict) ^ (1/alpha)) * sign(ridge_bin_predict)

    # Output predictions
    predictions       <- rbind(data.table(sample=test_samples, prediction = ridge_predict, type="frequency"),
    						   data.table(sample=test_samples, prediction = ridge_bin_predict, type="presence"))
    predictions       <- merge(predictions, test_target, 
    						   by.x="sample", by.y="cell_name")
    setnames(predictions, c("sample", "prediction", "type", "Compound", "target"))

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
saveRDS(main_auprc, "~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_mutation_auprc_log2.rds")
saveRDS(main_curve, "~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_mutation_auprc_curve.rds")

# Plot
pdf("~/Documents/FOLDER/LAB/GM/073017/cgp_mut/tcga_predictions.pdf", width=12, height=8)
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

