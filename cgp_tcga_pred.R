# cgp_tcga_pred.R
# Combines all cgp-based TCGA predictions for plotting
# Data obtained from sources such as cgp_mut_tcga_model.R and cgp_exp_tcga_model.R, along with ridge regression on cluster

library(data.table)
library(reshape2)
library(ROCR)
library(PRROC)
library(ggplot2)
library(cowplot)
library(viridis)

# Load files
mut_log2   <- readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_mutation_auprc_log2.rds")
exp_log2   <- readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_expression_auprc_log2.rds")
main_log2  <- rbind(mut_log2, exp_log2)

mut_curve  <- readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_mutation_auprc_curve.rds")
exp_curve  <- readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/tcga_expression_auprc_curve.rds")
main_curve <- rbind(mut_curve, exp_curve)

# Load additional files from gsea
gsea_ridge <- lapply(c("5-Fluorouracil", "Cisplatin", "Gemcitabine", "Temozolomide"), function(x) {
	file_in <- paste0("~/Documents/Rotation/PIPELINES/PREDICTIONS/c4_cancer_cgp_tcga_",x,"_gee_ridge_pt_IC50_gmv_train_F_gmv_test_F_kstest.rds")
	return(readRDS(file_in))
})
gsea_ridge <- do.call(rbind, gsea_ridge)
gsea_ridge <- gsea_ridge[p_val_th==0.05,]

gsea_ridge_auprc <- gsea_ridge[,list(auprc=pr.curve(-prediction[which(target==1)], -prediction[which(target==0)], curve=T)$auc.integral, 
							   baseline=sum(target)/length(target)),
						  by=c("Compound")]
gsea_ridge_auprc$log2_auprc <- with(gsea_ridge_auprc, log2(auprc/baseline))

gsea_ridge_curve <- gsea_ridge[,list(pr=performance(prediction(-prediction, target), "prec", "rec")@y.values[[1]],
                               rec=performance(prediction(-prediction, target), "prec", "rec")@x.values[[1]],
                               baseline=sum(target)/length(target)), 
						  by=c("Compound")]

main_log2  <- rbind(main_log2, data.table(gsea_ridge_auprc, type="GSEA ridge"))
main_curve <- rbind(main_curve, data.table(gsea_ridge_curve, type="GSEA ridge"))
main_log2$type <- factor(main_log2$type, levels=c("presence", "frequency", "raw expression", "Bn expression", "GSEA ridge"))

# Plot
main_curve[,baseline:=mean(baseline), by="Compound"]
print(main_log2)
pdf("~/Documents/FOLDER/LAB/GM/073017/predictions/tcga_predictions.pdf", width=12, height=8)
print(
	ggplot(main_log2, aes(Compound, log2_auprc, fill=type)) + 
		geom_bar(stat="identity", position="dodge") +
		scale_fill_viridis(option="D", discrete = T) +
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

main_log2 <- main_log2[Compound %in% c("5-Fluorouracil", "Temozolomide")][,c("Compound", "log2_auprc", "type"),with=F]
main_log2 <- rbind(main_log2, data.table(Compound=c("5-Fluorouracil", "Temozolomide"),
										 log2_auprc = c(0.3805289, 0.7862021),
										 type= "GSEA CNN"))

curve_update <- rbind(readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/log_5f.rds"),
					  readRDS("~/Documents/Rotation/PIPELINES/PREDICTIONS/log_temo.rds"))
print(curve_update)
curve_update <- curve_update[,list(pr=performance(prediction(-predicted_pt, actual), "prec", "rec")@y.values[[1]],
                               rec=performance(prediction(-predicted_pt, actual), "prec", "rec")@x.values[[1]],
                               baseline=sum(actual)/length(actual)), 
						  by=c("Compound")]
curve_update$type <- "GSEA CNN"
main_curve <- rbind(main_curve[Compound %in% c("5-Fluorouracil", "Temozolomide")],
					curve_update)
main_curve[,baseline:=mean(baseline), by="Compound"]

pdf("~/Documents/FOLDER/LAB/GM/073017/predictions/tcga_predictions_updated.pdf", width=12, height=8)
print(
	ggplot(main_log2, aes(Compound, log2_auprc, fill=type)) + 
		geom_bar(stat="identity", position="dodge") +
		scale_fill_viridis(option="D", discrete = T) +
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