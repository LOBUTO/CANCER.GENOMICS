# gs_encode_plot.R
# Find correlation of autoencoder features to expression features for samples
library(data.table)
library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)
library(parallel)

scriptsdir <- "/tigress/zamalloa/GIT/"
source(file.path(scriptsdir, "paper_1.R"))

# Arguments
targets       <- c("docetaxel", "cisplatin", "bortezomib_a", "bortezomib_b")
gsea     	  <- c("c4", "c4_cancer", "c4_cancer_c2.cp") #, "c2.cp")
norm     	  <- c("none" ,"sample", "gene", "batch")
gsea_m   	  <- "mean"
g_filter 	  <- seq(20, 200, 20) #c(seq(10, 40, 10), seq(50, 200, 50))
in_folder     <- "/tigress/zamalloa/"

# Load files
in_file_1     <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_logs.rds"
in_file_2     <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_pred_cor.rds"
out_plot_1    <- "GSEA_FILES/PLOTS/autoencoder_1_early_nrmse_best.pdf"
out_plot_2    <- "GSEA_FILES/PLOTS/autoencoder_1_early_nrmse_cgp_gee_cor_cor.pdf"
# best_table    <- readRDS(in_file_1)[["best_table"]]
cor_tables    <- readRDS(in_file_2)

main_train    <- data.table()
main_test     <- data.table()
cor_cor_train <- data.table()
cor_cor_test  <- data.table()
for (d in targets){

	train_files   <- Function_load_train(d, gee_exp=F, gee_target=T, geeproc_exp=T)
	cgp_old_exp   <- train_files[[2]]
	clin_exp      <- Function_load_test(d, gee_target_proc_exp=T)[[2]]

	for (nr in norm){
		# Normalizing expression
		all_exp    <- Function_wrap_normalization(cgp_old_exp, clin_exp, nr)
		print("Normalized")

		# Calculate correlation across samples with full expression
		train_cor_exp <- cor(all_exp[[1]], method="spearman")
		test_cor_exp  <- cor(all_exp[[2]], method="spearman")

		train_cor_exp <- data.table(melt(train_cor_exp))[Var1!=Var2]
		test_cor_exp  <- data.table(melt(test_cor_exp))[Var1!=Var2]
		for (g in gsea){
			for (f in g_filter){
				train_cor <- cor_tables$train_cor[Compound==d & gsea==g & norm==nr & g_filter==f,][Var1!=Var2,]
				if (nrow(train_cor)){
					print(c(d, g, nr, f))
					train_cor  <- merge(train_cor, train_cor_exp, by=c("Var1", "Var2"))	

					test_cor   <- cor_tables$test_cor[Compound==d & gsea==g & norm==nr & g_filter==f,][Var1!=Var2,]
					test_cor   <- merge(test_cor, test_cor_exp, by=c("Var1", "Var2"))

					main_train <- rbind(main_train, train_cor)
					main_test  <- rbind(main_test, test_cor)

					train_cor  <- train_cor[,list(COR=cor(value.x, value.y, method="spearman")), 
											 by=c("Compound", "gsea", "norm", "g_filter", "n_feat")]
					test_cor   <- test_cor[,list(COR=cor(value.x, value.y, method="spearman")), 
											 by=c("Compound", "gsea", "norm", "g_filter", "n_feat")]
					cor_cor_train <- rbind(cor_cor_train, train_cor)
					cor_cor_test  <- rbind(cor_cor_test, test_cor)
				}
			}
		}
	}
}

# Store
print("Storing files")
out_file_1    <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_cgp_gee_cor.rds"
out_file_2    <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_cgp_gee_cor_cor.rds"
saveRDS(list(main_train=main_train, main_test=main_test), out_file_1)
saveRDS(list(main_train=cor_cor_train, main_test=cor_cor_test), out_file_2)

# in_table      <- readRDS("GSEA_FILES/RESULTS/autoencoder_1_single_nrmse_cgp_gee_cor_cor.rds")
# cor_cor_train <- in_table[["main_train"]]
# cor_cor_test  <- in_table[["main_test"]]

print("Plotting")
# pdf(out_plot_1, width=12, height=8)
# for (d in targets){
# 	print(
# 		ggplot(best_table[Compound==d], aes(gsea, valid, fill=factor(g_filter))) +
# 			geom_bar(stat="identity", position="dodge") + facet_wrap(~norm, scales="free") +
# 			scale_fill_viridis(discrete = T, direction = -1) +
# 			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 			ggtitle(d)
# 		)
# }
# print(
# 	ggplot(best_table, aes(g_filter, gsea, fill=valid)) +
# 	facet_wrap(Compound~norm, scales="free") +
# 	geom_raster() + scale_fill_viridis(discrete = F, direction = -1)
# 	)
# dev.off()

# pdf(out_plot_2, width=12, height=8)
# for(d in targets){
# 	print(
# 		ggplot(rbind(data.table(cor_cor_train[Compound==d], source="Train"), 
# 					 data.table(cor_cor_test[Compound==d], source="Test")), 
# 			aes(gsea, COR, fill=factor(g_filter))) +
# 			geom_bar(stat="identity", position="dodge") + facet_grid(source~norm, scales="free") +
# 			scale_fill_viridis(discrete = T, direction = -1) +
# 			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 			ggtitle(d)
# 		)
# }

# print(
# 	ggplot(cor_cor_train, aes(factor(paste0(g_filter)), gsea, fill=COR)) +
# 		facet_grid(Compound~norm, scales="free") +
# 		geom_raster() + scale_fill_viridis(discrete = F, direction = -1)+
# 		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 		ggtitle("Train correlation of Autoencoded features to gene expression")
# 	)
# print(
# 	ggplot(cor_cor_test, aes(factor(g_filter), gsea, fill=COR)) +
# 		facet_grid(Compound~norm, scales="free") +
# 		geom_raster() + scale_fill_viridis(discrete = F, direction = -1) +
# 		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# 		ggtitle("Test correlation of Autoencoded features to gene expression")
# 		)
# dev.off()

cat("Done\n")