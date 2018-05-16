# gs_encode_results.R
library(data.table)
library(reshape2)
library(parallel)

# Load files
targets  <- c("docetaxel", "cisplatin", "bortezomib_a", "bortezomib_b")
gsea     <- c("c4", "c4_cancer", "c4_cancer_c2.cp") #, "c2.cp")
norm     <- c("none" , "sample", "gene", "batch")
gsea_m   <- "mean"
g_filter <- seq(20, 200, 20) #c(seq(10, 40, 10), seq(50, 200, 50)) #c(50, 100, 200, 500) #

best_table <- data.table()
all_table  <- data.table()

train_main_cor <- data.table()
test_main_cor  <- data.table()
for (d in targets){
	for (g in gsea){
		for (nr in norm){

			# Setup parallelization
			nodes          <- detectCores()
			cat(paste0("Number of available nodes detected: ", nodes, "\n"))
			cl   		   <- makeCluster(nodes)
			setDefaultCluster(cl)
			clusterExport(cl, varlist=c("data.table", "as.data.table", "acast",
										"g_filter", "d", "g", "nr", "gsea_m"),
						  envir=environment())
			cat("Done exporting to parallel\n")

			f_table <- parLapply(cl, g_fitler, function(f) {
				log_file    <- paste0("GSEA_FILES/RESULTS/",d,"_",g,"_",nr,"_",gsea_m,"_",f,"_autoencoder_1_early_nrmse_log")
				train_feat  <- paste0("GSEA_FILES/RESULTS/",d,"_",g,"_",nr,"_",gsea_m,"_",f,"_autoencoder_1_early_nrmse_train.txt")
				test_feat   <- paste0("GSEA_FILES/RESULTS/",d,"_",g,"_",nr,"_",gsea_m,"_",f,"_autoencoder_1_early_nrmse_test.txt")

				in_table    <- fread(log_file)

				if  (nrow(in_table)>0 & file.exists(train_feat)){ #(file.exists(train_feat)){
					print(c(d, g, nr, f))
					# print(nrow(in_table))
					best_index <- which.min(in_table$valid_error)
					best_train <- in_table$train_error[best_index]
					best_valid <- in_table$valid_error[best_index]

					best_table <- rbind(best_table,
										data.table(Compound=d, train=best_train, valid=best_valid, gsea=g, norm=nr, g_filter=f))
					all_table  <- rbind(all_table,
										data.table(in_table, gsea=g, norm=nr, g_filter=f))

					# Store feature correlation
					train_cor  <- data.frame(fread(train_feat, header=T), row.names=1)
					train_feat <- ncol(train_cor)

					if (train_feat>5){
						# print(train_cor[1:5,1:5])

						# Only store correlations if we have at least 5 features
						train_cor  <- data.table(melt(cor(t(train_cor), method="spearman")))
						print(train_cor)
						test_cor   <- data.frame(read.csv(test_feat, header=T, sep="\t"), row.names=1)
						test_cor   <- data.table(melt(cor(t(test_cor), method="spearman")))

						train_main_cor <- rbind(train_main_cor, data.table(train_cor, Compound=d, gsea=g, norm=nr, g_filter=f, n_feat=train_feat))
						test_main_cor  <- rbind(test_main_cor, data.table(test_cor, Compound=d, gsea=g, norm=nr, g_filter=f, n_feat=train_feat))
					}
				}
				return(list(b_table, a_table, tr_table, te_table))
			})
			stopCluster(cl)
		}
	}
}

out_file_1 <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_logs.rds"
out_file_2 <- "GSEA_FILES/RESULTS/autoencoder_1_early_nrmse_pred_cor.rds"

saveRDS(list(all_table=all_table, best_table=best_table), out_file_1)
saveRDS(list(train_cor=train_main_cor, test_cor=test_main_cor), out_file_2)

cat("Done\n")