# ridge_results.R
# Process predition results from ridge_predict.R

library(data.table)
library(reshape2)
library(ROCR)

# Load arguments
args       <- commandArgs(trailingOnly = TRUE)
drugs      <- c("bortezomib_a", "bortezomib_b", "cisplatin", "docetaxel")
gsea_type  <- args[1]
train_set  <- args[2]
train_gmv  <- args[3]
test_gmv   <- args[4]

# Process
main_table <- lapply(drugs, function(x) {

	file_name  <- paste0("GSEA_FILES/RESULTS/", gsea_type, "_", train_set, "_", x, "_gee_ridge_pt_IC50_gmv_train_",
				   		 train_gmv, "_gmv_test_", test_gmv, ".rds")
	temp_table <- readRDS(file_name)

	if (grepl("bor", x)) {
		temp_table <- temp_table[,list(AUC = as.numeric(performance(prediction(-prediction, Compound), "auc")@y.values),
									   cv_error = unique(cv_error)), by="p_val_th"]
	} else{
		temp_table <- temp_table[,list(AUC = as.numeric(performance(prediction(-prediction, target), "auc")@y.values),
									   cv_error = unique(cv_error)), by="p_val_th"]
	}

	temp_table$Compound <- x
	return(temp_table)
})

main_table <- do.call(rbind, main_table)

# Return
saveRDS(main_table, paste0("GSEA_FILES/RESULTS/", gsea_type, "_", train_set, "_gee_ridge_pt_IC50_gmv_train_", 
							train_gmv, "_gmv_test_", test_gmv, ".rds"))
print("Done")