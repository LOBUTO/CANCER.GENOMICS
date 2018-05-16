# gsea_encode.R
# Preps train and tests datasets for gs_econde.py

scriptsdir <- "/tigress/zamalloa/GIT/"
source(file.path(scriptsdir, "paper_1.R"))

# Load arguments
args         <- commandArgs(trailingOnly = TRUE)
target       <- args[1] #Done bortezomib_a, docetaxel, cisplatin, tcga_Cisplatin, tcga_Gemcitabine_gmv, tcga_Temozolomide_gmv
gsea         <- args[2] # "F" (not logical) or gsea choices ("c4", "c4_cancer"...)
norm         <- args[3] #batch, gene, sample
gsea_m       <- args[4] #mean/median
in_folder    <- "~/Documents/Rotation/PIPELINES/"
in_folder    <- "/tigress/zamalloa/" # For cluster
# out_folder   <- "~/Documents/FOLDER/LAB/GM/123017/RIDGE_EXP/"
gee_exp      <- F # cgp GEO used by Geeleher et al. 2014 processed by us
geeproc_exp  <- T # c(T,F) # cgp GEO used by Geeleher et al. 2014 pre-processed by them
gee_target   <- T # c(T,F) # target IC50 processed by Geeleher et al. 2014
gee_target_proc_exp <- c(T) # target GEO expression processed by Geeleher et al. 2014

# Execute
for (i in gee_target){
	for(j in gee_target_proc_exp){
		for(k in geeproc_exp){
			print(c(j,k,i))

			# Load files
			train_files <- Function_load_train(target, gee_exp=gee_exp, gee_target=i, geeproc_exp=k)
			test_files  <- Function_load_test(target, gee_target_proc_exp=j)

			# Choose normalization step
			tissue      <- Function_tissue(target, gee_target)
			bn_exp      <- Function_wrap_normalization(train_files[[2]], test_files[[2]], norm, 
													   train_files[[1]], tissue)

			# Define gsea and filter genes
			gsea_set     <- Function_load_gsea_updated(gsea)

			filtered     <- Function_gene_filter(bn_exp[[1]], bn_exp[[2]], gsea_set, 5)
			train_feat   <- filtered[[1]]
			test_feat    <- filtered[[2]]
			gsea_set     <- filtered[[3]]

			# Store
			cat("writing\n")
			out_file    <- paste0(in_folder, "GSEA_FILES/TRAIN_MATRICES/", target,
						"_gee_exp", gee_exp,"_geeprocexp",k,"_geetarget",i,"_geetargetprocexp",j,
						"_gsea", gsea, "_genenorm",norm,"_featmethod",gsea_m)

			write.table(x=train_feat, file=paste0(out_file, "_train"), sep="\t", row.names=T, col.names=T)
			write.table(x=test_feat,  file=paste0(out_file, "_test"), sep="\t", row.names=T, col.names=T)
			write.table(x=gsea_set,  file=paste0(out_file, "_gsea"), sep="\t", row.names=F, col.names=T)
		}
	}
}

cat("Done\n")