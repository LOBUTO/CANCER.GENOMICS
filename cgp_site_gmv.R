# cgp_site_gmv.R
# Add on script to perform Site-based (cancer tisse) gmv normalization on the cgp data set

library(data.table)
library(reshape2)
library(Hmisc)

# Load files
cgp_new <- readRDS("CGP_FILES/082916_cgp_new.rds")
cgp_exp <- readRDS("CGP_FILES/083016_cgp_exp.rds")

# Process
Function_process_gmv <- function(exp_table, target_table){

	common_cells <- intersect(colnames(exp_table), unique(target_table$cell_name))
  	target_table <- target_table[cell_name %in% common_cells,]

  	filter_cells <- unique(target_table[,c("cell_name", "Site"),with=F])
  	filter_cells[,N:=length(unique(cell_name)), by="Site"]

  	# Obtain gmv-site weighted expression
  	print("getting expression")
  	spec_mean    <- apply(exp_table[,filter_cells$cell_name], 1, function(x) wtd.mean(x, weights=(1/filter_cells$N), normwt = T))
  	spec_var     <- apply(exp_table[,filter_cells$cell_name], 1, function(x) wtd.var(x, weights=(1/filter_cells$N), normwt = T))
  	exp_gmv_spec <- (exp_table[,filter_cells$cell_name] - spec_mean) / spec_var

  	# Return
  	return(exp_gmv_spec)
}

gmv_exp <- Function_process_gmv(cgp_exp, cgp_new)

write.table(data.table(gmv_exp, keep.rownames = T), file="CGP_FILES/cgp_site_gmv_exp",  sep="\t", quote=F, row.names=F, col.names=T)
print("Done")
