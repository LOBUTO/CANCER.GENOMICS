#cgp_mut_inf.R
library(data.table)
library(reshape2)
library(parallel)

# This script will be to get essential mutations only

# Load file 
cgp_mut   <- fread("/tigress/zamalloa/OBJECTS/mutations.txt", select = c(1,4,7))
cgp_exp   <- readRDS("/tigress/zamalloa/CGP_FILES/083016_cgp_exp.rds")

# Process effect of each mutation
cgp_mut$mutation <- paste0(cgp_mut$Gene, "_", cgp_mut$AA)
cgp_mut[,n:=length(unique(SAMPLE)), by="mutation"]
cgp_mut <- cgp_mut[n>2,] #Filter for minimum number of samples affected by mutation
mutations        <- unique(cgp_mut$mutation)

total <- length(mutations)
cat(paste0("total mutations: ", total, "\n"))

nodes<-detectCores()
cat(paste0("nodes: ", nodes, "\n"))
cl<-makeCluster(nodes)
setDefaultCluster(cl)
clusterExport(cl, varlist=c("data.table", "as.data.table", "cgp_mut", "cgp_exp", "mutations", "total", "count"),envir=environment())

count <- 0
main_table <- lapply(mutations, function(i) {
	cat(paste0(count, "\n"))

	# Obtain cells
	pos_cells  <- unique(cgp_mut[mutation==i,]$SAMPLE)
	neg_cells  <- unique(cgp_mut[mutation!=i,]$SAMPLE)

	# Obtain influence of each gene in mutation by a wilcoxon test between mutated and non-mutated samples
	p_vals     <- parApply(cl, cgp_exp, 1, function(x) wilcox.test(x[pos_cells], x[neg_cells])$p.value) 
	p_vals     <- p.adjust(p_vals, method="fdr")

	# Keep genes that are differentially expressed
	diff_genes <- names(p_vals)[p_vals < 0.1]

	if (length(diff_genes)>1){
		return(data.table(mutation=i, diff_genes))
	}
	
	count <<- count + 1
	write.table("log", paste(count,"\n"), append=T, quote=F, row.names = F, col.names = F, sep="\t")
})

stopCluster(cl)
main_table <- do.call(rbind, main_table)

# Store results 
saveRDS(main_table, "/tigress/zamalloa/CGP_FILES/cgp_mut_inf.rds")
print("Done writing")