# gsea_exp.R
# Obtains gene expression binary scores based on pre-processed gene-set modules

library(data.table)
library(reshape2)
library(parallel)

# Load files
# cgp_exp     <- readRDS("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/121216_cgp_exp.rds")
cgp_exp <- readRDS("/tigress/zamalloa/CGP_FILES/083016_cgp_exp.rds")
# gsea <- fread("/Users/jzamalloa/Documents/Rotation/PIPELINES/GSEA_FILES/cancer_sets") #c6 only
# gsea <- fread("/tigress/zamalloa/GSEA_FILES/processed_sets") #All sets combines

#"h", "c1", 
for (g in c("c2", "c3", "c4", "c5", "c6", "c7")){
	print(g)
	gsea <- fread(paste0("/tigress/zamalloa/GSEA_FILES/",g,"_sets"))

	# Process
	print(gsea)
	gsea <- gsea[genes %in% rownames(cgp_exp)] #Actual count of genes present
	print(gsea)

	gsea[,N:=length(unique(genes)), by="Gene_set"]
	gene_sets <- unique(gsea[N<200][N>20]$Gene_set) # Filtering step
	print(length(unique(gene_sets)))
	
	gene_sets <- t(combn(gene_sets, 2))
	names     <- apply(gene_sets, 1, function(x) paste0(x[1], "$", x[2]))

	# Test for two patients
	samples <- colnames(cgp_exp)
	count   <- 0

	nodes<-detectCores()
	print(nodes)
	cl<-makeCluster(nodes)
	setDefaultCluster(cl)
	clusterExport(cl, varlist=c("data.table", "as.data.table", "count" ,"samples", "gene_sets", "cgp_exp", "gsea", "names"),envir=environment())

	main_table <- lapply(samples, function(i) {
		print(count/length(samples))

		p_vals <- parApply(cl, gene_sets, 1, function(x) {
			first_genes  <- as.vector(cgp_exp[gsea[Gene_set==x[1]]$genes, i])
			second_genes <- as.vector(cgp_exp[gsea[Gene_set==x[2]]$genes, i])
			return (wilcox.test(first_genes, second_genes, alternative="greater")$p.value)
		})
		p_vals <- p.adjust(p_vals, method="fdr")
		count <<-count+1

		return(data.table(sample = i, gs = names, p_vals))
	})

	stopCluster(cl)
	main_table <- do.call(rbind, main_table)
	print(paste0("Done writing ", g))

	# saveRDS(main_table, "/tigress/zamalloa/GSEA_FILES/all_gsea_cgp_pvals.rds")
	saveRDS(main_table, paste0("/tigress/zamalloa/GSEA_FILES/",g,"_gsea_cgp_pvals.rds"))
}
