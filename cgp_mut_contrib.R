# cgp_mut_contrib.R
library(data.table)
library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)

# Load files
cgp_mut <- fread("/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/mutations.txt")
cgp_new <- readRDS("~/Documents/Rotation/PIPELINES/CGP_FILES/082916_cgp_new.rds")

# Prep
# cgp_mut$gene_mut <- with(cgp_mut, paste0(Gene, "_", AA))
# cgp_mut[,mut_count:=length(unique(SAMPLE)), by="gene_mut"]
# cgp_mut <- cgp_mut[mut_count>=4,]

cgp_mut[,mut_count:=length(unique(SAMPLE)), by="Gene"]
cgp_mut <- cgp_mut[mut_count>=4]
cgp_mut[,sample_gene_mut:=length(cDNA), by=c("SAMPLE", "Gene")]
cgp_mut[,sample_n_mut:=length(cDNA), by="SAMPLE"]

# Process
drugs       <- c("Docetaxel", "Erlotinib", "Cisplatin", "Bortezomib", 
			 		"5-Fluorouracil", "Gemcitabine", "Temozolomide")

all_samples <- intersect(unique(cgp_new$cell_name), unique(cgp_mut$SAMPLE))
cgp_new     <- cgp_new[cell_name %in% all_samples,]
cgp_mut     <- cgp_mut[SAMPLE %in% all_samples,]

# cgp_mut$count <- 1
cgp_mut$count <- with(cgp_mut, sample_gene_mut/sample_n_mut)
cgp_mut     <- unique(cgp_mut[,c("SAMPLE", "Gene", "count"),with=F])
cgp_mut     <- acast(cgp_mut, Gene~SAMPLE, value.var="count", fill=0)

# Plot per compound
for (d in drugs){
	cat(paste0(d, "\n"))
	temp_table   <- cgp_new[Compound==d,][order(IC50)]

	# Classify sample cells
	sensitive    <- head(temp_table$cell_name, 15)
	resistant    <- tail(temp_table$cell_name, 55)
	temp_table   <- temp_table[cell_name %in% c(sensitive, resistant)]
	temp_mut     <- cgp_mut[,temp_table$cell_name]

	# Remove ZERO-features
	temp_mut     <- temp_mut[apply(temp_mut, 1, sum)!=0,]
	print(dim(temp_mut))

	# Obtain clustering
	sample_order <- hclust(dist(t(temp_mut)))[["order"]]
	temp_mut     <- temp_mut[,sample_order]
	gene_order 	 <- hclust(dist(temp_mut))[["order"]]
	temp_mut     <- temp_mut[gene_order,]

	# Table format
 	temp_mut     <- data.table(melt(temp_mut))
 	sample_cols  <- ifelse(temp_table$cell_name[sample_order] %in% sensitive, "red", "black")

 	# Plot
 	print(temp_mut[order(value)])
 	cat("Plotting\n")
 	pdf(paste0("~/Documents/FOLDER/LAB/GM/073017/cgp_mut/", d, "_mutation_phenotype_class_influence.pdf"), height=12, width=10)
 	print(
 		ggplot(temp_mut, aes(Var2, Var1, fill=value)) +
	 		geom_raster() + 
	 		theme(axis.text.x = element_text(colour=sample_cols, angle = 90, hjust = 1, size=8),
	            axis.text.y = element_text(size=4),
	            legend.position="top") +
	       	scale_fill_viridis(option = "C") +
	      	xlab("Cells") + ylab("Mutations") +
	      	ggtitle(paste0("Mutation map for: ", d))
 		)
 	
 	dev.off()
 }

 cat("Done")