# gs_plot_1.R

library(data.table)
library(reshape2)
library(ggplot2)
library(viridis)
library(cowplot)
library(gridExtra)

Function_clust_order <- function(in_feat){

	print("calculating")
	sample_order <- hclust(dist(t(in_feat)))[["order"]]
	in_feat      <- in_feat[,sample_order]
	# feat_order   <- hclust(dist(in_feat))[["order"]]
	# in_feat      <- in_feat[feat_order,]

	return(list(in_feat, sample_order))
}

# Load files
drug       <- "docetaxel"
bor_gs     <- fread(paste0("GSEA_FILES/c4_gsea_", drug, "_both_T_pvals"))
bor_gs_gmv <- fread(paste0("GSEA_FILES/c4_gsea_", drug, "_both_T_ext_gmv_T_pvals"))
# bor        <- readRDS("OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds")[["feat_table_a"]][Compound=="Bortezomib"]
bor        <- readRDS("OBJECTS/030217_GEE_DOCETAXEL.rds")[["feat_table"]][Compound=="Docetaxel"]

# Process
bor_gs[,p_vals:=p.adjust(pvals), by="sample"]
bor_gs_gmv[,p_vals:=p.adjust(pvals), by="sample"]
bor_gs$binary     <- ifelse(bor_gs$p_vals < 0.0001, 1, 0)
bor_gs_gmv$binary <- ifelse(bor_gs_gmv$p_vals < 0.0001, 1, 0)

# Obtain clustering
print("casting")
all_samples <- as.character(unique(bor$cell_name))
bor_gs      <- acast(bor_gs, gs~sample, value.var="binary")
bor_gs_gmv  <- acast(bor_gs_gmv, gs~sample, value.var="binary")
bor_gs      <- bor_gs[,all_samples]
bor_gs_gmv  <- bor_gs_gmv[,all_samples]

print("clustering")
bor_gs      <- Function_clust_order(bor_gs)
bor_gs_gmv  <- Function_clust_order(bor_gs_gmv)

# Post-process for plotting
gs_cols     <- ifelse(bor$target[bor_gs[[2]]] == 1, "red", "black")
gs_gmv_cols <- ifelse(bor$target[bor_gs_gmv[[2]]] == 1, "red", "black")

print("melting")
bor_gs      <- data.table(melt(bor_gs[[1]]))
bor_gs_gmv  <- data.table(melt(bor_gs_gmv[[1]]))

# Plot and save
print("plotting")
pdf(paste0("GSEA_FILES/PLOTS/Geeleher_",drug,"_GS_NON_GMV_vs_GMV.pdf"), , width = 16, height = 10)
grid.arrange(
	ggplot(bor_gs, aes(Var2, Var1, fill=value)) +
      geom_raster() + theme_bw() +
      theme(axis.text.x = element_text(colour=gs_cols, angle = 90, hjust = 1, size=4),
            axis.text.y = element_text(size=4),
            legend.position="top") +
      scale_fill_viridis(option = "D", discrete = F) +
      xlab("Samples") + ylab("Gene set features") +
      ggtitle("Original expression"),

    ggplot(bor_gs_gmv, aes(Var2, Var1, fill=value)) +
      geom_raster() + theme_bw() +
      theme(axis.text.x = element_text(colour=gs_gmv_cols, angle = 90, hjust = 1, size=4),
            axis.text.y = element_text(size=4),
            legend.position="top") +
      scale_fill_viridis(option = "D", discrete = F) +
      xlab("Samples") + ylab("Gene set features") +
      ggtitle("GMV-normalized expression"),
	ncol=2
	)
dev.off()

print("Done")





