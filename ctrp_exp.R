# ctrp_exp.R
# Parses ctrp expression data and activity data
# Expression for version 2.1, activity for all

library(data.table)
library(reshape2)

# EXPRESSION DATA #
# in_folder <- "/Users/jzamalloa/Documents/Rotation/DATABASES/CTRP/CTRPv2.1_2016_pub_NatChemBiol_12_109/"

# #Load files
# print("loading files")
# exp_file   <- fread(paste0(in_folder, "v21.data.gex_avg_log2.txt"), header=T)
# gene_file  <- fread(paste0(in_folder, "v21.meta.gex_features.txt"), header=T)[,c("idx_gene_feature", "gene_primary_name"),with=F]
# cell_file  <- fread(paste0(in_folder, "v21.meta.per_cell_line.txt"), header=T)[,c("master_ccl_id", "ccl_name"),with=F]

# #Parse
# print("parsing files")
# exp_file   <- merge(exp_file, cell_file, by="master_ccl_id")[,-c("master_ccl_id"),with=F]
# exp_file   <- merge(exp_file, gene_file, by="idx_gene_feature")[,-c("idx_gene_feature"),with=F]
# exp_file   <- exp_file[,list(exp=max(mrna_expression_avg_log2)), by=c("ccl_name","gene_primary_name")]

# exp_file   <- acast(exp_file, gene_primary_name~ccl_name, value.var = "exp")
# print(dim(exp_file))

# #Store
# saveRDS(exp_file, "/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/ctrp_exp_2.1.rds")
# print("Done")

# ACTIVITY DATA #
#2.0
in_folder <- "/Users/jzamalloa/Documents/Rotation/DATABASES/CTRP/"
exp_2.0   <- fread(paste0(in_folder, "CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_experiment.txt"))[,c("experiment_id", "master_ccl_id"),with=F]
cell_2.0  <- fread(paste0(in_folder, "CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_cell_line.txt"))[,c("master_ccl_id", "ccl_name"),with=F]
cpd_2.0   <- fread(paste0(in_folder, "CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt"))[,c("master_cpd_id", "cpd_name"),with=F]
act_2.0   <- fread(paste0(in_folder, "CTRPv2.0_2015_ctd2_ExpandedDataset/v20.data.curves_post_qc.txt"))[,c("experiment_id", "area_under_curve", "apparent_ec50_umol", "master_cpd_id"),with=F]

table_2.0 <- merge(exp_2.0, cell_2.0, by="master_ccl_id")[,-c("master_ccl_id"),with=F]
table_2.0 <- merge(act_2.0, table_2.0, by="experiment_id", allow.cartesian=TRUE)[,-c("experiment_id"),with=F]
table_2.0 <- merge(table_2.0, cpd_2.0, by="master_cpd_id", allow.cartesian=TRUE)
setkey(table_2.0)
table_2.0 <- unique(table_2.0)

#2.1
cell_2.1  <- fread(paste0(in_folder, "CTRPv2.1_2016_pub_NatChemBiol_12_109/v21.meta.per_cell_line.txt"))[,c("master_ccl_id", "ccl_name"),with=F]
cpd_2.1   <- fread(paste0(in_folder, "CTRPv2.1_2016_pub_NatChemBiol_12_109/v21.meta.per_compound.txt"))[,c("master_cpd_id", "cpd_name"),with=F]
act_2.1   <- fread(paste0(in_folder, "CTRPv2.1_2016_pub_NatChemBiol_12_109/v21.data.auc_sensitivities.txt"))[,c("area_under_curve", "master_cpd_id", "master_ccl_id"),with=F]

table_2.1 <- merge(act_2.1, cell_2.1, by=c("master_ccl_id"))[,-c("master_ccl_id"),with=F]
table_2.1 <- merge(table_2.1, cpd_2.1, by=c("master_cpd_id"))
setkey(table_2.1)
table_2.1 <- unique(table_2.1)

#2.2
cell_2.2  <- fread(paste0(in_folder, "CTRPv2.2_2015_pub_CancerDisc_5_1210/v22.meta.per_cell_line.txt"))[,c("index_ccl","ccl_name"),with=F]
cpd_2.2   <- fread(paste0(in_folder, "CTRPv2.2_2015_pub_CancerDisc_5_1210/v22.meta.per_compound.txt"))[,c("index_cpd", "cpd_name", "master_cpd_id"),with=F]
act_2.2   <- fread(paste0(in_folder, "CTRPv2.2_2015_pub_CancerDisc_5_1210/v22.data.auc_sensitivities.txt"))

table_2.2 <- merge(act_2.2, cell_2.2, by="index_ccl")[,-c("index_ccl"),with=F]
table_2.2 <- merge(table_2.2, cpd_2.2, by="index_cpd")[,-c("index_cpd"),with=F]
setkey(table_2.2)
table_2.2 <- unique(table_2.2)

#Return
saveRDS(list(ctrp_2.0 = table_2.0,
			 ctrp_2.1 = table_2.1,
			 ctrp_2.2 = table_2.2),
		file= "/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/ctrp_tables.rds")
print("Done")