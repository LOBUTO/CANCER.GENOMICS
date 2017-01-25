#cgp_new_prep_mf.R

library(data.table)
library(reshape2)

# FUNCTIONS
Function_top_cell_morgan_bits_features_extracted_mf <- function(feats, exp_table, morgan_table, max_cells=10, max_bits=10, pic50_scaled=T,
                                                                pic50_class=F, genes=F, rebalance=F, pca=F, common_genes,
                                                                nci_spiked, nci_morgan, nci_exp, nci_new){

  # Extract most variable cell features
  if (pca==F){
    if (genes==F){
      cell_feat <- cor(cgp_exp, method = "spearman")
      cell_var  <- data.table(cells = colnames(cell_feat),
                              VAR   = apply(cell_feat, 2, var))
      top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

      cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)

    } else {
      cgp_exp   <- t(log(cgp_exp))
      gene_var  <- data.table(genes = colnames(cgp_exp),
                              VAR   = apply(cgp_exp, 2, var))
      top_genes <- gene_var[order(-VAR),]$genes[1:max_cells]

      cell_feat <- data.table(cgp_exp[,top_genes], keep.rownames = T)
    }

  } else{
    cgp_pca     <- prcomp(t(cgp_exp[common_genes,]), center=T, scale. = T) # Scaling previously shown to be essential for out of set accuracy
    cell_feat   <- cgp_pca$x[,1:max_cells]
    cell_feat   <- data.table(cell_feat, keep.rownames = T)
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features
  if (max_bits > 0){
    if (pca == F){
      top_bits  <- morgan_table[,list(VAR = var(value)), by="bit_pos"][order(-VAR)]$bit_pos[1:max_bits]
      drug_feat <- morgan_table[bit_pos %in% top_bits, ]
      drug_feat <- acast(drug_feat, Compound~bit_pos, value.var="value")
      drug_feat <- data.table(drug_feat, keep.rownames = T)

    } else {
      drug_feat <- acast(morgan_table, Compound~bit_pos, value.var="value")
      drug_pca  <- prcomp(drug_feat, center=F, scale. = F)
      drug_feat <- data.table(drug_pca$x[, 1:max_bits], keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  } else {
    drug_feat <- data.table()
  }

  # Construct feature table based on classification/regression
  if(pic50_class==F){

    if(pic50_scaled==T){
      # feats$NORM.pIC50 <- scale(feats$pIC50) # WHOLE RE-SCALED - MODIFIED!!
      feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
    } else {
      feat_table <- feats[, c("Compound", "cell_name", "pIC50"), with=F]
    }

  } else{
    feat_table <- feats[, c("Compound", "cell_name", "pic50_class"),with=F]
    feat_table <- feat_table[pic50_class!=2,]
  }

  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  ############################## DO WE NEED TO NCI60 SPIKE IT??########################################
  if (nci_spiked==T){

    # Spike cell features first
    top_cell_feat <- colnames(cell_feat)[2:ncol(cell_feat)]
    common_genes  <- intersect(rownames(cgp_exp), rownames(nci_exp))
    nci_cells     <- colnames(nci_exp)
    cell_spiked   <- cbind(cgp_exp[common_genes, top_cell_feat], nci_exp[common_genes, nci_cells])
    cell_spiked   <- cor(cell_spiked, method="pearson")
    cell_spiked   <- data.table(cell_spiked[,top_cell_feat], keep.rownames=T)

    setnames(cell_spiked, c("cell_name", colnames(cell_spiked)[2:ncol(cell_spiked)]))

    # Spike drug features next
    if (max_bits > 0){
      top_drug_feat <- colnames(drug_feat)[2:ncol(drug_feat)]
      spiked_bits   <- nci_morgan[position %in% top_bits,]
      spiked_bits   <- acast(spiked_bits, Compound~position, value.var="value")
      spiked_bits   <- data.table(spiked_bits[,top_drug_feat], keep.rownames = T)

      setnames(spiked_bits, c("Compound", colnames(spiked_bits)[2:ncol(spiked_bits)]))

    } else {
      spiked_bits <- data.table()
    }

    # Lastly spike target features
    if(pic50_class==F){

      # NOTE: KEEP IN MIND these activity measures may not be analogous to IC50 and may need to be adjusted to be comparable
      if(pic50_scaled==T){
        nci_feat_table <- nci_new[, c("NSC", "CELL", "SCALE.ACT"), with=F]
      } else{
        nci_feat_table <- nci_new[, c("NSC", "CELL", "ACT"), with=F]
      }

    } else{
      nci_feat_table <- nci_new[, c("NSC", "CELL", "act_class"), with=F]
      nci_feat_table <- nci_feat_table[act_class!=2,]
    }
    setnames(nci_feat_table, c("Compound", "cell_name", "NORM_pIC50"))

    # Spiking
    cell_feat  <- rbind(cell_feat, cell_spiked)
    drug_feat  <- rbind(drug_feat, spiked_bits)
    feat_table <- rbind(feat_table, nci_feat_table)

    setkey(cell_feat)
    setkey(drug_feat)
    setkey(feat_table)
    cell_feat  <- unique(cell_feat)
    drug_feat  <- unique(drug_feat)
    feat_table <- unique(feat_table)
  }
  #####################################################################################################

  # Merge to obtain total number of possible combinations
  if (max_bits > 0){
    feat_table <- merge(feat_table, drug_feat[, 1:2, with=F], by="Compound", allow.cartesian=TRUE)
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name", allow.cartesian=TRUE)

  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Do we need to rebalance based on compound counts
  if (rebalance == T){
    feat_table[ , N:=length(cell_name), by="Compound"]
    max_count <- max(feat_table$N)

    feat_extra <- lapply(unique(feat_table[N < max_count, ]$Compound), function(d) {

      temp_table  <- feat_table[Compound == d,]
      n_needed    <- max_count - unique(temp_table$N)
      sample_rows <- sample(1:nrow(temp_table), n_needed, replace = T)
      temp_table  <- temp_table[sample_rows,]

      return(temp_table)
      })
    feat_extra <- do.call(rbind, feat_extra)

    feat_table <- rbind(feat_table, feat_extra)
    feat_table$N <- NULL

  }

  # Extract indices for combinations
  if (max_bits > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else {
    drug_index <- c()
  }
  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$NORM_pIC50

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))
}

Function_top_cell_morgan_counts_features_extracted_mf <- function(feats, exp_table, morgan_table, max_cells=10, max_counts=10, met_scaled=F, pic50_scaled=T,
                                                      pic50_class=F, genes=F, rebalance=F, pca=F, common_genes,
                                                      nci_spiked, nci_morgan, nci_exp, nci_new){
  # Constructs feature table using morgan bits and cell correlations as features, but limiting to most variable features
  # Both max_cells and max_counts are based in terms of decreasing variance
  # NOTE: Function can now be used for classification or regression
  #       If classification is chosen, pic50-based classes are chosen (1/0), while 2 is discarded

  # Extract most variable cell features
  if (pca==F){
    if (genes==F){
      cell_feat <- cor(cgp_exp, method = "spearman")
      cell_var  <- data.table(cells = colnames(cell_feat),
                              VAR   = apply(cell_feat, 2, var))
      top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]

      cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)

    } else {
      cgp_exp   <- t(log(cgp_exp))
      gene_var  <- data.table(genes = colnames(cgp_exp),
                              VAR   = apply(cgp_exp, 2, var))
      top_genes <- gene_var[order(-VAR),]$genes[1:max_cells]

      cell_feat <- data.table(cgp_exp[,top_genes], keep.rownames = T)
    }

  } else{
    cgp_pca     <- prcomp(t(cgp_exp[common_genes,]), center=T, scale. = T) # Scaling previously shown to be essential for out of set accuracy
    cell_feat   <- cgp_pca$x[,1:max_cells]
    cell_feat   <- data.table(cell_feat, keep.rownames = T)
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  # Extract most variable drug features (if needed)
  if (max_counts > 0 ){
    if (pca==F) {
      top_counts <- apply(acast(morgan_table, Compound~Substructure, value.var="Counts", fill=0), 2, var)
      top_counts <- names(sort(top_counts, decreasing = T))[1:max_counts]

      drug_feat  <- morgan_table[Substructure %in% top_counts, ]
      drug_feat  <- acast(drug_feat, Compound~Substructure, value.var="Counts", fill=0)
      drug_feat  <- data.table(drug_feat, keep.rownames = T)

    } else {
      drug_feat  <- acast(morgan_table, Compound~Substructure, value.var="Counts", fill=0)
      drug_pca   <- prcomp(drug_feat, center=F, scale. = F)
      drug_feat  <- data.table(drug_pca$x[, 1:max_counts], keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  } else {
    drug_feat  <- data.table()
  }

  # Construct feature table based on classification/regression
  if(pic50_class==F){

    if(pic50_scaled==T){
      # feats$NORM.pIC50 <- scale(feats$pIC50) # WHOLE RE-SCALED - MODIFIED!!
      feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
    } else {
      feat_table <- feats[, c("Compound", "cell_name", "pIC50"), with=F]
    }

  } else{
    feat_table <- feats[, c("Compound", "cell_name", "pic50_class"),with=F]
    feat_table <- feat_table[pic50_class!=2,]
  }

  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))

  ############################## DO WE NEED TO NCI60 SPIKE IT??########################################
  if (nci_spiked==T){

    # Spike cell features first
    top_cell_feat <- colnames(cell_feat)[2:ncol(cell_feat)]
    common_genes  <- intersect(rownames(cgp_exp), rownames(nci_exp))
    nci_cells     <- colnames(nci_exp)
    cell_spiked   <- cbind(cgp_exp[common_genes, top_cell_feat], nci_exp[common_genes, nci_cells])
    cell_spiked   <- cor(cell_spiked, method="pearson")
    cell_spiked   <- data.table(cell_spiked[,top_cell_feat], keep.rownames=T)

    setnames(cell_spiked, c("cell_name", colnames(cell_spiked)[2:ncol(cell_spiked)]))

    # Spike drug features next
    if (max_bits > 0){
      top_drug_feat <- colnames(drug_feat)[2:ncol(drug_feat)]
      spiked_counts   <- nci_morgan[substructure %in% top_counts,]
      spiked_counts   <- acast(spiked_counts, Compound~substructure, value.var="value", fill=0)
      spiked_counts   <- data.table(spiked_bits[,top_drug_feat], keep.rownames = T)

      setnames(spiked_counts, c("Compound", colnames(spiked_counts)[2:ncol(spiked_counts)]))

    } else {
      spiked_counts <- data.table()
    }

    # Lastly spike target features
    if(pic50_class==F){

      # NOTE: KEEP IN MIND these activity measures may not be analogous to IC50 and may need to be adjusted to be comparable
      if(pic50_scaled==T){
        nci_feat_table <- nci_new[, c("NSC", "CELL", "SCALE.ACT"), with=F]
      } else{
        nci_feat_table <- nci_new[, c("NSC", "CELL", "ACT"), with=F]
      }

    } else{
      nci_feat_table <- nci_feat_table[, c("NSC", "CELL", "act_class"), with=F]
      nci_feat_table <- nci_feat_table[act_class!=2,]
    }
    setnames(nci_feat_table, c("Compound", "cell_name", "NORM_pIC50"))

    # Spiking
    cell_feat  <- rbind(cell_feat, cell_spiked)
    drug_feat  <- rbind(drug_feat, spiked_bits)
    feat_table <- rbind(feat_table, nci_feat_table)
  }
  #####################################################################################################

  # Merge to obtain total number of possible combinations
  if (max_counts > 0){
    feat_table <- merge(feat_table, drug_feat[, 1:2, with=F], by="Compound")
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name")
  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Do we need to rebalance based on compound counts
  if (rebalance == T){
    feat_table[ , N:=length(cell_name), by="Compound"]
    max_count <- max(feat_table$N)

    feat_extra <- lapply(unique(feat_table[N < max_count, ]$Compound), function(d) {

      temp_table  <- feat_table[Compound == d,]
      n_needed    <- max_count - unique(temp_table$N)
      sample_rows <- sample(1:nrow(temp_table), n_needed, replace = T)
      temp_table  <- temp_table[sample_rows,]

      return(temp_table)
      })
    feat_extra <- do.call(rbind, feat_extra)

    feat_table <- rbind(feat_table, feat_extra)
    feat_table$N <- NULL

  }

  # Extract indices for combinations
  if (max_counts > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else{
    drug_index <- c()
  }

  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$NORM_pIC50

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))

}

Function.range.0.1 <- function(x){
  scale.1 <-  (x-min(x))/(max(x)-min(x))
  return(scale.1)
}

Function_load_morgan_bits <- function(morgan=T){
  if (morgan==T){
    morgan_nci_bits   <- fread(paste0(in_morgan, "NCI_MORGAN_BITS_r_16_b_2048.txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==16 & bits==2048,]
    morgan_nci_bits$position <- paste0("mcf_", morgan_nci_bits$position)
  } else{
    morgan_nci_bits   <- c()
  }
  return(morgan_nci_bits)
}

Function_load_morgan_counts <- function(morgan=T){
  if (morgan==T){
    morgan_nci_counts <- fread(paste0(in_morgan, "NCI_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "numeric", "numeric"))[radius==12,]
  } else{
    morgan_nci_counts <- c()
  }
}

# LOAD INPUT
args        <- commandArgs(trailingOnly = TRUE)

max_cells   <- as.numeric(args[1])
max_drugs   <- as.numeric(args[2])
file_name   <- args[3]
met_type    <- args[4]
class_mlp   <- as.logical(args[5])
samples     <- args[6]
genes       <- as.logical(args[7])
batch_norm  <- args[8]
pca         <- as.logical(args[9])
radii_set   <- as.numeric(args[10])
bit_set     <- as.numeric(args[11])
nci_spiked  <- F#as.logical(args[8])

# if (args[7] == "lab"){
#   in_folder   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/" #For lab
#   in_morgan   <- "/home/zamalloa/Documents/FOLDER/MORGAN_FILES/"
#   out_folder  <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"
# } else{
in_folder    <- "/tigress/zamalloa/CGP_FILES/" #For tigress
in_morgan    <- "/tigress/zamalloa/MORGAN_FILES/"
in_objects   <- "/tigress/zamalloa/OBJECTS/"
tcga_objects <- "/tigress/zamalloa/TCGA_FILES/"
out_folder   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/"
# }

# LOAD DATA
DRUGS.MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
cgp_new           <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
ccle_exp          <- readRDS(paste0(in_folder, "121116_ccle_exp.rds"))
nci_exp           <- readRDS(paste0(in_objects, "121216_nci60_exp.rds"))

if (batch_norm=="cgp_nci60"){
  print("batch_norm nci60")
  cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_nci60_exp_b_norm.rds"))[["EXP.1"]]

} else if (batch_norm=="cgp_ccle"){
  print("batch_norm ccle")
  cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_ccle_exp_b_norm.rds"))[["EXP.1"]]

} else if (batch_norm=="tcga_brca"){
  print("tcga_brca")
  cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.1"]]

} else if (batch_norm=="tcga_coad"){
  print("tcga_coad")
  cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.1"]]

} else if (batch_norm=="tcga_luad"){
  print("tcga_luad")
  cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.1"]]

} else if (batch_norm=="tcga_stad"){
  print("tcga_stad")
  cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.1"]]

} else{
  print("None")
  cgp_exp           <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))
}

common_genes      <- intersect(intersect(rownames(cgp_exp), rownames(nci_exp)), rownames(ccle_exp))
morgan_bits       <- fread(paste0(in_morgan, "CGP_MORGAN_BITS.txt"))[radius==radii_set & bits==bit_set,]  #8 and 2048 normal setting
morgan_counts     <- fread(paste0(in_morgan, "CGP_MORGAN_COUNTS.txt"),
                       colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,] #12 normal setting

nci_exp           <- c()#readRDS("/tigress/zamalloa/OBJECTS/061916.NCI60.EXP.rds")
nci_new           <- c()#readRDS("/tigress/zamalloa/OBJECTS/120916_nci60_new.rds")
#################################################### EXECUTE ####################################################
if (met_type=="drug_cor"){
  feat_table <- Function_top_cell_drug_features_extracted_mf(cgp_new, cgp_exp, DRUGS.MET.PROFILE,
                                                          max_cells = max_cells, max_drugs = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T)

  drug_sim   <- data.table(melt(cor(acast(DRUGS.MET.PROFILE, METABOLITE~DRUG, value.var = "TC"))))
  drug_sim   <- drug_sim[Var1!=Var2,]

} else if (met_type=="morgan_bits"){
  feat_table <- Function_top_cell_morgan_bits_features_extracted_mf(cgp_new, cgp_exp, morgan_bits,
                                                          max_cells = max_cells, max_bits = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T,
                                                          genes = genes, rebalance = F, nci_spiked = nci_spiked,
                                                          nci_morgan = Function_load_morgan_bits(nci_spiked),
                                                          nci_exp = nci_exp, nci_new = nci_new,
                                                          pca = pca, common_genes = common_genes)

  drug_sim   <- data.table(melt(cor(acast(morgan_bits, bit_pos~Compound, value.var = "value"))))
  drug_sim   <- drug_sim[Var1!=Var2,]

} else if (met_type=="morgan_counts"){
  feat_table <- Function_top_cell_morgan_counts_features_extracted_mf(cgp_new, cgp_exp, morgan_counts,
                                                          max_cells = max_cells, max_counts = max_drugs,
                                                          pic50_class = class_mlp, pic50_scaled = T,
                                                          genes = genes, rebalance = F, nci_spiked = nci_spiked,
                                                          nci_morgan = Function_load_morgan_counts(nci_spiked),
                                                          nci_exp = nci_exp, nci_new = nci_new,
                                                          pca = pca, common_genes = common_genes)

  drug_sim   <- data.table(melt(cor(acast(morgan_counts, Substructure~Compound, value.var = "Counts", fill = 0))))
  drug_sim   <- drug_sim[Var1!=Var2,]

}

# Split tables depending on all or drug samples
if (samples == "all"){

  set.seed(1234)
  train_rows       <- sample(1:length(feat_table$target), length(feat_table$target)*0.7)
  testing_rows     <- setdiff(1:length(feat_table$target), train_rows)
  set.seed(1234)
  valid_rows       <- sample(testing_rows, length(testing_rows)*0.5)
  test_rows        <- setdiff(testing_rows, valid_rows)

} else if (samples == "whole"){
  #AKA all data training

  set.seed(1234)
  train_rows       <- sample(1:length(feat_table$target), length(feat_table$target)*0.8)
  valid_rows       <- setdiff(1:length(feat_table$target), train_rows)
  test_rows        <- valid_rows

} else if (samples == "all_split"){
  print("all_split")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  train_drug       <- sample(all_drugs, length(all_drugs)*0.7)
  testing_drugs    <- setdiff(all_drugs, train_drug) #all but training
  set.seed(1234)
  valid_drug       <- sample(testing_drugs, length(testing_drugs)*0.5)
  test_drug        <- setdiff(testing_drugs, valid_drug)

  train_rows       <- which(feat_table$feat_table$Compound %in% train_drug)
  valid_rows       <- which(feat_table$feat_table$Compound %in% valid_drug)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else if (samples == "semi_split"){
  print ("semi_split")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  #test_drug        <- sample(all_drugs, 5)
  #test_drug        <- c("GDC0449", "MS-275", "PAC-1", "RDEA119", "TG101348")
  test_drug        <- sample(all_drugs, 0.2*length(all_drugs))
  # test_drug        <- c("17-AAG", "BHG712", "Bleomycin", "BX-912", "CEP-701", "CH5424802",
  #                       "CMK", "CP724714", "Docetaxel", "EHT 1864", "FR-180204", "FTI-277",
  #                       "GSK1070916", "GSK2126458", "IPA-3", "Ispinesib Mesylate", "JQ1",
  #                       "KU-55933", "Lapatinib", "LAQ824", "Midostaurin", "NSC-87877",
  #                       "NU-7441", "Nutlin-3a (-)", "Parthenolide", "PHA-793887", "Rapamycin",
  #                       "Sunitinib", "TAK-715", "Temozolomide", "Tipifarnib", "Trametinib",
  #                       "Tubastatin A", "YM201636")
  testing_drugs    <- setdiff(all_drugs, test_drug)

  testing_rows     <- which(feat_table$feat_table$Compound %in% testing_drugs)
  set.seed(1234)
  train_rows       <- sample(testing_rows, length(testing_rows)*0.80)
  valid_rows       <- setdiff(testing_rows, train_rows)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else if(samples == "all_split_balanced"){
  print("all_split_balanced")

  all_drugs        <- unique(feat_table$feat_table$Compound)
  set.seed(1234)
  test_drug        <- sample(all_drugs, length(all_drugs)*0.15)
  testing_drugs    <- setdiff(all_drugs, test_drug)

  # Sample based on similarity to testing for balanaced validation
  drug_sim         <- drug_sim[Var1 %in% test_drug,][Var2 %in% testing_drugs,]
  drug_sim         <- drug_sim[,list(mean_sim = mean(value)), by="Var2"] # Mean similarity to target
  valid_drug       <- sample(drug_sim$Var2, length(testing_drugs)*0.25 , prob= Function.range.0.1(drug_sim$mean_sim))
  train_drug       <- setdiff(testing_drugs, valid_drug)

  train_rows       <- which(feat_table$feat_table$Compound %in% train_drug)
  valid_rows       <- which(feat_table$feat_table$Compound %in% valid_drug)
  test_rows        <- which(feat_table$feat_table$Compound %in% test_drug)

} else {

  test_rows        <- which(feat_table$feat_table$Compound == samples)
  testing_rows     <- which(feat_table$feat_table$Compound != samples)

  set.seed(1234)
  train_rows       <- sample(testing_rows, length(testing_rows)*0.8)
  valid_rows       <- setdiff(testing_rows, train_rows)
}

train_drug       <- unique(names(feat_table$drug_index[train_rows]))
valid_drug       <- unique(names(feat_table$drug_index[valid_rows]))
test_drug        <- unique(names(feat_table$drug_index[test_rows]))

train_cell       <- unique(names(feat_table$cell_index[train_rows]))
valid_cell       <- unique(names(feat_table$cell_index[valid_rows]))
test_cell        <- unique(names(feat_table$cell_index[test_rows]))

if (max_drugs > 0){
  train_drug_table <- feat_table$drug_feat[Compound %in% train_drug,]
  valid_drug_table <- feat_table$drug_feat[Compound %in% valid_drug,]
  test_drug_table  <- feat_table$drug_feat[Compound %in% test_drug,]
} else{
  train_drug_table <- data.table()
  valid_drug_table <- data.table()
  test_drug_table  <- data.table()
}

train_cell_table <- feat_table$cell_feat[cell_name %in% train_cell,]
valid_cell_table <- feat_table$cell_feat[cell_name %in% valid_cell,]
test_cell_table  <- feat_table$cell_feat[cell_name %in% test_cell,]

# Obtain new index
train_feat_table <- feat_table$feat_table[train_rows,]
valid_feat_table <- feat_table$feat_table[valid_rows,]
test_feat_table  <- feat_table$feat_table[test_rows,]

if (max_drugs > 0){
  train_drug_index <- unlist(sapply(train_feat_table$Compound, function(x) which(x==train_drug_table$Compound)))
  valid_drug_index <- unlist(sapply(valid_feat_table$Compound, function(x) which(x==valid_drug_table$Compound)))
  test_drug_index  <- unlist(sapply(test_feat_table$Compound,  function(x) which(x==test_drug_table$Compound)))
} else {
  train_drug_index <- 0
  valid_drug_index <- 0
  test_drug_index  <- 0
}

train_cell_index <- unlist(sapply(train_feat_table$cell_name, function(x) which(x==train_cell_table$cell_name)))
valid_cell_index <- unlist(sapply(valid_feat_table$cell_name, function(x) which(x==valid_cell_table$cell_name)))
test_cell_index  <- unlist(sapply(test_feat_table$cell_name,  function(x) which(x==test_cell_table$cell_name)))

train_target     <- train_feat_table$NORM_pIC50
valid_target     <- valid_feat_table$NORM_pIC50
test_target      <- test_feat_table$NORM_pIC50

train_index      <- data.table(drug = train_drug_index - 1,
                               cell = train_cell_index - 1,
                               NORM_pIC50 = train_target)

valid_index      <- data.table(drug = valid_drug_index - 1,
                               cell = valid_cell_index - 1,
                               NORM_pIC50 = valid_target)

test_index       <- data.table(drug = test_drug_index - 1,
                               cell = test_cell_index - 1,
                               NORM_pIC50 = test_target)

# Scale tables with respect to training set (for Multiplicative_fusion, only scale drug features if they are continous variables)
if (met_type=="drug_cor"){
  not_scaling <- 1:3
  scaling     <- (max(not_scaling) + 1):ncol(feat_table)

} else if (met_type=="morgan_bits"){
  not_scaling <- 1
  scaling     <- 2:ncol(train_cell_table)

} else if (met_type=="morgan_counts"){
  not_scaling <- 1
  scaling     <- 2:ncol(train_cell_table)
}

scale_train_cell  <- scale(train_cell_table[, scaling, with=F])
train_cell_table  <- cbind(train_cell_table[, not_scaling, with=F], scale_train_cell)

train_cell_mean   <- attributes(scale_train_cell)$`scaled:center`
train_cell_sd     <- attributes(scale_train_cell)$`scaled:scale`

valid_cell_scale  <- sweep(valid_cell_table[, scaling, with=F], 2, train_cell_mean, "-")
valid_cell_scale  <- sweep(valid_cell_scale, 2, train_cell_sd, "/")
valid_cell_table  <- cbind(valid_cell_table[, not_scaling, with=F],
                     valid_cell_scale)

test_cell_scale  <- sweep(test_cell_table[, scaling, with=F], 2, train_cell_mean, "-")
test_cell_scale  <- sweep(test_cell_scale, 2, train_cell_sd, "/")
test_cell_table  <- cbind(test_cell_table[, not_scaling, with=F],
                     test_cell_scale)

# WRITE TABLES
write.table(train_index, paste0(out_folder, file_name, "_train_index"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(valid_index, paste0(out_folder, file_name, "_valid_index"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(test_index,  paste0(out_folder, file_name, "_test_index"),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(train_drug_table, paste0(out_folder, file_name, "_train_drug"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(valid_drug_table, paste0(out_folder, file_name, "_valid_drug"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(test_drug_table, paste0(out_folder, file_name,  "_test_drug"),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(train_cell_table, paste0(out_folder, file_name, "_train_cell"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(valid_cell_table, paste0(out_folder, file_name, "_valid_cell"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(test_cell_table, paste0(out_folder, file_name,  "_test_cell"),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing sel tables")
