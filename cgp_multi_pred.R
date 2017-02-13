# cgp_multi_pred.R
# Script to prep outer datasets to be predicted by CGP based models (mf-based models)
# i.e. NCI-60, CCLE, TCGA

library(data.table)
library(reshape2)

############################################################################################################
# FUNCTIONS
Function_prep_new <- function(target_new, type="", class = F, scaled=T){
  # Reformats input so that its ready for test table building
  # NOTE: For testing purposes we will only test 1000 compounds for nci_60 at the moment (Modified, see code)

  if (type=="nci_60"){

    # nci_filter <- target_new[act_class!=2,][,list(MEAN=mean(act_class), SUM=sum(act_class), COUNT=length(act_class)), by="NSC"]
    # nci_filter <- nci_filter[MEAN > 0.4 & MEAN < 0.7][COUNT>40,]$NSC

    ### CUSTOMIZED FILTER ###
    nci_cgp    <- c('740','3061','19893','26980','49842','83265','91874','94600','119875','122758','123127','124663','125066',
                    '125973','141540','157035','174939','180973','203821','226080','226080','252844','287459','330507','330507',
                    '339555','613327','628503','656576','673596','681239','683864','701554','701852','715055','718781','727989',
                    '732517','734950','743414','743444','745750','747599','747856','747971','748799','751249','753146','755880',
                    '756645','756714','757306','757441','757804','758184','758612','758645','759155','759850','759852','759856',
                    '759877')
  #  nci_filter  <- sample(setdiff(unique(target_new$NSC), nci_cgp), 50)
  #  nci_filter  <- c(nci_filter, nci_cgp)
   nci_filter  <- c('11963','122758','123127','125066','125973','141540','157035','166052','180973','19893','226080','252844',
                          '26980','287459','3061','330507','333329','339555','34864','353647','4307','49842','613327','621360','628503',
                          '634604','637921','639793','640681','641626','642728','644860','645374','648581','650914','653614','656251',
                          '656576','664331','666414','673596','674276','675554','675867','679743','681239','681446','682770','694490',
                          '697266','701554','701852','702687','711066','711886','715055','716139','716335','718388','718666','718781',
                          '720995','725063','727989','730248','732517','734950','740','743414','743444','745750','747599','747856','747971',
                          '750873','751249','755880','756645','757306','757441','757804','758612','758645','759155','759852','759856',
                          '759877','764320','764567','83265','87406','91874','94600','174939','617610','631137','683864','705297','709562',
                          '723897','732918')
   #nci_filter  <- setdiff(nci_filter, nci_cgp)

   ##########################

    target_new <- target_new[NSC %in% nci_filter,] # NECESSARY FILTER FOR OUT OF SET TESTING

    if (class == T){
      target_new <- target_new[,c("NSC", "cell_name", "act_class"), with=F]
      target_new <- target_new[act_class!=2,]
    } else {
      if (scaled == T){
        # target_new$SCALE.ACT <- scale(target_new$ACT) # WHOLE RE-SCALED - MODIFIED!!
        target_new <- target_new[,c("NSC", "cell_name", "SCALE.ACT"), with=F] #MAY NEED TO MODIFY!
      } else{
        target_new <- target_new[,c("NSC", "cell_name", "ACT"), with=F] #MAY NEED TO MODIFY!
      }
    }

  } else if (type=="ccle"){
    if (class == T){
      target_new  <- target_new[,c("Compound", "cell_name", "pic50_CLASS"), with=F]
      target_new  <- target_new[pic50_CLASS!=2,]

      ccle_filter <- target_new[,list(MEAN = mean(pic50_CLASS)), by="Compound"]
      ccle_filter <- ccle_filter[MEAN > 0.35 & MEAN < 0.65,]$Compound  # NECESSARY FILTER FOR OUT OF SET TESTING
      target_new  <- target_new[Compound %in% ccle_filter,]

    } else {
      if (scaled == T){
        # target_new$NORM.pIC50 <- scale(target_new$pIC50) # WHOLE RE-SCALED - MODIFIED!!
        target_new <- target_new[,c("Compound", "cell_name", "NORM.pIC50"), with=F]
      } else {
        target_new <- target_new[,c("Compound", "cell_name", "pIC50"), with=F]
      }
    }

  } else if (type=="cgp"){
    if (class == T){
      target_new <- target_new[,c("Compound", "cell_name", "pic50_class"),with=F]
      target_new <- target_new[pic50_class!=2,]
    } else {
      if (scaled == T){
        # target_new$NORM.pIC50 <- scale(target_new$pIC50) # WHOLE RE-SCALED - MODIFIED!!
        target_new <- target_new[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
      } else {
        target_new <- target_new[, c("Compound", "cell_name", "pIC50"), with=F]
      }
    }
  } else if (type=="tcga"){
    target_new    <- target_new[, c("bcr_patient_barcode", "drug_name", "Binary"), with=F]
    target_filter <- target_new[,list(N = length(bcr_patient_barcode)), by="drug_name"]
    target_filter <- target_filter[N>25,]$drug_name # MAY NEED TO INCREASE THRESHOLD
    target_new    <- target_new[drug_name %in% target_filter,] # NECESSARY FILTER FOR OUT OF SET TESTING

    target_new_a  <- copy(target_new)
    target_new_b  <- copy(target_new)
    target_new_a$cell_name <- paste0(gsub("-", ".", target_new_a$bcr_patient_barcode), ".01A")
    target_new_b$cell_name <- paste0(gsub("-", ".", target_new_b$bcr_patient_barcode), ".01B")
    target_new    <- rbind(target_new_a, target_new_b)
    target_new    <- target_new[, c("drug_name", "cell_name", "Binary"), with=F]

    target_new$Binary <- ifelse(target_new$Binary == "Effective", 1, 0)

  } else if (type == "tcga_all" | type == "tcga_multi"){

    if (type=="tcga_all"){
      tcga_compounds <- target_new[,list(samples=length(bcr_patient_barcode), cancers=length(unique(Cancer))), by="drug_name"]
      tcga_compounds <- tcga_compounds[samples>40,]$drug_name
      target_new     <- target_new[drug_name %in% tcga_compounds,]
    } else {
      target_new     <- target_new[,list(drug_name = paste0(sort(drug_name), collapse="+"),
                                         Binary = unique(Binary)), by="bcr_patient_barcode"]
    }

    target_new_a  <- copy(target_new)
    target_new_b  <- copy(target_new)
    target_new_a$cell_name <- paste0(gsub("-", ".", target_new_a$bcr_patient_barcode), ".01A")
    target_new_b$cell_name <- paste0(gsub("-", ".", target_new_b$bcr_patient_barcode), ".01B")
    target_new    <- rbind(target_new_a, target_new_b)
    target_new    <- target_new[, c("drug_name", "cell_name", "Binary"), with=F]

    target_new$Binary <- ifelse(target_new$Binary == "Effective", 1, 0)
  }

  setnames(target_new, c("Compound", "cell_name", "target"))

  return(target_new)
}

Function_load_morgan_bits <- function(morgan="nci_60", radii_set, bit_set){
  if (morgan=="nci_60"){
    morgan_nci_bits   <- fread(paste0(in_morgan, "NCI_MORGAN_r_", radii_set ,"_b_", bit_set ,".txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==radii_set & bits==bit_set,]
    morgan_nci_bits$position <- paste0("mcf_", morgan_nci_bits$position)
  } else if (morgan=="tcga"){
    morgan_nci_bits   <- fread(paste0(tcga_objects, "TCGA_MORGAN_BITS_r_", radii_set ,"_b_", bit_set ,".txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==radii_set & bits==bit_set,]
    morgan_nci_bits$position <- paste0("mcf_", morgan_nci_bits$position)

  } else if (morgan=="tcga_multi"){
    morgan_nci_bits   <- fread(paste0(tcga_objects, "TCGA_MULTI_MORGAN_BITS_r_", radii_set ,"_b_", bit_set ,".txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==radii_set & bits==bit_set,]
    morgan_add        <- fread(paste0(tcga_objects, "TCGA_MORGAN_BITS_r_", radii_set ,"_b_", bit_set ,".txt"),
                           colClasses = c("numeric", "numeric", "character", "numeric", "numeric"))[radius==radii_set & bits==bit_set,]
    morgan_nci_bits   <- rbind(morgan_nci_bits, morgan_add)

    morgan_nci_bits$position <- paste0("mcf_", morgan_nci_bits$position)

  } else{
    morgan_nci_bits   <- c()
  }
  return(morgan_nci_bits)
}

Function_load_morgan_counts <- function(morgan="nci_60", radii_set){

  if (morgan=="nci_60"){
    morgan_nci_counts <- fread(paste0(in_morgan, "NCI_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "character", "numeric"))[radius==radii_set,]

  } else if (morgan=="tcga") {
    morgan_nci_counts <- fread(paste0(tcga_objects, "TCGA_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "character", "numeric"))[radius==radii_set,]

  } else if (morgan=="tcga_multi"){
    morgan_nci_counts <- fread(paste0(tcga_objects, "TCGA_MULTI_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "character", "numeric"))[radius==radii_set,]
    morgan_add        <- fread(paste0(tcga_objects, "TCGA_MORGAN_COUNTS.txt"),
                           colClasses = c("numeric", "character", "character", "numeric"))[radius==radii_set,]
    morgan_nci_counts <- rbind(morgan_nci_counts, morgan_add)

  } else{
    morgan_nci_counts <- c()
  }

  return(morgan_nci_counts)
}

Function_prep_morgan_bits <- function(input_morgan, type=""){

  if (type == "nci_60" | type == "tcga" | type == "tcga_multi"){
    input_morgan <- input_morgan[,c("bits", "radius", "Compound", "position", "value"), with=F]
    setnames(input_morgan, c("bits", "radius", "Compound", "bit_pos", "value"))

    #MODIFIED!!!
    input_morgan <- input_morgan[!Compound %in% c("Gemcitabine","Doxorubicin","Carboplatin"),]
    input_morgan <- rbind(input_morgan,
                          morgan_bits[Compound %in% c("Gemcitabine","Doxorubicin","Carboplatin"),])

  } else{
    input_morgan <- input_morgan
  }

  return(input_morgan)
}

Function_prep_morgan_counts <- function(input_morgan, type=""){

  if (type == "nci_60" | type == "tcga" | type == "tcga_multi"){
    setnames(input_morgan, c("radius", "Compound", "Substructure", "Counts"))

    input_morgan <- rbind(input_morgan, morgan_counts[,c("radius", "Compound", "Substructure", "Counts"),with=F]) #Global morgan_counts
    setkey(input_morgan)
    input_morgan <- unique(input_morgan)

  } else {
    input_morgan <- input_morgan
  }

  return(input_morgan)
}

Function_file_out <- function(in_name){

  file_name <- strsplit(in_name, "/")[[1]]
  file_name <- file_name[length(file_name)]
  file_name <- strsplit(file_name, "_train_drug")[[1]][1]

  return(file_name)
}

Function_target_morgan_bits_features_extracted_mf <- function(target_new, exp_table, morgan_table, target_cells=c(), target_bits=c(),
                                                              cgp_exp, genes=F, scaling=T, pca=F, common_genes=c(),
                                                              original_exp = c(), original_bits=c()){
  # target_bits and target_cells have to be in the order that there trained in the input model
  # exp_table has to be as genesxsamples (rows=genes, columns=samples)

  ########## EXPRESSION FEATURES ###########
  if (pca==F){
    colnames(cgp_exp) <- paste0("CGP_", colnames(cgp_exp))
    target_cells      <- paste0("CGP_", target_cells)
    cgp_cell_cor      <- cor(cgp_exp, method = "spearman") # Original cor without filtering cells

    # Build based on common ordered cell features
    common_genes   <- intersect(rownames(cgp_exp), rownames(exp_table))
    cgp_exp        <- cgp_exp[common_genes, target_cells]
    exp_table      <- exp_table[common_genes,]
    target_samples <- colnames(exp_table)

    if (genes==F){
      cell_feat <- cbind(exp_table, cgp_exp)
      cell_feat <- cor(cell_feat, method = "spearman")
      cell_feat <- cell_feat[target_samples, target_cells]
      cell_feat <- scale(cell_feat) # HEAVILY MODIFIED!!!! MODIFIED!!!!!!!

      cell_feat <- data.table(cell_feat, keep.rownames = T)
    } else {
      exp_table <- t(log(exp_table))
      exp_table <- exp_table[,target_cells] #target_genes in this case

      cell_feat <- data.table(exp_table, keep.rownames = T)
    }
  } else {
    if (genes == F){
      print("Cor exp PCA pre-scaling")

      # Modify names
      colnames(cgp_exp) <- paste0("CGP_", colnames(cgp_exp))

      # Get original principal rotation from cgp_exp
      cgp_cor        <- cor(cgp_exp, method = "spearman")
      cgp_pca_exp    <- prcomp(cgp_cor, center = T, scale. = T)

      # Prepare target expression to apply PCA to it
      common_genes   <- intersect(rownames(cgp_exp), rownames(exp_table))
      feat_cells     <- colnames(cgp_exp) # NEED TO RENAME ORIGINAL
      target_samples <- colnames(exp_table) # NEED TO RENAME ORIGINAL IN CASE THEY HAVE SAME CELLS, THEN IT DOESN'T KNOW WHAT TO PICK
      cgp_exp        <- cgp_exp[common_genes, ]
      exp_table      <- exp_table[common_genes,]

      cell_feat      <- cbind(exp_table, cgp_exp)
      cell_feat      <- cor(cell_feat, method = "spearman")
      cell_feat      <- cell_feat[target_samples, feat_cells]
      cell_feat      <- scale(cell_feat) # COULD HEAVILY MODIFY HERE BEFORE APPLYING TRAIN SCALING

      # Apply scaling to cell_feat prior to rotation by PCA
      cgp_cor        <- scale(cgp_cor) # To obtain original pre-PCA scaling attributes
      cgp_cell_mean  <- attributes(cgp_cor)$`scaled:center`
      cgp_cell_sd    <- attributes(cgp_cor)$`scaled:scale`

      cell_feat_sc   <- sweep(cell_feat, 2, cgp_cell_mean, "-")
      cell_feat_sc   <- sweep(cell_feat_sc, 2, cgp_cell_sd, "/")

      # Apply rotation
      cell_feat      <- cell_feat_sc %*% cgp_pca_exp$rotation[,target_cells]
      cell_feat      <- data.table(cell_feat, keep.rownames = T)

    } else {
      #FIX!!! ORIGINAL EXPRESSION SHOULD HAVE ORIGINAL GENES TO OBTAIN ORIGINAL PCA ROTATION
      original_exp <- original_exp[common_genes, ]
      exp_table    <- exp_table[common_genes, ]
      cgp_pca_exp  <- prcomp(t(original_exp), center = T, scale. = T) # PCA on original to obtain rotation
      cell_feat    <- scale(t(exp_table)) %*% cgp_pca_exp$rotation[,target_cells] # Apply rotation from original data
      cell_feat    <- data.table(cell_feat, keep.rownames = T)
    }
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  ########## MORGAN FEATURES ###########
  if (length(target_bits) > 0){
    if (pca == F){
      drug_feat     <- morgan_table[bit_pos %in% target_bits,]
      drug_feat     <- acast(drug_feat, Compound~bit_pos, value.var="value")
      drug_feat     <- data.table(drug_feat[, target_bits], keep.rownames = T)

    } else {
      print("Cor morgan PCA")

      # Get original PCA rotation (since it is binary, there is no pre-scaling)
      original_bits <- acast(original_bits, Compound~bit_pos, value.var="value")
      original_pos  <- colnames(original_bits)
      drug_feat     <- acast(morgan_table, Compound~bit_pos, value.var="value")
      drug_feat     <- drug_feat[,original_pos]

      cgp_pca_bit   <- prcomp(original_bits, center=F, scale. = F)

      # Apply PCA rotation, no need to pre-scaling since original morgan bits were binary and unscaled
      drug_feat     <- drug_feat %*% cgp_pca_bit$rotation[,target_bits]
      drug_feat     <- data.table(drug_feat, keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))

  } else {
    drug_feat <- data.table()
  }

  ##################################################################################################
  # IS SCALING NECESSARY ?
  if (scaling==T){

    if (pca==T){
      if (genes==F){
        print("Cor exp PCA post-scaling")
        cgp_cell_feat <- data.table(cgp_pca_exp$x[, target_cells], keep.rownames = T)
      }

      if (length(target_bits) > 0){
        print("Cor morgan PCA post-scaling")
        cgp_drug_feat <- data.table(cgp_pca_bit$x[, target_bits], keep.rownames = T)

        # Do morgan scaling within loop since it is only applicable when PCA==T
        setnames(cgp_drug_feat, c("Compound", colnames(cgp_drug_feat)[2:ncol(cgp_drug_feat)]))

        not_scaling <- 1
        scaling     <- 2:ncol(cgp_drug_feat)

        scale_train_drug  <- scale(cgp_drug_feat[, scaling, with=F])

        train_drug_mean   <- attributes(scale_train_drug)$`scaled:center`
        train_drug_sd     <- attributes(scale_train_drug)$`scaled:scale`

        drug_feat_scaled  <- sweep(drug_feat[, scaling, with=F], 2, train_drug_mean, "-")
        drug_feat_scaled  <- sweep(drug_feat_scaled, 2, train_drug_sd, "/")
        drug_feat         <- cbind(drug_feat[, not_scaling, with=F],
                             drug_feat_scaled)
        print("morgan scaled")
      }
    } else {
      cgp_cell_feat <- data.table(cgp_cell_cor[, target_cells], keep.rownames = T)
    }

    setnames(cgp_cell_feat, c("cell_name", colnames(cgp_cell_feat)[2:ncol(cgp_cell_feat)]))

    not_scaling <- 1
    scaling     <- 2:ncol(cgp_cell_feat)

    scale_train_cell  <- scale(cgp_cell_feat[, scaling, with=F])

    train_cell_mean   <- attributes(scale_train_cell)$`scaled:center`
    train_cell_sd     <- attributes(scale_train_cell)$`scaled:scale`

    cell_feat_scaled  <- sweep(cell_feat[, scaling, with=F], 2, train_cell_mean, "-")
    cell_feat_scaled  <- sweep(cell_feat_scaled, 2, train_cell_sd, "/")
    cell_feat         <- cbind(cell_feat[, not_scaling, with=F],
                         cell_feat_scaled)

    print("cell scaled")
  }
  ##################################################################################################

  # Merge to obtain total number of possible combinations
  if (length(target_bits) > 0){
    feat_table <- merge(target_new, drug_feat[, 1:2, with=F], by="Compound", allow.cartesian=TRUE)

  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name", allow.cartesian=TRUE)

  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Extract indices for combinations
  if (length(target_bits) > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else {
    drug_index <- c()
  }

  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$target

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))
}

Function_target_morgan_counts_features_extracted_mf <- function(target_new, exp_table, morgan_table, target_cells=c(), target_counts=c(),
                                                                cgp_exp, genes=F, scaling=T){
  # target_counts and target_cells have to be in the order that there trained in the input model
  # exp_table has to be as genesxsamples (rows=genes, columns=samples)

  # Minor fix to account for common cells
  colnames(cgp_exp) <- paste0("CGP_", colnames(cgp_exp))
  target_cells      <- paste0("CGP_", target_cells)
  cgp_cell_cor      <- cor(cgp_exp, method = "spearman")

  # Build based on common ordered cell features
  common_genes   <- intersect(rownames(cgp_exp), rownames(exp_table))
  cgp_exp        <- cgp_exp[common_genes, target_cells]
  exp_table      <- exp_table[common_genes,]
  target_samples <- colnames(exp_table)

  if (genes==F){
    cell_feat <- cbind(exp_table, cgp_exp)
    cell_feat <- cor(cell_feat, method = "spearman")
    cell_feat <- cell_feat[target_samples, target_cells]

    cell_feat <- data.table(cell_feat, keep.rownames = T)
  } else {
    exp_table <- t(log(exp_table))
    exp_table <- exp_table[,target_cells] #target_genes in this case

    cell_feat <- data.table(exp_table, keep.rownames = T)
  }
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  ##################################################################################################
  # IS SCALING NECESSARY
  if (scaling==T){

    cgp_cell_feat <- data.table(cgp_cell_cor[, target_cells], keep.rownames = T)
    setnames(cgp_cell_feat, c("cell_name", colnames(cgp_cell_feat)[2:ncol(cgp_cell_feat)]))

    not_scaling <- 1
    scaling     <- 2:ncol(cgp_cell_feat)

    scale_train_cell  <- scale(cgp_cell_feat[, scaling, with=F])

    train_cell_mean   <- attributes(scale_train_cell)$`scaled:center`
    train_cell_sd     <- attributes(scale_train_cell)$`scaled:scale`

    cell_feat_scaled  <- sweep(cell_feat[, scaling, with=F], 2, train_cell_mean, "-")
    cell_feat_scaled  <- sweep(cell_feat_scaled, 2, train_cell_sd, "/")
    cell_feat         <- cbind(cell_feat[, not_scaling, with=F],
                         cell_feat_scaled)
    print("scaled")
  }
  ##################################################################################################

  # Build based on common ordered drug features
  if (length(target_counts) > 0 ){
    morgan_table$Substructure <- paste0("mc_", morgan_table$Substructure)
    target_counts <- paste0("mc_", target_counts)
    drug_feat     <- morgan_table[Substructure %in% target_counts, ]
    saveRDS(drug_feat, "x.rds")
    drug_feat     <- acast(drug_feat, Compound~Substructure, value.var="Counts", fill=0)

    drug_feat     <- data.table(drug_feat[, target_counts], keep.rownames = T)

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  } else {
    drug_feat  <- data.table()
  }

  # Merge to obtain total number of possible combinations
  if (length(target_counts) > 0){
    feat_table <- merge(target_new, drug_feat[, 1:2, with=F], by="Compound")
  }

  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name")
  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Extract indices for combinations
  if (length(target_counts) > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else{
    drug_index <- c()
  }

  cell_index <- sapply(feat_table$cell_name, function(x)  which(x==cell_feat$cell_name))
  drug_index <- unlist(drug_index)
  cell_index <- unlist(cell_index)
  target     <- feat_table$target

  # Return feature tables and indices
  return(list(drug_feat = drug_feat,
              cell_feat = cell_feat,
              drug_index = drug_index,
              cell_index = cell_index,
              target = target,
              feat_table = feat_table))

}

############################################################################################################
# LOAD INPUT
args      <- commandArgs(trailingOnly = TRUE)

target      <- args[1]
met_type    <- args[2]
cgp_drug    <- args[3] #_train_drug
cgp_cell    <- args[4] #_train_cell
class_mlp   <- as.logical(args[5])
batch_norm  <- args[6]
bn_external <- args[7]
pca         <- as.logical(args[8])
radii_set   <- as.numeric(args[9])
bit_set     <- as.numeric(args[10])
out_folder  <- "/tigress/zamalloa/PREDICTIONS/"
out_file    <- paste0(out_folder, Function_file_out(cgp_drug), "_bn_external_", bn_external)
print(out_file)

############################################################################################################
# LOAD DATA

#GENERAL
in_folder         <- "/tigress/zamalloa/CGP_FILES/"
in_morgan         <- "/tigress/zamalloa/MORGAN_FILES/"
in_objects        <- "/tigress/zamalloa/OBJECTS/"
tcga_objects      <- "/tigress/zamalloa/TCGA_FILES/"
fire_rnaseq       <- "/tigress/zamalloa/OBJECTS/FIREHOSE/RNASEQ/"
cgp_new           <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
nci_to_cgp_name   <- readRDS(paste0(in_objects, "080716.nci60_names_to_cgp.rds"))
tcga_ding         <- readRDS(paste0(in_objects, "121316_ding_tcga.rds"))
tcga_multi_ding   <- readRDS(paste0(in_objects, "020717_ding_multi.rds"))
morgan_bits       <- fread(paste0(in_morgan, "CGP_MORGAN_BITS.txt"))[radius==radii_set & bits==bit_set,]
morgan_counts     <- fread(paste0(in_morgan, "CGP_MORGAN_COUNTS.txt"),
                       colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,] #TO BE DECIDED

if (batch_norm=="cgp_nci60"){

 print("batch_norm nci60")
 cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_nci60_exp_b_norm.rds"))[["EXP.1"]]
 nci_exp         <- readRDS(paste0(in_objects, "121216_cgp_nci60_exp_b_norm.rds"))[["EXP.2"]]

} else if (batch_norm=="cgp_ccle"){

 print("batch_norm ccle")
 cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_ccle_exp_b_norm.rds"))[["EXP.1"]]
 ccle_exp        <- readRDS(paste0(in_objects, "121216_cgp_ccle_exp_b_norm.rds"))[["EXP.2"]]

} else if (batch_norm=="tcga_brca"){

 print("batch_norm tcga_brca")
 cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.1"]]
 tcga_brca_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.2"]]

} else if (batch_norm=="tcga_coad"){

 print("batch_norm tcga_coad")
 cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.1"]]
 tcga_coad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.2"]]

} else if (batch_norm=="tcga_stad"){

 print("batch_norm tcga_stad")
 cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.1"]]
 tcga_stad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.2"]]

} else if (batch_norm=="tcga_luad"){

 print("batch_norm tcga_luad")
 cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.1"]]
 tcga_luad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.2"]]

} else{
 print("None")
 cgp_exp         <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))
 ccle_exp        <- readRDS(paste0(in_folder, "121116_ccle_exp.rds"))
 nci_exp         <- readRDS(paste0(in_objects, "121216_nci60_exp.rds"))
 tcga_brca_exp   <- readRDS(paste0(tcga_objects, "041715.BRCA.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]
 tcga_luad_exp   <- readRDS(paste0(tcga_objects, "081915.LUAD.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]
 tcga_ucec_exp   <- readRDS(paste0(tcga_objects, "081915.UCEC.GA.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]
 tcga_stad_exp   <- readRDS(paste0(tcga_objects, "121416.STAD.RNASEQ.MATRICES.rds"))[["tumor"]]
}

#CCLE
ccle_new    <- readRDS(paste0(in_objects, "121116_ccle.rds"))

#NCI-60
nci_new     <- readRDS(paste0(in_objects, "120916_nci60_new.rds"))

#TCGA

#MODEL FEATURES
cell_features <- fread(cgp_cell, header=T, sep="\t")
cell_features <- colnames(cell_features)[2:ncol(cell_features)]

drug_features <- fread(cgp_drug, header=T, sep="\t")
drug_features <- colnames(drug_features)[2:ncol(drug_features)]

common_genes  <- intersect(intersect(rownames(cgp_exp), rownames(nci_exp)), rownames(ccle_exp))
############################################################################################################
# EXECUTE
#cgp_new <- cgp_new[Compound %in% c("GDC0449", "MS-275", "PAC-1", "RDEA119", "TG101348"),] #NOTE (Used during testing)

# Build feat tables
if (met_type == "morgan_bits"){

  if (target=="nci_60"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_nci60_exp_b_norm.rds"))[["EXP.1"]]
      nci_exp         <- readRDS(paste0(in_objects, "121216_cgp_nci60_exp_b_norm.rds"))[["EXP.2"]]
    }

    colnames(nci_exp) <- sapply(colnames(nci_exp), function(x)  unique(nci_new[,c("cell_name", "CELL"),with=F])[CELL==x,]$cell_name )

    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(nci_new, type="nci_60", class = class_mlp),
                                                                           nci_exp,
                                                                           Function_prep_morgan_bits(Function_load_morgan_bits("nci_60",
                                                                                                    radii_set, bit_set), type="nci_60"),
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=T,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="ccle"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(in_objects, "121216_cgp_ccle_exp_b_norm.rds"))[["EXP.1"]]
      ccle_exp        <- readRDS(paste0(in_objects, "121216_cgp_ccle_exp_b_norm.rds"))[["EXP.2"]]
    }
    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(ccle_new, type="ccle", class = class_mlp),
                                                                           ccle_exp,
                                                                           morgan_bits,
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=T,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="tcga_brca"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.1"]]
      tcga_brca_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.2"]]
    }

    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="BRCA",], type="tcga", class = class_mlp),
                                                                           tcga_brca_exp,
                                                                           morgan_bits,
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=F,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="tcga_coad"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.1"]]
      tcga_coad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.2"]]
    }
    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="COAD",], type="tcga", class = class_mlp),
                                                                           tcga_coad_exp,
                                                                           Function_prep_morgan_bits(Function_load_morgan_bits("tcga"), type="tcga",
                                                                                                     radii_set, bit_set),
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=F,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="tcga_stad"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.1"]]
      tcga_stad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.2"]]
    }
    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="STAD",], type="tcga", class = class_mlp),
                                                                           tcga_stad_exp,
                                                                           Function_prep_morgan_bits(Function_load_morgan_bits("tcga"), type="tcga",
                                                                                                     radii_set, bit_set),
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=F,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="tcga_luad"){
    if (as.logical(bn_external)==T){
      print("bn_external")
      cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.1"]]
      tcga_luad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.2"]]
    }
    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="LUAD",], type="tcga", class = class_mlp),
                                                                           tcga_luad_exp,
                                                                           Function_prep_morgan_bits(Function_load_morgan_bits("tcga"), type="tcga",
                                                                                                     radii_set, bit_set),
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=F,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)
  } else if (target=="self"){
    feat_table <- Function_target_morgan_bits_features_extracted_mf(Function_prep_new(cgp_new, type="cgp", class = class_mlp),
                                                                           cgp_exp,
                                                                           morgan_bits,
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=T,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)

  } else if (target=="tcga_all" | target=="tcga_multi") {

    if (target=="tcga_all"){
      tcga_target    <- Function_prep_new(tcga_ding, type=target)
      morgan_type    <- "tcga"
    } else {
      tcga_target    <- Function_prep_new(tcga_multi_ding, type=target)
      morgan_type    <- "tcga_multi"
    }

    cancer_types   <- c('BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC',
                      'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO',
                      'OV','PAAD','PCPG','PRAD','READ','SARC','SKCM',
                      'STAD','TGCT','THCA','UCEC','UCS')
    tcga_exp       <- lapply(cancer_types, function(x) {

      temp_exp <- readRDS(paste0(fire_rnaseq, x, ".rds"))[["tumor"]]
      return(temp_exp)
      })

    common_genes   <- Reduce(intersect, lapply(tcga_exp, function(x) rownames(x)))
    tcga_exp       <- do.call(cbind,    lapply(tcga_exp, function(x) x[common_genes,]))

    common_samples <- intersect(colnames(tcga_exp), tcga_target$cell_name)
    tcga_exp       <- tcga_exp[,common_samples]

    feat_table <- Function_target_morgan_bits_features_extracted_mf(tcga_target,
                                                                   tcga_exp,
                                                                   Function_prep_morgan_bits(Function_load_morgan_bits(morgan_type, radii_set, bit_set),
                                                                                              type=morgan_type),
                                                                   target_cells = cell_features, target_bits = drug_features,
                                                                   cgp_exp, genes=F, scaling=T,
                                                                   pca = pca, common_genes = common_genes,
                                                                   original_exp = cgp_exp, original_bits = morgan_bits)

  } else {
    # Assumes that in this scenario we have a specific compound belonging to a particular dataset

    drug_set    <- strsplit(target, "_")[[1]][1]
    drug        <- strsplit(target, "_")[[1]][2]

    set_type    <- ifelse(drug_set == "nci", "nci_60", drug_set)
    if (drug_set == "ccle") {
      morgan_bits <- morgan_bits
    } else{
      morgan_bits <- Function_prep_morgan_bits(Function_load_morgan_bits(set_type,
                               radii_set, bit_set), type=set_type)
    }

    target_new  <- Function_prep_new(get(paste0(drug_set,"_new")), type=set_type, class = class_mlp)

    feat_table  <- Function_target_morgan_bits_features_extracted_mf(target_new,
                                                                           get(paste0(drug_set,"_exp")),
                                                                           morgan_bits,
                                                                           target_cells = cell_features, target_bits = drug_features,
                                                                           cgp_exp, genes=F, scaling=T,
                                                                           pca = pca, common_genes = common_genes,
                                                                           original_exp = cgp_exp, original_bits = morgan_bits)

  }

} else if (met_type == "morgan_counts"){

  if (target=="nci_60"){
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(nci_new, type="nci_60", class = class_mlp),
                                                                             nci_exp,
                                                                             Function_prep_morgan_counts(Function_load_morgan_counts("nci_60"), type="nci_60",
                                                                                                         radii_set),
                                                                             target_cells = cell_features, target_counts = drug_features,
                                                                             cgp_exp, genes=F, scaling=T)
  } else if (target=="ccle"){
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(ccle_new, type="ccle", class = class_mlp),
                                                                             ccle_exp,
                                                                             morgan_counts,
                                                                             target_cells = cell_features, target_counts = drug_features,
                                                                             cgp_exp, genes=F, scaling=T)
  } else if (target=="tcga_brca"){
    cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.1"]]
    tcga_brca_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_brca_exp_b_norm.rds"))[["EXP.2"]]
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="BRCA",], type="tcga", class = class_mlp),
                                                                           tcga_brca_exp,
                                                                           morgan_counts,
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)
  } else if (target=="tcga_coad"){
    cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.1"]]
    tcga_coad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_coad_exp_b_norm.rds"))[["EXP.2"]]
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="COAD",], type="tcga", class = class_mlp),
                                                                           tcga_coad_exp,
                                                                           Function_prep_morgan_counts(Function_load_morgan_counts("tcga"), type="tcga",
                                                                                                       radii_set),
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)
  } else if (target=="tcga_stad"){
    cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.1"]]
    tcga_stad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_stad_exp_b_norm.rds"))[["EXP.2"]]
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="STAD",], type="tcga", class = class_mlp),
                                                                           tcga_stad_exp,
                                                                           Function_prep_morgan_counts(Function_load_morgan_counts("tcga"), type="tcga",
                                                                                                       radii_set),
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)
  } else if (target=="tcga_luad"){
    cgp_exp         <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.1"]]
    tcga_luad_exp   <- readRDS(paste0(tcga_objects, "121216_cgp_tcga_luad_exp_b_norm.rds"))[["EXP.2"]]
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(tcga_ding[Cancer=="LUAD",], type="tcga", class = class_mlp),
                                                                           tcga_luad_exp,
                                                                           Function_prep_morgan_counts(Function_load_morgan_counts("tcga"), type="tcga",
                                                                                                       radii_set),
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)
  } else if (target=="self"){
    feat_table <- Function_target_morgan_counts_features_extracted_mf(Function_prep_new(cgp_new, type="cgp", class = class_mlp),
                                                                           cgp_exp,
                                                                           morgan_counts,
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)
  }
}

# Prep to write out for prediction
index_table<- data.table(drug   = feat_table$drug_index - 1,
                         cell   = feat_table$cell_index - 1,
                         target = feat_table$target)

############################################################################################################
# WRITE
write.table(feat_table$drug_feat, paste0(out_file, "_", target, "_test_drug"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(feat_table$cell_feat, paste0(out_file, "_", target, "_test_cell"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(index_table, paste0(out_file, "_", target, "_test_index"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(feat_table$feat_table, paste0(out_file, "_", target, "_feat_table"),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing tables")
