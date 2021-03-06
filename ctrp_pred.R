# ctrp_pred.R
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
        target_new <- target_new[,c("Compound", "cell_name", "NORM.ACT"), with=F]
      } else {
        target_new <- target_new[,c("Compound", "cell_name", "ActArea"), with=F] #MODIFIED
      }
    }

  } else if (type=="cgp"){
    if (class == T){
      target_new <- target_new[,c("Compound", "cell_name", "pic50_class"),with=F]
      target_new <- target_new[pic50_class!=2,]
      cgp_filter <- target_new[,list(MEAN = mean(pic50_class)), by="Compound"]
      cgp_filter <- cgp_filter[MEAN > 0.35 & MEAN < 0.65,]$Compound  # NECESSARY FILTER FOR OUT OF SET TESTING
      target_new <- target_new[Compound %in% cgp_filter,]

    } else {
      if (scaled == T){
        target_new <- target_new[, c("Compound", "cell_name", "NORM.AUC"), with=F]
      } else {
        target_new <- target_new[, c("Compound", "cell_name", "AUC"), with=F]
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
      tcga_compounds <- tcga_compounds[samples>30,]$drug_name
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

    #MODIFIED!!! NOTE:cgp_bits needed for this (name match!)
    # input_morgan <- input_morgan[!Compound %in% c("Gemcitabine","Doxorubicin","Carboplatin"),]
    # input_morgan <- rbind(input_morgan,
    #                       ctrp_bits[Compound %in% c("Gemcitabine","Doxorubicin","Carboplatin"),])

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
                                                              genes=F, scaling=T, genespca=F, drugspca=F,
                                                              original_exp = c(), original_bits=c()){
  # target_bits and target_cells have to be in the order that there trained in the input model
  # exp_table has to be as genesxsamples (rows=genes, columns=samples)

  ########## EXPRESSION FEATURES ###########
  if (genespca==F){
    print("genespca==F")

    # Build based on common ordered cell features
    common_genes   <- intersect(rownames(original_exp), rownames(exp_table))
    exp_table      <- exp_table[common_genes,]
    target_samples <- colnames(exp_table) # Samples we want to predict
    colnames(original_exp) <- paste0("CTRP_", colnames(original_exp)) # Renaming cell names in original expression
    target_cells           <- paste0("CTRP_", target_cells) # Since name changed

    if (genes==F){
      print("genes==F")
      ctrp_cell_cor        <- cor(original_exp, method = "spearman") # Original cor without filtering cells

      cell_feat <- cbind(exp_table, original_exp[common_genes, target_cells])
      cell_feat <- cor(cell_feat, method = "spearman")
      cell_feat <- cell_feat[target_samples, target_cells]

      cell_feat <- data.table(cell_feat, keep.rownames = T)

    } else {
      print("genes==T")
      exp_table <- t(exp_table)
      exp_table <- exp_table[,target_cells] #target_genes in this case

      cell_feat <- data.table(exp_table, keep.rownames = T)
    }
  } else {
    print("genespca==T")
    if (genes == F){
      print("genes==F")

      # Modify names in case target exp has same cells
      colnames(original_exp) <- paste0("CTRP_", colnames(original_exp))

      # Get original principal rotation from cgp_exp
      ctrp_cell_cor  <- cor(original_exp, method = "spearman")
      ctrp_pca_exp   <- prcomp(ctrp_cell_cor, center = T, scale. = T)

      # Prepare target expression to apply PCA to it
      common_genes   <- intersect(rownames(original_exp), rownames(exp_table))
      feat_cells     <- colnames(original_exp)
      target_samples <- colnames(exp_table)
      exp_table      <- exp_table[common_genes,]

      cell_feat      <- cbind(exp_table, original_exp[common_genes, feat_cells])
      cell_feat      <- cor(cell_feat, method = "spearman")
      cell_feat      <- cell_feat[target_samples, feat_cells] #Complete feature matrix that needs to be scaled

      # Apply scaling to cell_feat prior to rotation by PCA
      ctrp_cell_cor_scale <- scale(ctrp_cell_cor) # To obtain original pre-PCA scaling attributes
      ctrp_cell_mean <- attributes(ctrp_cell_cor_scale)$`scaled:center`
      ctrp_cell_sd   <- attributes(ctrp_cell_cor_scale)$`scaled:scale`

      cell_feat_sc   <- sweep(cell_feat, 2, ctrp_cell_mean, "-")
      cell_feat_sc   <- sweep(cell_feat_sc, 2, ctrp_cell_sd, "/")

      # Apply rotation
      cell_feat      <- cell_feat_sc %*% ctrp_pca_exp$rotation[,target_cells]
      cell_feat      <- data.table(cell_feat, keep.rownames = T)

    } else {
      print("genes==T")
      #NOTE If PCA==T and genes==T, then it is assumed that expression matrices contain equal ordered sets of genes

      ctrp_pca_exp   <- prcomp(t(original_exp), center = T, scale. = T) # PCA on original to obtain rotation

      original_scale <- scale(t(original_exp))
      original_mean  <- attributes(original_scale)$`scaled:center`
      original_sd    <- attributes(original_scale)$`scaled:scale`

      cell_feat_sc   <- t(exp_table)
      cell_feat_sc   <- sweep(cell_feat_sc, 2, original_mean, "-")
      cell_feat_sc   <- sweep(cell_feat_sc, 2, original_sd, "/")

      cell_feat      <- cell_feat_sc %*% ctrp_pca_exp$rotation[,target_cells] # Apply rotation from original data
      cell_feat      <- data.table(cell_feat, keep.rownames = T)
    }
  }

  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))

  ########## MORGAN FEATURES ###########
  if (length(target_bits) > 0){
    if (drugspca == F){
      print("drugspca==F")
      drug_feat     <- morgan_table[bit_pos %in% target_bits,]
      drug_feat     <- acast(drug_feat, Compound~bit_pos, value.var="value")
      drug_feat     <- data.table(drug_feat[, target_bits], keep.rownames = T)

    } else {
      print("Morgan PCA")

      # Get original PCA rotation (since it is binary, there is no pre-scaling)
      original_bits <- acast(original_bits, Compound~bit_pos, value.var="value")
      original_pos  <- colnames(original_bits)
      drug_feat     <- acast(morgan_table, Compound~bit_pos, value.var="value")
      drug_feat     <- drug_feat[,original_pos]

      ctrp_pca_bit  <- prcomp(original_bits, center=F, scale. = F) #No PCA scaling with binary variables for original bits

      # Apply PCA rotation, no need to pre-scaling since original morgan bits were binary and unscaled
      drug_feat     <- drug_feat %*% ctrp_pca_bit$rotation[,target_bits]
      drug_feat     <- data.table(drug_feat, keep.rownames = T)
    }

    setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))

  } else {
    drug_feat <- data.table()
  }

  ##################################################################################################
  # IS SCALING NECESSARY? NOTE: ONLY SCALING IF PCA==F!!!
  # Since we don't scale binary/counts, then drug features are not scaled
  if (scaling==T){
    print("scaling==T")

    if (genespca==F){
      if (genes==F){
        ctrp_cell_feat <- data.table(ctrp_cell_cor[, target_cells], keep.rownames = T)
      } else {
        ctrp_cell_feat  <- t(original_exp)
        ctrp_cell_feat  <- data.table(ctrp_cell_feat[, target_cells], keep.rownames =T) # target_cells==target_genes
      }
      setnames(ctrp_cell_feat, c("cell_name", colnames(ctrp_cell_feat)[2:ncol(ctrp_cell_feat)]))

      not_scaling <- 1
      scaling     <- 2:ncol(ctrp_cell_feat)

      scale_train_cell  <- scale(ctrp_cell_feat[, scaling, with=F])

      train_cell_mean   <- attributes(scale_train_cell)$`scaled:center`
      train_cell_sd     <- attributes(scale_train_cell)$`scaled:scale`

      cell_feat_scaled  <- sweep(cell_feat[, scaling, with=F], 2, train_cell_mean, "-")
      cell_feat_scaled  <- sweep(cell_feat_scaled, 2, train_cell_sd, "/")
      cell_feat         <- cbind(cell_feat[, not_scaling, with=F],
                           cell_feat_scaled)

      print("cell scaled")
    }
  }
  ##################################################################################################

  # Merge to obtain total number of possible combinations
  if (length(target_bits) > 0){
    feat_table <- merge(target_new, drug_feat[, 1:2, with=F], by="Compound", allow.cartesian=TRUE)

  } else {
    feat_table <- target_new
  }
  feat_table <- merge(feat_table, cell_feat[, 1:2, with=F], by="cell_name", allow.cartesian=TRUE)

  setkey(feat_table)
  feat_table <- unique(feat_table)

  # Extract indices for combinations
  if (length(target_bits) > 0){
    drug_index <- sapply(feat_table$Compound,  function(x)  which(x==drug_feat$Compound))
  } else {
    drug_index <- 0
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
genespca    <- as.logical(args[7])
drugspca    <- as.logical(args[8])
radii_set   <- as.numeric(args[9])
bit_set     <- as.numeric(args[10])
genes       <- as.logical(args[11])
out_folder  <- "/tigress/zamalloa/PREDICTIONS/"
out_file    <- paste0(out_folder, Function_file_out(cgp_drug))
print(out_file)

############################################################################################################
# LOAD DATA

#GENERAL
in_folder         <- "/tigress/zamalloa/CGP_FILES/"
in_morgan         <- "/tigress/zamalloa/MORGAN_FILES/"
in_objects        <- "/tigress/zamalloa/OBJECTS/"
tcga_objects      <- "/tigress/zamalloa/TCGA_FILES/"
fire_rnaseq       <- "/tigress/zamalloa/OBJECTS/FIREHOSE/RNASEQ/"
gee_folder        <- "/tigress/zamalloa/OBJECTS/GEELEHER/"

# Load original data
ctrp_exp          <- readRDS(paste0(in_objects, "ctrp_exp_2.1.rds"))

ctrp_bits         <- fread(paste0(in_morgan, "CTRP_MORGAN_BITS_r_",radii_set,"_b_",bit_set,".txt"))
setnames(ctrp_bits, c("radius", "bits", "Compound", "bit_pos", "value"))
ctrp_bits$bit_pos <- paste0("mcf_", ctrp_bits$bit_pos)

ctrp_counts       <- fread(paste0(in_morgan, "CTRP_MORGAN_COUNTS.txt"),
                       colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,] #12 normal setting

# Choose target
if (target=="cgp"){
  feat_table    <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
  feat_table    <- Function_prep_new(feat_table, type="cgp", class=class_mlp, scaled=F) #NOTE: Manually inputted the class
  exp_table     <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))
  # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_cgp_exp_norm.rds"))[["EXP.2"]]
  morgan_bits   <- fread(paste0(in_morgan, "CGP_MORGAN_BITS.txt"))[radius==radii_set & bits==bit_set,]
  morgan_counts <- fread(paste0(in_morgan, "CGP_MORGAN_COUNTS.txt"),
                         colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,]

} else if (target=="ccle"){
  ccle_new      <- readRDS(paste0(in_objects, "121116_ccle.rds"))
  feat_table    <- Function_prep_new(ccle_new, type="ccle", class = class_mlp, scaled=F) #NOTE: Manually inputted the class
  exp_table     <- readRDS(paste0(in_folder, "121116_ccle_exp.rds"))
  # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_ccle_exp_norm.rds"))[["EXP.2"]]
  morgan_bits   <- fread(paste0(in_morgan, "CGP_MORGAN_BITS.txt"))[radius==radii_set & bits==bit_set,]
  morgan_counts <- fread(paste0(in_morgan, "CGP_MORGAN_COUNTS.txt"),
                         colClasses = c("numeric", "character", "numeric", "numeric"))[radius==radii_set,]

} else if (target=="nci_60"){
  nci_new             <- readRDS(paste0(in_objects, "120916_nci60_new.rds"))
  feat_table          <- Function_prep_new(nci_new, type="nci_60", class = class_mlp, scaled=F)
  exp_table           <- readRDS(paste0(in_objects, "121216_nci60_exp.rds"))
  colnames(exp_table) <- sapply(colnames(exp_table), function(x)  unique(nci_new[,c("cell_name", "CELL"),with=F])[CELL==x,]$cell_name )
  nci_to_cgp_name     <- readRDS(paste0(in_objects, "080716.nci60_names_to_cgp.rds"))
  morgan_bits         <- Function_prep_morgan_bits(Function_load_morgan_bits("nci_60", radii_set, bit_set), type="nci_60")

} else if (grepl("tcga", target)==T){
  tcga_ding       <- readRDS(paste0(in_objects, "121316_ding_tcga.rds"))
  tcga_multi_ding <- readRDS(paste0(in_objects, "020717_ding_multi.rds"))

  if(target=="tcga_brca"){
    feat_table    <- Function_prep_new(tcga_ding[Cancer==toupper(gsub("tcga_", "", target)),], type="tcga", class = class_mlp)
    exp_table     <- readRDS(paste0(tcga_objects, "041715.BRCA.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]

  } else if (target=="tcga_luad"){
    feat_table    <- Function_prep_new(tcga_ding[Cancer==toupper(gsub("tcga_", "", target)),], type="tcga", class = class_mlp)
    exp_table     <- readRDS(paste0(tcga_objects, "081915.LUAD.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]

  } else if (target=="tcga_ucec"){
    feat_table    <- Function_prep_new(tcga_ding[Cancer==toupper(gsub("tcga_", "", target)),], type="tcga", class = class_mlp)
    exp_table     <- readRDS(paste0(tcga_objects, "081915.UCEC.GA.RNASEQ.MATRICES.V2.RSEM.UQ.rds"))[["tumor"]]

  } else if (target=="tcga_stad"){
    feat_table    <- Function_prep_new(tcga_ding[Cancer==toupper(gsub("tcga_", "", target)),], type="tcga", class = class_mlp)
    exp_table     <- readRDS(paste0(tcga_objects, "121416.STAD.RNASEQ.MATRICES.rds"))[["tumor"]]

  } else if (target=="tcga_all" | target=="tcga_multi") {

    if (target=="tcga_all"){
      feat_table     <- Function_prep_new(tcga_ding, type=target)
    } else {
      feat_table     <- Function_prep_new(tcga_multi_ding, type=target)
    }

    cancer_types   <- c('BLCA','BRCA','CESC','COAD','ESCA','GBM','HNSC',
                      'KIRC','KIRP','LGG','LIHC','LUAD','LUSC','MESO',
                      'OV','PAAD','PCPG','PRAD','READ','SARC','SKCM',
                      'STAD','TGCT','THCA','UCEC','UCS')
    exp_table      <- lapply(cancer_types, function(x) {

      temp_exp <- readRDS(paste0(fire_rnaseq, x, ".rds"))[["tumor"]]
      return(temp_exp)
      })

    common_genes   <- Reduce(intersect, lapply(exp_table, function(x) rownames(x)))
    exp_table      <- do.call(cbind,    lapply(exp_table, function(x) x[common_genes,]))

    common_samples <- intersect(colnames(exp_table), feat_table$cell_name)
    exp_table      <- exp_table[,common_samples]
  }

  morgan_type  <- ifelse(target=="tcga_multi", "tcga_multi", "tcga")
  morgan_bits  <- Function_prep_morgan_bits(Function_load_morgan_bits(morgan_type, radii_set, bit_set),
                             type=morgan_type)

} else if (grepl("geeleher_", target)==T){

  if (target=="geeleher_cisplatin"){
    cis <- readRDS(paste0(gee_folder, "030217_GEE_CISPLATIN.rds"))
    feat_table <- cis$feat_table
    exp_table  <- cis$exp_table
    # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_cis_exp_norm.rds"))[["EXP.2"]]

  } else if (target=="geeleher_docetaxel"){
    doc <- readRDS(paste0(gee_folder, "030217_GEE_DOCETAXEL.rds"))
    feat_table <- doc$feat_table
    exp_table  <- doc$exp_table
    # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_doc_exp_norm.rds"))[["EXP.2"]]

  } else if (target=="geeleher_bortezomib_a"){
    bor <- readRDS(paste0(gee_folder, "030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
    feat_table <- bor$feat_table_a
    exp_table  <- bor$exp_table_a
    # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_bora_exp_norm.rds"))[["EXP.2"]]

  } else if (target=="geeleher_bortezomib_b"){
    bor <- readRDS(paste0(gee_folder, "030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds"))
    feat_table <- bor$feat_table_b
    exp_table  <- bor$exp_table_b
    # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_borb_exp_norm.rds"))[["EXP.2"]]

  } else if (target=="geeleher_erlotinib"){
    erl <- readRDS(paste0(gee_folder, "030417_GEE_ERLOTINIB.rds"))
    feat_table <- erl$feat_table
    feat_table$cell_name <- as.character(feat_table$cell_name)
    exp_table  <- erl$exp_table
    # exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_erlo_exp_norm.rds"))[["EXP.2"]]
  }

  exp_table     <- exp_table[,feat_table$cell_name]
  morgan_bits   <- fread(paste0(in_morgan, "CGP_MORGAN_BITS.txt"))[radius==radii_set & bits==bit_set,]
}

# Is batch normalization necessary?
if (batch_norm=="ctrp_cgp"){
  print("bn ctrp_cgp")
  # exp_table       <- readRDS(paste0(in_objects, "041817_cgp_ctrp_exp_norm.rds"))[["EXP.1"]]
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_cgp_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="ctrp_ccle"){
  print("bn ctrp_ccle")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_ccle_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="geeleher_cisplatin"){
  print("bn ctrp_cis")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_cis_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="geeleher_docetaxel"){
  print("bn ctrp_doc")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_doc_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="geeleher_erlotinib"){
  print("bn ctrp_erlo")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_erlo_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="geeleher_bortezomib_a"){
  print("bn ctrp_bora")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_bora_exp_norm.rds"))[["EXP.2"]]
} else if (batch_norm=="geeleher_bortezomib_b"){
  print("bn ctrp_borb")
  exp_table        <- readRDS(paste0(in_objects, "050917_ctrp_borb_exp_norm.rds"))[["EXP.2"]]
}

#MODEL FEATURES
cell_features  <- fread(cgp_cell, header=T, sep="\t")
cell_features  <- colnames(cell_features)[2:ncol(cell_features)]

info_drug_file <- file.info(cgp_drug)
if (info_drug_file$size > 1){
  drug_features  <- fread(cgp_drug, header=T, sep="\t")
  drug_features  <- colnames(drug_features)[2:ncol(drug_features)]
} else{
  drug_features  <- c()
}

############################################################################################################
# EXECUTE

# Build feat tables
if (met_type == "morgan_bits"){
  feat_table <- Function_target_morgan_bits_features_extracted_mf(feat_table, exp_table, morgan_bits,
                                                                         target_cells = cell_features, target_bits = drug_features,
                                                                         genes=genes, scaling=F,
                                                                         genespca = genespca, drugspca = drugspca,
                                                                         original_exp = ctrp_exp, original_bits = ctrp_bits)

} else if (met_type == "morgan_counts"){

  feat_table <- Function_target_morgan_counts_features_extracted_mf(feat_table, exp_table, morgan_counts,
                                                                           target_cells = cell_features, target_counts = drug_features,
                                                                           cgp_exp, genes=F, scaling=T)

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
write.table(feat_table$drug_feat, paste0(out_file, "_", target, "_td"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(feat_table$cell_feat, paste0(out_file, "_", target, "_tc"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(index_table, paste0(out_file, "_", target, "_ti"),
            quote=F, sep="\t", row.names=F, col.names=T)
write.table(feat_table$feat_table, paste0(out_file, "_", target, "_ft"),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing tables")
