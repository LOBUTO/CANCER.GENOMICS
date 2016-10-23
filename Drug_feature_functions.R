# Drug_feature_functions.R
# Set of functions designed to analyze best potential drug features that explained activity

##################################################################################################################

Function_top_cell_features_extracted <- function(feats, exp_table, max_cells=10){
  # Constructs feature table using ONLY cell correlations as features, but limiting to most variable features
  # Max_cells  decreasing variance
  
  # Extract most variable cell features
  cell_feat <- cor(t(scale(log(t(exp_table)))), method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]
  
  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))
  
  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, cell_feat, by="cell_name")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_drug_features_extracted <- function(feats, met_table, max_drugs=10, met_scaled=F){
  # Constructs feature table using drugs correlations ONLY as features, but limiting to most variable features
  # Max_drugs are based in terms of decreasing variance
  
  # Extract most variable drug features - is met scaling necessary?
  if (met_scaled==T){
    met_table[, SD:=sd(TC), by ="METABOLITE"]
    met_table <- met_table[SD!=0,]
    met_table <- met_table[, 1:3, with=F]
    
    met_table <- met_table[, TC := scale(TC), by="METABOLITE"]
  }
  
  drug_feat <- cor(acast(met_table, METABOLITE~DRUG, value.var = "TC"), method = "pearson")
  drug_var  <- data.table(drugs = colnames(drug_feat),
                          VAR   = apply(drug_feat, 2, var))
  top_drugs <- drug_var[order(-VAR),]$drugs[1:max_drugs]
  
  drug_feat <- data.table(drug_feat[, top_drugs], keep.rownames = T)
  setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  
  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, drug_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_metabo_features_extracted <- function(feats, met_table, max_mets=10, met_scaled=F){
  # Constructs feature table using METABOLIVE features, but limiting to most variable features
  # Max_mets are based in terms of decreasing variance
  
  # Extract most variable drug features - is met scaling necessary?
  if (met_scaled==T){
    met_table[, SD:=sd(TC), by ="METABOLITE"]
    met_table <- met_table[SD!=0,]
    met_table <- met_table[, 1:3, with=F]
    
    met_table <- met_table[, TC := scale(TC), by="METABOLITE"]
  }
  
  met_feat  <- acast(met_table, DRUG~METABOLITE, value.var = "TC")
  met_var   <- data.table(mets = colnames(met_feat),
                          VAR  = apply(met_feat, 2, var))
  top_mets  <- met_var[order(-VAR),]$mets[1:max_mets]
  
  met_feat <- data.table(met_feat[, top_mets], keep.rownames = T)
  setnames(met_feat, c("Compound", colnames(met_feat)[2:ncol(met_feat)]))
  
  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, met_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_morgan_features_extracted <- function(feats, morgan_table, max_morgan=10){
  # Constructs feature table using CIRCULAR MORGAN features, but limiting to most variable features
  # Max_morgan are based in terms of decreasing variance
  
  morgan_frame <- data.frame(morgan_table, row.names = 1)
  morgan_var   <- data.table(morgan = colnames(morgan_frame),
                             VAR    = apply(morgan_frame, 2, var))
  top_morgan   <- morgan_var[order(-VAR),]$morgan[1:max_morgan]
  
  morgan_feat  <- morgan_table[, c("Compound", top_morgan), with=F]
  
  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, morgan_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_gene_morgan_features_extracted <- function(feats, exp_table, morgan_table, max_genes=10, max_morgan=10){
  # Constructs feature table using gene expression and CIRCULAR MORGAN features, but limiting to most variable features
  # Max_morgan and max_genes are based in terms of decreasing variance
  
  # Extract top morgan features
  morgan_frame <- data.frame(morgan_table, row.names = 1)
  morgan_var   <- data.table(morgan = colnames(morgan_frame),
                             VAR    = apply(morgan_frame, 2, var))
  top_morgan   <- morgan_var[order(-VAR),]$morgan[1:max_morgan]
  
  morgan_feat  <- morgan_table[, c("Compound", top_morgan), with=F]
  
  # Extract top gene features
  top_genes  <- data.table(genes = rownames(exp_table), 
                           VAR = (apply(exp_table, 1, sd))^2 )
  top_genes  <- top_genes[order(-VAR),]$genes[1:max_genes]
  
  exp_table  <- t(exp_table[top_genes,])
  exp_table  <- data.table(exp_table, keep.rownames = T)
  setnames(exp_table, c("cell_name", colnames(exp_table)[2:ncol(exp_table)]))
  
  # Construct feature table
  feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, exp_table, by="cell_name")
  feat_table <- merge(feat_table, morgan_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_randomize_cell_drug_draw <- function(feat_table, target_drug, random_cells=c(), random_compounds=c()){
  # Function to randomly sample either cell_name or Compound attributes and return filtered down feature table
  
  target_table <- feat_table[Compound==target_drug,]
  feat_table   <- feat_table[Compound!=target_drug,]
  
  all_cells    <- unique(feat_table$cell_name)
  all_comp     <- unique(feat_table$Compound)
  
  # Select random samples
  if (length(random_cells)>0){
    new_cells  <- sample(all_cells, random_cells)
  } else {
    new_cells  <- all_cells
  }
  
  if (length(random_compounds)>0){
    new_compounds <- sample(all_comp, random_compounds)
  } else {
    new_compounds <- all_comp
  }
  
  # Obtain new feature table out of random samples
  feat_table <- feat_table[Compound %in% new_compounds,]
  feat_table <- feat_table[cell_name %in% new_cells,]
  
  # Clean up and Return
  feat_table <- rbind(feat_table, target_table)
  return(feat_table)
}

Function_top_cell_random_drug_features_extracted <- function(feats, exp_table, met_table, max_cells=10, max_drugs=10, pic50_scaled=T, 
                                                             rand="uniform"){
  # Produces a feature table based on cell-cell correlation features and random drug-drug features
  # Both max_cells and max_drugs are based in terms of decreasing variance
  # 'rand' can be "uniform" or "normal"
  
  # Extract most variable cell features
  cell_feat <- cor(t(scale(log(t(exp_table)))), method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]
  
  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))
  
  # Produce random drug features
  if (rand=="uniform"){
    met_table$TC <- runif(length(met_table$TC), min = -1, max = 1)
  } else if (rand=="normal"){
    met_table$TC <- Function.range.minus.plus.1(rnorm(length(met_table$TC)))  
  }
  
  drug_feat <- cor(acast(met_table, METABOLITE~DRUG, value.var = "TC"), method = "pearson")
  drug_var  <- data.table(drugs = colnames(drug_feat),
                          VAR   = apply(drug_feat, 2, var))
  top_drugs <- drug_var[order(-VAR),]$drugs[1:max_drugs]
  
  drug_feat <- data.table(drug_feat[, top_drugs], keep.rownames = T)
  setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  
  # Construct feature table
  if(pic50_scaled==T){
    feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  } else {
    feat_table <- feats[, c("Compound", "cell_name", "pIC50"), with=F]  
  }
  
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, cell_feat, by="cell_name")
  feat_table <- merge(feat_table, drug_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_cell_drug_features_extracted <- function(feats, exp_table, met_table, max_cells=10, max_drugs=10, met_scaled=F, pic50_scaled=T){
  # Constructs feature table using drugs and cell correlations as features, but limiting to most variable features
  # Both max_cells and max_drugs are based in terms of decreasing variance
  
  # Extract most variable cell features
  #cell_feat <- cor(t(scale(log(t(exp_table)))), method = "pearson") #Temporary change
  cell_feat <- cor(cgp_exp, method = "pearson")
  cell_var  <- data.table(cells = colnames(cell_feat),
                          VAR   = apply(cell_feat, 2, var))
  top_cells <- cell_var[order(-VAR),]$cells[1:max_cells]
  
  cell_feat <- data.table(cell_feat[, top_cells], keep.rownames = T)
  setnames(cell_feat, c("cell_name", colnames(cell_feat)[2:ncol(cell_feat)]))
  
  # Extract most variable drug features - is met scaling necessary?
  if (met_scaled==T){
    met_table[, SD:=sd(TC), by ="METABOLITE"]
    met_table <- met_table[SD!=0,]
    met_table <- met_table[, 1:3, with=F]
    
    met_table <- met_table[, TC := scale(TC), by="METABOLITE"]
  }
  
  drug_feat <- cor(acast(met_table, METABOLITE~DRUG, value.var = "TC"), method = "pearson")
  drug_var  <- data.table(drugs = colnames(drug_feat),
                          VAR   = apply(drug_feat, 2, var))
  top_drugs <- drug_var[order(-VAR),]$drugs[1:max_drugs]
  
  drug_feat <- data.table(drug_feat[, top_drugs], keep.rownames = T)
  setnames(drug_feat, c("Compound", colnames(drug_feat)[2:ncol(drug_feat)]))
  
  # Construct feature table
  if(pic50_scaled==T){
    feat_table <- feats[, c("Compound", "cell_name", "NORM.pIC50"), with=F]
  } else {
    feat_table <- feats[, c("Compound", "cell_name", "pIC50"), with=F]  
  }
  
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  
  feat_table <- merge(feat_table, cell_feat, by="cell_name")
  feat_table <- merge(feat_table, drug_feat, by="Compound")
  
  # Clean up and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_top_feat_pca_extracted <- function(feats, exp_table, met_table, target_drug, max_genes=10, max_mets=10){
  # Construct feat table by applying PCA to traininig set and applying training PCA rotation to test set
  # PCA is applied to whole constructed training feature table and then rotation is applied to testing table
  # Both max_genes and max_mets are based on total count of maximum features allowed prior to PCA ranked in order
  # of decreasing variance
  
  # Choose max items
  top_genes  <- data.table(genes = rownames(exp_table), 
                           VAR = (apply(exp_table, 1, sd))^2 )
  top_genes  <- top_genes[order(-VAR),]$genes[1:max_genes]
  
  top_mets   <- met_table[, list(VAR = var(TC)), by="METABOLITE"]
  top_mets   <- top_mets[order(-VAR),]$METABOLITE[1:max_mets]
  
  # Build feature table
  exp_table  <- exp_table[top_genes,]
  met_table  <- met_table[METABOLITE %in% top_mets,]
  
  exp_table  <- data.table(scale(t(log(exp_table))) , keep.rownames = T)
  met_table  <- data.table(acast(met_table, DRUG~METABOLITE, value.var = "TC"), keep.rownames = T)
  setnames(exp_table, c("cell_name", colnames(exp_table)[2:ncol(exp_table)]))
  setnames(met_table, c("Compound",  colnames(met_table)[2:ncol(met_table)]))
  
  feat_table <- feats[,c("Compound", "cell_name", "NORM.pIC50"),with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  feat_table <- merge(feat_table, exp_table, by="cell_name")
  feat_table <- merge(feat_table, met_table, by="Compound")
  
  # Split and apply PCA - use features that explain at least 95% of the variance
  train_set  <- feat_table[Compound==target_drug,]
  test_set   <- feat_table[Compound!=target_drug,]
  
  pca_train  <- prcomp(train_set[,4:ncol(train_set),with=F], center = F, scale. = F)
  train_var  <- pca_train$sdev^2
  train_cum  <- cumsum(train_var / sum(train_var))
  train_95   <- sum(train_cum < 0.95) + 1
  
  pca_test   <-  as.matrix(test_set[,4:ncol(test_set),with=F]) %*% pca_train$rotation[,1:train_95] #Var 95%
  
  pca_train  <- cbind(train_set[,1:3,with=F], data.table(pca_train$x[,1:train_95])) #Var 95%
  pca_test   <- cbind(test_set[,1:3,with=F],  data.table(pca_test))
  
  # Concatenate sets for return
  pca_feat   <- rbind(pca_train, pca_test)

  # Clean and return
  setkey(pca_feat)
  pca_feat   <- unique(pca_feat)
  return(pca_feat)
}

Function_top_feat_pca <- function(feats, exp_table, met_table, max_genes=10, max_mets=10) {
  # Construct feat table like "cgp_new_feat" based on independently extracted expression and structural features
  # feat as in cgp_new, exp_table like "cgp_exp" and met_table like "DRUGS.MET.PROFILE"
  # max_genes numbers chosen in terms for variance and max_mets chosen in terms of how many total mets to allow for PCA calculation
  
  # Returns PCAs that cover 95% of variance
  
  # Choose max items
  top_genes  <- data.table(genes = rownames(exp_table), 
                           VAR = (apply(exp_table, 1, sd))^2 )
  top_genes  <- top_genes[order(-VAR),]$genes[1:max_genes]
  
  top_mets   <- met_table[, list(VAR = var(TC)), by="METABOLITE"]
  top_mets   <- top_mets[order(-VAR),]$METABOLITE#[1:max_mets]
  
  # Apply PCA
  pca_genes  <- prcomp(scale(t(log(exp_table[top_genes,]))), center = F, scale. =F) 
  pca_mets   <- prcomp(acast(met_table[METABOLITE %in% top_mets,], DRUG~METABOLITE, value.var = "TC"), center = F, scale. = F)
  
  # Keep only those that cover at least 95% of variance
  gene_var   <- pca_genes$sdev^2
  gene_cum   <- cumsum(gene_var / sum(gene_var))
  gene_95    <- sum(gene_cum < 0.95) + 1
  
#   met_var    <- pca_mets$sdev^2
#   met_cum    <- cumsum(met_var / sum(met_var))
#   met_95     <- sum(met_cum < 0.95) + 1
  
  # Construct output
  pca_genes  <- pca_genes$x[, 1:gene_95]
  pca_mets   <- pca_mets$x[,  1:max_mets]
  colnames(pca_genes) <- paste0(colnames(pca_genes), ".G")
  colnames(pca_mets)  <- paste0(colnames(pca_mets),  ".M")
  
  pca_genes  <- data.table(pca_genes, keep.rownames = T)
  pca_mets   <- data.table(pca_mets , keep.rownames = T)
  setnames(pca_genes, c("cell_name", colnames(pca_genes)[2:ncol(pca_genes)]))
  setnames(pca_mets,  c("Compound",  colnames(pca_mets)[2:ncol(pca_mets)]))
  
  feat_table <- feats[,c("Compound", "cell_name", "NORM.pIC50"),with=F]
  setnames(feat_table, c("Compound", "cell_name", "NORM_pIC50"))
  feat_table <- merge(feat_table, pca_genes, by="cell_name")
  feat_table <- merge(feat_table, pca_mets, by="Compound")
  
  # Clean and return
  setkey(feat_table)
  feat_table <- unique(feat_table)
  return(feat_table)
}

Function_cell_weights <- function(cell_exp) {
  # Function to obtain cell-cell weights 
  # Expression with samples in rows and genes in columns
  
  w <- data.table(melt(cor(cell_exp, method = "pearson")))
  setnames(w, c("cell_1", "cell_2", "w"))
  
  w$cell_1 <- as.character(w$cell_1)
  w$cell_2 <- as.character(w$cell_2)
  
  setkey(w)
  w <- unique(w)
  
  return(w)
  
}

Function_drug_weights <- function(drug_feat) {
  # Function to obtain drug-drug weights
  # Tables such as DRUGS.MET.PROFILE
  
  w <- acast(drug_feat, METABOLITE~DRUG, value.var = "TC")
  w <- data.table(melt(cor(w, method="pearson")))
  setnames(w, c("drug_1", "drug_2", "w"))
  
  w$drug_1 <- as.character(w$drug_1)
  w$drug_2 <- as.character(w$drug_2)
  
  setkey(w)
  w <- unique(w)
  
  return(w)
}

Function_cgp_weight_separate <- function(cgp_table, target_drug, cell_w, drug_w, exponential=F) {
  
  # Make sure all data is present from weights
  cgp_table    <- cgp_table[cell_name %in% cell_w[[1]],]
  cgp_table    <- cgp_table[Compound %in% drug_w[[1]],]
  cgp_table    <- cgp_table[, c("Compound", "cell_name", "NORM_pIC50"), with=F]
  
  # Split into training and testing set
  train_table  <- cgp_table[Compound!=target_drug,]
  test_table   <- cgp_table[Compound==target_drug,]
  
  # Calculate combined weights
  setnames(drug_w, c("drug_1", "drug_2", "weight_d"))
  setnames(cell_w, c("cell_1", "cell_2", "weight_c"))
  
  drug_w       <- drug_w[drug_2 == target_drug,]
  cell_w       <- cell_w[cell_2  %in% test_table$cell_name,]
  
  weight_table <- merge(train_table,  drug_w, by.x="Compound",  by.y="drug_1", allow.cartesian=TRUE)
  weight_table <- merge(weight_table, cell_w, by.x="cell_name", by.y="cell_1", allow.cartesian=TRUE) 
  #[Compound, cell_name, NORM_pIC50, drug_2, weight_d, cell_2, weight_c]
  
  setkey(weight_table)
  weight_table <- unique(weight_table)
  
  # Predicting each target cell's IC50 in other drugs
  if (exponential == T) {
    
    weight_table$weight_d <- Function.range.0.1(exp(weight_table$weight_d))
    weight_table$weight_c <- Function.range.0.1(exp(weight_table$weight_c))
  } 
  print(weight_table)
  
  weight_table <- weight_table[, list(other_cell_pred = sum(weight_c * NORM_pIC50) / sum(weight_c),
                                      weight_d = weight_d), 
                               by = c("cell_2", "Compound")]
  setkey(weight_table)
  weight_table <- unique(weight_table)
  print(weight_table)
  
  # Predicting each target cell's IC50 in target drug
  weight_table <- weight_table[, list(Prediction = sum(other_cell_pred * weight_d) / sum(weight_d) ), 
                               by=c("cell_2")]
  
  # Outputing
  output       <- merge(test_table, weight_table, by.x="cell_name", by.y="cell_2")
  output$Cor   <- cor(output$NORM_pIC50, output$Prediction)
  output$NRMSE <- Function.NRMSE(output$Prediction, output$NORM_pIC50)
  
  # Return
  return(output)
}

Function_cgp_weight_compounded <- function(cgp_table, target_drug, exponential=F) {
  # Uses whole vector to calculate distances

  # It is numerically faster to calculate the distance per each sample being analyzed
  train_table  <- cgp_table[Compound!=target_drug,]
  target_table <- cgp_table[Compound==target_drug,]
  
  predictions  <- apply(target_table, 1, function(x) {
    print(x[2])
    target_vector <- as.numeric(x[4:ncol(target_table)])
    
    weights       <- apply(train_table, 1, function(y) {
      
      train_vector  <- as.numeric(y[4:ncol(train_table)])
      return(cor(target_vector, train_vector, method="pearson"))
    })
    
    if (exponential==T){
      weights <- Function.range.0.1(exp(weights))
    }
    
    prediction <- sum(weights * train_table$NORM_pIC50) / sum(weights)
    return(prediction)
  })
  
  #Return formatted predictions
  output       <- data.table(cell_name  = target_table$cell_name,
                             Compound   = target_drug,
                             NORM_pIC50 = target_table$NORM_pIC50,
                             Prediction = predictions)
  
  output$Cor   <- cor(output$NORM_pIC50, output$Prediction)
  output$NRMSE <- Function.NRMSE(output$Prediction, output$NORM_pIC50)
  
  return(output)
}

Function_cgp_weight_nearest_compounded <- function(cgp_table, target_drug, k=3) {
  # Uses whole vector to calculate distances (knn neighbors for regression)
  
  require(FNN)
  # It is numerically faster to calculate the distance per each sample being analyzed
  train_table  <- cgp_table[Compound!=target_drug,]
  target_table <- cgp_table[Compound==target_drug,]
  target_cells <- target_table$cell_name
  target_y     <- target_table$NORM_pIC50
  
  # Split for nearest neighbor modeling
  knn_x        <- as.vector(train_table$NORM_pIC50)
  knn_y        <- as.vector(target_table$NORM_pIC50)
  train_table  <- data.matrix(train_table[, 4:ncol(train_table), with=F])
  target_table <- data.matrix(target_table[, 4:ncol(target_table), with=F])
  
  # KNN Model
  knn_model    <- knn.reg(train = train_table, test = target_table, y = knn_x, k = k)
  
  #Return formatted predictions
  output       <- data.table(cell_name  = target_cells,
                             Compound   = target_drug,
                             NORM_pIC50 = target_y,
                             Prediction = knn_model$pred)
  
  output$Cor   <- cor(output$NORM_pIC50, output$Prediction)
  output$NRMSE <- Function.NRMSE(output$Prediction, output$NORM_pIC50)
  
  return(output)
}

Function_cgp_weight_linear_compounded <- function(cgp_table, target_drug, drug_w, lm_weight=F, scaled_features=F, lm_weight_exp="linear") {
  # Uses whole vector to calculate distances (lm for regression)
  
  # Pre-calculate weight vector in case it is needed later on
  train_table   <- cgp_table[Compound!=target_drug,]
  train_drugs   <- train_table$Compound
  
  lm_w_vector   <- drug_w[drug_1==target_drug,][, c("drug_2", "w"), with=F]
  train_table   <- merge(train_table, lm_w_vector, by.x="Compound", by.y="drug_2")
  lm_w_vector   <- train_table$w
  train_table$w <- NULL
  
  target_table  <- cgp_table[Compound==target_drug,]
  target_cells  <- target_table$cell_name
  target_y      <- target_table$NORM_pIC50
  
  # Split for nearest neighbor modeling
  train_table   <- train_table[, 3:ncol(train_table), with=F]
  target_table  <- target_table[, 4:ncol(target_table), with=F]
  
  # Do we need to scale features prior to modeling
  if (scaled_features==T){
    scaled_data  <- scale(train_table[,2:ncol(train_table),with=F])
    train_table  <- cbind(train_table[,c("NORM_pIC50"),with=F], scaled_data)
    
    train_mean   <- attributes(scaled_data)$`scaled:center`
    train_sd     <- attributes(scaled_data)$`scaled:scale`
    
    target_table <- sweep(target_table, 2, train_mean, "-")
    target_table <- sweep(target_table, 2, train_sd, "/")
  }
  
  # Linear Model (Apply weights if necessary)
  if (lm_weight==T){
    
    lm_w_vector <- lm_w_vector[which(lm_w_vector>0)]  # To make weights are between 0-1
    train_table <- train_table[which(lm_w_vector>0),]
    
    if (lm_weight_exp=="exp"){
      lm_w_vector <- exp(lm_w_vector)
      
    } else if (is.numeric(lm_weight_exp)){
      lm_w_vector <- (lm_w_vector)^lm_weight_exp
      
    } else {
      lm_w_vector <- lm_w_vector
    }
    
    lm_model    <- lm(formula = NORM_pIC50~.,  data = train_table, weights = lm_w_vector) 
    
  } else {
    lm_model    <- lm(formula = NORM_pIC50~.,  data = train_table)
  }
  
  lm_predict    <- predict(lm_model, target_table)
  
  # Return formatted predictions
  output        <- data.table(cell_name  = target_cells,
                             Compound   = target_drug,
                             NORM_pIC50 = target_y,
                             Prediction = lm_predict)
  
  output$Cor    <- cor(output$NORM_pIC50, output$Prediction)
  print(unique(output$Cor))
  output$NRMSE  <- Function.NRMSE(output$Prediction, output$NORM_pIC50)
  
  return(output)
}

Function_cgp_model_baseline <- function(cgp_table, target_drug, cell_w, drug_w, th, exponential=F, type="separate", k=3, scaled_features=F,
                                        lm_weight_exp="linear"){
  # Finds baseline prediction for target drug based on all the data and type of prediction
  # This will be in terms of weighted average
  # type can be:
    # "separate"   - using each cell_w and drug_w
    # "compounded" - using whole vector as distance in cgp_table rows [4:ncol]
  
  if (type == "separate"){
    
    output <- Function_cgp_weight_separate(cgp_table, target_drug, cell_w, drug_w, exponential=exponential)
    
  } else if (type == "compounded"){
    
    output <- Function_cgp_weight_compounded(cgp_table, target_drug, exponential=exponential)
    
  } else if (type == "knn_compounded"){
    
    output <- Function_cgp_weight_nearest_compounded(cgp_table, target_drug, k=k)
    
  } else if (type == "linear_compounded"){
    
    output <- Function_cgp_weight_linear_compounded(cgp_table, target_drug, drug_w, lm_weight = F, scaled_features = scaled_features)
    
  } else if (type == "linear_weighted"){
    print("linear_weighted")
    output <- Function_cgp_weight_linear_compounded(cgp_table, target_drug, drug_w, lm_weight = T, scaled_features = scaled_features,
                                                    lm_weight_exp = lm_weight_exp)
    
  } else if (type == "linear_compounded_thresholded"){
    
    cgp_table <- Function_feat_table_threshold_filtering(cgp_table, target_drug, drug_w, th)
    output    <- Function_cgp_weight_linear_compounded(cgp_table, target_drug, drug_w, lm_weight = F, scaled_features = scaled_features)
    
  } else if (type == "linear_weighted_thresholded"){
    
    cgp_table <- Function_feat_table_threshold_filtering(cgp_table, target_drug, drug_w, th)
    output    <- Function_cgp_weight_linear_compounded(cgp_table, target_drug, drug_w, lm_weight = T, scaled_features = scaled_features,
                                                       lm_weight_exp = lm_weight_exp)
  }
  
  # Return
  return(output)
}

Function_feat_table_threshold_filtering <- function(feat_table, target_drug, drug_w, th){
  # Filtering cgp_feat_table -like object by threshold of similarity according to drug_w table
  
  drug_w     <- drug_w[drug_1==target_drug,][w>=th,]
  filt_drugs <- unique(drug_w$drug_2)
  
  feat_table <- feat_table[Compound %in% filt_drugs,]
  
  return(feat_table)
}

Function_plot_comb_tables_plots <- function(table_list, labels=c(), target_index, x_lab, y_lab, title) {
  # Funciton to plot data from multiple tables
  # Assumes both tables have the same column indexes
  # target_index = column namesof index x row in plot
  
  # Add labels
  for (i in 1:length(table_list)){
    table_list[[i]][["Data"]] <- labels[i]
  }
  
  # Merge and assign label
  main_table <- do.call(rbind, table_list)
  
  # Plot relevant info
  main_table <- unique(main_table[,c("Compound", "Cor", target_index, "Data"),with=F])
  main_table[[target_index]] <- as.factor(main_table[[target_index]]) 
  
  ggplot(main_table, aes_string(target_index, "Cor", fill="Data")) + 
    geom_boxplot() + geom_jitter(size=0.4) + scale_fill_brewer(palette = "Set1") +
    theme_bw() + xlab(x_lab) + ylab(y_lab) + ggtitle(title)
}

Function_compare_drug_sim_performance <- function(act_table, sim_list, sim_types, sim_labels) {
  # Compares a list of similarity metrics for their linear cgp-drug predictive performance
  # sim_list is a list of tables/matrices similarities measures
  # sim_types expressed in what format the tables/matrices are submitted
  # sim_labels are the labels desired for each object in sim_list
  
  feat_table <- act_table[,c("cell_name", "Compound", "NORM.pIC50"),with=F]
  
  # Pre-format simarity tables
  all_sim <- data.table()
  for (i in 1:length(sim_list)){
    print(i)
    
    if (sim_types[i]=="feat_met_table"){
      
      sim_table <- cor(acast(sim_list[[i]], METABOLITE~DRUG, value.var="TC"), method="pearson")
      sim_table <- data.table(melt(sim_table))
      
    } else if (sim_types[i]=="feat_morgan_counts"){
      
      sim_table <- cor(acast(sim_list[[i]], Substructure~Compound, value.var = "Counts", fill = 0), method="pearson")
      sim_table <- data.table(melt(sim_table))
      
    } else if(sim_types[i]=="sim_table"){
      
      sim_table <- sim_list[[i]]
      
    } else if(sim_types[i]=="matrix_table"){
      
      sim_table <- cor(t(data.frame(sim_list[[i]], row.names = 1)), method = "pearson")
      sim_table <- data.table(melt(sim_table))
      
    } else if(sim_types[i]=="matrix_table_tc"){
      
      sim_table <- t(data.frame(sim_list[[i]], row.names = 1))
      sim_table <- Function_pairwise_matrix_column(sim_table, Function_tanimoto, c("drug_1", "drug_2", "sim"))
      
    }
    
    setnames(sim_table, c("drug_1", "drug_2", "sim"))
    sim_table$metric <- sim_labels[i]
    all_sim     <- rbind(all_sim, sim_table)
  }
  
  # Combine compound similarities and activity
  act_sim <- acast(feat_table, cell_name~Compound, value.var = "NORM.pIC50")
  act_sim <- cor(act_sim, use = "pairwise.complete.obs" ,method = "pearson")
  act_sim <- data.table(melt(act_sim))
  setnames(act_sim, c("drug_1", "drug_2", "Activity_correlation"))
  
  all_sim   <- merge(all_sim, act_sim, by=c("drug_1", "drug_2"))
  
  # Clean up and return
  return(all_sim)
}

Function_met_randomization_comparisson <- function(met_table) {
  # Function to compare variance in drug features against randomized features
  
  met_table$norm_rand <- Function.range.0.1(rnorm(length(met_table$TC)))
  met_table$unif_rand <- runif(length(met_table$TC))
  
  met_table <- melt(met_table, id.vars = c("DRUG", "METABOLITE"))
  feat_comp <- lapply(unique(met_table$variable), function(x) {
    
    temp_table <- met_table[variable==x,]
    temp_table <- acast(temp_table, METABOLITE~DRUG, value.var = "value")
    temp_table <- cor(temp_table, method="pearson")
    
    return(data.table(type = x, variances = apply(temp_table, 2, var) ))
  })
  
  feat_comp <- do.call(rbind, feat_comp)
  return(feat_comp)
  
}

Function_build_drug_drug_target_ml <- function(drug_table, target_table, target_name, comp_method="abs_diff"){
  # Function that builds feature table using vectorial difference between drugs
  # The target of this function is activiy similarity of compounds across common cell lines
  
  # Build target similarity
  sim_table  <- cor(acast(target_table, cell_name~Compound, value.var = target_name), method="pearson", use = "pairwise.complete.obs")
  sim_table  <- data.table(melt(sim_table))
  fix_name   <- paste0(strsplit(target_name, "[.]")[[1]], collapse = "_")
  setnames(sim_table, c("drug_1", "drug_2", fix_name))
  
  # Build feature table based on comp_method
  feat_table <- acast(drug_table, METABOLITE~DRUG, value.var = "TC")
  
  if (comp_method == "abs_diff"){
    feat_table <- Function_pairwise_matrix_column_for_vectors(feat_table, Function_vector_abs_diff, c("drug_1", "drug_2"))  
  } else if (comp_method == "product"){
    feat_table <- Function_pairwise_matrix_column_for_vectors(feat_table, Function_vector_product, c("drug_1", "drug_2"))  
  }
  
  # Combine and return
  feat_table <- merge(sim_table, feat_table, by=c("drug_1", "drug_2"))
  return(feat_table)
}

Function_drug_drug_target_prediction_lm_knn <- function(feat_table, target_drug, target_name, scaled_features=F, k=3){
  # Linear modeling and KNN prediction using target_name
  # Assumes first 2 columns are identifiers and 3rd columnd is target name
  
  require(FNN)
  
  # Clean up target name (filter for sd==0)
  setnames(feat_table, c("drug_1", "drug_2", "target", colnames(feat_table)[4:ncol(feat_table)]))
  feat_names    <- names(which(apply(feat_table[,4:ncol(feat_table),with=F], 2, sd) != 0))
  feat_table    <- feat_table[, c(colnames(feat_table)[1:3], feat_names), with=F]
  
  # Split tables
  train_table   <- feat_table[drug_1!=target_drug & drug_2!=target_drug, ]
  target_table  <- feat_table[drug_1==target_drug | drug_2==target_drug, ]
  
  target_labels <- target_table[, 1:3,with=F]
  target_table  <- target_table[, 4:ncol(target_table), with=F]
  
  # Perform prediction based on target_name
  if (scaled_features==T){
    scaled_data  <- scale(train_table[,4:ncol(train_table),with=F])
    train_table  <- cbind(train_table[,1:3,with=F], scaled_data)
    
    train_mean   <- attributes(scaled_data)$`scaled:center`
    train_sd     <- attributes(scaled_data)$`scaled:scale`
    
    target_table <- sweep(target_table, 2, train_mean, "-")
    target_table <- sweep(target_table, 2, train_sd, "/")
    
  } 
  
  # Linear Model
  lm_model       <- lm(target~., train_table[,3:ncol(train_table),with=F])
  lm_predict     <- predict(lm_model, target_table)
  
  # KNN Model
  knn_model      <- knn.reg(train = data.matrix(train_table[, 4:ncol(train_table), with=F]), test = target_table, 
                            y = train_table$target, k = k)
  knn_predict    <- knn_model$pred
  
  # Clean up to table
  lm_prediction  <- cbind(target_labels, data.table(Prediction = lm_predict))
  knn_prediction <- cbind(target_labels, data.table(Prediction = knn_predict))
  
  # Add labels and correlation stats
  lm_prediction$Cor   <- cor(lm_prediction$Prediction,  lm_prediction$target, method="pearson")
  knn_prediction$Cor  <- cor(knn_prediction$Prediction, knn_prediction$target, method="pearson")
  lm_prediction$type  <- "Linear Model"
  knn_prediction$type <- "KNN_Model"
  
  # Return
  return(rbind(lm_prediction, knn_prediction))
}
