# cgp_new_sel_predict_tcga_prep.R

library(data.table)
library(reshape2)

Function.exp.combat <- function(exp.matrix.1, exp.matrix.2) {
  #Function to normalize batch effects across two expression matrices
  #INPUT: 2 expression matrices with columns as samples and rows as genes
  #OUPOUT: 2 expression matrices, HOWEVER, genes that are not common across expression datasets will be removed

  require(sva)
  require(limma)

  common.genes <- intersect( unique(rownames(exp.matrix.1)) ,
                             unique(rownames(exp.matrix.2)))

  col.1 <- colnames(exp.matrix.1)
  col.2 <- colnames(exp.matrix.2)

  m.1 <- exp.matrix.1[common.genes,]
  m.2 <- exp.matrix.2[common.genes,]
  m.both <- cbind(m.1, m.2)

  batch <- c( rep("m1", length(col.1)) ,
              rep("m2", length(col.2)) )

  modcombat <- model.matrix(~1, data.frame(1:length(batch)) )

  combat.m <- ComBat(dat=m.both, batch=batch, mod=modcombat, par.prior = T, prior.plots = T)

  exp.matrix.1 <- combat.m[,col.1]
  exp.matrix.2 <- combat.m[,col.2]

  return(list(EXP.1=exp.matrix.1 , EXP.2=exp.matrix.2))
}

Function_new_cgp_sel_feat <- function(cgp_new, drug_met_cor, cgp_exp_cor, target_drug_feat, target_cell_feat,
                                      keep = c()){
  # Function to build total feature set using selected features
  # keep options = "none", "all_cells", "all_drugs", "all"

  # Build cell feature first
  if (length(keep)>0){
    keep_cols  <- intersect(colnames(cgp_exp_cor), keep)
    cell_table <- cgp_exp_cor[,keep_cols]

  } else {
    cell_table <- cgp_exp_cor[ , target_cell_feat]
  }
  cell_table <- data.table(cell_table, keep.rownames = T)
  setnames(cell_table, c("cell_name", colnames(cell_table)[2:ncol(cell_table)]))

  # Build drug feature next
  if (length(keep)>0){
    keep_cols  <- intersect(colnames(drug_met_cor), keep)
    drug_table <- drug_met_cor[,keep_cols]
  } else {
    drug_table <- drug_met_cor[ , target_drug_feat]
  }
  drug_table <- data.table(drug_table, keep.rownames = T)
  setnames(drug_table, c("Compound", colnames(drug_table)[2:ncol(drug_table)]))

  # Merge all with target
  main_table <- cgp_new[,c("cell_name", "Compound", "NORM.pIC50"),with=F]
  setnames(main_table, c("cell_name", "Compound", "NORM_pIC50"))
  main_table <- merge(main_table, cell_table, by=c("cell_name"))
  main_table <- merge(main_table, drug_table, by=c("Compound"))

  # Return
  return(main_table)
}

Function_prcomp_addon <- function(pr_result, original){
  # Function to calculate correlation and contributions of variables to PCs

  # Calculate correlation of variables to original table
  row_order   <- colnames(original)
  #cor_table   <- cor(cbind(pr_result$rotation[row_order,] ,t(original)[row_order,]), method="pearson") #NEED TO FIX!
  #cor_table   <- cor_table[colnames(original), colnames(pr_result$rotation)] #NEED TO FIX!!!

  #Calculate contributions
  #cos2_table <- (cor_table) ^ 2
  #contr_table <- (cos2_table * 100) / apply(cos2_table, 2, sum)
  #contr_table <- sweep(cos2_table, 2, colSums(cos2_table), "/")

  # Calculate contribution of each variable to each PC
  contr_table <- abs(pr_result$rotation)
  contr_table <- sweep(contr_table, 2, colSums(contr_table) , "/")

  # Calculate weighted contribution of each variable across all PCs
  all_contr   <- sweep(contr_table, 2, pr_result$sdev^2, "*")
  all_contr   <- apply(all_contr, 1, sum)

  # Extract those contributing features that pass expected threshold
  exp_contr   <- (pr_result$sdev^2) * (1/nrow(contr_table))
  exp_contr   <- sum(exp_contr)

  contr_feat  <- names(all_contr[all_contr>exp_contr])

  #Return
  #return(list(Cor = cor_table, Contr = contr_table, All_Contr = all_contr, Exp_Thr = exp_contr, Target_feat = contr_feat))
  return(list(Contr = contr_table, All_Contr = all_contr, Exp_Thr = exp_contr, Target_feat = contr_feat))
}

Function_build_feat_cgp <- function(cgp_table, DRUGS.MET.PROFILE, main_exp, side_exp, keep="none", mets=F, genes=F,
                                    chem_cor="pearson", cell_cor="pearson"){
  # Master function to get cgp features using best features derived form pca
  # Note:
  #   main_exp and side_exp should have genes in rows and cells in columns

  # Choose type of chemical features (drugs vs directly metabolites)
  DRUGS.MET.PROFILE  <- DRUGS.MET.PROFILE[DRUG!="CMK",] #Avoiding having same name as cell line

  # Remove duplicates, if any
  duplicates <- data.table(Var1 = c("MG-132", "ABT-263", "Nutlin-3", "Nutlin-3",
                                    "GDC-0449", "GSK-1904529A", "AZD-2281", "AZD6244",    "NVP-TAE684",
                                    "Crizotinib", "Lisitinib", "BEZ235",    "AZD-0530", "AZD-0530"),
                           Var2 = c("Z-LLNle-CHO", "Navitoclax", "Nutlin-3a", "Nutlin-3a (-)",
                                    "GDC0449", "GSK1904529A",   "Olaparib", "selumetinib","TAE684",
                                    "PF-02341066", "OSI-906",  "NVP-BEZ235","AZD0530",  "Saracatinib"))

  for (d in unique(duplicates$Var1)){

    com_drugs <- intersect(c(d, duplicates[Var1==d,]$Var2), unique(DRUGS.MET.PROFILE$DRUG))

    if (length(com_drugs)>1){

      keep_drug   <- com_drugs[1] #Keep the first one
      filter_drug <- setdiff(com_drugs, keep_drug)

      DRUGS.MET.PROFILE <- DRUGS.MET.PROFILE[!(DRUG %in% filter_drug),]
    }
  }

  if (mets==F){

    # Apply correlation
    drug_met_cor <- cor( acast(DRUGS.MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method=chem_cor)

  } else {

    # Cast to have mets as features
    drug_met_cor <- acast(DRUGS.MET.PROFILE, DRUG~METABOLITE, value.var = "TC")
  }

  # Obtain chemical features
  pr_drug_met_cor       <- prcomp(drug_met_cor, center = F, scale. = F)
  pr_drug_met_cor_addon <- Function_prcomp_addon(pr_drug_met_cor, drug_met_cor)

  # Obtain expression features
  if (genes==T){

    combat_scale_main_side         <- Function.exp.combat(t(scale(t(main_exp))), t(scale(t(side_exp))))
    main_exp_scale_combat_cor      <- t(combat_scale_main_side[["EXP.1"]])
    pr_main_exp_scale_combat       <- prcomp(main_exp_scale_combat_cor, center = F, scale. = F)
    pr_main_exp_scale_combat_addon <- Function_prcomp_addon(pr_main_exp_scale_combat, main_exp_scale_combat_cor)

  } else{

    combat_scale_main_side         <- Function.exp.combat(t(scale(t(main_exp))), t(scale(t(side_exp))))

    main_exp_scale_combat_cor      <- cor(combat_scale_main_side[["EXP.1"]], method=cell_cor)

    pr_main_exp_scale_combat       <- prcomp(main_exp_scale_combat_cor, center = F, scale. = F)
    pr_main_exp_scale_combat_addon <- Function_prcomp_addon(pr_main_exp_scale_combat, main_exp_scale_combat_cor)
  }

  # Finally, obtain feature table
  main_new_feat <- Function_new_cgp_sel_feat(cgp_table, drug_met_cor, main_exp_scale_combat_cor,
                                            pr_drug_met_cor_addon$Target_feat, pr_main_exp_scale_combat_addon$Target_feat, keep = keep)

  return(main_new_feat)
}

Function_build_feat_tcga <-function(feat_table, cancer_exp, cgp_exp, DRUGS.MET.PROFILE, tcga_resp, target_drug){

  # Prep expression first
  combat_scale_main_side         <- Function.exp.combat(t(scale(t(cgp_exp))), t(scale(t(cancer_exp))))
  main_exp_scale_combat_cor      <- cor(cbind(combat_scale_main_side[["EXP.1"]], combat_scale_main_side[["EXP.2"]]),
                                        method = "pearson")
  main_exp_scale_combat_cor      <- main_exp_scale_combat_cor[colnames(cancer_exp), colnames(cgp_exp)]

  cell_table   <- data.table(main_exp_scale_combat_cor, keep.rownames = T)
  setnames(cell_table, c("sample", colnames(cell_table)[2:ncol(cell_table)]))

  # Prep chemical features next - Remove duplicates, if any
  duplicates   <- data.table(Var1 = c("MG-132", "ABT-263", "Nutlin-3", "Nutlin-3",
                                    "GDC-0449", "GSK-1904529A", "AZD-2281", "AZD6244",    "NVP-TAE684",
                                    "Crizotinib", "Lisitinib", "BEZ235",    "AZD-0530", "AZD-0530"),
                           Var2 = c("Z-LLNle-CHO", "Navitoclax", "Nutlin-3a", "Nutlin-3a (-)",
                                    "GDC0449", "GSK1904529A",   "Olaparib", "selumetinib","TAE684",
                                    "PF-02341066", "OSI-906",  "NVP-BEZ235","AZD0530",  "Saracatinib"))

  for (d in unique(duplicates$Var1)){

    com_drugs <- intersect(c(d, duplicates[Var1==d,]$Var2), unique(DRUGS.MET.PROFILE$DRUG))

    if (length(com_drugs)>1){

      keep_drug   <- com_drugs[1] #Keep the first one
      filter_drug <- setdiff(com_drugs, keep_drug)

      DRUGS.MET.PROFILE <- DRUGS.MET.PROFILE[!(DRUG %in% filter_drug),]
    }
  }

  DRUGS.MET.PROFILE  <- DRUGS.MET.PROFILE[DRUG!="CMK",] #Avoiding having same name as cell line
  drug_met_cor <- cor( acast(DRUGS.MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
  drug_table   <- data.table(drug_met_cor, keep.rownames = T)
  setnames(drug_table, c("Compound", colnames(drug_table)[2:ncol(drug_table)]))

  # Merge all
  tcga_resp    <- tcga_resp[,c("sample", "Compound", "response"),with=F]
  tcga_resp$response <- as.character(tcga_resp$response)
  tcga_resp    <- merge(tcga_resp, cell_table, by="sample")
  tcga_resp    <- merge(tcga_resp, drug_table, by="Compound")

  # Select for target and filter for prepped feature columns in model training set
  tcga_resp    <- tcga_resp[Compound==target_drug,]
  tcga_resp    <- tcga_resp[ , c("sample", "Compound", "response", colnames(feat_table)[4:ncol(feat_table)]) , with=F]

  # Return
  return(tcga_resp)
}

######################################################################################################
# LOAD DATA
args         <- commandArgs(trailingOnly = TRUE)

cancer       <- args[1]
target_drug  <- args[2]
CANCER       <- cancer

in_folder    <- "/tigress/zamalloa/CGP_FILES/" #For tigress
MET.PROFILE  <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
cgp_new      <- readRDS(paste0(in_folder, "082916_cgp_new.rds"))
cgp_exp      <- readRDS(paste0(in_folder, "083016_cgp_exp.rds"))
cancer_exp   <- readRDS("/tigress/zamalloa/TCGA_FILES/090616_fireshose_all_exp.rds")[[cancer]]
tcga_resp    <- readRDS("/tigress/zamalloa/TCGA_FILES/090616_fireshose_all_response.rds")[cancer==CANCER,]

out_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress (same as in for now)
tcga_table   <- "/tigress/zamalloa/TCGA_FILES/TRAIN_TABLES/"
#####################################################################################################
# EXECUTE

# Minor clean up - removes invariants prior to expression standarization
cancer_exp$tumor <- cancer_exp$tumor[apply(cancer_exp$tumor, 1, sd)!=0,]

# Build cgp model data
feat_table   <- Function_build_feat_cgp(cgp_new, MET.PROFILE, cgp_exp, cancer_exp$tumor, keep=c(),
                                        mets=F, genes=F, chem_cor = "pearson", cell_cor="pearson")
print("Done building feat_table")

# Build tcga drug table for testing
target_table <- Function_build_feat_tcga(feat_table, cancer_exp$tumor, cgp_exp, MET.PROFILE, tcga_resp, target_drug)
print("Done building target_table")

######################################################################################################
# WRITE
saveRDS(feat_table, paste0(out_folder, cancer, "_all_cgp_new_", target_drug, ".rds"))

write.table(target_table, paste0(tcga_table, cancer, "_all_cgp_new_", target_drug) , col.names=T, row.names=F, quote=F, sep="\t")

print("Done prepping cgp table for tcga modeling and tcga target table for analysis")
