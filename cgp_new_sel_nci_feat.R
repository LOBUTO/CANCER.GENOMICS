# cgp_new_sel_nci_feat.R
# Function to predict nci60 drug activity using cgp_new on selected features

library(data.table)
library(reshape2)

Function.exp.combat <- function(exp.matrix.1, exp.matrix.2) {
  #Function to normalize batch effects across two expression matrices
  #INPUT: 2 expression matrices with columns as samples and rows as genes
  #OUPOUT: 2 expression matrices, HOWEVER, genes that are not common across expression datasets will be removed

  require(sva)
  require(pamr)
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

######################################################################################################
# LOAD DATA

drug_hmdb   <- fread("TABLES/TCGA.TRAINING/NCI60.TC.HMDB.FP4",
                    header=T, colClasses=c("character", "character", "numeric"))

#feat_table  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716_cgp_new_feat.rds")
#feat_table <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080816_cgp_new_feat_combat.rds")
feat_table <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_cgp_new_feat_combat.rds")
MET.PROFILE <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.DRUG.MET.PROFILE.rds")
#nci60.exp   <- readRDS("/home/zamalloa/Documents/FOLDER/OBJECTS/061916.NCI60.EXP.rds")
#cgp_exp     <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716_cgp_new_exp.rds")
nci60.exp   <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_NCI60_EXP_COMBAT.rds")
cgp_exp     <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/081016_CGP_EXP_COMBAT.rds")
nci60.gi50  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci.gi50.rds")
nci_to_cgp  <- readRDS("/home/zamalloa/Documents/FOLDER/CGP_FILES/080716.nci60_names_to_cgp.rds")

args        <- commandArgs(trailingOnly = TRUE)
target_drug <- args[1]
target_drug <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ")
out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

#####################################################################################################
# EXECUTE

# Prep drug table first
target_nsc   <- as.character(nci_to_cgp[Compound==target_drug,]$NSC)
drug_hmdb    <- drug_hmdb[DRUG %in% target_nsc,]
drug_hmdb$DRUG <- target_drug
MET.PROFILE  <- MET.PROFILE[DRUG!="CMK",] #Necessary to avoid confusion of also "CMK" cell name later
common_met   <- intersect(unique(drug_hmdb$METABOLITE), unique(MET.PROFILE$METABOLITE))
drug_table   <- rbind(drug_hmdb[METABOLITE %in% common_met,],
                    MET.PROFILE[METABOLITE %in% common_met,][!DRUG %in% target_drug,])

drug_table   <- cor( acast(drug_table, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
drug_table   <- drug_table[target_drug, unique(MET.PROFILE$DRUG), drop=F]
drug_table   <- data.table(drug_table, keep.rownames=T)

setnames(drug_table, c("Compound", colnames(drug_table)[2:ncol(drug_table)]))

# Then expression table
common_genes <- intersect(rownames(nci60.exp), rownames(cgp_exp))

colnames(nci60.exp) <- as.vector(sapply(colnames(nci60.exp), function(x) paste0(x, ".nci")))
colnames(cgp_exp)   <- as.vector(sapply(colnames(cgp_exp),   function(x) paste0(x, ".cgp")))

cell_table   <- cbind(nci60.exp[common_genes, ],
                      cgp_exp[common_genes,])
cell_table   <- cor(cell_table, method="pearson")
cell_table   <- cell_table[colnames(nci60.exp), colnames(cgp_exp)]
rownames(cell_table) <- sapply(rownames(cell_table), function(x)  strsplit(x, ".nci")[[1]][1])
colnames(cell_table) <- sapply(colnames(cell_table), function(x)  strsplit(x, ".cgp")[[1]][1])

cell_table   <- data.table(cell_table, keep.rownames = T)
setnames(cell_table, c("cell_name", colnames(cell_table)[2:ncol(cell_table)]))

# Then merge it all
main_table   <- merge(nci60.gi50, nci_to_cgp, by="NSC")
main_table   <- main_table[ ,c("CELL", "Compound", "SCALE.ACT"),with=F] # To be in accordance with nci60 names ("cell_name" is a cgp identifier)
setnames(main_table, c("cell_name", "Compound", "NORM_pIC50")) # To be in accordance with the pred.py script

main_table   <- merge(main_table, drug_table, by = "Compound")
main_table   <- merge(main_table, cell_table, by = "cell_name")

# Filter features and return
main_table   <- main_table[ , c("cell_name", "Compound", "NORM_pIC50", colnames(feat_table)[4:ncol(feat_table)] ), with= F]
setkey(main_table)
main_table   <- unique(main_table)

write.table(main_table, paste0(out_table, "NCI_CGP_SEL_FEAT.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

print("Done writing nci60 features")
