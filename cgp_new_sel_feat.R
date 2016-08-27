# cgp_new_sel_feat.R
# Function to split cgp table

library(data.table)
library(reshape2)

Function.cgp.select.best.feat <- function(cgp_new_feat, drug_met_cor, target_drug, weighted="length"){
  # As specified, finds potential best data to train for drug

  # Calculate weight contributions
  drug_index   <- which(colnames(cgp_new_feat) %in% setdiff(colnames(drug_met_cor), "CMK") ) #To avoid confusion with CMK cell
  cell_index   <- setdiff(1:ncol(cgp_new_feat), c(1:3, drug_index))

  if (weighted=="length"){
    cell_weight  <- length(cell_index) / (length(cell_index) + length(drug_index))
    drug_weight  <- length(drug_index) / (length(cell_index) + length(drug_index))
  } else if (weighted=="sd"){
    cell_sd      <- apply(cgp_new_feat[, cell_index, with=F], 2, sd)
    drug_sd      <- apply(cgp_new_feat[, drug_index, with=F], 2, sd)
    cell_weight  <- sum(cell_sd) / (sum(cell_sd) + sum(drug_sd))
    drug_weight  <- sum(drug_sd) / (sum(cell_sd) + sum(drug_sd))
  } else {
    cell_weight <- 0.5
    drug_weight <- 0.5
  }

  # Calculate average contribution of target table
  target_table <- cgp_new_feat[Compound==target_drug,]
  target_cell  <- cor(t(target_table[,cell_index,with=F]), method = "pearson")
  target_mean  <- mean(data.table(melt(target_cell))[Var1!=Var2, ]$value) * cell_weight + 1 * drug_weight

  #Obtain correlation tables for drug and cell line based on features
  cell_cor     <- cgp_new_feat[,c(2, cell_index),with=F]
  setkey(cell_cor)
  cell_cor     <- unique(cell_cor)
  cell_cor     <- as.matrix(data.frame(cell_cor, row.names = 1))
  cell_cor     <- cor(t(cell_cor), method = "pearson")
  cell_cor     <- melt(cell_cor)
  setnames(cell_cor, c("cell_name", "Cell", "cell_cor"))

  drug_cor     <- cgp_new_feat[,c(1, drug_index),with=F]
  setkey(drug_cor)
  drug_cor     <- unique(drug_cor)
  drug_cor     <- as.matrix(data.frame(drug_cor, row.names = 1))
  drug_cor     <- cor(t(drug_cor), method = "pearson")
  drug_cor     <- melt(drug_cor)
  setnames(drug_cor, c("Compound", "Drug", "drug_cor"))

  # Calculate correlation table between cell and drug features for target samples vs rests
  main_table   <- cgp_new_feat[Compound!=target_drug,][,c("Compound", "cell_name"),with=F]
  setkey(main_table)
  main_table   <- unique(main_table)
  main_table$c <- 1

  to_table     <- cgp_new_feat[Compound==target_drug,][,c("Compound", "cell_name"),with=F] # To target table
  setkey(to_table)
  to_table     <- unique(to_table)
  setnames(to_table, c("Drug", "Cell"))
  to_table$c   <- 1

  main_table   <- merge(main_table, to_table, by="c", allow.cartesian=T)
  main_table   <- merge(main_table, cell_cor, by=c("cell_name", "Cell"))
  main_table   <- merge(main_table, drug_cor, by=c("Compound", "Drug"))

  # Calculate filter table for combined weights
  main_table$weighted_score <- main_table$cell_cor * cell_weight + main_table$drug_cor * drug_weight
  main_table   <- main_table[,list(final_score = mean(weighted_score)), by=c("Compound", "cell_name")]
  main_table   <- main_table[final_score >= target_mean,]

  # Return filtered features table
  main_table   <- merge(cgp_new_feat, main_table[,c("Compound", "cell_name"),with=F], by=c("Compound", "cell_name"))
  main_table   <- rbind(main_table, target_table)
  return(main_table)
}

Function.scale.data.table<-function(d.t, col.protect=1:3, criteria=c(), scaling="SCALE"){
  #Scaling can be regular z-scoring ("SCALE") or 0-1 range ("ZERO")

  y.colnames <- colnames(d.t[,-col.protect,with=F])

  if (length(criteria)>0){
    all.criteria <- unique(d.t[[criteria]])

    all.y <- lapply(all.criteria, function(x)  {
      z <- d.t[d.t[[criteria]]==x,]
      z <- data.frame(z[,-col.protect,with=F])
      col.except <- apply(z, 2, sd)==0
      not.col.except <- apply(z, 2, sd)!=0

      if (scaling=="SCALE"){
        z <- cbind(z[,col.except], scale(z[,not.col.except]))
      } else if (scaling=="ZERO"){
        z <- cbind(z[,col.except], apply(z[,not.col.except], 2 ,function(zz) Function.range.0.1(zz)) )
      }

      return (z)
    })
    y <- do.call(rbind, all.y)

  } else {
    y<-data.frame(d.t[,-col.protect,with=F])

    if (scaling=="SCALE"){
      y<-scale(y)
    } else if (scaling=="ZERO"){
      y <- apply(y, 2, function(zz) Function.range.0.1(zz))
    }

  }

  x.colnames<-colnames(d.t[,col.protect,with=F])
  x<-cbind(d.t[,col.protect,with=F], y)
  setnames(x, c(x.colnames, y.colnames))

  return(x)
}


######################################################################################################
# LOAD DATA
args        <- commandArgs(trailingOnly = TRUE)

target_drug <- args[1]
target_drug <- paste0(strsplit(target_drug, "_")[[1]], collapse = " ")
usage       <- args[2]
modifier    <- args[3]

in_folder   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/" #For lab
in_folder   <- "/tigress/zamalloa/CGP_FILES/" #For tigress

MET.PROFILE <- readRDS(paste0(in_folder, "082316.DRUG.MET.PROFILE.rds"))
if (usage=="nci60"){
  feat_table <- readRDS(paste0(in_folder,"082316_cgp_new_feat_combat.rds"))
  #feat_table <- Function.scale.data.table(feat_table, col.protect=1:3)
  #feat_table <- readRDS(paste0(in_folder,"082116_cgp_new_feat_combat_all_drugs.rds"))
  #feat_table  <- readRDS(paste0(in_folder,"121615.CGP.TABLE.PIC50.CELL.SCALED.SCALED.rds"))
} else if (usage=="ccle"){
  feat_table <- readRDS(paste0(in_folder,"081616_cgp_new_feat_combat_ccle_based.rds"))
}

out_table   <- "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/" #For lab
out_table   <- "/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/" #For tigress
#####################################################################################################
# EXECUTE

test_table  <- feat_table[Compound==target_drug]
temp_table  <- feat_table[Compound!=target_drug]

if (modifier=="target_cells"){
  target_cells <- unique(test_table$cell_name)
  temp_table   <- temp_table[cell_name %in% target_cells,]
  temp_table[,NORM_pIC50:=scale(NORM_pIC50), by="Compound"]

} else if (modifier=="target_drugs"){
  drug_met_cor <- cor( acast(MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
  drug_met_cor <- data.table(melt(drug_met_cor))
  drug_met_cor <- drug_met_cor[Var1==target_drug,][value>0.95,]

  target_drugs <- unique(drug_met_cor$Var2)
  temp_table   <- temp_table[Compound %in% target_drugs,]

} else if (modifier=="both"){
  drug_met_cor <- cor( acast(MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
  drug_met_cor <- data.table(melt(drug_met_cor))
  drug_met_cor <- drug_met_cor[Var1==target_drug,][value>0.95,]

  target_drugs <- unique(drug_met_cor$Var2)
  target_cells <- unique(test_table$cell_name)
  temp_table   <- temp_table[Compound %in% target_drugs,][cell_name %in% target_cells,]

} else if (modifier=="target_feat"){
  drug_met_cor <- cor( acast(MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
  temp_table   <- Function.cgp.select.best.feat(feat_table, drug_met_cor, target_drug, weighted="length")
  temp_table   <- temp_table[Compound!=target_drug,]

} else if (modifier=="target_met"){
  feat_table <- readRDS(paste0(in_folder,"082216_cgp_new_feat_combat_met.rds"))
  test_table <- feat_table[Compound==target_drug]
  temp_table <- feat_table[Compound!=target_drug]
} else {
  print("no modifier")
}

train_rows  <- sample(1:nrow(temp_table), 0.8*nrow(temp_table))
valid_rows  <- setdiff(1:nrow(temp_table), train_rows)

train_table <- temp_table[train_rows, ]
valid_table <- temp_table[valid_rows, ]

# Obtain weigthed index
drug_met_cor <- cor( acast(MET.PROFILE, METABOLITE~DRUG, value.var = "TC")  , method="pearson")
drug_met_cor <- data.table(melt(drug_met_cor)) # [Var1, Var2, value]
drug_met_cor <- drug_met_cor[Var1==target_drug,]

train_w      <- sapply(train_table$Compound, function(x) unique(drug_met_cor[Var2==x]$value))
train_w      <- data.table(W = train_w)
valid_w      <- sapply(valid_table$Compound, function(x) unique(drug_met_cor[Var2==x]$value))
valid_w      <- data.table(W = valid_w)

######################################################################################################
# WRITE
write.table(train_table, paste0(out_table, usage , ".", modifier, "_TRAIN_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_table, paste0(out_table, usage , ".", modifier, "_VALID_CGP_SEL.", target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(test_table,  paste0(out_table, usage , ".", modifier, "_TEST_CGP_SEL.",  target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(train_w,     paste0(out_table, usage , ".", modifier, "_TRAIN_CGP_SEL.", target_drug, ".weights"),
            quote=F, sep="\t", row.names=F, col.names=T)

write.table(valid_w,     paste0(out_table, usage , ".", modifier, "_VALID_CGP_SEL.", target_drug, ".weights"),
            quote=F, sep="\t", row.names=F, col.names=T)


print("Done writing sel tables")
