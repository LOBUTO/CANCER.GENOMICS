#cgp_pic50_train_nci.R
library(data.table)
library(reshape2)

Function.cgp.all.data.cor <- function(cgp.table, nci.table, target_drug) {

  non_target_drugs <- unique(cgp.table[Compound!=target_drug,]$Compound)
  nci_target_table <- nci.table[Compound==target_drug,]

  count    <- 1
  total    <- length(non_target_drugs)
  cgp_list <- lapply(non_target_drugs, function(x)  {

    print(c(x, count/total))

    target_table <- rbind(nci_target_table,
                          cgp.table[Compound==x,])

    cgp_labels   <- paste0(target_table$Compound, "_", 1:nrow(target_table))
    cgp_ident    <- target_table[,c("cell_name", "Compound", "NORM_pIC50"),with=F]
    cgp_ident$LABELS <- cgp_labels

    cgp_calc     <- as.matrix(target_table[,4:ncol(target_table), with=F])
    rownames(cgp_calc) <- cgp_labels
    cgp_calc     <- cor(t(cgp_calc))

    cgp_calc     <- data.table(melt(cgp_calc))
    cgp_calc     <- cgp_calc[Var1!=Var2,]
    cgp_calc     <- merge(cgp_calc, cgp_ident, by.x="Var1", by.y="LABELS")
    cgp_calc     <- merge(cgp_calc, cgp_ident, by.x="Var2", by.y="LABELS")

    #Cut at threshold
    cgp_calc     <- cgp_calc[Compound.x==target_drug,] #Table of all distances from each target_drug-cell pair to all others
    print(dim(cgp_calc))

    cgp_min      <- cgp_calc[Compound.y==target_drug,]
    cgp_min[ , min_value := min(value), by = "cell_name.x"] #Minimun distance of each cell in target_drug set to all other cells within drug set
    cgp_min      <- cgp_min[ , c("cell_name.x" , "Compound.x", "min_value"),with=F]
    setkey(cgp_min)
    cgp_min      <- unique(cgp_min)

    cgp_calc     <- merge(cgp_calc, cgp_min, by=c("cell_name.x", "Compound.x"), allow.cartesian=TRUE)
    cgp_calc     <- cgp_calc[value >= min_value,]
      # That means we are only keeping those drug-cell pairs that are similar enough as a within set drug-pair
    print(dim(cgp_calc))

    #Clean up and return
    cgp_calc$target_abs_diff <- abs(cgp_calc$NORM_pIC50.x - cgp_calc$NORM_pIC50.y)
    cgp_calc     <- cgp_calc[,c("cell_name.x", "Compound.x", "NORM_pIC50.x",
                                "cell_name.y", "Compound.y", "NORM_pIC50.y",
                                "value", "target_abs_diff"),with=F]

    count <<- count + 1
    return(cgp_calc)

  })

  #Clean up and return all
  main_table <- do.call(rbind, cgp_list)
  setkey(main_table)
  main_table <- unique(main_table)

  return(main_table)
}


################################################################################
#LOAD FILES
CGP_FILES    <- "/tigress/zamalloa/CGP_FILES/"
TRAIN_TABLES <- "/tigress/zamalloa/CGP_FILES/CGP_TRAIN_TABLES/"

args         <- commandArgs(trailingOnly = TRUE)
target_drug  <- args[1]

cgp.cor.pIC50 <- readRDS(paste0(CGP_FILES, "cgp_cor_pIC50.rds"))
nci.cor.pIC50 <- readRDS(paste0(CGP_FILES, "nci60_cgpfeat_cancertable.rds"))
nci.cor.pIC50$Compound <- as.character(nci.cor.pIC50$Compound)

print("Done loading data")
################################################################################
#EXECUTE

cgp_anal1 <- Function.cgp.all.data.cor(cgp.cor.pIC50, nci.cor.pIC50, target_drug = target_drug)
cgp_anal1 <- unique(cgp_anal1[,c("cell_name.y", "Compound.y"),with=F])
cgp_anal1 <- cgp_anal1[Compound.y != target_drug, ]


sim_table <- merge(cgp.cor.pIC50, cgp_anal1,
              by.x=c("cell_name", "Compound"), by.y=c("cell_name.y", "Compound.y"))

target_table <- nci.cor.pIC50[Compound == target_drug,]

target_table <- rbind(target_table, sim_table)

###################################################################################################################
#STORE

write.table(target_table, paste0(TRAIN_TABLES, "pre_train_",target_drug),
            quote=F, sep="\t", row.names=F, col.names=T)


print("Done pre-treatment")
