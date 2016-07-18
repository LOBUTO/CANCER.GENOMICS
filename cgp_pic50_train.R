#Function to prep cgp pIC50 train tables
library(data.table)
library(reshape2)

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#Load files and parameteres
main.table <- readRDS("~/Documents/FOLDER/CGP_FILES/cgp_cor_pIC50.rds")
nci.data   <- readRDS("~/Documents/FOLDER/CGP_FILES/nci60_cgpfeat_cancerlist.rds")

target_drug <- args[1] #For testing
target_cancers <- c("Non-Small Cell Lung") #EXAMPLE for cells related to Erlotinib (nscl)
splits <- 0.8
nci_boost <- F

#Introduce tables
test_table <- main.table[Compound==target_drug,]
target_cells <- unique(test_table$cell_name)

main.table <- main.table[Compound!=target_drug,]
main.table <- main.table[cell_name %in% target_cells,]

if(nci_boost==T){
  nci.table <- do.call(rbind, lapply(target_cancers, function(x) nci.data[[x]]))
  main.table <- rbind(main.table, nci.table)
}

#Split tables
all_rows <- 1:nrow(main.table)
train_rows <- sample(all_rows, length(all_rows)*0.8)
valid_rows <- setdiff(all_rows, train_rows)

train_table <- main.table[train_rows,]
valid_table <- main.table[valid_rows,]

#Normalizing data with train scales
cols = ncol(train_table)
train_mean <- apply(train_table[,4:cols,with=F], 2, mean)
train_sd   <- apply(train_table[,4:cols,with=F], 2, sd)

train_table <- cbind(train_table[,1:3,with=F],
                     sweep(sweep(train_table[,4:cols,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                     )

valid_table <- cbind(valid_table[,1:3,with=F],
                    sweep(sweep(valid_table[,4:cols,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                    )

test_table  <- cbind(test_table[,1:3,with=F],
                    sweep(sweep(test_table[,4:cols,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                    )

#Write table
write.table(train_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TRAIN.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(valid_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/VALID.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(test_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TEST.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

print("Done writing tables")
