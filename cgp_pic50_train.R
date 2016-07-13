#Function to prep cgp pIC50 train tables
library(data.table)
library(reshape2)

main.table <- readRDS("~/Documents/FOLDER/CGP_FILES/cgp_cor_pIC50.rds")
target_drug <- "Erlotinib" #For testing
splits = 0.8

#Split tables
test_table <- main.table[Compound=="Erlotinib",]
main.table <- main.table[Compound!="Erlotinib",]

all_rows <- 1:nrow(main.table)
train_rows <- sample(all_rows, length(all_rows)*0.8)
valid_rows <- setdiff(all_rows, train_rows)

train_table <- main.table[train_rows,]
valid_table <- main.table[valid_rows,]

#Write table
write.table(train_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TRAIN.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(valid_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/VALID.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(test_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TEST.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

print("Done writing tables")
