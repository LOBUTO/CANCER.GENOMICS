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

#Normalizing data with train scales
train_mean <- apply(train_table[,4:ncol,with=F], 2, mean)
train_sd   <- apply(train_table[,4:ncol,with=F], 2, sd)

train_table <- cbind(train_table[,1:3,with=F],
                     sweep(sweep(train_table[,4:ncol,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                     )

valid_table <- cbind(valid_table[,1:3,with=F],
                    sweep(sweep(valid_table[,4:ncol,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                    )

test_table  <- cbind(test_table[,1:3,with=F],
                    sweep(sweep(test_table[,4:ncol,with=F], 2, train_mean, "-"), 2, train_sd, "/")
                    )

#Write table
write.table(train_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TRAIN.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(valid_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/VALID.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

write.table(test_table, paste0("~/Documents/FOLDER/TABLES/CGP.TRAINING/TEST.", target_drug, ".pIC50.csv"),
            sep="\t", quote=F, row.names=F, col.names=T)

print("Done writing tables")
