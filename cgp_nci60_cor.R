#cgp_nci60_cor.R
#Functions to calculate similarities between nci60 drugs and cgp compounds

library(data.table)
library(reshape2)

CGP_FILES <- "~/Documents/FOLDER/CGP_FILES/"

x <- readRDS(paste0(CGP_FILES, "nci60_cgpfeat_cancerlist.rds"))

main.table <- data.table()
for (names(cancer) in x){

  print(x)

  if (nrow(x[[cancer]])>1){

    y <- x[[cancer]][,c(2,4:142),with=F]
    setkey(y)
    y <- unique(y)

    y <- melt(y, id=c("Compound"))
    setkey(y)

    setnames(y, c("NSC", "CGP", "COR"))

    main.table <- rbind(main.table, y)

  }
}

setkey(main.table)
main.table <- unique(main.table)

#Return
saveRDS(main.table, paste0(CGP_FILES, "nci60_cgp_cor.rds"))
print("Done saving")
