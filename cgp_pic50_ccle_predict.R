#cgp_pic50_ccle_predict.R
#Graph for prediction
library(data.table)
library(reshape2)
library(ggplot2)

####################################################################################
# LOAD FILES
main_table <- fread("CGP_FILES/CCLE_RESULTS/ccle_specific_cgp_based_prediction.csv", sep ="\t",
                    colClasses=c("character", "character", "character", "numeric", "numeric"))

file_out <- paste0("FIGURES/CGP.MLP/", as.character(Sys.Date()), "_on_ccle.pdf")

####################################################################################
# EXECUTE
main_table <- main_table[, list(Cor = cor(ACTUAL, PREDICTED)), by=c("cgp_base", "ccle", "layers")]

main_table$based_model <- main_table$cgp_base == main_table$ccle

main_table$has_base    <- main_table$ccle %in% main_table$cgp_base

#Plot
main_table$layers <- factor(main_table$layers, levels = c("50.50", "100.100", "200.200", "400.400", "500.500"), ordered=T)

main_ccle <- main_table[has_base==T,]

pdf(file_out, width=12, height=8)

print (ggplot(main_table, aes(layers, Cor, fill=based_model)) + geom_boxplot() + geom_jitter(size=0.5) +
        facet_grid(~has_base) + theme_classic() + scale_fill_brewer(palette="Set1") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
        ggtitle("Prediction on ccle drugs using cgp based model") +
        xlab("MLP architecture") + ylab("Predicted/Actual correlation")
        )

print (ggplot(main_ccle, aes(layers, Cor, fill=based_model)) + geom_boxplot() + geom_jitter(size=0.5) +
        facet_grid(~ccle) + theme_classic() + scale_fill_brewer(palette="Set1") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
        )

dev.off()

print ("Done plotting")
