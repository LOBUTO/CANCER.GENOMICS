# Gee_Bor_Dex.R
library(data.table)
library(affy)
library(GEOquery)
library(hgu133a.db)
library(hgu133b.db)
library(ggplot2)

Function_exp_symbol <- function(rma, Symbols){
  
  rma <- rma[!is.na(Symbols),]
  rownames(rma) <- Symbols[!is.na(Symbols)]
  
  #Obtain max of gene for more than one probe per gene
  rma <- data.table(rma,keep.rownames = T  )
  rma <- rma[,lapply(.SD, mean), by=rn]
  rma <- data.frame(rma,row.names = 1)
  rma <- data.matrix(rma)
  
  #Return
  return(rma)
}

Function_process_class <- function(class_obj){
  
  class_table <- class_obj[,c("geo_accession", "characteristics_ch1.8", "characteristics_ch1.1"), with=F]
  class_table$target <- ifelse(class_table$characteristics_ch1.8 == "PGx_Responder = R", 1, #Sensitive (1), Resistant (0)
                           ifelse(class_table$characteristics_ch1.8 == "PGx_Responder = NR", 0, 2))
  class_table <- class_table[target!=2,]
  class_table$Compound <- sapply(as.character(class_table$characteristics_ch1.1), function(x) strsplit(x, " = ")[[1]][2])
  class_table <- class_table[,-c("characteristics_ch1.8", "characteristics_ch1.1"),with=F]
  setnames(class_table, c("cell_name", "target", "Compound"))

  class_table$Compound <- ifelse(class_table$Compound == "PS341", "Bortezomib", 
  								 ifelse(class_table$Compound == "Dex", "Dexamethasone", "None"))
  
  return(class_table)
}

# Load data
date      	    <- Sys.Date()
bortezomib_mas5 <- getGEO("GSE9782", getGPL = F)

# Process expression
a_exp <- exprs(bortezomib_mas5[[1]])
b_exp <- exprs(bortezomib_mas5[[2]])

a_probes <- row.names(a_exp)
b_probes <- row.names(b_exp)

a_Symbols = unlist(mget(a_probes, hgu133aSYMBOL, ifnotfound=NA))
b_Symbols = unlist(mget(b_probes, hgu133bSYMBOL, ifnotfound=NA))

a_exp <- Function_exp_symbol(a_exp, a_Symbols)
b_exp <- Function_exp_symbol(b_exp, b_Symbols)

# Get classes
a_obj <- data.table(pData(phenoData(bortezomib_mas5[[1]])))
b_obj <- data.table(pData(phenoData(bortezomib_mas5[[2]])))

a_class <- Function_process_class(a_obj)
b_class <- Function_process_class(b_obj)
print(a_class)
print(dim(a_exp))
print(b_class)
print(dim(b_exp))

print(a_class[Compound=="None",])
print(b_class[Compound=="None",])

# Save files
saveRDS(list(exp_table_a = a_exp, feat_table_a = a_class,
             exp_table_b = b_exp, feat_table_b = b_class), 
        file="/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030417_GEE_BORTEZOMIB_DEXAMETHASONE.rds")

# Make some nice plots
# a_obj_plot <- a_obj[,list(N=length(status), Response=characteristics_ch1.7, Responder_Class=characteristics_ch1.8, Treatment = characteristics_ch1.1), 
#                     by=c("characteristics_ch1.1", "characteristics_ch1.7", "characteristics_ch1.8")]
# a_obj_plot$Response <- factor(c("No change", "Partial response", "Minimal response", "IE", "Progressive disease", "Complete response",
#                                 "Progressive disease", "Partial response", "Minimal response", "No change", "IE", "Complete response"))

# pdf(file=paste0("~/Documents/FOLDER/LAB/GM/033017/",date,"_mulligan_set_1.pdf"), height=11, width=10)
# ggplot(a_obj_plot,
#        aes(factor(Responder_Class), N, fill=Response)) +
#   geom_bar(stat="identity", position="dodge") +
#   facet_wrap(~Treatment, ncol = 1) +
#   geom_text(aes(label=N), size=5, position=position_dodge(width=0.9), vjust=-0.25) +
#   theme_bw() + ylab("Number of samples") + xlab("Responder Class") +
#   ggtitle("Myeloma samples classified by response to treatment", subtitle = "Mulligan et al., 2007") +
#   scale_y_continuous(limits = c(0,65)) +
#   theme(axis.text  = element_text(size=12, colour="grey20"),
#         axis.title = element_text(size=14, face="bold"),
#         plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
#         plot.subtitle = element_text(hjust = 0.5, size = 13),
#         strip.text.x = element_text(size = 12, colour = "black")) +
#   scale_fill_brewer(palette = "Spectral")
# dev.off()

print("Done")