# Gee_Doc.R
# Process Docetaxel clinical drug data from Geeleher 

library(data.table)
library(reshape2)
library(affy)
library(oligo)
library(hgu95av2hsentrezgcdf)
library(hgu95av2hsentrezg.db)
library(ggplot2)

Function_cel_process <- function(geo.files, exp_proc){
  
  # raw.data=read.celfiles(verbose=FALSE, filenames=geo.files)
  raw.data=ReadAffy(verbose=TRUE, filenames=geo.files)
  data.rma.norm=affy::rma(raw.data)
  rma=exprs(data.rma.norm)
  
  rma <- merge(data.table(rma, keep.rownames=T), exp_proc, by.x="rn", by.y="ID_REF")
  rma <- rma[,-c("rn"), with=F]

  #Obtain max of gene for more than one probe per gene
  rma<-rma[,lapply(.SD, max), by=IDENTIFIER]
  rma<-data.frame(rma,row.names = 1)
  rma<-data.matrix(rma)
  
  #Clean headers
  colnames(rma)<-substr(colnames(rma), 1, 7)  
  
  #Return
  return(rma)
}

# Load data
date      <- Sys.Date()
resistant <- "GSM4901,GSM4902,GSM4904,GSM4905,GSM4906,GSM4909,GSM4910,GSM4911,GSM4912,GSM4913,GSM4916,GSM4918,GSM4922,GSM4924"
sensitive <- "GSM4903,GSM4907,GSM4908,GSM4914,GSM4915,GSM4917,GSM4919,GSM4920,GSM4921,GSM4923"

resistant <- strsplit(resistant, ",")[[1]]
sensitive <- strsplit(sensitive, ",")[[1]]

# x         <- fread("/Users/jzamalloa/Documents/Rotation/DATABASES/GEELEHER/DOCETAXEL/EXP_PROC.txt", header=T, sep="\t", drop=1)

# Process
# x <- aggregate(.~IDENTIFIER, x, mean)
# x <- data.frame(x, row.names = 1)
# x <- data.matrix(x)

# x <- x[,c(resistant, sensitive)]

cel_folder <- "/Users/jzamalloa/Documents/Rotation/DATABASES/GEELEHER/DOCETAXEL/"
files_349  <- paste0(cel_folder, "GSE349_RAW/")
files_349  <- list.files(files_349, full.names=T)
files_350  <- paste0(cel_folder, "GSE350_RAW/")
files_350  <- list.files(files_350, full.names=T)
cel_files  <- c(files_349, files_350)

exp_proc   <- fread("/Users/jzamalloa/Documents/Rotation/DATABASES/GEELEHER/DOCETAXEL/EXP_PROC.txt", header=T, sep="\t", select=c(1,2))

exp_table  <- Function_cel_process(cel_files, exp_proc)
print(dim(exp_table))
feat_table <- rbind(data.table(Compound="Docetaxel", cell_name = resistant, target=0),
					data.table(Compound="Docetaxel", cell_name = sensitive, target=1))
print(feat_table)
# Save files
saveRDS(list(exp_table = exp_table, feat_table = feat_table), file="/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030217_GEE_DOCETAXEL.rds")

# Make some nice plots
obj_plot <- feat_table[,list(N=length(cell_name)), by="target"]
obj_plot$target <- factor(obj_plot$target)

pdf(file=paste0("~/Documents/FOLDER/LAB/GM/033017/",date,"_chang_set_1.pdf"), height=8, width=8)
ggplot(obj_plot,
       aes(factor(target), N, fill=target)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label=N), size=5, position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() + ylab("Number of samples") + xlab("Responder Class") +
  ggtitle("Breast cancer samples classified by response to treatment", subtitle = "Docetaxel - Chang et al., 2003") +
  scale_y_continuous(limits = c(0,20)) +
  theme(axis.text  = element_text(size=12, colour="grey20"),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        strip.text.x = element_text(size = 12, colour = "black")) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

print("Done")