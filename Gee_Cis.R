# Gee_Cis.R
# Process Cisplatin clinical drug data from Geeleher 

library(data.table)
library(reshape2)
library(affy)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(ggplot2)

Function_cel_process <- function(geo.folder){
  geo.files<-paste(geo.folder, list.files(geo.folder), sep="/")
  
  raw.data=read.celfiles(verbose=FALSE, filenames=geo.files)
  #raw.data=ReadAffy(verbose=TRUE, filenames=geo.files)
  data.rma.norm=rma(raw.data)
  rma=exprs(data.rma.norm)
  
  #Extract probe ids, entrez symbols, and entrez ids
  probes=row.names(rma)
  # Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, ifnotfound=NA))
  Symbols = unlist(mget(probes, hgu133plus2SYMBOL, ifnotfound=NA))
  # Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, ifnotfound=NA))
  
  #Assign hugo symbols to matrix
  rma<-rma[!is.na(Symbols),]
  rownames(rma)<-Symbols[!is.na(Symbols)]
  
  #Obtain max of gene for more than one probe per gene
  rma<-data.table(rma,keep.rownames = T  )
  rma<-rma[,lapply(.SD, max), by=rn]
  rma<-data.frame(rma,row.names = 1)
  rma<-data.matrix(rma)
  
  #Clean headers
  colnames(rma)<-substr(colnames(rma), 1, 9)  
  
  #Return
  return(rma)
}

# Load data
date    <- Sys.Date()
samples <- "GSM467523 GSM467524 GSM467525 GSM467526 GSM467527 GSM467528 GSM467529 GSM467530 GSM467531 GSM467532 GSM467533 GSM467534 GSM467535 GSM467536 GSM467537 GSM467538 GSM467539 GSM467540 GSM467541 GSM467542 GSM467543 GSM467544 GSM467545 GSM467546 GSM467547 GSM467548 GSM467549 GSM467550 GSM467551 GSM467552 GSM467553 GSM467554 GSM467555 GSM467556 GSM467557 GSM467558 GSM467559 GSM467560 GSM467561 GSM467562 GSM467563 GSM467564 GSM467565 GSM467566 GSM467567 GSM467568 GSM467569 GSM467570 GSM467571 GSM467572 GSM467573 GSM467574 GSM467575 GSM467576 GSM467577 GSM467578 GSM467579 GSM467580 GSM467581 GSM467582 GSM467583 GSM467584 GSM467585 GSM467586 GSM467587 GSM467588 GSM467589 GSM467590 GSM467591 GSM467592 GSM467593 GSM467594 GSM467595 GSM467596 GSM467597 GSM467598 GSM467599 GSM467600 GSM467601 GSM467602 GSM467603 GSM467604 GSM467605 GSM467606"
samples <- strsplit(samples, " ")[[1]]

classes <- '"miller-payne response: 3"	"miller-payne response: 4"	"miller-payne response: 5"	"miller-payne response: 1"	"miller-payne response: 5"	"miller-payne response: 1"	"miller-payne response: 4"	"miller-payne response: 4"	"miller-payne response: 4"	"miller-payne response: 3"	"miller-payne response: 1"	"miller-payne response: 2"	"miller-payne response: 0"	"miller-payne response: 1"	"miller-payne response: 5"	"miller-payne response: 0"	"miller-payne response: 2"	"miller-payne response: 3"	"miller-payne response: 2"	"miller-payne response: 3"	"miller-payne response: 0"	"miller-payne response: 0"	"miller-payne response: 2"	"miller-payne response: 5"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"	"miller-payne response: n/a"'
classes <- sapply(strsplit(classes, split = '"	"')[[1]], function(x) strsplit(x,": ")[[1]][2])


# Process
feat_table <- data.table(Compound="Cisplatin",
						 cell_name=samples[1:24],
						 target=as.numeric(classes[1:24]))

feat_table$target <- ifelse(feat_table$target>=3, 1, 0) #Positive response is 1 and negative response is 0
print(feat_table)

exp_table  <- Function_cel_process("/Users/jzamalloa/Documents/Rotation/DATABASES/GEELEHER/CISPLATIN/GSE18864_RAW")
print(dim(exp_table))

# Save files
saveRDS(list(exp_table = exp_table, feat_table = feat_table), file="/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/030217_GEE_CISPLATIN.rds")

# Make some nice plots
obj_plot <- feat_table[,list(N=length(cell_name)), by="target"]
obj_plot$target <- factor(obj_plot$target)

pdf(file=paste0("~/Documents/FOLDER/LAB/GM/033017/",date,"_silver_set_1.pdf"), height=8, width=8)
ggplot(obj_plot,
       aes(factor(target), N, fill=target)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label=N), size=5, position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw() + ylab("Number of samples") + xlab("Responder Class") +
  ggtitle("Triple Negative Breast cancer samples\nclassified by response to treatment", subtitle = "Cisplatin - Silver et al., 2010") +
  scale_y_continuous(limits = c(0,20)) +
  theme(axis.text  = element_text(size=12, colour="grey20"),
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, size = 18, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        strip.text.x = element_text(size = 12, colour = "black")) +
  scale_fill_brewer(palette = "Spectral")
dev.off()

print("Done")