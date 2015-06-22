#Function.process.maf.files
#061414
#Takes .maf files and returns patient sample, gene and number of mutations per gene in 3 columns as data.table

args<-commandArgs(trailingOnly=T)

main.maf<-args[1]
output.file<-args[2]

if (length(args)>2) {
  additional.mafs<-args[3:length(args)]
} else
  additional.mafs<-c()

Function.process.maf.files<-function(main.file=c(), additional.files=c()) {
  #Process maf.files from TCGA to 2 objects:
  #First is a table containing:
  #   Tumor_Sample_Barcode
  #   Hugo_Symbol
  #   N.MUTATIONS per gene per patient sample
  #Second is a vector containing all theoretically sequenced genes
  #   This is regardless of type ("silent"), since this are all possible sequenced genes
  
  require(data.table)
  require(reshape2)

  #Process main.file
  dummy.a<-read.csv(main.file, header=T,  sep="\t", stringsAsFactors=F)
  dummy.a$Line_Number<-NULL #To remove duplicates
  dummy.a<-unique(dummy.a)
  SEQUENCED.GENES<-unique(as.vector(dummy.a$Hugo_Symbol)) #ALL POTENTIALLY SEQUENCED GENES
  
  #Filter for info we want
  dummy.a<-dummy.a[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Start_Position", "Variant_Classification")] #Only want sample names and mutations
  dummy.a<-as.data.table(dummy.a)
  
  for (files in additional.files) {
    a.1<-read.csv(files, header=T, sep="\t", stringsAsFactors=F)
    a.1$Line_Number<-NULL
    a.1<-unique(a.1)
    
    SEQUENCED.GENES<-unique(c(SEQUENCED.GENES, unique(as.vector(a.1$Hugo_Symbol))))
    
    a.1<-a.1[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Start_Position","Variant_Classification")] #Only want sample names and mutations
    a.1<-as.data.table(a.1)
    
    #Combine and unique
    dummy.a<-unique(rbind(dummy.a, a.1))
  }
  
  dummy.a$Start_Position<-NULL
  dummy.a$Variant_Classification[dummy.a$Variant_Classification!="Silent"]<-"Missense" 

  #To classify mutations in only two categories
  dummy.a$fill<-1
  dummy.a<-as.data.table(dcast(dummy.a, Tumor_Sample_Barcode + Hugo_Symbol ~ Variant_Classification, value.var="fill", fill=0 ,fun.aggregate=sum))
  dummy.a$Hugo_Symbol<-as.character(dummy.a$Hugo_Symbol)
  
  #Convert Tumor sample barcode to patient barcode
  dummy.a$PATIENT<-sapply(dummy.a$Tumor_Sample_Barcode, function(x) 
    paste0(strsplit(x, "-")[[1]][1:4], collapse="."))
  dummy.a$Tumor_Sample_Barcode<-NULL

  #NOTE!! - Remove duplicated elements due to different "Start_Position" by different agencies, will keep the maximum reported number of mutations
  dummy.a<-dummy.a[order(PATIENT,-Missense),]
  dummy.a<-dummy.a[!duplicated(dummy.a[,c("PATIENT", "Hugo_Symbol"), with=F]),]

  #Return
  dummy.list<-list(table.1=dummy.a, background.genes=SEQUENCED.GENES)
  return (dummy.list)
}

#Call function
dummy.output<-Function.process.maf.files(main.maf, additional.mafs)

#Write to output file
saveRDS(dummy.output, file=output.file)
