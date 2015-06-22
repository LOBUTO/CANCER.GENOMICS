#Function.process.maf.files.bg
#Difference from Function.process.maf.files() is that this is for "Silent" mutations only
#061814
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
  
  #Process main.file
  dummy.a<-read.csv(main.file, header=T,  sep="\t")
  dummy.a$Line_Number<-NULL #To remove duplicates
  dummy.a<-unique(dummy.a)
  SEQUENCED.GENES<-unique(as.vector(dummy.a$Hugo_Symbol)) #ALL POTENTIALLY SEQUENCED GENES
  dummy.a<-dummy.a[dummy.a$Variant_Classification=="Silent",] #Select for silent mutations
  
  #Filter for info we want
  dummy.a<-dummy.a[,c("Tumor_Sample_Barcode", "Hugo_Symbol")] #Only want sample names and mutations
  dummy.a<-as.data.table(dummy.a)
  dummy.a$fill<-1 #To count number of mutations per gene
  dummy.a<-dummy.a[,list(N.MUTATIONS=length(fill)), by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]
  
  #Process additional files:
  for (files in additional.files) {
    a.1<-read.csv(files, header=T, sep="\t")
    
    #Add to background genes
    SEQUENCED.GENES<-unique(c(SEQUENCED.GENES, unique(as.vector(a.1$Hugo_Symbol))))
    
    #Only keep samples not found in main.file
    a.1<-a.1[!(a.1$Tumor_Sample_Barcode %in% unique(as.vector(dummy.a$Tumor_Sample_Barcode))),]
    
    #Continue processing
    a.1$Line_Number<-NULL
    a.1<-unique(a.1)
    a.1<-a.1[a.1$Variant_Classification!="Silent",]
    a.1<-a.1[,c("Tumor_Sample_Barcode", "Hugo_Symbol")] #Only want sample names and mutations
    a.1<-as.data.table(a.1)
    a.1$fill<-1 #To count number of mutations per gene
    a.1<-a.1[,list(N.MUTATIONS=length(fill)), by=c("Tumor_Sample_Barcode", "Hugo_Symbol")]  
    
    #Add to main file
    dummy.a<-rbind(dummy.a,a.1)
  }
  
  #Return
  dummy.list<-list(table.1=dummy.a, background.genes=SEQUENCED.GENES)
  return (dummy.list)
}

#Call function
dummy.output<-Function.process.maf.files(main.maf, additional.mafs)

#Write to output file
saveRDS(dummy.output, file=output.file)