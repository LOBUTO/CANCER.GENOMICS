#031814
#Get position frequency of mature.mirna
source("http://bioconductor.org/biocLite.R")
biocLite("seqLogo")

mirna.pos.logo<-function(mirna.file, type) { #Takes files of the form 031814_BRCA_SEED_POS
  library(seqLogo)  
  SEED<-as.data.table(read.table(mirna.file, sep="\t", header=T,stringsAsFactors=F))
  
  #First make frequency
  GGPLOT<-ggplot(SEED, aes(x=factor(Position))) + geom_histogram() +
    theme(axis.text.y=element_text(size=rel(1.5)), plot.title=element_text(size=rel(2)), axis.text.x=element_text(size=rel(2.0)), 
          strip.text.x = element_text(size = 30) , legend.text = element_text(size = 22)) +
    labs(title=type, size=7)
  
  #Then information content logos
  SEED.COUNT<-SEED[,list(COUNT=length(Tumor_Sample_Barcode)), by=c("Position", "NT")]
  SEED.COUNT<-cast(SEED.COUNT, NT~Position,fill=0)
  rownames(SEED.COUNT)<-SEED.COUNT$NT
  SEED.COUNT$NT<-NULL
  SEED.COUNT<-sweep(SEED.COUNT, 2,colSums(SEED.COUNT), FUN="/")

  LOGOS<-seqLogo(makePWM(SEED.COUNT),xaxis=F)
  
  RETURN<-list(GGPLOT, LOGOS)
  return(RETURN)
}

#Run on cancers that have large population of matures - Due to logo, run one by one and save for now
TEST<-mirna.pos.logo("DATABASES/CANCER_DATA/miRNA/031814_UCEC_SEED_POS", "UCEC")
TEST[[1]]

