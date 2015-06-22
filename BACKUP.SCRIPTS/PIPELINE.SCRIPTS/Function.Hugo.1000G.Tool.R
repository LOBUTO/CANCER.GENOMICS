#Function.Hugo.1000G.Tool
#Assign Hugo information to processed thousand genome position
#111414

Function.prep<-function(thousand.files, exon.file){
  require(data.table)
  
  #Parse 1000G files
  thousand<-fread(thousand.files[1],header=T,sep="\t",stringsAsFactors=F,)
  for (record in thousand.files[1:length(thousand.files)]) {
    record.fread<-fread(record,header=T,sep="\t",stringsAsFactors=F)
    thousand<-cbind(thousand,record.fread)
  }
  
  #Parse Exon file
  EXONS<-fread(exon.file, header=T, sep="\t", stringsAsFactors=F, drop=c(3,6:8))  
  EXONS<-unique(EXONS)
  
  #Return
  return(list(thousand=thousand, exons=EXONS))
  
}

Function.main<-function(thousand.table, exons.table){
  require(data.table)
  require(IRanges)
  
  Chromosomes<-unique(as.vector(thousand.table$Chrom))
  
  #Apply function
  MAPPED.LISTS<-lapply(Chromosomes, function(x) {
    
    print (x)
    
    #Separate tables per chromosome
    thousand.chrm<-thousand.table[Chrom==x,]
    exon.chrm<-exons.table[Chrom==x,]
    
    #Create ranges objects
    thousand.ranges<-with(thousand.chrm, IRanges(Position, width=1,NULL))
    exon.ranges<-with(exon.chrm, IRanges(START, END, names=Hugo_Symbol))
    
    #Find coordinate overlaps to identify genes in 1000G
    coordinate.overlaps<-findOverlaps(thousand.ranges, exon.ranges)
    thousand.mapped<-cbind(thousand.chrm[queryHits(coordinate.overlaps),], exon.chrm[subjectHits(coordinate.overlaps),])
    
    #Return mapped chromosome
    return(thousand.mapped)
  })
  
  MAPPED<-do.call(rbind, MAPPED.LISTS)
  
  #Return
  return(MAPPED)
}

thousand.files<-c("PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.SNP.CNV", "PIPELINES/METABOLIC.DRIVERS/TABLES/111014.THOUSAND.Y.SNP.CNV")
exon.file<-"PIPELINES/METABOLIC.DRIVERS/TABLES/110914.EXON.COORDIANTES"

thousand.mapped<-Function.main(thousand, unique(EXONS[,c(1,2,4,5),with=F]) )
thousand.mapped[Hugo_Symbol=="TP53",]