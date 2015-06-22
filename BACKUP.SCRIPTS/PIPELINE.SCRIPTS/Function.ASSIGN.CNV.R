######Function.ASSIGN.CNV.R#######
#021615
#Assigns cnv-mapped genes through ANNOVAR to TCGA patients and constructs matrix
#PART OF: TCGA.CNV.MATRIX.CONSTRUCT.sh

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

internal.function<-function(hugos){
  
  if (grepl("(NM", hugos, fixed=T)==T){
    hugos<-as.vector(sapply(unlist(strsplit(hugos, "),")), function(y) unlist(strsplit(y, "\\(NM"))[1]))
  } else{
    hugos<-unlist(strsplit(hugos, ","))
  }
  
  return(list(Hugo_Symbol=hugos))
}

Function.avinput.Process<-function(avinput.file){
  
  #Load file
  avinput<-fread(avinput.file, sep="\t",header=F, stringsAsFactors=F)
  
  #Filter for regions that fall on genes
  avinput<-avinput[V1 %in% c("intronic","exonic", "ncRNA_exonic","UTR3","UTR5","splicing","ncRNA_intronic","UTR5;UTR3"),]
  
  #Clean table
  avinput$V6<-NULL
  avinput$V7<-NULL
  setkey(avinput)
  avinput<-unique(avinput)
  
  #Apply hugo parsing function
  avinput<-avinput[,internal.function(V2), by=c("V1","V2","V3","V4","V5")]
  
  #Filter
  avinput<-avinput[,c("V3","V4","V5","Hugo_Symbol"),with=F]
  setkey(avinput)
  avinput<-unique(avinput)
  
  #Clean and Return
  setnames(avinput, c("Chromosome","Start","End","Hugo_Symbol"))
  return(avinput)
}

Function.main<-function(cnv.hugomap, cnv.sample.map, cnv.folder){
  
  #Load cnv tcga sample map
  cnv.map<-fread(cnv.sample.map, header=T, sep="\t", stringsAsFactors=F)
  setnames(cnv.map, c("file", "sample"))
  
  #Clean tumor sample labels to match other data in project
  cnv.map$SAMPLE<-sapply(cnv.map$sample, function(x) paste(unlist(strsplit(x, "-"))[1:4], collapse=".") )
  cnv.map$sample<-NULL
  
  #Filter for nocnv_hg19 (normalized hg 19 data)
  cnv.map<-cnv.map[grepl("nocnv_hg19", file),]
  
  #Prepping parallelization
  print ("prepping for parallelization")
  
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("cnv.map", "cnv.hugomap","as.data.table","cnv.folder","data.table","fread","setnames") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  main.list<-parLapply(cl,1:nrow(cnv.map), function(x) {
    
    #Extract file and sample
    file.name<-cnv.map$file[x]
    SAMPLE<-cnv.map$SAMPLE[x]
    
    #add path to file name
    file.name<-paste(cnv.folder, file.name, sep="/")
    
    #load file
    cnv.file<-fread(file.name, header=T, sep="\t", stringsAsFactors=F, drop=c(1,5))
    
    #Assign hugo info
    cnv.file<-merge(cnv.file, cnv.hugomap, by=c("Chromosome","Start","End"))
    
    #There is a chance that we have multiple cnv reads per hugo, use median
    cnv.file<-cnv.file[,list(CNV=median(Segment_Mean)), by="Hugo_Symbol"]
    
    #Add patient info as column name (for later merging) [Hugo_Symbol, CNV, SAMPLE]
    cnv.file$SAMPLE<-SAMPLE
    
    #Clean and Return
    return(cnv.file)
    
  } )
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done integrating hugos to cnv data")
  
  #Table Constructiong
  print ("Constructing table")
  main.table<-do.call(rbind, main.list)
  print ("Done building table")
  
  #Obtain normal and cancer patient lists
  all.patients<-unique(as.vector(main.table$SAMPLE))
  normal.patients<-all.patients[grepl("11A", all.patients)]
  cancer.patients<-all.patients[grepl("01A", all.patients)]
  
  #Return as list
  return(list(main.table=main.table, normal.patients=normal.patients, cancer.patients=cancer.patients))
  
}
##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
avinput.file<-args[1]
cnv.sample.map<-args[2]
cnv.folder<-args[3]
output.file<-args[4]
print("Files loaded")
##########################################

##################EXECUTE#################
cnv.hugomap<-Function.avinput.Process(avinput.file)
print ("Processed avinput file")

cnv.mapped<-Function.main(cnv.hugomap, cnv.sample.map, cnv.folder)
print ("Hugo CNV mapped")

###############WRITING OUTPUT############
saveRDS(cnv.mapped, output.file)