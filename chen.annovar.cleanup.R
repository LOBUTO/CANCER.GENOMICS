#chen.annovar.cleanup.R
#010815
#Cleaning up of annovar processing of chen file
#STEP 3 of CHEN.REPLICATION.GENE.ASSIGN.sh

library(data.table)

args<-commandArgs(trailingOnly=T) 

file.in<-args[1] #This is the chen.rt.avinput.variant_function file
output.file<-args[2]


########FUNCTIONS#######
#Extraction function# - Also filters for far off intergenic (threshold at 500000)
internal.function<-function(V1,V2.part, threshold=500000){  
  
  if ((V1=="UTR5") | (V1=="UTR3") | (V1=="splicing") ){
    Hugos<-unlist(strsplit(V2.part,"\\("))[1]
    
  } else if ( (V1=="exonic;splicing") | (V1=="UTR5;UTR3") ){
    Hugos<-unlist(strsplit(V2.part,";"))
    Hugos<-unique(as.vector(sapply(Hugos, function(x) unlist(strsplit(x, "\\(NM"))[1])))
    
  } else if (V1=="ncRNA_splicing"){
    Hugos<-unlist(strsplit(V2.part,"),"))
    Hugos<-as.vector(sapply(Hugos, function(x)  unlist(strsplit(x, "\\("))[1] ))
    
  } else if (V1=="intergenic"){ #Thresholding at a particular distance
    dist.hugos<-unlist(strsplit(V2.part,",")) 
    dist.hugos<-sapply(dist.hugos, function(x) {
      dist<-unlist(strsplit(x, '[=)]'))[2]
      
      if (dist=="NONE"){ #Filter for "NONE" hugo/distance annotations
        Hugo<-"filter"
      }else if (as.numeric(dist)<=threshold){
        Hugo<-unlist(strsplit(x,"\\("))[1]
      } else{
        Hugo<-"filter" #If doesn't pass threshold then return "filter" for later filtering
      }
      return(Hugo)  
    })
    Hugos<-dist.hugos
    
  } else if (grepl("),", V2.part) ){
    Hugos<-unlist(strsplit(V2.part,","))
    Hugos<-as.vector(sapply(Hugos, function(x)  unlist(strsplit(x, "\\("))[1] ))
    
  } else if ( (grepl(",",V2.part)) | (grepl(";",V2.part)) ){
    Hugos<-unlist(strsplit(V2.part,"[,;]"))
    
  } else {
    Hugos<-V2.part
  }
  
  return(list(Hugos=Hugos))
}

#Filterinf function#
internal.rt.rank<-function(V1.part, V8.part){
  
  #Build table
  main<-data.table(TYPE=V1.part, RT=V8.part)
  
  #Check by rank
  if (nrow(main[TYPE %in% c("exonic","intronic"),])>0){
    main.rt<-main[TYPE %in% c("exonic","intronic"),]
    RT<-median(as.vector(main.rt$RT))
    
  } else if (nrow(main[TYPE %in% c("ncRNA_exonic","ncRNA_intronic"),])>0){
    main.rt<-main[TYPE %in% c("ncRNA_exonic","ncRNA_intronic"),]
    RT<-median(as.vector(main.rt$RT))
    
  } else if (nrow(main[TYPE %in% c("UTR5","UTR3","UTR5;UTR3"),])>0){
    main.rt<-main[TYPE %in% c("UTR5","UTR3"),]
    RT<-median(as.vector(main.rt$RT))
    
  } else if (nrow(main[TYPE %in% c("exonic;splicing","splicing"),])>0){
    main.rt<-main[TYPE %in% c("exonic;splicing","splicing"),]
    RT<-median(as.vector(main.rt$RT))
  
  } else if (nrow(main[TYPE %in% c("ncRNA_splicing"),])>0){
    main.rt<-main[TYPE %in% c("ncRNA_splicing"),]
    RT<-median(as.vector(main.rt$RT))
    
  } else if (nrow(main[TYPE %in% c("downstream","upstream", "upstream;downstream"),])>0){
    main.rt<-main[TYPE %in% c("downstream","upstream", "upstream;downstream"),]
    RT<-median(as.vector(main.rt$RT))
    
  } else if (nrow(main[TYPE %in% c("intergenic"),])>0){
    main.rt<-main[TYPE %in% c("intergenic"),]
    RT<-median(as.vector(main.rt$RT))
  }
  
  return(list(RT=RT))
}

########Load file#######
chen.annovar<-fread(file.in, header=F, sep="\t", stringsAsFactors=F, drop=5:7)

#########Find best replication time candidate per gene#######
#Obtain RT score per gene - NOTE: OR WE CAN ADJUST BY POSITION BY EXTRAPOLATING (ie. decay function, future) 
chen.annovar<-chen.annovar[,internal.function(V1,V2), by=c("V1","V2","V3","V4","V8")]
chen.annovar<-chen.annovar[Hugos!="filter",]

#Rank: exonic > intronic > ncRNA_exonic == ncRNA_intronice > UTR5 == UTR3 > splicing > upstream == downstream >> intergenic 
chen.annovar<-chen.annovar[,internal.rt.rank(V1,V8), by=c("Hugos", "V3")]

#Classify into replication category: high, medium, low
chen.annovar$cut<-cut(chen.annovar$RT, breaks=c(0, as.vector(quantile(chen.annovar$RT, c(0.33,0.66))),1),include.lowest=T,
                      labels=c("low", "medium","high"))

#Label
setnames(chen.annovar, c("Hugo_Symbol", "Chrom","REP.TIME", "REP.CLASS"))

########Write to output file########
write.table(file=output.file, chen.annovar, sep="\t", quote=F, row.names=F, col.names=T)
cat ("done annotating replication times")

#####################################TEST USAGE#######################################################################################
# test<-fread("FOLDER/SCRIPTS/desprat.rt.avinput.variant_function", header=F, sep="\t",stringsAsFactors=F,drop=5:7)
#   
# test<-test[,internal.function(V1,V2), by=c("V1","V2","V3","V4","V8")]
# test<-test[Hugos!="filter",]
# 
# test<-test[,internal.rt.rank(V1,V8), by=c("Hugos","V3")]
# 
# test$cuts<-cut(test$RT, breaks=c(0, as.vector(quantile(test$RT, c(0.33,0.66))),1),include.lowest=T, labels=c("low", "medium","high"))
# table(test$cuts)
# 
# ggplot(test, aes(RT, fill=V3))+geom_histogram() +theme.format + facet_wrap(~V3)
# ggplot(test, aes(RT, colour=cuts)) + geom_histogram() + theme.format + facet_wrap(~cuts)
# 
# rm(test)





