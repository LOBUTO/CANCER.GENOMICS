#031314 - Obtain coordinates for pri- and mature mirna

#Load mirna coordinates file
mirna.coor<-read.table("DATABASES/CANCER_DATA/miRNA/hsa._coordinatesgff3.gff3", header=F, sep="\t", skip=13,stringsAsFactors=F)
head(mirna.coor)
mirna.coor.prim<-mirna.coor[mirna.coor$V3=="miRNA_primary_transcript",]
mirna.coor.prim$Synonyms<-sapply(mirna.coor.prim$V9, function(x)  strsplit(strsplit(x,";")[[1]][3],"=")[[1]][2] )
head(mirna.coor.prim)

#Mature transcripts
mirna.coor.mature<-mirna.coor[mirna.coor$V3=="miRNA",]
mirna.coor.mature$ID<-sapply(mirna.coor.mature$V9, function(x)  strsplit(strsplit(x,";")[[1]][4],"=")[[1]][2] )
mirna.coor.mature$name<-sapply(mirna.coor.mature$V9, function(x)  strsplit(strsplit(x,";")[[1]][3],"=")[[1]][2] )
mirna.coor.mature<-mirna.coor.mature[,c(1,4,5,7,10,11)]
colnames(mirna.coor.mature)<-c("Chromosome","Start","End","Strand","Accession","Name")
head(mirna.coor.mature)
write.table(file="DATABASES/CANCER_DATA/miRNA/031614_mature_coor", mirna.coor.mature, quote=F, row.names=F,sep="\t")

#Load file to map name (i.e. hsa-mir-7-2) to Approved Symbol (MIR7-2)
mirna.map<-read.table("DATABASES/CANCER_DATA/miRNA/RNA_micro.txt.txt", header=T, sep="\t", stringsAsFactors=F)
mirna.map<-mirna.map[,c(2,7)]
head(mirna.map)
#Edited to complete non-identied entries in .maf - 03/17/14
mirna.map<-rbind(mirna.map, c("MIR124-1", "hsa-mir-124-1"))
mirna.map<-rbind(mirna.map, c("MIR1269A", "hsa-mir-1269a"))
mirna.map<-rbind(mirna.map, c("MIR128-2", "hsa-mir-128-2"))
mirna.map<-rbind(mirna.map, c("MIR27B", "hsa-mir-27b"))
mirna.map<-rbind(mirna.map, c("MIR301A", "hsa-mir-301a"))
mirna.map<-rbind(mirna.map, c("MIR3130-1", "hsa-mir-3130-1"))
mirna.map<-rbind(mirna.map, c("MIR450A1", "hsa-mir-450a-1"))
mirna.map<-rbind(mirna.map, c("MIR450A2", "hsa-mir-450a-2"))
mirna.map<-rbind(mirna.map, c("MIR490", "hsa-mir-490"))
mirna.map<-rbind(mirna.map, c("MIR513A1", "hsa-mir-513a-1"))
mirna.map<-rbind(mirna.map, c("MIR516A2", "hsa-mir-516a-2"))
mirna.map<-rbind(mirna.map, c("MIR550A1", "hsa-mir-550a-1"))
mirna.map<-rbind(mirna.map, c("MIR92A1", "hsa-mir-92a-1"))

#Merge coordinate to map
mirna.pri<-merge(mirna.map, mirna.coor.prim, by="Synonyms")
mirna.pri$ID<-sapply(mirna.pri$V9, function(x)  strsplit(strsplit(x,";")[[1]][1],"=")[[1]][2] )
mirna.pri<-mirna.pri[,c(1,2,3,6,7,9,12)]
colnames(mirna.pri)<-c("Synonyms", "Symbol","Chromosome","Start","End","Strand","Accession")
head(mirna.pri)
dim(mirna.pri)
write.table(file="DATABASES/CANCER_DATA/miRNA/031614_primirna_coor", mirna.pri, quote=F, row.names=F,sep="\t")

#Assign Symbol of Parent to mature miRs - SOME MAY NOT HAVE AN HGNC SYMBOL
head(mirna.coor.mature)
dim(mirna.coor.mature)
mirna.coor.mature.parent<-merge(mirna.coor.mature, mirna.pri[,c(7,2)], by="Accession")
colnames(mirna.coor.mature.parent)[7]<-"Parent"
head(mirna.coor.mature.parent)

