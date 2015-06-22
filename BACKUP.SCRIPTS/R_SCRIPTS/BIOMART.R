#LEARNING BIOMART - REVIEWED_UNREVIEWED UNIPROT DATASET
#Open file
uniprot.gene=read.table("uniprot-organism%3A%22Homo+sapiens+%5B9606%5D%22.tab",sep="\t",header=TRUE, quote="",comment.char="",)
head(uniprot.gene)
nrow(uniprot.gene)

#Filter for reviewed entries only - OPTIONAL
uniprot.gene.reviewed=uniprot.gene[grep("^reviewed", uniprot.gene$Status),]
nrow(uniprot.gene.reviewed)
uniprot.to.gene.reviewed=uniprot.gene.reviewed[, c(1,5)]
write.table(x=uniprot.to.gene.reviewed, "uniprot_to_gene_reviewed.text", sep="\t", row.names=F,quote=F, col.names=F)

#Get uniprot and gene columns
uniprot.to.gene=uniprot.gene[, c(1,5)]
head(uniprot.to.gene)
write.table(x=uniprot.to.gene, "uniprot_to_gene.text", sep="\t", row.names=F,quote=F, col.names=F)

#BIOMART
library("biomaRt")
listMarts()
#Query ensembl biomart database
ensembl=useMart("ensembl")
#list datasets in ensembl database
listDatasets(ensembl)
listDatasets(ensembl)[,1]
#select homo sapiens (hsapiens) dataset
ensembl=useDataset(dataset="hsapiens_gene_ensembl",mart=ensembl)
ensembl
#To build a biomart query (getBM), we need to introduce three argumentes: filters, attributes, values
#Filters are restrictions on the query. To see filters that can be used in selected dataset
filters=listFilters(ensembl)
head(filters)
filters[grep("hgnc", filters$name),]
#Attributes are what we are intereted in retrieving. To see attributes that can be used with dataset
attributes=listAttributes(ensembl)
head(attributes)
attributes[grep("name",attributes$name),]
#getBM has four main arguments: attributes, filters, values, mart

#Example - genes to uniprot (hgnc_symbol as filter, uniprot_swissprot_accession as attribute)
GENES=read.table("GENES")
RESULTS=getBM(attributes=c("hgnc_symbol", "uniprot_swissprot_accession"), filters="hgnc_symbol", 
              values=GENES[,1], mart=ensembl)
nrow(RESULTS)
nrow(RESULTS[grep("^$", RESULTS$uniprot_swissprot_accession),])
RESULTS[grep("^NA$", RESULTS$uniprot_swissprot_accession),]
RESULTS[grep("TRIM6-TRIM34", RESULTS$hgnc_symbol),]

#SEARCH FOR CATALYTIC SITE TERM IN ATTRIBUTES
attributes[grep("Protein",attributes$description),]

PDB=getBM(attributes=c("uniprot_swissprot_accession", "pdb"), filters="uniprot_swissprot_accession", 
              values=c("P69905"), mart=ensembl)
PDB

testing.table=read.table("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/TEST")
testing.table
UNIPROT=PDB=getBM(attributes=c("pdb","uniprot_swissprot_accession"),
                  filters="pdb", values=testing.table[,1], mart=ensembl)
UNIPROT

PDB.SEQ=getSequence(id=c("1A0N"), type="pdb", seqType="peptide", mart=ensembl)
PDB.SEQ
nrow(PDB.SEQ)
PDB.SEQ$peptide[1]