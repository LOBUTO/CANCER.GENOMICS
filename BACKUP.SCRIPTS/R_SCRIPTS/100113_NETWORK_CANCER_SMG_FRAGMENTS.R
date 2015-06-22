#NETWORK_CANCER_SMG_FRAGMENTS
#10/01/13

normalit<-function(m){
  A=(m - min(m))/(max(m)-min(m))
  return (A)
}

#LOAD FILE
FILE=read.table("NETWORK/R_ANALYSIS/100113_CANCER_SMG_BLOCKS", header=TRUE)

#Remove first column and assign to row labels
DATA=FILE[,-1]; DATA
row.names(DATA)<-FILE[,1]

#Convert to matrix for heatmap
DATA_MATRIX=as.matrix(DATA)
DATA_MATRIX=normalit(DATA_MATRIX); DATA_MATRIX

#Obtain heatmap
library(RColorBrewer)
hmcol<-brewer.pal(9,"Blues")
heatmap(DATA_MATRIX, col=hmcol)
