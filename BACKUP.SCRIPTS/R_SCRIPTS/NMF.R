#LEARNING NON_NEGATIVE MATRIX FACTORIZATION
#08/26/13

X<-matrix(1:12, 3, 4)
X

library(NMFN)
Z.MM=nnmf(X, 10,method=multiplicative update) #K is the number of factors/components in W and H matrices
Z.MM
Z.MM$W %*% Z.MM$H

library(NMF)
data(esGolub)
esGolub
esGolub<-esGolub[1:200,]
esGolub$Sample<-NULL

RES<-nmf(esGolub,3)
RES
fit(RES)