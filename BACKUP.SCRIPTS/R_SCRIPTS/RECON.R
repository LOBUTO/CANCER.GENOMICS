###RECON TESTING###
library(SBMLR)
mock.sbml<-readSBML("DATABASES/RECON/recon2_02/recon2.v02.xml")
mock.sbml$species[["_6523_1_g"]]
mock.results<-sim(mock.sbml, seq(1,100,1))
dim(mock.results)
mock.results[1:3,1:3]
plot(mock.results[,"time"], mock.results[,"_6523_1_g"])
plot(mock.results[,"time"], mock.results[,"M_acac_c"])
table(colMeans(mock.results))

mock.results[1:3,17230:17239]