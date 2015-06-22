read.expression.files = function(path=".", n=50) {
  #There are "null" in table that function reads as NA instead to deal with easier
  tables = lapply(list.files(path)[1:n], function(f) read.csv(paste0(path, "/", f), header=FALSE, sep="\t", skip=2, na.strings=c("null"),
                                                            col.names=c("gene", strsplit(f, "\\.")[[1]][1])))
  
  m = tables[[1]]
  for (i in 2:length(tables)) {
    m = merge(m, tables[[i]], all=FALSE) #So that only common genes in tables are merged
  }
  rownames(m) = m$gene
  m$gene = NULL
  return (as.matrix(m))
}


n = 50
breast.m = read.expression.files("BRCA/Level_3", n)
coad.m = read.expression.files("COAD/Level_3", n)

colnames(breast.m) = paste0(colnames(breast.m), ".breast")
colnames(coad.m) = paste0(colnames(coad.m), ".coad")

combined = cbind(breast.m, coad.m[rownames(breast.m), ])
#Check if there are NAs

co = cor(combined, use="complete.obs") #complete.obs is to delete cases where comparisson includes NAs
heatmap(co)
#m = merge(breast.m, coad.m)

library(limma)
d = data.frame(cancer=rep(c("breast", "coad"), each=n))

library(GSEAL)
library(org.Hs.eg.db)

combined.filtered = combined[rownames(combined) %in% names(as.list(org.Hs.egSYMBOL2EG)), ] #Keep only those that have an entrez gene identifier
rownames(combined.filtered) = as.character(org.Hs.egSYMBOL2EG)[rownames(combined.filtered)] #change names in actual matrix
combined.filtered = combined.filtered[complete.cases(combined.filtered), ] #Remove NAs

boxplot(colSums(combined.filtered) ~d$cancer) #each dot (column) is a patient, so in unnormalized data the expression data cummulative sum of the expression data per patient is higher in breast than colorectal, so need to be normalized

norm = normalizeBetweenArrays(combined.filtered, method="quantile") #normalizes expression intensities across set of arrays
fit = lmFit(norm, model.matrix(~ cancer, d)) # model.matrix(~cancer, d) inputs a design matrix so that each gene row in "norm" is split 50/50 between COAD and BRCA, keep in mind that this is because we originally let n=50, otherwise we have to create design matrix according to colnames(combined.filtered)
eb = eBayes(fit)
hist(eb$p.value[, 2])
volcanoplot(eb, highlight=20)

y = eb$coefficients[, 2]
hist(y)

mm = MembershipMatrix(organism = "Hs.eg", ontology = "BP", min.size = 5, max.size = 500, ancestors=FALSE)

wilcoxon.mm = TestEnrichment(mm, rownames(norm), y, method = "wilcoxon")
ttest.mm = TestEnrichment(mm, rownames(norm), y, method="t.test")

CompareTopSets(ttest.mm)
CompareTopSets(wilcoxon.mm)

