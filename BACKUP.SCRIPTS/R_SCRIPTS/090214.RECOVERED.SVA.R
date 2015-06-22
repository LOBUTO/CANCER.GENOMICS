
Function.RNAseq.Differential.Expression<-function(normal.matrix, cancer.matrix) {
  #Performs differential expression between the two matrices
  #Produces topTable
  
  require(plyr)
  require(limma)
  require(edgeR)
  
  #Build design matrix
  G1.n.samples<-length(colnames(cancer.matrix))
  G0.n.samples<-length(colnames(normal.matrix))
  G.design.matrix<-data.frame(G=c(rep("G1", G1.n.samples), rep("G0", G0.n.samples)))
  G.design.matrix<-model.matrix(~G, G.design.matrix)
  
  #Combine matrices
  cancer.matrix$rn<-rownames(cancer.matrix)
  normal.matrix$rn<-rownames(normal.matrix)
  dummy.expression.matrix<-join_all(list(as.data.frame(cancer.matrix), as.data.frame(normal.matrix)), by="rn", type="inner")
  rownames(dummy.expression.matrix)<-dummy.expression.matrix$rn
  dummy.expression.matrix$rn<-NULL #Remove column used to combine data frames
  dummy.expression.matrix<-dummy.expression.matrix[complete.cases(dummy.expression.matrix),] #Remove NAs
  
  #Convert RNAseq counts to log-count per million and normalize
  G.all<-as.matrix(dummy.expression.matrix)
  G.isexpr<- rowSums(cpm(G.all)>1) >= 20 #Keep genes with at least 1 count-per-million reads (cpm) in at least 20 samples
  G.all<-G.all[G.isexpr,]
  
  G.all<-DGEList(counts=G.all) #For scale normalization
  G.all<-calcNormFactors(G.all) #TMM - Scale normalization #KEEP IN MIND THAT THIS MAY NOT BE NECESSARY AS RNASEQ V2 files may already be
  # upper quantile normalized (TCGA)
  
  G.all<-voom(G.all, G.design.matrix,normalize.method="quantile") #Convert RNAseq expression values to log-count per million + quantile normalization
  
  #Do differential expression
  G.fit = lmFit(G.all, G.design.matrix) #fitting data 
  G.eb = eBayes(G.fit)
  
  #Get topTable
  all.G.fit<-topTable(G.eb, coef=2, n=Inf)
  
  #Return topTable
  return(all.G.fit)
}


######
cancer.matrices$normal
cancer.matrices$tumor

shared.genes = intersect(rownames(cancer.matrices$normal), rownames(cancer.matrices$tumor))
normal = cancer.matrices$normal[shared.genes, ]
tumor = cancer.matrices$tumor[shared.genes, ]

num.each = 100

colnames(normal) = paste("Norm", colnames(normal))
colnames(tumor) = paste("Tumor", colnames(tumor))
matrix.combined = as.matrix(cbind(normal[, sample(ncol(normal), num.each)], tumor[, sample(ncol(tumor), num.each)]))

matrix.full = as.matrix(cbind(normal, tumor))

matrix.combined = matrix.combined[rowSums(matrix.combined) >= 100, ]

design = data.frame(type=rep(c("Normal", "Tumor"), each=num.each))
mod = model.matrix(~ type, design)
mod0 = model.matrix(~ 1, design)

ss = svaseq(matrix.combined, mod, mod0, n.sv=16)

# design.full = data.frame(type=rep(c("Normal", "Tumor"), c(ncol(normal), ncol(tumor))))
# mod.full = model.matrix(~ type, design.full)
# mod0.full = model.matrix(~ 1, design.full)
# full.ss = svaseq(matrix.full, mod.full, mod0.full, n.sv=4)

cors = cor(matrix.combined, method="spearman")

v = voom(matrix.combined, normalize.method="quantile")
f = lmFit(v, model.matrix(~ type, design))
eb = eBayes(f)
hist(eb$p.value[, 2])

sum(p.adjust(eb$p.value[, 2]) < .05)

q = qvalue(eb$p.value[, 2])
sum(q$qvalue < .05)

f.ss = lmFit(v, model.matrix(~ type + ss$sv, design)) #Just like regular fit but added batch effect normalization
eb.ss = eBayes(f.ss)
head(topTable(eb.ss, coef=2, n=Inf))
hist(eb.ss$p.value[, 2])
sum(p.adjust(eb.ss$p.value[, 2]) < .05)

# remove effect of 3 surrogate variables
heatmap(cor(v$E))

ss.removed = v$E - f.ss$coefficients[, -c(1, 2)] %*% t(f.ss$design[, -c(1, 2)])

heatmap(cor(v$E), scale="none")
heatmap(cor(ss.removed), scale="none")

v.p.table