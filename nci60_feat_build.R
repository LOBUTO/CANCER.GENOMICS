#Function nci60_feat_build.R
#Function to develop nci60 feat table for classfication methods

library(data.table)
library(reshape2)

Function.nci60.hmdb.table <- function(nci.file, hmdb, pre.table=data.table()){

  nci60 <- fread(nci.file, header=T)
  nci60$DRUG <- as.character(nci60$DRUG)

  if(nrow(pre.table)>0){
    nci60 <- nci60[!DRUG %in% unique(pre.table$NSC),]
  }

  mets <- intersect(unique(hmdb$METABOLITE), unique(nci60$METABOLITE))

  nci60 <- nci60[METABOLITE %in% mets,]
  nci60[,SD:=sd(TC), by="DRUG"]
  nci60 <- nci60[SD!=0,]

  hmdb <- hmdb[METABOLITE %in% mets,]

  met.matrix <- cbind(acast(nci60, METABOLITE~DRUG, value.var="TC")[mets, ],
                      acast(hmdb,  METABOLITE~DRUG, value.var="TC")[mets, ])
  print (dim(met.matrix))
  nci.drugs <- unique(nci60$DRUG)
  cgp.drugs <- unique(hmdb$DRUG)

  main.table <- data.table()
  for (i in nci.drugs){
    print(i)
    cors <- sapply(cgp.drugs, function(x) cor(met.matrix[,i], met.matrix[,x]) )
    cors.table <- data.table(NSC=i, CGP=cgp.drugs, COR=cors)

    main.table <- rbind(main.table, cors.table)
  }

  #Clean up and return
  if (nrow(pre.table)>0){
    main.table <- rbind(main.table, pre.table)
  }

  return(main.table)
}

Function.build.nci60.feat.class <- function(nci60.exp, nci.gi50, nci60.cgp.drugs, cgp.cor.AUC, cgp.exp.matrix, feat.scaling=T){

  #Build co-expression table
  cgp.exp.matrix <- cgp.exp.matrix[, intersect(colnames(cgp.exp.matrix), colnames(cgp.cor.AUC))]
  common.genes <- intersect(rownames(cgp.exp.matrix), rownames(nci60.exp))
  main.exp <- cbind(cgp.exp.matrix[common.genes,], nci60.exp[common.genes,])
  main.exp <- cor(main.exp, method = "pearson")

  main.exp <- main.exp[colnames(nci60.exp), colnames(cgp.exp.matrix)]
  main.exp <- data.table(main.exp, keep.rownames = T)
  setnames(main.exp, c("cell_name", colnames(main.exp)[2:ncol(main.exp)]))

  #Build co-chemical feature table
  main.chem <- nci60.cgp.drugs[CGP %in% colnames(cgp.cor.AUC),]
  main.chem <- acast(main.chem, NSC~CGP, value.var="COR")
  main.chem <- data.table(main.chem, keep.rownames = T)
  setnames(main.chem, c("Compound", colnames(main.chem)[2:ncol(main.chem)]))

  #Build activity table and classify
  main.act <- nci.gi50[,c("CELL", "NSC", "SCALE.ACT"),with=F]
  main.act$CLASS <- ifelse(main.act$SCALE.ACT>=0.8, 1,
                           ifelse(main.act$SCALE.ACT<=-0.8, 0, 2))
  main.act <- main.act[CLASS!=2,]
  main.act$SCALE.ACT <- NULL
  setnames(main.act, c("cell_name", "Compound", "NORM_AUC"))

  #Merge all
  main.table <- merge(main.act, main.exp, by="cell_name")#, allow.cartesian=TRUE)
  main.table <- merge(main.table, main.chem, by="Compound")#, allow.cartesian=TRUE)

  #Clean up and return
  main.table <- main.table[,colnames(cgp.cor.AUC), with=F]
  if (feat.scaling==T){
    main.table <- Function.scale.data.table(main.table, col.protect=1:3)
  }
  return(main.table)
}

########################################################################################################################################

OBJ_FOLDER <- "~/Documents/FOLDER/OBJECTS/"
TABLES_FOLDER <- "~/Documents/FOLDER/TABLES/TCGA.TRAINING/"

nci60.cgp.drugs <- readRDS(paste0(OBJ_FOLDER, "052316.NCI60.HMDB.TABLE.643711.rds")) #UPDATE as you process
nci60.exp <- readRDS(paste0(OBJ_FOLDER, "061916.NCI60.EXP.rds"))
nci.gi50 <- readRDS(paste0(OBJ_FOLDER, "061916.NCI.GI50.rds"))

DRUGS.MET.PROFILE <- readRDS(paste0(OBJ_FOLDER, "061916.DRUGS.MET.PROFILE.rds"))

cgp.cor.AUC <- readRDS(paste0(OBJ_FOLDER, "052516.CGP.COR.AUC.rds"))
cgp.exp.matrix <- readRDS(paste0(OBJ_FOLDER, "061716.CGP.EXP.MATRIX.rds"))

nci60.cgp.drugs <- Function.nci60.hmdb.table(paste0(TABLES_FOLDER,"NCI60.TC.HMDB.FP4"),
                                             DRUGS.MET.PROFILE[DRUG %in% colnames(cgp.cor.AUC),],
                                             nci60.cgp.drugs)
saveRDS(nci60.cgp.drugs, paste0(OBJ_FOLDER, "052316.NCI60.HMDB.TABLE.PRELAST.rds"))

#Unscaled features
nci.cgp.feat.class <- Function.build.nci60.feat.class(nci60.exp, nci.gi50, nci60.cgp.drugs, cgp.cor.AUC, cgp.exp.matrix, feat.scaling=F)

########################################################################################################################################
#Store table
saveRDS(nci.cgp.feat.class, paste0(OBJ_FOLDER, "061816.NCI.CGP.FEAT.CLASS.UNSCALED_FEAT.rds"))

print ("Done!!")
