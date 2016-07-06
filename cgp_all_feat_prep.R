#cgp_all_feat.R
#Function to build feature matrix for cgp based on AUC values using all available features

#####################################################################################################################
#Functions
library(data.table)
library(reshape2)

Function.scale.data.table<-function(d.t, col.protect=1:3, criteria=c(), scaling="SCALE"){
  #Scaling can be regular z-scoring ("SCALE") or 0-1 range ("ZERO")

  y.colnames <- colnames(d.t[,-col.protect,with=F])

  if (length(criteria)>0){
    all.criteria <- unique(d.t[[criteria]])

    all.y <- lapply(all.criteria, function(x)  {
      z <- d.t[d.t[[criteria]]==x,]
      z <- data.frame(z[,-col.protect,with=F])
      col.except <- apply(z, 2, sd)==0
      not.col.except <- apply(z, 2, sd)!=0

      if (scaling=="SCALE"){
        z <- cbind(z[,col.except], scale(z[,not.col.except]))
      } else if (scaling=="ZERO"){
        z <- cbind(z[,col.except], apply(z[,not.col.except], 2 ,function(zz) Function.range.0.1(zz)) )
      }

      return (z)
    })
    y <- do.call(rbind, all.y)

  } else {
    y<-data.frame(d.t[,-col.protect,with=F])

    if (scaling=="SCALE"){
      y<-scale(y)
    } else if (scaling=="ZERO"){
      y <- apply(y, 2, function(zz) Function.range.0.1(zz))
    }

  }

  x.colnames<-colnames(d.t[,col.protect,with=F])
  x<-cbind(d.t[,col.protect,with=F], y)
  setnames(x, c(x.colnames, y.colnames))

  return(x)
}

Function.pre.cgp.3 <- function (cgp.file, norm.resp="Compound") {
  #Function to pre-process cgp files
  #Obtains AUC and NORM.AUC of compounds across each cell

  #Load drug file
  #cgp.file <- "DATABASES/CANCERRXGENE/gdsc_manova_input_w5.csv"
  cgp.drug <- fread(cgp.file, header=T, sep=",")

  #Keep only AUC information
  cgp.ic50.col <- colnames(cgp.drug)[ (grepl("_AUC", colnames(cgp.drug), ignore.case = T)) ]

  cgp.drug <- cgp.drug[,c("Cell Line", cgp.ic50.col),with=F]
  cgp.drug <- cgp.drug[6:nrow(cgp.drug),]
  cgp.drug <- cgp.drug[, !duplicated(colnames(cgp.drug)),with=F]

  #Melt to long format and remove empty cells
  cgp.drug <- melt(cgp.drug, id.vars = "Cell Line")[value!="",]

  #Obtain normalized AUC values per cells
  setnames(cgp.drug, c("cell_name", "Compound", "AUC"))

  #cgp.drug[,COUNT:=length(unique(Compound)), by="cell_name"] #To remove low counting cell experiments
  #cgp.drug <- cgp.drug[COUNT>10,]
  #cgp.drug$COUNT <- NULL

  cgp.drug$AUC <- as.numeric(cgp.drug$AUC)
  if (norm.resp!="NONE"){
    cgp.drug <- cgp.drug[,NORM.AUC:=scale(AUC), by=norm.resp]
  } else if (norm.resp=="NONE") {
    cgp.drug$NORM.AUC <- cgp.drug$AUC
  }

  #Clean up
  cgp.drug$Compound <- as.character(cgp.drug$Compound)
  cgp.drug$Compound <- sapply(cgp.drug$Compound, function(d) strsplit(d, split = "_")[[1]][1])

  #Return
  return(cgp.drug)
}

Function.pre.cgp.2 <- function (cgp.file, norm.resp="Compound") {
  #Function to pre-process cgp files
  #Obtains pIC50 and normalized pIC50 of compounds across each cell or compounds across cells
  #norm.resp can either be "Compound" or "cell_name"

  #Load drug file
  cgp.file <- "DATABASES/CANCERRXGENE/gdsc_manova_input_w5.csv"
  cgp.drug <- fread(cgp.file, header=T, sep=",")

  #Keep only IC50 information
  cgp.ic50.col <- colnames(cgp.drug)[ (grepl("IC_50", colnames(cgp.drug), ignore.case = T)) &
                                        !grepl("IC_50_LOW", colnames(cgp.drug), ignore.case = T) &
                                        !grepl("IC_50_HIGH", colnames(cgp.drug), ignore.case = T)]

  cgp.drug <- cgp.drug[,c("Cell Line", cgp.ic50.col),with=F]
  cgp.drug <- cgp.drug[6:nrow(cgp.drug),]
  cgp.drug <- cgp.drug[, !duplicated(colnames(cgp.drug)),with=F]

  #Melt to long format and remove empty cells
  cgp.drug <- melt(cgp.drug, id.vars = "Cell Line")[value!="",]

  #Obtain normalized IC50 values per cells
  setnames(cgp.drug, c("cell_name", "Compound", "IC50"))
  cgp.drug$IC50 <- as.numeric(cgp.drug$IC50)

  #cgp.drug$IC50 <- exp(as.numeric(cgp.drug$IC50)) #since values are obtained as natural log micromolar

  #cgp.drug[,COUNT:=length(unique(Compound)), by="cell_name"] #To remove low counting cell experiments
  #cgp.drug <- cgp.drug[COUNT>10,]
  #cgp.drug$COUNT <- NULL

  cgp.drug$pIC50 <- -log10(exp(cgp.drug$IC50)) #since values are obtained as natural log micromolar

  #Normalize with respect to choice
  cgp.drug <- cgp.drug[,NORM.pIC50:=scale(pIC50), by=norm.resp] #Normalizing by choice now (change on 011516)
  cgp.drug <- cgp.drug[,NORM.IC50:=scale(IC50), by=norm.resp]

  #Clean up
  cgp.drug$Compound <- as.character(cgp.drug$Compound)
  cgp.drug$Compound <- sapply(cgp.drug$Compound, function(d) strsplit(d, split = "_")[[1]][1])

  #Return
  return(cgp.drug)
}

Function_cgp_all_feat <- function(cgp.file,  DRUGS.MET.PROFILE, cosmic.expression, norm.resp="Compound",
                                  response="pIC50", feat.scaling=T, classify=T) {

  #Process target variable
  if (response=="pIC50" | response=="IC50"){
    act.table <- Function.pre.cgp.2(cgp.file, norm.resp = norm.resp)
  } else if (response=="AUC"){
    act.table <- Function.pre.cgp.3(cgp.file, norm.resp = norm.resp)
  }
  print("Done processing activity")

  #Process expression using cosmic data
  cosmic.cast <- acast(cosmic.expression, SAMPLE_NAME~GENE_NAME, value.var="Z_SCORE")
  cosmic.cast <- data.table(cosmic.cast, keep.rownames = T)
  setnames(cosmic.cast, c("cell_name", colnames(cosmic.cast)[2:ncol(cosmic.cast)]))
  print("Done casting expression")

  #Process chemical Tanimoto data
  DRUGS.MET.PROFILE.CGP <- DRUGS.MET.PROFILE[,1:3,with=F]
  DRUGS.MET.PROFILE.CGP <- DRUGS.MET.PROFILE.CGP[DRUG %in% unique(act.table$Compound),]

  drugs.prof <- acast(DRUGS.MET.PROFILE.CGP, DRUG~METABOLITE, value.var="TC")
  drugs.prof <- data.table(drugs.prof, keep.rownames = T)
  setnames(drugs.prof, c("DRUG", colnames(drugs.prof)[2:ncol(drugs.prof)]))
  print ("Done casting drug chem")

  #Is this for a classification task
  act.table <- act.table[,c("cell_name", "Compound", "NORM.AUC"),with=F]
  if (classify==T){
    print ("scaling")
    act.table$CLASS <- ifelse(act.table$NORM.AUC <= -0.8, 1,
                              ifelse(act.table$NORM.AUC >= 0.8, 0, 2))
    act.table <- act.table[CLASS!=2,]
    act.table$NORM.AUC <- act.table$CLASS
    act.table <- act.table[,c("cell_name", "Compound", "NORM.AUC"),with=F]
  }

  #Combine all
  setnames(act.table, c("cell_name", "DRUG", "NORM.AUC"))

  print ("Building main table")
  main.table <- merge(act.table, drugs.prof, by="DRUG")
  main.table <- merge(main.table, cosmic.cast, by = "cell_name")

  #Do we need to scale features
  if (feat.scaling==T){
    main.table <- Function.scale.data.table(main.table, col.protect = 1:3)
  }

  #Return
  return(main.table)
}
#####################################################################################################################
#Load files
FOLDER = "/tigress/zamalloa/CGP_FILES/" #For cluster
FOLDER = "/home/zamalloa/Documents/FOLDER/CGP_FILES/"

cgp.file = paste0(FOLDER, "gdsc_manova_input_w5.csv")
DRUGS.MET.PROFILE = readRDS(paste0(FOLDER, "DRUGS.MET.PROFILE"))
cosmic.exp <- fread(paste0(FOLDER, "CosmicCLP_CompleteGeneExpression.tsv"), sep="\t", header = T, drop=c(1,4))

#####################################################################################################################
#Execute
cgp_all_feat <- Function_cgp_all_feat(cgp.file, DRUGS.MET.PROFILE,
                                      cosmic.exp, norm.resp="Compound",
                                      response="AUC", feat.scaling=T, classify = F)

#####################################################################################################################
#Store
write.table(cgp_all_feat, FOLDER + "cgp_all_feat", sep="\t", col.names=T, row.names=F, quote=F)

print("Done writing")
