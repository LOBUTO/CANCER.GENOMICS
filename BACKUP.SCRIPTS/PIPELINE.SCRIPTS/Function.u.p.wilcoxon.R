#Function.u.p.wilcoxon.R
#100714
#Calculates the wilcoxon statistic per gene based on background mutation rate across patients
#Produces wilcoxon.table and mapped back genes to patients of wilcoxon results

#table.1$table.1
#table.1.obj<-copy(table.1)
#uniprot.file<-"/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/DATABASES/UNIPROT/091714_UNIPROT_MAPPING.tab"

Function.BMR<-function(bmr.main.table, background.length) {
    #Calculates mutation significance of each gene based on the background mutation rate per patient
    #This is done for missense mutations only

    require(data.table)

    #Get total mutations per patients [PATIENT, BM]
    BM.table<-bmr.main.table[,list(BM=sum(Missense)), by="PATIENT"]

    #Merge to main table
    main.table<-as.data.table(merge(as.data.frame(bmr.main.table),
        as.data.frame(BM.table), by="PATIENT"))

    #Calculate significance
    sig.table<-main.table[,list(P.VAL=phyper(q=Missense-1, m=Length, 
        n=background.length-Length, k=BM, lower.tail=F)), 
    by=c("PATIENT", "Hugo_Symbol", "Missense", "Silent", "Length")]
    sig.table$P.VAL.ADJ<-p.adjust(sig.table$P.VAL, method="fdr")

    #Clean up and Return
    sig.table<-sig.table[order(P.VAL.ADJ),]
    return(sig.table)
}

Function.pre.process<-function(table.1, background.genes, uniprot.file, bmr=TRUE){
    #Pre process table.1 to give appropriate protein length and obtain total background amino acid length
    #If bmr=TRUE, will apply background mutation function before wilcoxon analysis
    
    require(data.table)

    #Use uniprot gene names and synonyms to map lengths for table.1 [Hugo_Symbol, Length]
    uniprot.table<-as.data.table(read.csv(uniprot.file, header=T, sep="\t", stringsAsFactors=F))
    uniprot.table<-uniprot.table[,c("Gene.names", "Length"),with=F]
    uniprot.table$ID<-1:nrow(uniprot.table)
    uniprot.table<-uniprot.table[,list(Hugo_Symbol= sub(";", "",strsplit(Gene.names," ")[[1]]), Length=Length ), by="ID"]    
    uniprot.table$ID<-NULL
    
    #Filter for duplicated gene records based on largest gene size
    uniprot.table<-uniprot.table[order(Length, decreasing=T),]
    uniprot.table<-uniprot.table[!duplicated(Hugo_Symbol),]

    #Merge with table.1 to get length information 
    #[Hugo_Symbol, Missense, Silent, PATIENT, Length]
    #main.table<-table.1.pval[P.VAL.ADJ<0.05,][,c("Hugo_Symbol", "Missense","Silent", "PATIENT"), with=F]
    #main.table<-as.data.table(merge(as.data.frame(main.table), as.data.frame(uniprot.table), by="Hugo_Symbol"))
    main.table<-as.data.table(merge(as.data.frame(table.1), 
        as.data.frame(uniprot.table), by="Hugo_Symbol"))

    #Get total hypothetical background mutation length
    background.length<-uniprot.table[Hugo_Symbol %in% as.vector(background.genes),]
    background.length<-sum(as.vector(background.length$Length))
    background.length<-9000000
    #background.length<-10000000
    #background.length<-11317118

    #Are we pre-filtering by bmr.function()?
    if (bmr==TRUE){
        main.table<-Function.BMR(main.table, background.length)
        main.table<-main.table[,c("Hugo_Symbol", "Missense","Silent", "PATIENT","Length","P.VAL.ADJ"), with=F]
    }

    #Return
    return(list(main.table=main.table, background.length=background.length))
}

Function.WILCOXON.BMR<-function(main.table, background.length, bmr=TRUE, lm=FALSE) {
    #Calculates mutation significance of each gene based on the background mutation rate across patients using a wilcoxon measure
    #This is done for missense mutations only

    require(data.table)

    #Get total mutations per patients [PATIENT, BMR] 
    #MODIFY BASED ON ALL MISSENSE!!!! THAT MEANS NOT FILTERING ABOVE main.table BY P.VAL.ADJ!!!

    #Do we want to apply linear modeling to fit expected silent mutations?
    if (lm==TRUE){

        #Get actual count of mutations per patient
        patient.count<-main.table[,list(Patient.Missense=sum(Missense), 
            Patient.Silent=sum(Silent), Patient.TOTAL=sum(Missense) + sum(Silent)), by="PATIENT"]

        #Fit to linear model Total~Silent
        model.lm<-lm(Patient.Silent~Patient.TOTAL, data=patient.count)
        model.coeff<-coef(model.lm)
        patient.count$BM<-patient.count$Patient.TOTAL*model.coeff[[2]] + model.coeff[[1]]
        patient.count$BMR<-patient.count$BM/background.length
        BMR.table<-patient.count[,c("PATIENT", "BMR"),with=F]
        
    } else {
        BMR.table<-main.table[,list(BMR=sum(Silent)/background.length), by="PATIENT"]
    }

    #Then filter to get only significant genes from bmr test if applied
    if (bmr==TRUE){
        main.table<-main.table[P.VAL.ADJ<0.05,]
    }   

    #Merge to main table
    main.table<-as.data.table(merge(as.data.frame(main.table),
        as.data.frame(BMR.table), by="PATIENT"))

    #Filter out non-missense records
    main.table<-main.table[Missense!=0,]

    #Get gene mutation rate - MAY NEED TO BE MODELED A PRIORI (PER PATIENT, PER GENE??)
    main.table$GMR<-main.table$Missense/main.table$Length

    #Calculate wilcoxon per gene
    wilcoxon.table<-main.table[,list(P.VAL=wilcox.test(GMR, BMR, alternative="greater", paired=TRUE)$p.value, 
        W=wilcox.test(GMR, BMR, alternative="greater", paired=TRUE)$statistic), by="Hugo_Symbol"]
    wilcoxon.table$P.VAL.ADJ<-p.adjust(wilcoxon.table$P.VAL, method="fdr")

    #Clean up and Return
    wilcoxon.table<-wilcoxon.table[order(P.VAL.ADJ),]
    return(wilcoxon.table)
}

Function.WILCOXON.MAP<-function(wilcoxon.table, main.table, threshold){
    #Maps results of wilcoxon table to patients
    require(data.table)

    #Map back to significant wilcoxon genes
    target.map<-merge(main.table[,c("Hugo_Symbol","Missense","Silent","PATIENT"), with=F][Missense!=0,],
        wilcoxon.table[,c("Hugo_Symbol","W","P.VAL.ADJ"),with=F][P.VAL.ADJ<threshold,], by="Hugo_Symbol")

    #Return
    return(target.map)
}

Function.SUBTYPE.SPLIT<-function(table.1.obj, subtype.table){
    #Splits table.1.object by molecular subtype found in subtype table
    #Make sure subtype table is from same cancer type

    require(data.table)
    
    #Obtain subtypes per patient
    target.table<-merge(table.1.obj$table.1, subtype.table[,c("PATIENT", "TYPE"), with=F],
        by="PATIENT")

    #Split by subtypes
    target.tables<-split(target.table, target.table$TYPE)

    #Return lists of table.1s
    return(target.tables)
}

Function.WILCOXON.MAX<-function(wilcoxon.table, table.1){
    #Gets maximum coverage based on wilcoxon genes by W and minimum number of genes needed to achieve such coverage by W

    require(data.table)

    #Get restricted wilcoxon.table 
    wilcoxon.table<-wilcoxon.table[W>1,]
    wilcoxon.table<-wilcoxon.table[order(W, decreasing=T),]
    wilcoxon.genes<-unique(as.vector(wilcoxon.table$Hugo_Symbol))

    #Get maximum number of patients covered by restricted wilcoxon genes
    max.patients<-length(unique(as.vector(table.1[Hugo_Symbol %in% wilcoxon.genes,][Missense!=0,]$PATIENT)))

    #Get number of patients covered by ordered list of wilcoxon genes
    patient.by.w<-sapply(1:length(wilcoxon.genes), 
        function(x) 
        length(unique(as.vector(table.1[Hugo_Symbol %in% wilcoxon.genes[1:x],][Missense!=0,]$PATIENT))))
    patient.by.w<-data.table(Hugo_Symbol=wilcoxon.genes, COVERAGE=patient.by.w)

    #Get maximum attainable number of patients by minimum gene coverage
    min.genes<-which.max(patient.by.w$COVERAGE)

    #Return
    return(list(max.patients=max.patients, min.genes=min.genes))
}

args<-commandArgs(trailingOnly=T)

table.1.obj<-readRDS(args[1]) 
uniprot.file<-args[2]
output.file<-args[3]

if (length(args)==3){
    pre.process<-Function.pre.process(table.1.obj$table.1, table.1.obj$background.genes, uniprot.file, bmr=TRUE)
    print ("Done with pre-processing")

    wilcoxon.table<-Function.WILCOXON.BMR(pre.process$main.table, pre.process$background.length, lm=F)
    print ("Wilcoxon table build")

    wilcoxon.map<-Function.WILCOXON.MAP(wilcoxon.table, pre.process$main.table, 0.05)
    print ("Done Wilcoxon mapping")

    coverage.table<-data.frame(TOTAL=length(unique(as.vector(pre.process$main.table$PATIENT))),
        COVERED=length(unique(as.vector(wilcoxon.map$PATIENT))))
    print ("Done Coverage")

    saveRDS(object=list(wilcoxon.table=wilcoxon.table, wilcoxon.map=wilcoxon.map, 
        coverage.table=coverage.table), file=output.file)
    print ("done saving")

#subtype.table<-readRDS("/Users/jzamalloa/Desktop/FOLDER/ECLIPSE/workspace/Rotation/PIPELINES/METABOLIC.DRIVERS/OBJECTS/BRCA/100414.BRCA.PATIENT.SUBTYPES.rds")
} else if (length(args)==4){
    subtype.table<-readRDS(args[4])

    #Split into table.1s per subtype
    split.table.1.obj<-Function.SUBTYPE.SPLIT(table.1.obj, subtype.table)

    #Initialize coverage table and result list
    bmr=TRUE
    if (bmr==TRUE){
        COVERAGE.TABLE<-data.frame(SUBTYPE=c(), ALL.PATIENTS=c(), PRE.WILCOXON=c(),COVERED=c(), 
        N.STARTING.GENES=c(),N.GENES.PRE=c(), N.GENES=c(), MAX.W.PATIENTS=c(), MIN.W.GENES=c())    
    } else {
        COVERAGE.TABLE<-data.frame(SUBTYPE=c(), ALL.PATIENTS=c(),COVERED=c(), 
        N.STARTING.GENES=c(), N.GENES=c(),MAX.W.PATIENTS=c(), MIN.W.GENES=c())
    }
    SUBTYPE.TABLE.1.OBJ<-list()

    #Run functions per each subtype
    for (subtype in names(split.table.1.obj)){

        pre.process<-Function.pre.process(split.table.1.obj[[subtype]], table.1.obj$background.genes, 
            uniprot.file, bmr=bmr)
        print ("Done with pre-processing")

        wilcoxon.table<-Function.WILCOXON.BMR(pre.process$main.table, pre.process$background.length, lm=T, bmr=bmr)
        print ("Wilcoxon table build")

        wilcoxon.map<-Function.WILCOXON.MAP(wilcoxon.table, pre.process$main.table, 0.05)
        print ("Done Wilcoxon mapping")

        maximum.coverage<-Function.WILCOXON.MAX(wilcoxon.table, split.table.1.obj[[subtype]])
        print ("Done maximum coverage")

        #Append to coverage table and main list
        require(data.table)
        
        SUBTYPE.TABLE.1.OBJ[[subtype]]<-list(wilcoxon.table=wilcoxon.table, wilcoxon.map=wilcoxon.map)
        
        if (bmr==TRUE){
            COVERED.TABLE<-data.frame(SUBTYPE=subtype, 
            ALL.PATIENTS=length(unique(as.vector(split.table.1.obj[[subtype]]$PATIENT))), 
            PRE.WILCOXON=length(unique(as.vector(pre.process$main.table[P.VAL.ADJ<0.05,]$PATIENT))),
            COVERED=length(unique(as.vector(wilcoxon.map$PATIENT))),
            N.STARTING.GENES=length(unique(as.vector(split.table.1.obj[[subtype]]$Hugo_Symbol))),
            N.GENES.PRE=length(unique(as.vector(pre.process$main.table[P.VAL.ADJ<0.05,]$Hugo_Symbol))),
            N.GENES=length(unique(as.vector(wilcoxon.map$Hugo_Symbol))),
            MAX.W.PATIENTS=maximum.coverage$max.patients,
            MIN.W.GENES=maximum.coverage$min.genes)    
        } else {
            COVERED.TABLE<-data.frame(SUBTYPE=subtype, 
            ALL.PATIENTS=length(unique(as.vector(split.table.1.obj[[subtype]]$PATIENT))), 
            COVERED=length(unique(as.vector(wilcoxon.map$PATIENT))),
            N.STARTING.GENES=length(unique(as.vector(split.table.1.obj[[subtype]]$Hugo_Symbol))),
            N.GENES=length(unique(as.vector(wilcoxon.map$Hugo_Symbol))),
            MAX.W.PATIENTS=maximum.coverage$max.patients,
            MIN.W.GENES=maximum.coverage$min.genes)
        }
        
        COVERAGE.TABLE<-rbind(COVERAGE.TABLE, COVERED.TABLE)
        print ("Done Coverage")
    }

    #Clean up and return
    saveRDS(object=list(SUBTYPES=SUBTYPE.TABLE.1.OBJ, COVERAGE.TABLE=COVERAGE.TABLE), file=output.file)
}

#length(unique(as.vector(table.1.obj$table.1$PATIENT)))
#length(unique(as.vector(table.1.obj$table.1[Hugo_Symbol %in% as.vector(wilcoxon.table[P.VAL.ADJ<0.01,]$Hugo_Symbol),]$PATIENT)))

# wilcoxon.table[P.VAL.ADJ<0.05,][order(P.VAL.ADJ),]
# wilcoxon.table[Hugo_Symbol %in% as.vector(COSMIC.BRCA$Symbol),]

# BRCA.CANCER.DIST.1.5
# test.1<-Function.pre.process(table.1.obj, uniprot.file, bmr=TRUE)
# test.2<-Function.WILCOXON.BMR(test.1$main.table, test.1$background.length, lm=TRUE)
# test.2[P.VAL.ADJ<0.05,]
# test.2[Hugo_Symbol %in% COSMIC.BRCA$Symbol,]
# length(unique(as.vector(table.1$table.1$PATIENT)))
# length(unique(as.vector(table.1$table.1[Hugo_Symbol %in% as.vector(test.2[P.VAL.ADJ<0.05,]$Hugo_Symbol),]$PATIENT)))

# GBM.table.1<-readRDS("GBM/100814.GBM.Table.1.rds")
# test.1.GBM<-Function.pre.process(GBM.table.1, uniprot.file, bmr=TRUE)
# test.2.GBM<-Function.WILCOXON.BMR(test.1.GBM$main.table, test.1.GBM$background.length, lm=TRUE)
# test.2.GBM[P.VAL.ADJ<0.05,]
# test.2.GBM[Hugo_Symbol %in% COSMIC.GBM$Symbol,]

# length(unique(as.vector(GBM.table.1$table.1$PATIENT)))
# length(unique(as.vector(GBM.table.1$table.1[Hugo_Symbol %in% as.vector(test.2.GBM[P.VAL.ADJ<0.05,]$Hugo_Symbol),]$PATIENT)))

# COAD.table.1<-readRDS("COAD/100914.COAD.Table.1.rds")
# COAD.GA.table.1<-readRDS("COAD/100914.COAD.GA.Table.1.rds")
# COAD.SOLID.table.1<-readRDS("COAD/100914.COAD.SOLID.Table.1.rds")

# test.1.COAD<-Function.pre.process(COAD.table.1, uniprot.file, bmr=TRUE)
# test.2.COAD<-Function.WILCOXON.BMR(test.1.COAD$main.table, test.1.COAD$background.length, lm=TRUE)
# test.2.COAD[P.VAL.ADJ<0.05,]
# test.2.COAD[Hugo_Symbol %in% COSMIC.COAD$Symbol,]

# length(unique(as.vector(COAD.table.1$table.1$PATIENT)))
# length(unique(as.vector(COAD.table.1$table.1[Hugo_Symbol %in% as.vector(test.2.COAD[P.VAL.ADJ<0.05,]$Hugo_Symbol),]$PATIENT)))

# OV.table.1<-readRDS("OV/101014.OV.Table.1.rds")
# test.1.OV<-Function.pre.process(OV.table.1, uniprot.file, bmr=TRUE)
# test.2.OV<-Function.WILCOXON.BMR(test.1.OV$main.table, test.1.OV$background.length, lm=TRUE)
# test.2.OV[P.VAL.ADJ<0.05,]
# test.2.OV[Hugo_Symbol %in% COSMIC.OV$Symbol,]

# length(unique(as.vector(OV.table.1$table.1$PATIENT)))
# length(unique(as.vector(OV.table.1$table.1[Hugo_Symbol %in% as.vector(test.2.OV[P.VAL.ADJ<0.05,]$Hugo_Symbol),]$PATIENT)))

# READ.table.1<-readRDS("READ/101014.READ.Table.1.rds")
# test.1.READ<-Function.pre.process(READ.table.1, uniprot.file, bmr=TRUE)
# test.2.READ<-Function.WILCOXON.BMR(test.1.READ$main.table, test.1.READ$background.length, lm=TRUE)
# test.2.READ[P.VAL.ADJ<0.05,]
# test.2.READ[Hugo_Symbol %in% COSMIC.READ$Symbol,]

# length(unique(as.vector(READ.table.1$table.1$PATIENT)))
# length(unique(as.vector(READ.table.1$table.1[Hugo_Symbol %in% as.vector(test.2.READ[P.VAL.ADJ<0.05,]$Hugo_Symbol),]$PATIENT)))
