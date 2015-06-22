#Function.Filter.CNV.R
#092214
#Function to filter CNV data obtained from GISTIC (Function.process.GISTIC.TH.R()) based on their influence on gene expression

Function.pre.process<-function(exp.obj, cnv.table, threshold){
    #Pre-process expression matrix and CNV table

    require(data.table)
    require(reshape2)

    #Threshold CNV Table before merging
    cnv.table<<-cnv.table[abs(CNV.TH)>=threshold,]

    #Convert whole expression matrix to short format
    EXP.MELTED<-as.data.table(melt(exp.obj$combined.matrices))
    setnames(EXP.MELTED, c("Hugo_Symbol", "PATIENT", "EXP"))

    #Obtain normal expression table and filter out genes not found in CNV.table
    NORMAL.EXP<-EXP.MELTED[grepl("11A", PATIENT),]
    NORMAL.EXP$CNV.TH<-0
    NORMAL.EXP$TYPE<-"NORMAL"
    NORMAL.EXP<-NORMAL.EXP[Hugo_Symbol %in% as.vector(cnv.table$Hugo_Symbol),]

    #Merge to main table - Since merging with CNV Table, we are removing first normal patients
    MERGED.TABLE<-as.data.table(merge(as.data.frame(EXP.MELTED), as.data.frame(cnv.table), by=c("PATIENT", "Hugo_Symbol")))
    MERGED.TABLE$TYPE<-"TUMOR"

    #Bind to normals
    MERGED.TABLE<-rbind(MERGED.TABLE, NORMAL.EXP)

    #Return 
    return(MERGED.TABLE)
}

Function.CNV.SIG<-function(MERGED.TABLE) {
    #Find significantly dysregulated CNV genes

    require(data.table)

    internal.function<-function(exp, cnv.th, type, hugo) {
        internal.table<-as.data.table(data.frame(EXP=exp, CNV.TH=cnv.th, TYPE=type))
        #print (unique(as.vector(MERGED.TABLE$Hugo_Symbol))[hugo])
        #print (internal.table)
        #[EXP, CNV.TH, TYPE]
        normals<-internal.table[TYPE=="NORMAL",]
        tumors<-internal.table[TYPE=="TUMOR",]

        #print (length(as.vector(tumors[CNV.TH<0,]$EXP)))
        #print (as.vector(tumors[CNV.TH<0,]$EXP))
        #print (length(as.vector(tumors[CNV.TH>0,]$EXP)))
        #print (as.vector(tumors[CNV.TH>0,]$EXP))

        if (length(as.vector(tumors[CNV.TH<0,]$EXP))>1){
            NEG.PVAL=wilcox.test(as.vector(tumors[CNV.TH<0,]$EXP), as.vector(normals$EXP), alternative="less")$p.value
        } else{
            NEG.PVAL=1
        }

        if (length(as.vector(tumors[CNV.TH>0,]$EXP))>1){
            POS.PVAL=wilcox.test(as.vector(tumors[CNV.TH>0,]$EXP), as.vector(normals$EXP), alternative="greater")$p.value
        } else{
            POS.PVAL=1
        }

        #Keep count
        count<<-count+1
        print (count/max.count)
                   
        return (list(NEG.PVAL=NEG.PVAL, POS.PVAL=POS.PVAL))
    }

    #To keep count
    count<-0
    max.count<-length(unique(as.vector(MERGED.TABLE$Hugo_Symbol)))

    #Calculate pvalues for both positive and negative CNV 
    SIG.TEST<-MERGED.TABLE[, internal.function(EXP, CNV.TH, TYPE), by="Hugo_Symbol"]

    #Correct for multiple hypothesis
    SIG.TEST$NEG.PVAL.ADJ<-p.adjust(SIG.TEST$NEG.PVAL, method="fdr")
    SIG.TEST$POS.PVAL.ADJ<-p.adjust(SIG.TEST$POS.PVAL, method="fdr")

    #Return
    return(SIG.TEST)

}

Function.Apply.Filter<-function(SIG.TEST){
    #Applies significant structural rearregenments found to CNV table per patient

    require(data.table)

    #Filters genes for fdr in pos and neg cases
    POS.SIG<-SIG.TEST[,c("Hugo_Symbol", "POS.PVAL.ADJ"), with=F]
    POS.SIG<-POS.SIG[POS.PVAL.ADJ<0.05,]

    NEG.SIG<-SIG.TEST[,c("Hugo_Symbol", "NEG.PVAL.ADJ"), with=F]
    NEG.SIG<-NEG.SIG[NEG.PVAL.ADJ<0.05,]    

    #Splits cnv.table into neg (deletions) and pos (insertions)
    POS.CNV<-cnv.table[sign(CNV.TH)==1,]
    NEG.CNV<-cnv.table[sign(CNV.TH)==-1,]

    #Apply filter
    POS.CNV<-as.data.table(merge(as.data.frame(POS.CNV), as.data.frame(POS.SIG), by="Hugo_Symbol"))
    NEG.CNV<-as.data.table(merge(as.data.frame(NEG.CNV), as.data.frame(NEG.SIG), by="Hugo_Symbol"))

    #Clean up and return
    setnames(POS.CNV, c("Hugo_Symbol", "PATIENT", "CNV.TH", "P.VAL.ADJ"))
    setnames(NEG.CNV, c("Hugo_Symbol", "PATIENT", "CNV.TH", "P.VAL.ADJ"))
    FILTERED.CNV<-unique(rbind(POS.CNV, NEG.CNV))
    return(FILTERED.CNV)
}

#ggplot(FILTERED.CNV[,list(CNV.COUNT=length(Hugo_Symbol)), by="PATIENT"], aes(CNV.COUNT)) + geom_histogram()

args<-commandArgs(trailingOnly=T)

exp.obj<-readRDS(args[1]) 
cnv.table<-readRDS(args[2])
output.file<-args[3]

MERGED.TABLE<-Function.pre.process(exp.obj, cnv.table, 2)
print ("done pre-processing")

SIG.TEST<-Function.CNV.SIG(MERGED.TABLE)
print ("done with p-values")

main.result<-Function.Apply.Filter(SIG.TEST)
print ("done applying filter")

saveRDS(object=main.result, file=output.file)
print ("done saving")