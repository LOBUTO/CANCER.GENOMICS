#Function.table.1.bmr.R
#081714
#Takes table.1 and uniprot_mapping.tab file to obtain p-values for each gene in each patient for background mutation rates
#NOTE: Results will contain records of silent genes!!  

Function.pre.process<-function(table.1, uniprot.file){
    #Pre process table.1 to give appropriate protein length and obtain total background amino acid length
    
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
    main.table<-as.data.table(merge(as.data.frame(table.1$table.1), 
        as.data.frame(uniprot.table), by="Hugo_Symbol"))

    #Get total hypothetical background mutation length
    background.length<-uniprot.table[Hugo_Symbol %in% as.vector(table.1$background.genes),]
    background.length<-sum(as.vector(background.length$Length))

    #Return
    return(list(main.table=main.table, background.length=background.length))
}

Function.BMR<-function(main.table, background.length) {
    #Calculates mutation significance of each gene based on the background mutation rate per patient
    #This is done for missense mutations only

    require(data.table)

    #Get total mutations per patients [PATIENT, BMR]
    BM.table<-main.table[,list(BM=sum(Missense)), by="PATIENT"]

    #Merge to main table
    main.table<-as.data.table(merge(as.data.frame(main.table),
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

args<-commandArgs(trailingOnly=T)

table.1<-readRDS(args[1]) 
uniprot.file<-args[2] #Such as 091714_UNIPROT_MAPPING.tab
output.file<-args[3]

pre.process<-Function.pre.process(table.1, uniprot.file)
print ("done pre-processing")

main.result<-Function.BMR(pre.process$main.table, pre.process$background.length)
print ("done with BMR p-values")

saveRDS(object=main.result, file=output.file)
print ("done saving ")