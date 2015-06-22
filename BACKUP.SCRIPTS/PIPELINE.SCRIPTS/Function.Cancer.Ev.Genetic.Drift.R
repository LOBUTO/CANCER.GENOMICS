#Function.Cancer.Ev.Genetic.Drift.R - IN PROGRESS
#091814
#Finds an optimal gene path across patients going from non-cancer to cancer based on significantly mutated genes and copy number variations
#Results in table giving type of driver gene found (Mutation of CNV)
#Requires:
#   exp.matrix.obj
#   table.1.pval (from Function.table.1.bmr.R())

Function.dist.matrix<-function(exp.matrix){
    #Pre-process exp matrix to produced median distance matrix to normal
    require(data.table)

    #First get pairwise distance matrix
    CORR.MATRIX<-cor(exp.matrix$combined.matrices[,c(exp.matrix$normal.patients,
    exp.matrix$cancer.patients)], method="spearman")
    DISS.MATRIX<-1-CORR.MATRIX

    #Get median distance of cancer patients to normal patients
    CANCER.DIST<-sapply(exp.matrix$cancer.patients, function(x) median(DISS.MATRIX[x, exp.matrix$normal.patients]))
    CANCER.DIST<-as.data.table(CANCER.DIST, keep.rownames=T)
    setnames(CANCER.DIST, c("PATIENT", "DIST.TO.NORMAL"))

    #Return
    return(CANCER.DIST)
}

Function.pre.main<-function(table.1.pval, cnv.table, threshold=0.05) {
    #Creates main table 
    require(data.table)

    #Filter cnv.table
    cnv.table<-cnv.table[,c("PATIENT", "Hugo_Symbol", "P.VAL.ADJ"),with=F] #TRACK
    cnv.table$TYPE<-"CNV" #TRACK

    #Filter for significant genes in table.1.pval based on threshold - MAYBE
    #main.table<-table.1.pval[P.VAL.ADJ<threshold, ]
    main.table<-table.1.pval[,c("PATIENT", "Hugo_Symbol", "P.VAL.ADJ"),with=F]
    main.table$TYPE<-"MUTATION" #TRACK

    #Merge mutation and cnv information per patient
    main.table<-rbind(main.table, cnv.table)

    #Obtain diss information per patient in main table
    main.table<-as.data.table(merge(as.data.frame(main.table), as.data.frame(CANCER.DIST), by="PATIENT"))

    #Order by shortest distance to normal - Beginning with closest to normal
    #[PATIENT, Hugo_Symbol, P.VAL.ADJ, TYPE, DIST.TO.NORMAL]
    main.table<-main.table[order(DIST.TO.NORMAL),]

    #Return
    return(main.table)
}

Function.CEGD<-function(main.table, seed=0, vector=c(), threshold=0.05) {
    #Uncovers genetic path from normal to cancer based on the following rules:
    #   First check if any gene in patient is already present in pool of chosen genes, if present, then pick it first. 
    #   If gene is not present in pool, then choose one with the lowest p-value.
    #   Dealing with ties at gene level: If more than one is present, choose one with lowest P.VAL.ADJ. If both have the same P.VAL.ADJ, choose one with the highest frequency in grand frequency table
    #   Dealing with ties at the patient level: Do all patients simultaneously per dist.to.normal rank

    require(data.table)

    #Create grand frequency table of genes across patients for those that pass threshold
    #main.table.freq<-main.table[P.VAL.ADJ<threshold, ]
    #main.table.freq<-as.data.table(table(as.vector(main.table.freq$Hugo_Symbol)), keep.rownames=T)
    #setnames(main.table.freq, c("Hugo_Symbol", "Freq"))
    #main.table.freq<-main.table.freq[order(Freq, decreasing=T),]
    
    #Create 2 grand frequency tables, one for CNV and one for mut
    #Both will be for values that pass threshold
    mut.table.freq<-main.table[P.VAL.ADJ<0.05 & TYPE=="MUTATION",]
    mut.table.freq<-as.data.table(table(as.vector(mut.table.freq$Hugo_Symbol)), keep.rownames=T)
    setnames(mut.table.freq, c("Hugo_Symbol", "Freq"))
    mut.table.freq<-mut.table.freq[order(Freq, decreasing=T),]

    cnv.table.freq<-main.table[P.VAL.ADJ<0.05 & TYPE=="CNV",]
    cnv.table.freq<-as.data.table(table(as.vector(cnv.table.freq$Hugo_Symbol)), keep.rownames=T)
    setnames(cnv.table.freq, c("Hugo_Symbol", "Freq"))
    cnv.table.freq<-cnv.table.freq[order(Freq, decreasing=T),]

    #Spit ordered main table to recursively applied method
    main.table.split<-split(main.table, main.table$DIST.TO.NORMAL)

    #Initialize gene pools - Two vectors corresponding to MUTATION and CNV pool
    if (length(vector)>0){
        ordered.candidate.pool<-vector
    } else if (seed==0){
        #ordered.candidate.pool<-c() #TRACK!
        mut.candidate.pool<-c()
        cnv.candidate.pool<-c()
    } else {
        ordered.candidate.pool<-as.vector(main.table.freq$Hugo_Symbol)[1:seed]
    }

    #Initialize count
    current.count<-0
    max.count<-length(names(main.table.split))

    #ALTERNATIVE
    main.table.1<-main.table[P.VAL.ADJ<0.05,]
    main.table.1.split<-split(main.table.1, main.table.1$DIST.TO.NORMAL)

    ####BUILD MODEL####
    main.table[TYPE=="MUTATION",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")][order(N.GENES, decreasing=T),]
    MUT.LM.EXP<-lm(log(N.GENES)~DIST.TO.NORMAL, main.table[TYPE=="MUTATION",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")] )
    MUT.LM.EXP<-MUT.LM.EXP$coefficients
    ggplot(main.table[TYPE=="MUTATION",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")], aes(DIST.TO.NORMAL, N.GENES)) + geom_point() + scale_y_log10()

    main.table[TYPE=="CNV",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")][order(N.GENES, decreasing=T),]
    CNV.LM<-lm(N.GENES~DIST.TO.NORMAL, main.table[TYPE=="CNV",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")] )
    CNV.LM<-CNV.LM$coefficients
    ggplot(main.table[TYPE=="CNV",][P.VAL.ADJ<0.05,][,list(N.GENES=length(Hugo_Symbol)),by=c("PATIENT","DIST.TO.NORMAL")], aes(DIST.TO.NORMAL, N.GENES)) + geom_point()

    ###################

    #Get coverage
    main.table.1.coverage<-unique(main.table.1[,c("PATIENT", "Hugo_Symbol"), with=F])

    #Apply method - Cannot do parallelization since it is recursive
    CEGD.Table<-lapply(main.table.1.split, function(x) {

        #Split per patient if multiple patients present
        target.table.split<-split(x, x$PATIENT)

        #Temporal candidate pool - Initialize with empty vectors
        temporal.mut.pool<-c()
        temporal.cnv.pool<-c()

        #Model Threshold
        MUT.TH<-MUT.LM.EXP[1]+unique(as.vector(x$DIST.TO.NORMAL)) * exp(MUT.LM.EXP[2])
        CNV.TH<-CNV.LM[1]+unique(as.vector(x$DIST.TO.NORMAL))*CNV.LM[2]

        #Analyse each patient table independently and simultanesouly
        to.pool<-sapply(target.table.split, function(y) {

            #target.table<-y[order(P.VAL.ADJ),] #TRACK!
            #candidate.genes<-as.vector(target.table$Hugo_Symbol) #TRACK!

            #Filter candidate genes by coverage overlap - TRACK!!!
            #Add genes only if at least 20% of their patients are not already on list
            CURRENT.POOL<-c(mut.candidate.pool, cnv.candidate.pool)
            CURRENT.POOL<-unique(as.vector(main.table.1.coverage[Hugo_Symbol %in% CURRENT.POOL,]$PATIENT))
            OVERLAP.TH<-main.table.1.coverage[Hugo_Symbol %in% y$Hugo_Symbol,]
            OVERLAP.TH<-OVERLAP.TH[,list(OVERLAP.DIFF=length(setdiff(PATIENT, CURRENT.POOL))/length(PATIENT)), by="Hugo_Symbol"]
            OVERLAP.TH<-OVERLAP.TH[OVERLAP.DIFF>0.2,]
            OVERLAP.CANDIDATES<-as.vector(OVERLAP.TH$Hugo_Symbol)
            y<-y[Hugo_Symbol %in% OVERLAP.CANDIDATES,]

            #Obtain candidate genes for cnv and mut by thresholding to model
            target.mut.table<-y[TYPE=="MUTATION",][order(P.VAL.ADJ),]
            if (nrow(target.mut.table)>MUT.TH){
                target.mut.table<-target.mut.table[1:MUT.TH,]
            }
            candidate.mut.genes<-as.vector(target.mut.table$Hugo_Symbol)
            
            target.cnv.table<-y[TYPE=="CNV",][order(P.VAL.ADJ),]
            if (nrow(target.cnv.table)>CNV.TH){
                target.cnv.table<-target.cnv.table[1:CNV.TH,]
            }
            candidate.cnv.genes<-as.vector(target.cnv.table$Hugo_Symbol)

            #If any candidate genes are already present, don't need to add more to pool
            if ( (sum(candidate.mut.genes %in% mut.candidate.pool)>0) | 
                (sum(candidate.cnv.genes %in% cnv.candidate.pool)>0) ) {
                a<-0

            #Else if candidates not in pool, add based on lowest p-val
            } else {

                #Add to corresponding pools if patient has both type (MUT and CNV) of aberrations, just MUT, or just CNV

                #If patient actually has both MUT and CNV aberrations - TRACK!!!
                if ( (length(candidate.mut.genes)>0) & (length(candidate.cnv.genes)>0) ){

                    #Chosen genes for MUTATION
                    top.p.val.mut<-min(as.vector(target.mut.table$P.VAL.ADJ))
                    chosen.mut.genes<-as.vector(target.mut.table[P.VAL.ADJ==top.p.val.mut,]$Hugo_Symbol)
                    chosen.mut.genes<-chosen.mut.genes[!(chosen.mut.genes %in% mut.candidate.pool)]

                    #If more than one MUTATION chose gene based on freq
                    if (length(chosen.mut.genes)>1){
                        chosen.mut.gene<-mut.table.freq[Hugo_Symbol %in% chosen.mut.genes,]$Hugo_Symbol[1]
                        temporal.mut.pool<<-c(temporal.mut.pool, chosen.mut.gene)
                    } else{
                        temporal.mut.pool<<-c(temporal.mut.pool, chosen.mut.genes)
                    }
                    
                    #Chosen genes for CNV
                    top.p.val.cnv<-min(as.vector(target.cnv.table$P.VAL.ADJ))
                    chosen.cnv.genes<-as.vector(target.cnv.table[P.VAL.ADJ==top.p.val.cnv,]$Hugo_Symbol)
                    chosen.cnv.genes<-chosen.cnv.genes[!(chosen.cnv.genes %in% cnv.candidate.pool)]

                    #If more than one CNV chose gene based on freq
                    if (length(chosen.cnv.genes)>1){
                        chosen.cnv.gene<-cnv.table.freq[Hugo_Symbol %in% chosen.cnv.genes,]$Hugo_Symbol[1]
                        temporal.cnv.pool<<-c(temporal.cnv.pool, chosen.cnv.gene)
                    } else{
                        temporal.cnv.pool<<-c(temporal.cnv.pool, chosen.cnv.genes)
                    }

                #OR If patient only has mutation information
                } else if ((length(candidate.mut.genes)>0) & (length(candidate.cnv.genes)==0)) {

                    #Chosen genes for MUTATION
                    top.p.val.mut<-min(as.vector(target.mut.table$P.VAL.ADJ))
                    chosen.mut.genes<-as.vector(target.mut.table[P.VAL.ADJ==top.p.val.mut,]$Hugo_Symbol)
                    chosen.mut.genes<-chosen.mut.genes[!(chosen.mut.genes %in% mut.candidate.pool)]   

                    #If more than one MUTATION chose gene based on freq
                    if (length(chosen.mut.genes)>1){
                        chosen.mut.gene<-mut.table.freq[Hugo_Symbol %in% chosen.mut.genes,]$Hugo_Symbol[1]
                        temporal.mut.pool<<-c(temporal.mut.pool, chosen.mut.gene)
                    } else{
                        temporal.mut.pool<<-c(temporal.mut.pool, chosen.mut.genes)
                    }

                #OR if patient only has cnv information
                } else {
                    #Chosen genes for CNV
                    top.p.val.cnv<-min(as.vector(target.cnv.table$P.VAL.ADJ))
                    chosen.cnv.genes<-as.vector(target.cnv.table[P.VAL.ADJ==top.p.val.cnv,]$Hugo_Symbol)
                    chosen.cnv.genes<-chosen.cnv.genes[!(chosen.cnv.genes %in% cnv.candidate.pool)]

                    #If more than one CNV chose gene based on freq
                    if (length(chosen.cnv.genes)>1){
                        chosen.cnv.gene<-cnv.table.freq[Hugo_Symbol %in% chosen.cnv.genes,]$Hugo_Symbol[1]
                        temporal.cnv.pool<<-c(temporal.cnv.pool, chosen.cnv.gene)
                    } else{
                        temporal.cnv.pool<<-c(temporal.cnv.pool, chosen.cnv.genes)
                    }
                }

                #top.p.val<-min(as.vector(target.table$P.VAL.ADJ)) #TRACK!
                #chosen.genes<-as.vector(target.table[P.VAL.ADJ==top.p.val,]$Hugo_Symbol) #TRACK!!

                #If more than one top chosen genes, add based on total frequency found in the grand frequency table
                #if (length(chosen.genes)>1){
                #    chosen.gene<-main.table.freq[Hugo_Symbol %in% chosen.genes,]$Hugo_Symbol[1]
                #    return (chosen.gene)
                #} else{
                #    return (chosen.genes)
                #}
            }
            return(c(temporal.mut.pool, temporal.cnv.pool))
        })
    
        print (to.pool)
        #Add chosen candidates to pool
        mut.candidate.pool<<-c(mut.candidate.pool, temporal.mut.pool)
        cnv.candidate.pool<<-c(cnv.candidate.pool, temporal.cnv.pool)
        #ordered.candidate.pool<<-c(ordered.candidate.pool, to.pool) #TRACK!!
        
        #Keep count
        current.count<<-current.count+1
        print (current.count/max.count)

    })
    print ("CEGD gene vector built")

    #Assign frequency to candidates
    CEGD.Table.MUT<-data.frame(Hugo_Symbol=mut.candidate.pool)
    CEGD.Table.MUT$TYPE<-"MUT"
    CEGD.Table.MUT<-as.data.table(merge(as.data.frame(CEGD.Table.MUT), as.data.frame(mut.table.freq)))
    CEGD.Table.CNV<-data.frame(Hugo_Symbol=cnv.candidate.pool)
    CEGD.Table.CNV$TYPE<-"CNV"
    CEGD.Table.CNV<-as.data.table(merge(as.data.frame(CEGD.Table.CNV), as.data.frame(cnv.table.freq)))
    #CEGD.Table<-data.frame(Hugo_Symbol=ordered.candidate.pool)
    CEGD.Table<-as.data.table(rbind(CEGD.Table.MUT, CEGD.Table.CNV))
    print (CEGD.Table)
    CEGD.Table<-as.data.table(merge(CEGD.Table, as.data.frame(main.table.freq), by="Hugo_Symbol"))   

    #Return
    return(CEGD.Table)
}

args<-commandArgs(trailingOnly=T)

table.1.pval<-readRDS(args[1]) 
exp.matrix<-readRDS(args[2])
cnv.filtered.table<-readRDS(args[3])
output.file<-args[4]

if (length(args)>4) {
  if (grepl(pattern="rds", as.character(args[5]))==TRUE)  {
    personal.vector<-readRDS(args[5])
    print (personal.vector)
    seed<-0
  } else {
    seed<-as.numeric(args[5])
    personal.vector<-c()
  }
} else {
    personal.vector<-c()
    seed<-0
}

CANCER.DIST<-Function.dist.matrix(exp.matrix)
print ("distance matrix built")

main.table<-Function.pre.main(table.1.pval, cnv.filtered.table)
print ("main table built")

CEGD<-Function.CEGD(main.table, seed=seed, vector=personal.vector)
print ("CEGD built")

saveRDS(object=CEGD, file=output.file)
print ("done saving ")