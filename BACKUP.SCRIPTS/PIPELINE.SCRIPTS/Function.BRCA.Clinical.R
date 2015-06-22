#Function.BRCA.Clinical.R
#093014
#Process files of the form:
# nationwidechildrens.org_clinical_patient_brca.txt
# nationwidechildrens.org_biospecimen_cqcf_brca.txt

Function.BRCA.Clinical<-function(clinical_patient.file, clinical_cqcf.file){

    library(data.table)

    #Open first clinical file
    clinical<-as.data.table(read.csv(clinical_patient.file, sep="\t", header=T, skip=2))

    #Filter for lines of interest (CHANGE IF NECESSARY!!)
    clinical<-clinical[,c(1,7,9,14,15,16,34,35,38:40,41,44,50,56), with=F]
    setnames(clinical, c("PATIENT", "GENDER", "RACE", "VITAL.STATUS", "DAYS.TO.LAST.CONTACT", "DAYS.TO.DEATH", "EXAMINED.LYMPH.NODES", "POSITIVE.LYMPH.NODES",
                         "T.STATUS","N.STATUS","M.STATUS", "TUMOR.STAGE", "ER.STATUS", "PR.STATUS","HER2.STATUS"))
    
    #Get percent positive lymph nodes
    clinical$PERCENT.LYMPH.NODES<-ifelse( (clinical$EXAMINED.LYMPH.NODES!="[Not Available]" & clinical$POSITIVE.LYMPH.NODES!="[Not Available]"), as.numeric(as.character(clinical$POSITIVE.LYMPH.NODES))/as.numeric(as.character(clinical$EXAMINED.LYMPH.NODES)), NA)

    #Open file with procurement information
    clinical.procurement<-as.data.table(read.csv(clinical_cqcf.file, sep="\t", header=T, skip=1))
    clinical.procurement<-clinical.procurement[,c(1,5), with=F]
    setnames(clinical.procurement, c("PATIENT","DAYS.TO.PROCUREMENT"))

    #Clean procurement info
    clinical.procurement<-clinical.procurement[DAYS.TO.PROCUREMENT!="[Not Available]",]
    clinical.procurement$DAYS.TO.PROCUREMENT<-as.numeric(as.character(clinical.procurement$DAYS.TO.PROCUREMENT))

    #Merge all clinical
    main.clinical<-as.data.table(merge(as.data.frame(clinical), as.data.frame(clinical.procurement)))

    #Obtain overall survival - Longest alive time (DEATH or CONTACT)
    main.clinical$OVERALL.SURVIVAL<-ifelse(main.clinical$DAYS.TO.LAST.CONTACT!="[Not Available]",
        as.numeric(as.character(main.clinical$DAYS.TO.LAST.CONTACT)), 
        as.numeric(as.character(main.clinical$DAYS.TO.DEATH)))

    #Obtain TCGA surival - From procurement to overall
    main.clinical$TCGA.SURVIVAL<-main.clinical$OVERALL.SURVIVAL - main.clinical$DAYS.TO.PROCUREMENT

    #Clean up patient nomenclature
    main.clinical$PATIENT<-as.character(main.clinical$PATIENT)
    main.clinical$PATIENT<-sapply(main.clinical$PATIENT, function(x) paste0(strsplit(x, "-")[[1]],collapse="."))

    #Clean up sample types - Add 01A and 01B to patient labels as well as 11A and 11B normals
    main.clinical.B<-copy(main.clinical)
    main.clinical.C<-copy(main.clinical)
    main.clinical.D<-copy(main.clinical)
    main.clinical$PATIENT<-paste0(main.clinical$PATIENT, sep=".01A")
    main.clinical.B$PATIENT<-paste0(main.clinical.B$PATIENT, sep=".01B")
    main.clinical.C$PATIENT<-paste0(main.clinical.C$PATIENT, sep=".11A")
    main.clinical.D$PATIENT<-paste0(main.clinical.D$PATIENT, sep=".11B")
    main.clinical<-rbind(main.clinical, main.clinical.B, main.clinical.C, main.clinical.D)
    
    #Summarize Clinical sets (i.e I=c(I,IA, IB))
    main.clinical$TUMOR.STAGE.SET<-ifelse(main.clinical$TUMOR.STAGE %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC"), "Stage I", 
        ifelse(main.clinical$TUMOR.STAGE %in% c("Stage II", "Stage IIA", "Stage IIB"), "Stage II",
        ifelse(main.clinical$TUMOR.STAGE %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"), "Stage III", 
        ifelse(main.clinical$TUMOR.STAGE %in% c("Stage IV", "Stage IVA", "Stage IVB"), "Stage IV",
        ifelse(main.clinical$TUMOR.STAGE %in% c("Stage Tis"), "Stage 0", as.character(main.clinical$TUMOR.STAGE) )))))
  
    #Group tumor stages into early or late
    main.clinical$TUMOR.STAGING<-ifelse(main.clinical$TUMOR.STAGE.SET %in% c("Stage 0", "Stage I", "Stage II"), "Early",
                                 ifelse(main.clinical$TUMOR.STAGE.SET %in% c("Stage III", "Stage IV"), "Late",
                                 ifelse(main.clinical$T.STATUS        %in% c("T4a", "T4b", "T4c","T4d"), "Late",
                                 ifelse( (main.clinical$T.STATUS=="T2" & main.clinical$N.STATUS %in% c("N2","N2a","N2b","N2c","N2d") & main.clinical$M.STATUS=="M0"),"Late",
                                 ifelse( (main.clinical$T.STATUS %in% c("T1a", "T1b","T1c","T1d") & main.clinical$N.STATUS %in% c("N1","N1a","N1b","N1c","N1d") & 
                                                                                                        main.clinical$M.STATUS=="M0"), "Early",
                                 ifelse( (main.clinical$T.STATUS=="T2" & main.clinical$N.STATUS=="N0" & main.clinical$M.STATUS=="M0"), "Early",
                                        as.character(main.clinical$TUMOR.STAGE.SET) ))))))
    
    #Return
    return(main.clinical)
}

args<-commandArgs(trailingOnly=T)

clinical_patient.file<-args[1]
clinical_cqcf.file<-args[2] 
output.file<-args[3]

main.result<-Function.BRCA.Clinical(clinical_patient.file, clinical_cqcf.file)

saveRDS(object=main.result, file=output.file)
print ("done saving ")