#Function.Parse.Biogrid.R
#092114
#Parses Biogrid files to obtain pairwise interaction of Hugo identifiers
#NOTE: Table may contain synonyms, so it will not work for enrichment calculations

BIOGRID.short<-unique(BIOGRID[,c("Alt.IDs.Interactor.A", "Alt.IDs.Interactor.B"), with=F])

BIOGRID.short<-unique(BIOGRID.short[,list(Hugo.A=strsplit(Alt.IDs.Interactor.A, "\\|")[[1]]), by=c("Alt.IDs.Interactor.A", "Alt.IDs.Interactor.B")])

BIOGRID.short$Alt.IDs.Interactor.A<-NULL

BIOGRID.short<-unique(BIOGRID.short[,list(Hugo.B=strsplit(Alt.IDs.Interactor.B, "\\|")[[1]]), by=c("Hugo.A", "Alt.IDs.Interactor.B")])

BIOGRID.short$Alt.IDs.Interactor.B<-NULL

BIOGRID.short<-unique(BIOGRID.short)

#Clean names
BIOGRID.short$Hugo.A<-sapply(BIOGRID.short$Hugo.A, function(x) 
    strsplit(x, ":")[[1]][2])
BIOGRID.short$Hugo.B<-sapply(BIOGRID.short$Hugo.B, function(x) 
    strsplit(x, ":")[[1]][2])
BIOGRID.short<-unique(BIOGRID.short)

#Get rid of "self interactors"
BIOGRID.short<-BIOGRID.short[Hugo.A!=Hugo.B,]

#Convert so that all unique interactions can be called from first column
BIOGRID.short<-unique(rbind(BIOGRID.short, setnames(rev(BIOGRID.short), names(BIOGRID.short))))

