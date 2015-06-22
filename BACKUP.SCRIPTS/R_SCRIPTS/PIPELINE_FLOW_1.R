#PIPELINE_FLOW1

library(diagram)

#Build matrix template
PIPELINE.1<-matrix(ncol=8,nrow=8,byrow=T, data=0)

#Build labels
LABELS<-c(".maf file","Processed .maf file \n (i.e. BRCA.table.1)",
          "...from PC", "Processed HMDB table \n (i.e. BRCA.table.2)", "Filtered table.2 against KEGG \n (i.e. BRCA.table.2.processed)",
          "...from PC", "Processed KEGG table \n (i.e BRCA.table.3)", "Filtered table.3 out of product \n (i.2. BRCA.table.3.processed)"
          )

#Build connectors (Column is origin and row is destination)
PIPELINE.1[2,1]<-"Function.process.maf.files"
PIPELINE.1[4,3]<-"Function..."
PIPELINE.1[7,6]<-"Function..."
PIPELINE.1[5,4]<-"Function.post.process.table.2"
PIPELINE.1[8,7]<-"Function.post.process.table.3"
PIPELINE.1[5,8]<-""

plotmat(PIPELINE.1, name=LABELS, pos=c(2,3,3),
        cex.txt=0.8, arr.pos=0.4)

openplotmat()

#Get hirerchical labels
LABELS<-c(".maf file","...from PC", "...from PC",
          "Processed .maf file \n (i.e. BRCA.table.1)", "Processed HMDB table \n (i.e. BRCA.table.2)", "Processed KEGG table \n (i.e BRCA.table.3)",
          "Filtered table.2 against KEGG \n (i.e. BRCA.table.2.processed)","Filtered table.3 out of product \n (i.e. BRCA.table.3.processed)"
          )

#Get hierchical coordinates - This is what builds the position matrix of labels
LABELS.POS<-coordinates(c(3, 3,2))

#Connect hirerchical labels with arrows
treearrow(from=LABELS.POS[1,],to=LABELS.POS[4,],lwd=6)  
treearrow(from=LABELS.POS[2,],to=LABELS.POS[5,],lwd=6)  
treearrow(from=LABELS.POS[3,],to=LABELS.POS[6,],lwd=6)  
treearrow(from=LABELS.POS[5,],to=LABELS.POS[7,],lwd=6)  
treearrow(from=LABELS.POS[6,],to=LABELS.POS[7,],lwd=6)  
treearrow(from=LABELS.POS[6,],to=LABELS.POS[8,],lwd=6)  

#Construct labels to diagrama
for (label in 1:length(LABELS)) {
  textround(LABELS.POS[label,],radx=0.08,rady=0.05,lab=LABELS[label])
}
