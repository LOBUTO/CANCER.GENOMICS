######Function.METABOLIC.COR.PVALUE.R#######
#040615
#Calculate enzymatic P-values from Function.METABOLIC.COR.R and Function.METABOLIC.COR.BACKGROUND.R output
#NOTE: Both ouputs must be computed from identical enzyme lists

##################FUNCTIONS################
library(data.table)
library(reshape2)
library(parallel)

Function.MAIN<-function(main.cor, back.cor){
  
  #Load files
  main.cor<-readRDS(main.cor)
  back.cor<-readRDS(back.cor)
  print ("done loading files")
  
  #Prepping parallelization
  print ("prepping for parallelization")
  nodes<-detectCores()
  cl<-makeCluster(nodes)
  setDefaultCluster(cl)
  clusterExport(cl, varlist=c("as.data.table","data.table", "main.cor", "back.cor") ,envir=environment())
  print ("Done exporting values")
  
  #Execute
  print ("Calculating p-values")
  main.list<-parLapply(cl, names(main.cor), function(x) {
    
    #Check that mut.class table is not NA due to empty correlation
    cor.table<-main.cor[[x]]$MAIN.TABLE
    cor.table<-cor.table[complete.cases(cor.table),]
    
    if (length(cor.table)>0){
      
      #Convert back to matrix
      cor.table<-as.matrix(cor.table, rownames.force=T)
      
      #Get genes of interest from cor.table
      genes<-rownames(cor.table)
      
      #Obtain background correlation for sample siz
      n.samples<-main.cor[[x]]$N.SAMPLES
      back.table<-back.cor[[as.character(n.samples)]]
      
      #Calculate empirical p-value
      p.vals<-sapply(genes, function(y) mean(as.vector(cor.table[y,])<=back.table[y,]))
      
      #Build results
      main.table<-data.table(Hugo_Symbol=genes, P.VAL=p.vals)
      main.table$CLASS<-x
      
      #Correct for multiple hypothesis testing
      main.table$P.VAL.ADJ<-p.adjust(main.table$P.VAL, method="fdr")
      
      #Return
      return(main.table)  
    }
  })
  
  #Stop parallelization
  stopCluster(cl)
  print ("Done parallelizing")
  
  #Converge into main table
  main.table<-do.call(rbind, main.list)
  
  #Clean up and return
  main.table<-main.table[order(P.VAL.ADJ),]
  return(main.table)
}

##########################################

################LOAD FILES################
args<-commandArgs(trailingOnly=T)
main.cor<-args[1]
back.cor<-args[2]
output.file<-args[3]
print("opened files")
##########################################

##################EXECUTE#################
main.function<-Function.MAIN(main.cor, back.cor)
print ("done with main function")

###############WRITING OUTPUT############
write.table(main.function, output.file, sep="\t", quote=F, row.names=F, col.names=T)
print ("Done writing to file")