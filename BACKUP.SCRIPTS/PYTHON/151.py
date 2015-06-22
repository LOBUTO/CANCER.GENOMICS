#Uses the GRCh37 chromosome gbs files to construct table of mRNA exons at gene level

from Bio import SeqIO
import subprocess

#Load gbs files
GBS_FILES=subprocess.check_output("ls", cwd="SOFTWARE/pipeline/data/gbs").splitlines()
GBS_FOLDER="SOFTWARE/pipeline/data/gbs"

#Open output file
#Write to file based on priority of feature
FILE_OUT1=open("DATABASES/REFSEQ/070114_GENE_COORDIANTES","w")
FILE_OUT1.write("Hugo_Symbol"+"\t"+"Chrom"+"\t"+"Strand"+"\t"+
                "START"+"\t"+"END"+"\t"+"FEATURE"+"\t"+"FEAT_START"+"\t"+"FEAT_END")

#Loop through files
for gbs in GBS_FILES:
    CHROMOSOME=gbs.split("chr")[1].split(".")[0]
    FILE_IN=GBS_FOLDER+"/"+gbs

    #Open file
    RECORD=SeqIO.read(open(r"%s"%FILE_IN), "genbank")
    
    #Filter for gene, exon, mRNA or CDS
    RECORD=filter(lambda x: x.type in ["gene", "exon","mRNA", "CDS"], RECORD.features)
    
    #Create gene dictionary template
    DICT=dict((y.qualifiers["gene"][0],{}) for y in RECORD if y.type=="gene")
    
    #Loop to obtain gene position and feature position per gene as dictionary
    for record in RECORD:
        print record
        
        #Get gene indices and strand first
        if record.type=="gene":
            DICT[record.qualifiers["gene"][0]]["START"]=record.location.nofuzzy_start
            DICT[record.qualifiers["gene"][0]]["END"]=record.location.nofuzzy_end
        
            DICT[record.qualifiers["gene"][0]]["STRAND"]=record.strand
            
        #Given that we have filtered for exon, mRNA or CDS, get features next 
        else: 
            TYPE=str(record.type)
            DICT[record.qualifiers["gene"][0]][TYPE]={}
            
            #Given that it could have multiple exons
            if len(record.sub_features)>0:
                POSITION_PAIRS=[]
                for sub in record.sub_features:
                    POSITION_PAIRS=POSITION_PAIRS+[[sub.location.nofuzzy_start, sub.location.nofuzzy_end]]
                DICT[record.qualifiers["gene"][0]][TYPE]["POSITION_PAIRS"]=POSITION_PAIRS
                
            else:
                POSITION_PAIRS=[[record.location.nofuzzy_start, record.location.nofuzzy_end]]
                DICT[record.qualifiers["gene"][0]][TYPE]["POSITION_PAIRS"]=POSITION_PAIRS
    
    for gene in DICT.keys():
        
        #Make sure to add chromosomal information    
        #Write feature available based on priority mRNA>CDS>exon per gene
        if "mRNA" in DICT[gene].keys():
            
            for pair in DICT[gene]["mRNA"]["POSITION_PAIRS"]:
                FILE_OUT1.write("\n"+gene+"\t"+ CHROMOSOME +"\t"+
                        str(DICT[gene]["STRAND"])+"\t"+str(DICT[gene]["START"])+"\t"+str(DICT[gene]["END"]))
            
                FILE_OUT1.write("\t"+"mRNA")
                
                FILE_OUT1.write("\t"+str(pair[0])+
                            "\t"+str(pair[1]))
        
        elif "CDS" in DICT[gene].keys():
            
            for pair in DICT[gene]["CDS"]["POSITION_PAIRS"]:
                FILE_OUT1.write("\n"+gene+"\t"+ CHROMOSOME +"\t"+
                        str(DICT[gene]["STRAND"])+"\t"+str(DICT[gene]["START"])+"\t"+str(DICT[gene]["END"]))
            
                FILE_OUT1.write("\t"+"CDS")
                
                FILE_OUT1.write("\t"+str(pair[0])+
                            "\t"+str(pair[1]))
        
        elif "exon" in DICT[gene].keys():
            
            for pair in DICT[gene]["exon"]["POSITION_PAIRS"]:
                FILE_OUT1.write("\n"+gene+"\t"+ CHROMOSOME +"\t"+
                        str(DICT[gene]["STRAND"])+"\t"+str(DICT[gene]["START"])+"\t"+str(DICT[gene]["END"]))
            
                FILE_OUT1.write("\t"+"exon")
                
                FILE_OUT1.write("\t"+str(pair[0])+
                            "\t"+str(pair[1]))

FILE_OUT1.close()
