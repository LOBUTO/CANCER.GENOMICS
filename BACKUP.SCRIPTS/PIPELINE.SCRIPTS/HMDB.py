#PIPELINE FOR PARSING HMDB XML FILES
#Requires:
	#HMDB folder where XML files are contained
	#Text file containing list of drugs to be filtered out of XML files
	#LOGS folder where to deposit temporal files created by BABEL
	#Name of output file
	#TC filter level (0.0 - 1.0)

import subprocess, pickle, sys, os
import xml.etree.ElementTree as ET

#Folder input where XML files are at
FOLDER=sys.argv[1]

#Text file where drugs to be filtered against are at - File contains one drug per line
FILE_IN1=open(sys.argv[2])
DRUG_FILTER=[DRUG.upper() for DRUG in FILE_IN1.read().splitlines()]
FILE_IN1.close()

#Logs folder
LOGS=sys.argv[3]

#Output file
FILE_OUT=open(sys.argv[4], "w")

#Get XML file names
FILES=subprocess.check_output("ls", cwd=FOLDER).splitlines()

#Store in Dictionary
DICT={}

#STEP 1 - PARSE XML FILES FOR METABOLITES
for record in FILES:
	print record
	INFILE=os.path.join(FOLDER, record)

	#Process XML
	ROOT=ET.parse(INFILE).getroot()

	ORIGIN=[X.text for X in ROOT.iter("origin")]
	ORIGIN_CHECK=set(["Endogenous", "Food"]) | set(ORIGIN)
	
	#Check if it proceeds from Food Source or if it is Endogenous
	if (len(ORIGIN_CHECK)<=2) and (len(ORIGIN)>0):
	
		#Get Name
		NAME=ROOT.find("name").text.strip()		
		
		#Get Synonyms
		SYNONYMS=[SYN.text for SYN in ROOT.findall("synonyms/synonym")]
			
		#Get Smile - KEEP IN MIND THAT SOME MAY HAVE SMILES="None"
		SMILES=str(ROOT.find("smiles").text).strip()
		
		#Get KEGG id - KEEP IN MIND THAT SOME MAY HAVE kegg_id="None"
		KEGG_ID=str(ROOT.find("kegg_id").text).strip()
		
		#Get UNIPROTS
		UNIPROTS=[X.text.strip() for X in ROOT.findall("protein_associations/protein/uniprot_id")]

		#Get Gene names
		GENES=[str(Y.text).strip() for Y in ROOT.findall("protein_associations/protein/gene_name")]	
		GENES=list(set(filter(lambda x: x!="None", GENES)))
		
		#Append to dictionary per each depending on wether it has a gene name or not, if not, then append uniprot
		#ONLY APPEND TO DICTIONARY IF IT HAS A PROTEIN ASSOCIATION RECORD
		RESULTS=zip(UNIPROTS,GENES)
		if len(RESULTS)>0:
			DICT[NAME]={"SMILES":"", "UNIPROTS":[], "GENES":[], "SYNONYMS":[]}
			DICT[NAME]["SMILES"]=SMILES		#String
			DICT[NAME]["KEGG_ID"]=KEGG_ID #Keep in mind this is a string and not a list
			DICT[NAME]["GENES"]=[gene[1] for gene in RESULTS if gene[1]!="None"]
			DICT[NAME]["UNIPROTS"]=[uniprot[0] for uniprot in RESULTS if uniprot[1]=="None"]
			DICT[NAME]["SYNONYMS"]=SYNONYMS

print len(DICT)
#STEP 2 - FILTER FOR DRUGS - BE CAREFUL WITH THIS STEP, MAKE SURE YOU ARE AWARE OF WHAT YOU ARE FILTERING OUT

for METABOLITE in DICT.keys():
	if len(set([MET.upper() for MET in [METABOLITE]+DICT[METABOLITE]["SYNONYMS"]]) & set(DRUG_FILTER))>0:
		print METABOLITE
		del DICT[METABOLITE]
print len(DICT)

#STEP 3 - FILTER FOR SYNONYMS BY SMILES
#This step filters out over represented molecules such as triglycerides that have equal tanimoto coefficients
for metabolite in DICT.keys():
	print len(DICT)
	#Check that it is still in dictionary after filtering due to iteration and that it has a smiles code
	if DICT.has_key(metabolite) and DICT[metabolite]["SMILES"]!="None":
		
		#Prep files
		QUERY_FILE=open(os.path.join(LOGS, "TEMP1.smi"), "w")
		QUERY_FILE.write(DICT[metabolite]["SMILES"])
		QUERY_FILE.close()
		
		LIST_FILE=open(os.path.join(LOGS, "TEMP2.smi"), "w")
		OTHER_SMILES=[SMI for SMI in DICT.keys() if SMI!=metabolite and DICT[SMI]["SMILES"]!="None"]
		for record in OTHER_SMILES:
			LIST_FILE.write(DICT[record]["SMILES"] + "\t" + record + "\n")
		LIST_FILE.close()
		
		#Run BABEL - Verify Babel path to be correct, FP2 used
		BABEL=subprocess.Popen(["-c", "/usr/local/bin/babel %s %s -ofpt -xfFP2"%(os.path.join(LOGS,"TEMP1.smi"), os.path .join(LOGS,"TEMP2.smi"))],
		stdout=subprocess.PIPE, shell=True)
				
		BABEL=BABEL.communicate()[0].splitlines()
		
		#Filter out troublesome lines
		BABEL=filter(lambda x:len(x)>3 and x!="Possible superstructure of first mol", BABEL)
		
		#Rank BABEL output in descending order
		RANK=[[float(X.split("=")[1].strip()), X.split("Tanimoto")[0].strip().strip(">")] for X in BABEL]
		RANK=sorted(RANK, reverse=True)
		
		#Return metabolites that have a Tanimoto coefficient equal or greater to desired FILTER, 
		#molecules that do not pass filter (have greater or equal TC when compared to key) will be removed
		FILTER=[TC[1] for TC in RANK if TC[0]>=float(sys.argv[5])] #TC FILTER TO BE APPLIED
		
		#Remove TC-synonyms from dictionary
		for TC in FILTER:
			if TC in DICT:
				del DICT[TC]

#STEP 4 - REMOVE WATER MOLECULE
#Its appearance in the database could be an artifact of pdb structures produced displaying water molecules
if ("Water" in DICT)==True:
	del DICT["Water"]

#STEP 5 - WRITE TO TABLE (Rather than pickle out)
#Only write to table if it has genes associated with it
#KEEP IN MIND, we are only writing those with KEGG identifiers!!!!!
FILE_OUT.write("METABOLITE"+"\t"+ "KEGG_ID"+"\t"+ "Hugo_Symbol")
for metabolite in DICT.keys():
	
	if (len(DICT[metabolite]["GENES"])>0  and DICT[metabolite]["KEGG_ID"]!="None") :
		MET_GENES=DICT[metabolite]["GENES"] 
		for gene in MET_GENES:
			
			FILE_OUT.write("\n"+metabolite+"\t"+ DICT[metabolite]["KEGG_ID"] + "\t" + gene)

FILE_OUT.close()
		
#Pickle out
#pickle.dump(DICT, PICKLE_OUT)
#PICKLE_OUT.close()


		
