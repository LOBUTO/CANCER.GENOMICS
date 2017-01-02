# cgp_morgan.py
# Extract morgan features for cgp drugs

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem

# Load cgp drugs
file_in = open("/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/pubchem_id_2016_filter_smiles.txt")
drugs_pubmed = dict((x.split("\t")[0], x.split("\t")[1]) for x in file_in.read().splitlines()[1:])
file_in.close()

file_in = open("/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/DRUG_SMILES_ALL.txt")
drugs_all = dict((x.split("\t")[0], x.split("\t")[1]) for x in file_in.read().splitlines()[1:])
file_in.close()

for i in drugs_all.keys():
	if i in drugs_pubmed.keys():
		del drugs_all[i]

drugs_pubmed.update(drugs_all)

############################################################################################################################################
# OBTAIN MORGAN BITS
# Across radii and bit count
# for r in [1, 2, 4, 8, 12, 16]:
# 	for b in [256, 512, 1024, 2048]:

# 		FILE_OUT = open( ("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CGP_MORGAN_BITS_r_%i_b_%i.txt") % (r, b), "w")
FILE_OUT = open("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CGP_MORGAN_BITS.txt", "w")
FILE_OUT.write("radius" + "\t" + "bits" + "\t" + "Compound" + "\t" + "bit_pos" + "\t" + "value")
FILE_OUT.close()

for d in drugs_pubmed.keys():
	
	smiles = Chem.MolFromSmiles(drugs_pubmed[d])

	if smiles is not None:

		for r in [1, 2, 4, 8, 12, 16]:
			for b in [256, 512, 1024, 2048]:
			    
			    mcf     = AllChem.GetMorganFingerprintAsBitVect(smiles, r, nBits=b)
			    mcf_bit = mcf.ToBitString()
			    
			    # Write out
			    with open("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CGP_MORGAN_BITS.txt", "a") as FILE_OUT:
				    for j in xrange(len(mcf_bit)):
				        	FILE_OUT.write("\n" + str(r) + "\t" + str(b) + "\t" + d + "\t" + "mcf_"+str(j+1) + "\t" + mcf_bit[j])

print("Done writing morgan bits")

############################################################################################################################################
# OBTAIN MORGAN COUNTS
FILE_OUT = open("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CGP_MORGAN_COUNTS.txt", "w")
FILE_OUT.write("radius" + "\t" + "Compound" + "\t" + "Substructure" + "\t" + "Counts")

for d in drugs_pubmed.keys():
	
	smiles = Chem.MolFromSmiles(drugs_pubmed[d])

	if smiles is not None:

		for r in [1, 2, 4, 8, 12, 16]:
			    
		    mcf        = AllChem.GetMorganFingerprint(smiles, r)
		    mcf_counts = mcf.GetNonzeroElements()
		    
		    # Write out
		    for sub in mcf_counts.keys():
        		FILE_OUT.write("\n" + str(r) + "\t" + d + "\t" + str(sub) + "\t" + str(mcf_counts[sub]))		    

FILE_OUT.close()

print("Done writing morgan counts")
