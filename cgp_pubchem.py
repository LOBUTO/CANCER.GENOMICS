# cgp_pubchem.py
# Extracts cgp pubchem smiles using the pubchem IDs

import urllib
import pandas as pd

file_out = "/Users/jzamalloa/Documents/Rotation/DATABASES/CANCERRXGENE/pubchem_id_2016_filter_smiles.txt"

# Load pubchem cgp file
file_in  = pd.read_csv("~/Documents/Rotation/DATABASES/CANCERRXGENE/pubchem_id_2016_filter.txt", sep="\t")
p_dict   = dict((file_in.Compound[i], file_in.Pub_ID[i]) for i in xrange(file_in.shape[0]))

# Extract SMILES through API
out = open(file_out, "w")
out.write("Compound" + "\t" + "SMILES")
for i in p_dict.keys():

	pub_id = p_dict[i]

	smiles       = urllib.urlopen(("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/%s/property/CanonicalSMILES/txt")%pub_id ).read()
	smiles_split = smiles.split("\n")

	if smiles != "\n" and len(smiles_split)==2:
		print(smiles_split[0])
		out.write("\n" + i + "\t" + smiles_split[0])

out.close()

print("Done writing")
