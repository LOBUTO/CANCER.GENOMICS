# ctrp_morgan.py
# Extract morgan features for ctrp compounds

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit import Chem

# Load files
file_in = open("/Users/jzamalloa/Documents/Rotation/DATABASES/CTRP/CTRPv2.0_2015_ctd2_ExpandedDataset/v20.meta.per_compound.txt")
drugs   = dict((x.split("\t")[1], x.split("\t")[10]) for x in file_in.read().splitlines()[1:])
file_in.close()

############################################################################################################################################
# OBTAIN MORGAN BITS
# Across radii and bit count
for r in [1, 2, 4, 8, 12, 16]:
	for b in [256, 512, 1024, 2048]:

		FILE_OUT = open( ("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CTRP_MORGAN_BITS_r_%i_b_%i.txt") % (r, b), "w")
		FILE_OUT.write("radius" + "\t" + "bits" + "\t" + "Compound" + "\t" + "position" + "\t" + "value")
		FILE_OUT.close()

for d in drugs.keys():
	print(d)
	smiles = Chem.MolFromSmiles(drugs[d])

	if smiles is not None:

		for r in [1, 2, 4, 8, 12, 16]:
			for b in [256, 512, 1024, 2048]:
			    
			    mcf     = AllChem.GetMorganFingerprintAsBitVect(smiles, r, nBits=b)
			    mcf_bit = mcf.ToBitString()
			    
			    # Write out
			    with open( ("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CTRP_MORGAN_BITS_r_%i_b_%i.txt") % (r, b), "a") as FILE_OUT:
				    for j in xrange(len(mcf_bit)):
				        	FILE_OUT.write("\n" + str(r) + "\t" + str(b) + "\t" + d + "\t" + str(j+1) + "\t" + mcf_bit[j])

print("Done writing morgan bits")

############################################################################################################################################
# OBTAIN MORGAN COUNTS
FILE_OUT = open("/Users/jzamalloa/Documents/Rotation/PIPELINES/OBJECTS/CTRP_MORGAN_COUNTS.txt", "w")
FILE_OUT.write("radius" + "\t" + "Compound" + "\t" + "substructure" + "\t" + "value")

for d in drugs.keys():
	print(d)
	smiles = Chem.MolFromSmiles(drugs[d])

	if smiles is not None:

		for r in [1, 2, 4, 8, 12, 16]:
			    
		    mcf        = AllChem.GetMorganFingerprint(smiles, r)
		    mcf_counts = mcf.GetNonzeroElements()
		    
		    # Write out
		    for sub in mcf_counts.keys():
        		FILE_OUT.write("\n" + str(r) + "\t" + d + "\t" + str(sub) + "\t" + str(mcf_counts[sub]))		    

FILE_OUT.close()

print("Done writing morgan counts")