#Run cythons

import cython159

TEST=cython159.vpCUMPATIENTDEL("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p.SNS.0",
                      "PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081914.BRCA.Table.1.GENE.PATIENT")

print len(TEST)

FILE_OUT=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/081814.BRCA.Table.v.p.filtered.GENE.PATIENT.SNS.0","w")
FILE_OUT.write("Hugo_Symbol"+"\t"+"v.PROTEIN")
for i in TEST:
    FILE_OUT.write("\n"+i[0]+"\t"+str(i[1]))
FILE_OUT.close()
