####081814 
#Filter genes in v(p) for subset of negatively differentiated genes
import timeit

#Load files
FILE_IN1=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p.neg")
BRCA_NEGS=[X.split() for X in FILE_IN1.read().splitlines()]
FILE_IN1.close()

FILE_IN2=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p")
BRCA_vp=[x.split("\t")[:2] for x in FILE_IN2.read().splitlines()[1:]]
FILE_IN2.close()

for i in BRCA_vp[:10]:
    print i

print "SORTED"
#Sort BRCA_vp by v(p)
BRCA_vp=sorted(BRCA_vp, key=lambda x: x[1], reverse=True)
for i in BRCA_vp[:10]:
    print i
    
#Make dict per gene of negatively affected genes
NEG_DICT=dict((BRCA_vp[x][0], BRCA_NEGS[x]) for x in range(len(BRCA_vp)))
print len(NEG_DICT)

print NEG_DICT["A1CF"][:10]
print NEG_DICT["hsa-let-7f-2"][:10]

#Check speed for first 100
tic=timeit.default_timer()
X=0
SORTED_GENES=[x[0] for x in BRCA_vp]

for count in range(len(SORTED_GENES)):
    
    gene=SORTED_GENES[count]
    if gene in NEG_DICT:
        
        SELF_GENES=set(NEG_DICT[gene])
        OTHER_GENES=SORTED_GENES[count+1:]
        OTHER_GENES=filter(lambda x:x in NEG_DICT, OTHER_GENES)
        
        #Filter out other genes if affected genes are all contained within gene being studied with 0.25% error
        OTHER_CONTAINED=filter(lambda x: 
                               len(set.intersection(set(NEG_DICT[x]), SELF_GENES)) in range(int(len(NEG_DICT[x])*0.975), len(NEG_DICT[x])+1), 
                               OTHER_GENES)
        
        #Delete genes in they are contatined by top gene
        for other in OTHER_CONTAINED:
            print other
            del NEG_DICT[other]
            
    X=X+1
    print X

toc=timeit.default_timer()        
print toc-tic
print len(NEG_DICT)
