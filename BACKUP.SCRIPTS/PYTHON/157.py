#081514 
#GET CUMMULATIVE GENE COUNT FOR NEGS
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
BRCA_vp=sorted(BRCA_vp, key=lambda x: x[1])
for i in BRCA_vp[:10]:
    print i
    
#Make dict per gene of negatively affected genes
NEG_DICT=dict((BRCA_vp[x][0], BRCA_NEGS[x]) for x in range(len(BRCA_vp)))
print len(NEG_DICT)

print NEG_DICT["A1CF"][:10]
print NEG_DICT["hsa-let-7f-2"][:10]

#Set up parallelization
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 
pool = ThreadPool() 

#Get cummulative unions
def neg_cum(v_count):
    print v_count
    CAUSAL_GENES=[X[0] for X in BRCA_vp[:v_count+1]]
    NEG_GENS=set([X for Y in CAUSAL_GENES for X in NEG_DICT[Y]])
    NEG_GENS_LEN=len(NEG_GENS)
    return [CAUSAL_GENES[-1], NEG_GENS_LEN]

tic=timeit.default_timer()
NEG_CUM_COUNTS=map(neg_cum,range(1000))
toc=timeit.default_timer()
pool.close()
pool.join()

print toc-tic
print len(NEG_CUM_COUNTS)
for i in NEG_CUM_COUNTS[:10]:
    print i