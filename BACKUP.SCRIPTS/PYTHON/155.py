#####081314#####

import itertools
import timeit

FILE_IN1=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p")
BRCAvp=[x.split("\t") for x in FILE_IN1.read().splitlines()[1:]]
FILE_IN1.close()
print len(BRCAvp)
for i in BRCAvp[:4]:
    print i

FILE_IN2=open("PIPELINES/METABOLIC.DRIVERS/TABLES/BRCA/080914.BRCA.Table.v.p.neg")
BRCAneg=[x.split() for x in FILE_IN2.read().splitlines()]
FILE_IN2.close()

print len(BRCAneg)
for i in BRCAneg[:9]:
    print i[:10]

#Make dict per gene of negatively affected genes
NEG_DICT=dict((BRCAvp[x][0], BRCAneg[x]) for x in range(len(BRCAvp)))
print len(NEG_DICT)

print NEG_DICT["A1CF"][:10]
print NEG_DICT["hsa-let-7f-2"][:10]

#Get pairwise combinations of all genes in dict
NEG_PAIRS=list(itertools.combinations(NEG_DICT.keys(),2))
print len(NEG_PAIRS)
print NEG_PAIRS[:5]

#Set up parallelization
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool 
pool = ThreadPool() 

#Get intersection of pairwises
def neg_int(PAIR):
    CALLS=[PAIR[0], PAIR[1], len(set.intersection(set(NEG_DICT[PAIR[0]]), set(NEG_DICT[PAIR[1]])))]
    return CALLS

tic=timeit.default_timer()
NEG_INTERSECTS_20=map( neg_int, NEG_PAIRS[1:20000000])
toc=timeit.default_timer()
pool.close()
pool.join()

print toc-tic
print NEG_INTERSECTS[1:10]
