import numpy as np
import pylab as pl
import scipy.stats

DATA_LOGP=np.loadtxt("LIGAND_OUTPUT/LOGP_SINGLE_SET")
DATA_MW=np.loadtxt("LIGAND_OUTPUT/MW_SINGLE_SET")

print scipy.stats.spearmanr(DATA_LOGP[:,0], DATA_LOGP[:,1])
pl.plot(DATA_LOGP[:,0], DATA_LOGP[:,1], "ro")
pl.show()

print scipy.stats.spearmanr(DATA_MW[:,0], DATA_MW[:,1])
pl.plot(DATA_MW[:,0], DATA_MW[:,1], "g*")
pl.show()