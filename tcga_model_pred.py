import pandas as pd
import os
import sys
import timeit
import numpy
import theano
import cPickle
import theano.tensor as T
from theano.sandbox.rng_mrg import MRG_RandomStreams as RandomStreams
from scipy.stats import pearsonr

def shared_drug_dataset_IC50(drug_file, integers=True, target="AUC"):

    drug_data=pd.read_csv(drug_file, sep="\t")

    data_x=drug_data.iloc[:,3:]

    if target=="AUC":
        data_y=list(drug_data.NORM_AUC)
    elif target=="pIC50":
        data_y=list(drug_data.NORM_pIC50)
    elif target=="IC50":
        data_y=list(drug_data.NORM_IC50)
    elif target=="LIVED":
        data_y=list(drug_data.LIVED)
    elif target=="CLASS":
        data_y=list(drug_data.CLASS)

    shared_x = theano.shared(numpy.asarray(data_x, dtype=theano.config.floatX), borrow=True)
    shared_y = theano.shared(numpy.asarray(data_y, dtype=theano.config.floatX), borrow=True)

    if integers==True:
        return shared_x, T.cast(shared_y, 'int32')
    else:
        return shared_x, shared_y

#Load file to predict
IN_FILE = sys.argv[1]
METRIC = "AUC" #assumed!!!

FEAT_MATRIX = shared_drug_dataset_IC50(IN_FILE, integers=True, target=METRIC)

#Load model
MODEL_FILE = sys.argv[2]
