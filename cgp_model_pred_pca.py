#cgp.cor.AUC_class_pred.py
#Function to predict class from cgp.cor.AUC.class

##################################################################################################################

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
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

class LogisticRegression(object):

    def __init__(self, input, n_in, n_out, rng, W=None, b=None):

        if W is None:
            W_values = numpy.asarray(
                    rng.uniform(
                        low=-numpy.sqrt(6. / (n_in + n_out)),
                        high=numpy.sqrt(6. / (n_in + n_out)),
                        size=(n_in, n_out)
                    ),
                    dtype=theano.config.floatX
                )
            W = theano.shared(value=W_values, name='W', borrow=True)

        if b is None :
            b = theano.shared(
                value=numpy.zeros(
                    (n_out,),
                    dtype=theano.config.floatX
                ),
                name='b',
                borrow=True
            )

        self.W = W
        self.b = b

        self.p_y_given_x = T.nnet.softmax(T.dot(input, self.W) + self.b)

        self.y_pred = T.argmax(self.p_y_given_x, axis=1)

        self.params = [self.W, self.b]

        #self.input = input

    def pred(self, y):
        """Returns prediction only
        """
        return(self.y_pred)

    def negative_log_likelihood(self, y):

        return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]), y])
        #return T.nnet.categorical_crossentropy(self.p_y_given_x, y).mean()

    def errors(self, y):

        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                'y should have the same shape as self.y_pred',
                ('y', y.type, 'y_pred', self.y_pred.type)
            )
        # check if y is of the correct datatype

        return T.mean(T.neq(self.y_pred, y))
        # if y.dtype.startswith('int'):
        #     # the T.neq operator returns a vector of 0s and 1s, where 1
        #     # represents a mistake in prediction
        #     return T.mean(T.neq(self.y_pred, y))
        # else:
        #     raise NotImplementedError()

    def diff(self, y):
        setdiff = self.y_pred - y
        return T.mean(abs(setdiff))

def drop(input, rng, p=0.5):
    """
    :type input: numpy.array
    :param input: layer or weight matrix on which dropout resp. dropconnect is applied

    :type p: float or double between 0. and 1.
    :param p: p probability of NOT dropping out a unit or connection, therefore (1.-p) is the drop rate.

    """
    srng = RandomStreams(rng.randint(999999))

    mask = srng.binomial(n=1, p=1.-p, size=input.shape)
    return input * T.cast(mask, theano.config.floatX) / (1.-p)

def prelu(x, alpha):
    return theano.tensor.switch(x<0, alpha*x, x)

class HiddenLayer(object):
    def __init__(self, rng, is_train, input, n_in, n_out, W=None, b=None, alpha=None,
                 activation=T.tanh, p=0.5, dropout=False):

        self.input = input

        if W is None:
            W_values = numpy.asarray(
                rng.uniform(
                    low=-numpy.sqrt(6. / (n_in + n_out)),
                    high=numpy.sqrt(6. / (n_in + n_out)),
                    size=(n_in, n_out)
                ),
                dtype=theano.config.floatX
            )

            if activation == theano.tensor.nnet.sigmoid:
                W_values *= 4

            W = theano.shared(value=W_values, name='W', borrow=True)

        if b is None:
            b_values = numpy.zeros((n_out,), dtype=theano.config.floatX)
            b = theano.shared(value=b_values, name='b', borrow=True)

        if alpha is None:
            alpha_value = numpy.full((n_out), .1,  dtype=theano.config.floatX)
            alpha = theano.shared(value=alpha_value, name='alpha', borrow=True)

        self.W = W
        self.b = b
        self.alpha = alpha

        lin_output = T.dot(input, self.W) + self.b
        output = (
            lin_output if activation is None
            else activation(lin_output, self.alpha)
        )

        #Is droput necessary?
        if dropout==True:
            train_output = drop(output, rng=rng, p=p)
            self.output = T.switch(T.neq(is_train, 0), train_output, output)
        else:
            self.output = output

        # parameters of the model
        self.params = [self.W, self.b, self.alpha]

class MLP(object):
    def __init__(self, rng, input, is_train, n_in, n_hidden, n_out, p=0.5, dropout=False, input_p=0.1): #, batch_size=20):

        #Need input dropout layer
        if input_p!=None:
            self.input_layer = drop(input, rng=rng, p=input_p)
            self.input_layer = T.switch(T.neq(is_train, 0), self.input_layer, input)
        else:
            self.input_layer=input

        self.layer_0 = HiddenLayer(
            rng=rng,
            input=self.input_layer,
            n_in=n_in,
            n_out=n_hidden[0],
            activation=prelu,
            is_train=is_train,
            p=p,
            dropout=dropout
        )

        self.params = self.layer_0.params

        #Add more layers accordingly
        layer_number = 1
        if len(n_hidden)>1:

            for layer in n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "layer_" + str(layer_number-1)).output,
                                                    n_in=n_hidden[layer_number-1],
                                                    n_out=n_hidden[layer_number],
                                                    activation=prelu,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
                                                )

                setattr(self, "layer_" + str(layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "layer_" + str(layer_number)).params

                layer_number = layer_number + 1


        # The logistic regression layer gets as input the hidden units
        # of the hidden layer
        self.logRegressionLayer = LogisticRegression(
            input=getattr(self, "layer_" + str(layer_number-1)).output,
            n_in=n_hidden[layer_number-1],
            n_out=n_out,
            rng=rng #,batch_size=batch_size
        )
        self.params = self.params + self.logRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.layer_0.W).sum() + abs(self.logRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.layer_0.W ** 2).sum() + (self.logRegressionLayer.W ** 2).sum()
        )
        #

        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )

        self.errors = self.logRegressionLayer.errors
        self.pred = self.logRegressionLayer.pred
        self.diff = self.logRegressionLayer.diff
        # self.errors = self.linearRegressionLayer.errors
        # self.loss = self.linearRegressionLayer.loss
        # self.NRMSE = self.linearRegressionLayer.NRMSE
        # self.pred = self.linearRegressionLayer.pred

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

def shared_drug_dataset_pred(drug_data, integers=True, target="AUC"):

    data_x=drug_data.iloc[:,1:]

    # if target=="AUC":
    #     data_y=list(drug_data.NORM_AUC)
    # elif target=="pIC50":
    #     data_y=list(drug_data.NORM_pIC50)
    # elif target=="IC50":
    #     data_y=list(drug_data.NORM_IC50)
    # elif target=="LIVED":
    #     data_y=list(drug_data.LIVED)
    # elif target=="CLASS":
    #     data_y=list(drug_data.CLASS)
    data_y = list(drug_data[0])

    shared_x = theano.shared(numpy.asarray(data_x, dtype=theano.config.floatX), borrow=True)
    shared_y = theano.shared(numpy.asarray(data_y, dtype=theano.config.floatX), borrow=True)

    if integers==True:
        return shared_x, T.cast(shared_y, 'int32')
    else:
        return shared_x, shared_y

def model_prediction(MODEL_FILE, test_drug_x):

    model_obj = []

    for l in  cPickle.load(open(MODEL_FILE, "rb")):
        model_obj = model_obj + [l]

    #Perform prediction
    n_layers = len(model_obj)-2 #minus loglayer and rng

    input = test_drug_x
    for l in xrange(1, n_layers+1):

        input = prelu(
                        T.dot(input, model_obj[l].W) + model_obj[l].b ,
                        model_obj[l].alpha
                      )

    log_layer = T.nnet.softmax(
                                T.dot(input, model_obj[0].W) + model_obj[0].b
                              )
    prediction = T.argmax(log_layer, axis=1)

    return prediction.eval()

##################################################################################################################
#Load class file to make predictions for each drugs
IN_FILE = sys.argv[1]
BASE_FILE = sys.argv[2]
PCA_FILE = sys.argv[3]
MODEL_FILE = sys.argv[4]
METRIC = "AUC"
n_pcas = 700

all_cgp = pd.read_csv(IN_FILE, sep="\t")
all_drugs = list(set(all_cgp.DRUG))

used_feat = pd.read_csv(BASE_FILE, sep="\t", nrows=100)
used_feat = [c for c in used_feat]

with open(PCA_FILE, "rb") as md:
    pca_model = cPickle.load(md)

#Obtain predictions per cancer type
rotation = pca_model.components_[:n_pcas]
rotation = rotation.transpose()

#Predict and write tables for each drug
OUT_FOLDER = "/tigress/zamalloa/RESULTS/TCGA.TRAINING/" #For tigress
OUT_FOLDER = "/home/zamalloa/Documents/FOLDER/RESULTS/TCGA.TRAINING/" #For Lab

FILE_OUT_val = open(OUT_FOLDER + "cgp_auc_class_predction", "w")
FILE_OUT_val.write("DRUG" + "\t" + "ACTUAL" + "\t" + "PREDICTED")

COUNT = 0.0
for drug in all_drugs:

    target_table = all_cgp[all_cgp.DRUG==drug]
    target_labels = target_table["NORM.AUC"]
    target_labels = pd.DataFrame({"NORM_AUC": list(target_labels)})

    target_table = target_table[used_feat]
    print(target_table.shape)
    #target_table = scale(target_table) #Supposedly scaled

    target_table = numpy.dot(target_table, rotation)
    target_table = scale(target_table)
    target_table = pd.DataFrame(target_table)

    target_table = pd.concat([target_labels, target_table], axis=1)
    print(target_table.shape)

    test_drug_x, test_drug_y = shared_drug_dataset_pred(target_table, integers=True, target=METRIC)

    prediction = model_prediction(MODEL_FILE, test_drug_x)
    actual = test_drug_y.get_value()
    print(actual[:5])

    for n in xrange(len(actual)):
        FILE_OUT_val.write("\n" + drug + "\t" + str(actual[n]) + "\t" + str(prediction[n]) )

    COUNT = COUNT+1
    print(COUNT/len(all_drugs))

FILE_OUT_val.close()

print("Done predicting!!")
