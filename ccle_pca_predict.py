# ccle_pca_predict.py

import pandas as pd
import os
import sys
import timeit
import numpy as np
import theano
import cPickle
import theano.tensor as T
from theano.sandbox.rng_mrg import MRG_RandomStreams as RandomStreams
#from theano.tensor.shared_randomstreams import RandomStreams
#from theano.sandbox.cuda.rng_curand import CURAND_RandomStreams as RandomStreams

class LogisticRegression(object):

    def __init__(self, input, n_in, n_out, rng, W=None, b=None):

        if W is None:
            W_values = np.asarray(
                    rng.uniform(
                        low=-np.sqrt(6. / (n_in + n_out)),
                        high=np.sqrt(6. / (n_in + n_out)),
                        size=(n_in, n_out)
                    ),
                    dtype=theano.config.floatX
                )
            W = theano.shared(value=W_values, name='W', borrow=True)

        if b is None :
            b = theano.shared(
                value=np.zeros(
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

        #return -T.mean(T.log(self.p_y_given_x)[T.arange(y.shape[0]), y])
        return T.nnet.categorical_crossentropy(self.p_y_given_x, y).mean()

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

def drop(input, rng, p=0.5):
    """
    :type input: np.array
    :param input: layer or weight matrix on which dropout resp. dropconnect is applied

    :type p: float or double between 0. and 1.
    :param p: p probability of NOT dropping out a unit or connection, therefore (1.-p) is the drop rate.

    """
    srng = RandomStreams(rng.randint(999999))

    mask = srng.binomial(n=1, p=1.-p, size=input.shape)
    return input * T.cast(mask, theano.config.floatX) / (1.-p)

def relu(x):
    return theano.tensor.switch(x<0, 0, x)

def prelu(x, alpha):
    return theano.tensor.switch(x<0, alpha*x, x)

class HiddenLayer(object):
    def __init__(self, rng, is_train, input, n_in, n_out, W=None, b=None, alpha=None,
                 activation=T.tanh, p=0.5, dropout=False):

        self.input = input

        if W is None:
            W_values = np.asarray(
                rng.uniform(
                    low=-np.sqrt(6. / (n_in + n_out)),
                    high=np.sqrt(6. / (n_in + n_out)),
                    size=(n_in, n_out)
                ),
                dtype=theano.config.floatX
            )

            if activation == theano.tensor.nnet.sigmoid:
                W_values *= 4

            W = theano.shared(value=W_values, name='W', borrow=True)

        if b is None:
            b_values = np.zeros((n_out,), dtype=theano.config.floatX)
            b = theano.shared(value=b_values, name='b', borrow=True)

        if alpha is None:
            alpha_value = np.full((n_out), .1,  dtype=theano.config.floatX)
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

def rescale_weights(params, incoming_max):
    incoming_max = np.cast[theano.config.floatX](incoming_max)
    for p in params:
        w = p.get_value()
        w_sum = (w**2).sum(axis=0)
        w[:, w_sum>incoming_max] = w[:, w_sum>incoming_max] * np.sqrt(incoming_max) / w_sum[w_sum>incoming_max]
        p.set_value(w)

class MLP(object):
    def __init__(self, rng, input, is_train, n_in, n_hidden, n_out, p=0.5, dropout=False, input_p=0.1): #, batch_size=20):

        #Need input dropout layer
        if input_p!=None:
            self.input_layer = drop(input, rng=rng, p=input_p)
            self.input_layer = T.switch(T.neq(is_train, 0), self.input_layer, input)
        else:
            self.input_layer=input

        param_to_scale = [] #To scale weights to square length of 15

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
        param_to_scale = param_to_scale + [self.layer_0.params[0]]

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

                param_to_scale = param_to_scale + [getattr(self, "layer_" + str(layer_number)).params[0]]

                layer_number = layer_number + 1


        # The logistic regression layer gets as input the hidden units
        # of the hidden layer

        self.logRegressionLayer = LogisticRegression(
            input=getattr(self, "layer_" + str(layer_number-1)).output,
            n_in=n_hidden[layer_number-1],
            n_out=n_out,
            rng=rng
        )
        self.params = self.params + self.logRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.layer_0.W).sum() + abs(self.logRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.layer_0.W ** 2).sum() + (self.logRegressionLayer.W ** 2).sum()
        )

        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )

        self.errors = self.logRegressionLayer.errors
        self.pred = self.logRegressionLayer.pred
        self.param_to_scale = param_to_scale
        # self.diff = self.logRegressionLayer.diff
        # self.errors = self.linearRegressionLayer.errors
        # self.loss = self.linearRegressionLayer.loss
        # self.NRMSE = self.linearRegressionLayer.NRMSE
        # self.pred = self.linearRegressionLayer.pred

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

def shared_drug_dataset_IC50(drug_data, integers=True):

    data_x=drug_data.iloc[:,3:]

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
    # elif target=="PCA":
    #     data_y=list(drug_data.PCA)

    data_y = list(drug_data.pic50_class)

    shared_x = theano.shared(np.asarray(data_x, dtype=theano.config.floatX), borrow=True)
    shared_y = theano.shared(np.asarray(data_y, dtype=theano.config.floatX), borrow=True )

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

####################################################################################################################################################################################################
# LOAD FILES

drug  = sys.argv[1]
m_pca = sys.argv[2]
g_pca = sys.argv[3]

file_to_model = "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/" + "PRE_CCLE_PCA." + drug + "_" + m_pca + "_" + g_pca
model_file    = "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_RESULTS/new_cgp_pca_class_model_" + drug + "_" + g_pca + "_" + m_pca + ".pkl"
####################################################################################################################################################################################################
# EXECUTE

ccle_data      = pd.read_csv(file_to_model, sep="\t")

ccle_x, ccle_y = shared_drug_dataset_IC50(ccle_data, integers = True)

predictions    = model_prediction(file_to_model, ccle_x)

actual         = ccle_y.eval()

####################################################################################################################################################################################################
# WRITE OUT
file_out = open("/home/zamalloa/Documents/FOLDER/CGP_FILES/CCLE_PRED/cgp_pca_model_M" + m_pca + "_G_" + g_pca + "_D_" + drug, "w")
file_out.write("M_PCA" + "\t" + "G_PCA" + "\t" + "DRUG" + "\t" + "ACTUAL" + "\t" + "PREDICTION")

for n in xrange(len(actual)):

    file_out.write("\n" + m_pca + "\t" + g_pca + "\t" + drug + "\t" + actual[n] + "\t" + predictions[n])

file_out.close()

print("Done with predictions")
