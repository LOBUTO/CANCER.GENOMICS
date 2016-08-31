# cgp_new_sel_pred_tcga_post.py

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

def pear (x, y):

    #Calculate the pearson correlation between 2 vectors
    pear = pear = ( T.sum(x*y) - (T.sum(x)*T.sum(y))/x.shape[0] )  / (T.sqrt(  ( T.sum(T.sqr(x)) - (T.sqr(T.sum(x)))/x.shape[0] ) * ( T.sum(T.sqr(y)) - (T.sqr(T.sum(y)))/y.shape[0]  ) ))
    return pear

class LinearRegression(object):

    def __init__(self, input, n_in, n_out, rng, W=None, b=None):
        """ Initialize the parameters of the logistic regression

        :type input: theano.tensor.TensorType
        :param input: symbolic variable that describes the input of the
                      architecture (one minibatch)

        :type n_in: int
        :param n_in: number of input units, the dimension of the space in
                     which the datapoints lie

        :type n_out: int
        :param n_out: number of output units, the dimension of the space in
                      which the labels lie

        """
        # start-snippet-1
        # initialize with 0 the weights W as a matrix of shape (n_in, n_out)
        #rng = np.random.RandomState(23455)
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


        self.p_y_given_x = T.dot(input, self.W) + self.b

        # symbolic description of how to compute prediction as class whose
        # probability is maximal
        #self.y_pred = T.argmax(self.p_y_given_x, axis=1)
        self.y_pred = self.p_y_given_x[:,0]
        # end-snippet-1

        # parameters of the model
        self.params = [self.W, self.b]

        #self.batch_size = batch_size #Only for purposes of calculating error

    def pred(self, y):
        """Returns prediction only
        """
        return(self.y_pred)

    def errors(self, y):
        """Return a float representing the number of errors in the minibatch
        over the total number of examples of the minibatch ; zero one
        loss over the size of the minibatch

        :type y: theano.tensor.TensorType
        :param y: corresponds to a vector that gives for each example the
                  correct label
        """

        # check if y has same dimension of y_pred
        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                'y should have the same shape as self.y_pred',
                ('y', y.type, 'y_pred', self.y_pred.type)
            )
        # check if y is of the correct datatype
        if y.dtype.startswith('flo'): #CHANGED!!!!!
            # the T.neq operator returns a vector of 0s and 1s, where 1
            # represents a mistake in prediction

            return T.sqrt(T.mean(T.sqr(y-self.y_pred))) #/ (T.max(y) - T.min(y)) #NRMSE
            #return (1/ (2. * batch_size ) ) * T.sum(T.sqr(y-self.y_pred))

        else:
            raise NotImplementedError()

    def loss(self, y):
        """Return a float representing the number of errors in the minibatch
        over the total number of examples of the minibatch ; zero one
        loss over the size of the minibatch

        :type y: theano.tensor.TensorType
        :param y: corresponds to a vector that gives for each example the
                  correct label
        """

        # check if y has same dimension of y_pred
        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                'y should have the same shape as self.y_pred',
                ('y', y.type, 'y_pred', self.y_pred.type)
            )
        # check if y is of the correct datatype
        if y.dtype.startswith('flo'): #CHANGED!!!!!
            # the T.neq operator returns a vector of 0s and 1s, where 1
            # represents a mistake in prediction

            #return T.sqrt(T.mean(T.sqr(y-self.y_pred))) / (T.max(y) - T.min(y)) #NRMSE
            return pear(y, self.y_pred)
        else:
            raise NotImplementedError()

    def NRMSE(self, y):
        """Return a float representing the number of errors in the minibatch
        over the total number of examples of the minibatch ; zero one
        loss over the size of the minibatch

        :type y: theano.tensor.TensorType
        :param y: corresponds to a vector that gives for each example the
                  correct label
        """

        # check if y has same dimension of y_pred
        if y.ndim != self.y_pred.ndim:
            raise TypeError(
                'y should have the same shape as self.y_pred',
                ('y', y.type, 'y_pred', self.y_pred.type)
            )
        # check if y is of the correct datatype
        if y.dtype.startswith('flo'): #CHANGED!!!!!
            # the T.neq operator returns a vector of 0s and 1s, where 1
            # represents a mistake in prediction

            return T.sqrt(T.mean(T.sqr(y-self.y_pred))) / (T.max(y) - T.min(y)) #NRMSE

        else:
            raise NotImplementedError()


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

        self.linearRegressionLayer = LinearRegression(
            input=getattr(self, "layer_" + str(layer_number-1)).output,
            n_in=n_hidden[layer_number-1],
            n_out=n_out,
            rng=rng #,batch_size=batch_size
        )
        self.params = self.params + self.linearRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.layer_0.W).sum() + abs(self.linearRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.layer_0.W ** 2).sum() + (self.linearRegressionLayer.W ** 2).sum()
        )

        # self.negative_log_likelihood = (
        #     self.logRegressionLayer.negative_log_likelihood
        # )
        #
        # self.errors = self.logRegressionLayer.errors
        # self.pred = self.logRegressionLayer.pred
        self.param_to_scale = param_to_scale
        # self.diff = self.logRegressionLayer.diff
        self.errors = self.linearRegressionLayer.errors
        self.loss = self.linearRegressionLayer.loss
        self.NRMSE = self.linearRegressionLayer.NRMSE
        self.pred = self.linearRegressionLayer.pred

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

    data_y = list(drug_data.NORM_pIC50)

    shared_x = theano.shared(np.asarray(data_x, dtype=theano.config.floatX), borrow=True)
    shared_y = theano.shared(np.asarray(data_y, dtype=theano.config.floatX), borrow=True )

    if integers==True:
        return shared_x, T.cast(shared_y, 'int32')
    else:
        return shared_x, shared_y

def model_prediction(MODEL_FILE, test_drug_x):

    model_obj = []

    for l in cPickle.load(open(MODEL_FILE, "rb")):
        model_obj = model_obj + [l]

    #Perform prediction
    n_layers = len(model_obj)-2 #minus loglayer and rng

    input = test_drug_x
    for l in xrange(1, n_layers+1):

        input = prelu(
                        T.dot(input, model_obj[l].W) + model_obj[l].b ,
                        model_obj[l].alpha
                      )

    p_y_given_x = T.dot(input, model_obj[0].W) + model_obj[0].b

    prediction = p_y_given_x[:,0]

    return prediction.eval()

####################################################################################################################################################################################################
# Load data
drug_target = sys.argv[1]
drug_target = (" ").join(drug_target.split("_"))
cancer      = sys.argv[2]
modifier    = sys.argv[3]
extra       = sys.argv[4]

model_folder   = "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_RESULTS/"
model_file     = model_folder + "tcga_" + cancer + "." + modifier + "_new_cgp_sel_model_" + extra + "_" + drug_target + ".pkl"

target_folder  = "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TRAIN_TABLES/"
target_file    = target_folder + cancer + "_all_cgp_new_" + drug_target

out_folder     = "/home/zamalloa/Documents/FOLDER/TCGA_FILES/TCGA_NEW_RESULTS/"
file_out       = out_folder + "cgp_new_modeling_tcga_" + cancer + "_" + drug_target
file_out_2     = out_folder + "cgp_new_modeling_tcga_" + drug_target

####################################################################################################################################################################################################
# Execute for tcga

target_table = pd.read_csv(target_file, sep="\t")
data_x       = target_table.iloc[:,3:]
test_drug_x  = theano.shared(np.asarray(data_x, dtype=theano.config.floatX), borrow=True)

prediction   = model_prediction(model_file, test_drug_x)
actual       = list(target_table.response)

with open(file_out, "w") as f:
    f.write("Cancer" + "\t" + "Compound" + "\t" + "sample" + "\t" + "Actual" + "\t" + "Predicted")

for n in xrange(len(actual)):
    with open(file_out, "a") as dd:
        dd.write("\n" + cancer + "\t" + drug_target + "\t" + target_table.sample[n] + "\t" + str(actual[n]) + "\t" + str(prediction[n]) )

####################################################################################################################################################################################################
# Execute for self
cgp_folder = "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

test_table = pd.read_csv(cgp_folder + "tcga_" + cancer + "." + modifier +  "_TEST_CGP_SEL."  +  drug_target, sep="\t")

test_drug_x, test_drug_y   = shared_drug_dataset_IC50(test_table,  integers=False)

prediction   = model_prediction(model_file, test_drug_x)
actual       = test_drug_y.get_value()

with open(file_out_2, "w") as f:
    f.write("Cancer" + "\t" + "Compound" + "\t" + "cell_name" + "\t" + "Actual" + "\t" + "Predicted")

for n in xrange(len(actual)):
    with open(file_out_2, "a") as dd:
        dd.write("\n" + cancer + "\t" + drug_target + "\t" + test_table.cell_name[n] + "\t" + str(actual[n]) + "\t" + str(prediction[n]) )

####################################################################################################################################################################################################
# Done
print("Done predicting")
