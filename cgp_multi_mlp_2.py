# cgp_multi_mlp.py

import pandas as pd
import os
import sys
import timeit
import numpy as np
import theano
import cPickle
import itertools
import theano.tensor as T
from theano.tensor.nnet import relu
from theano.sandbox.rng_mrg import MRG_RandomStreams as RandomStreams
#from theano.tensor.shared_randomstreams import RandomStreams
#from theano.sandbox.cuda.rng_curand import CURAND_RandomStreams as RandomStreams

def pear (x, y):

    #Calculate the pearson correlation between 2 vectors
    #pear = ( T.sum(x*y) - (T.sum(x)*T.sum(y))/x.shape[0] )  / (T.sqrt(  ( T.sum(T.sqr(x)) - (T.sqr(T.sum(x)))/x.shape[0] ) * ( T.sum(T.sqr(y)) - (T.sqr(T.sum(y)))/y.shape[0]  ) ))
    pear  = T.sum(((x - T.mean(x)) / T.std(x)) * ((y - T.mean(y)) / T.std(y)))  / x.shape[0]
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
        # if W is None:
        #     W_values = np.asarray(
        #             rng.uniform(
        #                 low=-np.sqrt(6. / (n_in + n_out)),
        #                 high=np.sqrt(6. / (n_in + n_out)),
        #                 size=(n_in, n_out)
        #             ),
        #             dtype=theano.config.floatX
        #         )
        #     W = theano.shared(value=W_values, name='W', borrow=True)
        #
        # if b is None :
        #     b = theano.shared(
        #         value=np.zeros(
        #             (n_out,),
        #             dtype=theano.config.floatX
        #         ),
        #         name='b',
        #         borrow=True
        #     )
        #
        # self.W = W
        # self.b = b
        self.W = theano.shared(
            value=np.zeros(
                (n_in, n_out),
                dtype=theano.config.floatX
            ),
            name='W',
            borrow=True
        )
        # initialize the biases b as a vector of n_out 0s
        self.b = theano.shared(
            value=np.zeros(
                (n_out,),
                dtype=theano.config.floatX
            ),
            name='b',
            borrow=True
        )

        self.p_y_given_x = T.dot(input, self.W) + self.b

        self.y_pred = self.p_y_given_x[:,0]

        # parameters of the model
        self.params = [self.W, self.b]

        #self.batch_size = batch_size #Only for purposes of calculating error

    def pred(self, y):
        """Returns prediction only
        """
        return(self.y_pred)

    def errors(self, y):

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

    # def weighted_errors(self, y, y_weights):
    #
    #     # check if y has same dimension of y_pred
    #     if y.ndim != self.y_pred.ndim:
    #         raise TypeError(
    #             'y should have the same shape as self.y_pred',
    #             ('y', y.type, 'y_pred', self.y_pred.type)
    #         )
    #     # check if y is of the correct datatype
    #     if y.dtype.startswith('flo'): #CHANGED!!!!!
    #
    #         return T.sqrt(T.sum(T.sqr(y - self.y_pred) * y_weights) / T.sum(y_weights))
    #
    #     else:
    #         raise NotImplementedError()

    def pear_loss(self, y):
        #return -pear(y, self.y_pred)
        return -T.sum(((y - T.mean(y)) / T.std(y)) * ((self.y_pred - T.mean(self.y_pred)) / T.std(self.y_pred)))  / y.shape[0]

    def pear_check(self, y):

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
            return T.sum(((y - T.mean(y)) / T.std(y)) * ((self.y_pred - T.mean(self.y_pred)) / T.std(self.y_pred)))  / y.shape[0]
        else:
            raise NotImplementedError()

    def NRMSE(self, y):

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

class LogisticRegression(object):

    def __init__(self, input, n_in, n_out):

        # initialize with 0 the weights W as a matrix of shape (n_in, n_out)
        self.W = theano.shared(
            value=np.zeros(
                (n_in, n_out),
                dtype=theano.config.floatX
            ),
            name='W',
            borrow=True
        )
        # initialize the biases b as a vector of n_out 0s
        self.b = theano.shared(
            value=np.zeros(
                (n_out,),
                dtype=theano.config.floatX
            ),
            name='b',
            borrow=True
        )

        self.p_y_given_x = T.nnet.softmax(T.dot(input, self.W) + self.b)

        self.y_pred = T.argmax(self.p_y_given_x, axis=1)

        self.params = [self.W, self.b]

        self.input = input

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
        if y.dtype.startswith('int'):
            # the T.neq operator returns a vector of 0s and 1s, where 1
            # represents a mistake in prediction
            return T.mean(T.neq(self.y_pred, y))
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
            alpha_value = np.full((n_out), .3,  dtype=theano.config.floatX)
            alpha = theano.shared(value=alpha_value, name='alpha', borrow=True)

        self.W = W
        self.b = b
        self.alpha = alpha

        lin_output = T.dot(input, self.W) + self.b
        # if (activation == T.tanh) or (activation == relu):
        #     output = activation(lin_output)
        #     self.params = [self.W, self.b]
        #
        # elif activation == prelu:
        #     output = activation(lin_output, self.alpha)
        #     self.params = [self.W, self.b, self.alpha]
        # else:
        #     output = lin_output

        output = activation(lin_output, self.alpha)
        self.params = [self.W, self.b, self.alpha]

        # output = (
        #     lin_output if activation is None
        #     else activation(lin_output, self.alpha)
        # )

        #Is droput necessary?
        if dropout==True:
            train_output = drop(output, rng=rng, p=p)
            self.output = T.switch(T.neq(is_train, 0), train_output, output)
        else:
            self.output = output

        # parameters of the model
        # if activation == prelu:
        #     self.params = [self.W, self.b, self.alpha]
        # else:
        #     self.params = [self.W, self.b]

class Multiplicative_fusion(object):
    # Performs combinations of two neural layers (drug_input, cell_input) from
    # different sources and output the Multiplicative_fusion layer along with
    # 4 parameters to be learned.

    def __init__(self, rng, drug_input, cell_input, drug_in, cell_in, neural_range,
                 cell_alpha=None, drug_alpha=None, cell_beta=None, drug_beta=None):

        self.drug_input = drug_input
        self.cell_input = cell_input

        if cell_alpha is None:

            cell_alpha_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(cell_in)),
                                                high =  np.sqrt(6. /(cell_in)),
                                                size = (cell_in)
                                          ), dtype=theano.config.floatX)
            cell_alpha = theano.shared(value=cell_alpha_value, name='cell_alpha', borrow=True)

        if drug_alpha is None:

            #drug_alpha_value = np.full((drug_in), .1,  dtype=theano.config.floatX)
            drug_alpha_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(drug_in)),
                                                high =  np.sqrt(6. /(drug_in)),
                                                size = (drug_in)
                                          ), dtype=theano.config.floatX)
            drug_alpha = theano.shared(value=drug_alpha_value, name='drug_alpha', borrow=True)

        if cell_beta is None:

            #cell_beta_value = np.full((cell_in), .1,  dtype=theano.config.floatX)
            cell_beta_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(cell_in)),
                                                high =  np.sqrt(6. /(cell_in)),
                                                size = (cell_in)
                                          ), dtype=theano.config.floatX)
            cell_beta = theano.shared(value=cell_beta_value, name='cell_beta', borrow=True)

        if drug_beta is None:

            #drug_beta_value = np.full((drug_in), .1,  dtype=theano.config.floatX)
            drug_beta_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(drug_in)),
                                                high =  np.sqrt(6. /(drug_in)),
                                                size = (drug_in)
                                          ), dtype=theano.config.floatX)
            drug_beta = theano.shared(value=drug_beta_value, name='drug_beta', borrow=True)

        self.cell_alpha = cell_alpha
        self.drug_alpha = drug_alpha
        self.cell_beta  = cell_beta
        self.drug_beta  = drug_beta

        # Apply linear modifications to both input based on parameters
        drug_output = (drug_input * self.drug_alpha) + self.drug_beta
        cell_output = (cell_input * self.cell_alpha) + self.cell_beta

        # Apply pairwise neuron multiplication
        # output      = T.concatenate(
        #                               [drug_output[:, [i]] * cell_output for i in neural_range],
        #                                axis = 1
        #                             )
        # output      = drug_output * cell_output # Keep in mind that both have to have the same number of neurons
        output       = T.concatenate([drug_output, cell_output], axis = 1) #NOTE: Modified to "true concatenation"
        # Output
        self.output = output

        # Parameters of the fusion
        self.params = [self.cell_alpha, self.drug_alpha, self.cell_beta, self.drug_beta]

class Multiplicative_fusion_zero_drug(object):
    # Performs synthetic combination cell input layer
    # and output the Multiplicative_fusion layer along with
    # 2 parameters to be learned.

    def __init__(self, rng, cell_input, cell_in,
                 cell_alpha=None, cell_beta=None):

        self.cell_input = cell_input

        if cell_alpha is None:

            cell_alpha_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(cell_in)),
                                                high =  np.sqrt(6. /(cell_in)),
                                                size = (cell_in)
                                          ), dtype=theano.config.floatX)
            cell_alpha = theano.shared(value=cell_alpha_value, name='cell_alpha', borrow=True)

        if cell_beta is None:

            #cell_beta_value = np.full((cell_in), .1,  dtype=theano.config.floatX)
            cell_beta_value = np.asarray(rng.uniform(
                                                low  = -np.sqrt(6. /(cell_in)),
                                                high =  np.sqrt(6. /(cell_in)),
                                                size = (cell_in)
                                          ), dtype=theano.config.floatX)
            cell_beta = theano.shared(value=cell_beta_value, name='cell_beta', borrow=True)

        self.cell_alpha = cell_alpha
        self.cell_beta  = cell_beta

        # Apply linear modifications to both input based on parameters
        cell_output = (cell_input * self.cell_alpha) + self.cell_beta

        output      = cell_output
        # Output
        self.output = output

        # Parameters of the fusion
        self.params = [self.cell_alpha, self.cell_beta]

def rescale_weights(params, incoming_max):
    incoming_max = np.cast[theano.config.floatX](incoming_max)
    for p in params:
        w = p.get_value()
        w_sum = (w**2).sum(axis=0)
        w[:, w_sum>incoming_max] = w[:, w_sum>incoming_max] * np.sqrt(incoming_max) / w_sum[w_sum>incoming_max]
        p.set_value(w)

class Multi_MLP_Regression(object):
    def __init__(self, rng, cell_input, drug_input, is_train,
                 cell_n_in, drug_n_in, cell_n_hidden, drug_n_hidden, fusion_n_hidden, neural_range,
                 n_out,
                 cell_p=0.5, cell_dropout=False, cell_input_p=0.1,
                 drug_p=0.5, drug_dropout=False, drug_input_p=0.1):

        # PROCESS CELL INPUT FIRST
        if cell_input_p!=None:
            self.cell_input_layer = drop(cell_input, rng=rng, p=cell_input_p)
            self.cell_input_layer = T.switch(T.neq(is_train, 0), self.cell_input_layer, cell_input)
        else:
            self.cell_input_layer = cell_input

        param_to_scale = [] #To scale weights to square length of 15

        self.cell_layer_0 = HiddenLayer(
            rng=rng,
            input=self.cell_input_layer,
            n_in=cell_n_in,
            n_out=cell_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=cell_p,
            dropout=cell_dropout
        )

        self.params = self.cell_layer_0.params
        param_to_scale = param_to_scale + [self.cell_layer_0.params[0]]

        cell_layer_number = 1
        if len(cell_n_hidden)>1:

            for layer in cell_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
                                                    n_in=cell_n_hidden[cell_layer_number-1],
                                                    n_out=cell_n_hidden[cell_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=cell_p,
                                                    dropout=cell_dropout
                                                )

                setattr(self, "cell_layer_" + str(cell_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "cell_layer_" + str(cell_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "cell_layer_" + str(cell_layer_number)).params[0]]

                cell_layer_number = cell_layer_number + 1

        # PROCESS DRUG INPUT NEXT
        if drug_input_p!=None:
            self.drug_input_layer = drop(drug_input, rng=rng, p=drug_input_p)
            self.drug_input_layer = T.switch(T.neq(is_train, 0), self.drug_input_layer, drug_input)
        else:
            self.drug_input_layer = drug_input

        self.drug_layer_0 = HiddenLayer(
            rng=rng,
            input=self.drug_input_layer,
            n_in=drug_n_in,
            n_out=drug_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=drug_p,
            dropout=drug_dropout
        )

        self.params = self.params + self.drug_layer_0.params # Adding to previous cell params
        param_to_scale = param_to_scale + [self.drug_layer_0.params[0]]

        drug_layer_number = 1
        if len(drug_n_hidden)>1:

            for layer in drug_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "drug_layer_" + str(drug_layer_number-1)).output,
                                                    n_in=drug_n_hidden[drug_layer_number-1],
                                                    n_out=drug_n_hidden[drug_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=drug_p,
                                                    dropout=drug_dropout
                                                )

                setattr(self, "drug_layer_" + str(drug_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "drug_layer_" + str(drug_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "drug_layer_" + str(drug_layer_number)).params[0]]

                drug_layer_number = drug_layer_number + 1

        # APPLY FUSION
        # Combine both outputs to obtain Multiplicative fusion layer
        self.multiplicative_input = Multiplicative_fusion(
            drug_input = getattr(self, "drug_layer_" + str(drug_layer_number-1)).output,
            cell_input = getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
            drug_in  = drug_n_hidden[drug_layer_number-1],
            cell_in  = cell_n_hidden[cell_layer_number-1],
            neural_range = neural_range, #NOTE: Integer representing total number of summed features between drug and cell previous layer
            rng =  rng
        )
        self.params = self.params + self.multiplicative_input.params

        # PROCESS FUSION LAYERS
        self.fusion_layer_0 = HiddenLayer(
            rng=rng,
            input=self.multiplicative_input.output,
            n_in=neural_range, #NOTE: Integer representing total number of summed features between drug and cell previous layer
            n_out=fusion_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=cell_p,
            dropout=cell_dropout
        )

        self.params = self.params + self.fusion_layer_0.params # Plus previous separate layer params (drug + cell + mf)
        param_to_scale = param_to_scale + [self.fusion_layer_0.params[0]]

        fusion_layer_number = 1
        if len(fusion_n_hidden)>1:

            for layer in fusion_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
                                                    n_in=fusion_n_hidden[fusion_layer_number-1],
                                                    n_out=fusion_n_hidden[fusion_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=cell_p, #MAY CHANGE TO OWN VARIABLE
                                                    dropout=cell_dropout #MAY CHANGE TO OWN VARIABLE
                                                )

                setattr(self, "fusion_layer_" + str(fusion_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "fusion_layer_" + str(fusion_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "fusion_layer_" + str(fusion_layer_number)).params[0]]

                fusion_layer_number = fusion_layer_number + 1


        # APPLY LINEAR REGRESSION
        self.linearRegressionLayer = LinearRegression(
            input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
            n_in=fusion_n_hidden[fusion_layer_number-1],
            n_out=n_out,
            rng = rng
        )

        self.params = self.params + self.linearRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.drug_layer_0.W).sum() + abs(self.linearRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.drug_layer_0.W ** 2).sum() + (self.linearRegressionLayer.W ** 2).sum()
        )

        self.errors = self.linearRegressionLayer.errors
        self.pear_loss = self.linearRegressionLayer.pear_loss
        self.pear_check = self.linearRegressionLayer.pear_check
        self.NRMSE = self.linearRegressionLayer.NRMSE
        self.pred = self.linearRegressionLayer.pred

        self.param_to_scale = param_to_scale

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

class Multi_MLP_Regression_Zero_Drug(object):
    def __init__(self, rng, cell_input, is_train,
                 cell_n_in, cell_n_hidden,
                 n_out, p=0.5, dropout=False, input_p=0.1):

        # PROCESS CELL INPUT FIRST
        if input_p!=None:
            self.cell_input_layer = drop(cell_input, rng=rng, p=input_p)
            self.cell_input_layer = T.switch(T.neq(is_train, 0), self.cell_input_layer, cell_input)
        else:
            self.cell_input_layer = cell_input

        param_to_scale = [] #To scale weights to square length of 15

        self.cell_layer_0 = HiddenLayer(
            rng=rng,
            input=self.cell_input_layer,
            n_in=cell_n_in,
            n_out=cell_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=p,
            dropout=dropout
        )

        self.params = self.cell_layer_0.params
        param_to_scale = param_to_scale + [self.cell_layer_0.params[0]]

        cell_layer_number = 1
        if len(cell_n_hidden)>1:

            for layer in cell_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
                                                    n_in=cell_n_hidden[cell_layer_number-1],
                                                    n_out=cell_n_hidden[cell_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
                                                )

                setattr(self, "cell_layer_" + str(cell_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "cell_layer_" + str(cell_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "cell_layer_" + str(cell_layer_number)).params[0]]

                cell_layer_number = cell_layer_number + 1

        # NO FUSION #
        # APPLY LINEAR REGRESSION
        self.linearRegressionLayer = LinearRegression(
            input=getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
            n_in=cell_n_hidden[cell_layer_number-1],
            n_out=n_out,
            rng = rng
        )

        self.params = self.params + self.linearRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.cell_layer_0.W).sum() + abs(self.linearRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.cell_layer_0.W ** 2).sum() + (self.linearRegressionLayer.W ** 2).sum()
        )

        self.errors = self.linearRegressionLayer.errors
        self.pear_loss = self.linearRegressionLayer.pear_loss
        self.pear_check = self.linearRegressionLayer.pear_check
        self.NRMSE = self.linearRegressionLayer.NRMSE
        self.pred = self.linearRegressionLayer.pred

        self.param_to_scale = param_to_scale

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

class Multi_MLP_Class(object):
    def __init__(self, rng, cell_input, drug_input, is_train,
                 cell_n_in, drug_n_in, cell_n_hidden, drug_n_hidden, fusion_n_hidden, neural_range,
                 n_out,
                 cell_p=0.5, cell_dropout=False, cell_input_p=0.1,
                 drug_p=0.5, drug_dropout=False, drug_input_p=0.1):

        # PROCESS CELL INPUT FIRST
        if cell_input_p!=None:
            self.cell_input_layer = drop(cell_input, rng=rng, p=cell_input_p)
            self.cell_input_layer = T.switch(T.neq(is_train, 0), self.cell_input_layer, cell_input)
        else:
            self.cell_input_layer = cell_input

        param_to_scale = [] #To scale weights to square length of 15

        self.cell_layer_0 = HiddenLayer(
            rng=rng,
            input=self.cell_input_layer,
            n_in=cell_n_in,
            n_out=cell_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=cell_p,
            dropout=cell_dropout
        )

        self.params = self.cell_layer_0.params
        param_to_scale = param_to_scale + [self.cell_layer_0.params[0]]

        cell_layer_number = 1
        if len(cell_n_hidden)>1:

            for layer in cell_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
                                                    n_in=cell_n_hidden[cell_layer_number-1],
                                                    n_out=cell_n_hidden[cell_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=cell_p,
                                                    dropout=cell_dropout
                                                )

                setattr(self, "cell_layer_" + str(cell_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "cell_layer_" + str(cell_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "cell_layer_" + str(cell_layer_number)).params[0]]

                cell_layer_number = cell_layer_number + 1

        # PROCESS DRUG INPUT NEXT
        if drug_input_p!=None:
            self.drug_input_layer = drop(drug_input, rng=rng, p=drug_input_p)
            self.drug_input_layer = T.switch(T.neq(is_train, 0), self.drug_input_layer, drug_input)
        else:
            self.drug_input_layer = drug_input

        self.drug_layer_0 = HiddenLayer(
            rng=rng,
            input=self.drug_input_layer,
            n_in=drug_n_in,
            n_out=drug_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=drug_p,
            dropout=drug_dropout
        )

        self.params = self.params + self.drug_layer_0.params # Adding to previous cell params
        param_to_scale = param_to_scale + [self.drug_layer_0.params[0]]

        drug_layer_number = 1
        if len(drug_n_hidden)>1:

            for layer in drug_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "drug_layer_" + str(drug_layer_number-1)).output,
                                                    n_in=drug_n_hidden[drug_layer_number-1],
                                                    n_out=drug_n_hidden[drug_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=drug_p,
                                                    dropout=drug_dropout
                                                )

                setattr(self, "drug_layer_" + str(drug_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "drug_layer_" + str(drug_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "drug_layer_" + str(drug_layer_number)).params[0]]

                drug_layer_number = drug_layer_number + 1

        # APPLY FUSION
        # Combine both outputs to obtain Multiplicative fusion layer
        self.multiplicative_input = Multiplicative_fusion(
            drug_input = getattr(self, "drug_layer_" + str(drug_layer_number-1)).output,
            cell_input = getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
            drug_in  = drug_n_hidden[drug_layer_number-1],
            cell_in  = cell_n_hidden[cell_layer_number-1],
            neural_range = neural_range,
            rng = rng
        )
        self.params = self.params + self.multiplicative_input.params

        # PROCESS FUSION LAYERS
        self.fusion_layer_0 = HiddenLayer(
            rng=rng,
            input=self.multiplicative_input.output,
            n_in=neural_range,#drug_n_hidden[drug_layer_number-1] * cell_n_hidden[cell_layer_number-1],
            n_out=fusion_n_hidden[0],
            activation=relu,
            is_train=is_train,
            p=cell_p, #MAY CHANGE TO OWN VARIABLE
            dropout=cell_dropout #MAY CHANGE TO OWN VARIABLE
        )

        self.params = self.params + self.fusion_layer_0.params # Plus previous separate layer params (drug + cell + mf)
        param_to_scale = param_to_scale + [self.fusion_layer_0.params[0]]

        fusion_layer_number = 1
        if len(fusion_n_hidden)>1:

            for layer in fusion_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
                                                    n_in=fusion_n_hidden[fusion_layer_number-1],
                                                    n_out=fusion_n_hidden[fusion_layer_number],
                                                    activation=relu,
                                                    is_train=is_train,
                                                    p=cell_p, #MAY CHANGE TO OWN VARIABLE
                                                    dropout=cell_dropout #MAY CHANGE TO OWN VARIABLE
                                                )

                setattr(self, "fusion_layer_" + str(fusion_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "fusion_layer_" + str(fusion_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "fusion_layer_" + str(fusion_layer_number)).params[0]]

                fusion_layer_number = fusion_layer_number + 1


        # APPLY LOGISTIC REGRESSION
        # The logistic regression layer gets as input the fused multiplicate_input
        self.logRegressionLayer = LogisticRegression(
            input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
            n_in=fusion_n_hidden[fusion_layer_number-1],
            n_out=n_out
        )

        self.params = self.params + self.logRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.drug_layer_0.W).sum() + abs(self.logRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.drug_layer_0.W ** 2).sum() + (self.logRegressionLayer.W ** 2).sum()
        )

        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )

        self.errors = self.logRegressionLayer.errors
        self.pred = self.logRegressionLayer.pred
        self.param_to_scale = param_to_scale

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

class Multi_MLP_Class_Zero_Drug(object):
    def __init__(self, rng, cell_input, is_train,
                 cell_n_in, cell_n_hidden, fusion_n_hidden,
                 n_out, p=0.5, dropout=False, input_p=0.1):

        # PROCESS CELL INPUT FIRST
        if input_p!=None:
            self.cell_input_layer = drop(cell_input, rng=rng, p=input_p)
            self.cell_input_layer = T.switch(T.neq(is_train, 0), self.cell_input_layer, cell_input)
        else:
            self.cell_input_layer = cell_input

        param_to_scale = [] #To scale weights to square length of 15

        self.cell_layer_0 = HiddenLayer(
            rng=rng,
            input=self.cell_input_layer,
            n_in=cell_n_in,
            n_out=cell_n_hidden[0],
            activation=prelu,
            is_train=is_train,
            p=p,
            dropout=dropout
        )

        self.params = self.cell_layer_0.params
        param_to_scale = param_to_scale + [self.cell_layer_0.params[0]]

        cell_layer_number = 1
        if len(cell_n_hidden)>1:

            for layer in cell_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
                                                    n_in=cell_n_hidden[cell_layer_number-1],
                                                    n_out=cell_n_hidden[cell_layer_number],
                                                    activation=prelu,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
                                                )

                setattr(self, "cell_layer_" + str(cell_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "cell_layer_" + str(cell_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "cell_layer_" + str(cell_layer_number)).params[0]]

                cell_layer_number = cell_layer_number + 1


        # APPLY FUSION
        # Combine single output to obtain Multiplicative fusion layer
        self.multiplicative_input = Multiplicative_fusion_zero_drug(
            cell_input = getattr(self, "cell_layer_" + str(cell_layer_number-1)).output,
            cell_in  = cell_n_hidden[cell_layer_number-1],
            rng = rng
        )
        self.params = self.params + self.multiplicative_input.params

        # PROCESS FUSION LAYERS
        self.fusion_layer_0 = HiddenLayer(
            rng=rng,
            input=self.multiplicative_input.output,
            n_in=cell_n_hidden[cell_layer_number-1],
            n_out=fusion_n_hidden[0],
            activation=prelu,
            is_train=is_train,
            p=p,
            dropout=dropout
        )

        self.params = self.params + self.fusion_layer_0.params # Plus previous separate layer params (drug + cell + mf)
        param_to_scale = param_to_scale + [self.fusion_layer_0.params[0]]

        fusion_layer_number = 1
        if len(fusion_n_hidden)>1:

            for layer in fusion_n_hidden[1:]:

                current_hidden_layer = HiddenLayer(
                                                    rng=rng,
                                                    input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
                                                    n_in=fusion_n_hidden[fusion_layer_number-1],
                                                    n_out=fusion_n_hidden[fusion_layer_number],
                                                    activation=prelu,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
                                                )

                setattr(self, "fusion_layer_" + str(fusion_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "fusion_layer_" + str(fusion_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "fusion_layer_" + str(fusion_layer_number)).params[0]]

                fusion_layer_number = fusion_layer_number + 1


        # APPLY LOGISTIC REGRESSION
        # The logistic regression layer gets as input the fused multiplicate_input
        self.logRegressionLayer = LogisticRegression(
            input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
            n_in=fusion_n_hidden[fusion_layer_number-1],
            n_out=n_out,
        )

        self.params = self.params + self.logRegressionLayer.params

        #L1 and L2 regularization
        self.L1 = (
            abs(self.cell_layer_0.W).sum() + abs(self.logRegressionLayer.W).sum()
        )

        self.L2_sqr = (
            (self.cell_layer_0.W ** 2).sum() + (self.logRegressionLayer.W ** 2).sum()
        )

        self.negative_log_likelihood = (
            self.logRegressionLayer.negative_log_likelihood
        )

        self.errors = self.logRegressionLayer.errors
        self.pred = self.logRegressionLayer.pred
        self.param_to_scale = param_to_scale

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!


def regression_mlp_mf(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=20,
             cell_n_hidden=[500,200,100], drug_n_hidden=[500,200,100], mf_manual = "None", fusion_n_hidden = [500,200,100],
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):

    #Demonstrate stochastic gradient descent optimization for a multilayer
    #perceptron for parallel drug and cell layers to Multiplicative fusion layer
    train_drug_x, train_cell_x, train_drug_index_x, train_cell_index_x, train_set_y = datasets[0]
    valid_drug_x, valid_cell_x, valid_drug_index_x, valid_cell_index_x, valid_set_y = datasets[1]
    test_drug_x,  test_cell_x,  test_drug_index_x,  test_cell_index_x,  test_set_y = datasets[2]

    # valid_batch_size = valid_drug_index_x.eval().shape[0] # Could have been valid_cell_index_x since they are of equal length
    valid_batch_size = train_batch_size #MODIFIED
    test_batch_size  = train_batch_size #MODIFIED
    train_samples    = train_drug_index_x.eval().shape[0] # Could have been train_cell_index_x since they are of equal length
    valid_samples    = valid_drug_index_x.eval().shape[0] #MODIFIED
    test_samples     = test_drug_index_x.eval().shape[0] #MODIFIED

    # Compute input layer for both drug and cell networks
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]
    DRUG_N_IN = train_drug_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size
    n_valid_batches = valid_samples / valid_batch_size #MODIFIED
    n_test_batches  = test_samples / test_batch_size

    # compute fusion neural range
    # if mf_manual=="None":
    #
    #     if cell_n_hidden[-1] != drug_n_hidden[-1]:
    #
    #         neural_range  = min([cell_n_hidden[-1], drug_n_hidden[-1]])
    #         if neural_range == cell_n_hidden[-1]:
    #             drug_n_hidden = drug_n_hidden + [neural_range] # In place modification!!!
    #         else :
    #             cell_n_hidden = cell_n_hidden + [neural_range] # In place modification!!!
    #     else:
    #         neural_range = cell_n_hidden[-1]
    # else:
    #     neural_range  = mf_manual
    #     drug_n_hidden = drug_n_hidden + [mf_manual]
    #     cell_n_hidden = cell_n_hidden + [mf_manual]
    neural_range = cell_n_hidden[-1] + drug_n_hidden[-1] #NOTE: Addition of last two layers

    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')
    x_d = T.matrix('x_d')

    y = T.fvector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Regression(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        drug_input = x_d, #needed
        cell_n_in = CELL_N_IN, #calculated
        drug_n_in = DRUG_N_IN, #calculated
        cell_n_hidden = cell_n_hidden, #calculated
        drug_n_hidden = drug_n_hidden, #calculated
        fusion_n_hidden = fusion_n_hidden,
        neural_range = neural_range, #calculated
        n_out=1,
        cell_p=p,
        cell_dropout=dropout,
        cell_input_p=input_p,
        drug_p=0.7,
        drug_dropout=True,
        drug_input_p=0.2
    )

    cost = (
        classifier.NRMSE(y) #pear_loss, NRMSE is more accurate if we are randomly sampling every mini-batch
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )

    validate_model = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            # x_c: valid_cell_x[valid_cell_index_x,],
            # x_d: valid_drug_x[valid_drug_index_x,],
            # y: valid_set_y,
            x_c: valid_cell_x[valid_cell_index_x[index * valid_batch_size:(index + 1) * valid_batch_size],],
            x_d: valid_drug_x[valid_drug_index_x[index * valid_batch_size:(index + 1) * valid_batch_size],],
            y: valid_set_y[index * valid_batch_size:(index + 1) * valid_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    #NOTE: Temporarily modified
    # test_cor = theano.function(
    #     inputs=[index],
    #     outputs=classifier.pear_check(y),
    #     givens={
    #         x_c: test_cell_x[test_cell_index_x[index * test_batch_size:(index + 1) * test_batch_size],],
    #         x_d: test_drug_x[test_drug_index_x[index * test_batch_size:(index + 1) * test_batch_size],],
    #         y: test_set_y[index * test_batch_size:(index + 1) * test_batch_size],
    #         is_train: np.cast['int32'](0)
    #     },
    #     on_unused_input='warn',
    # )

    test_nrmse = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: test_cell_x[test_cell_index_x[index * test_batch_size:(index + 1) * test_batch_size],],
            x_d: test_drug_x[test_drug_index_x[index * test_batch_size:(index + 1) * test_batch_size],],
            y: test_set_y[index * test_batch_size:(index + 1) * test_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    #NOTE: Temporarily modified
    # test_pred = theano.function(
    #     inputs=[index],
    #     outputs=classifier.pred(y),
    #     givens={
    #         x_c: test_cell_x[test_cell_index_x,],
    #         x_d: test_drug_x[test_drug_index_x,],
    #         y: test_set_y,
    #         is_train: np.cast['int32'](0)
    #     },
    #     on_unused_input='warn',
    # )

    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            x_d: train_drug_x[train_drug_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            x_d: train_drug_x[train_drug_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf #MODIFIED
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    test_pear = 0
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN"+ "\t"+"VALID.ERROR" + "\t" + "TEST.NRMSE")# + "\t" + "TEST.COR")#NOTE:Temporarily modified
    FILE_OUT.close()

    FILE_OUT_val = open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "w")
    FILE_OUT_val.write("EPOCH" +"\t" + "ACTUAL" +"\t"+"PREDICTED")
    FILE_OUT_val.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            # iteration number
            #iter = (epoch - 1) * n_train_batches + minibatch_index

            #if (iter + 1) % validation_frequency == 0:
            if (minibatch_index + 1) % EPOCH_SIZE == 0:
                # compute zero-one loss on validation set

                validation_losses = [validate_model(i) for i in xrange(n_valid_batches)]
                this_validation_loss = np.mean(validation_losses)

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)


                log = ('epoch %i, minibatch %i/%i, train error %f ,validation error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error ,
                        this_validation_loss
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error) + "\t"+ str(this_validation_loss) \
                                   +"\t" +str(test_loss))# + "\t" + str(test_pear)) #NOTE:Temporarily modified

                # if we got the best validation score until now
                if this_validation_loss < best_validation_loss: #MODIFIED
                    LR_COUNT = 0

                    #improve patience if loss improvement is good enough
                    # if (
                    #     this_validation_loss < best_validation_loss *
                    #     improvement_threshold
                    # ):
                    #     patience = max(patience, iter * patience_increase)

                    best_validation_loss = this_validation_loss
                    best_iter = iter

                    # test it on the test set
                    test_losses = [test_nrmse(i) for i in xrange(n_test_batches)]
                    test_loss = np.mean(test_losses)

                    #NOTE: Temporarily modified
                    # test_pears = [test_cor(i) for i in xrange(n_test_batches)]
                    # test_pear = np.mean(test_pears)

                    log = ((' epoch %i, minibatch %i/%i, test error of '
                        'best nrmse and pear %f %%') %
                        (epoch, minibatch_index + 1, EPOCH_SIZE, test_loss))#, test_pear))#NOTE: Temporarily modified
                    # print(log)
                    with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                        logfile.write(log + "\n")

                    #NOTE:ONLY SAVE MODEL if validation improves
                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["drug_n_hidden"]   = [getattr(classifier, "drug_layer_" + str(e)) for e in xrange(len(drug_n_hidden))]
                    MODEL["multiplicative"]  = classifier.multiplicative_input
                    MODEL["fusion_n_hidden"] = [getattr(classifier, "fusion_layer_" + str(e)) for e in xrange(len(fusion_n_hidden))]
                    MODEL["linear"]          = classifier.linearRegressionLayer

                    with open(OUT_FOLDER + "/" + drug_name + "_" + str(epoch) + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                    #NOTE:ONLY write if validation improvement
                    # ACTUAL = test_set_y.get_value()
                    # PREDICTED = [test_pred(i) for i in xrange(n_test_batches)][0]
                    #
                    # with open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "a") as FILE_OUT_val:
                    #     for l in xrange(len(ACTUAL)):
                    #         FILE_OUT_val.write("\n" + str(epoch) + "\t" + str(ACTUAL[l]) + "\t" + str(PREDICTED[l]))
                else:
                    LR_COUNT = LR_COUNT+1

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))

def regression_mlp_mf_zero_drug(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=10,
             cell_n_hidden=[500,200,100], mf_manual = "None",
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):

    #Demonstrate stochastic gradient descent optimization for a multilayer
    #perceptron for parallel drug and cell layers to Multiplicative fusion layer
    train_cell_x, train_cell_index_x, train_set_y = datasets[0]
    valid_cell_x, valid_cell_index_x, valid_set_y = datasets[1]
    test_cell_x,  test_cell_index_x,  test_set_y = datasets[2]

    # valid_batch_size = valid_cell_index_x.eval().shape[0]
    valid_batch_size = train_batch_size #MODIFIED
    test_batch_size  = test_cell_index_x.eval().shape[0]
    train_samples    = train_cell_index_x.eval().shape[0]
    valid_samples    = valid_cell_index_x.eval().shape[0] #MODIFIED

    # Compute input layer
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size
    n_valid_batches = valid_samples / valid_batch_size #MODIFIED
    n_test_batches  = test_batch_size / test_batch_size #1

    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')

    y = T.fvector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Regression_Zero_Drug(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        cell_n_in = CELL_N_IN, #calculated
        cell_n_hidden = cell_n_hidden, #calculated
        n_out=1,
        p=p,
        dropout=dropout,
        input_p=input_p
    )

    cost = (
        classifier.NRMSE(y)
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )

    validate_model = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: valid_cell_x[valid_cell_index_x[index * valid_batch_size:(index + 1) * valid_batch_size],],
            y: valid_set_y[index * valid_batch_size:(index + 1) * valid_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_cor = theano.function(
        inputs=[index],
        outputs=classifier.pear_check(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_nrmse = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_pred = theano.function(
        inputs=[index],
        outputs=classifier.pred(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.pear_check(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    test_pear = 0
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN"+ "\t"+"VALID.ERROR" + "\t" + "TEST.COR" + "\t" + "TEST.ERROR")
    FILE_OUT.close()

    FILE_OUT_val = open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "w")
    FILE_OUT_val.write("EPOCH" +"\t" + "ACTUAL" +"\t"+"PREDICTED")
    FILE_OUT_val.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            # iteration number
            #iter = (epoch - 1) * n_train_batches + minibatch_index

            #if (iter + 1) % validation_frequency == 0:
            if (minibatch_index + 1) % EPOCH_SIZE == 0:
                # compute zero-one loss on validation set

                validation_losses = [validate_model(i) for i in xrange(n_valid_batches)]
                this_validation_loss = np.mean(validation_losses)

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)


                log = ('epoch %i, minibatch %i/%i, train error %f ,validation error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error ,
                        this_validation_loss
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error) + "\t"+ str(this_validation_loss) \
                                   +"\t" +str(test_pear) + "\t" + str(test_loss))

                # if we got the best validation score until now
                if this_validation_loss < best_validation_loss:
                    LR_COUNT = 0

                    #improve patience if loss improvement is good enough
                    # if (
                    #     this_validation_loss < best_validation_loss *
                    #     improvement_threshold
                    # ):
                    #     patience = max(patience, iter * patience_increase)

                    best_validation_loss = this_validation_loss
                    best_iter = iter

                    # test it on the test set
                    test_losses = [test_nrmse(i) for i in xrange(n_test_batches)]
                    test_loss = np.mean(test_losses)

                    test_pears = [test_cor(i) for i in xrange(n_test_batches)]
                    test_pear = np.mean(test_pears)

                    log = ((' epoch %i, minibatch %i/%i, test error of '
                        'best nrmse and pear %f,%f %%') %
                        (epoch, minibatch_index + 1, EPOCH_SIZE, test_loss, test_pear))
                    # print(log)
                    with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                        logfile.write(log + "\n")

                    #ONLY SAVE MODEL if validation improves
                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["linear"]          = classifier.linearRegressionLayer

                    # MODEL = []
                    # for e in xrange(len(cell_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "cell_layer_" + str(e))]
                    # for e in xrange(len(drug_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "drug_layer_" + str(e))]
                    # for e in xrange(len(fusion_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "fusion_layer_" + str(e))]
                    #
                    # MODEL = MODEL + [classifier.logRegressionLayer]
                    with open(OUT_FOLDER + "/" + drug_name + "_" + str(epoch) + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                    #Only write if validation improvement
                    ACTUAL = test_set_y.get_value()
                    PREDICTED = [test_pred(i) for i in xrange(n_test_batches)][0]

                    with open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "a") as FILE_OUT_val:
                        for l in xrange(len(ACTUAL)):
                            FILE_OUT_val.write("\n" + str(epoch) + "\t" + str(ACTUAL[l]) + "\t" + str(PREDICTED[l]))
                else:
                    LR_COUNT = LR_COUNT+1

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))

def regression_mlp_train_mf_zero_drug(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=10,
             cell_n_hidden=[500,200,100], mf_manual = "None",
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):

    train_cell_x, train_cell_index_x, train_set_y = datasets[0]
    train_samples    = train_cell_index_x.eval().shape[0]

    # Compute input layer
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size

    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')

    y = T.fvector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Regression_Zero_Drug(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        cell_n_in = CELL_N_IN, #calculated
        cell_n_hidden = cell_n_hidden,
        n_out=1,
        p=p,
        dropout=dropout,
        input_p=input_p
    )

    cost = (
        classifier.NRMSE(y)
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )
    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    test_pear = 0
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN")
    FILE_OUT.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            # iteration number
            #iter = (epoch - 1) * n_train_batches + minibatch_index

            #if (iter + 1) % validation_frequency == 0:
            if (minibatch_index + 1) % EPOCH_SIZE == 0:

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)


                log = ('epoch %i, minibatch %i/%i, train error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error))

                if (epoch % 100)==0: #Write model every 200 epochs
                    LR_COUNT = 0

                    #ONLY SAVE MODEL if validation improves
                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["linear"]          = classifier.linearRegressionLayer

                    with open(OUT_FOLDER + "/" + drug_name + "_" + str(epoch) + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                else:
                    LR_COUNT = LR_COUNT+1

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))


def class_mlp_mf_zero_drug(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=20,
             cell_n_hidden=[500,200,100], fusion_n_hidden = [500,200,100],
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):

    #Demonstrate stochastic gradient descent optimization for a multilayer
    #perceptron for parallel drug and cell layers to Multiplicative fusion layer
    train_cell_x, train_cell_index_x, train_set_y = datasets[0]
    valid_cell_x, valid_cell_index_x, valid_set_y = datasets[1]
    test_cell_x,  test_cell_index_x,  test_set_y = datasets[2]

    valid_batch_size = valid_cell_index_x.eval().shape[0]
    test_batch_size  = test_cell_index_x.eval().shape[0]
    train_samples    = train_cell_index_x.eval().shape[0]

    # Compute input layer
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size
    n_valid_batches = valid_batch_size / valid_batch_size #1
    n_test_batches  = test_batch_size / test_batch_size #1

    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')

    y = T.ivector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Class_Zero_Drug(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        cell_n_in = CELL_N_IN, #calculated
        cell_n_hidden = cell_n_hidden, #calculated
        fusion_n_hidden = fusion_n_hidden,
        n_out=2,
        p=p,
        dropout=dropout,
        input_p=input_p
    )

    cost = (
        classifier.negative_log_likelihood(y)
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )

    validate_model = theano.function(
        inputs=[index],
        outputs=classifier.negative_log_likelihood(y),
        givens={
            x_c: valid_cell_x[valid_cell_index_x,],
            y: valid_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_model = theano.function(
        inputs=[index],
        outputs=classifier.errors(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_pred = theano.function(
        inputs=[index],
        outputs=classifier.pred(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.errors(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN"+ "\t"+"VALID.ERROR" + "\t" + "TEST.ERROR")
    FILE_OUT.close()

    FILE_OUT_val = open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "w")
    FILE_OUT_val.write("EPOCH" +"\t" + "ACTUAL" +"\t"+"PREDICTED")
    FILE_OUT_val.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            # iteration number
            #iter = (epoch - 1) * n_train_batches + minibatch_index

            #if (iter + 1) % validation_frequency == 0:
            if (minibatch_index + 1) % EPOCH_SIZE == 0:
                # compute zero-one loss on validation set

                validation_losses = [validate_model(i) for i in xrange(n_valid_batches)]
                this_validation_loss = np.mean(validation_losses)

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)


                log = ('epoch %i, minibatch %i/%i, train error %f ,validation error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error ,
                        this_validation_loss
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error) + "\t"+ str(this_validation_loss) \
                                   + "\t" + str(test_loss))

                # if we got the best validation score until now
                if this_validation_loss < best_validation_loss:
                    LR_COUNT = 0

                    #improve patience if loss improvement is good enough
                    # if (
                    #     this_validation_loss < best_validation_loss *
                    #     improvement_threshold
                    # ):
                    #     patience = max(patience, iter * patience_increase)

                    best_validation_loss = this_validation_loss
                    best_iter = iter

                    # test it on the test set
                    test_losses = [test_model(i) for i in xrange(n_test_batches)]
                    test_loss = np.mean(test_losses)

                    log = ((' epoch %i, minibatch %i/%i, test error of '
                        'best loss %f %%') %
                        (epoch, minibatch_index + 1, EPOCH_SIZE, test_loss))
                    # print(log)
                    with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                        logfile.write(log + "\n")

                    #ONLY SAVE MODEL if validation improves
                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["multiplicative"]  = classifier.multiplicative_input
                    MODEL["fusion_n_hidden"] = [getattr(classifier, "fusion_layer_" + str(e)) for e in xrange(len(fusion_n_hidden))]
                    MODEL["logistic"]        = classifier.logRegressionLayer

                    # MODEL = []
                    # for e in xrange(len(cell_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "cell_layer_" + str(e))]
                    # for e in xrange(len(drug_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "drug_layer_" + str(e))]
                    # for e in xrange(len(fusion_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "fusion_layer_" + str(e))]
                    #
                    # MODEL = MODEL + [classifier.logRegressionLayer]
                    with open(OUT_FOLDER + "/" + drug_name + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                    #Only write if validation improvement
                    ACTUAL = test_set_y.eval()
                    PREDICTED = [test_pred(i) for i in xrange(n_test_batches)][0]

                    with open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "a") as FILE_OUT_val:
                        for l in xrange(len(ACTUAL)):
                            FILE_OUT_val.write("\n" + str(epoch) + "\t" + str(ACTUAL[l]) + "\t" + str(PREDICTED[l]))
                else:
                    LR_COUNT = LR_COUNT+1

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))

def class_mlp_mf(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=20,
             cell_n_hidden=[500,200,100], drug_n_hidden=[500,200,100], mf_manual="None", fusion_n_hidden = [500,200,100],
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):

    #Demonstrate stochastic gradient descent optimization for a multilayer
    #perceptron for parallel drug and cell layers to Multiplicative fusion layer
    train_drug_x, train_cell_x, train_drug_index_x, train_cell_index_x, train_set_y = datasets[0]
    valid_drug_x, valid_cell_x, valid_drug_index_x, valid_cell_index_x, valid_set_y = datasets[1]
    test_drug_x,  test_cell_x,  test_drug_index_x,  test_cell_index_x,  test_set_y = datasets[2]

    valid_batch_size = valid_drug_index_x.eval().shape[0] # Could have been valid_cell_index_x since they are of equal length
    test_batch_size  = test_drug_index_x.eval().shape[0]
    train_samples    = train_drug_index_x.eval().shape[0] # Could have been train_cell_index_x since they are of equal length

    # Compute input layer for both drug and cell networks
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]
    DRUG_N_IN = train_drug_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size
    n_valid_batches = valid_batch_size / valid_batch_size #1
    n_test_batches  = test_batch_size / test_batch_size #1

    # compute fusion neural range
    #neural_range = range(drug_n_hidden[-1])
    if mf_manual=="None":

        if cell_n_hidden[-1] != drug_n_hidden[-1]:

            neural_range  = min([cell_n_hidden[-1], drug_n_hidden[-1]])
            if neural_range == cell_n_hidden[-1]:
                drug_n_hidden = drug_n_hidden + [neural_range] # In place modification!!!
            else :
                cell_n_hidden = cell_n_hidden + [neural_range] # In place modification!!!
        else:
            neural_range = cell_n_hidden[-1]
    else:
        neural_range  = mf_manual
        drug_n_hidden = drug_n_hidden + [mf_manual]
        cell_n_hidden = cell_n_hidden + [mf_manual]
    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')
    x_d = T.matrix('x_d')

    y = T.ivector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Class(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        drug_input = x_d, #needed
        cell_n_in = CELL_N_IN, #calculated
        drug_n_in = DRUG_N_IN, #calculated
        cell_n_hidden = cell_n_hidden, #calculated
        drug_n_hidden = drug_n_hidden, #calculated
        fusion_n_hidden = fusion_n_hidden,
        neural_range = neural_range, #calculated
        n_out=2,
        cell_p=p,
        cell_dropout=dropout,
        cell_input_p=input_p,
        drug_p=0.7,
        drug_dropout=True,
        drug_input_p=0.2
    )

    cost = (
        classifier.negative_log_likelihood(y)
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )

    validate_model = theano.function(
        inputs=[index],
        outputs=classifier.negative_log_likelihood(y),
        givens={
            x_c: valid_cell_x[valid_cell_index_x,],
            x_d: valid_drug_x[valid_drug_index_x,],
            y: valid_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_model = theano.function(
        inputs=[index],
        outputs=classifier.errors(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            x_d: test_drug_x[test_drug_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_pred = theano.function(
        inputs=[index],
        outputs=classifier.pred(y),
        givens={
            x_c: test_cell_x[test_cell_index_x,],
            x_d: test_drug_x[test_drug_index_x,],
            y: test_set_y,
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            x_d: train_drug_x[train_drug_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.errors(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            x_d: train_drug_x[train_drug_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN"+ "\t"+"VALID.ERROR" + "\t" + "TEST.ERROR")
    FILE_OUT.close()

    FILE_OUT_val = open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "w")
    FILE_OUT_val.write("EPOCH" +"\t" + "ACTUAL" +"\t"+"PREDICTED")
    FILE_OUT_val.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            # iteration number
            #iter = (epoch - 1) * n_train_batches + minibatch_index

            #if (iter + 1) % validation_frequency == 0:
            if (minibatch_index + 1) % EPOCH_SIZE == 0:
                # compute zero-one loss on validation set

                validation_losses = [validate_model(i) for i in xrange(n_valid_batches)]
                this_validation_loss = np.mean(validation_losses)

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)


                log = ('epoch %i, minibatch %i/%i, train error %f ,validation error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error ,
                        this_validation_loss
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error) + "\t"+ str(this_validation_loss) \
                                   + "\t" + str(test_loss))

                # if we got the best validation score until now
                if this_validation_loss < best_validation_loss:
                    LR_COUNT = 0

                    #improve patience if loss improvement is good enough
                    # if (
                    #     this_validation_loss < best_validation_loss *
                    #     improvement_threshold
                    # ):
                    #     patience = max(patience, iter * patience_increase)

                    best_validation_loss = this_validation_loss
                    best_iter = iter

                    # test it on the test set
                    test_losses = [test_model(i) for i in xrange(n_test_batches)]
                    test_loss = np.mean(test_losses)

                    log = ((' epoch %i, minibatch %i/%i, test error of '
                        'best loss %f %%') %
                        (epoch, minibatch_index + 1, EPOCH_SIZE, test_loss))
                    # print(log)
                    with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                        logfile.write(log + "\n")

                    #ONLY SAVE MODEL if validation improves
                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["drug_n_hidden"]   = [getattr(classifier, "drug_layer_" + str(e)) for e in xrange(len(drug_n_hidden))]
                    MODEL["multiplicative"]  = classifier.multiplicative_input
                    MODEL["fusion_n_hidden"] = [getattr(classifier, "fusion_layer_" + str(e)) for e in xrange(len(fusion_n_hidden))]
                    MODEL["logistic"]        = classifier.logRegressionLayer

                    # MODEL = []
                    # for e in xrange(len(cell_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "cell_layer_" + str(e))]
                    # for e in xrange(len(drug_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "drug_layer_" + str(e))]
                    # for e in xrange(len(fusion_n_hidden)):
                    #     MODEL = MODEL + [getattr(classifier, "fusion_layer_" + str(e))]
                    #
                    # MODEL = MODEL + [classifier.logRegressionLayer]
                    with open(OUT_FOLDER + "/" + str(epoch) + "_" + drug_name + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                    #Only write if validation improvement
                    ACTUAL = test_set_y.eval()
                    PREDICTED = [test_pred(i) for i in xrange(n_test_batches)][0]

                    with open(OUT_FOLDER + "/combined_D_values." + drug_name + ".txt", "a") as FILE_OUT_val:
                        for l in xrange(len(ACTUAL)):
                            FILE_OUT_val.write("\n" + str(epoch) + "\t" + str(ACTUAL[l]) + "\t" + str(PREDICTED[l]))
                else:
                    LR_COUNT = LR_COUNT+1

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))

def regression_mlp_train_mf(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=20,
             cell_n_hidden=[500,200,100], drug_n_hidden=[500,200,100], mf_manual = "None", fusion_n_hidden = [500,200,100],
             p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"
             ):
    #NOTE: TRAIN SET ONLY TRAINING

    train_drug_x, train_cell_x, train_drug_index_x, train_cell_index_x, train_set_y = datasets[0]
    train_samples    = train_drug_index_x.eval().shape[0] # Could have been train_cell_index_x since they are of equal length

    # Compute input layer for both drug and cell networks
    CELL_N_IN = train_cell_x.get_value(borrow=True).shape[1]
    DRUG_N_IN = train_drug_x.get_value(borrow=True).shape[1]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_samples / train_batch_size

    # compute fusion neural range
    #neural_range = range(drug_n_hidden[-1])
    # if mf_manual=="None":
    #
    #     if cell_n_hidden[-1] != drug_n_hidden[-1]:
    #
    #         neural_range  = min([cell_n_hidden[-1], drug_n_hidden[-1]])
    #         if neural_range == cell_n_hidden[-1]:
    #             drug_n_hidden = drug_n_hidden + [neural_range] # In place modification!!!
    #         else :
    #             cell_n_hidden = cell_n_hidden + [neural_range] # In place modification!!!
    #     else:
    #         neural_range = cell_n_hidden[-1]
    # else:
    #     neural_range  = mf_manual
    #     drug_n_hidden = drug_n_hidden + [mf_manual]
    #     cell_n_hidden = cell_n_hidden + [mf_manual]
    neural_range = cell_n_hidden[-1] + drug_n_hidden[-1] #NOTE: Addition of last two layers
    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x_c = T.matrix('x_c')
    x_d = T.matrix('x_d')

    y = T.fvector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    classifier = Multi_MLP_Regression(
        rng=rng,
        is_train = is_train, #needed
        cell_input = x_c, #needed
        drug_input = x_d, #needed
        cell_n_in = CELL_N_IN, #calculated
        drug_n_in = DRUG_N_IN, #calculated
        cell_n_hidden = cell_n_hidden, #calculated
        drug_n_hidden = drug_n_hidden, #calculated
        fusion_n_hidden = fusion_n_hidden,
        neural_range = neural_range, #calculated
        n_out=1,
        cell_p=p,
        cell_dropout=dropout,
        cell_input_p=input_p,
        drug_p=0.7,
        drug_dropout=True,
        drug_input_p=0.2
    )

    cost = (
        classifier.NRMSE(y) #pear_loss, NRMSE is more accurate if we are randomly sampling every mini-batch
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )
    ###################################

    #learning rate to shared
    learning_rate = theano.shared(np.cast[theano.config.floatX](learning_rate) )

    # momentum implementation stolen from
    # http://nbviewer.ipython.org/github/craffel/theano-tutorial/blob/master/Theano%20Tutorial.ipynb
    assert initial_momentum >= 0. and initial_momentum < 1.
    momentum =theano.shared(np.cast[theano.config.floatX](initial_momentum), name='momentum', borrow=True)

    # List of update steps for each parameter
    updates = []
    #Just gradient descent on cost
    for param in classifier.params:
        # For each parameter, we'll create a param_update shared variable.
        # This variable will keep track of the parameter's update step across iterations.
        # We initialize it to 0
        param_update = theano.shared(param.get_value()*0., broadcastable=param.broadcastable, borrow=True)
        # Each parameter is updated by taking a step in the direction of the gradient.
        # However, we also "mix in" the previous step according to the given momentum value.
        # Note that when updating param_update, we are using its old value and also the new gradient step.
        updates.append((param, param - learning_rate*param_update))
        # Note that we don't need to derive backpropagation to compute updates - just use T.grad!
        updates.append((param_update, momentum*param_update + (1. - momentum)*T.grad(cost, param)/(2*train_batch_size) ))

    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x_c: train_cell_x[train_cell_index_x[vector],],
            x_d: train_drug_x[train_drug_index_x[vector],],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x_c: train_cell_x[train_cell_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            x_d: train_drug_x[train_drug_index_x[index * train_batch_size:(index + 1) * train_batch_size],],
            y: train_set_y[index * train_batch_size:(index + 1) * train_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    ###############
    # TRAIN MODEL #
    ###############
    print '... training'

    # early-stopping parameters
    patience = 18000000 # look as this many examples regardless

    patience_increase = 2 # wait this much longer when a new best is found
    improvement_threshold = 0.995 # a relative improvement of this much is considered significant (default = 0.995)
    validation_frequency = min(n_train_batches, patience / 2)

    best_validation_loss = np.inf #MODIFIED
    best_iter = 0
    start_time = timeit.default_timer()

    epoch = 0
    done_looping = False
    test_loss = 1
    test_pear = 0
    LR_COUNT = 1

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN")
    FILE_OUT.close()

    EPOCH_SIZE = n_train_batches
    while (epoch < n_epochs) and (not done_looping):
        epoch = epoch + 1
        # print "momentum: ", momentum.get_value()
        # print "learning rate: ", learning_rate.get_value()
        log = "momentum: " + str(momentum.get_value()) + "; learning_rate: " + str(learning_rate.get_value())
        with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
            logfile.write(log + "\n")

        # if LR_COUNT==1000:
        #     new_learning_rate = learning_rate.get_value() * 0.2
        #     print new_learning_rate
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

        #for minibatch_index in xrange(n_train_batches):
        for minibatch_index in xrange(EPOCH_SIZE):

            ran_index = list(np.random.randint(low=0, high=train_samples-1, size=train_batch_size))
            minibatch_avg_cost = train_model(ran_index)

            rescale_weights(classifier.param_to_scale, 15.)

            if (minibatch_index + 1) % EPOCH_SIZE == 0:

                this_train_error = [train_error(i) for i in xrange(n_train_batches)]
                this_train_error = np.mean(this_train_error)

                log = ('epoch %i, minibatch %i/%i, train error %f %%' %
                    (
                        epoch,
                        minibatch_index + 1,
                        EPOCH_SIZE,
                        this_train_error
                    ))
                # print(log)
                with open(OUT_FOLDER + "/log." + drug_name + ".txt", "a") as logfile:
                    logfile.write(log + "\n")

                with open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "a") as FILE_OUT:
                    FILE_OUT.write("\n"+ str(epoch) + "\t" + str(this_train_error))


                if (epoch % 100)==0: #Write model every 200 epochs
                    LR_COUNT = 0

                    MODEL = {}
                    MODEL["cell_n_hidden"]   = [getattr(classifier, "cell_layer_" + str(e)) for e in xrange(len(cell_n_hidden))]
                    MODEL["drug_n_hidden"]   = [getattr(classifier, "drug_layer_" + str(e)) for e in xrange(len(drug_n_hidden))]
                    MODEL["multiplicative"]  = classifier.multiplicative_input
                    MODEL["fusion_n_hidden"] = [getattr(classifier, "fusion_layer_" + str(e)) for e in xrange(len(fusion_n_hidden))]
                    MODEL["linear"]          = classifier.linearRegressionLayer

                    with open(OUT_FOLDER + "/" + drug_name + "_" + str(epoch) + ".pkl", "wb") as f:
                        cPickle.dump(MODEL, f)

                else:
                    LR_COUNT = LR_COUNT+1 #MEANINGLESS

            # if patience <= iter:
            #     done_looping = True
            #     break
            # if LR_COUNT==100:
            #     done_looping = True
            #     break

        # adaption of momentum
        if momentum.get_value() < 0.99:
            new_momentum = 1. - (1. - momentum.get_value()) * 0.999
            momentum.set_value(np.cast[theano.config.floatX](new_momentum))
        # adaption of learning rate
        new_learning_rate = learning_rate.get_value() * 0.998
        learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))
        # if epoch%500 == 0:
        #     new_learning_rate = learning_rate.get_value() * 0.1
        #     learning_rate.set_value(np.cast[theano.config.floatX](new_learning_rate))

    end_time = timeit.default_timer()

    print(('Optimization complete. Best validation score of %f %% '
            'obtained at iteration %i, with test performance %f %%') %
            (best_validation_loss, best_iter + 1, test_loss ))

    print >> sys.stderr, ('The code for file ' +
                                os.path.split("__file__")[1] +
                                ' ran for %.2fm' % ((end_time - start_time) / 60.))

def shared_drug_dataset_IC50_mf(drug_data, cell_data, index_data, integers=True):

    cell_data     = cell_data.iloc[:,1:]
    cell_index    = list(index_data.cell)
    shared_cell_x = theano.shared(np.asarray(cell_data, dtype=theano.config.floatX), borrow=True)
    shared_cell_i = theano.shared(np.asarray(cell_index, dtype=theano.config.floatX), borrow=True)
    shared_cell_i = T.cast(shared_cell_i, 'int32')

    data_y     = list(index_data.NORM_pIC50)
    shared_y   = theano.shared(np.asarray(data_y, dtype=theano.config.floatX), borrow=True)

    if drug_data is not None:
        drug_data     = drug_data.iloc[:,1:]
        drug_index    = list(index_data.drug)
        shared_drug_x = theano.shared(np.asarray(drug_data, dtype=theano.config.floatX), borrow=True)
        shared_drug_i = theano.shared(np.asarray(drug_index, dtype=theano.config.floatX), borrow=True)
        shared_drug_i = T.cast(shared_drug_i, 'int32')

        if integers==True:
            return shared_drug_x, shared_cell_x, shared_drug_i, shared_cell_i, T.cast(shared_y, 'int32')
        else:
            return shared_drug_x, shared_cell_x, shared_drug_i, shared_cell_i, shared_y

    else:
        if integers==True:
            return shared_cell_x, shared_cell_i, T.cast(shared_y, 'int32')
        else:
            return shared_cell_x, shared_cell_i, shared_y

####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################

#OBTAIN SETTINGS (if available)
# if (len(sys.argv)==2):
#     n_epochs = 500
#     out_file = sys.argv[1]
# else:
#     n_epochs = int(sys.argv[2])
#     out_file = sys.argv[1] + "n_epoch_" + sys.argv[2]

n_epochs = 3000
out_file = sys.argv[1]

if sys.argv[6] != "0":

    DRUG_ALL     = int(sys.argv[6])
    DRUG_STD     = int(sys.argv[6]) * 2/3
    DRUG_THIRD   = int(sys.argv[6]) * 1/3
    DRUG_SQRT    = int(round(np.sqrt(int(sys.argv[6]))))
    DRUG_50      = int(sys.argv[6]) * 3/2

    if sys.argv[2][:6] == "manual":
        d_neurons    = [int(x) for x in sys.argv[2].split("_")[1:]]
    elif sys.argv[2]   == "sqrt":
        d_neurons    = [DRUG_SQRT]*2
    elif sys.argv[2] == "std_sqrt":
        d_neurons    = [DRUG_STD, DRUG_SQRT]
    elif sys.argv[2] == "std_std":
        d_neurons    = [DRUG_STD, DRUG_STD]
    elif sys.argv[2] == "std_third_sqrt":
        d_neurons    = [DRUG_STD, DRUG_THIRD, DRUG_SQRT]
    elif sys.argv[2] == "std_third":
        d_neurons    = [DRUG_STD, DRUG_THIRD]
    elif sys.argv[2] == "all":
        d_neurons    = [DRUG_ALL, DRUG_ALL]
    elif sys.argv[2] == "1_50":
        d_neurons    = [DRUG_50, DRUG_50]
    else :
        DRUG_NEURONS = int(sys.argv[2])
        d_neurons    = [DRUG_NEURONS]*1

    print "d_neurons " + str(d_neurons)

CELL_ALL     = int(sys.argv[7])
CELL_STD     = int(sys.argv[7]) * 2/3
CELL_THIRD   = int(sys.argv[7]) * 1/3
CELL_SQRT    = int(round(np.sqrt(int(sys.argv[7]))))
CELL_50      = int(sys.argv[7]) * 3/2

if sys.argv[3][:6] == "manual":
    c_neurons    = [int(x) for x in sys.argv[3].split("_")[1:]]
elif sys.argv[3] == "sqrt":
    c_neurons    = [CELL_SQRT]*2
elif sys.argv[3] == "std_sqrt":
    c_neurons    = [CELL_STD, CELL_SQRT]
elif sys.argv[3] == "std_std":
    c_neurons    = [CELL_STD, CELL_STD]
elif sys.argv[3] == "std_third_sqrt":
    c_neurons    = [CELL_STD, CELL_THIRD, CELL_SQRT]
elif sys.argv[3] == "std_third":
    c_neurons    = [CELL_STD, CELL_THIRD]
elif sys.argv[3] == "all":
    c_neurons    = [CELL_ALL, CELL_ALL]
elif sys.argv[3] == "1_50":
    c_neurons    = [CELL_50, CELL_50]
else :
    CELL_NEURONS = int(sys.argv[3])
    c_neurons    = [CELL_NEURONS]*1
print "c_neurons " + str(c_neurons)


FUSION_NEURONS = [int(x) for x in sys.argv[4].split("_")[1:]] #NOTE: Assumes in the form "manual_x_x_..."
print "FUSION_NEURONS " + str(FUSION_NEURONS)

if sys.argv[5] == "T":
    class_mlp = True
else:
    class_mlp = False

fold = sys.argv[8] #NOTE: Adding fold dependent feature. Determines how we train

if len(sys.argv)==10:
    mf_manual = int(sys.argv[9])
else:
    mf_manual = "None"
print "mf_manual " + str(mf_manual)

#OBTAIN FILES
file_name     = sys.argv[1]

IN_FOLDER="/tigress/zamalloa/CGP_FILES/TRAIN_TABLES/" #For tigress
# IN_FOLDER="/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/"

if class_mlp is True:
    OUT_FOLDER="/tigress/zamalloa/CGP_FILES/CLASS_RESULTS/" #For tigress
    # OUT_FOLDER="/home/zamalloa/Documents/FOLDER/CGP_FILES/CLASS_RESULTS/"
else:
    OUT_FOLDER="/tigress/zamalloa/CGP_FILES/REGRESSION_RESULTS/" #For tigress
    # OUT_FOLDER="/home/zamalloa/Documents/FOLDER/CGP_FILES/REGRESSION_RESULTS/"

# EXECUTE LEARNING - Do we have any drug features?
if sys.argv[6] != "0":

    train_drug  = pd.read_csv(IN_FOLDER + file_name + "_train_drug", sep="\t")
    train_cell  = pd.read_csv(IN_FOLDER + file_name + "_train_cell", sep="\t")
    train_index = pd.read_csv(IN_FOLDER + file_name + "_train_index", sep="\t")

    train_drug, train_cell, train_drug_index, train_cell_index, train_set_y = shared_drug_dataset_IC50_mf(train_drug, train_cell, train_index, integers=class_mlp)
    train_samples = train_drug_index.eval().shape[0]

    if fold=="fold_none":
        # That is, if we are training with all the data, no early dropout, so there only exists a train set, no valid or test sets
        drugval = [(train_drug, train_cell, train_drug_index, train_cell_index, train_set_y)]

        if class_mlp is True:

            class_mlp_train_mf(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=50,
                         cell_n_hidden=c_neurons, drug_n_hidden= d_neurons, mf_manual=mf_manual, fusion_n_hidden = FUSION_NEURONS,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)
        else:

            regression_mlp_train_mf(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=50,
                         cell_n_hidden=c_neurons, drug_n_hidden= d_neurons, mf_manual=mf_manual, fusion_n_hidden = FUSION_NEURONS,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)

    else:

        valid_drug  = pd.read_csv(IN_FOLDER + file_name + "_valid_drug", sep="\t")
        valid_cell  = pd.read_csv(IN_FOLDER + file_name + "_valid_cell", sep="\t")
        valid_index = pd.read_csv(IN_FOLDER + file_name + "_valid_index", sep="\t")

        test_drug   = pd.read_csv(IN_FOLDER + file_name + "_test_drug", sep="\t")
        test_cell   = pd.read_csv(IN_FOLDER + file_name + "_test_cell", sep="\t")
        test_index  = pd.read_csv(IN_FOLDER + file_name + "_test_index", sep="\t")

        valid_drug, valid_cell, valid_drug_index, valid_cell_index, valid_set_y = shared_drug_dataset_IC50_mf(valid_drug, valid_cell, valid_index, integers=class_mlp)
        test_drug,  test_cell,  test_drug_index,  test_cell_index,  test_set_y  = shared_drug_dataset_IC50_mf(test_drug,  test_cell,  test_index, integers=class_mlp)

        drugval = [(train_drug, train_cell, train_drug_index, train_cell_index, train_set_y),
                   (valid_drug, valid_cell, valid_drug_index, valid_cell_index, valid_set_y),
                   (test_drug,  test_cell,  test_drug_index,  test_cell_index,  test_set_y)]

        if class_mlp is True:

            class_mlp_mf(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=50,
                         cell_n_hidden=c_neurons, drug_n_hidden= d_neurons, mf_manual=mf_manual, fusion_n_hidden = FUSION_NEURONS,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)
        else:

            regression_mlp_mf(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=50,
                         cell_n_hidden=c_neurons, drug_n_hidden= d_neurons, mf_manual=mf_manual, fusion_n_hidden = FUSION_NEURONS,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)
else:
    #If no drug features

    print("zero drug neurons")
    train_cell  = pd.read_csv(IN_FOLDER + file_name + "_train_cell", sep="\t")
    train_index = pd.read_csv(IN_FOLDER + file_name + "_train_index", sep="\t")
    train_cell, train_cell_index, train_set_y = shared_drug_dataset_IC50_mf(None, train_cell, train_index, integers=class_mlp)

    if fold=="fold_none":
        # That is, if we are training with all the data, no early dropout, so there only exists a train set, no valid or test sets
        drugval = [(train_cell, train_cell_index, train_set_y)]
        print("fold_none")
        if class_mlp is True:

            class_mlp_train_mf_zero_drug(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=10,
                         cell_n_hidden=c_neurons, mf_manual=mf_manual,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)
        else:

            regression_mlp_train_mf_zero_drug(learning_rate=1.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=20,
                         cell_n_hidden=c_neurons, mf_manual=mf_manual,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)

    else:
        valid_cell  = pd.read_csv(IN_FOLDER + file_name + "_valid_cell", sep="\t")
        valid_index = pd.read_csv(IN_FOLDER + file_name + "_valid_index", sep="\t")

        test_cell   = pd.read_csv(IN_FOLDER + file_name + "_test_cell", sep="\t")
        test_index  = pd.read_csv(IN_FOLDER + file_name + "_test_index", sep="\t")

        valid_cell, valid_cell_index, valid_set_y = shared_drug_dataset_IC50_mf(None, valid_cell, valid_index, integers=class_mlp)
        test_cell,  test_cell_index,  test_set_y  = shared_drug_dataset_IC50_mf(None,  test_cell,  test_index, integers=class_mlp)

        drugval = [(train_cell, train_cell_index, train_set_y),
                   (valid_cell, valid_cell_index, valid_set_y),
                   (test_cell,  test_cell_index,  test_set_y)]

        #DEEP LEARNING WITHOUT DROPOUT
        if class_mlp is True:

            class_mlp_mf_zero_drug(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=10,
                         cell_n_hidden=c_neurons, mf_manual=mf_manual,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)
        else:

            regression_mlp_mf_zero_drug(learning_rate=0.1, L1_reg=0, L2_reg=0.0000000, n_epochs=n_epochs, initial_momentum=0.5, input_p=0.2,
                         datasets=drugval, train_batch_size=10,
                         cell_n_hidden=c_neurons, mf_manual=mf_manual,
                         p=0.7, dropout=True,
                         drug_name=out_file,
                         OUT_FOLDER = OUT_FOLDER)


print "DONE"
