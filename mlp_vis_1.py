# mlp_vis_1.py

import pandas as pd
import os
import sys
import timeit
import numpy as np
import theano
import cPickle
import itertools
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
        if (activation == T.tanh) or (activation == relu):
            output = activation(lin_output)
            self.params = [self.W, self.b]

        elif activation == prelu:
            output = activation(lin_output, self.alpha)
            self.params = [self.W, self.b, self.alpha]
        else:
            output = lin_output

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

        # neural_range is the neuron range from drug_out
        # (ie. precalculated range of last drug hidden layer)

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
        output      = drug_output * cell_output # Keep in mind that both have to have the same number of neurons
        #output       = T.concatenate([drug_output, cell_output], axis = 1)
        # Output
        self.output = output

        # Parameters of the fusion
        self.params = [self.cell_alpha, self.drug_alpha, self.cell_beta, self.drug_beta]

class Multiplicative_fusion_zero_drug(object):
    # Performs combinations of two neural layers (drug_input, cell_input) from
    # different sources and output the Multiplicative_fusion layer along with
    # 4 parameters to be learned.

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

        # PROCESS DRUG INPUT NEXT
        if input_p!=None:
            self.drug_input_layer = drop(drug_input, rng=rng, p=input_p)
            self.drug_input_layer = T.switch(T.neq(is_train, 0), self.drug_input_layer, drug_input)
        else:
            self.drug_input_layer = drug_input

        self.drug_layer_0 = HiddenLayer(
            rng=rng,
            input=self.drug_input_layer,
            n_in=drug_n_in,
            n_out=drug_n_hidden[0],
            activation=prelu,
            is_train=is_train,
            p=p,
            dropout=dropout
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
                                                    activation=prelu,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
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
            neural_range = neural_range
        )
        self.params = self.params + self.multiplicative_input.params

        # PROCESS FUSION LAYERS
        self.fusion_layer_0 = HiddenLayer(
            rng=rng,
            input=self.multiplicative_input.output,
            n_in=neural_range, #drug_n_hidden[drug_layer_number-1],
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
        self.linearRegressionLayer = LinearRegression(
            input=getattr(self, "fusion_layer_" + str(fusion_layer_number-1)).output,
            n_in=fusion_n_hidden[fusion_layer_number-1],
            n_out=n_out,
            rng=rng
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
        self.loss = self.linearRegressionLayer.loss
        self.NRMSE = self.linearRegressionLayer.NRMSE
        self.pred = self.linearRegressionLayer.pred

        self.param_to_scale = param_to_scale

        self.input = input #KEEP IN MIND THIS IS DIFFERENT THAN self.input_layer!!!

class Multi_MLP_Class(object):
    def __init__(self, rng, cell_input, drug_input, is_train,
                 cell_n_in, drug_n_in, cell_n_hidden, drug_n_hidden, fusion_n_hidden, neural_range,
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
            activation=T.tanh,
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
                                                    activation=T.tanh,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
                                                )

                setattr(self, "cell_layer_" + str(cell_layer_number), current_hidden_layer)

                self.params = self.params + getattr(self, "cell_layer_" + str(cell_layer_number)).params

                param_to_scale = param_to_scale + [getattr(self, "cell_layer_" + str(cell_layer_number)).params[0]]

                cell_layer_number = cell_layer_number + 1

        # PROCESS DRUG INPUT NEXT
        if input_p!=None:
            self.drug_input_layer = drop(drug_input, rng=rng, p=input_p)
            self.drug_input_layer = T.switch(T.neq(is_train, 0), self.drug_input_layer, drug_input)
        else:
            self.drug_input_layer = drug_input

        self.drug_layer_0 = HiddenLayer(
            rng=rng,
            input=self.drug_input_layer,
            n_in=drug_n_in,
            n_out=drug_n_hidden[0],
            activation=T.tanh,
            is_train=is_train,
            p=p,
            dropout=dropout
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
                                                    activation=T.tanh,
                                                    is_train=is_train,
                                                    p=p,
                                                    dropout=dropout
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
            activation=T.tanh,
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
                                                    activation=T.tanh,
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

####### LOAD DATA #######
file_in    = sys.argv[1]
drug_layer = sys.argv[2]=="T"
out_folder = "/tigress/zamalloa/CGP_FILES/MLP_LAYERS/"
file_out   = out_folder + file_in.split("/")[-1].rstrip(".pkl")
model      = cPickle.load(open(file_in, "rb"))

# Obtain weights per layer
cell_layers = [[a.W.get_value(), a.b.get_value()]  for a in model["cell_n_hidden"]]
mf_params   = [model["multiplicative"].cell_alpha.get_value(), model["multiplicative"].cell_beta.get_value()]
mf_layers   = [[a.W.get_value(), a.b.get_value()] for a in model["fusion_n_hidden"]]
log_layer   = [model["logistic"].W.get_value(), model["logistic"].b.get_value()]

if drug_layer==True:
    drug_layers = [[a.W.get_value(), a.b.get_value()]  for a in model["drug_n_hidden"]]
    mf_params   = mf_params + [model["multiplicative"].drug_alpha.get_value(), model["multiplicative"].drug_beta.get_value()]


# print([a[0].shape for a in cell_layers])
# print([a[0].shape for a in drug_layers])
# print([a.shape for a in mf_params])
# print([a[0].shape for a in mf_layers])
# print(log_layer)
# print([a.shape for a in log_layer])

####### SAVE DATA #######
layer=1
for i in cell_layers:
    #Save W and b
    np.savetxt(file_out + "_cell_W" +str(layer)  ,i[0], delimiter = "\t")
    np.savetxt(file_out + "_cell_b" +str(layer)  ,i[1], delimiter = "\t")
    layer+=1

if drug_layer==True:
    layer=1
    for i in drug_layers:
        #Save W and b
        np.savetxt(file_out + "_drug_W" +str(layer)  ,i[0], delimiter = "\t")
        np.savetxt(file_out + "_drug_b" +str(layer)  ,i[1], delimiter = "\t")
        layer+=1

np.savetxt(file_out + "_cell_alpha", mf_params[0], delimiter = "\t")
np.savetxt(file_out + "_cell_beta",  mf_params[1],  delimiter = "\t")
if drug_layer==True:
    np.savetxt(file_out + "_drug_alpha", mf_params[2], delimiter = "\t")
    np.savetxt(file_out + "_drug_beta",  mf_params[3],  delimiter = "\t")
    
layer=1
for i in mf_layers:
    np.savetxt(file_out + "_mf_W" +str(layer)  ,i[0], delimiter = "\t")
    np.savetxt(file_out + "_mf_b" +str(layer)  ,i[1], delimiter = "\t")
    layer+=1

np.savetxt(file_out + "_log_W", log_layer[0], delimiter = "\t")
np.savetxt(file_out + "_log_b", log_layer[1], delimiter = "\t")

print("Done writing layer parameters")