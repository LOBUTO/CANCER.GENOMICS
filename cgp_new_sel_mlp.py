# cgp_new_sel_mlp.py

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

def test_mlp(learning_rate=0.01, L1_reg=0.00, L2_reg=0.0001, n_epochs=1000, initial_momentum = 0.5,
             datasets="datasets", train_batch_size=20,
             n_hidden=[500,200,100], p=0.5, dropout=False, input_p=None, drug_name=None, OUT_FOLDER="OUT_FOLDER"):

    #Demonstrate stochastic gradient descent optimization for a multilayer
    #perceptron

    train_set_x, train_set_y = datasets[0]
    valid_set_x, valid_set_y = datasets[1]
    test_set_x, test_set_y = datasets[2]
    #erlo_x, erlo_y = datasets[3] #MODIFIED

    valid_batch_size = valid_set_x.get_value(borrow=True).shape[0]
    test_batch_size= test_set_x.get_value(borrow=True).shape[0]
    N_IN=valid_set_x.get_value(borrow=True).shape[1]
    train_samples = train_set_x.get_value(borrow=True).shape[0]

    # compute number of minibatches for training, validation and testing
    n_train_batches = train_set_x.get_value(borrow=True).shape[0] / train_batch_size
    n_valid_batches = valid_set_x.get_value(borrow=True).shape[0] / valid_batch_size
    n_test_batches = test_set_x.get_value(borrow=True).shape[0] / test_batch_size


    ######################
    # BUILD ACTUAL MODEL #
    ######################
    print '... building the model'

    # allocate symbolic variables for the data
    index = T.lscalar("i") # index to a [mini]batch
    vector = T.vector("v", dtype='int32')
    x = T.matrix('x')
    y = T.vector('y')

    is_train = T.iscalar('is_train') # pseudo boolean for switching between training and prediction

    rng = np.random.RandomState(1234)

    # construct the MLP class
    N_HIDDEN = ".".join([str(NN) for NN in n_hidden])
    classifier = MLP(
        rng=rng,
        is_train = is_train,
        input=x,
        n_in=N_IN,   #FIXED !!!!!!
        n_hidden=n_hidden,
        n_out=1,
        p=p,
        dropout=dropout,
        input_p=input_p #, batch_size=batch_size
    )

    #classifier.negative_log_likelihood(y)
    cost = (
        classifier.errors(y)
        + L1_reg * classifier.L1
        + L2_reg * classifier.L2_sqr
    )

    validate_model = theano.function(
        inputs=[index],
        outputs=classifier.errors(y), #negative_log_likelihood(y)
        givens={
            x: valid_set_x[index * valid_batch_size:(index + 1) * valid_batch_size],
            y: valid_set_y[index * valid_batch_size:(index + 1) * valid_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_cor = theano.function(
        inputs=[index],
        outputs=classifier.loss(y),
        givens={
            x: test_set_x[index * test_batch_size:(index + 1) * test_batch_size],
            y: test_set_y[index * test_batch_size:(index + 1) * test_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_nrmse = theano.function(
        inputs=[index],
        outputs=classifier.NRMSE(y),
        givens={
            x: test_set_x[index * test_batch_size:(index + 1) * test_batch_size],
            y: test_set_y[index * test_batch_size:(index + 1) * test_batch_size],
            is_train: np.cast['int32'](0)
        },
        on_unused_input='warn',
    )

    test_pred = theano.function(
        inputs=[index],
        outputs=classifier.pred(y),
        givens={
            x: test_set_x[index * test_batch_size:(index + 1) * test_batch_size],
            y: test_set_y[index * test_batch_size:(index + 1) * test_batch_size],
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

    """
    gparams = [T.grad(cost, param) for param in classifier.params]

    updates = [
        (param, param - learning_rate * gparam)
        for param, gparam in zip(classifier.params, gparams)
    ]
    """
    train_model = theano.function(
        inputs=[vector],
        outputs=cost,
        updates=updates,
        givens={
            x: train_set_x[vector,],
            y: train_set_y[vector,],
            is_train: np.cast['int32'](1)
        },
        on_unused_input='warn',
    )

    train_error = theano.function(
        inputs=[index],
        outputs=classifier.errors(y),
        givens={
            x: train_set_x[index * train_batch_size:(index + 1) * train_batch_size],
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

    # STORE_FILE="_LR"+str(learning_rate)+"_EPOCHS"+str(n_epochs) + "_BATCH_SIZE"+str(train_batch_size) + \
    # "_N_HIDDEN"+str(N_HIDDEN)+"_DROPOUT"+str(dropout)+"_P"+str(p)+"_IP"+str(input_p)
    #
    # STORE_RESULTS=open(OUT_FOLDER +"/"+ drug_name + STORE_FILE, "w")
    # STORE_RESULTS.write("LR"+"\t"+"EPOCHS"+"\t"+"BATCH_SIZE"+"\t"+
    #                     "L1"+"\t"+"L2"+"\t"+"N_HIDDEN"+"\t"+"P_HIDDEN"+"\t"+"DROPOUT"+"\t"+ "INPUT_DROPOUT"+"\t"+
    #                     "EPOCH_N"+"\t"+"BATCH_TYPE" + "\t" +"LOSS")

    FILE_OUT =  open(OUT_FOLDER + "/combined_D." + drug_name + ".txt", "w")
    FILE_OUT.write("EPOCH" + "\t" + "TRAIN"+ "\t"+"VALID.ERROR" + "\t" + "TEST.COR" + "\t" + "TEST.NRMSE")
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
                    MODEL = [classifier.linearRegressionLayer]
                    for e in xrange(len(n_hidden)):
                        MODEL = MODEL + [getattr(classifier, "layer_" + str(e))]
                    MODEL = MODEL + [rng]
                    with open(OUT_FOLDER + "/" + drug_name + ".pkl", "wb") as f:
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



import pandas as pd
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

    data_y = list(drug_data.NORM.pIC50)

    shared_x = theano.shared(np.asarray(data_x, dtype=theano.config.floatX), borrow=True)
    shared_y = theano.shared(np.asarray(data_y, dtype=theano.config.floatX), borrow=True )

    if integers==True:
        return shared_x, T.cast(shared_y, 'int32')
    else:
        return shared_x, shared_y

####################################################################################################################################################################################################
####################################################################################################################################################################################################
####################################################################################################################################################################################################

#OBTAIN FILES
# input_layers = sys.argv[1].split("f")[0].split("_")
# input_layers = [int(l) for l in input_layers]
# print(input_layers)
#input_name = ".".join([str(c) for c in input_layers])

drug   = sys.argv[1]

OUT_FOLDER="/tigress/zamalloa/CGP_FILES/CGP_RESULTS" #For tigress
OUT_FOLDER = "/home/zamalloa/Documents/FOLDER/CGP_FILES/CGP_NEW_RESULTS/" #For Lab

IN_FOLDER="/tigress/zamalloa/CGP_FILES/CGP_TRAIN_TABLES" #For tigress
IN_FOLDER = "/home/zamalloa/Documents/FOLDER/CGP_FILES/TRAIN_TABLES/" #For Lab

train_table = pd.read_csv(IN_FOLDER + "TRAIN_CGP_SEL." +  drug, sep="\t")
valid_table = pd.read_csv(IN_FOLDER + "VALID_CGP_SEL." +  drug, sep="\t")
test_table  = pd.read_csv(IN_FOLDER + "TEST_CGP_SEL."  +  drug, sep="\t")

train_drug_x, train_drug_y = shared_drug_dataset_IC50(train_table, integers=False)
valid_drug_x, valid_drug_y = shared_drug_dataset_IC50(valid_table, integers=False)
test_drug_x, test_drug_y   = shared_drug_dataset_IC50(test_table,  integers=False)

drugval= [(train_drug_x, train_drug_y), (valid_drug_x, valid_drug_y),(test_drug_x, test_drug_y)]

#DEEP LEARNING WITHOUT DROPOUT
NEURONS = (valid_drug_x.get_value(borrow=True).shape[1] +1)*2/3
#NEURONS = int(NEURONS * 2)
print NEURONS

for drop_out in [0.5]:

    for l in [2]:
        test_mlp(learning_rate=10.0, L1_reg=0, L2_reg=0.0000000, n_epochs=3000, initial_momentum=0.5, input_p=0.2,
                     datasets=drugval, train_batch_size=50,
                     n_hidden=[NEURONS]*l, p=drop_out, dropout=True,
                     drug_name="new_cgp_sel_model_" + drug,
                     OUT_FOLDER = OUT_FOLDER)

print "DONE"
