#############################################################Testing H20 Package for Neural Networks########################################################
###Example
## Start a local cluster with 4GB RAM
localH2O = h2o.init(ip = "localhost", port = 54321, startH2O = TRUE, Xmx = '4g') 

#Load data in package mlbench
data(BreastCancer)

## Convert Breast Cancer into H2O
dat <- BreastCancer[, -1] # remove the ID column
head(dat)
dim(dat)
dat_h2o <- as.h2o(localH2O, dat, key = 'dat') 

## Split 60/40 for Training/Test Dataset
set.seed(1234)
y_all <- as.matrix(dat_h2o$Class)
rand_folds <- createFolds(as.factor(y_all), k = 5)
row_train <- as.integer(unlist(rand_folds[1:3]))
row_test <- as.integer(unlist(rand_folds[4:5]))
y_train <- as.factor(y_all[row_train])
y_test <- as.factor(y_all[row_test])

###### Train the model without dropout######
model <- h2o.deeplearning(x = 1:9, # column numbers for predictors
                          y = 10, # column number for label
                          data = dat_h2o[row_train, ],
                          activation = "Tanh", #tanh without dropout 
                          balance_classes = TRUE,
                          hidden = c(100,100), ## three hidden layers
                          epochs = 500) #number of epochs or iterations of training data

## Evaluate 
yhat_train <- h2o.predict(model, dat_h2o[row_train, ])$predict
yhat_train <- as.factor(as.matrix(yhat_train))
yhat_test <- h2o.predict(model, dat_h2o[row_test, ])$predict
yhat_test <- as.factor(as.matrix(yhat_test))

confusionMatrix(yhat_train, y_train)$overall[1]
confusionMatrix(yhat_test, y_test)$overall[1]

###### Train the model with dropout######
model <- h2o.deeplearning(x = 1:9, # column numbers for predictors
                          y = 10, # column number for label
                          data = dat_h2o[row_train, ],
                          activation = "TanhWithDropout", #as name states, tanh with dropout
                          input_dropout_ratio = 0.2,
                          hidden_dropout_ratios = c(0.5,0.5, 0.5), #for each layer
                          balance_classes = TRUE,
                          hidden = c(50,50,50), ## three hidden layers
                          epochs = 500)

###############################################################################################################################################################