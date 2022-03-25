# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' Applies a discriminant to perform a classification
#'
#' This function performs a classification using a discriminant (linear or quadratic)
#' from different inputs and returns the results of the classification performed by the discriminant.
#'
#' @seealso
#' * \code{\link{testModel.WaveAnalisys}}
#' * \code{\link{testModel.lda}}
#' * \code{\link{testModel.qda}}
#' @export
testModel <- function(x,...) {
  UseMethod("testModel")
}

#' Computes a classification from WaveAnalisys objects
#'
#' This function trains a discriminant (linear or quadratic) using two objects
#' of class WaveAnalysis used as training and test set.
#'
#' @param train WaveAnalisys class object to be used as training set for the selected model.
#' @param test WaveAnalisys class object to be used as test set.
#' @param grps Vector that determines the class to which each of the observations provided in the test set belongs.
#' @param method Select the discriminant used. The possible options are "linear" and "quadratic".
#' @param returnClassification Allows to select if the raw result classification is returned.
#'
#' @return if returnClassification is false return a object of class ConfunsionMatrix
#'         if returnClassification is true, it returns a list containing an object of the ConfusionMatrix class and a vector with the classification result.
#' @examples CM <- testModel(MWATrain,MWATest,c(1,2,2,1),"linear")
#'
#' @seealso
#' * \code{\link{testModel}}
#' * \code{\link{testModel.lda}}
#' * \code{\link{testModel.qda}}
#' * \code{\link{MultiVaweAnalisys}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' @export
testModel.WaveAnalisys <- function (train,test,grps,method, returnClassification = FALSE) {
  if (missing(train)) stop("The argument \"train\" must be provided. \"train\" must be an object of class WaveAnalisys")
  if (missing(test)) stop("The argument \"test\" must be provided. \"test\" must be an object of class WaveAnalisys")
  if (missing(grps)) stop("The arguemnt \"grps\" must be provided.")
  if (missing(method)) stop("The argument \"method\" must be provided. The avaiable options are \"linear\" or \"quadratic\"")
  stopifnot(class(train) == "WaveAnalisys",class(test) == "WaveAnalisys")
  if (test$Observations != length(grps)) stop("The number of observations in the test set does not correspond to the classes provided in the grps parameter.")


  model <- trainModel(train,grps,method)
  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(grps))

  if (returnClassification) {
    return(list("CM" = CM, "Clasification" = prediction[[1]]))
  }

  return(CM)
}


#' Computes a classification from a pretrained discriminant
#'
#' This function uses a pre-trained linear discriminant to classify a set of test data
#' . As output it returns a confusion matrix and optionally the raw classification result.
#'
#' @param model Trained linear discriminant. see \code{\link{trainModel}}
#' @param test WaveAnalisys class object to be used as test set.
#' @param grps Vector that determines the class to which each of the observations provided in the test set belongs.
#' @param returnClassification Allows to select if the raw result classification is returned.
#'
#' @return if returnClassification is false return a object of class ConfunsionMatrix
#'         if returnClassification is true, it returns a list containing an object of the ConfusionMatrix class and a vector with the classification result.
#'
#' @examples CM <- testModel(ldaDiscriminant,MWAtest,c(1,1,2,2))
#'
#' @seealso
#' \code{\link{testModel}}
#' \code{\link{testModel.WaveAnalisys}}
#'
#' @export
testModel.lda <- function (model,test,grps,returnClassification = FALSE) {
  if (missing(model)) stop("The argument \"model\" must be provided")
  if (missing(test)) stop("The argument \"test\" must be provided")
  if (missing(grps)) stop("The argument \"grps\" must be provided")
  stopifnot(class(model) == "lda" || class(model) == "qda")
  if (test$Observations != length(grps)) stop("The number of observations in the test set does not correspond to the classes provided in the grps parameter.")


  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(grps))

  if (returnClassification) {
    return(list("CM" = CM, "Clasification" = prediction[[1]]))
  }

  return(CM)
}

#' testModel.qda
#'
#' This function is similar to \code{\link{testModel.lda}} except that it uses a quadratic discriminant instead of a linear one.
#' @param model
#'
#' @export
testModel.qda <- function (model,test,grps,returnClassification = FALSE) {
  testModel.lda(model,test,grps,returnClassification)
}

#' Leave-One-Out Cross Validation
#'
#' This function performs the Leave-One-Out Cross Validation (LOOCV) process with different types of input parameters.
#' @seealso
#' \code{\link{LOOCV.array}}
#' \code{\link{LOOCV.WaveAnalisys}}
#' @export
LOOCV <- function(x,...){
  UseMethod("LOOCV")
}

#' Generates and validates a discriminant model generated directly from the data.
#'
#' It generates and validates a discriminant model starting from the data, which
#'  must be provided in 2 groups according to their classification. First, a
#'  WaveAnalisys object is obtained according to the selected characteristics,
#'  filter and levels. Then, the most important features are selected using a
#'  stepwise discriminant that allows to select a maximum number of variables
#'  (maxVars) or a minimum enhancement step (VStep). Finally, the model is
#'  trained using the subset of features and validated using
#'  Leave-One-Out Cross Validation (LOOCV).
#'
#' @param XSeries1 Sample from the population 1 (dim x length x cases)
#' @param XSeries2 Sample from the population 2 (dim x length x cases)
#' @param f Selected filter for the MODWT (to see the available filters use the function \code{\link{availableFilters}}
#' @param method Selected method for the discriminant. Valid values "Linear" "Quadratic"
#' @param maxvars maximun number of variables included by the stepDiscrim algorithm (Note that if you defined this, can not define VStep)
#' @param Vstep Minimum value of V above which all other variables are considered
#'        irrelevant and therefore will not be included. (Note that if you defined
#'         this, can not defined maxvars) For more information see StepDiscrim documentation
#' @param lev Determines the number of decomposition levels for MODWT (by default the optimum is calculated).
#' @param features A list of characteristics that will be used for the classification process. To see the available features see \code{\link{availableFeatures}}
#' @param returnClassification Allows to select if the raw result classification is returned.
#' @param nCores Determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return if returnClassification is false return a object of class ConfunsionMatrix
#'         if returnClassification is true, it returns a list containing an object of the ConfusionMatrix class and a vector with the classification result.
#'
#' @examples CM <- LOOCVC(XSeries1,XSeries2,"haar","linear")
#'
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{LOOCV.WaveAnalisys}}
#' * \code{\link{availableFilters}}
#' * \code{\link{availableFeatures}}
#' @export
LOOCV.array <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0,features = c("Var","Cor","IQR","PE","DM"),returnClassification = FALSE,nCores = 0) {
  if (missing(XSeries1)) stop("The argument \"XSeries1\" must be provided")
  if (missing(XSeries2)) stop("The argument \"XSeries2\" must be provided")
  if(missing(f)) stop("The argument \"f\" (filter) must be provided. To see the avaiable filters use availableFilters()")
  if (missing(method)) stop("The argument \"method\" must be defined. The available  methods are \"linear\" and \"quadratic\"")
  if (length(features) == 0) stop("At least one feature must be provided. To see the available filters use the availableFeatures()")

  MWA <- generateStepDiscrim(XSeries1,XSeries2,f,method,maxvars,Vstep,lev,features,nCores)
  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))
  return(LOOCV(MWA,grps,method,returnClassification))
}

#' LOOCV
#'
#' Performs a leave-one-cross-validation (LOOCV) method on a WaveAnalysis object.
#'  It is advisable to have selected a subset of all features (\code{\link{StepDiscrim}},\code{\link{StepDiscrimV}})
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalisys function
#' @param grps labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options "linear" "quadratic"
#' @param returnClassification Allows to select if the raw result classification is returned.
#'
#' @return if returnClassification is false return a object of class ConfunsionMatrix
#'         if returnClassification is true, it returns a list containing an object of the ConfusionMatrix class and a vector with the classification result.
#'
#' @examples CM <- LOOCV(MWA,c(1,1,2,2),"linear")
#'
#' @export
#'
#' @importFrom caret confusionMatrix
LOOCV.WaveAnalisys <- function(MWA,grps,method, returnClassification = FALSE) {
  if (missing(MWA)) stop("The argument \"MWA\" must be provided")
  if (missing(grps)) stop("The argument \"grps\" must be provided")
  if (missing(method)) stop("The argument \"method\" must be provided")
  stopifnot(class(MWA) == "WaveAnalisys")


  n <- MWA$Observations
  class <- vector("numeric",n)
  for (i in 1:n) {
    aux <- extractSubset(MWA,c(i))
    MWATest <- aux[[1]]
    MWATrain <- aux[[2]]
    grpsT <-grps[-i]
    model <- trainModel(MWATrain,grpsT,method)
    class[i] <- classify(model, MWATest)[[1]]
  }
  CM <- confusionMatrix(as.factor(class),as.factor(grps))

  if (returnClassification) {
    return(list("CM" = CM, "classification" = class))
  }

  return (CM)
}

#' KFCV
#'
#' Performs k-fold cross-validation where groups are chosen randomly.
#'  In case the value k is not divisor of the number of observations the last
#'  group will have nobs mod k observations.
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalisys function
#' @param grps labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options "Linear" "Quadratic"
#' @param k k parameter for the K-fold method.
#' @param seed seed for random selection of k groups
#'
#' @return A vector containing one confusion matrix object for each group
#' @export
#'
#'
# TODO--------- MWA$Values to values.
KFCV <- function(MWA,grps,method,k, seed = 10){
  stopifnot(class(MWA) == "WaveAnalisys")
  set.seed(seed)
  nobs <- MWA$Observations
  reorderMWA <- MWA
  cols <- sample(nobs)
  reorderMWA$Values <- reorderMWA$Values[,cols]
  reorderGrps <- grps[cols]
  CMvector <- vector(mode = "list")
  n <- 1
  m <- k
  while (m <= nobs){
    aux <- extractSubset(reorderMWA,n:m)
    MWATest <- aux[[1]]
    MWATrain <- aux[[2]]
    grpsT <- reorderGrps[-c(n:m)]
    model <- trainModel(MWATrain,grpsT,method)
    class <- classify(model,MWATest)[[1]]
    CMvector <- append(CMvector,confusionMatrix(as.factor(class), as.factor(reorderGrps[n:m])))
    n <- n + k
    m <- m + k
    m <- max(m,nobs)
  }
  return(CMvector)
}

#' Generate a Discriminant Model
#'
#' This function allows training of a discriminant model using different inputs
#'
#' @return A trained discriminant model
#'
#' @seealso
#' * \code{\link{trainModel.array}}
#' * \code{\link{trainModel.WaveAnalysis}}
#'
#' @export
#'
trainModel <- function(x, ...) {
  UseMethod("trainModel")
}

#' Generates a discriminant model from training data.
#'
#' It generates a discriminant model starting from the training data, which must
#'  be provided in 2 groups depending on their classification. The method first
#'  obtains the variances and correlations using MODWT, the f filter is applied
#'  with a number of levels lev. Then a subset of all the generated features
#'  will be obtained by means of a stepwise discriminant, which can be driven
#'  by a maximum number of features or by a minimum metric to be met. Finally,
#'  the selected discriminant model is trained with the subset obtained.
#'
#' @param XSeries1 Sample from the population 1 (dim x length x cases)
#' @param XSeries2 Sample from the population 2 (dim x length x cases)
#' @param f Selected filter for the MODWT (to see the available filters use the function avaibleFilters)
#' @param method Selected method for the discriminant. Valid values "Linear" "Quadratic"
#' @param maxvars maximun number of variables included by the stepDiscrim algorithm (Note that if you defined this, can not define VStep)
#' @param Vstep Minimum value of V above which all other variables are considered
#'        irrelevant and therefore will not be included. (Note that if you defined
#'         this, can not defined maxvars) For more information see StepDiscrim documentation
#' @param lev Determines the number of decomposition levels for MODWT (by default the optimum is calculated).
#' @param returnClassification Allows to select if the raw result classification is returned.
#' @param nCores Determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return a discriminant model object (lda or qda)
#'
#' @examples model <- trainModel(XSeries1,XSeries2,"d6","linear", maxvars = 20)
#'
#' @export
#' @seealso
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' * \code{\link{trainModel}}
#' @md

trainModel.array <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0, features = c("Var","Cor","IQR","PE","DM"),nCores = 0) {
  if (missing(XSeries1)) stop("The argument \"XSeries1\" must be provided")
  if (missing(XSeries2)) stop("The argument \"XSeries2\" must be provided")
  if(missing(f)) stop("The argument \"f\" (filter) must be provided. To see the avaiable filters use availableFilters()")
  if (missing(method)) stop("The argument \"method\" must be defined. The available  methods are \"linear\" and \"quadratic\"")
  if (length(features) == 0) stop("At least one feature must be provided. To see the available filters use the availableFeatures()")

  MWA <- generateStepDiscrim(XSeries1,XSeries2,f,method,maxvars,Vstep,lev,features,nCores)
  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))
  return (trainModel(MWA,grps,method))

}
#' Generates a discriminant model from an already generated "WaveAnalysis".
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalisys function
#' @param groups labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options are"linear" and "quadratic"
#'
#' @return a discriminant model based on selected method.
#'
#' @example model <- trainModel(MWA,c(1,1,2,2),"quadratic")
#'
#' @seealso
#' * \code{\link{MultiVaweAnalisys}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' @export
#'
#' @importFrom MASS lda qda
trainModel.WaveAnalisys <- function (MWA,grps,method) {
  if (missing(MWA)) stop("The argument \"MWA\" must be provided")
  if (missing(grps)) stop("The argument \"grps\" must be provided")
  if (missing(method)) stop("The argument \"method\" must be defined. The available  methods are \"linear\" and \"quadratic\"")

  stopifnot(class(MWA) == "WaveAnalisys")
  values = values(MWA)
  if (method == "linear"){
    model <- lda(t(values),grps)
  } else if (method == "quadratic"){
    model <- qda(t(values),grps)
  } else{
    stop (paste(c("Method",method ,"not supported")))
  }
  return(model)
}

#' Classifies observations based on a pre-trained model.
#'
#' This function allows to qualify observations based on a pre-trained model
#' that could have been obtained in several ways (such as using the train model
#' function).
#'
#' @param model Pre-trained discriminant model (lda or qda)
#' @param data Data to be classified by the model. Remember that it must be an
#' object of type WaveAnalisys. Note that it should have the same characteristics
#'  as those used to generate the model.
#'
#' @return A vector with predicted class of each observation
#'
#' @example prediction <- classify(model,data)
#'
#' @seealso
#' * \code{\link{trainModel}}
#'
#' @importFrom stats predict
#' @importFrom magrittr %>%
classify <- function (model,data) {
  if (missing(model)) stop("The argument \"model\" must be provided")
  if (missing(data)) stop("The argument  \"data\" must be provided")
  stopifnot(class(data) == "WaveAnalisys")
  stopifnot(class(model) == "lda" || class(model) =="qda")

  values = values(data)
  return(predict(model,t(values)))
}
