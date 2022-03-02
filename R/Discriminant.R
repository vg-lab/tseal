# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' @export
testModel <- function(x,...) {
  UseMethod("testModel")
}
#' @export
testModel.WaveAnalisys <- function (train,test,grps,method, returnClassification = FALSE) {
  if (missing(train)) {
    stop("\"train\" must be provided. \"train\" must be an object of class WaveAnalisys")
  }
  if (missing(test)) {
    stop("\"test\" must be provided. \"test\" must be an object of class WaveAnalisys")
  }
  if (missing(grps)) {
    stop("\"grps\" must be provided.")
  }
  if (missing(method)) {
    stop("\"method\" must be provided. The avaiable options are \"linear\" or \"quadratic\"")
  }

  stopifnot(class(train) == "WaveAnalisys",
            class(test) == "WaveAnalisys")


  model <- trainModel(train,grps,method)
  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(grps))

  if (returnClassification) {
    return(list("CM" = CM, "Clasification" = prediction[[1]]))
  }

  return(CM)
}
#' @export
testModel.lda <- function (model,test,grps,returnClassification = FALSE) {
  stopifnot(!missing(model),
            !missing(test),
            !missing(grps),
            class(model) == "lda" || class(model) == "qda"
            )

  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(grps))

  if (returnClassification) {
    return(list("CM" = CM, "Clasification" = prediction[[1]]))
  }

  return(CM)
}

#' @export
testModel.qda <- function (model,test,grps,returnClassification = FALSE) {
  testModel.lda(model,test,grps,returnClassification)
}

#' @export
LOOCV <- function(x,...){
  UseMethod("LOOCV")
}

#' @export
LOOCV.array <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0,features = c("Var","Cor","IQR","PE","DM"),returnClassification = FALSE,nCores = 0) {
  if (missing(XSeries1)) {
    stop("XSeries1 must be provided")
  }

  if (missing(XSeries2)) {
    stop("XSeries2 must be provided")
  }

  if(missing(f)) {
    stop("The parameter \"f\" (filter) must be defined. To see the avaiable filters use avaibleFilters function")
  }

  if (missing(method)){
    stop("The parameter \"method\" must be defined. The avaiable  methods are \"linear\" and \"quadratic\"")
  }

  MWA <- generateStepDiscrim(XSeries1,XSeries2,f,method,maxvars,Vstep,lev,features,nCores)
  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))
  return(LOOCV(MWA,grps,method,returnClassification))
}

#' LOOCV
#'
#' Performs a leave-one-cross-validation (LOOCV) method
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalisys function
#' @param grps labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options "Linear" "Quadratic"
#'
#' @export
#'
#' @importFrom caret confusionMatrix
LOOCV.WaveAnalisys <- function(MWA,grps,method, returnClassification = FALSE) {
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
#' @return
#' A vector containing one confusion matrix object for each group
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

trainModel <- function(x, ...) {
  UseMethod("trainModel")
}

#' Train model
#' Generates a prediction model from training data.
#'
#' It generates a prediction model starting from the training data, which must
#'  be provided in 2 groups depending on their classification. The method first
#'  obtains the variances and correlations using MODWT, the f filter is applied
#'  with a number of levels lev. Then a subset of all the generated features
#'  will be obtained by means of a stepwise discriminant, which can be driven
#'  by a maximum number of features or by a minimum metric to be met. Finally,
#'  the selected prediction model is trained with the subset obtained.
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
#' @param Var Determines whether variances are to be taken into account for the classifier.
#' @param Cor Determines whether correlations are to be taken into account for the classifier.
#'
#' @return a prediction model object
#'
#' @export
#' @seealso
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' @md

trainModel.array <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0, features = c("Var","Cor","IQR","PE","DM"),nCores = 0) {
  MWA <- generateStepDiscrim(XSeries1,XSeries2,f,method,maxvars,Vstep,lev,features,nCores)
  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))
  return (trainModel(MWA,grps,method))

}
#' train Model
#' Generates a prediction model from an already generated "WaveAnalysis".
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalisys function
#' @param groups labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options "Linear" "Quadratic"
#'
#' @return a prediction model based on selected method.
#' @export
#'
#' @importFrom MASS lda qda
trainModel.WaveAnalisys <- function (MWA,groups,method) {
  stopifnot(class(MWA) == "WaveAnalisys")
  values = values(MWA)
  if (method == "linear"){
    model <- lda(t(values),groups)
  } else if (method == "quadratic"){
    model <- qda(t(values),groups)
  } else{
    stop (paste(c("Method",method ,"not supported")))
  }
  return(model)
}

#' @importFrom stats predict
#' @importFrom magrittr %>%
classify <- function (model,data) {
  stopifnot(class(data) == "WaveAnalisys")
  stopifnot(class(model) == "lda" || class(model) =="qda")

  values = values(data)
  return(predict(model,t(values)))
}
