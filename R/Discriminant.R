# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21


testModel <- function (train,test,Tgroups,method) {
  model <- trainModel(train,Tgroups,method)
  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(Tgroups))
  return(CM)
}

testModel <- function (model,test,Tgroups) {
  prediction <- classify(model,test)
  CM <- confusionMatrix(as.factor(prediction[[1]]),as.factor(Tgroups))
  return(CM)
}

LOOCV <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0, Var = TRUE, Cor = TRUE) {
  if (missing(maxvars) && missing(Vstep)){
    simpleError("maxvars o Vstep must be defined")
  }
  if (!missing(maxvars) && !missing(Vstep)){
    simpleError("only maxvars or Vstep can be defined.")
  }
  if (!missing(maxvars) && maxVars < 0) {
    SimpleError("maxVars must be >0")
  }
  if (!missing(Vstep && Vstep <0)){
    simpleError("Vstep must be >0")
  }

  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,XSeries2[[3]],1))
  MWA <- MultiVaweAnalisys(Xeries1,Xseries2,f,lev,VAR,COR)
  if (!missing(maxvars)){
    MWA <- StepDiscrim(MWA,grps,maxvars,Var,Cor)
  } else {
    MWA <- StepDiscrimV(MWA,grps,VStep,Var,Cor)
  }

  return(LOOCV(MWA,grps,method))
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
LOOCV <- function(MWA,grps,method) {
  stopifnot(class(MWA) == "WaveAnalisys")
  n <- MWA$Observations
  class <- vector("numeric",n)
  for (i in 1:n) {
    aux <- extractSubset(MWA,c(i))
    MWATest <- aux[[1]]
    MWATrain <- aux[[2]]
    grpsT <-grps[-i]
    model <- trainModelMWA(MWATrain,grpsT,method)
    class[i] <- classify(model, MWATest)[[1]]
  }

  return (confusionMatrix(as.factor(class),as.factor(grps)))
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
    model <- trainModelMWA(MWATrain,grpsT,method)
    class <- classify(model,MWATest)[[1]]
    CMvector <- append(CMvector,confusionMatrix(as.factor(class), as.factor(reorderGrps[n:m])))
    n <- n + k
    m <- m + k
    m <- max(m,nobs)
  }
  return(CMvector)
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

trainModelData <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0, Var = TRUE, Cor = TRUE) {
  if (missing(maxvars) && missing(Vstep)){
    simpleError("maxvars o Vstep must be defined")
  }
  if (!missing(maxvars) && !missing(Vstep)){
    simpleError("only maxvars or Vstep can be defined.")
  }
  if (!missing(maxvars) && maxvars < 0) {
    simpleError("maxVars must be >0")
  }
  if (!missing(Vstep && Vstep <0)){
    simpleError("Vstep must be >0")
  }

  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,XSeries2[[3]],1))
  MWA <- MultiVaweAnalisys(XSeries1,XSeries2,f,lev,Var,Cor)
  if (!missing(maxvars)){
    MWA <- StepDiscrim(MWA,grps,maxvars,Var,Cor)
  } else {
    MWA <- StepDiscrimV(MWA,grps,Vstep,Var,Cor)
  }

  return (trainModelMWA(MWA,grps,method))

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
trainModelMWA <- function (MWA,groups,method) {
  stopifnot(class(MWA) == "WaveAnalisys")
  if (method == "linear"){
    model <- lda(t(MWA$Values),groups)
  } else if (method == "quadratic"){
    model <- qda(t(MWA$Values),groups)
  } else{
    simpleError ("Method not supported")
  }
  return(model)
}

#' @importFrom stats predict
#' @importFrom magrittr %>%
classify <- function (model,data) {
  stopifnot(class(data) == "WaveAnalisys")
  stopifnot(class(model) == "lda" || class(model) =="qda")

  return(model %>% predict(data$Values))
}
