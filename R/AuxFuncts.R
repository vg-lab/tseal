# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' @importFrom utils read.csv
#' @importFrom abind abind
read_data <- function (infile,header = TRUE) {
  if (missing(infile)) stop("\"infile must be provided\"")
  files <- list.files(infile, full.names = TRUE)

  data <- read.csv(files[1], header = header)
  files <- files[-1]
  for (currentFile in files) {
    aux <- read.csv(currentFile,header = header)
    data <- abind(data,aux,along = 3)
  }

  saveRDS(data, file = "totalExperimentsECG.rds")
  return(data)
}

D3toD2 <- function(i,j,k,nRows,nCols,nPages) {

  if (missing(i) && missing(j)) {
    return (list(k,list(0,nRows * nCols)))
  }

  if (missing(i)) {
    return (k,list)
  }

  if (missing(j)) {
    return (list(k, list((i-1)*nCols+1, (i) * nCols)))
  }

  return (list( k,list(i*nCols + j,i*nCols + j)))
}

#'  Generate StepDiscrim from raw data
#'
#'  This function allows to obtain in a single step the complete WaveAnalysis
#'  and the selection of the most discriminating variables of the WaveAnalysis.
#'
#' @param XSeries Sample from the population (dim x length x cases)
#' @param grps labeled vector that classify the observations
#' @param f Selected filter for the MODWT (to see the available filters use the function \code{\link{availableFilters}}
#' @param maxvars maximun number of variables included by the stepDiscrim algorithm (Note that if you defined this, can not define VStep)
#' @param Vstep Minimum value of V above which all other variables are considered
#'        irrelevant and therefore will not be included. (Note that if you defined
#'         this, can not defined maxvars) For more information see StepDiscrim documentation
#' @param lev Determines the number of decomposition levels for MODWT (by default the optimum is calculated).
#' @param features A list of characteristics that will be used for the classification process. To see the available features see \code{\link{availableFeatures}}
#' @param nCores Determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return A MultiWaveAnalysis with the most discriminant variables based on the features indicated.
#'
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#'
#' @examples
#'           ECGExample <- loadECGExample()
#'           # The dataset has the first 5 elements of class 1
#'           # and the last 5 of class 2.
#'           grps <- c(rep(1,5),rep(2,5))
#'           MWADiscrim <- generateStepDiscrim(ECGExample, grps, "haar", maxvars = 20)
#' @export
generateStepDiscrim <- function(XSeries,grps,f, maxvars = 0, Vstep = 0, lev = 0, features = c("Var","Cor","IQR","PE","DM"),nCores = 0) {
  if (missing(XSeries)) stop("The argument \"XSeries\" must be provided.")
  if (missing(f)) stop("The argument \"f\" (filter) must be provided. To see available filter use availableFilters()")

  if (missing(maxvars) && missing(Vstep)){
    stop("maxvars o Vstep must be defined")
  }
  if ((!missing(maxvars) && maxvars != 0) && (!missing(Vstep) && Vstep != 0)){
    stop("only maxvars or Vstep can be defined.")
  }
  if (!missing(maxvars) && maxvars < 0) {
    stop("maxVars must be >0")
  }
  if (!missing(Vstep) && Vstep <0){
    stop("Vstep must be >0")
  }

  MWA <- MultiWaveAnalysis(XSeries,f,lev,features,nCores)
  if (!missing(maxvars)){
    MWA <- StepDiscrim(MWA,grps,maxvars,features,nCores)
  } else {
    MWA <- StepDiscrimV(MWA,grps,Vstep,features,nCores)
  }

  return(MWA)
}




#' testFilters
#'
#' This function performs a test with a series of filters defined by the user,
#'  for the maximum number of variables determined. This function can be used to
#'  compare the performance of different filters with a different number of
#'  variables to be considered and the differences between a linear and a quadratic discriminant.
#'
#' @param XSeries Samples from the population (dim x length x cases)
#' @param grps labeled vector that classify the observations.
#' @param maxvars maximun number of variables included by the stepDiscrim algorithm
#' @param filters Vector indicating the filters to be tested. To see the available filters use the function \code{\link{availableFilters}}
#' @param features A list of characteristics that will be used for the classification process. To see the available features see \code{\link{availableFeatures}}
#' @param lev Wavelet descomposition level, by default is selcted using the "conservative" strategy. See \code{\link{chooseLevel}} function.
#'
#' @return A list that each element contains:
#'   * CM: confusion matrix with a particular configuration using LOOCV
#'   * Classification: a vector with the raw classification result. "1" if the observation belongs to the population 1 and "2" if belongs to the population 2.
#'   * NVars: the total numbers of variables have been taken into account in the classification process
#'   * Method: Type of classifier used.
#'   * Filter: filter used in the MultiWave analysis process
#'   * Features: Vector containing the features taken into account
#'
#' @examples
#' \dontrun{
#'           ECGExample <- loadECGExample()
#'           # The dataset has the first 5 elements of class 1
#'           # and the last 5 of class 2.
#'           grps <- c(rep(1,5),rep(2,5))
#'           result <- testFilters(ECGExample, grps, maxvars = 3)
#'           }
#'
#' @export
#'
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{availableFilters}}
#' * \code{\link{availableFeatures}}
#'
#' @md
testFilters <- function(XSeries, grps ,maxvars,filters = c("haar","d4","d6","d8","la8"),features = c("Var","Cor","IQR","PE","DM"),lev = 0) {
  if (missing(XSeries)) stop("XSeries must be provided")
  if (missing(grps)) stop("grps must be provided")
  if (missing(maxvars)) stop("maxvars must be provided")
  if (length(filters) == 0) stop("At least one filter must be provided. To see the available filters use the availableFilters()")
  if (length(features) == 0) stop("At least one feature must be provided. To see the available filters use the availableFeatures()")

  data <- list()
  nFeatures <- length(features)

  for (f in filters) {
    print(f)
    MWA <- MultiWaveAnalysis(XSeries,f)
    for ( i in 1:nFeatures) {
      comFeatures <- combn(features,i)
      listFeatures <- split(comFeatures, rep(1:ncol(comFeatures), each = nrow(comFeatures)))
      for (cFeatures in listFeatures){
        print(cFeatures)
        aux <- StepDiscrim(MWA,grps,maxvars,cFeatures,pos=TRUE)
        Tr <- aux[[1]]
        incl <- aux[[2]]
        print("Validation")
        maxVar <- min(maxvars,length(incl))
        for (v in 2:maxVar){
          print(v)
          filterValues <- Tr[incl[1:v],]
          MWAAux <- list(Features = list(Var = filterValues, Cor = NA, IQR = NA, DM = NA, PE = NA),
                         StepSelection = list(Var = NA,Cor = NA,IQR = NA,DM = NA, PE = NA),
                         Observations = MWA$Observations,
                         NLevels = MWA$NLevels,
                         filter = MWA$filter)
          attr(MWAAux,"class") <- "WaveAnalysis"
          aux <- LOOCV(MWAAux,grps,"linear",TRUE)
          data <- append(data,list(list(CM = aux[[1]],classification = aux[[2]],NVars = v,Method = "linear",filter = f, Features = cFeatures)))
          aux <- LOOCV(MWAAux,grps,"quadratic",TRUE)
          data <- append(data,list(list(CM = aux[[1]],classification = aux[[2]],NVars = v,Method = "quadratic",filter = f, Features = cFeatures)))
        }
      }
    }
  }
  return(data)
}

#' Loads a small example DataSet
#'
#'
#'
#' @return A matrix
#' @export
#'
#' @examples
#'          ECGExample <- loadECGExample()
#'          summary(ECGExample)
#' @importFrom utils data
loadECGExample <- function() {
  .myDataEnv <- new.env(parent=emptyenv())
  data(ECGExample,envir = .myDataEnv)
  return(.myDataEnv$ECGExample)
}

