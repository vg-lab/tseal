# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' @importFrom utils read.csv
#' @importFrom abind abind
read_data <- function (infile,header) {
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

generateStepDiscrimWMA <- function(XSeries1,XSeries2,f,method, maxvars = 0, Vstep = 0, lev = 0, features = c("Var","Cor","IQR","PE","DM"),nCores = 0) {
  if (missing(maxvars) && missing(Vstep)){
    stop("maxvars o Vstep must be defined")
  }
  if (!missing(maxvars) && !missing(Vstep)){
    stop("only maxvars or Vstep can be defined.")
  }
  if (!missing(maxvars) && maxvars < 0) {
    stop("maxVars must be >0")
  }
  if (!missing(Vstep && Vstep <0)){
    stop("Vstep must be >0")
  }

  grps = rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,XSeries2[[3]],1))
  MWA <- MultiVaweAnalisys(XSeries1,XSeries2,f,lev,Var,Cor)
  if (!missing(maxvars)){
    MWA <- StepDiscrim(MWA,grps,maxvars,Var,Cor)
  } else {
    MWA <- StepDiscrimV(MWA,grps,Vstep,Var,Cor)
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
#' @param XSeries1 Samples from the population 1 (dim x length x cases)
#' @param XSeries2 Samples from the population 2 (dim x length x cases)
#' @param var Maximum number of selectable variables by stepwise discriminant
#' @param filters Vector que contiene los filtros a probar. Para ver los filtros disponibles ejecuta avaibleFilters()
#' @param lev Wavelet descomposition level, by default is selcted using the "conservative" strategy. See \code{\link{chooseLevel}} function.
#'
#' @return A list that each element contains:
#'   * CM: confusion matrix with a particular configuration using LOOCV
#'   * Vars: a boolean that determines if the variances have been taken into account in the variable selection process.
#'   * Cors: a boolean that determines if the Correlations have been taken into account in the variable selection process.
#'   * NVars: the total numbers of variables have been taken into account in the clasification process
#'   * Method: Type of classifier used.
#'   * Filter: filter used in the MultiWave analisys process
#'
#' @export
#'
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{MultiVaweAnalisys}}
#' * \code{\link{StepDiscrim}}
#'
#' @md
testFilters <- function(XSeries1,XSeries2,var,filters = c("haar","d4","d6","d8","la8"),features = c("Var","Cor","IQR","PE","DM"),lev = 0) { #change maxvars for vector, and add a vector of filters.
  data <- list()
  grps <- rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))
  nFeatures <- length(features)

  for (f in filters) {
    print(f)
    MWA <- MultiVaweAnalisys(XSeries1,XSeries2,f)
    for ( i in 1:nFeatures) {
      comFeatures <- combn(features,i)
      listFeatures <- split(comFeatures, rep(1:ncol(comFeatures), each = nrow(comFeatures)))
      for (cFeatures in listFeatures){
        print(cFeatures)
        aux <- StepDiscrim(MWA,grps,var,cFeatures,pos=TRUE)
        Tr <- aux[[1]]
        incl <- aux[[2]]
        print("Validation")
        maxVar <- min(var,length(incl))
        for (v in 2:maxVar){
          print(v)
          filterValues <- Tr[incl[1:v],]
          MWAAux <- list(Var = filterValues, Cor = NA, IQR = NA, DM = NA, PE = NA, Observations = MWA$Observations, NLevels = MWA$NLevels, filter = MWA$filter)
          attr(MWAAux,"class") <- "WaveAnalisys"
          data <- append(data,list(list(CM = LOOCV(MWAAux,grps,"linear"),NVars = v,Method = "linear",filter = f, Features = cFeatures)))
          data <- append(data,list(list(CM = LOOCV(MWAAux,grps,"quadratic"),NVars = v,Method = "quadratic",filter = f, Features = cFeatures)))
        }
      }
    }
  }
  return(data)
}
