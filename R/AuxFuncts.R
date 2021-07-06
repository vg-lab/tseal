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
testFilters <- function(XSeries1,XSeries2,var,filters = c("haar","d4","d6","d8","la8"),lev = 0) { #change maxvars for vector, and add a vector of filters.
  data <- list()
  grps <- rbind(matrix(1,dim(XSeries1)[[3]],1),matrix(2,dim(XSeries2)[[3]],1))

  for (f in filters) {
    print(f)
    MWA <- MultiVaweAnalisys(XSeries1,XSeries2,f)
    print("Var&Cor")
    MWAVarCOR <- StepDiscrim(MWA,grps,maxvars = var)
    print("Var")
    MWAVar <- StepDiscrim(MWA,grps,maxvars = var,Var = TRUE, Cor = FALSE)
    print("Cor")
    MWACor <- StepDiscrim(MWA,grps,maxvars = var,Var = FALSE, Cor = TRUE)
    print("Discrim")
    for (v in 2:var) {
      print(v)

      aux <- MWAVarCOR
      aux$Values <- aux$Values[sort(aux$importance[1:v]),]
      data <- append(data,list(list(CM = LOOCV(aux,grps,"linear"),Vars = TRUE,Cors = TRUE, NVars = v,Method = "linear",Filter = f)))
      data <- append(data,list(list(CM = LOOCV(aux,grps,"quadratic"),Vars = TRUE,Cors = TRUE, NVars = v,Method = "quadratic",Filter = f)))

      aux <- MWAVar
      aux$Values <- aux$Values[sort(aux$importance[1:v]),]
      data <- append(data,list(list(CM = LOOCV(aux,grps,"linear"),Vars = TRUE,Cors = FALSE, NVars = v,Method = "linear",Filter = f)))
      data <- append(data,list(list(CM = LOOCV(aux,grps,"quadratic"),Vars = TRUE,Cors = FALSE, NVars = v,Method = "quadratic",Filter = f)))

      aux <- MWACor
      aux$Values <- aux$Values[sort(aux$importance[1:v]),]
      data <- append(data,list(list(CM = LOOCV(aux,grps,"linear"),Vars = FALSE,Cors = TRUE, NVars = v,Method = "linear",Filter = f)))
      data <- append(data,list(list(CM = LOOCV(aux,grps,"quadratic"),Vars = FALSE,Cors = TRUE, NVars = v,Method = "quadratic",Filter = f)))

    }
  }

  return(data)
}
