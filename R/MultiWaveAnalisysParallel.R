# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' Generate a MultiWave analisys
#'
#' Generates a multivariate analysis by calculating a series of features from
#'  the result of applying MODWT to the input data.
#'
#' @param Xseries Sample from the population (dim x length x cases)
#' @param f Selected wavelet filter for the analysis. To see the available filters use the function \code{\link{availableFilters}}
#' @param lev Wavelet decomposition level by default is selected using the "conservative" strategy. See \code{\link{chooseLevel}} function.
#' @param features It allows to select the characteristics to be calculated for the analysis. To see the available features use the function \code{\link{availableFeatures}}
#' @param nCores determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return A multivariate analysis with the characteristics indicated in the parameter features. This is an object of class WaveAnalysis
#'
#' @examples
#'    ECGExample <- loadECGExample()
#'    MWA <- MultiWaveAnalysis(ECGExample, "haar")
#'
#' @seealso
#' * \code{\link{availableFilters}}
#' * \code{\link{availableFeatures}}
#' @export
#'
#' @importFrom bigmemory as.matrix GetMatrixSize
#' @importFrom waveslim modwt wave.variance wave.correlation wave.filter
#' @importFrom parallel makeCluster clusterExport clusterEvalQ detectCores
#' @importFrom pryr object_size
#' @importFrom stats sd
#' @importFrom utils combn
MultiWaveAnalysis <- function(Xseries,f,lev = 0,features = c("Var","Cor","IQR","PE","DM"), nCores = 0) {
  if (missing(Xseries)) stop("The argument \"Xseries\" must be provided.")
  if (missing(f)) stop("The argument \"f\" (filter) must be provided. To see available filter use availableFilters()")
  if (length(features) == 0) stop("At least one feature must be provided. To see the available filters use availableFeatures()")
  if (length(dim(Xseries)) == 2) stop("It seems that a dimension is missing, in case your series contains only one case, make sure that you have activated the option \"drop = FALSE\" as in the following example Series1 = Series2 [,,1, drop = FALSE].")


  if (nCores == 0) {
    nCores = detectCores() - 1
  }

  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    nCores <- 2L
  }

  HVar <- "Var" %in% features
  HCor <- "Cor" %in% features
  HIQR <- "IQR" %in% features
  HPE  <- "PE"  %in% features
  HDM  <- "DM"  %in% features

  dim1 <- dim(Xseries)

  nv1 <- dim1[1]  #Time series Dimension
  nr1 <- dim1[2]  #Time series Lenght
  nc1 <- dim1[3]  #Sample Size

  nCores <- min(nCores, nc1)

  if (lev == 0) {
    lev <- chooseLevel("conservative",f,nr1)
  }

  #set parallel enviorement
  c <- parallel::makeCluster(nCores)
  mgrinit(c, boost = TRUE)
  result = tryCatch({
  makebarr(c, boost = TRUE)


  # Standarizing the data
  mgrmakevar(c,"YSeries1",nc1,nv1*nr1)

  for (i in 1:nc1) {
    YSeries1[i,] <- as.vector(t(Xseries[,,i]))
  }


  clusterExport(c,c("nv1","nr1","nc1"), envir= environment())

  clusterExport(c,c("D3toD2","computeIQR","computePermutationEntropy",
                    "computeDMeasure"), envir = loadNamespace("MTSC"))

  clusterEvalQ(c,{
    ids <- getidxs(nv1)
    for (i in ids) {
      for (j in 1:nc1) {
        aux <- D3toD2(i,,j,nv1,nr1,nc1)
        Vtemporal <- YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
        YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]] <- (Vtemporal - mean(Vtemporal)) / sd(Vtemporal)
      }
    }
    barr()
  })

  NVar <- nv1 * lev
  if (HCor || HDM) {
    NbK <- combn(1:nv1,2)
    NumberNbK <- dim(NbK)[2]
    NCor <- NumberNbK * lev
  } else {
    NCor <- 0
    NumberNbK <- 0
    NbK <- 0
  }


  clusterExport(c,c("NumberNbK","NbK","lev","f","HVar","HCor","HIQR","HPE","HDM"),envir = environment())
  if (HVar) {
    mgrmakevar(c,"Var",NVar,nc1)
  } else {
    Var = NA
  }
  if (HCor) {
    mgrmakevar(c,"Cor",NCor,nc1)
  } else {
    Cor = NA
  }
  if (HIQR) {
    mgrmakevar(c,"IQR",NVar,nc1)
  } else {
    IQR = NA
  }
  if (HPE) {
    mgrmakevar(c,"PE",NVar,nc1)
  } else {
    PE = NA
  }
  if (HDM) {
    mgrmakevar(c,"DM",NCor,nc1)
  } else {
    DM = NA
  }

  clusterEvalQ (c, {
    ids <- getidxs(nc1)
    print(ids)
    for (i in ids){
      WJ <- list()
      for (j in 1:nv1){
        aux = D3toD2(j,,i,nv1,nr1,nc1)
        Vtemporal <- YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
        aux.modwt <- waveslim::modwt(Vtemporal,f,lev)
        WJ <- append(WJ,list(aux.modwt))

        if (HVar) {
          WVARAux <- waveslim::wave.variance(aux.modwt)
          Var[((j-1) * lev + 1):(j*lev),i] <- WVARAux[1:lev,1]
        }
        if (HIQR) {
          WIQRAux <- computeIQR(aux.modwt)
          IQR[((j-1) * lev + 1):(j*lev),i] <- WIQRAux
        }
        if (HPE) {
          WPEAux <- computePermutationEntropy(aux.modwt)
          PE[((j-1) * lev + 1):(j*lev),i] <- WPEAux
        }
      }

      print("Cors")
      if (HCor || HDM) {
        for (k in 1:NumberNbK) {
          if (HCor){
            WCOR <- suppressWarnings(waveslim::wave.correlation(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]],nr1))
            Cor[((k-1)*lev + 1):(k * lev),i] <- WCOR[1:lev,1]
          }
          if (HDM) {
            WDM <- computeDMeasure(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]])
            DM[((k-1)*lev + 1):(k * lev),i] <- WDM
          }
        }
      }
      print("finishCor")
    }
    print("barr")
    barr()
  })
  }, finally = stoprdsm(c))

  aVar <- as.matrix(Var)
  aCor <- as.matrix(Cor)
  aIQR <- as.matrix(IQR)
  aDM <- as.matrix(DM)
  aPE <- as.matrix(PE)

  if(nc1 == 1) {
    aVar <- matrix(aVar)
    aCor <- matrix(aCor)
    aIQR <- matrix(aIQR)
    aDM <- matrix(aDM)
    aPE <- matrix(aPE)
  }

  x <- list(Features = list(Var = aVar, Cor = aCor, IQR = aIQR, DM = aDM, PE = aPE),
            StepSelection = list(Var = NA,Cor = NA,IQR = NA,DM = NA, PE = NA),
            Observations = nc1,
            NLevels = lev,
            filter = f)
  attr(x,"class") <- "WaveAnalysis"
  return(x)


}

#' Select the dwt level of decomposition based on wavelet filter, data series length and a user choice
#'
#' @param choice
#'         Valid values:
#'              - "Conservative" : L < log2 ( N / (L - 1) + 1)
#'              - "max" : L <= log2(N)
#'              - "superMax" <= log2(1.5 * N)
#' @param f Wavelet transform filter name. To see the available filters use the function \code{\link{availableFilters}}
#' @param N Number of observations.
#'
#' @return Number of level of decomposition based in selection criteria
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
#'   Time Series Analysis. Cambridge: Cambridge University Press.
#'
#' @examples
#'  lev <- chooseLevel("conservative","haar",8)
#' @export
chooseLevel <- function (choice,f,N) {
  if (missing(choice)) stop("The argument \"choice\" must be provided. The available options are \"Conservative\", \"max\" and \"superMax\"")

  L <- wave.filter(f)[[1]]
  if (choice == "conservative") {
    J0 <- floor(log2( (N / (L - 1)) -1))
    return(J0 - 1)
  }
  if (choice == "max") {
    J0 <- floor(log2(N))
    return(J0 - 1)
  }
  if (choice == "supermax"){
    J0 <- floor(log2(1.5 * N))
    return(J0 - 1)
  }
  stop(paste(c("selected choice",choice,"its not valid. The available options are \"Conservative\", \"max\" and \"superMax\"")))
}

#' Extract observations from a WaveAnalysis
#'
#' This function permits to extract certain observations from a WaveAnalysis
#'
#' @param MWA WaveAnalysis from which the desired observations will be extracted
#' @param indices Indices that will indicate which observations will be extracted
#'
#' @return A list with two elements:
#'         * MWA: The WaveAnalysis provided minus the extracted observations.
#'         * MWAExtracted: A new WaveAnalysis with the extracted observations
#' @export
#'
#' @examples
#'           ECGExample <- loadECGExample()
#'           MWA <- MultiWaveAnalysis(ECGExample,"haar")
#'           aux <- extractSubset(MWA, c(1,2,3))
#'           MWATrain <- aux[[1]]
#'           MWATest <- aux[[2]]
extractSubset <- function(MWA,indices){
  if (missing(MWA)) stop("The argument \"MWA\" must be provided.")
  if (missing(indices)) stop("The argument \"indices\" must be provided.")

  n <- length(indices)

  MWA1 <- MWA
  MWA1$Observations <- n

  MWA2 <- MWA
  MWA2$Observations <- MWA2$Observations - n

  for (feature in names(MWA$Features)) {
    if ( !(is.na(MWA$Features[[feature]][1])) ){
        MWA1$Features[[feature]] <- MWA1$Features[[feature]][,indices, drop = FALSE]
        MWA2$Features[[feature]] <- MWA2$Features[[feature]][,-indices, drop = FALSE]
      }
    }
  return(list(MWA1,MWA2))
}

#' availableFeatures
#'
#' Print the available features for the \code{\link{MultiWaveAnalysis}} and \code{\link{StepDiscrim}}
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#'
#' @export
availableFeatures <- function() {
  writeLines("available features:
          - Variance (Var)
          - Correlation (Cor)
          - Interquartile range (IQR)
          - Permutation Entropy (PE)
          - Hoefflin s D measure (D)

          Example:
          c(\"Var\",\"PE\")
        ")
}

#' availableFilters
#'
#' Print the available filters for the wave analysis
#'
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#'
#' @export
availableFilters <- function(){
  writeLines("Avaible filters:
          - harr
          - d4
          - mb4
          - w4
          - bs3
          - fk4
          - d6
          - fk6
          - d8
          - fk8
          - la8
          - mb8
          - bl14
          - fk14
          - d16
          - la16
          - mb16
          - la20
          - bl20
          - fk22
          - mb24")
}

#' @importFrom stats IQR
computeIQR <- function(X) {
  aux <- head(X,-1)
  return (unlist(lapply(aux, function(x) stats::IQR(x))))
}

#' @importFrom utils head
#' @importFrom statcomp permutation_entropy ordinal_pattern_distribution
computePermutationEntropy <- function(X) {
  aux <- head(X,-1)
  return (unlist(lapply(aux, function(x) statcomp::permutation_entropy(statcomp::ordinal_pattern_distribution(x,6)))))
}

#' @importFrom utils head
#' @importFrom wdm wdm
computeDMeasure <- function(X,Y) {
  X <- head(X,-1)
  Y <- head(Y,-1)
  return (mapply(function(x,y) wdm::wdm(x,y,"hoeffding"),X,Y))
}



#' @export
print.WaveAnalysis <- function (x, ...) {
  print(paste("The number of correlation variables are",x$Cors))
  print(paste("The number of varianze variables are",x$Vars))
  print(paste("The Filter used is", x$filter))
  print(paste("The number of levels are", x$NLevels))
}

values <- function(MWA){
  stopifnot(class(MWA) == "WaveAnalysis")
  values <- matrix(0,nrow = 0,ncol = MWA$Observations)
  for (feature in MWA$Features) {
    if(!(is.na(feature[1]))) {
      values <- rbind(values,as.matrix(feature))
    }
  }
  return(values)
}

getAllFeatures <- function(){
  return (c("Var","Cor","IQR","PE","DM"))
}
