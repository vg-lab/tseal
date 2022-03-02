# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' Generate a MultiWave analisys
#'
#' @param XSeries1 Sample from the population 1 (dim x length x cases)
#' @param XSeries2 Sample from the population 2 (dim x length x cases)
#' @param f Selected wavelet filter for the analisys. To view all the avaivle filters type ...
#' @param lev Wavelet descomposition level by default is selcted using the "conservative" strategy. See chooseLevel function.
#' @param Var a boolean that determines if the analissy contains the variances.
#' @param Cor a boolean that determines if the analissy contains the correlations.
#' @param nCores determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return An object of class WaveAnalisys
#'
#' @export
#'
#' @importFrom bigmemory as.matrix GetMatrixSize
#' @importFrom waveslim modwt wave.variance wave.correlation wave.filter
#' @importFrom parallel makeCluster clusterExport clusterEvalQ detectCores
#' @importFrom pryr object_size
#' @importFrom stats sd
#' @importFrom utils combn
#' @importFrom statcomp permutation_entropy ordinal_pattern_distribution
#' @importFrom wdm wdm
MultiVaweAnalisys <- function (XSeries1,XSeries2,f,lev = 0,features = c("Var","Cor","IQR","PE","DM"), nCores = 0){
  if (nCores == 0) {
    nCores = detectCores() - 1
  }
  HVar <- "Var" %in% features
  HCor <- "Cor" %in% features
  HIQR <- "IQR" %in% features
  HPE  <- "PE"  %in% features
  HDM  <- "DM"  %in% features

  dim1 <- dim(XSeries1)
  dim2 <- dim(XSeries2)

  nv1 <- dim1[1]  #Time series Dimension
  nr1 <- dim1[2]  #Time series Lenght
  nc1 <- dim1[3]  #Sample Size

  nv2 <- dim2[1]
  nr2 <- dim2[2]
  nc2 <- dim2[3]

  if (nv1 != nv2){
    simpleWarning("Dimension Mismatch")
    nv1 <- min(nv1,nv2)
    nv2 <- nv1
  }

  if (nr1 != nr2) {
    simpleWarning("Length Mismatch")
    nr1 <- min(nr1,nr2)
    nr2 <- nr1
  }

  if (lev == 0) {
    lev <- chooseLevel("conservative",f,nr1)
  }

  #set parallel enviorement
  c <- makeCluster(nCores)
  mgrinit(c)
  makebarr(c)


  # Standarizing the data
  mgrmakevar(c,"YSeries1",nc1,nv1*nr1)
  mgrmakevar(c,"YSeries2",nc2,nv1*nr1)

  for (i in 1:nc1) {
    YSeries1[i,] <- as.vector(t(XSeries1[,,i]))
  }

  for (i in 1:nc2) {
    YSeries2[i,] <- as.vector(t(XSeries2[,,i]))
  }


  clusterExport(c,c("nv1","nr1","nc1","nc2"), envir= environment())
  clusterExport(c,c("modwt","wave.variance","wave.correlation","D3toD2","computeIQR",
                    "computePermutationEntropy","permutation_entropy", "wdm",
                    "ordinal_pattern_distribution","computeDMeasure"))
  clusterEvalQ(c,{
  ids <- getidxs(nv1)
  for (i in ids) {
    for (j in 1:nc1) {
      aux <- D3toD2(i,,j,nv1,nr1,nc1)
      Vtemporal <- YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
      YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]] <- (Vtemporal - mean(Vtemporal)) / sd(Vtemporal)
    }
    for (j in 1:nc2) {
      aux <- D3toD2(i,,j,nv1,nr1,nc2)
      Vtemporal <- YSeries2[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
      YSeries2[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]] <- (Vtemporal - mean(Vtemporal)) / sd(Vtemporal)
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
    mgrmakevar(c,"Var",NVar,nc1 + nc2)
  } else {
    Var = NA
  }
  if (HCor) {
    mgrmakevar(c,"Cor",NCor,nc1 + nc2)
  } else {
    Cor = NA
  }
  if (HIQR) {
    mgrmakevar(c,"IQR",NVar,nc1 + nc2)
  } else {
    IQR = NA
  }
  if (HPE) {
    mgrmakevar(c,"PE",NVar,nc1 + nc2)
  } else {
    PE = NA
  }
  if (HDM) {
    mgrmakevar(c,"DM",NCor,nc1 + nc2)
  } else {
    DM = NA
  }

  clusterEvalQ (c, {
    ids <- getidxs(nc1)
    for (i in ids){
      WJ <- list()
      for (j in 1:nv1){
        aux = D3toD2(j,,i,nv1,nr1,nc1)
        Vtemporal <- YSeries1[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
        aux.modwt <- modwt(Vtemporal,f,lev)
        WJ <- append(WJ,list(aux.modwt))

        if (HVar) {
          WVARAux <- wave.variance(aux.modwt)
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

      if (HCor || HDM) {
        for (k in 1:NumberNbK) {
          if (HCor){
            e1 <- NbK[1,k]
            e2 <- NbK[2,k]
              WCOR <- suppressWarnings(wave.correlation(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]],nr1))
              Cor[((k-1)*lev + 1):(k * lev),i] <- WCOR[1:lev,1]
            }
          if (HDM) {
            WDM <- computeDMeasure(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]])
            DM[((k-1)*lev + 1):(k * lev),i] <- WDM
          }
        }
      }
    }

    ids <- getidxs(nc2)
    ids <- sapply(ids, function(x) x+nc1)
    for (i in ids) {
      WJ <- list()
      for (j in 1:nv1){
        aux = D3toD2(j,,i-nc1,nv1,nr1,nc2)
        Vtemporal <- YSeries2[aux[[1]],aux[[2]][[1]]:aux[[2]][[2]]]
        aux.modwt <- modwt(Vtemporal,f,lev)
        WJ <- append(WJ,list(aux.modwt))


        if (HVar) {
          WVARAux <- wave.variance(aux.modwt)
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

      if (HCor || HDM) {
        for (k in 1:NumberNbK) {
          if (HCor){
            e1 <- NbK[1,k]
            e2 <- NbK[2,k]
            WCOR <- suppressWarnings(wave.correlation(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]],nr1))
            Cor[((k-1)*lev + 1):(k * lev),i] <- WCOR[1:lev,1]
          }
          if (HDM) {
            WDM <- computeDMeasure(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]])
            DM[((k-1)*lev + 1):(k * lev),i] <- WDM
          }
        }
      }
    }
    barr()
  })

  stoprdsm(c)

  x <- list(Features = list(Var = as.matrix(Var), Cor = as.matrix(Cor), IQR = as.matrix(IQR), DM = as.matrix(DM), PE = as.matrix(PE)), Observations = nc1+nc2, NLevels = lev, filter = f)
  attr(x,"class") <- "WaveAnalisys"
  return(x)
}

#' Select the dwt level of decomposition based on wavelet filter, data series length and a user choice
#'
#' @param choice
#'         Valid values:
#'              - "Conservative" : L < log2 ( N / (L - 1) + 1)
#'              - "max" : L <= log2(N)
#'              - "superMax" <= log2(1.5 * N)
#' @param f Wavelet transform filter name. For avaible filters use ...
#' @param N Number of observations.
#'
#' @return Number of level of descomposition based in selection criteria
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
#'   Time Series Analysis. Cambridge: Cambridge University Press.
#'
#' @examples
#'  lev <- chooseLevel("conservative","haar",8)
#' @export
chooseLevel <- function (choice,f,N) {
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
  stop(paste(c("selected choice",choice,"its not valid")))
}

extractSubset <- function(MWA,x){
  n <- length(x)

  MWA1 <- MWA
  MWA1$Observations <- n

  MWA2 <- MWA
  MWA2$Observations <- MWA2$Observations - n

  for (feature in names(MWA$Features)) {
    if ( !(is.na(MWA$Features[[feature]][1])) ){
      if ( dim(MWA$Features[[feature]])[1] == 1){
        MWA1$Features[[feature]] <- t(as.matrix(MWA1$Features[[feature]][,x]))
        MWA2$Features[[feature]] <- t(as.matrix(MWA2$Features[[feature]][,-x]))
      } else {
        MWA1$Features[[feature]] <- MWA1$Features[[feature]][,x]
        MWA2$Features[[feature]] <- MWA2$Features[[feature]][,-x]
      }
    }
  }

# -------------------------------------------------------------------------



  return(list(MWA1,MWA2))

}

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

computeIQR <- function(X) {
  aux <- head(X,-1)
  return (unlist(lapply(aux, function(x) IQR(x))))
}

computePermutationEntropy <- function(X) {
  aux <- head(X,-1)
  return (unlist(lapply(aux, function(x) permutation_entropy(ordinal_pattern_distribution(x,6)))))
}

computeDMeasure <- function(X,Y) {
  X <- head(X,-1)
  Y <- head(Y,-1)
  return (mapply(function(x,y) wdm(x,y,"hoeffding"),X,Y))
}



#' @export
print.WaveAnalisys <- function (x) {
  print(paste("The number of correlation variables are",x$Cors))
  print(paste("The number of varianze variables are",x$Vars))
  print(paste("The Filter used is", x$filter))
  print(paste("The number of levels are", x$NLevels))
}

#' @export
values <- function(MWA){
  stopifnot(class(MWA) == "WaveAnalisys")
  values <- matrix(0,nrow = 0,ncol = MWA$Observations)
  for (feature in MWA$Features) {
    if(!(is.na(feature[1]))) {
      values <- rbind(values,as.matrix(feature))
    }
  }
  return(values)
}

#' @export
getAllFeatures <- function(){
  return (c("Var","Cor","IQR","PE","DM"))
}
