#' # Title     : TODO
#' # Objective : TODO
#' # Created by: ivan
#' # Created on: 8/4/21
#'
#' #' Generate a MultiWave analisys
#' #'
#' #' @param XSeries1 Sample from the population 1 ()
#' #' @param XSeries2 Sample from the population 2 ()
#' #' @param f Selected wavelet filter for the analisys. To view all the avaivle filters type ...
#' #' @param lev Wavelet descomposition level by default is selcted using the "conservative" strategy. See chooseLevel function.
#' #' @param Var a boolean that determines if the analissy contains the variances.
#' #' @param Cor a boolean that determines if the analissy contains the correlations.
#' #'
#' #' @return An object of class WaveAnalisys
#' #'
#' #' @export
#' #'
#' #' @examples
#' # MultiVaweAnalisys <- function (XSeries1,XSeries2,f,lev = 0,Var = TRUE, Cor = TRUE){
#' #   dim1 <- dim(XSeries1)
#' #   dim2 <- dim(XSeries2)
#' #
#' #   nv1 <- dim1[1]  #Time series Dimension
#' #   nr1 <- dim1[2]  #Time series Lenght
#' #   nc1 <- dim1[3]  #Sample Size
#' #
#' #   nv2 <- dim2[1]
#' #   nr2 <- dim2[2]
#' #   nc2 <- dim2[3]
#' #
#' #   if (nv1 != nv2){
#' #     simpleWarning("Dimension Mismatch")
#' #     nv1 <- min(nv1,nv2)
#' #     nv2 <- nv1
#' #   }
#' #
#' #   if (nr1 != nr2) {
#' #     simpleWarning("Length Mismatch")
#' #     nr1 <- min(nr1,nr2)
#' #     nr2 <- nr1
#' #   }
#' #
#' #   if (lev == 0) {
#' #     lev <- chooseLevel("conservative",f,nr1)
#' #   }
#' #
#' #
#' #   # Standarizing the data
#' #   YSeries1 <- XSeries1[1:nv1,1:nr1,1:nc1]
#' #   YSeries2 <- XSeries2[1:nv1,1:nr1,1:nc2]
#' #   for (i in 1:nv1) {
#' #     for (j in 1:nc1) {
#' #       Vtemporal <- YSeries1[i,,j]
#' #       YSeries1[i,,j] <- (Vtemporal - mean(Vtemporal)) / sd(Vtemporal)
#' #     }
#' #     for (j in 1:nc2) {
#' #       Vtemporal <- YSeries2[i,,j]
#' #       YSeries2[i,,j] <- (Vtemporal - mean(Vtemporal)) / sd(Vtemporal)
#' #     }
#' #   }
#' #
#' #   NVar <- if(Var) nv1 * lev else 0
#' #   if (Cor) {
#' #     NbK <- combn(1:nv1,2)
#' #     NumberNbK <- dim(NbK)[2]
#' #     NCor <- NumberNbK * lev
#' #   } else {
#' #     NCor <- 0
#' #   }
#' #
#' #   WVC <- matrix(0,NVar + NCor,nc1 + nc2 )
#' #
#' #   for (i in 1:nc1){
#' #     WJ <- list()
#' #       for (j in 1:nv1){
#' #         Vtemporal <- YSeries1[j,,i]
#' #         aux.modwt <- modwt(Vtemporal,f,lev)
#' #         WJ <- append(WJ,list(aux.modwt))
#' #         if (Var) {
#' #           WVARAux <- wave.variance(aux.modwt)
#' #           WVC[((j-1) * lev + 1):(j*lev),i] <- WVARAux[1:lev,1]
#' #         }
#' #       }
#' #
#' #     if (Cor) {
#' #       for (k in 1:NumberNbK) {
#' #         WCOR <- suppressWarnings(wave.correlation(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]],nr1))
#' #         WVC[(NVar+(k-1)*lev + 1):(NVar + k * lev),i] <- WCOR[1:lev,1]
#' #       }
#' #     }
#' #   }
#' #
#' #   for (i in (nc1 + 1):(nc1 + nc2)) {
#' #     WJ <- list()
#' #       for (j in 1:nv1){
#' #         Vtemporal <- YSeries2[j,,i-nc1]
#' #         aux.modwt <- modwt(Vtemporal,f,lev)
#' #         WJ <- append(WJ,list(aux.modwt))
#' #         if (Var) {
#' #           WVARAux <- wave.variance(aux.modwt)
#' #           WVC[((j-1) * lev + 1):(j*lev),i] <- WVARAux[1:lev,1]
#' #         }
#' #       }
#' #
#' #     if (Cor) {
#' #       for (k in 1:NumberNbK) {
#' #         WCOR <- suppressWarnings(wave.correlation(WJ[[NbK[1,k]]],WJ[[NbK[2,k]]],nr1))
#' #         WVC[(NVar+(k-1)*lev + 1):(NVar + k * lev),i] <- WCOR[1:lev,1]
#' #       }
#' #     }
#' #   }
#' #
#' #   x <- list(Values = WVC,Observations = nc1+nc2, Vars = NVar, Cors = NCor, NLevels = lev, filter = f, importance = vector("numeric"))
#' #   attr(x,"class") <- "WaveAnalisys"
#' #   return(x)
#' # }
#'
#' #' Select the dwt level of decomposition based on wavelet filter, data series length and a user choice
#' #'
#' #' @param choice
#' #'         Valid values:
#' #'              - "Conservative" : L < log2 ( N / (L - 1) + 1)
#' #'              - "max" : L <= log2(N)
#' #'              - "superMax" <= log2(1.5 * N)
#' #' @param f Wavelet transform filter name. For avaible filters use ...
#' #' @param N Number of observations.
#' #'
#' #' @return Number of level of descomposition based in selection criteria
#' #' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for
#' #'   Time Series Analysis. Cambridge: Cambridge University Press.
#' #'
#' #' @examples
#' #'  ChooseLevel("conservative","haar",8)
#' chooseLevel <- function (choice,f,N) {
#'   L <- wave.filter(f)[[1]]
#'   if (choice == "conservative") {
#'     J0 <- floor(log2( (N / (L - 1)) -1))
#'     return(J0 - 1)
#'   }
#'   if (choice == "max") {
#'     J0 <- floor(log2(N))
#'     return(J0 - 1)
#'   }
#'   if (choice == "supermax"){
#'     J0 <- floor(log2(1.5 * N))
#'     return(J0 - 1)
#'   }
#'   simpleError("selected choice its not valid")
#' }
#'
#' extractSubset <- function(WMA,x){
#'   n <- length(x)
#'   WMA1 <- WMA
#'   WMA1$Values <- WMA1$Values[,x]
#'   WMA1$Observations <- n
#'
#'   WMA2 <- WMA
#'   WMA2$Values <- WMA2$Values[,-x]
#'   WMA2$Observations <- WMA2$Observations - n
#'   return(list(WMA1,WMA2))
#'
#' }
#'
#' avaibleFilters <- function(){
#'   print("Avaible filters:
#'           harr
#'           d4
#'           mb4
#'           w4
#'           bs3
#'           fk4
#'           d6
#'           fk6
#'           d8
#'           fk8
#'           la8
#'           mb8
#'           bl14
#'           fk14
#'           d16
#'           la16
#'           mb16
#'           la20
#'           bl20
#'           fk22
#'           mb24")
#' }
#'
#' print.WaveAnalisys <- function (x) {
#'   print(paste("The number of correlation variables are",x$Cors))
#'   print(paste("The number of varianze variables are",x$Vars))
#'   print(paste("The Filter used is", x$filter))
#'   print(paste("The number of levels are", x$NLevels))
#' }
#'


