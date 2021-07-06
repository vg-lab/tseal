# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

Desing <- function (grp){
  nobs <- length(grp)
  index <- unique(grp)
  ngrps <- length(index)
  G <- matrix(0,nobs,ngrps)
  for (i in 1:nobs){
    for (g in 1:ngrps){
      if (grp[i] == index[g]){
        G[i,g] <- 1
      }
    }
  }
  return (G)
}


#' @importFrom MASS ginv
Lawley <- function (X,grps) {
  print(X)
  if (is.vector(X)){
    nobs<- length(X)
    p <- 1
    means <- mean(X)
  } else {
    nobs <- dim(X)[1]
    p <- dim(X)[2]
    means <- colMeans(X)
    print(X)
  }

  G <- Desing(grps)
  ngrps <- dim(G)[2]

  mean_W <- ginv(t(G) %*% G) %*% t(G) %*% X
  devs_W <- X - G %*% mean_W
  W <- t(devs_W) %*% devs_W

  devs_T <- X - matrix(1,nobs,1) %*% means
  T <- t(devs_T) %*% devs_T
  B <- T - W

  rcond_W <- rcond(W)

  if (rcond_W < 1e-8){
    V <- NaN
    Vp <- NaN
  } else {
    invW <- solve(W)
    V <- sum(diag(B %*% invW))
    Vp <- V/p
  }

  return (list(V,Vp))
}

StepDiscrim_ <- function (X,grps,maxvars,nCores) {
  #set parallel enviorement
  c <- makeCluster(nCores)
  mgrinit(c)
  makebarr(c)

  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- length(grps)

  mgrmakevar(c,"Xs",n,p)
  Xs[,] <- X[,]

  if (r == 1){
    grps <- t(grps)
  }

  vars_incl <- matrix(0,1,p)
  Vcum <- matrix(0,1,maxvars)
  Vpcum <- matrix(0,1,maxvars)

  for (step in 1:maxvars) {
    vi <- which(vars_incl>0)
    vni <- which(vars_incl==0)
    nvni <- length(vni)


    #V <- vector("numeric",nvni)
    #Vp <- vector("numeric",nvni)
    mgrmakevar(c,"Vs",nvni,1)
    mgrmakevar(c,"Vps",nvni,1)
    clusterExport(c,c("vi","vni","grps","nvni"), envir = environment())
    clusterExport(c,c("Lawley","ginv"))

    clusterEvalQ(c,{
      ids <-getidxs(nvni)
      for (v in ids) {
        aux <- Lawley(Xs[,c(vi,vni[v])],grps)
        Vs[v] <- aux[[1]]
        Vps[v] <- aux[[2]]
      }
      barr()
    })

    V <- as.vector(as.matrix(Vs))
    Vp<- as.vector(as.matrix(Vps))

    if (is.finite(sum(V))){
      i <- which.max(Vp)
      Vpmax <- Vp[i]
      vars_incl[vni[i]] <- step
      Vpcum[step] <- Vpmax
      Vcum[step] <- V[i]
    } else {
      break
    }
  }

  stoprdsm(c)

  aux <- sort(vars_incl,index.return = TRUE)
  y <- aux[[1]]
  incl <- aux[[2]]
  i <- which(y>0)
  incl <- incl[i]
  len_inc <- length(incl)
  Vcum <- Vcum[1:len_inc]
  Vpcum <-Vpcum[1:len_inc]
  return(list(incl,Vcum,Vpcum))
}

StepDiscrimV_ <- function(X,grps,VStep,nCores) {
  #set parallel enviorement
  c <- makeCluster(nCores)
  mgrinit(c)
  makebarr(c)

  n <- dim(X)[1]
  p <- dim(X)[2]
  r <- length(grps)

  mgrmakevar(X,"Xs",n,p)
  Xs[,] <- X[,]

  if (r == 1){
    grps <- t(grps)
  }

  vars_incl <- vector("numeric")
  vni <- 1:p
  Vcum <- vector("numeric")
  Vpcum <- vector("numeric")
  Vpmax = 100000
  while ( Vpmax > VStep){
    vi <- vars_incl
    nvni <- length(vni)

    mgrmakevar(c,"Vs",nvni,1)
    mgrmakevar(c,"Vps",nvni,1)
    clusterExport(c,c("vi","vni","grps","nvni"), envir = environment())
    clusterExport(c,c("Lawley","ginv"))

    clusterEvalQ(c,{
      ids <-getidxs(nvni)
      for (v in ids) {
        aux <- Lawley(Xs[,c(vi,vni[v])],grps)
        Vs[v] <- aux[[1]]
        Vps[v] <- aux[[2]]
      }
      barr()
    })

    V <- as.vector(as.matrix(Vs))
    Vp <- as.vector(as.matrix(Vps))

    if (is.finite(sum(V))){
      i <- which.max(Vp)
      Vpmax <- Vp[i]
      vars_incl <- append(vars_incl,vni[i])
      Vpcum <- append(Vpcum,Vpmax)
      Vcum <- append(Vcum,V[i])
      vni <- vni[-i]
    } else {
      break
    }
  }

  aux <- sort(vars_incl,index.return = TRUE)
  y <- aux[[1]]
  incl <- aux[[2]]
  i <- which(y>0)
  incl <- incl[i]
  len_inc <- length(incl)
  Vcum <- Vcum[1:len_inc]
  Vpcum <-Vpcum[1:len_inc]
  return(list(incl,Vcum,Vpcum))
}

#' Stepwise discriminant
#'
#' Stepwise discriminant analysis to determine the best subset of variables.
#' Introduces variables so as to maximize at each step the Lawley-Hotelling
#' trace (=Rao's V).  This measure is proportional to the mean Mahalanobis distance.
#'
#' Based on StepDiscrim of R.E. Strauss
#'
#' @param WVC WaveAnalisys object obtained with MultiVaweAnalisys function
#' @param grps labeled vector that classify the observations.
#' @param maxvars The number of desired values.
#' @param Var Determines if the algorithm take in account the variances.
#' For use this option, the provided WaveAnalisys has to have variances
#' @param Cor Determines if the algorithm take in account the correlations.
#'  For use this option, the provided WaveAnalisys has to have correlations
#' @param nCores determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return A WaveAnalisys object with the maxvars most discriminant variables,
#'  and a importance vector that determines the importance order of the selected variables.
#'
#' @export
#'
StepDiscrim <- function (WVC,grps,maxvars,Var = TRUE, Cor = TRUE, nCores = 0) {
  stopifnot(class(WVC) == "WaveAnalisys")
  NVar <- WVC$Vars
  NCor <- WVC$Cors

  if (Var & NVar<=0){
    simpleError("The provided analysis does not contain the variances")
  }

  if (Cor & NCor<=0){
    simpleError("The analysis provided does not contain correlations.")
  }

  if (nCores == 0) {
    nCores <- detectCores() - 1
  }

  rv <- WVC[[1]]
  if (Var & Cor) {
    Tr <- rv[1:(NVar + NCor),]
  } else if (Var) {
    Tr <- rv[1:NVar,]
  } else if (Cor){
    if (NVar>0){
      Tr <- rv[(NVar + 1):(NVar + NCor),]
      NVar <- 0
    } else {
      Tr <- rv[1:NCor,]
    }
  }

  incl <- StepDiscrim_(t(Tr),grps,maxvars,nCores)[[1]]
  inclSorted <- sort(incl, index.return = TRUE)
  WVC[[1]] <- Tr[inclSorted$x,]
  WVC$Vars <- length(which(incl<NVar))
  WVC$Cors <- length(incl) - WVC$Vars

  importance <- rep(0,length(incl))
  for (i in 1:length(incl)) {
    importance[inclSorted$ix[i]] =  i
  }


  WVC$importance <- importance
  return(WVC)
}

#' Stepwise discriminant
#'
#' Stepwise discriminant analysis to determine the best subset of variables.
#' Introduces variables so as to maximize at each step the Lawley-Hotelling
#' trace (=Rao's V).  This measure is proportional to the mean Mahalanobis distance.
#'
#' Based on StepDiscrim of R.E. Strauss
#'
#' @param WVC WaveAnalisys object obtained with MultiVaweAnalisys function
#' @param grps labeled vector that classify the observations-
#' @param VStep Determine the minimum value of V to continue adding new variables. Ex if an determinate step the maximum V is 0.2 but VStep is 0.3 the algorithm end.
#' @param Var Determines if the algorithm take in account the variances. For use this option, the provided WaveAnalisys has to have variances
#' @param Cor Determines if the algorithm take in account the correlations. For use this option, the provided WaveAnalisys has to have correlations
#' @param nCores determines the number of processes that will be used in the function, by default it uses all but one of the system cores.
#'
#' @return A WaveAnalisys object with the maxvars most discriminant variables, and a importance vector that determines the importance order of the selected variables.
#'
#' @export
#'
StepDiscrimV <- function (WVC,grps,VStep,Var = TRUE, Cor = TRUE, nCores = 0) {
  stopifnot(class(WVC) == "WaveAnalisys")
  NVar <- WVC$Vars
  NCor <- WVC$Cors

  # TODO TESTING
  if (Var & NVar<=0){
    simpleError("The provided analysis does not contain the variances")
  }

  if (Cor & NCor<=0){
    simpleError("The analysis provided does not contain correlations.")
  }

  if (nCores == 0) {
    nCores <- detectCores() -1
  }

  rv <- WVC[[1]]
  if (Var & Cor) {
    Tr <- rv[1:(NVar + NCor),]
  } else if (Var) {
    Tr <- rv[1:NVar,]
  } else if (Cor){
    if (NVar>0){
      Tr <- rv[NVar:(NVar + NCor),]
    } else {
      Tr <- rv[1:NCor,]
    }
  }

  incl <- StepDiscrimV_(t(Tr),grps,VStep,nCores)[[1]]
  inclSorted <- sort(incl, index.return = TRUE)
  WVC[[1]] <- Tr[inclSorted$x,]
  WVC$Vars <- length(which(incl<NVar))
  WVC$Cors <- length(incl) - WVC$Vars

  importance <- rep(0,length(incl))
  for (i in 1:length(incl)) {
    importance[inclSorted$ix[i]] =  i
  }

  WVC$importace <- importance
  return(WVC)
}
