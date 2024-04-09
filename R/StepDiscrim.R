# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

Desing <- function(grp) {
    nobs <- length(grp)
    index <- unique(grp)
    ngrps <- length(index)
    G <- matrix(0, nobs, ngrps)
    for (i in seq_len(nobs)) {
        for (g in seq_len(ngrps)) {
            if (grp[i] == index[g]) {
                G[i, g] <- 1
            }
        }
    }
    return(G)
}


#' @importFrom MASS ginv
#' @noRd
Lawley <- function(X, grps) {
    if (is.vector(X)) {
        nobs <- length(X)
        p <- 1
        means <- mean(X)
    } else {
        nobs <- dim(X)[1]
        p <- dim(X)[2]
        means <- colMeans(X)
    }

    G <- Desing(grps)
    ngrps <- dim(G)[2]

    mean_W <- MASS::ginv(t(G) %*% G) %*% t(G) %*% X
    devs_W <- X - G %*% mean_W
    W <- t(devs_W) %*% devs_W

    devs_T <- X - matrix(1, nobs, 1) %*% means
    TMatrix <- t(devs_T) %*% devs_T
    B <- TMatrix - W

    rcond_W <- rcond(W)

    if (rcond_W < 1e-8) {
        V <- 0
        Vp <- 0
    } else {
        invW <- solve(W)
        V <- sum(diag(B %*% invW))
        Vp <- V / p
    }

    return(list(V, Vp))
}

geTSEALectedFeatures <- function(MWA, features) {
    features <- tolower(features)
    features <- unique(features)
    MWAFeatures <- tolower(names(MWA$Features))
    idx <-
        vapply(features, function(x) {
            match(x, MWAFeatures)
        }, FUN.VALUE = 0L)

    Tr <- matrix(0, nrow = 0, ncol = MWA$Observations)
    for (ix in seq_along(idx)) {
        id <- idx[ix]
        if (is.na(id) || all(is.na(MWA$Features[[id]]))) {
            stop("The provided analysis does not contain",
                 as.character(features[ix]))
        } else {
            Tr <- rbind(Tr, MWA$Features[[id]])
        }
    }
    return(list(Tr, idx))
}

GenerateMWADiscrim <- function(MWA, incl, idx) {
    MWAAux <- list(
        Features = list(
            Var = NA,
            Cor = NA,
            IQR = NA,
            DM = NA,
            PE = NA
        ),
        StepSelection = list(
            Var = NA,
            Cor = NA,
            IQR = NA,
            DM = NA,
            PE = NA
        ),
        Observations = MWA$Observations,
        NLevels = MWA$NLevels,
        filter = MWA$filter
    )
    attr(MWAAux, "class") <- "WaveAnalysis"
    acc <- 1
    for (id in idx) {
        size <- dim(MWA$Features[[id]])[1]
        upperLimit <- acc + size - 1
        aux <- which(incl >= acc & incl <= upperLimit)
        index <-
            vapply(incl[aux], function(x) {
                x - (acc - 1)
            }, FUN.VALUE = 0)

        if (length(index) > 0) {
            MWAAux$StepSelection[[id]] <- index
            if (length(index) == 1) {
                MWAAux$Features[[id]] <- t(as.matrix(MWA$Features[[id]][index, ]))
            } else {
                MWAAux$Features[[id]] <- MWA$Features[[id]][index, ]
            }
        }
        acc <- upperLimit + 1
    }

    return(MWAAux)
}

#' @importFrom methods is
#' @noRd
StepDiscrimRaw_ <-
    function(MWA, grps, maxvars, features, nCores = 0) {
        checkmate::anyMissing(c(MWA, grps, maxvars, features))

        if (length(features) == 0) {
            stop(
                "At least one feature must be provided. To see the available
                features use availableFeatures()"
            )
        }

        if (length(grps) != MWA$Observations) {
            stop("The \"grps\" length mismatches with the observations of \"MWA\"")
        }

        maxvars <- checkmate::asCount(maxvars)


        stopifnot(is(MWA, "WaveAnalysis"))

        Tr <- geTSEALectedFeatures(MWA, features)[[1]]

        if (nCores == 0) {
            nCores <- parallelly::availableCores(omit = 1)
        }

        incl <- StepDiscrim_(t(Tr), grps, maxvars, nCores)[[1]]
        return(list(Tr, incl))
    }

StepDiscrim_ <- function(X, grps, maxvars, nCores) {
    # maxVars to maximum numbers of vars available
    maxvars <- min(maxvars, dim(X)[2])
    # set parallel enviorement
    c <- makeCluster(nCores)
    mgrinit(c)
    tryCatch({
        n <- dim(X)[1]
        p <- dim(X)[2]
        r <- length(grps)

        mgrmakevar(c, "Xs", n, p)
        Xs[, ] <- X[, ]

        if (r == 1) {
            grps <- t(grps)
        }

        vars_incl <- matrix(0, 1, p)
        Vcum <- matrix(0, 1, maxvars)
        Vpcum <- matrix(0, 1, maxvars)

        for (step in seq_len(maxvars)) {
            vi <- which(vars_incl > 0)
            vni <- which(vars_incl == 0)
            nvni <- length(vni)


            # V <- vector("numeric",nvni)
            # Vp <- vector("numeric",nvni)
            mgrmakevar(c, "Vs", nvni, 1)
            mgrmakevar(c, "Vps", nvni, 1)
            clusterExport(c, c("vi", "vni", "grps", "nvni"), envir = environment())
            clusterExport(c, c("Lawley", "Desing"), envir = loadNamespace("TSEAL"))

            clusterEvalQ(c, {
                ids <- getidxs(nvni)
                if (length(ids) != 0) {
                    for (v in ids) {
                        aux <- Lawley(Xs[, c(vi, vni[v])], grps)
                        Vs[v] <- aux[[1]]
                        Vps[v] <- aux[[2]]
                    }
                }
            })

            V <- as.vector(as.matrix(Vs))
            Vp <- as.vector(as.matrix(Vps))

            if (is.finite(sum(V))) {
                i <- which.max(Vp)
                Vpmax <- Vp[i]
                vars_incl[vni[i]] <- step
                Vpcum[step] <- Vpmax
                Vcum[step] <- V[i]
            } else {
                break
            }
        }
    },
    finally = stoprdsm(c))

    aux <- sort(vars_incl, index.return = TRUE)
    y <- aux[[1]]
    incl <- aux[[2]]
    i <- which(y > 0)
    incl <- incl[i]
    len_inc <- length(incl)
    Vcum <- Vcum[seq_len(len_inc)]
    Vpcum <- Vpcum[seq_len(len_inc)]
    return(list(incl, Vcum, Vpcum))
}

StepDiscrimV_ <- function(X, grps, VStep, nCores) {
    # set parallel enviorement
    c <- makeCluster(nCores)
    mgrinit(c)
    tryCatch({
        n <- dim(X)[1]
        p <- dim(X)[2]
        r <- length(grps)

        mgrmakevar(c, "Xs", n, p)
        Xs[, ] <- X[, ]

        if (r == 1) {
            grps <- t(grps)
        }

        vars_incl <- vector("numeric")
        vni <- seq_len(p)
        nvni <- length(vni)
        Vcum <- vector("numeric")
        Vpcum <- vector("numeric")
        Vpmax <- 100000
        while (Vpmax > VStep && nvni > 0) {
            vi <- vars_incl
            nvni <- length(vni)

            mgrmakevar(c, "Vs", nvni, 1)
            mgrmakevar(c, "Vps", nvni, 1)
            clusterExport(c, c("vi", "vni", "grps", "nvni"), envir = environment())
            clusterExport(c, c("Lawley"), envir = loadNamespace("TSEAL"))

            clusterEvalQ(c, {
                ids <- getidxs(nvni)
                if (length(ids) != 0) {
                    for (v in ids) {
                        aux <- Lawley(Xs[, c(vi, vni[v])], grps)
                        Vs[v] <- aux[[1]]
                        Vps[v] <- aux[[2]]
                    }
                }
            })

            V <- as.vector(as.matrix(Vs))
            Vp <- as.vector(as.matrix(Vps))

            if (is.finite(sum(V))) {
                i <- which.max(Vp)
                Vpmax <- Vp[i]
                vars_incl <- append(vars_incl, vni[i])
                Vpcum <- append(Vpcum, Vpmax)
                Vcum <- append(Vcum, V[i])
                vni <- vni[-i]
            } else {
                break
            }
        }
    },
    finally = stoprdsm(c))

    return(list(vars_incl, Vcum, Vpcum))
}

#' Select the most discriminating variables
#'
#' Stepwise discriminant analysis to determine the best subset of variables.
#' Introduces variables so as to maximize at each step the Lawley-Hotelling
#' trace (=Rao's V).  This measure is proportional to the mean Mahalanobis
#'  distance.
#'
#' Based on StepDiscrim of R.E. Strauss
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalysis function
#' @param grps Labeled vector that classify the observations.
#' @param maxvars The number of desired values. Must be a positive integer
#' @param features A list of characteristics that will be used for the
#'        classification process. To see the available features see
#'        \code{\link{availableFeatures}}
#' @param nCores Determines the number of processes that will be used in the
#'        function, by default it uses all but one of the system cores. Must be
#'        a positive integer, where 0 corresponds to the default behavior
#'
#' @return A WaveAnalysis object with the maxvars most discriminant variables
#'
#' @examples
#' ECGExample <- loadECGExample()
#' MWA <- MultiWaveAnalysis(ECGExample, "haar")
#' MWADiscrim <- StepDiscrim(
#'   MWA, c(rep(1, 5), rep(2, 5)), 20,
#'   c("Var", "Cor", "IQR")
#' )
#'
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrimV}}
#' @export
#'
StepDiscrim <- function(MWA,
                        grps,
                        maxvars,
                        features = c("Var", "Cor", "IQR", "PE", "DM"),
                        nCores = 0) {
    checkmate::anyMissing(c(MWA, grps, maxvars))


    if (length(features) == 0) {
        stop(
            "At least one feature must be provided. To see the available filters
         use the availableFeatures()"
        )
    }

    if (length(grps) != MWA$Observations) {
        stop("The \"grps\" length mismatches with the observations of \"MWA\"")
    }

    maxvars <- checkmate::asCount(maxvars)
    nCores <- checkmate::asCount(nCores)


    aux <- geTSEALectedFeatures(MWA, features)
    Tr <- aux[[1]]
    idx <- aux[[2]]

    if (nCores == 0) {
        nCores <- parallelly::availableCores(omit = 1)
    }


    incl <- StepDiscrim_(t(Tr), grps, maxvars, nCores)[[1]]

    return(GenerateMWADiscrim(MWA, incl, idx))
}

#' Select the most discriminating variables
#'
#' Stepwise discriminant analysis to determine the best subset of variables.
#' Introduces variables so as to maximize at each step the Lawley-Hotelling
#' trace (=Rao's V).  This measure is proportional to the mean Mahalanobis
#' distance.The process ends when in one step the value of the Lawley-Hotelling
#' trace is less than a given value.
#'
#' Based on StepDiscrim of R.E. Strauss
#'
#' @param MWA WaveAnalysis object obtained with MultiWaveAnalysis function
#' @param grps Labeled vector that classify the observations.
#' @param VStep Determine the minimum value of V to continue adding new
#'        variables. Ex if an determinate step the maximum V is 0.2 but VStep is
#'         0.3 the algorithm end. Must be greater than 0.
#' @param features A list of characteristics that will be used for the
#'        classification process. To see the available features see
#'        \code{\link{availableFeatures}}
#' @param nCores Determines the number of processes that will be used in the
#'        function, by default it uses all but one of the system cores. Must be
#'        a positive integer, where 0 corresponds to the default behavior
#'
#' @return A WaveAnalysis object with the most discriminant variables
#'
#' @examples
#' ECGExample <- loadECGExample()
#' MWA <- MultiWaveAnalysis(ECGExample, "haar")
#' MWADiscrim <- StepDiscrimV(
#'   MWA, c(rep(1, 5), rep(2, 5)), 0.2,
#'   c("Var", "Cor", "IQR")
#' )
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrim}}
#'
#' @export
#'
StepDiscrimV <- function(MWA,
                         grps,
                         VStep,
                         features = c("Var", "Cor", "IQR", "PE", "DM"),
                         nCores = 0) {
    checkmate::anyMissing(c(MWA, grps, VStep))

    if (VStep <= 0) {
        stop("The argument \"VStep\" must be provided and must be grater than 0")
    }

    if (length(features) == 0) {
        stop(
            "At least one feature must be provided. To see the available filters
         use the availableFeatures()"
        )
    }

    if (length(grps) != MWA$Observations) {
        stop("The \"grps\" length mismatches with the observations of \"MWA\"")
    }

    if (!is.numeric(VStep) || length(VStep) != 1 || VStep <= 0) {
        stop("The argument \"VStep\" must be a number greater than 0")
    }

    nCores <- checkmate::asCount(nCores)


    aux <- geTSEALectedFeatures(MWA, features)
    Tr <- aux[[1]]
    idx <- aux[[2]]

    if (nCores == 0) {
        nCores <- parallelly::availableCores(omit = 1)
    }


    incl <- StepDiscrimV_(t(Tr), grps, VStep, nCores)[[1]]

    return(GenerateMWADiscrim(MWA, incl, idx))
}

#' Allows to select the same variables for a given StepDiscrim
#'
#' Allows to perform the same variable selection in a new MWA object starting
#' from a MWA object with the variables already selected (it is advisable that
#' the parameters of the MWA and of the selection are the same).
#'
#' @param MWA WaveAnalysis object on which variables are to be selected.
#' @param MWADiscrim WaveAnalysis object on which certain variables have been
#'  previously selected, using \code{\link{StepDiscrim}} or
#'   \code{\link{StepDiscrimV}}
#'
#' @return An object of class WaveAnalysis with the same variables selected as in the
#'         MWADiscrim object.
#' @export
#'
#' @examples
#' ECGExample <- loadECGExample()
#' # We simulate that the second series has been obtained after
#' Series1 <- ECGExample[, , 1:9]
#' Series2 <- ECGExample[, , 10, drop = FALSE]
#' MWA <- MultiWaveAnalysis(Series1, "haar", features = c("Cor", "PE"))
#' MWADiscrim <- StepDiscrimV(MWA, c(rep(1, 5), rep(2, 4)), 0.2,
#'   features = c("Cor", "PE")
#' )
#'
#' MWA2 <- MultiWaveAnalysis(Series2, "haar", features = c("Cor", "PE"))
#' MWA2Discrim <- SameDiscrim(MWA2, MWADiscrim)
#' # At this point MWA2Discrim has the same variables that MWADiscrim
#' # and can be used in a pretrained model with MWADiscrim
#'
#' @seealso
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
SameDiscrim <- function(MWA, MWADiscrim) {

    checkmate::anyMissing(c(MWA, MWADiscrim))

    if (MWA$NLevels != MWADiscrim$NLevels) {
        stop("The number of levels of descomposition must be the same")
    }

    if (all(is.na(MWADiscrim$StepSelection))) {
        stop(
            "The \"MWADiscrimination\" provided has no variables selected, probably
          because they have not been selected. See MWADiscrim function for more
         information."
        )
    }

    MWAAux <- list(
        Features = list(
            Var = NA,
            Cor = NA,
            IQR = NA,
            DM = NA,
            PE = NA
        ),
        StepSelection = list(
            Var = NA,
            Cor = NA,
            IQR = NA,
            DM = NA,
            PE = NA
        ),
        Observations = MWA$Observations,
        NLevels = MWA$NLevels,
        filter = MWA$filter
    )
    attr(MWAAux, "class") <- "WaveAnalysis"

    for (feature in names(MWA$Features)) {
        if (!all(is.na(MWADiscrim$StepSelection[[feature]]))) {
            selection <- MWADiscrim$StepSelection[[feature]]
            MWAAux$Features[[feature]] <-
                MWA$Features[[feature]][selection, ,
                                        drop = FALSE]
            MWAAux$StepSelection[[feature]] <- selection
        }
    }

    return(MWAAux)
}
