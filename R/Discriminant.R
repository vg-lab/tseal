# Title     : TODO
# Objective : TODO
# Created by: ivan
# Created on: 8/4/21

#' Applies a discriminant to perform a classification
#'
#' This function performs a classification using a pretrained discriminant
#' (linear or quadratic), and return a confusion matrix with the results
#' obtained
#'
#' @param model Previously trained model (linear or quadratic) that will be
#'        used for classification.
#' @param ... Additional arguments
#'
#'
#' @seealso
#' * \code{\link{testModel.lda}}
#' * \code{\link{testModel.qda}}
#' @export
testModel <- function(model, ...) {
    UseMethod("testModel")
}

#' Computes a classification from a pretrained discriminant
#'
#' This function uses a pretrained linear discriminant to classify a set of
#' test data. As output it returns a confusion matrix and optionally the raw
#' classification result.
#'
#' @param model Trained linear discriminant (lda object).
#'              see \code{\link{trainModel}}
#' @param test WaveAnalysis class object to be used as test set.
#' @param grps Vector that determines the class to which each of the
#'             observations provided in the test set belongs.
#' @param returnClassification Allows to select if the raw result classification
#'        is returned.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'  object of the confusionMatrix class and a vector with the
#'  classification result.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' # The dataset has the first 5 elements of class 1
#' # and the last 5 of class 2.
#' grps <- c(rep(1, 5), rep(2, 5))
#' MWA <- generateStepDiscrim(ECGExample, grps, "haar", maxvars = 5, features = c("var"))
#' aux <- extractSubset(MWA, c(1, 2, 9, 10))
#' MWATest <- aux[[1]]
#' MWATrain <- aux[[2]]
#' ldaDiscriminant <- trainModel(MWATrain, grps[3:8], "linear")
#' CM <- testModel(ldaDiscriminant, MWATest, grps[c(1, 2, 9, 10)])
#' }
#'
#' @seealso
#' \code{\link{testModel}}
#'
#' @importFrom checkmate anyMissing assertFlag
#' @export
testModel.lda <- function(model,
                          test,
                          grps,
                          returnClassification = FALSE,
                          ...) {
    checkmate::anyMissing(c(model, test, grps))
    checkmate::assertFlag(returnClassification)
    stopifnot(is(model, "lda") || is(model, "qda"))
    if (test$Observations != length(grps)) {
        stop(
            "The number of observations in the test set does not correspond to the
         classes provided in the grps parameter."
        )
    }

    if (length(grps) < 2) {
        stop(
            "The minimun numer of observations is 2. If you want to classify only
         one observation use \"classify\" function"
        )
    }


    prediction <- classify(model, test)
    CM <- confusionMatrix(as.factor(prediction), as.factor(grps))

    if (returnClassification) {
        return(list("CM" = CM, "Clasification" = prediction[[1]]))
    }

    return(CM)
}

#' testModel.qda
#'
#' This function uses a pretrained quadratic discriminant to classify a set of
#' test data. As output it returns a confusion matrix and optionally the raw
#' classification result.
#'
#' @param model Trained quadratic discriminant (instance of class qda).
#'        see \code{\link{trainModel}}
#' @param test WaveAnalysis class object to be used as test set.
#' @param grps Vector that determines the class to which each of the
#'             observations provided in the test set belongs.
#' @param returnClassification Allows to select if the raw result classification
#'        is returned.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'    object of the confusionMatrix class and a vector with the
#'    classification result.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' # The dataset has the first 5 elements of class 1
#' # and the last 5 of class 2.
#' grps <- c(rep(1, 5), rep(2, 5))
#' MWA <- generateStepDiscrim(ECGExample, grps, "haar", maxvars = 2, features = c("var"))
#' aux <- extractSubset(MWA, c(1, 2, 9, 10))
#' MWATest <- aux[[1]]
#' MWATrain <- aux[[2]]
#' qdaDiscriminant <- trainModel(MWATrain, grps[3:8], "quadratic")
#' CM <- testModel(qdaDiscriminant, MWATest, grps[c(1, 2, 9, 10)])
#' }
#'
#' @export
testModel.qda <- function(model,
                          test,
                          grps,
                          returnClassification = FALSE,
                          ...) {
    testModel.lda(model, test, grps, returnClassification)
}

#' Leave-One-Out Cross Validation
#'
#' This function performs the Leave-One-Out Cross Validation (LOOCV) process
#' with different types of input parameters.
#'
#' @param data Starting data to generate the validation. It can be either the
#'        raw data, or a previously generated WaveAnalysis object.
#' @param ... Additional arguments
#'
#' @seealso
#' * \code{\link{LOOCV.array}}
#' * \code{\link{LOOCV.WaveAnalysis}}
#'
#' @importFrom checkmate anyMissing assertFlag
#' @export
LOOCV <- function(data, ...) {
    UseMethod("LOOCV")
}

#' Generates and validates a discriminant model generated directly from the
#' data.
#'
#' It generates and validates a discriminant model starting from the data.
#' First, a WaveAnalysis object is obtained according to the selected
#' characteristics, filter and levels. Then, the most important features are
#' selected using a stepwise discriminant that allows to select a maximum number
#' of variables (maxvars) or a minimum enhancement step (VStep). Finally, the
#' model is trained using the subset of features and validated using
#' Leave-One-Out Cross Validation (LOOCV).
#'
#' @param data Sample from the population (dim x length x cases)
#' @param grps Labeled vector that classify the observations
#' @param f Selected filter for the MODWT (to see the available filters use the
#'        function \code{\link{availableFilters}}
#' @param method Selected method for the discriminant.
#'        Valid values "linear" "quadratic"
#' @param maxvars Maximum number of variables included by the StepDiscrim
#'        algorithm (Note that if you defined this, can not define VStep).
#'        Must be a positive integer greater than 0.
#' @param VStep Minimum value of V above which all other variables are
#'        considered irrelevant and therefore will not be included. (Note that
#'        if you defined this, can not defined maxvars). Must be a positive
#'        number and greater than 0. For more information see StepDiscrim
#'        documentation
#' @param lev Determines the number of decomposition levels for MODWT
#'        (by default the optimum is calculated using the "conservative"
#'        strategy). Must be a positive integer (including 0 to auto-select
#'        the level)
#' @param features A list of characteristics that will be used for the
#'        classification process. To see the available features see
#'        \code{\link{availableFeatures}}
#' @param returnClassification Allows to select if the raw result classification
#'        is returned.
#' @param nCores Determines the number of processes that will be used in the
#'        function, by default it uses all but one of the system cores. Must be
#'        a positive integer, where 0 corresponds to the default behavior.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'    object of the confusionMatrix class and a vector with the
#'    classification result.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' grps <- c(rep(1, 5), rep(2, 5))
#' CM <- LOOCV(ECGExample, grps, "haar", "linear",
#'   maxvars = 5,
#'   features = c("Var"), returnClassification = FALSE
#' )
#' # or with VStep
#' CMV <- LOOCV(ECGExample, grps, "haar", "linear",
#'  VStep = 5,
#'  features = c("Var", "Cor"), returnClassification = FALSE
#' )
#' }
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{LOOCV.WaveAnalysis}}
#' * \code{\link{availableFilters}}
#' * \code{\link{availableFeatures}}
#'
#' @importFrom checkmate anyMissing
#' @export
LOOCV.array <-
    function(data,
             grps,
             f,
             method,
             maxvars,
             VStep,
             lev = 0,
             features = c("Var", "Cor", "IQR", "PE", "DM"),
             returnClassification = FALSE,
             nCores = 0,
             ...) {
        checkmate::anyMissing(c(data, grps, f, method))


        if (length(features) == 0) {
            stop(
                "At least one feature must be provided. To see the available filters
         use the availableFeatures()"
            )
        }

        if (length(dim(data)) != 3) {
            stop(
                "It seems that a dimension is missing, in case your series contains
         only one case, make sure that you have activated the option
         \"drop = FALSE\" as in the following example
         Series1 = Series2 [,,1, drop = FALSE]."
            )
        }

        checkmate::assertFlag(returnClassification)

        if (missing(maxvars) && missing(VStep)) {
            stop("maxvars o VStep must be defined")
        }

        if (!missing(maxvars) && !missing(VStep)) {
            stop("only maxvars or VStep can be defined.")
        }

        if (!missing(maxvars)) {
            maxvars <- asCount(maxvars, positive = TRUE)
        } else {
            if (!is.numeric(VStep) || length(VStep) != 1 || VStep <= 0) {
                stop("The argument \"VStep\" must be a number greater than 0")
            }
        }

        XSeries <- data

        f <- tolower(f)
        method <- tolower(method)
        features <- tolower(features)

        MWA <-
            generateStepDiscrim(XSeries, grps, f, maxvars, VStep, lev, features,
                                nCores)
        return(LOOCV(MWA, grps, method, returnClassification))
    }

#' LOOCV
#'
#' Performs a leave-one-cross-validation (LOOCV) method on a WaveAnalysis
#' object. It is advisable to have selected a subset of all features
#' (\code{\link{StepDiscrim}},\code{\link{StepDiscrimV}})
#'
#' @param data WaveAnalysis object obtained with MultiWaveAnalysis and
#'        preferably obtained a subset of its characteristics
#'        (\code{\link{StepDiscrim}}, \code{\link{StepDiscrimV}})
#' @param grps Labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options
#'        "linear" "quadratic"
#' @param returnClassification Allows to select if the raw result classification
#'       is returned.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'    object of the confusionMatrix class and a vector with the
#'    classification result.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' MWA <- MultiWaveAnalysis(ECGExample, "haar", features = c("var"))
#' MWADiscrim <- StepDiscrim(MWA, c(rep(1, 5), rep(2, 5)), 5)
#' CM <- LOOCV(MWADiscrim, c(rep(1, 5), rep(2, 5)), "linear")
#' }
#'
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{LOOCV.array}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#'
#' @export
#'
#' @importFrom caret confusionMatrix
LOOCV.WaveAnalysis <-
    function(data,
             grps,
             method,
             returnClassification = FALSE,
             ...) {
        if (missing(data)) {
            stop("The argument \"data\" must be provided")
        }
        if (missing(grps)) {
            stop("The argument \"grps\" must be provided")
        }
        if (missing(method)) {
            stop("The argument \"method\" must be provided")
        }
        if (is.numeric(returnClassification)) {
            stop("The argument \"returnClassification\" must be a logical value")
        }

        if (length(grps) != data$Observations) {
            stop("The \"grps\" length mismatches with the observations of \"MWA\"")
        }

        stopifnot(is(data, "WaveAnalysis"))

        MWA <- data

        method <- tolower(method)

        n <- MWA$Observations
        class <- vector("numeric", n)
        for (i in seq_len(n)) {
            aux <- extractSubset(MWA, c(i))
            MWATest <- aux[[1]]
            MWATrain <- aux[[2]]
            grpsT <- grps[-i]
            model <- trainModel(MWATrain, grpsT, method)
            class[i] <- classify(model, MWATest)
        }
        CM <- confusionMatrix(as.factor(class), as.factor(grps))

        if (returnClassification) {
            return(list("CM" = CM, "classification" = class))
        }

        return(CM)
    }

#' K-Fold Cross Validation (KFCV)
#'
#' This function performs the K-Fold Cross Validation (KFCV) process
#' with different types of input parameters.
#'
#' @param data Starting data to generate the validation. It can be either the
#'        raw data, or a previously generated WaveAnalysis object.
#' @param ... Additional arguments
#'
#' @seealso
#' * \code{\link{KFCV.array}}
#' * \code{\link{KFCV.WaveAnalysis}}
#' @export
KFCV <- function(data, ...) {
    UseMethod("KFCV")
}


#' Generates and validates a discriminant model generated directly from the
#' data.
#'
#' It generates and validates a discriminant model starting from the data. First
#' , a WaveAnalysis object is obtained according to the selected characteristics
#' ,filter and levels. Then, the most important features are selected using a
#' stepwise discriminant that allows to select a maximum number of variables
#' (maxvars) or a minimum enhancement step (VStep). Finally, the model is
#' trained using the subset of features and validated using
#' K-Fold Cross Validation (KFCV).
#'
#' @param data Sample from the population (dim x length x cases)
#' @param grps Labeled vector that classify the observations
#' @param f Selected filter for the MODWT (to see the available filters use the
#'        function \code{\link{availableFilters}}
#' @param method Selected method for the discriminant.
#'        Valid values "linear" "quadratic"
#' @param maxvars Maximum number of variables included by the StepDiscrim
#'        algorithm (Note that if you defined this, can not define VStep).
#'        Must be a positive integer greater than 0.
#' @param VStep Minimum value of V above which all other variables are
#'        considered irrelevant and therefore will not be included. (Note that
#'        if you defined this, can not defined maxvars). Must be a positive
#'        number and greater than 0. For more information see StepDiscrim
#'        documentation
#' @param k The number of folds in KFCV. Must be a positive integer lower or
#'        equal than the number of observations
#' @param lev Determines the number of decomposition levels for MODWT
#'        (by default the optimum is calculated using the "conservative"
#'        strategy). Must be a positive integer (including 0 to auto-select
#'        the level)
#' @param features A list of characteristics that will be used for the
#'        classification process. To see the available features see
#'        \code{\link{availableFeatures}}
#' @param returnClassification Allows to select if the raw result classification
#'        is returned.
#' @param nCores Determines the number of processes that will be used in the
#'        function, by default it uses all but one of the system cores. Must be
#'        a positive integer, where 0 corresponds to the default behavior.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'    object of the confusionMatrix class and a vector with the
#'    classification result.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' grps <- c(rep(1, 5), rep(2, 5))
#' CM <- KFCV(ECGExample, grps, "haar", "linear",
#'   maxvars = 5,
#'   features = c("Var"), returnClassification = FALSE
#' )
#' # or with VStep
#' CMV <- KFCV(ECGExample, grps, "haar", "linear",
#'  k = 5,
#'  VStep = 5,
#'  features = c("Var", "Cor"), returnClassification = FALSE
#' )
#' }
#' @seealso
#' * \code{\link{LOOCV}}
#' * \code{\link{LOOCV.WaveAnalysis}}
#' * \code{\link{availableFilters}}
#' * \code{\link{availableFeatures}}
#' @export
KFCV.array <-
    function(data,
             grps,
             f,
             method,
             maxvars,
             VStep,
             k = 5L,
             lev = 0L,
             features = c("Var", "Cor", "IQR", "PE", "DM"),
             returnClassification = FALSE,
             nCores = 0,
             ...) {
        checkmate::anyMissing(c(data, grps, f, method, features))

        if (length(features) == 0) {
            stop(
                "At least one feature must be provided. To see the available filters
         use the availableFeatures()"
            )
        }

        if (length(dim(data)) != 3) {
            stop(
                "It seems that a dimension is missing, in case your series contains
         only one case, make sure that you have activated the option
         \"drop = FALSE\" as in the following example
         Series1 = Series2 [,,1, drop = FALSE]."
            )
        }


        if (is.numeric(returnClassification)) {
            stop("The argument \"returnClassification\" must be a logical value")
        }

        if (missing(maxvars) && missing(VStep)) {
            stop("maxvars o VStep must be defined")
        }

        if (!missing(maxvars) && !missing(VStep)) {
            stop("only maxvars or VStep can be defined.")
        }

        if (!missing(maxvars)) {
            if (!is.numeric(maxvars) || length(maxvars) != 1 || maxvars <= 0) {
                stop("The argument \"maxvars\" must be an integer greater than 0")
            }
        } else {
            if (!is.numeric(VStep) || length(VStep) != 1 || VStep <= 0) {
                stop("The argument \"VStep\" must be a number greater than 0")
            }
        }

        k <- checkmate::asCount(k, positive = TRUE)
        lev <- checkmate::asCount(lev)
        nCores <- checkmate::asCount(nCores)


        XSeries <- data

        f <- tolower(f)
        method <- tolower(method)
        features <- tolower(features)

        MWA <-
            generateStepDiscrim(XSeries, grps, f, maxvars, VStep, lev, features,
                                nCores)
        return(KFCV(MWA, grps, method, k, returnClassification))
    }

#' KFCV
#'
#' Performs k-fold cross-validation where groups are chosen randomly.
#' In case the value k is not divisor of the number of observations the last
#' group will have nobs mod k observations.
#'
#' @param data WaveAnalysis (MWA) object obtained with MultiWaveAnalysis and
#'        preferably obtained a subset of its characteristics
#'        (\code{\link{StepDiscrim}},\code{\link{StepDiscrimV}})
#' @param grps labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options
#'       "linear" "quadratic"
#' @param k the number of folds in KFCV. Must be a positive integer and lower or
#'        equal than the number of observations
#' @param returnClassification Allows to select if the raw result classification
#'        is returned.
#' @param ... Additional arguments
#'
#' @return * if returnClassification is false return a object of class
#'            confusionMatrix
#'  * if returnClassification is true, it returns a list containing an
#'    object of the confusionMatrix class and a vector with the
#'    classification result.
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' MWA <- MultiWaveAnalysis(ECGExample, "haar", features = c("var"))
#' MWADiscrim <- StepDiscrim(MWA, c(rep(1, 5), rep(2, 5)), 5)
#' CM <- KFCV(MWADiscrim, c(rep(1, 5), rep(2, 5)), "linear", 5,
#'   returnClassification = FALSE
#' )
#' }
#'
#' @importFrom checkmate anyMissing assertFlag
#'
#' @export
KFCV.WaveAnalysis <- function(data,
                              grps,
                              method,
                              k = 5L,
                              returnClassification = FALSE,
                              ...) {
    checkmate::anyMissing(c(data, grps, method))
    checkmate::assertFlag(returnClassification)

    if (length(grps) != data$Observations) {
        stop("The \"grps\" length mismatches with the observations of \"MWA\"")
    }

    k <- asCount(k)

    MWA <- data

    nobs <- MWA$Observations

    if (k > nobs) {
        k <- nobs
        warning(
            "The number of groups is greater than the number of
                  observations. The value of k:",
            k,
            "will be used"
        )
    }

    method <- tolower(method)

    reorderMWA <- MWA
    cols <- sample(nobs)

    for (feature in names(reorderMWA$Features)) {
        if (!all(is.na(reorderMWA$Features[feature]))) {
            reorderMWA$Features[[feature]] <-
                reorderMWA$Features[[feature]][, cols, drop = FALSE]
        }
    }

    reorderGrps <- grps[cols]
    class <- vector("numeric", nobs)
    inc <- floor(nobs / k)
    n <- 1
    m <- inc
    while (m >= n) {
        aux <- extractSubset(reorderMWA, n:m)
        MWATest <- aux[[1]]
        MWATrain <- aux[[2]]
        grpsT <- reorderGrps[-c(n:m)]
        model <- trainModel(MWATrain, grpsT, method)
        class[n:m] <- classify(model, MWATest)
        n <- n + inc
        m <- m + inc
        m <- min(m, nobs)
    }

    CM <- confusionMatrix(as.factor(class), as.factor(reorderGrps))

    if (returnClassification) {
        return(list("CM" = CM, "classification" = class))
    }

    return(CM)
}

#' Generate a Discriminant Model
#'
#' This function allows training of a discriminant model using different inputs
#'
#' @param data Starting data to generate a discriminator (linear or quadratic).
#'        This starting data can be either the raw data, or a WaveAnalysis
#'        object generated earlier.
#' @param ... Additional arguments
#'
#' @return A trained discriminant model
#'
#' @seealso
#' * \code{\link{trainModel.array}}
#' * \code{\link{trainModel.WaveAnalysis}}
#'
#' @export
#'
trainModel <- function(data, ...) {
    UseMethod("trainModel")
}

#' Generates a discriminant model from training data.
#'
#' It generates a discriminant model starting from the training data, which must
#' be provided in 2 groups depending on their classification. The method first
#' obtains the variances and correlations using MODWT, the f filter is applied
#' with a number of levels lev. Then a subset of all the generated features
#' will be obtained by means of a stepwise discriminant, which can be driven
#' by a maximum number of features or by a minimum metric to be met. Finally,
#' the selected discriminant model is trained with the subset obtained.
#'
#' @param data Sample from the population (dim x length x cases)
#' @param grps Labeled vector that classify the observations
#' @param f Selected filter for the MODWT (to see the available filters use the
#'        function availableFilters)
#' @param method Selected method for the discriminant. Valid values
#'        "linear" "quadratic"
#' @param maxvars Maximum number of variables included by the StepDiscrim
#'        algorithm (Note that if you defined this, can not define VStep). Must
#'        be a positive integer greater than 0.
#' @param VStep Minimum value of V above which all other variables are
#'        considered irrelevant and therefore will not be included. (Note that
#'        if you defined this, can not defined maxvars).Must be a positive
#'        number greater than 0. For more information see StepDiscrim
#'        documentation
#' @param lev Determines the number of decomposition levels for MODWT
#'        (by default the optimum is calculated). Must be a positive integer
#' @param features A list of characteristics that will be used for the
#'        classification process. To see the available features
#'        see \code{\link{availableFeatures}}
#' @param nCores Determines the number of processes that will be used in the
#'        function, by default it uses all but one of the system cores. Must be
#'        a positive integer, where 0 corresponds to the default behavior.
#' @param ... Additional arguments
#'
#' @return A discriminant model object (lda or qda)
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' # The dataset has the first 5 elements of class 1 and the last 5 of class 2.
#' grps <- c(rep(1, 5), rep(2, 5))
#' model <- trainModel(ECGExample, grps, "d6", "linear",
#'   maxvars = 5, features = c("Var")
#' )
#' # or using VStep
#' modelV <- trainModel(ECGExample, grps, "d6", "linear",
#'     VStep = 14.5, features = c("Var")
#' )
#' }
#' @export
#' @seealso
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' * \code{\link{trainModel}}
#' @md
trainModel.array <-
    function(data,
             grps,
             f,
             method,
             maxvars,
             VStep,
             lev = 0,
             features = c("Var", "Cor", "IQR", "PE", "DM"),
             nCores = 0,
             ...) {
        if (missing(data)) {
            stop("The argument \"XSeries\" must be provided")
        }
        if (missing(grps)) {
            stop("The argument \"grps\" must be provided")
        }
        if (missing(f)) {
            stop(
                "The argument \"f\" (filter) must be provided. To see the avaiable
         filters use availableFilters()"
            )
        }

        if (missing(method)) {
            stop(
                "The argument \"method\" must be defined. The available  methods are
         \"linear\" and \"quadratic\""
            )
        }

        if (length(features) == 0) {
            stop(
                "At least one feature must be provided. To see the available filters
         use the availableFeatures()"
            )
        }

        if (length(dim(data)) != 3) {
            stop(
                "It seems that a dimension is missing, in case your series contains
         only one case, make sure that you have activated the option
         \"drop = FALSE\" as in the following example
         Series1 = Series2 [,,1, drop = FALSE]."
            )
        }

        if (length(grps)) {
            if (missing(maxvars) && missing(VStep)) {
                stop("maxvars o VStep must be defined")
            }
        }
        if (!missing(maxvars) && !missing(VStep)) {
            stop("only maxvars or VStep can be defined.")
        }

        if (!missing(maxvars)) {
            if (!is.numeric(maxvars) || length(maxvars) != 1 || maxvars <= 0) {
                stop("The argument \"maxvars\" must be an integer greater than 0")
            }
        } else {
            if (!is.numeric(VStep) || length(VStep) != 1 || VStep <= 0) {
                stop("The argument \"VStep\" must be a number greater than 0")
            }
        }

        XSeries <- data

        f <- tolower(f)
        method <- tolower(method)
        features <- tolower(features)

        MWA <-
            generateStepDiscrim(XSeries, grps, f, maxvars, VStep, lev,
                                features, nCores)
        return(trainModel(MWA, grps, method))
    }
#' Generates a discriminant model from an already generated "WaveAnalysis".
#'
#' @param data A WaveAnalysis object obtained with MultiWaveAnalysis function
#' @param grps Labeled vector that classify the observations.
#' @param method Selected method for discrimination. Valid options are
#'        "linear" and "quadratic"
#' @param ... Additional arguments
#'
#' @return A discriminant model based on selected method. It can be an object of
#'         the class lda or qda.
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' MWA <- MultiWaveAnalysis(ECGExample, "d6", features = c("Var"))
#' MWADiscrim <- StepDiscrim(MWA, c(rep(1, 5), rep(2, 5)), 5,
#'   features = c("Var")
#' )
#' model <- trainModel(MWADiscrim, c(rep(1, 5), rep(2, 5)), "linear")
#' }
#'
#' @seealso
#' * \code{\link{MultiWaveAnalysis}}
#' * \code{\link{StepDiscrim}}
#' * \code{\link{StepDiscrimV}}
#' @export
#'
#' @importFrom MASS lda qda
trainModel.WaveAnalysis <- function(data, grps, method, ...) {
    MWA <- data
    if (missing(MWA)) {
        stop("The argument \"MWA\" must be provided")
    }
    if (missing(grps)) {
        stop("The argument \"grps\" must be provided")
    }

    if (missing(method)) {
        stop(
            "The argument \"method\" must be defined. The available  methods are
         \"linear\" and \"quadratic\""
        )
    }

    if (length(grps) != data$Observations) {
        stop("The \"grps\" length mismatches with the observations of \"MWA\"")
    }

    if (length(method) > 1) {
        stop("The argument must be a single string
         (\"linear\" or \"quadratic\") not a list")
    }

    stopifnot(is(MWA, "WaveAnalysis"))

    method <- tolower(method)

    values <- values(MWA)
    if (method == "linear") {
        model <- lda(t(values), grps)
    } else if (method == "quadratic") {
        model <- qda(t(values), grps)
    } else {
        stop("Method", as.character(method), "not supported")
    }
    return(model)
}

#' Classifies observations based on a pretrained model.
#'
#' This function allows to classify observations based on a pretrained model
#' that could have been obtained in several ways (such as using the train model
#' function).
#'
#' @param model pretrained discriminant model (lda or qda)
#' @param data Data to be classified by the model. Remember that it must be an
#'        object of type WaveAnalysis. Note that it should have the same
#'        variables selected as those used to generate the model.
#'
#' @return A factor with predicted class of each observation
#'
#' @examples
#' \donttest{
#' ECGExample <- loadECGExample()
#' # We simulate that the second series has been obtained after
#' Series1 <- ECGExample[, , 1:9]
#' Series2 <- ECGExample[, , 10, drop = FALSE]
#'
#' # Training a discriminant model
#' MWA <- MultiWaveAnalysis(Series1, "haar", features = c("var"))
#' MWADiscrim <- StepDiscrim(MWA, c(rep(1, 5), rep(2, 4)), maxvars = 5)
#' model <- trainModel(MWADiscrim, c(rep(1, 5), rep(2, 4)), "linear")
#'
#' # Using the discriminant trained on new data
#' MWA2 <- MultiWaveAnalysis(Series2, "haar", features = c("var"))
#' MWA2Discrim <- SameDiscrim(MWA2, MWADiscrim)
#' prediction <- classify(model, MWA2Discrim)
#' }
#'
#' @seealso
#' * \code{\link{trainModel}}
#'
#' @export
#'
#' @importFrom stats predict
#' @importFrom magrittr %>%
classify <- function(model, data) {
    if (missing(model)) {
        stop("The argument \"model\" must be provided")
    }
    stopifnot(is(data, "WaveAnalysis"))
    stopifnot(is(model, "lda") || is(model, "qda"))

    values <- values(data)
    return(predict(model, t(values))[[1]])
}
