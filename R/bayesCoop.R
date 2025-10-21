#' @title bayesCoop
#'
#' @description This function implements the BayesCOOP methodology for supervised multimodal integration. It combines jittered group spike-and-slab LASSO regularization with intermediate fusion to enable integrative learning across multiple data modalities. For uncertainty quantification, BayesCOOP applies the Bayesian bootstrap to generate approximate posterior samples by performing maximum a posteriori (MAP) estimation on jittered and resampled multimodal datasets. Currently, only continuous outcomes are supported.
#'
#' @param data_train a list of feature_table, sample_metadata and feature_metadata from training data(unstandardized). See \emph{IntegratedLearner} for more details.
#' @param data_test a list of feature_table, sample_metadata and feature_metadata from testing data(unstandardized). See \emph{IntegratedLearner} for more details..
#' @param family currently supports only Gaussian family. Default value: "gaussian".
#' @param ss a two-dimensional vector representing the scale parameters of the spike and slab parts, namely, s0 and s1 respectively, with s0 < s1. Default value: ss = c(s0, s1) = c(0.05, 1).
#' @param group Logical. If \code{TRUE}, \code{group} is a list of variable names where \code{group[[k]]} includes variables in the \eqn{k}-th group. For ungrouped predictors, the prior is double-exponential with scale \code{ss[2]} and mean 0. If \code{FALSE}, \code{group = NULL}. Default is \code{TRUE}.
#' @param bb logical. If TRUE, Bayesian bootstrap is implemented, else, maximum a posteriori (MAP) estimation is implemented which requires additional parameters to be passed to \code{control}. Default value: TRUE.
#' @param alpha_dirich shape parameter for Dirichlet distribution. Default value: 1.
#' @param bbiters number of Bayesian bootstrap iterations. Default value: 1100.
#' @param bbburn number of burn-in samples. Default value: 100.
#' @param maxit integer giving the maximal number of EM iterations. Default value: 100.
#' @param filter logical. If TRUE, features are filtered as per abundance and prevalence thresholds. Default value: TRUE.
#' @param abd_thresh threshold for abundance of features. Default value: 0.
#' @param prev_thresh threshold for prevalence of features. Default value: 0.1.
#' @param Warning logical. If TRUE, show the error messages of not convergence. Default value: TRUE.
#' @param verbose logical. If TRUE, print out number of iterations and computational time. Default value: TRUE.
#' @param control a named list with one element `rho`, which is a numeric vector of user-specified rho values.
#' @return a list containing the following components is returned if `bb` is TRUE.
#' \item{y_samples}{approximate posterior predictive samples}
#' \item{y_pred}{predicted values of response}
#' \item{mspe}{mean square prediction error}
#' \item{beta_samples}{approximate posterior samples of beta}
#' \item{beta_postmed}{posterior median as an estimate for beta}
#' \item{rho_samples}{approximate posterior samples of rho}
#' \item{rho_postmed}{posterior median as an estimate for rho}
#' \item{errVar_samples}{approximate posterior samples of error variance}
#' \item{time}{time taken for BayesCoop to run}
#' a list containing the following components if `bb` is FALSE.
#' \item{y_pred}{predicted values of response}
#' \item{mspe}{mean square prediction error}
#' \item{beta_MAP}{MAP estimate of beta}
#' \item{rho_MAP}{optimal rho among the specified values chosen corresponding to minimum mean square prediction error, returns the same value if only one value for rho is specified in `control`}
#' \item{time}{time taken for MAP estimation to run}
#' @import BhGLM
#' @import caret
#' @import dplyr
#' @import glmnet
#' @import MCMCpack
#' @import rmutil
#' @import stats
#' @import truncnorm
#' @import survival
#' @export

bayesCoop <- function(data_train, data_test, family = "gaussian",
                      ss = c(0.05, 1), group = TRUE,
                      bb = TRUE, alpha_dirich = 1,
                      bbiters = 1100, bbburn = 100, maxit = 100,
                      filter = TRUE, abd_thresh = 0, prev_thresh = 0.1,
                      Warning = TRUE, verbose = TRUE, control = list()) {

    # Some basic checks
    if(length(data_train) != 3 | any(names(data_train) != c("feature_table", "sample_metadata", "feature_metadata")) == TRUE)
        stop("data_train must be a list of feature_table, sample_metadata and feature_metadata!")
    if(length(data_test) != 3 | any(names(data_test) != c("feature_table", "sample_metadata", "feature_metadata")) == TRUE)
        stop("data_test must be a list of feature_table, sample_metadata and feature_metadata!")
    if(length(ss) != 2 | ss[1] >= ss[2])
        stop("ss must be a two-dimensional vector c(s0, s1) with s0 < s1!")
    if(family != "gaussian")
        stop("Non-gaussian families are currently not supported!")
    if(bbiters <= bbburn)
        stop("(bbiters - bbburn) must be a positive integer!")

    output = list()

    #########################################################################
    ######################### data pre-processing ###########################
    #########################################################################

    y_train <- gen_datalist(data_train)$y
    y_test <- gen_datalist(data_test)$y
    xList_train <- gen_datalist(data_train)$xList
    xList_test <- gen_datalist(data_test)$xList

    ## Filtering features
    if(filter == TRUE){
        xList_train <- lapply(xList_train, function(foo) filter_features(foo, abd_thresh, prev_thresh))
        xList_test <- lapply(xList_test, function(foo) filter_features(foo, abd_thresh, prev_thresh))
    }

    ## Considering the common features between train and test set
    keep_features <- vector("list", length = length(xList_train))
    for(i in 1:length(keep_features)) {
        keep_features[[i]] <- intersect(names(xList_train[[i]]), names(xList_test[[i]]))
        xList_train[[i]] <- as.matrix(xList_train[[i]][, keep_features[[i]], drop = FALSE])
        xList_test[[i]] <- as.matrix(xList_test[[i]][, keep_features[[i]], drop = FALSE])
    }

    ## Scaling the predictors for bmlasso
    xList_train_std <- lapply(xList_train, function(foo) as.data.frame(scale(foo)))
    xList_test_std <- lapply(xList_test, function(foo) as.data.frame(scale(foo)))

    ## centering the response
    y_train <- y_train - mean(y_train)
    y_test <- y_test - mean(y_test)

    # Check if group information is provided or not
    if(group == TRUE)
        group <- lapply(xList_train_std, names) else group <- NULL

    #########################################################################
    ##################### BayesCoop Implementation ##########################
    #########################################################################

    if(bb) {
        message("Implementing BayesCoop...")

        ## this is the default implementation, i.e, implement jittered bb with sampling rho and error variance
        ## Initialization ##
        rho <- 0.5
        pis <- unlist(lapply(xList_train_std, ncol))
        beta0 <- runif(sum(pis), min = -5, max = 5)
        picums <- cumsum(pis)
        theta <- vector("list", length(xList_train_std))
        theta[[1]] <- beta0[1:picums[1]]
        for(ii in 2:length(xList_train_std)){
            theta[[ii]] <- beta0[(picums[ii-1]+1):picums[ii]]
        }
        errVar <- 1

        cts <- 0
        beta_samples <- matrix(NA, bbiters - bbburn, length(beta0))
        rho_samples <- numeric(length = bbiters - bbburn)
        errVar_samples <- numeric(length = bbiters - bbburn)
        y_samples <- matrix(NA, bbiters - bbburn, length(y_test))

        start.time <- Sys.time()
        for(its in 1:bbiters) {

            ###### Update rho ######
            denom <- 0
            for(i in 1:(length(xList_train_std) -1)) {
                for(j in (i + 1):length(xList_train_std)) {
                    denom <- denom + sum(drop(as.matrix(xList_train_std[[i]]) %*% theta[[i]] - as.matrix(xList_train_std[[j]]) %*% theta[[j]]) ^ 2) + 1e-5
                }
            }
            rho <- (truncnorm::rtruncnorm(1, a = 0, b = 1, mean = 0, sd = sqrt(errVar / denom))) ^ 2

            ## Augmented data for training set using updated rho
            dataAug_train <- data.Augment(y_train, xList_train_std, rho)
            y_aug_train <- dataAug_train$y_aug; x_aug_train <- as.matrix(dataAug_train$x_aug)


            ###### Update beta ######
            likelihood_weights <- nrow(x_aug_train) * MCMCpack::rdirichlet(1, rep(alpha_dirich, nrow(x_aug_train)))


            # draw mu_t
            jitter <- rmutil::rlaplace(n = (ncol(x_aug_train)+1), m = 0, s = ss[1]) ## uses the smaller scale, i.e., the scale parameter of spike density

            fit_bmlasso_bb <- bmlasso.weighted(x = x_aug_train, y = y_aug_train, family = family,
                                               maxit = maxit, alpha = 1, ss = ss,
                                               lhood_weights = likelihood_weights, jitter = jitter,
                                               group = group, Warning = Warning, verbose = verbose)

            beta_new <- fit_bmlasso_bb$coefficients[-1]

            theta[[1]] <- beta_new[1:picums[1]]
            for(kk in 2:length(xList_train_std)){
                theta[[kk]] <- beta_new[(picums[kk-1]+1):picums[kk]]
            }

            ###### Update error variance #######
            rss <- sum(drop(y_aug_train - x_aug_train %*% beta_new) ^ 2)
            errVar <- 1 / rgamma(1, shape = (nrow(x_aug_train) + 1) / 2, rate = (rss + 1) / 2)

            # store samples
            if(its > bbburn) {
                cts <- cts + 1
                beta_samples[cts, ] <- beta_new
                rho_samples[cts] <- rho
                errVar_samples[cts] <- errVar
            }

            message("iter = ", its)
        }
        stop.time <- Sys.time()

        # Estimated beta and rho
        beta_postmed <- apply(beta_samples, 2, function(foo) quantile(foo, 0.5))
        rho_postmed <- median(rho_samples)
        time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")

        #########################################################################
        #################### Prediction performance #############################
        #########################################################################

        ## Augmented test data
        dataAug_test <- data.Augment(y_test, xList_test_std, rho_postmed)
        y_aug_test <- dataAug_test$y_aug; x_aug_test <- as.matrix(dataAug_test$x_aug)

        y_aug_samples <- do.call("rbind", lapply(1:nrow(y_samples),
                                function(tt) {drop(x_aug_test %*% beta_samples[tt, ]) + rnorm(nrow(x_aug_test), 0, sqrt(errVar_samples[tt]))}))
        y_samples <- y_aug_samples[, 1:length(y_test)]
        y_pred <- apply(y_samples, 2, median)
        mspe <- mean((y_pred - y_test)^2)
        
        output$y_samples <- y_samples
        output$y_pred <- y_pred
        output$mspe <- mspe
        output$beta_samples <- beta_samples
        output$beta_postmed <- beta_postmed
        output$rho_samples <- rho_samples
        output$rho_postmed <- rho_postmed
        output$errVar_samples <- errVar_samples
        output$time <- time

    } else {

        if(is.null(control$rho)){
            stop("At least one non-null value of rho is needed for MAP estimation!")
        }

        rho_grid <- control$rho
        message("Implementing MAP estimation...")

        start.time <- Sys.time()
        fit <- bmlasso_cv(y = y_train, xList = xList_train_std,
                       rho_grid = rho_grid, group = group, ss = ss, family = family,
                       maxit = maxit, Warning = Warning, verbose = verbose)
        stop.time <- Sys.time()

        # Estimated beta and rho
        betahat_MAP <- fit$beta_hat_bmlasso
        rho_MAP <- fit$rho_bmlasso
        time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")

        #########################################################################
        #################### Prediction performance #############################
        #########################################################################

        ## Augmented test data
        dataAug_test <- data.Augment(y_test, xList_test_std, rho_MAP)
        y_aug_test <- dataAug_test$y_aug; x_aug_test <- as.matrix(dataAug_test$x_aug)
        y_aug_hat <- x_aug_test %*% betahat_MAP
        y_pred <- y_aug_hat[(1:length(y_test))]
        mspe <- mean((y_pred - y_test)^2)

        output$y_pred <- y_pred
        output$mspe <- mspe
        output$beta_MAP <- betahat_MAP
        output$rho_MAP <- rho_MAP
        output$time <- time

    }

    #########################################################################
    ################################ Return #################################
    #########################################################################
    output
}