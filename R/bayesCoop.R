library(truncnorm)
library(MCMCpack)
library(rmutil)

#### source helper functions #### (please change the path accordingly)
source("https://github.com/himelmallick/BayesCOOP/blob/master/R/helperFunctions.R")

# Function to implement bayesCoop
# data_train: a list of feature_table, sample_metadata and feature_metadata from training data. {Cite integrated learner for more references}.
# data_test: a list of feature_table, sample_metadata and feature_metadata from testing data. {Cite integrated learner for more references}.
# family: currently supports only Gaussian family. Default value: "gaussian".
# ss = (s0, s1): s0 and s1 respectively stands for the scale parameters of the spike and slab parts, s0 < s1. Default value: s0 = 0.05, s1 = 1.
# group: logical. If TRUE, group is a list of variable names, group[[k]] includes variables in the k-th group. For ungrouped predictors, the prior is double-exponential with scale ss[2] and mean 0. If FALSE, group = NULL. Default value: TRUE.
# alpha_dirich: shape parameter for Dirichlet distribution. Default value: 1.
# maxit: integer giving the maximal number of EM iterations.
# bbiters: number of Bayesian bootstrap iterations. Default value: 1100.
# bbburn: number of burn-in samples. Default value: 100.
# abd_thresh: threshold for abundance of features. Default value: 0.
# prev_thresh: threshold for prevalence of features. Default value: 0.1.

bayesCoop <- function(data_train, data_test, family = "gaussian", 
                     ss = c(0.05, 1), group = TRUE,
                     alpha_dirich = 1, maxit = 100, 
                     bbiters = 1100, bbburn = 100,
                     abd_thresh = 0, prev_thresh = 0.1,
                     Warning = TRUE, verbose = TRUE) {
  
  #########################################################################
  ######################### data pre-processing ###########################
  #########################################################################
  
  y_train <- gen_datalist(data_train)$y
  y_test <- gen_datalist(data_test)$y
  xList_train <- gen_datalist(data_train)$xList
  xList_test <- gen_datalist(data_test)$xList
  
  ## Filtering features
  xList_train <- lapply(xList_train, function(foo) filter_features(foo, abd_thresh, prev_thresh))
  xList_test <- lapply(xList_test, function(foo) filter_features(foo, abd_thresh, prev_thresh))
  
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
  
  #########################################################################
  
  #########################################################################
  ##################### BayesCoop Implementation ##########################
  #########################################################################
  
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
  
  # Check if group information is provided or not
  if(group == TRUE) 
    group <- lapply(xList_train_std, names) else group <- NULL
  
  cts <- 0
  beta_samples <- matrix(NA, bbiters - bbburn, length(beta0))
  rho_samples <- numeric(length = bbiters - bbburn)
  errVar_samples <- numeric(length = bbiters - bbburn)
  
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
    if(length(xList_train_std) == 2) {
      dataAug_train <- data.Augment.2views(y_train, xList_train_std, rho)
    } else {
      dataAug_train <- data.Augment.mviews(y_train, xList_train_std, rho)
    }
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
  
  # posterior median of beta and rho
  beta_postmed <- apply(beta_samples, 2, function(foo) quantile(foo, 0.5))
  rho_postmed <- median(rho_samples)
  time <- as.numeric(round(difftime(stop.time, start.time, units="min"), 3), units = "mins")
  
  #########################################################################
  #################### Prediction performance #############################
  #########################################################################
  
  ## Augmented test data
  if(length(xList_train) == 2) {
    dataAug_test <- data.Augment.2views(y_test, xList_test_std, rho_postmed)
  } else {
    dataAug_test <- data.Augment.mviews(y_test, xList_test_std, rho_postmed)
  }
  y_aug_test <- dataAug_test$y_aug; x_aug_test <- as.matrix(dataAug_test$x_aug)
  y_aug_hat <- x_aug_test %*% beta_postmed
  y_pred <- y_aug_hat[(1:length(y_test))]
  mspe <- mean((y_pred - y_test)^2)
  
  #########################################################################
  ################################ Return #################################
  #########################################################################
  return(list(beta_postmed = beta_postmed,
              rho_postmed = rho_postmed,
              beta_samples = beta_samples,
              y_pred = y_pred,
              mspe = mspe,
              time = time))
}
