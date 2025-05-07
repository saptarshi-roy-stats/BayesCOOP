# BayesCOOP

The repository contains the analysis codes from the Bayesian coopertive
project.

## Dependencies

`BayesCOOP` requires the following `R` package: `devtools` (for
installation only). Please install it before installing `BayesCOOP`,
which can be done as follows (execute from within a fresh R session):

    install.packages("devtools")
    library(devtools)

## Installation

Once the dependencies are installed, `BayesCOOP` can be loaded using the
following command:

    devtools::install_github("himelmallick/BayesCOOP")
    library(BayesCOOP)

## Load libraries

    library(BhGLM)
    library(tidyverse)
    library(caret)

## sourcing required functions

    rm(list = ls())
    source("~/bayesCOOP/BC/bayesCoop.R")

    ## Loading required package: coda

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## ##
    ## ## Markov Chain Monte Carlo Package (MCMCpack)

    ## ## Copyright (C) 2003-2025 Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park

    ## ##
    ## ## Support provided by the U.S. National Science Foundation

    ## ## (Grants SES-0350646 and SES-0350613)
    ## ##

    ## 
    ## Attaching package: 'rmutil'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     nesting

    ## The following object is masked from 'package:BhGLM':
    ## 
    ##     covariates

    ## The following object is masked from 'package:stats':
    ## 
    ##     nobs

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.data.frame, units

## Loading the StelzerDOS real dataset

    data_train = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS.RData"))); 
    rm(pcl)

    data_test = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS_valid.RData"))); 
    rm(pcl)

## Pre-processing the longitudinal data by considering only baseline observations

### Remove metabolomics from the train set to match with validation

    data_train$feature_metadata = data_train$feature_metadata %>% dplyr::filter(featureType!='Metabolomics')
    data_train$feature_table = data_train$feature_table[rownames(data_train$feature_metadata),]

### Consider only baseline observations for the train set

    positions = grep("A", colnames(data_train$feature_table), ignore.case = TRUE)
    data_train$feature_table = data_train$feature_table[, positions]
    data_train$sample_metadata = data_train$sample_metadata[positions, ]
    rm(positions)

### Consider only baseline observations for the validation set

    positions = grep("G1", colnames(data_test$feature_table))
    data_test$feature_table = data_test$feature_table[, positions]
    data_test$sample_metadata = data_test$sample_metadata[positions, ]
    rm(positions)

# BayesCoop Implementation

    set.seed(1)
    result = bayesCoop(data_train, data_test, family = "gaussian", 
                       ss = c(0.05, 1), group = TRUE,
                       alpha_dirich = 1, maxit = 1, 
                       bbiters = 11, bbburn = 10,
                       abd_thresh = 0, prev_thresh = 0.1,
                       Warning = TRUE, verbose = TRUE)

    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.023 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.022 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.035 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.039 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.027 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.024 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.015 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.024 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.013 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.023 minutes 
    ## EM Coordinate Decent Iterations: 1 
    ## Computational time: 0.019 minutes

    #result$beta_postmed

    result$rho_postmed

    ## [1] 0.0005647181

    #result$beta_samples

    result$y_pred

    ## [1]  25.5238395  -6.8344675 -10.3726248   9.7306132   0.8182733 -16.7161457
    ## [7] -10.9832246   8.8337367

    result$mspe

    ## [1] 483.148

    result$time

    ## [1] 0.331
