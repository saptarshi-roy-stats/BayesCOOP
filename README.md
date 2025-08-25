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

    source("https://github.com/himelmallick/BayesCOOP/blob/master/R/helperFunctions.R")

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

### Using top features for the train set

    ########################
    # Prevalence filtering #
    ########################
     
    abd_threshold = 0
    prev_threshold = 0.1
    data_train$feature_table <- data_train$feature_table[rowMeans(data_train$feature_table > abd_threshold) > prev_threshold,]
    data_train$feature_metadata<-data_train$feature_metadata[rownames(data_train$feature_table),]
     
    ######################
    # Variance filtering #
    ######################
     
    data_train$feature_metadata$featureType<-as.factor(data_train$feature_metadata$featureType)
    name_layers<-with(droplevels(data_train$feature_metadata), list(levels = levels(data_train$feature_metadata$featureType)), nlevels = nlevels(data_train$featureType))$levels
    keep_features<-c()
     
    for (i in seq_along(name_layers)){
     
      ##################
      # Extract layers #
      ##################
     
      include_list<-data_train$feature_metadata %>% filter(featureType == name_layers[i])
      dat_slice<-data_train$feature_table[rownames(data_train$feature_table) %in% include_list$featureID, ]
     
      ################################
      # Per-layer variance filtering #
      ################################
     
      sds<-apply(dat_slice, 1, sd)
      dat_slice_filtered<-dat_slice[which(sds>quantile(sds, .95)),]
      keep_features<-c(keep_features, rownames(dat_slice_filtered))
      rm(dat_slice); rm(dat_slice_filtered)
    }
     
    #########################
    # Final filtered tables #
    #########################
     
    data_train$feature_table<-data_train$feature_table[keep_features, ]
    data_train$feature_metadata<-data_train$feature_metadata[keep_features,]
     
     
    ################
    # Sanity check #
    ################
     
    all(rownames(data_train$feature_table)==rownames(data_train$feature_metadata))

    ## [1] TRUE

    all(colnames(data_train$feature_table)==rownames(data_train$sample_metadata))

    ## [1] TRUE

### Using top features for the test set

    ########################
    # Prevalence filtering #
    ########################

    abd_threshold = 0
    prev_threshold = 0.1
    data_test$feature_table <- data_test$feature_table[rowMeans(data_test$feature_table > abd_threshold) > prev_threshold,]
    data_test$feature_metadata<-data_test$feature_metadata[rownames(data_test$feature_table),]

    ######################
    # Variance filtering #
    ######################

    data_test$feature_metadata$featureType<-as.factor(data_test$feature_metadata$featureType)
    name_layers<-with(droplevels(data_test$feature_metadata), list(levels = levels(data_test$feature_metadata$featureType)), nlevels = nlevels(data_test$featureType))$levels
    keep_features<-c()

    for (i in seq_along(name_layers)){
      
      ##################
      # Extract layers #
      ##################
      
      include_list<-data_test$feature_metadata %>% filter(featureType == name_layers[i])
      dat_slice<-data_test$feature_table[rownames(data_test$feature_table) %in% include_list$featureID, ]
      
      ################################
      # Per-layer variance filtering #
      ################################
      
      sds<-apply(dat_slice, 1, sd)
      dat_slice_filtered<-dat_slice[which(sds>quantile(sds, .95)),]
      keep_features<-c(keep_features, rownames(dat_slice_filtered))
      rm(dat_slice); rm(dat_slice_filtered)
    }

    #########################
    # Final filtered tables #
    #########################

    data_test$feature_table<-data_test$feature_table[keep_features, ]
    data_test$feature_metadata<-data_test$feature_metadata[keep_features,]


    ################
    # Sanity check #
    ################

    all(rownames(data_test$feature_table)==rownames(data_test$feature_metadata))

    ## [1] TRUE

    all(colnames(data_test$feature_table)==rownames(data_test$sample_metadata))

    ## [1] TRUE

# BayesCoop Implementation

    set.seed(1)
    result = bayesCoop(data_train, data_test, family = "gaussian", 
                       ss = c(0.05, 1), group = TRUE,
                       alpha_dirich = 1, maxit = 100, 
                       bbiters = 1100, bbburn = 1000,
                       abd_thresh = 0, prev_thresh = 0.1,
                       Warning = TRUE, verbose = TRUE)

    result$mspe

    ## [1] 672.6159

    result$time

    ## [1] 0.061

Let’s check MCMC convergence of the reciprocal Bayesian LASSO estimator
through two visualizations: trace plots and histograms.

    ######################################
    # Visualization of Posterior Samples #
    ######################################

    ##############
    # Trace Plot #
    ##############

    library(coda)
    par(mar=c(2,2,2,2))
    plot(mcmc(result$beta_samples[,1:9]),density=FALSE,smooth=TRUE)
```
![](https://github.com/himelmallick/BayesCOOP/blob/master/misc/unnamed-chunk-13-1.png)

```         
#############
# Histogram #
#############

library(psych)
multi.hist(result$beta_samples[,1:9],density=TRUE,main="")
```

![](https://github.com/himelmallick/BayesCOOP/blob/master/misc/unnamed-chunk-14-1.png)
