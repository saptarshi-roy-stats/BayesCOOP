# BayesMultiview

The repository contains the analysis and simulation codes from the Bayesian multiview project.

## Dependencies

`BayesCOOP` requires the following `R` package: `devtools` (for installation only). Please install it before installing `BayesCOOP`, which can be done as follows (execute from within a fresh R session):

install.packages("devtools")
library(devtools)

## Installation

Once the dependencies are installed, `BayesCOOP` can be loaded using the following command:

devtools::install_github("himelmallick/BayesCOOP")
library(BayesCOOP)

## Load libraries 

library(BhGLM)
library(tidyverse)
library(caret)

## sourcing required functions 
rm(list = ls())
source("https://github.com/himelmallick/BayesCOOP/tree/master/R/bayesCoop.R")

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
                   alpha_dirich = 1, maxit = 100, 
                   bbiters = 1100, bbburn = 100,
                   abd_thresh = 0, prev_thresh = 0.1,
                   Warning = TRUE, verbose = TRUE)

