# BayesCOOP - Bayesian Cooperative Learning for Multimodal Integration

The repository houses the **`BayesCOOP`** R package for multimodal
integration for prediction of continuous outcomes.

## Dependencies

`BayesCOOP` requires the following `R` package: `devtools` (for
installation only). Please install it before installing `BayesCOOP`,
which can be done as follows (execute from within a fresh R session):

    install.packages("devtools")
    #> 
    #> The downloaded binary packages are in
    #>  /var/folders/c8/vl4b3y8s66g_s23mlmvpzvj80000gn/T//RtmpTQcgpe/downloaded_packages
    library(devtools)
    #> Loading required package: usethis

## Installation

Once the dependencies are installed, `BayesCOOP` can be loaded using the
following command:

    # If the package isn't installed, load it from source (package root) for knitting
    if (!requireNamespace("BayesCOOP", quietly = TRUE)) {
      if (requireNamespace("pkgload", quietly = TRUE)) {
        # load the current package source without attaching everything to the search path
        pkgload::load_all(export_all = FALSE, helpers = FALSE, attach_testthat = FALSE)
      }
    } else {
      devtools::install_github("himelmallick/BayesCOOP")
      library(BayesCOOP)
    }
    #> ℹ Loading BayesCOOP

## Example Implementation

## Loading the StelzerDOS real dataset

    data_train = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS.RData"))); rm(pcl)
    data_test = get(load(url("https://raw.githubusercontent.com/himelmallick/IntegratedLearner/master/data/StelzerDOS_valid.RData"))); rm(pcl)

## Pre-processing the longitudinal data by considering only baseline observations

### Remove metabolomics from the train set to match with validation

    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr", repos = "https://cloud.r-project.org")
    }
    library(dplyr)

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
