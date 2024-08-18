# README

# TMBocelot

### TMBocelot: an Omnibus statistical Control modeL Optimizing the TMB Thresholds with systematic measurement errors

## Introduction

This codebase  focuses on validating the performance of statistical model frameworks under various error scenarios.

## Key Components

This code repository primarily includes four simulation experiments for Bayesian estimators, where each simulation experiment encompasses the parameter estimation process as well as the threshold calculation process through TMB-cat and the evaluation of result. In addition, there are also some other codes to assist with inference.

### Bayesian Joint Error Estimator

- **Code**: `pairError/joint.R`
- **Description**: This code simulates the performance of a Bayesian joint error estimator, which considers pair-wise error.

### Bayesian Single Error Estimator

- **Code**: `OnlyError/OnlyTMBError`.R
- **Description**: This code simulates the performance of a Bayesian single error estimator, only considering TMB errors.
- **Code**: `OnlyMis/OnlyMis.R`
- **Description**: This code simulates the performance of a Bayesian single error estimator, only considering the misclassification of objective tumor response.

### Bayesian Naive Estimator

- **Code**: `Naive/Naive.R`
- **Description**: This code simulates the performance of a Bayesian naïve estimator ignoring errors. This serves as a baseline for the joint and single error estimators.

### Other

- **Code**: `pairError/TMBDensity.R` and `pairError/ErrorDensity.R`
- **Description**: These two codes present the approaches based on Metropolis-Hasting algorithm for estimating Dirichlet mixture models.
- **Code**: `pairError/TMB-cat.R`
- **Description**: The code for calculating thresholds.
- **Code**: `pairError/Median_thres.R`
- **Description**: The code for calculating median-thresholds.
- **Code**: `Generate.R`
- **Description**: The code for generating simulation data.

## Note

Users can modify the data parameters to suit their specific research needs or to test the robustness of the models under different conditions.

