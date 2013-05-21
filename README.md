regularisation-path-following
=============================

Regularisation-path algorithm using pattern of input variables as covariates to fit a Cox Proportional Hazards regression model.

To build necessary components, please run **build.sh** in a command shell.

Usage examples can be found in **test.R**. Functions to load test data files are located in **load.data.R**.

Important functions include: `itemset.coxpath` (single iteration of the path-following algorithm) and `run.cv` (cross-validated iterations). 
