regularisation-path-following
=============================

Regularisation-path algorithm using patterns of input variables (itemsets) as covariates to fit a Cox Proportional Hazards regression model.

Full description of method and algorithm published in [Oxford Journal's Bioinformatics](http://bioinformatics.oxfordjournals.org/content/early/2013/10/04/bioinformatics.btt532.short).

Based on work by [Park & Hastie (2006)](http://www.stanford.edu/~hastie/Papers/glmpath.jrssb.pdf): *L1 Regularization Path Algorithm for Generalized Linear Models*.

To build necessary components, please run **build.sh** in a command shell.

Usage examples can be found in **test.R**. Functions to load test data files are located in **load.data.R**.

Important functions include: `itemset.coxpath` (single iteration of the path-following algorithm) and `run.cv` (cross-validated iterations). 
