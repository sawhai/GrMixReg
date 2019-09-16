# Grouped Mixture of Regressions (GMR)

This repository contains the code to run the simulation study used in the paper "Performance Evaluation of Automotive
Dealerships using Grouped Mixture of Regressions".

The paper introduces Grouped Mixture of Regressions (GMR) as an extension of the regular Finite Mixture of Regressions (FMR). GMR allows for known group structure among observations, in addition to their posibly unknown heterogeneity. The fitting algorithm assigns groups (instead of individual observations) to clusters. An Expectation Maximization (EM) algorithm is developed  to fit this model. In addition, the Maximum a Posterior (MAP) prediction for new observations is derived and implemented. 

The code can be used to recreate the simulation study in the paper, as detailed below. The easiest way is perhaps to open  the project file `GrMixReg.Rproj` in RStudio and then run each of the simulation files below.

## Simulation files
The following are the main simulation files:

* `var_sig_level_sims.R` runs the simulations with variable noise level.
* `var_grp_size_sims.R` runs the simulations with variable numer of groups (G).
* `find_optimal_k.R` runs the simulations to find the optimal K (number of clusters) using cross-validation.

The following use results form simulations ran elsewhere to plot comparison plots:

* `plot_predict_comparison.R`plots the comparison between MMCL++ and GMR, using the 

The .RData files in `data/` folder are the results saved and used for plotting.

## Modules
The following modules are used:

* `modules/sims.R` contains functions executing the main simulation loops
* `modules/fit_GMR.R` contains functions for fitting and predication based on GMR.
* `modules/GMR_data_gen.R` contains functions to generate data from the GMR model.
* `modeuls/network_commons.R` contains functions for computing confusion matrix and NMI for cluster comparisons.




